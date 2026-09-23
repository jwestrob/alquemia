#!/usr/bin/env python3
"""Fixed native-GFN solvent derivatives after two own-seed continuations."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import (HA_TO_KCAL, InvalidArtifact, read_json, record, verify,
                               write_new, xyz, paired)
from mace_site_kinematics import Kinematics
from structure_informed_starts import sources, role_for_mode
from accommodation_torsion_profiles import geometry_check
from compact_solvation import completed, diagnostics
from run_orca_task_manifest import load_manifest_tasks

PROTOCOL='native_GFN2_self_continued_solvent_force_check_v1'
CASES=('4MAE','q88jh5-pqq-la_model')
MODES={'4MAE':'A/172/chi3','q88jh5-pqq-la_model':'A/221/chi3'}
OFFSETS=(-.001,-.0005,.0005,.001)
SETTINGS={'step_radian':.001,'maximum_heavy_displacement_A':.8,
          'derivative_absolute_tolerance_kcal_mol_radian':.2,'derivative_relative_tolerance':.05,
          'cell_energy_tolerance_kcal_mol':.1,'solvent_contrast_tolerance_kcal_mol':.2,
          'moving_heavy_jacobian_threshold_A_radian':1e-10,
          'MaxIter':500,'temperature_K':300,'mpi_ranks':8,'workers':8,'maxcore_MiB':2000,
          'reported_stage':2,'gradient_quantity':'gradient_not_force'}


def recipe(charge,multiplicity,medium,stage):
    from native_pool_continuation import recipe as restart_recipe
    value=restart_recipe(charge,multiplicity,medium,gradient=False)
    if stage==0:value=value.replace('! Native-GFN2-xTB','! Native-GFN2-xTB NoAutostart',1)
    return value


def pinned_sources(inputs,seeds):
    _,srcs,_=sources(inputs);selected={s['case_id']:s for s in srcs if s['case_id'] in CASES}
    sd=read_json(seeds);rows=[r for r in sd['rows'] if r['case_id'] in CASES and r['candidate']=='origin']
    if len(rows)!=8:raise InvalidArtifact('exact eight original seed-source cells required')
    for cid in CASES:
        s=selected[cid];ca,la=s['tasks'];paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
        for t in s['tasks']:
            for medium in ('vacuum','alpb'):
                r=next(r for r in rows if (r['case_id'],r['metal'],r['medium'])==(cid,t['metal'],medium))
                if xyz(verify(r['xyz']))!=xyz(verify(t['xyz'])) or (r['charge'],r['multiplicity'])!=(t['charge'],t['multiplicity']):
                    raise InvalidArtifact('current common-pool origin differs from source mapping/state')
    return selected,rows


def prepare(inputs,seeds,agreement,output):
    srcs,seedrows=pinned_sources(inputs,seeds);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        f=impl/p.name;shutil.copyfile(p,f);pins[p.name]=record(f)
    cells=[];maps={}
    for cid in CASES:
        src=srcs[cid]
        for t in src['tasks']:
            kin=Kinematics(read_json(verify(t['mapping']))['context']);idx=next(i for i,m in enumerate(kin.modes) if m['id']==MODES[cid])
            role,res,chi,_=role_for_mode(src,kin,kin.modes[idx])
            if (role,res,chi)!=('anchor_glutamate','GLU','chi3'):raise InvalidArtifact('terminal anchor Glu identity differs')
            zeros=np.zeros(len(kin.modes));p0,c0,j0,cj0=kin.evaluate(zeros)
            original=xyz(verify(t['xyz']));symbols=[a[0] for a in original]
            if not np.allclose(c0,[a[1:] for a in original],atol=1e-12,rtol=0):raise InvalidArtifact('mapping q0 replay differs')
            heavy=np.array(kin.heavy);speed=np.linalg.norm(j0[idx,heavy],axis=1);moving=heavy[speed>SETTINGS['moving_heavy_jacobian_threshold_A_radian']]
            if not len(moving):raise InvalidArtifact('no physical heavy movement')
            meta={'source_task':t,'physical_mapping':t['mapping'],'mode_id':MODES[cid],'mode_index':idx,
                  'source_role':src['roles'][role],'adaptive_selected':idx in t['active_indices'],
                  'moving_heavy_indices':moving.tolist(),'moving_heavy_ids':[kin.data['physical_ids'][i] for i in moving],
                  'RMS_moving_heavy_speed_A_radian':float(np.sqrt(np.mean(np.sum(j0[idx,moving]**2,axis=1)))),
                  'RMS_all_heavy_speed_A_radian':float(np.sqrt(np.mean(np.sum(j0[idx,heavy]**2,axis=1))))}
            maps[cid+'__'+t['metal']]=meta
            for offset in OFFSETS:
                q=zeros.copy();q[idx]=offset;physical,coords,_,_=kin.evaluate(q);checks=geometry_check(kin,q,symbols)
                if not checks['pass'] or checks['maximum_heavy_displacement_A']>SETTINGS['maximum_heavy_displacement_A']:
                    raise InvalidArtifact('physical displacement check failed')
                coords=[(a[0],*p) for a,p in zip(original,coords)]
                ident=cid+'__'+t['metal']+'__'+str(offset).replace('-','m').replace('.','p')
                xp=out/'coordinates'/(ident+'.xyz');xp.parent.mkdir(exist_ok=True)
                xp.write_text(str(len(coords))+'\n'+PROTOCOL+'\n'+''.join(a[0]+' '+' '.join(format(v,'.16g') for v in a[1:])+'\n' for a in coords))
                for medium in ('vacuum','alpb'):
                    seed=next(r for r in seedrows if (r['case_id'],r['metal'],r['medium'])==(cid,t['metal'],medium))
                    cells.append({'cell_id':ident+'__'+medium,'case_id':cid,'metal':t['metal'],'medium':medium,
                        'offset_radian':offset,'full_q':q.tolist(),'mode_id':MODES[cid],'physical_mapping':t['mapping'],
                        'charge':t['charge'],'multiplicity':t['multiplicity'],'xyz':record(xp),
                        'geometry_checks':checks,'source':seed,'source_task':t})
    receipt=read_json(verify(seedrows[0]['receipt']))
    d={'protocol_id':PROTOCOL,'settings':SETTINGS,'inputs':record(inputs),'seeds':record(seeds),'agreement':record(agreement),
       'implementation':pins,'orca':receipt['orca_executable'],'cells':cells,'maps':maps,'new_molecular_calls':96,
       'external_q0_cells':8,'external_q0_owner':'native_pool_continuation_20260923','new_DFT_MACE_calls':0,
       'baseline_changed':False}
    write_new(out/'design.json',d)
    return prepare_stage(out/'design.json',0,None)


def prepare_stage(design,stage,previous):
    d=read_json(design)
    if stage not in (0,1,2) or d['settings']!=SETTINGS:raise InvalidArtifact('fixed stage/settings differ')
    prior=read_json(previous) if previous else None
    if bool(prior)!=(stage>0):raise InvalidArtifact('preceding own stage required')
    if prior and (prior['design']!=record(design) or prior['stage']!=stage-1):raise InvalidArtifact('stage sequence differs')
    directory=Path(design).parent/('stage'+str(stage));directory.mkdir(exist_ok=False);tasks=[];blocked=[]
    for cell in d['cells']:
        source=cell['source'];seed=None
        if prior:
            row=next(r for r in prior['rows'] if r['cell_id']==cell['cell_id'])
            if row['status']!='complete':blocked.append({'cell_id':cell['cell_id'],'reason':'preceding cell unavailable','previous':row});continue
            seed=row['seed_after'];source={**source,'energy_hartree':row['actual']['energy_hartree'],
                'xtbw':seed['preserved_after'],'gbw':seed['preserved_gbw_after']}
        tid=cell['cell_id']+'__stage'+str(stage);td=directory/'tasks'/tid;td.mkdir(parents=True)
        shutil.copyfile(verify(cell['xyz']),td/'core.xyz');(td/'endpoint.inp').write_text(recipe(cell['charge'],cell['multiplicity'],cell['medium'],stage))
        task={**cell,'task_id':tid,'stage':stage,'candidate':'displacement','source':source,
              'input':record(td/'endpoint.inp'),'xyz':record(td/'core.xyz'),'output_path':str(td/'endpoint.out'),
              'gradient_requested':False,'restart':bool(stage),'active_seed_path':str(td/'endpoint.runtime.xtbw'),
              'active_gbw_path':str(td/'endpoint.runtime.gbw')}
        if seed:
            shutil.copyfile(verify(seed['preserved_after']),td/'seed.immutable.xtbw');shutil.copyfile(verify(seed['preserved_gbw_after']),td/'seed.immutable.gbw')
            task.update(immutable_seed=record(td/'seed.immutable.xtbw'),immutable_gbw=record(td/'seed.immutable.gbw'))
            task['seed_source']={'xtbw':seed['preserved_after'],'gbw':seed['preserved_gbw_after']}
        tasks.append(task)
    m={'protocol_id':PROTOCOL,'settings':SETTINGS,'stage':stage,'design':record(design),'previous':record(previous) if previous else None,
       'tasks':tasks,'blocked':blocked,'orca':d['orca'],'implementation':d['implementation'],
       'execution_resources':{'mpi_ranks':8,'concurrent_tasks':8},
       'execution_policy':{'task_runner':d['implementation']['run_orca_task_manifest.py'],
                           'runtime_renderer':d['implementation']['render_orca_runtime_input.py']}}
    mp=directory/'manifest.json';write_new(mp,m);result=validate(mp,True);write_new(directory/'PREFLIGHT.json',result);return result


def validate(manifest,fresh=False):
    m=read_json(manifest);d=read_json(verify(m['design']))
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or len(d['cells'])!=32:raise InvalidArtifact('force scope changed')
    for k in ('inputs','seeds','agreement','orca'):verify(d[k])
    for pin in d['implementation'].values():verify(pin)
    cells={c['cell_id']:c for c in d['cells']}
    if len(m['tasks'])+len(m['blocked'])!=32:raise InvalidArtifact('fixed denominator changed')
    if {t['cell_id'] for t in m['tasks']}|{t['cell_id'] for t in m['blocked']}!=set(cells):raise InvalidArtifact('cell membership differs')
    for t in m['tasks']:
        c=cells[t['cell_id']]
        for key in ('case_id','metal','medium','offset_radian','full_q','mode_id','physical_mapping','charge','multiplicity'):
            if t[key]!=c[key]:raise InvalidArtifact('state or coordinate changed')
        if xyz(verify(t['xyz']))!=xyz(verify(c['xyz'])) or verify(t['input']).read_text()!=recipe(t['charge'],t['multiplicity'],t['medium'],m['stage']):
            raise InvalidArtifact('input recipe or coordinates changed')
        kin=Kinematics(read_json(verify(t['physical_mapping']))['context'])
        if not np.allclose(kin.evaluate(t['full_q'])[1],[a[1:] for a in xyz(verify(t['xyz']))],atol=1e-12,rtol=0):raise InvalidArtifact('physical map changed')
        expected={'core.xyz','endpoint.inp'}
        if m['stage']:
            prior=read_json(verify(m['previous']));r=next(r for r in prior['rows'] if r['cell_id']==t['cell_id'])
            if r['status']!='complete' or verify(t['immutable_seed']).read_bytes()!=verify(r['seed_after']['preserved_after']).read_bytes() or verify(t['immutable_gbw']).read_bytes()!=verify(r['seed_after']['preserved_gbw_after']).read_bytes():raise InvalidArtifact('own preceding seed differs')
            expected|={'seed.immutable.xtbw','seed.immutable.gbw'}
        if fresh and {p.name for p in Path(t['output_path']).parent.iterdir()}!=expected:raise InvalidArtifact('nonfresh task folder')
    _,tasks=load_manifest_tasks(Path(manifest))
    return {'status':'validated','manifest':record(manifest),'stage':m['stage'],'tasks':len(tasks),'blocked':len(m['blocked']),'new_calls_in_validation':0}


def collect(manifest,output):
    from native_pool_continuation import collect_task
    validate(manifest);m=read_json(manifest);mp=Path(manifest);rows=[]
    beforep=mp.parent/'SEEDS_BEFORE.json';afterp=mp.parent/'SEEDS_AFTER.json'
    before={r['task_id']:r for r in read_json(beforep)['rows']} if beforep.exists() else {}
    after={r['task_id']:r for r in read_json(afterp)['rows']} if afterp.exists() else {}
    for t in m['tasks']:
        if m['stage']:
            r=collect_task(manifest,t,before.get(t['task_id']),after.get(t['task_id']))
            # Shared continuation parser uses confirmed_restart; normalize availability only.
            if r['status']=='confirmed_restart':r['status']='complete'
        else:
            r={'status':'unavailable','actual':None}
            pin=completed(manifest,t['task_id'])
            if pin:
                try:
                    from structure_informed_starts import scf_details
                    audit=diagnostics(pin,t);text=verify(pin['output']).read_text();detail=scf_details(text)
                    if ('INITIAL GUESS: SAD' not in text or not detail['native_mixer_observed'] or
                        audit['charge_sanity_status']!='pass' or audit['parameter_export']['sha256']!=t['source']['parameter_export']['sha256']):raise InvalidArtifact('initial native/state audit failed')
                    a=after[t['task_id']];verify(a['preserved_after']);verify(a['preserved_gbw_after'])
                    r.update(status='complete',actual=pin,audit=audit,details=detail,seed_after=a)
                except Exception as exc:r.update(status='audit_failed',reason=str(exc),actual=pin)
        r.update({k:t[k] for k in ('cell_id','task_id','case_id','metal','medium','offset_radian','stage')});rows.append(r)
    for b in m['blocked']:rows.append({'cell_id':b['cell_id'],'status':'blocked','reason':b['reason'],'actual':None})
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'design':m['design'],'stage':m['stage'],'rows':rows,
            'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete'}
    write_new(output,result);return result


def execute(manifest):
    from native_pool_continuation import execute_tasks
    return execute_tasks(manifest,validate,seeded=read_json(manifest)['stage']>0)


def run(design):
    design=Path(design).resolve();stages=[]
    for stage in range(3):
        mp=design.parent/('stage'+str(stage))/'manifest.json'
        if stage:prepare_stage(design,stage,stages[-1])
        error=None
        try:execute(mp)
        except Exception as exc:error=f'{type(exc).__name__}: {exc}'
        cp=mp.parent/'collection.json';collect(mp,cp);stages.append(cp)
        if error:write_new(mp.parent/'EXECUTION_ERROR.json',{'error':error})
    result={'design':record(design),'collections':[record(p) for p in stages],'new_DFT_MACE_calls':0}
    write_new(design.parent/'RUN.json',result);return result


def derivative(energies):
    return {'h':(energies[.001]-energies[-.001])*HA_TO_KCAL/.002,
            'half':(energies[.0005]-energies[-.0005])*HA_TO_KCAL/.001}


def derivative_check(coarse,fine,analytic,speed):
    tolerance=max(.2,.05*abs(fine));refinement=abs(coarse-fine);error=abs(analytic-fine)
    return {'centered_h_kcal_mol_radian':coarse,'centered_half_h_kcal_mol_radian':fine,
            'analytic_kcal_mol_radian':analytic,'refinement_error_kcal_mol_radian':refinement,
            'analytic_error_kcal_mol_radian':error,'acceptance_tolerance_kcal_mol_radian':tolerance,
            'refinement_pass':refinement<=tolerance,'analytic_pass':error<=tolerance,
            'pass':refinement<=tolerance and error<=tolerance,
            'physical_RMS_speed_A_radian':speed,
            'centered_half_h_kcal_mol_A':fine/speed,'analytic_kcal_mol_A':analytic/speed,
            'analytic_error_kcal_mol_A':error/speed}


def compare(design,centers,output):
    """Compare fixed stage2 energies with external q0 gradients; no engine calls."""
    from native_pool_continuation import gradient_result, recipe as center_recipe
    d=read_json(design);root=Path(design).parent
    stages=[read_json(root/('stage'+str(i))/'collection.json') for i in (0,1,2)]
    for i,c in enumerate(stages):
        if c['design']!=record(design) or c['stage']!=i:raise InvalidArtifact('displacement stage identity differs')
        verify(c['manifest'])
    center2=read_json(centers);cm2=read_json(verify(center2['manifest']))
    if center2['stage']!=2:raise InvalidArtifact('external center must use stage2')
    center1=read_json(verify(cm2['previous']));by_stage=[];origins={};origin_rows=[]
    if center1['stage']!=1:raise InvalidArtifact('external center previous stage differs')
    for c in stages:by_stage.append({r['cell_id']:r for r in c['rows']})
    for cid in CASES:
        for metal in ('Ca','La'):
            meta=d['maps'][cid+'__'+metal];kin=Kinematics(read_json(verify(meta['physical_mapping']))['context'])
            jac=kin.evaluate(np.zeros(len(kin.modes)))[3][meta['mode_index']]
            for medium in ('vacuum','alpb'):
                rows=[next(r for r in c['rows'] if (r['case_id'],r['candidate'],r['metal'],r['medium'])==(cid,'origin',metal,medium)) for c in (center1,center2)]
                key=(cid,metal,medium);origin_rows.append({'case_id':cid,'metal':metal,'medium':medium,'stage1':rows[0],'stage2':rows[1]})
                if any(r['status']!='confirmed_restart' for r in rows):continue
                t=next(t for t in cm2['tasks'] if t['task_id']==rows[1]['task_id'])
                if xyz(verify(t['xyz']))!=xyz(verify(meta['source_task']['xyz'])) or (t['charge'],t['multiplicity'])!=(meta['source_task']['charge'],meta['source_task']['multiplicity']):
                    raise InvalidArtifact('external center source/state mismatch')
                if verify(t['input']).read_text()!=center_recipe(t['charge'],t['multiplicity'],medium,True):raise InvalidArtifact('q0 native EnGrad recipe differs')
                pin=completed(verify(center2['manifest']),t['task_id'])
                if pin!=rows[1]['actual']:raise InvalidArtifact('external center receipt mismatch')
                g=gradient_result(t,pin)
                origins[key]={'energy_hartree':pin['energy_hartree'],
                    'analytic_kcal_mol_radian':float(np.sum(jac*np.asarray(g['gradient_kcal_mol_A']))),
                    'gradient':g,'stage1_energy_hartree':rows[0]['energy_hartree'],
                    'energy_change_kcal_mol':(pin['energy_hartree']-rows[0]['energy_hartree'])*HA_TO_KCAL}
    cells=[];case_results=[]
    for c in d['cells']:
        rr=[idx[c['cell_id']] for idx in by_stage]
        good=all(r['status']=='complete' for r in rr)
        values=[r['actual']['energy_hartree'] if r['status']=='complete' else None for r in rr]
        delta=(values[2]-values[1])*HA_TO_KCAL if values[1] is not None and values[2] is not None else None
        cells.append({**{k:c[k] for k in ('cell_id','case_id','metal','medium','offset_radian')},
            'status':'complete' if good else 'unavailable','energies_hartree_by_stage':values,
            'stage2_minus_stage1_kcal_mol':delta,'numerical_pass':delta is not None and abs(delta)<=.1})
    for cid in CASES:
        measures=[];quantities={};settling=[];center_checks=[]
        for metal in ('Ca','La'):
            for medium in ('vacuum','alpb'):
                key=(cid,metal,medium);cc=[c for c in cells if (c['case_id'],c['metal'],c['medium'])==key]
                source=origins.get(key)
                center_checks.append({'metal':metal,'medium':medium,'status':'complete' if source else 'unavailable',
                    'stage2_minus_stage1_kcal_mol':source['energy_change_kcal_mol'] if source else None,
                    'pass':source is not None and abs(source['energy_change_kcal_mol'])<=.1})
                name=metal+'_'+medium
                if source is None or len(cc)!=4 or any(c['status']!='complete' for c in cc):
                    quantities[name]=None;continue
                e={c['offset_radian']:c['energies_hartree_by_stage'][2] for c in cc}
                quantities[name]={'energies_hartree':e,'q0_energy_hartree':source['energy_hartree'],
                    'gradient':source['gradient'],'analytic':source['analytic_kcal_mol_radian'],**derivative(e)}
        for metal in ('Ca','La'):
            a,b=[quantities[metal+'_'+s] for s in ('alpb','vacuum')]
            quantities[metal+'_solvent']=({field:a[field]-b[field] for field in ('h','half','analytic')} if a and b else None)
        a,b=[quantities[z+'_solvent'] for z in ('Ca','La')]
        quantities['Ca_minus_La_solvent']=({field:a[field]-b[field] for field in ('h','half','analytic')} if a and b else None)
        for offset in (*OFFSETS,0.):
            deltas={}
            for metal in ('Ca','La'):
                for medium in ('vacuum','alpb'):
                    if offset==0.:
                        source=origins.get((cid,metal,medium));delta=source['energy_change_kcal_mol'] if source else None
                    else:delta=next(c['stage2_minus_stage1_kcal_mol'] for c in cells if (c['case_id'],c['metal'],c['medium'],c['offset_radian'])==(cid,metal,medium,offset))
                    deltas[metal+'_'+medium]=delta
            complete=all(v is not None for v in deltas.values())
            delta=(deltas['Ca_alpb']-deltas['Ca_vacuum']-deltas['La_alpb']+deltas['La_vacuum']) if complete else None
            settling.append({'offset_radian':offset,'component_changes_kcal_mol':deltas,
                'solvent_contrast_change_kcal_mol':delta,'pass':delta is not None and abs(delta)<=.2})
        speed=d['maps'][cid+'__Ca']['RMS_moving_heavy_speed_A_radian']
        if speed!=d['maps'][cid+'__La']['RMS_moving_heavy_speed_A_radian']:raise InvalidArtifact('paired physical coordinate measure differs')
        for name,q in quantities.items():
            row={'quantity':name,'status':'available' if q else 'unavailable','raw':q}
            if q:row.update(derivative_check(q['h'],q['half'],q['analytic'],speed))
            else:row['pass']=False
            measures.append(row)
        numeric=all(c['numerical_pass'] for c in cells if c['case_id']==cid) and all(c['pass'] for c in center_checks+settling)
        force=all(q['pass'] for q in measures)
        case_results.append({'case_id':cid,'numerical_pass':numeric,'force_pass':force,
            'qualified_for_this_diagnostic':numeric and force,
            'qualified_solvent_contrast_gradient_kcal_mol_radian':quantities['Ca_minus_La_solvent']['analytic'] if numeric and force else None,
            'derivatives':measures,
            'center_settling':center_checks,'solvent_contrast_settling':settling})
    result={'protocol_id':PROTOCOL,'design':record(design),'analysis_implementation':record(__file__),
        'collections':[record(root/('stage'+str(i))/'collection.json') for i in (0,1,2)],'external_centers':record(centers),
        'case_denominator':2,'displaced_cell_denominator':32,'external_center_denominator':8,
        'derivative_denominator':14,'cases':case_results,'cells':cells,'center_rows':origin_rows,
        'qualified_cases':sum(c['qualified_for_this_diagnostic'] for c in case_results),
        'passed_derivative_quantities':sum(r['pass'] for c in case_results for r in c['derivatives']),
        'new_molecular_calls_in_comparison':0,'classifier_gain':None,'entropy_correction':None,'production_changed':False}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('inputs','seeds','agreement','output'):a.add_argument('--'+k,required=True)
    a=s.add_parser('stage')
    a.add_argument('--design',required=True);a.add_argument('--stage',type=int,required=True);a.add_argument('--previous')
    for name in ('validate','execute','collect'):
        a=s.add_parser(name);a.add_argument('--manifest',required=True)
        if name=='collect':a.add_argument('--output',required=True)
    a=s.add_parser('run');a.add_argument('--design',required=True)
    a=s.add_parser('compare')
    for name in ('design','centers','output'):a.add_argument('--'+name,required=True)
    args=vars(p.parse_args());op=args.pop('op');fn=prepare_stage if op=='stage' else globals()[op]
    result=fn(**args);print(json.dumps({k:v for k,v in result.items() if k not in ('rows','cases','cells','center_rows')},indent=2))
if __name__=='__main__':main()
