"""Physical derivatives of the frozen OMOL + native GFN2 solvent descriptor."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import re
import shutil
import numpy as np
from affordable_common import BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new, xyz
from affordable_response import read_engrad
from affordable_workflow import dry_run, execute as execute_manifest
from compact_solvation import completed, diagnostics
from mace_hybrid import EV_TO_KCAL
from mace_site_kinematics import Kinematics

PROTOCOL='compact_solvent_physical_response_v1'
CASES=('1H4I','4MAE','1F6S','6IP9','1GLG','2FW0','2FVY')
MODES={'1H4I':'A/177/chi3','1GLG':'A/140/peptide_crankshaft'}
STEPS=(-.01,-.005,.005,.01)


def recipe(charge,medium,gradient):
    return ('! Native-GFN2-xTB NoAutostart'+(' EnGrad' if gradient else '')
            +(' ALPB(Water)' if medium=='alpb' else '')
            +'\n%maxcore 2000\n%method\n WriteXTBParam true\n ReadXTBParam false\nend\n'
            +'%scf\n SmearTemp 300\n Convergence Tight\n UseXTBMixer true\nend\n'
            +f'* xyzfile {charge} 1 core.xyz\n')


def prepare(source,agreement,output):
    from second_shell_context import parent_state
    from coordination_preparation_context import geometry
    from hydration_network import write_xyz
    src=read_json(source);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    sources={};maps={};all_tasks=[]
    for case in CASES:
        originals=[t for t in src['all_tasks'] if t['case_id']==case and t['representation']=='context']
        if len(originals)!=4:raise InvalidArtifact('missing source endpoint quartet')
        mm=verify(originals[0]['source_endpoint']['source_manifest'])
        physical_manifest=mm.parent/'manifest.json';sm=read_json(physical_manifest)
        cfg=read_json(verify(sm['configuration']));item=next(s for s in sm['states'] if s['case']==case)
        st=parent_state(item['source'],cfg['topology']);prep=read_json(verify(item['preparation']))
        for metal in ('Ca','La'):
            t=next(t for t in originals if t['metal']==metal and t['medium']=='vacuum')
            source_ep=t['source_endpoint'];r=read_json(verify(source_ep['native_MACE_receipt']))
            fm=read_json(verify(r['manifest']));ft=next(x for x in fm['tasks'] if x['task_id']==r['task_id'])
            rows=xyz(verify(t['xyz']));forces=np.load(verify(r['forces']),allow_pickle=False)
            if xyz(verify(ft['xyz']))!=rows or forces.shape!=(len(rows),3) or not np.isfinite(forces).all():
                raise InvalidArtifact('archived force geometry/order differs')
            if r['status']!='computed' or ft['charge']!=t['charge'] or ft['spin_multiplicity']!=1:
                raise InvalidArtifact('archived force state differs')
            gp=geometry(st,prep,rows,xyz(verify(item['source']['endpoints'][metal]['xyz'])))
            kin=Kinematics(gp['context']);zero=np.zeros(len(kin.modes));physical,core,_,jac=kin.evaluate(zero)
            if not np.allclose(core,[a[1:] for a in rows],atol=1e-10,rtol=0):raise InvalidArtifact('mapping origin differs')
            name=case+'__'+metal;mp=out/'maps'/(name+'.json');write_new(mp,gp)
            maps[name]=record(mp);sources[name]={'case_id':case,'metal':metal,'xyz':t['xyz'],'charge':t['charge'],
                'native_MACE_receipt':source_ep['native_MACE_receipt'],'forces':r['forces'],
                'physical_preparation':item['preparation'],'physical_manifest':record(physical_manifest),
                'configuration':sm['configuration'],'mode_ids':[x['id'] for x in kin.modes],
                'mode_units':[x['unit'] for x in kin.modes],
                'MACE_projected_gradient_kcal_per_unit':(-EV_TO_KCAL*np.einsum('mij,ij->m',jac,forces)).tolist()}
            candidates=[('center',0.,rows,None)]
            if case in MODES:
                idx=next(i for i,m in enumerate(kin.modes) if m['id']==MODES[case])
                for step in STEPS:
                    q=zero.copy();q[idx]=step;_,positions,_,_=kin.evaluate(q);checks=kin.check(q)
                    if not checks['pass'] or kin.displacement(q)>.05:raise InvalidArtifact('displacement geometry check failed')
                    candidates.append(('displacement',step,[(a[0],*p) for a,p in zip(rows,positions)],MODES[case]))
            for medium in ('vacuum','alpb'):
                parent=next(t for t in originals if t['metal']==metal and t['medium']==medium)
                for kind,step,coords,mode_id in candidates:
                    tid=name+'__'+medium+'__'+kind+(('_'+str(step).replace('-','m').replace('.','p')) if step else '')
                    stage='gate' if case=='1H4I' and kind=='center' else 'followup'
                    directory=out/stage/'tasks'/tid;directory.mkdir(parents=True);xp=directory/'core.xyz';ip=directory/'endpoint.inp'
                    if kind=='center':shutil.copyfile(verify(parent['xyz']),xp)
                    else:write_xyz(xp,coords,PROTOCOL)
                    ip.write_text(recipe(parent['charge'],medium,kind=='center'))
                    all_tasks.append({'task_id':tid,'case_id':case,'case':case,'metal':metal,'medium':medium,'kind':kind,
                        'mode_id':mode_id,'step_radian':step,'charge':parent['charge'],'multiplicity':1,'xyz':record(xp),
                        'input':record(ip),'output_path':str(directory/'endpoint.out'),'source_task':parent,
                        'physical_mapping':maps[name]})
    gate=[t for t in all_tasks if t['case_id']=='1H4I' and t['kind']=='center']
    follow=[t for t in all_tasks if t not in gate]
    design={'protocol_id':PROTOCOL,'agreement':record(agreement),'source_manifest':record(source),'sources':sources,
        'maps':maps,'all_tasks':all_tasks,'gate_tasks':[t['task_id'] for t in gate],
        'validation_modes':MODES,'validation_steps_radian':list(STEPS),'new_gradient_calls':28,'new_SP_calls':32,
        'new_MACE_calls':0,'new_DFT_calls':0,'baseline_changed':False,'calibrated_score':None}
    write_new(out/'design.json',design)
    for stage,tasks in [('gate',gate),('followup',follow)]:
        dest=out/stage;impl=dest/'implementation';impl.mkdir(parents=True);pins={}
        for p in Path(__file__).parent.glob('*.py'):
            f=impl/p.name;shutil.copyfile(p,f);pins[p.name]=record(f)
        manifest={'protocol_id':PROTOCOL,'stage':stage,'design':record(out/'design.json'),'agreement':record(agreement),
            'tasks':tasks,'implementation':pins,'orca':src['orca'],
            'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
            'execution_resources':{'mpi_ranks':8,'concurrent_tasks':8},'baseline_changed':False}
        if stage=='followup':manifest['gate_manifest']=record(out/'gate/manifest.json')
        write_new(dest/'manifest.json',manifest)
    return {'status':'prepared','gate':validate(out/'gate/manifest.json'),'followup':validate(out/'followup/manifest.json')}


def validate(manifest):
    m=read_json(manifest);d=read_json(verify(m['design']))
    if m['protocol_id']!=PROTOCOL or d['validation_modes']!=MODES or tuple(d['validation_steps_radian'])!=STEPS:
        raise InvalidArtifact('declared physical experiment changed')
    verify(d['source_manifest']);verify(m['agreement'])
    for p in m['implementation'].values():verify(p)
    expected=[t for t in d['all_tasks'] if (t['task_id'] in d['gate_tasks'])==(m['stage']=='gate')]
    if m['tasks']!=expected:raise InvalidArtifact('fixed task denominator changed')
    for t in m['tasks']:
        parent=t['source_task'];pin=d['sources'][t['case_id']+'__'+t['metal']]
        verify(pin['native_MACE_receipt']);verify(pin['forces']);verify(pin['configuration']);verify(pin['physical_preparation'])
        if (t['charge'],t['multiplicity'])!=(parent['charge'],parent['multiplicity']):raise InvalidArtifact('state changed')
        if verify(t['input']).read_text()!=recipe(t['charge'],t['medium'],t['kind']=='center'):raise InvalidArtifact('method changed')
        kin=Kinematics(read_json(verify(t['physical_mapping']))['context']);q=np.zeros(len(kin.modes))
        if t['kind']=='center':
            if verify(t['xyz']).read_bytes()!=verify(parent['xyz']).read_bytes():raise InvalidArtifact('center geometry changed')
        else:
            q[next(i for i,x in enumerate(kin.modes) if x['id']==t['mode_id'])]=t['step_radian']
            if not kin.check(q)['pass'] or kin.displacement(q)>.05:raise InvalidArtifact('physical displacement invalid')
        expected_xyz=kin.evaluate(q)[1];actual=xyz(verify(t['xyz']))
        if not np.allclose(expected_xyz,[a[1:] for a in actual],atol=1e-9,rtol=0):raise InvalidArtifact('coordinate mapping mismatch')
    return dry_run(manifest)


def parse_gradient(t,pin):
    import gemmi
    op=verify(pin['output']);text=op.read_text();path=op.parent/'endpoint.runtime.engrad'
    if re.search(r'numerical gradient|numerical differentiation',text,re.I):raise InvalidArtifact('numerical gradient fallback')
    if 'XTB SCF gradient' not in text or 'XTB CN gradient' not in text:
        raise InvalidArtifact('native analytic gradient driver not confirmed')
    if ('ALPB gradient' in text)!=(t['medium']=='alpb'):
        raise InvalidArtifact('analytic solvent derivative does not match endpoint medium')
    g=read_engrad(path);atoms=xyz(verify(t['xyz']))
    if abs(g['energy_Ha']-pin['energy_hartree'])>1e-8 or len(atoms)!=g['atom_count']:raise InvalidArtifact('gradient energy/count mismatch')
    if not np.array_equal(g['atomic_numbers'],[gemmi.Element(a[0]).atomic_number for a in atoms]):raise InvalidArtifact('gradient atom order mismatch')
    if not np.allclose(g['coordinates_bohr']*BOHR_TO_A,[a[1:] for a in atoms],atol=1e-6,rtol=0):raise InvalidArtifact('gradient coordinates mismatch')
    gradient=g['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A
    kin=Kinematics(read_json(verify(t['physical_mapping']))['context']);jac=kin.evaluate(np.zeros(len(kin.modes)))[3]
    return {'engrad':record(path),'gradient_kcal_mol_A':gradient.tolist(),
            'projected_gradient_kcal_per_unit':np.einsum('mij,ij->m',jac,gradient).tolist(),
            'numerical_gradient_used':False,'quantity':'gradient_not_force'}


def collect(manifest):
    validate(manifest);m=read_json(manifest);rows=[];pairs=[]
    for t in m['tasks']:
        row={k:t[k] for k in ('task_id','case_id','metal','medium','kind','mode_id','step_radian')}
        row.update(status='unavailable',energy_hartree=None,reason=None)
        pin=completed(manifest,t['task_id']);op=Path(t['output_path'])
        if op.exists():row['output']=record(op)
        if pin:
            try:
                audit=diagnostics(pin,t);text=op.read_text()
                if audit['charge_sanity_status']!='pass':raise InvalidArtifact('atomic charge sanity failed')
                if 'INFO: Using special xTB SCF mixer' not in text:raise InvalidArtifact('native mixer not confirmed')
                tol=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text)
                if len(tol)!=1 or float(tol[0])>1e-8:raise InvalidArtifact('tight convergence not confirmed')
                source=Path(t['source_task']['output_path']).parent/'endpoint.runtime.xtb.json'
                # Reused primary tasks have their actual output pinned in the original
                # collection, so fall back to its source task receipt when necessary.
                if not source.exists():
                    d=read_json(verify(m['design']));parent=read_json(verify(d['source_manifest']))
                    reused=parent['reused'][t['source_task']['task_id']];source=verify(reused['output']).parent/'endpoint.runtime.xtb.json'
                if read_json(source)!=read_json(verify(audit['parameter_export'])):raise InvalidArtifact('GFN2 parameters changed')
                row.update(pin);row.update(audit)
                if t['kind']=='center':row.update(parse_gradient(t,pin))
                row['status']='complete'
            except (InvalidArtifact,OSError,ValueError,KeyError,IndexError) as exc:row.update(status='unsupported',reason=str(exc))
        rows.append(row)
    by={(r['case_id'],r['metal'],r['medium']):r for r in rows if r['kind']=='center'}
    for case in dict.fromkeys(r['case_id'] for r in rows if r['kind']=='center'):
        delta={};status='unavailable'
        for metal in ('Ca','La'):
            a,b=[by[(case,metal,x)] for x in ('alpb','vacuum')]
            if any(r['status']!='complete' for r in (a,b)):continue
            old=[]
            for medium in ('alpb','vacuum'):
                task=next(t for t in m['tasks'] if (t['case_id'],t['metal'],t['medium'],t['kind'])==(case,metal,medium,'center'))
                d=read_json(verify(m['design']));pm=read_json(verify(d['source_manifest']))
                original=pm['reused'].get(task['source_task']['task_id']) or completed(verify(d['source_manifest']),task['source_task']['task_id'])
                if original is None:raise InvalidArtifact('primary center result unavailable')
                old.append(original['energy_hartree'])
            delta[metal]=((a['energy_hartree']-b['energy_hartree'])-(old[0]-old[1]))*HA_TO_KCAL
        if len(delta)==2:status='pass' if max(abs(x) for x in delta.values())<=.1 and abs(delta['Ca']-delta['La'])<=.2 else 'fail'
        pairs.append({'case_id':case,'status':status,'transfer_change_kcal_mol':delta,
                      'score_change_kcal_mol':delta['Ca']-delta['La'] if len(delta)==2 else None})
    return {'manifest':record(manifest),'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete',
            'rows':rows,'center_energy_checks':pairs,'analytic_support_gate':bool(pairs) and all(r['status']=='pass' for r in pairs),
            'new_DFT_MACE_calls':0,'affinity_correction':None,'entropy':None,'baseline_changed':False}


def execute(manifest):
    validate(manifest);m=read_json(manifest)
    if m['stage']=='followup':
        gate=collect(verify(m['gate_manifest']))
        if not gate['analytic_support_gate'] or gate['status']!='complete':raise InvalidArtifact('analytic support gate has not passed')
    return execute_manifest(manifest)


def analyze(gate_collection,followup_collection,output):
    saved=[read_json(p) for p in (gate_collection,followup_collection)]
    for c in saved:
        replay=collect(verify(c['manifest']))
        reference=copy.deepcopy(c)
        # BLAS/vector arithmetic differs by <1e-12 across the actual compute/login
        # hosts. Permit 1e-10 only for derived mode projections; raw artifacts,
        # Cartesian gradients, energies and all other fields remain exact.
        for a,b in zip(replay['rows'],reference['rows']):
            if 'projected_gradient_kcal_per_unit' in a:
                ga=a.pop('projected_gradient_kcal_per_unit');gb=b.pop('projected_gradient_kcal_per_unit')
                if not np.allclose(ga,gb,atol=1e-10,rtol=0):raise InvalidArtifact('projected gradient replay differs')
        if replay!=reference:raise InvalidArtifact('actual collection replay differs')
    gate,follow=saved;gm=read_json(verify(gate['manifest']));fm=read_json(verify(follow['manifest']))
    if gm['design']!=fm['design']:raise InvalidArtifact('different physical experiments')
    d=read_json(verify(gm['design']));allrows=gate['rows']+follow['rows'];centers={}
    for r in allrows:
        if r['kind']=='center':centers[(r['case_id'],r['metal'],r['medium'])]=r
    cases=[]
    for case in CASES:
        sources={m:d['sources'][case+'__'+m] for m in ('Ca','La')};ep={};missing=[]
        for metal in ('Ca','La'):
            a,b=[centers[(case,metal,medium)] for medium in ('alpb','vacuum')]
            if a['status']!='complete' or b['status']!='complete':missing.append(metal);continue
            gv=np.asarray(sources[metal]['MACE_projected_gradient_kcal_per_unit'])
            gs=np.asarray(a['projected_gradient_kcal_per_unit'])-b['projected_gradient_kcal_per_unit']
            ep[metal]={'vacuum_MACE':gv.tolist(),'solvent_transfer':gs.tolist(),'composite':(gv+gs).tolist()}
        row={'case_id':case,'status':'unavailable' if missing else 'complete','missing_metals':missing,
             'mode_ids':sources['Ca']['mode_ids'],'mode_units':sources['Ca']['mode_units'],'endpoints':ep,
             'derivative_scope':'shared physical donor coordinates; endpoint-specific fixed water H retained',
             'affinity_correction':None}
        if not missing:
            diff={term:(np.asarray(ep['Ca'][term])-ep['La'][term]).tolist() for term in ('vacuum_MACE','solvent_transfer','composite')}
            row['Ca_minus_La_projected_derivatives']=diff
            angular=np.array([u=='radian' for u in row['mode_units']]);a,b=[np.array(diff[k])[angular] for k in ('vacuum_MACE','composite')]
            flips=(a*b<0)&(np.abs(a)>.02)&(np.abs(b)>.02)
            row['angular_sign_changes_above_0p02_floor']=int(flips.sum());row['angular_modes']=int(angular.sum())
            row['angular_vector_cosine_vacuum_vs_composite']=float(a@b/(np.linalg.norm(a)*np.linalg.norm(b))) if np.linalg.norm(a)*np.linalg.norm(b)>0 else None
        cases.append(row)
    checks=[]
    for case,mode in MODES.items():
        for metal in ('Ca','La'):
            source=d['sources'][case+'__'+metal];idx=source['mode_ids'].index(mode)
            for medium in ('vacuum','alpb','transfer'):
                media=('alpb','vacuum') if medium=='transfer' else (medium,)
                relevant=[r for r in allrows if r['case_id']==case and r['metal']==metal and r['medium'] in media]
                check={'case_id':case,'metal':metal,'mode_id':mode,'medium':medium,'status':'unavailable'}
                if any(r['status']!='complete' for r in relevant):checks.append(check);continue
                def value(step):
                    vals=[next(r for r in relevant if r['kind']=='displacement' and r['step_radian']==step and r['medium']==s)['energy_hartree'] for s in media]
                    return vals[0]-(vals[1] if len(vals)==2 else 0)
                gradients=[centers[(case,metal,s)]['projected_gradient_kcal_per_unit'][idx] for s in media]
                exact=gradients[0]-(gradients[1] if len(gradients)==2 else 0)
                numerical={str(h):(value(h)-value(-h))/(2*h)*HA_TO_KCAL for h in (.005,.01)}
                tolerance=max(.02,.005*abs(exact));error=abs(numerical['0.005']-exact);refinement=abs(numerical['0.005']-numerical['0.01'])
                check.update(status='pass' if error<=tolerance and refinement<=tolerance else 'fail',
                    analytic_kcal_mol_rad=exact,central_energy_derivative_kcal_mol_rad=numerical,
                    analytic_error_kcal_mol_rad=error,refinement_difference_kcal_mol_rad=refinement,tolerance_kcal_mol_rad=tolerance)
                checks.append(check)
    result={'protocol_id':PROTOCOL,'design':gm['design'],'collections':[record(gate_collection),record(followup_collection)],
        'analyzer_implementation':record(__file__),'cases':cases,'derivative_checks':checks,
        'derivative_check_status':'pass' if all(r['status']=='pass' for r in checks) else 'failed_or_unavailable',
        'center_energy_checks':gate['center_energy_checks']+follow['center_energy_checks'],
        'complete_tasks':sum(r['status']=='complete' for r in allrows),'task_denominator':len(allrows),
        'affinity_correction':None,'entropy':None,'new_fits':0,'new_DFT_MACE_calls':0,'baseline_changed':False}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='operation',required=True)
    q=s.add_parser('prepare')
    for k in ('source','agreement','output'):q.add_argument('--'+k,required=True)
    for op in ('validate','execute','collect'):
        q=s.add_parser(op);q.add_argument('--manifest',required=True)
        if op=='collect':q.add_argument('--output',required=True)
    q=s.add_parser('analyze')
    for k in ('gate_collection','followup_collection','output'):q.add_argument('--'+k.replace('_','-'),required=True)
    a=vars(p.parse_args());op=a.pop('operation')
    if op=='collect':dest=a.pop('output');r=collect(**a);write_new(dest,r)
    else:r=globals()[op](**a)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','cases','derivative_checks')},indent=2))


if __name__=='__main__':main()
