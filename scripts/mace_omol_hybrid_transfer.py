"""Source-backed structural transfer of the unchanged matched vacuum response."""
from __future__ import annotations
import argparse,copy,json,shutil,time
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,cache_key,paired,read_json,record,verify,write_new,xyz
from affordable_response import source_key
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL,accepted_attempt,check_atoms,write_xyz
from mace_omol_gradient_run import CONFIG,check_parent,model
from mace_omol_hybrid_response import task as response_task, POLICY
from mace_omol_ablation import ADAPTER as MASK,COMPONENT
from mace_omol_ablation_run import SEMANTICS
from mace_omol_vacuum import METHOD,scientific_input,parse_endpoint
from mace_metal_response import grid,matrix
from mace_bounded_response import ball,endpoint_checks
STAGE='matched_hybrid_GGR_transfer'
STRUCTURES=('GGR_2FW0','GGR_2FVY')
REPRESENTATIONS={'extended':'alpha_caps','connected':'connected_segment'}
COUNTS={'extended':58,'connected':111}


def normalized_core(repair,whole):
    """Use actual source identities and cap anchors; do not re-protonate."""
    from mace_mechanics import physical_model,original_full,physical_id
    from mace_omol_matched_h import normalized_physical
    graph,coords,pos,mk,_,_=physical_model(repair)
    original=original_full(whole,graph,coords,pos,mk)
    replay=normalized_physical(whole,original)
    byid={a['id']:a for a in whole['physical_atoms']};ends={};mapping=[]
    if whole['source']!=repair['source_structure'] or whole['explicit_waters'] or repair['explicit_water_inventory']:
        raise InvalidArtifact('transfer must retain identical zero-water source')
    for metal in ('Ca','La'):
        e=repair['outputs'][metal];old=xyz(verify(e['xyz']));new=list(old);seen={0}
        for item in repair['atom_graph']['source_to_qm']:
            i=item['qm_index']
            if i in seen:raise InvalidArtifact('duplicate core source index')
            seen.add(i)
            if item['kind']=='source':
                aid=physical_id(item['source']);target=byid[aid];src=coords[source_key(item['source'])]
                if old[i][0]!=target['element'] or not np.allclose(old[i][1:],src,atol=1e-9,rtol=0):raise InvalidArtifact('core source identity differs')
                if old[i][0]=='H':new[i]=('H',*target['xyz_A'])
                elif not np.allclose(old[i][1:],target['xyz_A'],atol=1e-9,rtol=0):raise InvalidArtifact('physical heavy atom changed')
                if metal=='Ca':mapping.append({**item,'physical_id':aid})
            elif item['kind']=='sigma_link_H':
                x,y=(np.array(byid[physical_id(item[k])]['xyz_A']) for k in ('retained','omitted'))
                expected=x+item['length_A']*(y-x)/np.linalg.norm(y-x)
                if old[i][0]!='H' or not np.allclose(old[i][1:],expected,atol=1e-9,rtol=0):raise InvalidArtifact('cap does not match fixed source anchors')
                if metal=='Ca':mapping.append(copy.deepcopy(item))
            else:raise InvalidArtifact('unsupported source mapping')
        if seen!=set(range(len(old))) or e['charge']!=(-1 if metal=='Ca' else 0):raise InvalidArtifact('core mapping/charge does not close')
        if any(old[i]!=new[i] for i in range(len(old)) if old[i][0]!='H'):raise InvalidArtifact('heavy coordinate changed')
        ends[metal]={'atoms':new,'charge':e['charge'],'spin_multiplicity':1,'state':check_atoms(new,e['charge']),'metal_index':0}
    if ends['Ca']['atoms'][1:]!=ends['La']['atoms'][1:] or ends['Ca']['atoms'][0][1:]!=ends['La']['atoms'][0][1:]:raise InvalidArtifact('paired coordinates differ')
    return ends,mapping,replay


@cached_file_checks
def prepare(config,output):
    from ggr_preparation import prepare as expand
    from mace_omol_prepared import audit_preparation
    start=time.monotonic();cpu=time.process_time();cfg=read_json(config)
    if set(cfg['sources'])!=set(STRUCTURES):raise InvalidArtifact('fixed two-structure inventory differs')
    for key in ('agreement','inventory','software','core_report','reference_report','quantum_parent'):verify(cfg[key])
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);cases={};fulls={}
    for gid,refs in cfg['sources'].items():
        wp=verify(refs['whole']);w=read_json(wp);audit=audit_preparation(wp);old=read_json(verify(refs['scalar_report']))
        mp=verify(old['manifest']);pm=read_json(mp)
        if w['case_id']!=gid or old['preparation']!=refs['whole'] or not old['numerical_gate_pass']:raise InvalidArtifact('whole archive source differs')
        fulls[gid]={}
        for metal,e in w['endpoints'].items():
            t=next(t for t in pm['tasks'] if t['metal']==metal and t['position']=='bound');r=old['rows'][t['task_id']]
            if xyz(verify(t['xyz']))!=xyz(verify(e['xyz'])) or t['charge']!=e['charge'] or not any(accepted_attempt(a,t,mp)==r for a in (mp.parent/'execution'/t['task_id']).glob('attempt_*')):
                raise InvalidArtifact('whole scalar reuse lacks matching physical input/receipt')
            fulls[gid][metal]={**e,'metal_index':next(i for i,a in enumerate(w['physical_atoms']) if a['id']=='metal'),
                'preparation':refs['whole'],'assembly':w['assembly'],'microstate':w['microstate'],'explicit_waters':w['explicit_waters'],'evidence':w['evidence'],
                'archived_center_energy_eV':r['energy_eV'],'scalar_archive':refs['scalar_report']}
        row=read_json(verify(refs['source_row']))
        if row['preparation_manifest']!=w['source_preparation']:raise InvalidArtifact('source row and whole core parent differ')
        for representation,policy in REPRESENTATIONS.items():
            name=gid+'_'+representation;d=out/name;d.mkdir();expand(verify(row['preparation_manifest']),d/'source_core',policy)
            rp=d/'source_core/preparation_manifest.json'
            repair=read_json(rp);ends,mapping,replay=normalized_core(repair,w);stored={}
            for metal,e in ends.items():
                if len(e['atoms'])!=COUNTS[representation]:raise InvalidArtifact('transfer core composition differs')
                xp=d/(metal+'.xyz');write_xyz(xp,e['atoms']);stored[metal]={k:v for k,v in e.items() if k!='atoms'};stored[metal]['xyz']=record(xp)
            value={'case_id':name,'global_id':gid,'representation':representation,'source_preparation':record(rp),'whole_preparation':refs['whole'],
                   'endpoints':stored,'mapping':mapping,'normalization_replay':replay,'whole_audit':audit,'evidence':w['evidence'],'explicit_waters':w['explicit_waters']}
            write_new(d/'mapping.json',value);cases[name]=record(d/'mapping.json')
    result={'protocol_id':POLICY,'transfer_stage':STAGE,'config':record(config),'cases':cases,'fulls':fulls,
            'agreement':cfg['agreement'],'reference_report':cfg['reference_report'],'inventory':cfg['inventory'],'software':cfg['software'],
            'core_report':cfg['core_report'],'quantum_parent':cfg['quantum_parent'],'preparation_wall_seconds':time.monotonic()-start,'preparation_CPU_seconds':time.process_time()-cpu}
    write_new(out/'preparation.json',result);return {'status':'prepared','preparation':record(out/'preparation.json'),'cases':list(cases)}


def prepared_states(preparation):
    p=read_json(preparation);cfg=read_json(verify(p['config']));cases={name:read_json(verify(ref)) for name,ref in p['cases'].items()}
    if p['protocol_id']!=POLICY or p['transfer_stage']!=STAGE or set(cases)!={g+'_'+r for g in STRUCTURES for r in REPRESENTATIONS}:raise InvalidArtifact('transfer preparation inventory differs')
    for k in ('agreement','reference_report','inventory','software','core_report','quantum_parent'):
        if p[k]!=cfg[k]:raise InvalidArtifact('preparation source/config differs')
        verify(p[k])
    for name,c in cases.items():
        if c['case_id']!=name or c['global_id'] not in STRUCTURES or c['whole_preparation']!=cfg['sources'][c['global_id']]['whole']:raise InvalidArtifact('case source/identity differs')
        verify(c['source_preparation']);w=read_json(verify(c['whole_preparation']));repair=read_json(verify(c['source_preparation']))
        byid={a['id']:a for a in w['physical_atoms']}
        for metal,e in c['endpoints'].items():
            rows=xyz(verify(e['xyz']));original=xyz(verify(repair['outputs'][metal]['xyz']));want=list(original)
            for item in c['mapping']:
                i=item['qm_index']
                if item['kind']=='source' and original[i][0]=='H':want[i]=('H',*byid[item['physical_id']]['xyz_A'])
            if [x[0] for x in rows]!=[x[0] for x in want] or not np.allclose([x[1:] for x in rows],[x[1:] for x in want],atol=5e-10,rtol=0):raise InvalidArtifact('normalized coordinates changed')
            if e['charge']!=repair['outputs'][metal]['charge'] or e['state']!=check_atoms(rows,e['charge']):raise InvalidArtifact('state charge/parity differs')
            full=p['fulls'][c['global_id']][metal]
            if any(full[k]!=w['endpoints'][metal][k] for k in ('xyz','charge','spin_multiplicity','state')):raise InvalidArtifact('whole state differs')
        paired(verify(c['endpoints']['La']['xyz']),verify(c['endpoints']['Ca']['xyz']),0,-1)
    return p,cases


def quantum_manifest(p,preparation,states,out,phase,native_preparation=None):
    from mace_global_benchmark import snapshot
    from ggr_sensitivity import ORCA
    out.mkdir();parent=read_json(verify(p['quantum_parent']));pins=snapshot(out,(*parent['implementation'],'mace_omol_hybrid_transfer.py','mace_omol_hybrid_transfer_minimum.py'))
    tasks=[]
    for key,s in states.items():
        name,metal=s['case_id'],s['metal'];e=s['core'];d=out/key;d.mkdir();xp=d/'core.xyz';shutil.copyfile(verify(e['xyz']),xp);ip=d/'endpoint.inp';ip.write_text(scientific_input(e['charge']))
        t={'task_id':key,'case_id':name,'metal':metal,'charge':e['charge'],'multiplicity':1,'input':record(ip),'xyz':record(xp),
           'source_mapping':p['cases'][name],'output_path':str(d/'endpoint.out'),'engrad_path':str(d/'endpoint.engrad'),'task_type':'analytic_gradient'}
        t['cache_key']=cache_key({'task':t,'method':METHOD,'protocol':POLICY,'implementation':pins});tasks.append(t)
    tasks.sort(key=lambda t:(-len(xyz(verify(t['xyz']))),t['task_id']))
    m={'protocol_id':POLICY,'stage':STAGE,'phase':phase,'preparation':record(preparation),'agreement':p['agreement'],'implementation':pins,
       'tasks':tasks,'orca':record(ORCA),'method':METHOD,'energy_scope':'isolated_vacuum_endpoint',
       'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},'cost_tracking':{'compute_budget':None,'wall_time_limit':None}}
    if native_preparation is not None:m['native_preparation']=native_preparation
    write_new(out/'manifest.json',m);return record(out/'manifest.json')


@cached_file_checks
def initial(preparation,output):
    from mace_omol import common,seal
    p,cases=prepared_states(preparation);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);jobs=[]
    centers=[(e,name,metal,'core','center',[0.,0.,0.]) for name,c in cases.items() for metal,e in c['endpoints'].items()]
    centers += [(e,gid,metal,'full','center',[0.,0.,0.]) for gid,ends in p['fulls'].items() for metal,e in ends.items()]
    groups={'centers':centers}
    for gid,ends in p['fulls'].items():groups['grid_'+gid]=[(e,gid,metal,'full',point,u) for metal,e in ends.items() for point,u in grid().items() if point!='center']
    for group,items in groups.items():
        _,d,m=common(verify(p['inventory']),verify(p['software']),verify(p['agreement']),out/group,STAGE)
        m.update(protocol_id=POLICY,phase='initial',group=group,preparation=record(preparation),core_report=p['core_report'],gradient_settings=CONFIG,model=model(verify(p['software'])),tasks=[response_task(e,name,metal,kind,point,u,d) for e,name,metal,kind,point,u in items])
        check_parent(m);seal(d,m);jobs.append(record(d/'manifest.json'))
    states={name+'_'+metal:{'case_id':name,'metal':metal,'core':e} for name,c in cases.items() for metal,e in c['endpoints'].items()}
    q=quantum_manifest(p,preparation,states,out/'quantum','initial')
    result={'preparation':record(preparation),'MACE_manifests':jobs,'quantum_manifest':q,'new_MACE_calls':156,'new_DFT_calls':8}
    write_new(out/'initial.json',result);validate_quantum(verify(q));return result


def expected_tasks(m):
    if m['phase']=='native':
        from mace_omol_hybrid_transfer_minimum import expected_native
        return expected_native(m)
    p,cases=prepared_states(verify(m['preparation']))
    if m['phase']!='initial':raise InvalidArtifact('unsupported transfer phase')
    items=[]
    if m['group']=='centers':
        items=[(e,name,metal,'core','center',[0.,0.,0.]) for name,c in cases.items() for metal,e in c['endpoints'].items()]
        items += [(e,gid,metal,'full','center',[0.,0.,0.]) for gid,ends in p['fulls'].items() for metal,e in ends.items()]
    elif m['group'].startswith('grid_') and m['group'][5:] in STRUCTURES:
        gid=m['group'][5:];items=[(e,gid,metal,'full',point,u) for metal,e in p['fulls'][gid].items() for point,u in grid().items() if point!='center']
    else:raise InvalidArtifact('undeclared transfer group')
    return p,{kind+'__'+name+'_'+metal+'__'+point:(e,name,metal,kind,point,u) for e,name,metal,kind,point,u in items}


@cached_file_checks
def validate(manifest):
    m=read_json(manifest);p,expected=expected_tasks(m)
    if m['protocol_id']!=POLICY or m['stage']!=STAGE or m['software']!=p['software'] or m['model']!=model(verify(p['software'])) or m['gradient_settings']!=CONFIG or m['core_report']!=p['core_report']:raise InvalidArtifact('transfer method differs')
    for pin in [m['agreement'],*m['implementation'].values()]:verify(pin)
    check_parent(m)
    if len(m['tasks'])!=len(expected) or {t['task_id'] for t in m['tasks']}!=set(expected):raise InvalidArtifact('finite task inventory differs')
    for t in m['tasks']:
        e,name,metal,kind,point,u=expected[t['task_id']];source=xyz(verify(e['xyz']));wanted=np.array([a[1:] for a in source]);wanted[e['metal_index']]+=u;actual=xyz(verify(t['xyz']))
        fields={**e,'case_id':name,'metal':metal,'kind':kind,'point':point,'source_xyz':e['xyz'],'displacement_A':u,'energy_only':m['phase']=='initial' and point!='center','derivative_backend':'checkpointed','descriptor_gradient_experiment':POLICY,'energy_component':COMPONENT,'charge_feature_adapter':MASK,'output_semantics':SEMANTICS,'capture_native_readout':True}
        if any(t.get(k)!=v for k,v in fields.items() if k!='xyz'):raise InvalidArtifact('scientific task setting differs')
        if [a[0] for a in actual]!=[a[0] for a in source] or not np.allclose([a[1:] for a in actual],wanted,atol=5e-10,rtol=0) or check_atoms(actual,t['charge'])!=t['state']:raise InvalidArtifact('paired/source coordinate or charge differs')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):raise InvalidArtifact('transfer cache differs')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


@cached_file_checks
def validate_quantum(manifest):
    if read_json(manifest)['phase']=='native':
        from mace_omol_hybrid_transfer_minimum import validate_quantum as validate_native
        return validate_native(manifest)
    from affordable_workflow import dry_run
    m=read_json(manifest);p,cases=prepared_states(verify(m['preparation']))
    if m['phase']!='initial' or m['stage']!=STAGE or m['protocol_id']!=POLICY or m['method']!=METHOD:raise InvalidArtifact('transfer quantum method/phase differs')
    expected={name+'_'+metal:(name,metal,e) for name,c in cases.items() for metal,e in c['endpoints'].items()}
    if len(m['tasks'])!=len(expected) or {t['task_id'] for t in m['tasks']}!=set(expected):raise InvalidArtifact('quantum inventory differs')
    for pin in [m['agreement'],m['orca'],*m['implementation'].values()]:verify(pin)
    for t in m['tasks']:
        name,metal,e=expected[t['task_id']];payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['source_mapping']!=p['cases'][name] or t['xyz']['sha256']!=e['xyz']['sha256'] or t['charge']!=e['charge'] or t['multiplicity']!=1 or verify(t['input']).read_text()!=scientific_input(e['charge']):raise InvalidArtifact('quantum source coordinates/charge/recipe differs')
        if t['cache_key']!=cache_key({'task':payload,'method':METHOD,'protocol':POLICY,'implementation':m['implementation']}):raise InvalidArtifact('quantum cache differs')
    return dry_run(manifest)


@cached_file_checks
def collect_quantum(manifest):
    from ggr_sensitivity import executed
    validate_quantum(manifest);m,rows=executed(manifest)
    for t in m['tasks']:
        r=rows[t['task_id']]
        if r['status']!='complete':continue
        try:r.update(parse_endpoint(t,verify(r['output']),t['engrad_path']))
        except (ValueError,OSError) as exc:r.update(status='invalid',reason=str(exc),energy_hartree=None)
    return {'status':'complete' if all(r['status']=='complete' for r in rows.values()) else 'incomplete','manifest':record(manifest),'rows':rows}


@cached_file_checks
def collect(manifest):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[];checks=[]
    for t in m['tasks']:
        good=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp);attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None})
            if r is not None:good.append(r)
        r=good[-1] if good else {'status':'unavailable','energy_eV':None};rows[t['task_id']]=r
        if r['status']=='computed':
            error=r['native_readout']['component_sum_error_kcal_mol'];checks.append({'task_id':t['task_id'],'name':'readout_closure','error_kcal':error,'pass':abs(error)<=.01})
            if t['point']=='center' and 'archived_center_energy_eV' in t:
                error=(r['energy_eV']-t['archived_center_energy_eV'])*EV_TO_KCAL;checks.append({'task_id':t['task_id'],'name':'archived_scalar_agreement','error_kcal':error,'pass':abs(error)<=.01})
    return {'status':'complete' if all(r['status']=='computed' for r in rows.values()) else 'incomplete','manifest':record(mp),'rows':rows,'attempts':attempts,'checks':checks,'numerical_checks_pass':bool(checks) and all(c['pass'] for c in checks)}


@cached_file_checks
def assess(initial,output):
    start=time.monotonic();cpu=time.process_time();job=read_json(initial);p,cases=prepared_states(verify(job['preparation']))
    qc=collect_quantum(verify(job['quantum_manifest']));rows={};collections=[]
    if qc['status']!='complete':raise InvalidArtifact('initial quantum endpoints incomplete')
    for ref in job['MACE_manifests']:
        c=collect(verify(ref))
        if c['status']!='complete' or not c['numerical_checks_pass']:raise InvalidArtifact('initial learned endpoints incomplete or fail checks')
        if set(rows)&set(c['rows']):raise InvalidArtifact('duplicate learned task')
        rows.update(c['rows']);collections.append(c)
    if len(rows)!=156:raise InvalidArtifact('initial collection inventory differs')
    matrices={};odd={};result_rows={};scores={}
    for gid,ends in p['fulls'].items():
        for metal,e in ends.items():
            key=gid+'_'+metal;center=rows['full__'+key+'__center'];values={point:(rows['full__'+key+'__'+point]['energy_eV']-center['energy_eV'])*EV_TO_KCAL for point in grid()}
            matrices[key]={scale:matrix(values,scale) for scale in ('coarse','fine')};odd[key]=[]
            g=np.load(verify(center['gradient']))[e['metal_index']]*EV_TO_KCAL
            for scale,h in [('coarse',.02),('fine',.01)]:
                for i,axis in enumerate('xyz'):
                    observed=(values[scale+'_'+axis+'p']-values[scale+'_'+axis+'m'])/2;pred=h*g[i];tol=max(.01,.01*abs(pred))
                    odd[key].append({'scale':scale,'axis':axis,'error_kcal':float(observed-pred),'tolerance_kcal':float(tol),'pass':bool(abs(observed-pred)<=tol)})
    for name,c in cases.items():
        energies={};result_rows[name]={}
        for metal,e in c['endpoints'].items():
            key=name+'_'+metal;gid=c['global_id'];d=qc['rows'][key];core=rows['core__'+key+'__center'];full=rows['full__'+gid+'_'+metal+'__center'];fe=p['fulls'][gid][metal]
            gd=np.array(d['gradient_kcal_mol_per_A'][0]);gj=(np.load(verify(full['gradient']))[fe['metal_index']]-np.load(verify(core['gradient']))[0])*EV_TO_KCAL;g=gd+gj;k=matrices[gid+'_'+metal]
            fine=ball(k['fine'],g);coarse=ball(k['coarse'],g);pred=copy.deepcopy(fine)
            if fine['status']=='solved' and coarse['status']=='solved':
                u=np.array(fine['displacement_A']);err=float(.5*u@(k['fine']-k['coarse'])@u);delta=float(np.linalg.norm(u-coarse['displacement_A']))
                pred.update(refinement_error_kcal_mol=err,refinement_displacement_A=delta,status='eligible_for_native_check' if abs(err)<=.02 and delta<=.01 and all(x['pass'] for x in odd[gid+'_'+metal]) else 'numerical_response_not_qualified')
            energies[metal]={'DFT_hartree':d['energy_hartree'],'full_eV':full['energy_eV'],'core_eV':core['energy_eV']}
            result_rows[name][metal]={'prediction':pred,'DFT_gradient_kcal_mol_A':gd.tolist(),'J_gradient_kcal_mol_A':gj.tolist(),'gradient_kcal_mol_A':g.tolist(),'matrices_kcal_mol_A2':{key:value.tolist() for key,value in k.items()}}
        ca,la=energies['Ca'],energies['La'];dr=(ca['DFT_hartree']-la['DFT_hartree'])*HA_TO_KCAL;fr=(ca['full_eV']-la['full_eV'])*EV_TO_KCAL;cr=(ca['core_eV']-la['core_eV'])*EV_TO_KCAL
        scores[name]={'static_R_kcal_scale':dr+fr-cr,'DFT_R_kcal_mol':dr,'full_R_model_kcal':fr,'core_R_model_kcal':cr,'endpoints':energies,'evidence':c['evidence']}
    reference=read_json(verify(p['reference_report']));contrasts=[]
    for alpha in ('ALPHA_1F6S','ALPHA_6IP9'):
        for name,s in scores.items():
            margin=reference['scores'][alpha]['static_R_kcal_scale']-s['static_R_kcal_scale'];contrasts.append({'alpha':alpha,'GGR':name,'static_margin_kcal':margin,'static_direction_pass':margin>.02})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);sp=out/'assessment_source.py';sp.write_bytes(Path(__file__).read_bytes())
    result={'status':'complete','protocol_id':POLICY,'initial':record(initial),'preparation':job['preparation'],'quantum':qc,'MACE_collections':collections,'rows':result_rows,'static_scores':scores,'static_contrasts':contrasts,'odd_checks':odd,'implementation':record(sp),'baseline_changed':False,'predictions_are_not_validated_scores':True,'report_wall_seconds':time.monotonic()-start,'report_CPU_seconds':time.process_time()-cpu}
    write_new(out/'result.json',result);return {'status':'complete','static_contrasts':contrasts,'predictions':{n:{metal:x['prediction']['status'] for metal,x in row.items()} for n,row in result_rows.items()}}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    for op,keys in [('assess',('initial','output')),('prepare',('config','output')),('initial',('preparation','output')),('validate',('manifest',)),('validate_quantum',('manifest',)),('collect',('manifest',)),('collect_quantum',('manifest',))]:
        q=sub.add_parser(op)
        for k in keys:q.add_argument('--'+k,required=True)
    a=vars(p.parse_args());op=a.pop('command');print(json.dumps(globals()[op](**a),indent=2))
