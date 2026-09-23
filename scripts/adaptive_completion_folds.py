"""Primary225 source adapter for the unchanged scaled four-angle optimizer."""
from __future__ import annotations
import argparse
import copy
import fcntl
import json
import os
from pathlib import Path
import shutil
import time
import numpy as np
import scipy
import scipy.optimize._slsqp_py as scipy_slsqp
import scipy.optimize._slsqplib as scipy_kernel
from affordable_common import InvalidArtifact,paired,read_json,record,verify,write_new,xyz
import adaptive_completion as completion
import adaptive_angular_proposals as angular
from adaptive_force_diagnostic import endpoint,preview,project
from accommodation_nonlinear import WarmGPU
from mace_site_kinematics import Kinematics

PROTOCOL=completion.PROTOCOL
SCALED_POOL='nikasha_scaled_angular_common_geometry_native_OMOL_GFN2_ALPB_v1'


def source_pair(source,parent,cid,diagnostic):
    ts=[next(t for t in parent['tasks'] if t['case_id']==cid and t['metal']==z) for z in ('Ca','La')]
    pairs=[endpoint(source,parent,t,'origin') for t in ts]
    if any(data is None for _,data in pairs):raise InvalidArtifact('actual paired origin forces unavailable')
    ca,la=[data for _,data in pairs]
    if ca['mapping']!=la['mapping'] or not np.allclose(ca['coordinates'],la['coordinates'],atol=1e-12,rtol=0):
        raise InvalidArtifact('source pair geometry/map differs')
    choice=preview(ca['mapping']['modes'],ca['vectors'],ca['normalized'],la['normalized'])
    old=next(r for r in diagnostic['pairs'] if r['case_id']==cid)
    if choice!=old['selection_preview'] or len(choice['selected'])!=4:
        raise InvalidArtifact('original common four-mode selector does not replay')
    names=[v['id'] for v in ca['mapping']['modes']];ids=[v['id'] for v in choice['selected']]
    selected=[names.index(n) for n in ids];tasks=[]
    for t,(row,data) in zip(ts,pairs):
        if t['q0_status']!='available':raise InvalidArtifact('source q0 components unavailable')
        oldpoint=read_json(verify(row['proposal_receipt']))['origin']
        if any(oldpoint['full_q']):raise InvalidArtifact('source origin is not q0')
        f=np.load(verify(oldpoint['forces']),allow_pickle=False)
        _,_,raw,_,_=project(data['mapping'],oldpoint['full_q'],f)
        point={**oldpoint,'active_q_radian':[0.]*4,'gradient_kcal_mol_rad':raw[selected].tolist()}
        new=copy.deepcopy(t)
        new.update(active_indices=selected,active_mode_ids=ids,active_roles=['adaptive_physical_angular']*4,
                   original_active_mode_ids=t['active_mode_ids'],selector=choice,
                   prior_proposal_manifest=record(source),
                   origin_reuse={'proposal_receipt':row['proposal_receipt'],'source_point':oldpoint,'point':point,
                                 'projection':'actual_same_q0_cartesian_forces_through_selected_physical_J'})
        tasks.append(new)
    return tasks


def prepare(source,base_collection,diagnostic,agreement,output,execution_shards=1):
    if execution_shards not in (1,4):raise InvalidArtifact('only declared serial/four-shard execution supported')
    parent=read_json(source);base=read_json(base_collection);bm=read_json(verify(base['manifest']));diag=read_json(diagnostic)
    if bm['source_manifest']!=record(source) or len(base['cases'])!=225:
        raise InvalidArtifact('exact primary225 base/source union required')
    if record(source) not in [a['manifest'] for a in diag['archives']]:raise InvalidArtifact('force diagnostic source differs')
    sourcecases={c['case_id']:c for c in parent['cases']+parent['unavailable_cases']}
    if set(sourcecases)!={c['case_id'] for c in base['cases']}:raise InvalidArtifact('primary225 membership differs')
    tasks=[];cases=[]
    for b in base['cases']:
        cid=b['case_id'];c={'case_id':cid,'source_case':sourcecases[cid],
                          'inherited_base_pool_status':b['pool']['status'],'status':'unavailable',
                          'reason':b['pool'].get('reason',b.get('reason'))}
        if b['pool']['status']=='available':
            try:
                pair=source_pair(source,parent,cid,diag)
                tasks.extend(pair);c.update(status='prepared',reason=None)
            except (InvalidArtifact,KeyError,ValueError,OSError) as exc:c['reason']=str(exc)
        cases.append(c)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        dst=impl/p.name;shutil.copyfile(p,dst);pins[p.name]=record(dst)
    m={**{k:parent[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources')},
       'protocol_id':PROTOCOL,'settings':completion.SETTINGS,'stage':'primary225','source_manifest':record(source),
       'base_collection':record(base_collection),'diagnostic':record(diagnostic),'agreement':record(agreement),
       'declared_case_ids':[c['case_id'] for c in cases],'cases':cases,'tasks':tasks,'implementation':pins,
       'optimizer_software':{'version':scipy.__version__,'wrapper':record(scipy_slsqp.__file__),'kernel':record(scipy_kernel.__file__)},
       'new_optimizer_starts':len(tasks),'endpoint_denominator':450,'new_DFT_calls':0,'new_GFN2_calls_in_this_adapter':0,
       'source_states_changed':False,'reference':None,'production_changed':False,'all_evidence_consumed':True,
       'execution_requires_separate_frozen_reference':True,
       'execution_partition':{'rule':'task_index_modulo','shards':execution_shards}}
    write_new(out/'manifest.json',m);r=validate(out/'manifest.json');write_new(out/'PREFLIGHT.json',r);return r


def validate(manifest):
    m=read_json(manifest);parent=read_json(verify(m['source_manifest']));base=read_json(verify(m['base_collection']))
    bm=read_json(verify(base['manifest']));diag=read_json(verify(m['diagnostic']))
    if m['protocol_id']!=PROTOCOL or m['settings']!=completion.SETTINGS or m['stage']!='primary225':
        raise InvalidArtifact('scaled numerical protocol differs')
    partition=m.get('execution_partition',{'rule':'task_index_modulo','shards':1})
    if partition.get('rule')!='task_index_modulo' or partition.get('shards') not in (1,4):
        raise InvalidArtifact('unsupported execution partition')
    if bm['source_manifest']!=m['source_manifest'] or m['declared_case_ids']!=[b['case_id'] for b in base['cases']] or len(m['cases'])!=225:
        raise InvalidArtifact('source population/order changed')
    verify(m['agreement'])
    for k in ('software','orca','cpu_executable','gpu_executable'):verify(m[k])
    for k in ('model','software','orca','cpu_python','gpu_python'):
        if m[k]!=parent[k]:raise InvalidArtifact('scientific execution settings changed')
    for p in m['implementation'].values():verify(p)
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('SciPy version differs')
    tasks={t['task_id']:t for t in m['tasks']}
    if len(tasks)!=len(m['tasks']) or m['endpoint_denominator']!=450:raise InvalidArtifact('endpoint denominator differs')
    for c,b in zip(m['cases'],base['cases']):
        if c['case_id']!=b['case_id'] or c['inherited_base_pool_status']!=b['pool']['status']:
            raise InvalidArtifact('base availability changed')
        pair=[tasks.get(c['case_id']+'__'+z) for z in ('Ca','La')]
        if c['status']!='prepared':
            if any(t is not None for t in pair):raise InvalidArtifact('unavailable source gained executable task')
            continue
        if b['pool']['status']!='available' or any(t is None for t in pair):
            raise InvalidArtifact('prepared case lacks compatible base/paired tasks')
        for t in pair:
            old=next(o for o in parent['tasks'] if o['task_id']==t['task_id'])
            for k in ('xyz','mapping','source_preparation','charge','multiplicity','q0','mode_count'):
                if t[k]!=old[k]:raise InvalidArtifact('physical/electronic source changed: '+k)
            choice=next(r['selection_preview'] for r in diag['pairs'] if r['case_id']==c['case_id'])
            if t['selector']!=choice or t['active_mode_ids']!=[r['id'] for r in choice['selected']]:
                raise InvalidArtifact('common selector differs')
            origin=t['origin_reuse'];recorded=read_json(verify(origin['proposal_receipt']))
            if recorded['origin']!=origin['source_point'] or recorded['manifest']!=m['source_manifest']:
                raise InvalidArtifact('archived original force record differs')
            point=origin['point'];kin=Kinematics(read_json(verify(t['mapping']))['context'])
            forces=np.load(verify(point['forces']),allow_pickle=False)
            coords,_,raw,_,_=project(kin.data,point['full_q'],forces)
            expected={**origin['source_point'],'active_q_radian':[0.]*4,
                      'gradient_kcal_mol_rad':raw[t['active_indices']].tolist()}
            if point!=expected or any(point['full_q']):raise InvalidArtifact('actual-force reprojection differs')
            if [kin.modes[i]['id'] for i in t['active_indices']]!=t['active_mode_ids']:
                raise InvalidArtifact('selected physical map indices differ')
            if not np.allclose(coords,[a[1:] for a in xyz(verify(t['xyz']))],atol=1e-12,rtol=0):
                raise InvalidArtifact('source coordinate reconciliation exceeds existing1e-12 policy')
            native=read_json(verify(point['MACE']));request=read_json(verify(native['request']))
            if (native['model']!=m['model'] or native['forces']!=point['forces'] or
                    request['charge']!=t['charge'] or request['multiplicity']!=t['multiplicity'] or
                    xyz(verify(request['xyz']))!=xyz(verify(t['xyz']))):
                raise InvalidArtifact('native source state/model differs')
        if pair[0]['active_mode_ids']!=pair[1]['active_mode_ids']:raise InvalidArtifact('paired selected IDs differ')
        paired(verify(pair[1]['xyz']),verify(pair[0]['xyz']),pair[1]['charge'],pair[0]['charge'])
    return {'status':'prepared_dry_run_pass','manifest':record(manifest),'cases':225,
            'prepared_cases':sum(c['status']=='prepared' for c in m['cases']),
            'tasks':len(m['tasks']),'endpoint_denominator':450,'new_molecular_calls':0}


def shard_tasks(manifest_data,shard_index):
    partition=manifest_data.get('execution_partition',{'rule':'task_index_modulo','shards':1})
    count=partition['shards']
    if partition['rule']!='task_index_modulo' or count not in (1,4) or not 0<=shard_index<count:
        raise InvalidArtifact('invalid fixed execution shard')
    return [(i,t) for i,t in enumerate(manifest_data['tasks']) if i%count==shard_index]


def execute(manifest,reference,shard_index=0):
    validate(manifest);ref=read_json(reference)
    if (ref.get('protocol_id')!=SCALED_POOL or ref.get('proposal_protocol_id')!=PROTOCOL or
            ref.get('proposal_settings')!=completion.SETTINGS or not ref.get('frozen_at_UTC') or
            ref.get('noncanonical_folds_used_for_calibration') is not False):
        raise InvalidArtifact('matching separately frozen canonical reference required')
    if (not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or
            int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=200000):raise InvalidArtifact('existing32CPU/200000MiB/GPU allocation required')
    m=read_json(manifest);root=Path(manifest).parent;selected=shard_tasks(m,shard_index)
    partition=m.get('execution_partition',{'rule':'task_index_modulo','shards':1})
    lock_name='proposal.lock' if partition['shards']==1 else 'proposal_shard_'+str(shard_index)+'.lock'
    lock=(root/lock_name).open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);began=time.monotonic();gpu=None;error=None;results=[]
    try:
        gpu=WarmGPU(manifest)
        for _,t in selected:
            pin=completion.optimize(manifest,t,gpu);results.append(pin)
            print(json.dumps({'task_id':t['task_id'],'status':read_json(verify(pin))['status']}),flush=True)
    except Exception as exc:error=str(exc)
    finally:
        try:
            if gpu is not None:gpu.close()
        except Exception as exc:error=(error+'; ' if error else '')+str(exc)
        elapsed=time.monotonic()-began
        write_new(root/('execution_'+os.environ['SLURM_JOB_ID']+'_shard_'+str(shard_index)+'.json'),
                  {'manifest':record(manifest),'frozen_reference':record(reference),'results':results,'error':error,
                   'execution_partition':partition,'shard_index':shard_index,'lock_path':str(root/lock_name),
                   'selected_task_indices':[i for i,_ in selected],
                   'selected_task_ids':[t['task_id'] for _,t in selected],
                   'wall_seconds':elapsed,'allocated_core_seconds':elapsed*32,'allocated_GPU_seconds':elapsed,
                   'job_id':os.environ['SLURM_JOB_ID'],'new_DFT_calls':0,'new_GFN2_calls':0})
    if error:raise InvalidArtifact(error)
    return {'status':'proposal_shard_complete','shard_index':shard_index,
            'executed_tasks':len(selected),'endpoint_denominator':450}


def collect(manifest,output):
    validate(manifest);m=read_json(manifest);root=Path(manifest).parent;tasks={t['task_id']:t for t in m['tasks']};rows=[];cases=[]
    for c in m['cases']:
        pair=[]
        for z in ('Ca','La'):
            tid=c['case_id']+'__'+z;t=tasks.get(tid)
            row={'task_id':tid,'case_id':c['case_id'],'metal':z,'status':'unavailable',
                 'reason':c['reason'] if c['status']!='prepared' else 'not_run','candidate':None,
                 'proposal_receipt':None,'candidate_GFN2_vacuum_hartree':None,'candidate_GFN2_ALPB_hartree':None,
                 'composite_candidate_energy':None,'boundary_flag':None,'residual':None}
            if t is not None:
                row.update(source_preparation=t['source_preparation'],mapping=t['mapping'],charge=t['charge'],multiplicity=t['multiplicity'],
                           common_active_mode_ids=t['active_mode_ids'],origin_xyz=t['xyz'],origin_components=t['q0']['components'])
                rp=root/'proposals'/tid/'result.json';r=read_json(rp) if rp.exists() else None
                if r:
                    if r['manifest']!=record(manifest):raise InvalidArtifact('proposal manifest differs')
                    row.update(proposal_receipt=record(rp),reason=r.get('reason'),boundary_flag=r.get('boundary_flag'),residual=r.get('residual'))
                    if r['status']=='proposal_available':
                        kin=Kinematics(read_json(verify(t['mapping']))['context']);q=np.asarray(r['proposal']['full_q'])
                        atoms=xyz(verify(r['proposal']['coordinate']))
                        if any(q[i]!=0 for i in set(range(len(q)))-set(t['active_indices'])):raise InvalidArtifact('unselected mode moved')
                        angular.final_geometry(kin,t,q[t['active_indices']],[a[0] for a in atoms])
                        if not np.allclose(kin.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):raise InvalidArtifact('candidate map differs')
                        row.update(status='candidate_available',candidate=r['proposal'])
            pair.append(row);rows.append(row)
        cases.append({'case_id':c['case_id'],'status':'ready_for_common_pool' if len(pair)==2 and all(r['status']=='candidate_available' for r in pair) else 'unavailable',
                      'reason':c['reason'],'selected_geometry':None,'score':None,'source_case':c['source_case']})
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'source_manifest':m['source_manifest'],
            'cases':cases,'endpoints':rows,'case_denominator':225,'endpoint_denominator':450,
            'available_candidates':sum(r['status']=='candidate_available' for r in rows),'new_DFT_calls':0,
            'new_GFN2_calls_in_adapter':0,'common_pool_scoring_performed':False,'production_changed':False,'reference':None}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('cases','endpoints')}


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('prepare')
    for k in ('source','base_collection','diagnostic','agreement','output'):q.add_argument('--'+k.replace('_','-'),required=True)
    q.add_argument('--execution-shards',type=int,choices=(1,4),default=1)
    q=s.add_parser('dry-run');q.add_argument('--manifest',required=True)
    q=s.add_parser('execute');q.add_argument('--manifest',required=True);q.add_argument('--reference',required=True)
    q.add_argument('--shard-index',type=int,default=0)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('op');fn={'prepare':prepare,'dry-run':validate,'execute':execute,'collect':collect}[op]
    print(json.dumps(fn(**a),indent=2))


if __name__=='__main__':main()
