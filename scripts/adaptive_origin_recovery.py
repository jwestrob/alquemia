"""Two explicit origin-only adaptive pools; reuse the established physical engine."""
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
from adaptive_completion_folds import source_pair
from adaptive_force_diagnostic import project
from accommodation_nonlinear import WarmGPU
from mace_site_kinematics import Kinematics
from mace_hybrid import check_atoms,write_xyz
import nikasha_pool as pool

PROTOCOL='Nikasha_direct_origin_minimal_adaptive_recovery_v1'
CASES=('p12293-pqq-la_model__conditioned_Ca__seed-1_sample-4',
       'q60ar6-pqq-la_model__conditioned_Ca__seed-1_sample-1')
CANDIDATES=('origin','adaptive_Ca','adaptive_La')


def snapshot(source,dest):
    dest.mkdir();pins={}
    for p in Path(source).glob('*.py'):
        q=dest/p.name;shutil.copyfile(p,q);pins[p.name]=record(q)
    return pins


def check_reference(reference,parent):
    r=read_json(reference)
    if (r['protocol_id']!='Nikasha_origin_and_two_adaptive_candidates_replay_v1' or
        r['candidate_ids']!=list(CANDIDATES) or r['noncanonical_folds_used_for_calibration'] or
        r['signature']['proposal_protocol_id']!=completion.PROTOCOL or
        r['signature']['proposal_settings']!=completion.SETTINGS):raise InvalidArtifact('minimal reference method differs')
    for k in ('model','software','orca'):
        if r['signature'][k]!=parent[k]:raise InvalidArtifact('reference electronic method differs')
    verify(r['collection']);verify(r['agreement']);verify(r['implementation'])
    return r


def origin_cell(t,parent):
    # Reuse the existing exact source/recipe/receipt audit, without any old-candidate gate.
    return pool.reused_cell(t,{'origin_xyz':t['xyz'],'origin_components':t['q0']['components']},
                            {'origin':t['origin_reuse']['source_point']},'origin',parent)


def prepare(source,diagnostic,reference,ledger,agreement,output):
    parent=read_json(source);diag=read_json(diagnostic);check_reference(reference,parent)
    prior=read_json(ledger)
    if len(prior['cases'])!=225:raise InvalidArtifact('225-source ledger required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];cases=[]
    for cid in CASES:
        old=next(c for c in prior['cases'] if c['case_id']==cid)
        c={'case_id':cid,'status':'unavailable','reason':None,'old_row':old,'origins':{}}
        try:
            pair=source_pair(Path(source).resolve(),parent,cid,diag)
            for t in pair:
                cell=origin_cell(t,parent);c['origins'][t['metal']]=cell
            tasks.extend(pair);c.update(status='prepared')
        except (InvalidArtifact,OSError,ValueError,KeyError,StopIteration) as exc:c['reason']=str(exc)
        cases.append(c)
    m={k:parent[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources')}
    m.update(protocol_id=PROTOCOL,settings=completion.SETTINGS,proposal_protocol_id=completion.PROTOCOL,
             source_manifest=record(source),diagnostic=record(diagnostic),reference=record(reference),ledger=record(ledger),
             agreement=record(agreement),declared_case_ids=list(CASES),cases=cases,tasks=tasks,
             implementation=snapshot(Path(__file__).parent,out/'implementation'),
             optimizer_software={'version':scipy.__version__,'wrapper':record(scipy_slsqp.__file__),'kernel':record(scipy_kernel.__file__)},
             optimizer_start_limit=4,new_cross_MACE_limit=4,new_GFN2_limit=16,new_DFT_calls=0,production_changed=False)
    path=out/'manifest.json';write_new(path,m);v=validate(path);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest):
    m=read_json(manifest);parent=read_json(verify(m['source_manifest']));diag=read_json(verify(m['diagnostic']))
    check_reference(verify(m['reference']),parent);verify(m['agreement']);verify(m['ledger'])
    if (m['protocol_id']!=PROTOCOL or m['settings']!=completion.SETTINGS or m['declared_case_ids']!=list(CASES)
        or [c['case_id'] for c in m['cases']]!=list(CASES) or m['optimizer_start_limit']!=4
        or m['new_cross_MACE_limit']!=4 or m['new_GFN2_limit']!=16):raise InvalidArtifact('frozen recovery scope changed')
    for pin in m['implementation'].values():verify(pin)
    for k in ('software','orca','cpu_executable','gpu_executable'):verify(m[k])
    for k in ('model','software','orca','cpu_python','gpu_python'):
        if m[k]!=parent[k]:raise InvalidArtifact('source method changed')
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer software differs')
    tasks={t['task_id']:t for t in m['tasks']};expected=[]
    for c in m['cases']:
        if c['status']!='prepared':continue
        pair=source_pair(verify(m['source_manifest']),parent,c['case_id'],diag)
        for t in pair:
            expected.append(t['task_id'])
            if tasks.get(t['task_id'])!=t or c['origins'][t['metal']]!=origin_cell(t,parent):raise InvalidArtifact('exact source/selector/origin changed')
            kin=Kinematics(read_json(verify(t['mapping']))['context'])
            angular.final_geometry(kin,t,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))])
        paired(verify(pair[1]['xyz']),verify(pair[0]['xyz']),pair[1]['charge'],pair[0]['charge'])
    if [t['task_id'] for t in m['tasks']]!=expected:raise InvalidArtifact('undeclared or missing task')
    return {'status':'validated','manifest':record(manifest),'cases':2,'prepared':sum(c['status']=='prepared' for c in m['cases']),
            'optimizer_starts':len(m['tasks']),'new_calls_in_validation':0}


def allocation():
    if (not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or
        int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=200000):raise InvalidArtifact('32CPU/200000MiB/oneGPU allocation required')


def execute(manifest):
    validate(manifest);allocation();m=read_json(manifest);root=Path(manifest).parent
    with (root/'proposal.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);start=time.monotonic();gpu=None;error=None;results=[]
        try:
            if m['tasks']:gpu=WarmGPU(manifest)
            for t in m['tasks']:
                p=completion.optimize(manifest,t,gpu);results.append(p)
                print(json.dumps({'task_id':t['task_id'],'status':read_json(verify(p))['status']}),flush=True)
        except Exception as exc:error=str(exc)
        finally:
            try:
                if gpu:gpu.close()
            except Exception as exc:error=(error+'; ' if error else '')+str(exc)
            elapsed=time.monotonic()-start
            write_new(root/('EXECUTION_'+os.environ['SLURM_JOB_ID']+'.json'),{'manifest':record(manifest),'results':results,
                'error':error,'wall_seconds':elapsed,'allocated_core_seconds':32*elapsed,'GPU_seconds':elapsed,
                'job_id':os.environ['SLURM_JOB_ID']})
        if error:raise InvalidArtifact(error)
    return {'completed_search_records':len(results)}


def prepare_pool(proposals,output):
    validate(proposals);m=read_json(proposals);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    bytask={(t['case_id'],t['metal']):t for t in m['tasks']};cases=[];tasks=[]
    for source in m['cases']:
        cid=source['case_id'];c={'case_id':cid,'status':'unavailable','reason':source['reason'],
            'old_result':source['old_row']['old_result'],'prior_pool':source['old_row']['pool'],
            'candidates':[{'id':'origin','xyz':bytask[cid,'Ca']['xyz']}] if source['status']=='prepared' else [],
            'matrix':{z:{'origin':source['origins'][z]} if z in source['origins'] else {} for z in ('Ca','La')}}
        pending=[]
        try:
            if source['status']!='prepared':raise InvalidArtifact('source origin unsupported')
            points={}
            for z in ('Ca','La'):
                t=bytask[cid,z];rp=Path(proposals).parent/'proposals'/t['task_id']/'result.json';r=read_json(rp)
                if r['manifest']!=record(proposals) or r['status']!='proposal_available':raise InvalidArtifact('required adaptive candidate unavailable: '+z)
                point=r['proposal'];kin=Kinematics(read_json(verify(t['mapping']))['context']);q=np.asarray(point['full_q'])
                if any(q[i] for i in set(range(len(q)))-set(t['active_indices'])):raise InvalidArtifact('inactive mode moved')
                atoms=xyz(verify(point['coordinate']));angular.final_geometry(kin,t,q[t['active_indices']],[a[0] for a in atoms])
                if not np.allclose(kin.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):raise InvalidArtifact('candidate mapping mismatch')
                points[z]=point;c['candidates'].append({'id':'adaptive_'+z,'xyz':point['coordinate'],'proposal_receipt':record(rp)})
            for z in ('Ca','La'):
                for maker in ('Ca','La'):
                    t=bytask[cid,z];name='adaptive_'+maker;point=points[maker];a=xyz(verify(point['coordinate']));atoms=[(z,*a[0][1:]),*a[1:]]
                    check_atoms(atoms,t['charge']);tid=cid+'__'+z+'__at_'+name;d=out/'cells'/tid;d.mkdir(parents=True)
                    xp=d/'context.xyz';write_xyz(xp,atoms)
                    if xyz(xp)!=atoms:raise InvalidArtifact('coordinate serialization changed')
                    nt={'task_id':tid,'case_id':cid,'metal':z,'candidate':name,'xyz':record(xp),
                        'charge':t['charge'],'multiplicity':t['multiplicity'],'source_mapping':t['mapping'],'source_preparation':t['source_preparation']}
                    if z==maker:nt['native_reuse']=point['MACE'];pool.native_reuse(nt,m)
                    pending.append(nt);c['matrix'][z][name]={'status':'pending','task_id':tid,'xyz':record(xp),'reused':False}
            c.update(status='prepared',reason=None);tasks.extend(pending)
        except (InvalidArtifact,OSError,KeyError,ValueError) as exc:c['reason']=str(exc)
        cases.append(c)
    pm={k:m[k] for k in ('model','software','orca','cpu_python','gpu_python','resources','agreement','reference','source_manifest','ledger')}
    pm.update(protocol_id=PROTOCOL,proposal_manifest=record(proposals),proposal_settings=completion.SETTINGS,
        settings={**pool.SETTINGS,'candidate_order':list(CANDIDATES)},tasks=tasks,cases=cases,
        implementation=snapshot(verify(m['implementation']['adaptive_origin_recovery.py']).parent,out/'implementation'),
        shard_count=1,GFN2_maxiter=500,numerical_policy_id='native_GFN2_MaxIter500_unchanged_convergence_v1',
        new_MACE_cells=sum(not t.get('native_reuse') for t in tasks),maximum_new_GFN2_calls=2*len(tasks))
    pm['numerical_qualification']=read_json(verify(m['reference']))['numerical_provenance'][0]['qualification']
    qualification=read_json(verify(pm['numerical_qualification']))
    if not qualification['both_controls_pass'] or not qualification['formerly_failed_cell_complete']:
        raise InvalidArtifact('MaxIter500 qualification absent')
    if pm['new_MACE_cells']>4 or pm['maximum_new_GFN2_calls']>16:raise InvalidArtifact('finite recovery calls exceeded')
    mp=out/'manifest.json';write_new(mp,pm);pool.low_prepare(mp);v=validate_pool(mp);write_new(out/'PREFLIGHT.json',v);return v


def validate_pool(manifest):
    m=read_json(manifest);parent=read_json(verify(m['proposal_manifest']));validate(verify(m['proposal_manifest']))
    if (m['protocol_id']!=PROTOCOL or [c['case_id'] for c in m['cases']]!=list(CASES) or
        m['settings']!={**pool.SETTINGS,'candidate_order':list(CANDIDATES)} or m['GFN2_maxiter']!=500):raise InvalidArtifact('pool policy changed')
    for k in ('model','software','orca','reference','agreement'):
        if m[k]!=parent[k]:raise InvalidArtifact('pool source method differs')
    for pin in m['implementation'].values():verify(pin)
    qualification=read_json(verify(m['numerical_qualification']))
    if not qualification['both_controls_pass'] or not qualification['formerly_failed_cell_complete']:
        raise InvalidArtifact('MaxIter500 qualification absent')
    tasks={t['task_id']:t for t in m['tasks']};expected=[]
    for c in m['cases']:
        origin=next(v for v in parent['cases'] if v['case_id']==c['case_id'])
        if c['status']!='prepared':continue
        if [v['id'] for v in c['candidates']]!=list(CANDIDATES):raise InvalidArtifact('candidate set changed')
        for candidate in c['candidates'][1:]:
            maker=candidate['id'].split('_')[-1]
            ot=next(t for t in parent['tasks'] if t['case_id']==c['case_id'] and t['metal']==maker)
            rp=Path(m['proposal_manifest']['path']).parent/'proposals'/ot['task_id']/'result.json'
            if record(rp)!=candidate['proposal_receipt']:raise InvalidArtifact('candidate receipt differs')
            r=read_json(rp)
            if r['manifest']!=m['proposal_manifest'] or r['status']!='proposal_available' or r['proposal']['coordinate']!=candidate['xyz']:
                raise InvalidArtifact('candidate not produced by declared physical search')
            kin=Kinematics(read_json(verify(ot['mapping']))['context']);q=np.asarray(r['proposal']['full_q'])
            if any(q[i] for i in set(range(len(q)))-set(ot['active_indices'])):raise InvalidArtifact('inactive coordinate moved')
            angular.final_geometry(kin,ot,q[ot['active_indices']],[a[0] for a in xyz(verify(candidate['xyz']))])
        for z in ('Ca','La'):
            if c['matrix'][z]['origin']!=origin['origins'][z]:raise InvalidArtifact('origin component changes')
            ot=next(t for t in parent['tasks'] if t['case_id']==c['case_id'] and t['metal']==z)
            for candidate in c['candidates'][1:]:
                t=tasks[c['matrix'][z][candidate['id']]['task_id']];expected.append(t['task_id']);atoms=xyz(verify(candidate['xyz']))
                if xyz(verify(t['xyz']))!=[(z,*atoms[0][1:]),*atoms[1:]] or (t['charge'],t['multiplicity'])!=(ot['charge'],ot['multiplicity']):raise InvalidArtifact('cross-cell geometry/state differs')
                if t.get('native_reuse'):pool.native_reuse(t,m)
    if len(expected)!=len(tasks) or set(expected)!=set(tasks) or len(tasks)>8:raise InvalidArtifact('finite cross-cell population differs')
    return {'status':'validated','manifest':record(manifest),'prepared':sum(c['status']=='prepared' for c in m['cases']),
            'new_MACE_cells':m['new_MACE_cells'],'new_GFN2_calls':m['maximum_new_GFN2_calls']}


def execute_cross(manifest):
    validate_pool(manifest);allocation();m=read_json(manifest);root=Path(manifest).parent
    with (root/'cross.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);start=time.monotonic();gpu=None;error=None;receipts=[]
        try:
            for t in m['tasks']:
                if t.get('native_reuse'):pool.native_reuse(t,m);continue
                if gpu is None:gpu=WarmGPU(manifest)
                d=root/'cells'/t['task_id'];request=d/'request.json'
                write_new(request,{k:t[k] for k in ('task_id','xyz','charge','multiplicity')}|{'manifest':record(manifest)})
                _,pin=gpu.evaluate(record(request));receipts.append(pin)
        except Exception as exc:error=str(exc)
        finally:
            try:
                if gpu:gpu.close()
            except Exception as exc:error=(error+'; ' if error else '')+str(exc)
            elapsed=time.monotonic()-start
            write_new(root/('CROSS_EXECUTION_'+os.environ['SLURM_JOB_ID']+'.json'),{'manifest':record(manifest),'receipts':receipts,
                 'error':error,'wall_seconds':elapsed,'allocated_core_seconds':32*elapsed,'GPU_seconds':elapsed,'job_id':os.environ['SLURM_JOB_ID']})
        if error:raise InvalidArtifact(error)
    return {'new_cross_calls':len(receipts)}


def collect(manifest,output):
    validate_pool(manifest)
    return pool.collect(manifest,output)


def score(collection,output):
    from nikasha_pool_compare import method_record
    result=read_json(collection);mp=verify(result['manifest']);validate_pool(mp);m=read_json(mp)
    ref=check_reference(verify(m['reference']),m);rows=[]
    if [c['case_id'] for c in result['cases']]!=list(CASES):raise InvalidArtifact('result source population differs')
    for c in result['cases']:
        p=c['pool'];expected=c['old_result']['expected_class'];values={}
        if c['status']=='prepared' and pool.choose_rows(c['matrix'],CANDIDATES)!=p:raise InvalidArtifact('pool algebra differs')
        for mode in ('mathematical','operational'):
            value=p[mode]['composite_R_model_kcal_mol'] if p['status']=='available' else None
            values[mode]=method_record(value,ref['variants'][mode]['bands'],expected)
        rows.append({'case_id':c['case_id'],'expected_class':expected,'pool_status':p['status'],
                     'reason':c.get('reason') or p.get('reason'),'decisions':values,'prior_pool':c['prior_pool'],
                     'old_static_R':c['old_result']['R0'],'selected_pool':p})
    out={'protocol_id':PROTOCOL,'collection':record(collection),'reference':m['reference'],'implementation':record(__file__),
         'rows':rows,'denominator':2,'available':sum(r['pool_status']=='available' for r in rows),
         'full225_overlay_owned_by_parent':True,'new_reference_fitted':False,'production_changed':False}
    write_new(output,out);return {k:v for k,v in out.items() if k!='rows'}


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='operation',required=True)
    q=sub.add_parser('prepare')
    for k in ('source','diagnostic','reference','ledger','agreement','output'):q.add_argument('--'+k,required=True,type=Path)
    for name in ('validate','execute','validate-pool','execute-cross','collect'):
        q=sub.add_parser(name);q.add_argument('--manifest',required=True,type=Path)
        if name=='collect':q.add_argument('--output',required=True,type=Path)
    q=sub.add_parser('prepare-pool');q.add_argument('--proposals',required=True,type=Path);q.add_argument('--output',required=True,type=Path)
    q=sub.add_parser('score');q.add_argument('--collection',required=True,type=Path);q.add_argument('--output',required=True,type=Path)
    args=vars(p.parse_args());op=args.pop('operation').replace('-','_');r=globals()[op](**args)
    print(json.dumps(r,indent=2))


if __name__=='__main__':main()
