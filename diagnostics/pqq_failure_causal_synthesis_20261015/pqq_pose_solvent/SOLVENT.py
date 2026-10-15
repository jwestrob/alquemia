#!/usr/bin/env python3
"""Finite pose continuation: reuse qualified SCC/native ALPB and full-pool policy."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/pqq_pose_solvent'
TEMPLATE=ROOT/'workspaces/exaa_metal_midpoint_20261014/scoring_v1/CONFIG.json'
BASES=('origin','adaptive_Ca','adaptive_La');METALS=('Ca','La');POLICIES=('operational','mathematical','origin')
VARIANTS=tuple(f'{q}__PQQ_{sign}1deg' for q in BASES for sign in ('minus','plus'));ORDER=BASES+VARIANTS
PROTOCOL='ExaA_intact_PQQ_e3_one_degree_frozen_ALPB_common9_v1'
KEYS=('matrix','runtime','settings','factorial','component','native_manifest','native_build','native_binary','native_adapter','MACE_model','corrected_ORCA_parameter_provenance','ORCA_correction_policy')

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])!=p:raise ValueError('Changed artifact: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def require(value,reason):
    if not value:raise ValueError(reason)

def prepare(score_config,score_collection,output):
    """Actual terminal artifact binding only; scientific gate replay is Slurm-only."""
    sc=read(score_config);sr=read(score_collection);old=read(TEMPLATE)
    require(sr['config']==pin(score_config) and sr['terminal'] and sr['complete_new_cells']==24,'Exact complete24-cell MACE collection required')
    require(sr['stage_gate']['decision']=='qualifies_for_root_solvent_review' and sr['stage_gate']['threshold_kcal_mol']==1.0,'Recorded declared solvent continuation gate not passed')
    require(len(sc['cells'])==24 and sc['new_order']==list(ORDER),'Declared common9 scope differs')
    execution=Path(score_collection).with_name('EXECUTION.json');ex=read(execution)
    require(ex['status']=='complete' and ex['config']==pin(score_config) and ex['actual_model_forwards']==24,'MACE execution is not complete')
    cells=[]
    for t in sc['cells']:
        source,=[r for r in sr['rows'] if r['sample']==t['sample']];v=source['matrix'][t['metal']][t['candidate']]
        require(v['status']=='complete' and v['xyz']==t['xyz'],'Completed MACE geometry missing')
        verify(v['receipt']);verify(t['xyz'])
        cells.append({**t,'MACE':v['receipt'],'MACE_eV':v['energy_eV']})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'SOLVENT.py');shutil.copyfile(HERE/'RUN.sbatch',out/'RUN.sbatch')
    c={k:old[k] for k in KEYS}
    c.update(protocol_id=PROTOCOL,implementation=pin(out/'SOLVENT.py'),wrapper=pin(out/'RUN.sbatch'),
       plan=pin(HERE.parent/'PQQ_POSE_PLAN.md'),adapter_plan=pin(HERE/'PLAN.md'),template=pin(TEMPLATE),
       score_config=pin(score_config),score_collection=pin(score_collection),score_execution=pin(execution),score_adapter=sc['implementation'],
       preparation=sc['preparation'],preparation_manifest=sc['preparation_manifest'],inventory=sc['inventory'],rows=sc['rows'],cells=cells,
       solvent_adapter=old['solvent_adapter'],pool_adapter=old['pool_adapter'],scorer=old['scorer'],
       baseline_config=old['previous_pool_config'],baseline_collection=old['previous_pool_collection'],
       numerical_python=sc['python'],CPUs=7,declared_sources=2,declared_new_cells=24,reused_original_cells=12,
       original_order=list(BASES),new_order=list(ORDER),maximum_SCC_calls=96,maximum_native_component_calls=48,
       new_MACE_calls=0,new_DFT_calls=0,new_optimizations=0,classification=None,
       parameter_role='Published standalone GFN2 parameter file and qualified internal native ALPB; corrected full ORCA JSON is audited provenance only, never default ORCA GFN2')
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def preflight(c,f,r,n,a,h,pose):
    for k in ('plan','adapter_plan','wrapper','corrected_ORCA_parameter_provenance','ORCA_correction_policy'):verify(c[k])
    require(pin(sys.executable)==c['numerical_python'],'Qualified numerical interpreter differs')
    old=read(c['template'])
    for k in KEYS:require(c[k]==old[k],'Qualified solvent kernel/settings changed: '+k)
    require(read(c['ORCA_correction_policy'])['candidate_sha256']==c['corrected_ORCA_parameter_provenance']['sha256'],'Corrected full La parameter lineage differs')
    require(c['protocol_id']==PROTOCOL and c['declared_new_cells']==24 and c['reused_original_cells']==12,'Finite scope differs')
    require(c['original_order']==list(BASES) and c['new_order']==list(ORDER),'Pool order differs')
    sc=read(c['score_config']);sr=read(c['score_collection']);ex=read(c['score_execution']);prep=read(c['preparation'])
    require(sr['config']==c['score_config'] and sr['terminal'] and sr['complete_new_cells']==24 and ex['status']=='complete' and ex['actual_model_forwards']==24,'Actual complete MACE matrix required')
    require(ex['config']==c['score_config'] and sr['preparation']==c['preparation']==sc['preparation'] and prep['rows']==c['rows'],'Actual preparation/MACE joins differ')
    require(sc['implementation']==c['score_adapter'] and sc['plan']==c['plan'] and sc['inventory']==c['inventory'],'MACE protocol provenance differs')
    worker=read(sc['worker_manifest']);b=load('_pqq_pose_saved_worker_metadata',worker['implementation']);s=load('_pqq_pose_original_runtime_checks',sc['previous_score_adapter'])
    pf=pose.preflight(sc,b,s);mr={};by_sample={row['sample']:row for row in sr['rows']};saved_results={}
    require(len(c['cells'])==24 and {row['sample'] for row in c['rows']}=={0,3},'Finite24/two-source scope differs')
    for cell in c['cells']:
        t,=[t for t in sc['cells'] if t['task_id']==cell['task_id']]
        require(all(cell[k]==t[k] for k in t),'Completed MACE input changed')
        v=by_sample[t['sample']]['matrix'][t['metal']][t['candidate']]
        require(v['status']=='complete' and v['receipt']==cell['MACE'] and v['energy_eV']==cell['MACE_eV'],'Completed MACE receipt changed')
        saved_results[t['task_id']]=v;mr[t['task_id']]=a.mace(cell,c)
        symbols,_=f.coordinates(cell['xyz'])
        require(len(symbols)==192 and sum(f.VALENCE[z] for z in symbols)-cell['charge']==cell['valence_electrons'],'Native valence inventory differs')
    replay=pose.collect(sc,pf,saved_results)
    require(replay['rows']==sr['rows'] and replay['stage_gate']==sr['stage_gate'],'Saved MACE full-pool/gate replay differs')
    require(all(row['common3_to9']['mathematical'] is not None for row in replay['rows']),'Full MACE common9 unavailable')
    require(any(abs(row['common3_to9']['mathematical']['paired_contrast_change_kcal_mol'])>=1.0 for row in replay['rows']),'Declared sign-blind1kcal stage gate failed')
    runtime=read(c['runtime']);nm=read(c['native_manifest']);build=read(c['native_build'])
    for k in ('xtb','xtb_parameter_file','tblite_python'):verify(runtime[k])
    for p in runtime['xtb_dynamic']['libraries']+runtime['tblite_distribution']['distribution_files']:verify(p)
    for d in runtime['tblite_dynamic'].values():
        for p in d['libraries']:verify(p)
    for p in list(nm['implementation'].values())+nm['build']['pins']+build['linked_libraries']:verify(p)
    verify(build['binary']);require(build['manifest']==c['native_manifest'] and build['binary']==c['native_binary'] and nm['MACE_model']==c['MACE_model'],'Qualified native build/model differs')
    baseline=read(c['baseline_collection']);bc=read(c['baseline_config']);inv=read(c['inventory']);old_rows={};count=0
    require(baseline['config']==c['baseline_config'] and baseline['available_cells']==20,'Actual completed historical baseline differs')
    for source in c['rows']:
        stored,=[row for row in baseline['rows'] if row['source_case_id']==source['source_case_id']]
        bound,=[row for row in bc['rows'] if row['source_case_id']==source['source_case_id']]
        src,=[row for row in inv['rows'] if row['case_id']==source['source_case_id']]
        evidence={z:{} for z in METALS}
        for cell in (cell for cell in bound['cells'] if cell['candidate'] in BASES):
            z,q=cell['metal'],cell['candidate'];original=src['matrix'][z][q]
            require(cell['original']==original,'Exact original inventory cell differs')
            for k in ('maps','state_signature','source_preparation'):require(cell[k]==source[k],'Original source state/map differs')
            a.mace(original,c);v=a.reusable(c,cell,n,r);saved=stored['evidence'][z][q]
            for k in ('MACE','MACE_eV','xyz','native_result','G_qvac_hartree','components_hartree'):require(v[k]==saved[k],'Saved baseline evidence changed: '+k)
            require(v['MACE']==pf['baseline'][str(source['sample'])][z][q]['receipt'],'MACE and solvent baseline sources differ')
            evidence[z][q]=v;count+=1
        require(all(set(evidence[z])==set(BASES) for z in METALS),'Incomplete original common3')
        old_rows[str(source['sample'])]={'evidence':evidence,'historical3':stored['historical3']}
    require(count==12,'Exactly twelve existing solvent cells must be reused')
    return {'status':'passed','stage_gate':replay['stage_gate'],'MACE_collection':c['score_collection'],'reused_MACE':mr,
            'old_rows':old_rows,'reused_original_cells':12,'new_molecular_calls':0,
            'qualification':'unchanged component-v2, native water diagonal, cavity equality and vacuum charge-precision gates; absolute v1 results retained'}

def summarize(c,pf,sr,a,h,n,r,selector):
    native={v['task_id']:v for v in sr};rows=[];residuals=[];available=0
    for source in c['rows']:
        old=pf['old_rows'][str(source['sample'])];evidence=copy.deepcopy(old['evidence']);fail=[]
        for cell in (cell for cell in c['cells'] if cell['sample']==source['sample']):
            z,q=cell['metal'],cell['candidate'];v=native.get(cell['task_id']);evidence[z][q]={'status':'unavailable','reason':'not executed'}
            try:
                require(v is not None and v['status']=='available','Native/component endpoint unavailable: '+str(v.get('reason') if v else 'missing'))
                original={**cell,'native_frozen_ALPB_result':v['result'],'native_frozen_ALPB_hartree':v['G_qvac_hartree']}
                evidence[z][q]=a.reusable(c,{'original':original,**{k:cell[k] for k in ('state_signature','maps','source_preparation')}},n,r);available+=1
            except Exception as e:
                evidence[z][q]={'status':'unavailable','reason':f'{type(e).__name__}: {e}'};fail.append({'task_id':cell['task_id'],'metal':z,'candidate':q,'native_record':v,'reason':str(e)})
        matrix={z:{q:({'MACE_eV':v['MACE_eV'],'solvent_hartree':v['G_qvac_hartree'],'xyz':v['xyz'],'MACE':v['MACE'],'native_result':v['native_result']} if v['status']=='available' else None) for q,v in evidence[z].items()} for z in METALS}
        pools={name:h.choose(selector,matrix,order) for name,order in (('original3',BASES),('expanded9',ORDER))}
        h.reproduce(pools['original3'],old['historical3'],residuals)
        works={z:{q:a.work(v,evidence[z]['origin']) for q,v in evidence[z].items()} for z in METALS};changes={}
        for policy in POLICIES:
            before,after=pools['original3'][policy],pools['expanded9'][policy]
            if before is None or after is None:changes[policy]=None;continue
            metals={z:a.difference(works[z][after[z]['candidate']],works[z][before[z]['candidate']]) for z in METALS}
            delta=a.difference(metals['Ca'],metals['La']);raw=after['R_model_kcal_mol']-before['R_model_kcal_mol'];h.close(raw,delta['total'],residuals)
            changes[policy]={'per_metal_change_kcal_mol':metals,'differential_change_kcal_mol':delta,'R_model_change_kcal_mol':raw,
                'before':{z:before[z]['candidate'] for z in METALS},'after':{z:after[z]['candidate'] for z in METALS}}
        finite={}
        for q in VARIANTS:
            base=q.split('__PQQ_',1)[0];metal={z:a.difference(works[z][q],works[z][base]) for z in METALS}
            finite[q]={'base':base,'per_metal_pose_work_kcal_mol':metal,'Ca_minus_La_pose_work_kcal_mol':a.difference(metal['Ca'],metal['La'])}
        rows.append({'sample':source['sample'],'role':source['role'],'source_case_id':source['source_case_id'],'pools':pools,
            'common3_to9':changes,'all_origin_relative_works_kcal_mol':works,'finite_pose_responses':finite,
            'evidence':evidence,'failed_new_cells':fail,'classification':None})
    by={row['sample']:row for row in rows};contrasts={}
    for policy in POLICIES:
        target,control=by[3]['common3_to9'][policy],by[0]['common3_to9'][policy]
        contrasts[policy]=a.difference(target['differential_change_kcal_mol'],control['differential_change_kcal_mol']) if target and control else None
    return {'protocol_id':PROTOCOL,'terminal':True,'declared_sources':2,'declared_new_cells':24,'available_new_cells':available,'reused_original_cells':12,
        'rows':rows,'target_minus_control_change_kcal_mol':contrasts,'stage_gate':pf['stage_gate'],
        'pool_orders':{'original3':list(BASES),'expanded9':list(ORDER)},
        'saved_algebra_checks':{'count':len(residuals),'maximum_residual_kcal_mol':max(residuals,default=None),'tolerance_kcal_mol':h.TOL},
        'new_MACE_calls':0,'new_DFT_calls':0,'new_optimizations':0,'classification':None,'new_calibration':False,'production_changed':False,
        'interpretation':'Finite fixed-pose common-pool response; no affinity, physical minimum or predictive-validation claim'}

def execute(path,action):
    require(bool(os.environ.get('SLURM_JOB_ID')),'Actual Slurm required before scientific imports')
    c=read(path);require(c['implementation']==pin(__file__),'Wrong adapter snapshot')
    matrix=load('_pqq_pose_unchanged_solvent_dispatch',c['matrix']);placement=matrix.placement()
    out=Path(path).parent/('execution' if action=='execute' else 'preflight_'+os.environ['SLURM_JOB_ID']);start=time.monotonic();pf=None;sr=[]
    receipt={'status':'failed','config':pin(path),'placement':placement,'action':action}
    with (Path(path).parent/'SOLVENT.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(exist_ok=False)
        try:
            sys.path.insert(0,str(ROOT/'scripts'))
            f=load('_pqq_pose_factorial',c['factorial']);r=load('_pqq_pose_component',c['component']);n=load('_pqq_pose_native',c['native_adapter'])
            a=load('_pqq_pose_prior_solvent_checks',c['solvent_adapter']);h=load('_pqq_pose_pool_checks',c['pool_adapter']);selector=load('_pqq_pose_selection',c['scorer']);pose=load('_pqq_pose_vacuum_gate',c['score_adapter'])
            pf=preflight(c,f,r,n,a,h,pose);write(out/'PREFLIGHT.json',pf)
            if action=='execute':sr=matrix.solvent(c,out,f,r,n,pf['reused_MACE'])
            receipt['status']='complete'
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
        finally:
            try:
                if pf and action=='execute':
                    if not sr:
                        for cell in c['cells']:
                            p=out/'solvent/replay/tasks'/cell['task_id']/'RESULT.json'
                            if p.exists():
                                v=read(p);sr.append({k:v[k] for k in ('task_id','status','reason','G_qvac_hartree')}|{'result':pin(p)})
                    result={'config':pin(path),**summarize(c,pf,sr,a,h,n,r,selector)};write(out/'FINAL_COLLECTION.json',result)
                    if result['available_new_cells']!=24:receipt['status']='partial' if receipt['status']=='complete' else receipt['status']
            except BaseException as e:receipt.update(status='failed',collection_failure={'reason':str(e),'traceback':traceback.format_exc()});raise
            finally:
                elapsed=time.monotonic()-start;receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,
                    SCC_started=sum(1 for _ in out.glob('solvent/points/*/*/START.json')),SCC_terminal=sum(1 for _ in out.glob('solvent/points/*/*/EXECUTION.json')),
                    native_component_calls=sum(1 for _ in out.glob('solvent/replay/tasks/*/*/native_EXECUTION.json')),
                    measured_resource_receipts={'SCC':[pin(p) for p in sorted(out.glob('solvent/points/*/*/EXECUTION.json'))],
                      'native_ALPB':[pin(p) for p in sorted(out.glob('solvent/replay/tasks/*/*/native_EXECUTION.json'))]},
                    new_MACE_calls=0,new_DFT_calls=0,new_optimizations=0)
                write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('prepare');q.add_argument('--score-config',required=True);q.add_argument('--score-collection',required=True);q.add_argument('--output',required=True)
    for name in ('preflight','execute'):
        q=sub.add_parser(name);q.add_argument('--config',required=True);q.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.score_config,a.score_collection,a.output)
    else:execute(verify({'path':str(Path(a.config).resolve()),'sha256':a.config_sha}),a.action)
