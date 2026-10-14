#!/usr/bin/env python3
"""One finite midpoint: existing warm MACE, component-v2/native ALPB, saved pools."""
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
import signal
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE=ROOT/'diagnostics/exaa_metal_midpoint_20261014/scoring'
PREVIOUS=ROOT/'workspaces/exaa_metal_position_20261014'
THREE=('origin','adaptive_Ca','adaptive_La')
FIVE=THREE+('transplanted_origin','transplanted_adaptive_La')
SEVEN=FIVE+('midpoint_adaptive_Ca','midpoint_adaptive_La')
DONORS=('adaptive_Ca','adaptive_La');METALS=('Ca','La');POLICIES=('operational','mathematical','origin')
PROTOCOL='ExaA_one_half_metal_midpoint_common7_v1'

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])!=p:raise ValueError('Changed pin: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,data):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def prepare(preparation,output):
    p=read(preparation);scorepath=PREVIOUS/'scoring_v1/CONFIG.json';score=read(scorepath)
    solventpath=PREVIOUS/'solvent_v1/CONFIG.json';solvent=read(solventpath)
    poolpath=PREVIOUS/'full_pool_check_v1/CONFIG.json';pool=read(poolpath)
    oldpath=poolpath.parent/'execution/FINAL_COLLECTION.json';old=read(oldpath)
    inventory=read(score['inventory']);rows=p['rows'];cells=[]
    if len(rows)!=2 or {r['sample'] for r in rows}!={0,3}:raise ValueError('Both declared source rows required')
    if old['config']!=pin(poolpath) or not old['terminal'] or old['available_cells']!=20:raise ValueError('Exact completed20-cell baseline required')
    for row in rows:
        source,=[r for r in inventory['rows'] if r['case_id']==row['source_case_id']]
        if row['status']!='prepared':continue
        for z in METALS:
            if set(row['matrix'][z])!=set(DONORS):raise ValueError('Wrong midpoint donor set')
            for donor in DONORS:
                v=row['matrix'][z][donor]
                cells.append({'task_id':f"sample{row['sample']}_{z}_midpoint_{donor}",'case_id':row['case_id'],
                  'sample':row['sample'],'metal':z,'candidate':'midpoint_'+donor,'donor_candidate':donor,
                  **{k:v[k] for k in ('xyz','charge','multiplicity','atom_count')},
                  'valence_electrons':source['matrix'][z][donor]['GFN2_valence_electrons'],
                  **{k:row[k] for k in ('maps','state_signature','source_preparation','geometry_admission')}})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'RUN.py');shutil.copyfile(HERE/'RUN.sbatch',out/'RUN.sbatch')
    worker=copy.deepcopy(read(score['worker_manifest']))
    worker.update(protocol_id=PROTOCOL,tasks=cells,origin_task_count=len(cells),cases=[],fixtures=[],preparation=pin(preparation),
      finite_scope={'declared_midpoint_cells':8,'admitted_cells':len(cells),'maximum_MACE_calls':8,'maximum_SCC_calls':32,'maximum_component_calls':16,'optimizer_starts':0},
      authorization='Root-authorized single-fraction midpoint diagnostic; actual Slurm only')
    write(out/'worker_manifest.json',worker)
    c={k:solvent[k] for k in ('matrix','runtime','settings','factorial','component','native_manifest','native_build','native_binary','native_adapter','MACE_model','corrected_ORCA_parameter_provenance','ORCA_correction_policy')}
    c.update(protocol_id=PROTOCOL,implementation=pin(out/'RUN.py'),wrapper=pin(out/'RUN.sbatch'),plan=pin(HERE.parent/'PLAN.md'),
      preparation=pin(preparation),rows=rows,cells=cells,worker_manifest=pin(out/'worker_manifest.json'),
      previous_scoring_config=pin(scorepath),previous_score_adapter=score['implementation'],
      previous_solvent_config=pin(solventpath),solvent_adapter=solvent['implementation'],
      previous_pool_config=pin(poolpath),pool_adapter=pool['implementation'],previous_pool_collection=pin(oldpath),
      scorer=pool['scorer'],inventory=score['inventory'],previous_preparation=score['preparation'],
      force_motivation=pin(PREVIOUS/'force_review_v1/execution/RESULT.json'),fraction=0.5,CPUs=7,
      declared_new_cells=8,reused_old_cells=20,maximum_MACE_calls=8,maximum_SCC_calls=32,maximum_native_component_calls=16,
      new_optimizations=0,new_DFT_calls=0,classification=None)
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def preflight(c,b,s,a,n,r,h):
    for k in ('plan','wrapper','force_motivation','corrected_ORCA_parameter_provenance','ORCA_correction_policy'):verify(c[k])
    oldsc=read(c['previous_scoring_config']);pf=s.preflight(oldsc,b);worker=read(c['worker_manifest']);oldworker=read(oldsc['worker_manifest'])
    for k in ('implementation','model','settings','pool_settings','source_modules','source_snapshot','calculator_source','model_python','model_executable','optimizer_python','optimizer_executable','gates'):
        if worker[k]!=oldworker[k]:raise ValueError('Warm model/runtime changed: '+k)
    if worker['tasks']!=c['cells'] or worker['fixtures'] or worker['cases'] or worker['threads_per_model']!=7 or worker['lanes']!=1:raise ValueError('Finite worker scope differs')
    if c['fraction']!=0.5 or c['declared_new_cells']!=8 or len(c['cells'])>8:raise ValueError('Declared midpoint scope differs')
    cfg=read(c['previous_solvent_config']);runtime=read(c['runtime']);nm=read(c['native_manifest']);build=read(c['native_build'])
    for key in ('runtime','settings','factorial','component','native_manifest','native_build','native_binary','native_adapter','MACE_model','corrected_ORCA_parameter_provenance','ORCA_correction_policy'):
        if c[key]!=cfg[key]:raise ValueError('Qualified solvent protocol changed: '+key)
    for k in ('xtb','xtb_parameter_file','tblite_python'):verify(runtime[k])
    for p in runtime['xtb_dynamic']['libraries']+runtime['tblite_distribution']['distribution_files']:verify(p)
    for d in runtime['tblite_dynamic'].values():
        for p in d['libraries']:verify(p)
    for p in list(nm['implementation'].values())+nm['build']['pins']+build['linked_libraries']:verify(p)
    verify(build['binary'])
    if build['manifest']!=c['native_manifest'] or build['binary']!=c['native_binary']:raise ValueError('Native build differs')
    inv=read(c['inventory']);prep=read(c['preparation']);prior=read(c['previous_preparation']);old=read(c['previous_pool_collection']);pc=read(c['previous_pool_config'])
    if prep['rows']!=c['rows'] or old['config']!=c['previous_pool_config'] or old['available_cells']!=20:raise ValueError('Actual preparation/baseline join differs')
    reused={};checks=[]
    for oldrow in old['rows']:
        bound,=[q for q in pc['rows'] if q['source_case_id']==oldrow['source_case_id']]
        for cell in bound['cells']:
            a.mace(cell['original'],cfg);v=a.reusable(cfg,cell,n,r);saved=oldrow['evidence'][cell['metal']][cell['candidate']]
            for k in ('MACE','MACE_eV','xyz','native_result','G_qvac_hartree','components_hartree'):
                if v[k]!=saved[k]:raise ValueError('Saved baseline evidence changed: '+k)
        reused[str(oldrow['sample'])]=oldrow
    for row in c['rows']:
        source,=[q for q in inv['rows'] if q['case_id']==row['source_case_id']]
        previous,=[q for q in prior['rows'] if q['source_case_id']==row['source_case_id']]
        if row['status']!='prepared':
            checks.append({'sample':row['sample'],'status':'unavailable','reason':row.get('reason')});continue
        admission=read(row['geometry_admission'])
        if not admission['pass'] or admission['status']!='prepared':raise ValueError('Midpoint geometry admission failed')
        for key in ('maps','state_signature','source_preparation'):
            if row[key]!=source[key]:raise ValueError('Physical state/map provenance differs')
        origin=b.xyz(verify(source['matrix']['Ca']['origin']['xyz']))[0][1:]
        transferred=b.xyz(verify(previous['matrix']['Ca']['origin']['xyz']))[0][1:]
        expected=tuple(x+0.5*(y-x) for x,y in zip(origin,transferred));paired={}
        for t in (t for t in c['cells'] if t['sample']==row['sample']):
            z,donor=t['metal'],t['donor_candidate'];v=row['matrix'][z][donor]
            if t['xyz']!=v['xyz'] or v['geometry_admission']!=row['geometry_admission'] or not v['geometry_checks']['pass']:raise ValueError('Actual cell geometry admission differs')
            new=verify(t['xyz']);oldxyz=verify(source['matrix'][z][donor]['xyz']);coords=b.xyz(new)
            if new.read_text().splitlines()[3:]!=oldxyz.read_text().splitlines()[3:]:raise ValueError('Nonmetal coordinate lines changed')
            if len(coords)!=192 or coords[0][0]!=z or tuple(coords[0][1:])!=expected:raise ValueError('Exact declared midpoint/composition differs')
            state=read(row['state_signature'])
            if t['charge']!=state['endpoint_charges'][z] or t['multiplicity']!=1 or state['water_inventory']:raise ValueError('Midpoint physical state differs')
            if not all(math.isfinite(x) for atom in coords for x in atom[1:]):raise ValueError('Nonfinite midpoint coordinates')
            paired[z,donor]=coords
        for donor in DONORS:
            ca,la=paired['Ca',donor],paired['La',donor]
            if ca[1:]!=la[1:] or ca[0][1:]!=la[0][1:]:raise ValueError('Paired midpoint geometry differs')
        checks.append({'sample':row['sample'],'status':'passed','geometry_admission':row['geometry_admission'],'exact_fraction':0.5,'unchanged_nonmetal_lines':True})
    return {'status':'passed','versions':pf['versions'],'checks':checks,'old_rows':reused,'reused_old_cells':20,'new_model_calls':0}

def run_mace(c,m,b,s,base,versions):
    out=base/'mace';(out/'lane0').mkdir(parents=True);transport=None;results={};handlers=b.handlers();error=None;start=time.monotonic()
    try:
        if c['cells']:
            transport=b.WarmCPU(verify(c['worker_manifest']),m,out,0)
            if transport.ready['threads']!=7 or transport.ready['versions']!=versions:raise ValueError('Loaded model runtime differs')
        for t in c['cells']:
            request={'manifest':c['worker_manifest'],'task_id':t['task_id'],**{k:t[k] for k in ('xyz','charge','multiplicity')}}
            result,rp=b.cell(out/'cells'/t['task_id'],request,transport,m['model'])
            row={'status':'unavailable','receipt':rp,'energy_eV':None,'reason':result.get('reason')}
            try:
                s.admitted_result(result,read(result['request']),t,m['model'],versions)
                row.update(status='complete',energy_eV=result['energy_eV'])
            except Exception as e:row['reason']=str(e)
            results[t['task_id']]=row
        if transport:transport.close();transport=None
    except BaseException as e:error=str(e);raise
    finally:
        if transport:transport.abort()
        for sig,handler in handlers.items():signal.signal(sig,handler)
        write(out/'PHASE.json',{'rows':results,'reason':error,'elapsed_seconds':time.monotonic()-start,'declared_cells':8,
           'dispatched_calls':sum(1 for _ in out.glob('cells/*/DISPATCH_STARTED.json')),
           'actual_model_forwards':sum(read(p).get('actual_model_forwards',0) for p in out.glob('cells/*/mace_result.json'))})
    return results

def summarize(c,pf,mr,sr,a,h,selector):
    native={r['task_id']:r for r in sr};rows=[];residuals=[]
    for source in c['rows']:
        old=pf['old_rows'][str(source['sample'])];evidence=copy.deepcopy(old['evidence']);fail=[]
        for z in METALS:
            for donor in DONORS:evidence[z]['midpoint_'+donor]={'status':'unavailable','reason':'midpoint preparation or evaluation unavailable'}
        for t in (t for t in c['cells'] if t['sample']==source['sample']):
            mm=mr.get(t['task_id']);ss=native.get(t['task_id']);z,q=t['metal'],t['candidate']
            if not mm or mm['status']!='complete' or not ss or ss['status']!='available':
                fail.append({'metal':z,'candidate':q,'MACE':mm,'solvent':ss});continue
            evidence[z][q]=a.energy_cell(c,{**t,'MACE':mm['receipt'],'MACE_eV':mm['energy_eV']},ss)
        matrix={z:{q:({'MACE_eV':v['MACE_eV'],'solvent_hartree':v['G_qvac_hartree'],'xyz':v['xyz'],'MACE':v['MACE'],'native_result':v['native_result']} if v['status']=='available' else None) for q,v in evidence[z].items()} for z in METALS}
        pools={name:h.choose(selector,matrix,order) for name,order in (('historical3',THREE),('previous5',FIVE),('new7',SEVEN))}
        h.reproduce(pools['historical3'],old['historical3'],residuals)
        for policy in POLICIES:
            now,prior=pools['previous5'][policy],old['extended5'][policy]
            h.close(now['R_model_kcal_mol'],prior['R_model_kcal_mol'],residuals)
            for z in METALS:
                if now[z]['candidate']!=prior[z]['candidate']:raise ValueError('Saved previous5 selection changed')
                for q in FIVE:h.close(now[z]['all_works_kcal_mol'][q],prior[z]['all_works_kcal_mol'][q],residuals)
        works={z:{q:a.work(v,evidence[z]['origin']) for q,v in evidence[z].items()} for z in METALS};changes={}
        for baseline in ('historical3','previous5'):
            changes[baseline]={}
            for policy in POLICIES:
                oldpool,newpool=pools[baseline][policy],pools['new7'][policy]
                if oldpool is None or newpool is None:changes[baseline][policy]=None;continue
                metals={z:a.difference(works[z][newpool[z]['candidate']],works[z][oldpool[z]['candidate']]) for z in METALS}
                diff=a.difference(metals['Ca'],metals['La']);delta=newpool['R_model_kcal_mol']-oldpool['R_model_kcal_mol'];h.close(delta,diff['total'],residuals)
                changes[baseline][policy]={'per_metal_change_kcal_mol':metals,'differential_change_kcal_mol':diff,'R_model_change_kcal_mol':delta,
                    'choices':{z:{'before':oldpool[z]['candidate'],'after':newpool[z]['candidate'],'before_xyz':oldpool[z]['components']['xyz'],'after_xyz':newpool[z]['components']['xyz']} for z in METALS}}
        rows.append({'sample':source['sample'],'role':source['role'],'source_case_id':source['source_case_id'],'preparation_status':source['status'],
            'pools':pools,'new7_minus_baseline':changes,'all_origin_relative_works_kcal_mol':works,'evidence':evidence,'failed_new_cells':fail,'classification':None})
    by={r['sample']:r for r in rows};contrasts={}
    for baseline in ('historical3','previous5'):
        contrasts[baseline]={}
        for policy in POLICIES:
            a0=by[0]['new7_minus_baseline'][baseline][policy];a3=by[3]['new7_minus_baseline'][baseline][policy]
            contrasts[baseline][policy]={'target_minus_control_change_kcal_mol':a.difference(a3['differential_change_kcal_mol'],a0['differential_change_kcal_mol'])} if a0 and a3 else None
    return {'protocol_id':PROTOCOL,'terminal':True,'declared_sources':2,'declared_new_cells':8,'reused_old_cells':20,
       'available_new_cells':sum(v['status']=='available' for r in rows for z in METALS for q,v in r['evidence'][z].items() if q.startswith('midpoint_')),
       'rows':rows,'target_minus_control_changes':contrasts,'actual_saved_algebra_checks':{'count':len(residuals),'maximum_residual':max(residuals,default=None),'tolerance_kcal_mol':h.TOL},
       'pool_orders':{'historical3':THREE,'previous5':FIVE,'new7':SEVEN},'classification':None,'new_optimizations':0,'new_DFT_calls':0,
       'interpretation':'One declared midpoint; lower energy alone does not demonstrate better classification; no automatic further search'}

def execute(path,action):
    if not os.environ.get('SLURM_JOB_ID'):raise ValueError('Actual Slurm required')
    c=read(path)
    if c['implementation']!=pin(__file__):raise ValueError('Wrong adapter snapshot')
    matrix=load('_midpoint_existing_matrix',c['matrix']);place=matrix.placement()
    base=Path(path).parent/('execution' if action=='execute' else 'preflight_'+os.environ['SLURM_JOB_ID'])
    base.mkdir(exist_ok=False);start=time.monotonic();pf=None;mr={};sr=[];receipt={'status':'failed','config':pin(path),'placement':place,'action':action}
    with (Path(path).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            sys.path.insert(0,str(ROOT/'scripts'))
            m=read(c['worker_manifest']);b=load('_midpoint_warm_cpu',m['implementation']);s=load('_midpoint_prior_mace_checks',c['previous_score_adapter'])
            a=load('_midpoint_prior_solvent',c['solvent_adapter']);h=load('_midpoint_prior_pool',c['pool_adapter']);selector=load('_midpoint_selector',c['scorer'])
            f=load('_midpoint_factorial',c['factorial']);r=load('_midpoint_component',c['component']);n=load('_midpoint_native',c['native_adapter'])
            pf=preflight(c,b,s,a,n,r,h);write(base/'PREFLIGHT.json',pf)
            if action=='execute':
                mr=run_mace(c,m,b,s,base,pf['versions'])
                # Admission, never energy/sign, controls dispatch. Failed MACE cells remain unavailable.
                selected=[t for t in c['cells'] if mr.get(t['task_id'],{}).get('status')=='complete']
                solvent_start=time.monotonic()
                if selected:sr=matrix.solvent({**c,'cells':selected},base,f,r,n,mr)
                write(base/'SOLVENT_PHASE.json',{'elapsed_seconds':time.monotonic()-solvent_start,'declared_cells':8,'eligible_cells':len(selected),'rows':sr})
            receipt['status']='complete'
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
        finally:
            try:
                if pf and action=='execute':
                    if not mr and (base/'mace/PHASE.json').exists():mr=read(base/'mace/PHASE.json')['rows']
                    if not sr:
                        for t in c['cells']:
                            p=base/'solvent/replay/tasks'/t['task_id']/'RESULT.json'
                            if p.exists():
                                v=read(p);sr.append({k:v[k] for k in ('task_id','status','reason','G_qvac_hartree')}|{'result':pin(p)})
                    write(base/'FINAL_COLLECTION.json',{'config':pin(path),**summarize(c,pf,mr,sr,a,h,selector)})
            except BaseException as e:
                receipt.update(status='failed',collection_failure={'reason':str(e),'traceback':traceback.format_exc()});raise
            finally:
                elapsed=time.monotonic()-start
                receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,
                  MACE_forwards=sum(read(p).get('actual_model_forwards',0) for p in base.glob('mace/cells/*/mace_result.json')),
                  SCC_started=sum(1 for _ in base.glob('solvent/points/*/*/START.json')),SCC_terminal=sum(1 for _ in base.glob('solvent/points/*/*/EXECUTION.json')),
                  native_component_calls=sum(1 for _ in base.glob('solvent/replay/tasks/*/*/native_EXECUTION.json')),
                  measured_resource_receipts={
                    'MACE':[pin(p) for p in sorted(base.glob('mace/lane*/CPU_WORKER_SUMMARY.json'))],
                    'SCC':[pin(p) for p in sorted(base.glob('solvent/points/*/*/EXECUTION.json'))],
                    'native_ALPB':[pin(p) for p in sorted(base.glob('solvent/replay/tasks/*/*/native_EXECUTION.json'))]})
                write(base/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('prepare');q.add_argument('--preparation',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    for name in ('preflight','execute'):
        q=sub.add_parser(name);q.add_argument('--config',type=Path,required=True);q.add_argument('--config-sha',required=True)
    args=p.parse_args()
    if args.action=='prepare':prepare(args.preparation,args.output)
    else:execute(verify({'path':str(args.config.resolve()),'sha256':args.config_sha}),args.action)
