#!/usr/bin/env python3
"""Sixteen fixed native updates; reuse the completed twelve-update origin scout."""
import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import signal
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE=ROOT/'diagnostics/native_cavity_accommodation_20261014/native'
ORIGIN=ROOT/'workspaces/native_cavity_shielding_20261014/native_v1'
DONORS=('adaptive_Ca','adaptive_La');METALS=('Ca','La');SAMPLES=(0,3)
IDS=tuple(f'sample{s}_{d}_{z}' for s in SAMPLES for d in DONORS for z in METALS)
PROTOCOL='native_capped_union_accommodation_v1'
TOL=1e-9;CEILING=1e-6

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])!=p:raise ValueError('Changed artifact: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def require(ok,reason):
    if not ok:raise ValueError(reason)
def load(p):
    s=importlib.util.spec_from_file_location('_qualified_cavity_driver',verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def prepare(preparation,output):
    p=read(preparation);oc=read(ORIGIN/'CONFIG.json');old=read(ORIGIN/'execution/FINAL_COLLECTION.json')
    require(p['status']=='prepared' and len(p['rows'])==8 and {r['id'] for r in p['rows']}==set(IDS),'Exact eight prepared donor endpoints required')
    require(old['status']=='complete' and old['config']==pin(ORIGIN/'CONFIG.json'),'Completed original scout required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('RUN.py','RUN.sbatch'):shutil.copyfile(HERE/name,out/name)
    c={k:oc[k] for k in ('build','binary','native_manifest','native_adapter','native_source','numerical_allowance_Eh','informative_charge_bound_Eh','HA_TO_KCAL')}
    c.update(protocol=PROTOCOL,implementation=pin(out/'RUN.py'),wrapper=pin(out/'RUN.sbatch'),plan=pin(HERE.parent/'PLAN.md'),
      preparation=pin(preparation),origin_config=pin(ORIGIN/'CONFIG.json'),origin_driver=oc['implementation'],
      origin_collection=pin(ORIGIN/'execution/FINAL_COLLECTION.json'),origin_execution=pin(ORIGIN/'execution/EXECUTION.json'),
      origin_qualification=pin(ORIGIN/'execution/SCOUT_QUALIFICATION.json'),task_order=IDS,
      maximum_cavity_updates=16,reused_origin_updates=12,native_threads=1,reserved_CPUs=4,new_SCC_calls=0,new_MACE_calls=0,new_DFT_calls=0)
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def preflight(c,h,n,place):
    require(c['protocol']==PROTOCOL and c['numerical_allowance_Eh']==TOL and c['informative_charge_bound_Eh']==CEILING,'Frozen protocol differs')
    for key in ('implementation','wrapper','plan'):verify(c[key])
    oc=read(c['origin_config'])
    for key in ('build','binary','native_manifest','native_adapter','native_source','numerical_allowance_Eh','informative_charge_bound_Eh','HA_TO_KCAL'):
        require(c[key]==oc[key],'Qualified runtime or policy differs: '+key)
    require(oc['implementation']==c['origin_driver'],'Origin runtime implementation differs')
    origins,origin_pf=h.preflight(oc,n,place)
    old=read(c['origin_collection']);execution=read(c['origin_execution']);gate=read(c['origin_qualification'])
    require(old['terminal'] and old['status']=='complete' and old['config']==c['origin_config'],'Origin collection unavailable')
    require(gate['status']=='passed' and all(v['passed'] for v in gate['checks']) and gate['continuation_memory_gate']['passed'],'Origin scout gate failed')
    require(execution['status']=='complete' and execution['config']==c['origin_config'] and execution['successful_native_updates']==12,'Actual twelve-update origin execution unavailable')
    for p in execution['resource_receipts']:verify(p)
    oldrows={r['id']:r for r in old['rows']}
    for key,row in oldrows.items():
        require(row['status']=='available' and row['archived_local']==origins[key]['source']['archived_native_result'],'Origin identity differs')
        for result in row['updates'].values():
            saved=read(result['result']);require(saved['status']=='available' and all(x['passed'] for x in saved['checks']),'Origin native update failed')
            require(all(saved[k]==result[k] for k in ('input','components_Eh','cap_zero_components_Eh','rounding_bounds_Eh','shape_Eh')),'Origin result collection mismatch')
            for p in ('execution','cavity'):verify(saved[p])
    prep=read(c['preparation']);read(prep['manifest']);rows={r['id']:r for r in prep['rows']}
    require(prep['status']=='prepared' and prep['terminal'] and set(rows)==set(IDS),'Eight prepared endpoints unavailable')
    nm=read(c['native_manifest']);records={};paired={};memory=[]
    for key in IDS:
        row=rows[key];origin_id=f"sample{row['sample']}_{row['metal']}";donor=row['donor_candidate']
        require(key==f"sample{row['sample']}_{donor}_{row['metal']}" and donor in DONORS and row['origin_row_id']==origin_id,'Source/candidate identity differs')
        origin=origins[origin_id]['source']
        require(row['origin_preparation']==oc['preparation'] and row['status']=='prepared' and row['core_atom_count']==192 and row['cap_count']==14,'Origin preparation/core linkage differs')
        for field in ('state_signature','source_preparation','protonated_source'):
            require(row[field]==origin[field],'Physical source/state changed: '+field);verify(row[field])
        old_native=read(row['archived_native_result'])
        require(old_native['status']=='available' and old_native['manifest']==c['native_manifest'],'Saved donor native result unavailable')
        archived=n.task_preflight(nm,old_native['input'])
        require(all(x['passed'] for x in n.checks(nm,old_native['water'],archived)) and all(x['informative'] for x in old_native['vacuum_precision']),'Saved donor native component gate failed')
        require(row['saved_solvent']['component_v2_acceptance']['accepted'],'Saved donor component-v2 acceptance missing')
        for medium in ('water','vacuum'):
            for field in ('execution','xyz','charges','cavity','born_matrix_with_HB','result'):verify(old_native[medium][field])
        require(set(row['variants'])=={'local','union','cap_removed'},'Unexpected new variant')
        local,union,removed=[h.variant(row['variants'][name]) for name in ('local','union','cap_removed')]
        original_union=h.variant(origin['variants']['union'])
        require(row['variants']['local']['xyz']['sha256']==old_native['input']['xyz']['sha256'],'Saved donor geometry changed')
        require(local[1]==list(map(float,verify(old_native['vacuum']['charges']).read_text().split())),'Saved donor charges changed')
        require(len(local[0])==192 and sum(f[0] for f in local[2])==14 and all(f[1] for f in local[2]),'Local cap/uncertainty mapping differs')
        require(local[2]==original_union[2][:192],'Cap/real identities changed from origin')
        require(all(union[k][:192]==local[k] for k in range(3)),'Union changed the saved donor core')
        require(all(union[k][192:]==original_union[k][192:] for k in range(3)),'Appended physical source centers or exact-zero charges changed')
        keep=[i for i,f in enumerate(union[2]) if not f[0]]
        require(all(removed[k]==[union[k][i] for i in keep] for k in range(3)),'Cap removal changed a physical center/charge')
        require(len(union[0])==len(original_union[0]),'Whole-source cavity size changed')
        paired[row['sample'],donor,row['metal']]=union[0]
        records[key]={'source':row,'archived':old_native,'origin':oldrows[origin_id]}
        memory.append(len(union[0]))
    for sample in SAMPLES:
        for donor in DONORS:
            ca,la=[paired[sample,donor,z] for z in METALS]
            require(ca[0][0]=='Ca' and la[0][0]=='La' and ca[0][1]==la[0][1] and ca[1:]==la[1:],'Paired metal physical geometry differs')
    require(max(memory)==origin_pf['largest_atom_count'],'New allocation differs from qualified origin size')
    measured=execution['peak_child_RSS_KiB']
    require(measured*1024*1.25+512*1024**2<=place['allocated_memory_bytes'],'Actual origin memory plus headroom exceeds allocation')
    return records,{'status':'passed','declared_endpoints':8,'cavity_updates':16,'origin_qualification':c['origin_qualification'],
      'origin_execution':c['origin_execution'],'qualified_memory_estimate':origin_pf,'actual_origin_peak_RSS_KiB':measured,
      'unchanged_appended_physical_source':True,'paired_geometry_identical':True}

def collect(c,h,records,results):
    rows=[]
    for key in IDS:
        row={'id':key,'status':'unavailable','reason':'native update incomplete','updates':results.get(key,{}),'classification':None}
        if key in records:
            rec=records[key];source=rec['source'];old=rec['archived']['vacuum'];origin=rec['origin']
            row.update(sample=source['sample'],metal=source['metal'],donor_candidate=source['donor_candidate'],candidate=source['donor_candidate'],
              archived_local=source['archived_native_result'],origin_row_id=source['origin_row_id'])
            if rec.get('qualification'):row['endpoint_qualification']=rec['qualification']
            r=results.get(key,{})
            if all(name in r and r[name]['status']=='available' for name in ('union','cap_removed')):
                union,removed=[r[x] for x in ('union','cap_removed')];bounds=h.difference_bounds(r,old)
                T=union['components_Eh']['Gelec']-old['components_Eh']['Gelec'];shape=union['shape_Eh']['full']-h.old_shape(old)
                cap=union['cap_zero_components_Eh']['Gelec']-removed['components_Eh']['Gelec'];capshape=union['shape_Eh']['cap_zero']-removed['shape_Eh']['full']
                informative=all(v<CEILING for v in bounds.values())
                row.update(status='available' if informative else 'unavailable',reason=None if informative else 'Propagated endpoint bound uninformative',
                  transfer_Eh={'Gelec':T,'Gshape':shape,'Gelec_without_shape':T-shape},
                  cap_volume_Eh={'Gelec':cap,'Gshape':capshape,'Gelec_without_shape':cap-capshape},endpoint_charge_bounds_Eh=bounds)
                row['accommodation_work_change_Eh']={term:row['transfer_Eh'][term]-origin['transfer_Eh'][term] for term in row['transfer_Eh']}
                row['cap_volume_work_change_Eh']={term:row['cap_volume_Eh'][term]-origin['cap_volume_Eh'][term] for term in row['cap_volume_Eh']}
                row['work_change_charge_bounds_Eh']={'Gelec':bounds['union_minus_local']+origin['conservative_transfer_charge_bound_Eh'],
                  'cap_volume':bounds['cap_zero_union_minus_cap_removed']+origin['conservative_cap_volume_charge_bound_Eh']}
                row['work_change_precision_informative']=all(v<CEILING for v in row['work_change_charge_bounds_Eh'].values())
        rows.append(row)
    by={r['id']:r for r in rows};paired={};contrasts={}
    for donor in DONORS:
        paired[donor]={}
        for sample in SAMPLES:
            ca,la=[by[f'sample{sample}_{donor}_{z}'] for z in METALS]
            if ca['status']==la['status']=='available':
                p={name:{term:ca[name][term]-la[name][term] for term in ca[name]} for name in ('accommodation_work_change_Eh','cap_volume_work_change_Eh')}
                p['charge_bounds_Eh']={term:ca['work_change_charge_bounds_Eh'][term]+la['work_change_charge_bounds_Eh'][term] for term in ('Gelec','cap_volume')}
                p['precision_informative']=all(v<CEILING for v in p['charge_bounds_Eh'].values());paired[donor][str(sample)]=p
            else:paired[donor][str(sample)]=None
        a,b=paired[donor]['0'],paired[donor]['3']
        if a and b:
            contrasts[donor]={name:{term:b[name][term]-a[name][term] for term in a[name]} for name in ('accommodation_work_change_Eh','cap_volume_work_change_Eh')}
            contrasts[donor]['charge_bounds_Eh']={term:a['charge_bounds_Eh'][term]+b['charge_bounds_Eh'][term] for term in ('Gelec','cap_volume')}
            contrasts[donor]['precision_informative']=all(v<CEILING for v in contrasts[donor]['charge_bounds_Eh'].values())
        else:contrasts[donor]=None
    return {'protocol':PROTOCOL,'terminal':True,'status':'complete' if all(r['status']=='available' for r in rows) else 'incomplete',
      'preparation':c['preparation'],'origin_collection':c['origin_collection'],'build':c['build'],'rows':rows,
      'Ca_minus_La_accommodation_changes':paired,'sample3_minus_sample0':contrasts,'hartree_to_kcal_mol':c['HA_TO_KCAL'],
      'interpretation':'Matched electrostatic exclusion changes only; cap controls are algebraic, not physical electronic states. No classifier or full-solvation correction applied.',
      'classification':None}

def execute(path):
    require(bool(os.environ.get('SLURM_JOB_ID')),'Actual Slurm required')
    c=read(path);require(c['implementation']==pin(__file__),'Runner snapshot differs');h=load(c['origin_driver']);place=h.allocation();n=h.load(c['native_adapter'])
    out=Path(path).parent/'execution';out.mkdir(exist_ok=False);records={};results={};start=time.monotonic();previous={}
    receipt={'status':'failed','config':pin(path),'placement':place}
    with (Path(path).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for sig in (signal.SIGTERM,signal.SIGHUP,signal.SIGINT):previous[sig]=signal.signal(sig,lambda *_:n.STOP.set())
        try:
            records,pf=preflight(c,h,n,place);write(out/'PREFLIGHT.json',pf)
            env=dict(os.environ);env['LD_LIBRARY_PATH']=read(c['native_manifest'])['build']['library_directory']
            for key in IDS:
                results[key]={}
                for name in ('union','cap_removed'):results[key][name]=h.evaluate(c,n,records[key]['source'],name,out,env)
                bounds=h.difference_bounds(results[key],records[key]['archived']['vacuum'])
                peak=max(read(v['execution'])['peak_RSS_KiB'] for v in results[key].values())
                memory_ok=peak*1024*1.25+512*1024**2<=place['allocated_memory_bytes']
                qualified=all(b<CEILING for b in bounds.values()) and memory_ok
                write(out/'tasks'/key/'QUALIFICATION.json',{'status':'passed' if qualified else 'failed',
                  'difference_bounds_Eh':bounds,'informative_ceiling_Eh':CEILING,'measured_peak_RSS_KiB':peak,
                  'memory_gate':{'passed':memory_ok,'allocated_memory_bytes':place['allocated_memory_bytes'],'headroom':'25percent of measured peak plus512MiB'},
                  'first_new_endpoint':key==IDS[0],'native_checks_passed':True})
                records[key]['qualification']=pin(out/'tasks'/key/'QUALIFICATION.json')
                require(qualified,'Endpoint precision/memory qualification failed; further updates blocked')
            receipt['status']='complete'
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
        finally:
            for sig,handler in previous.items():signal.signal(sig,handler)
            try:write(out/'FINAL_COLLECTION.json',{'config':pin(path),**collect(c,h,records,results)})
            except BaseException as e:receipt.update(status='failed',collection_failure={'reason':str(e),'traceback':traceback.format_exc()});raise
            finally:
                pins=[pin(p) for p in sorted(out.glob('tasks/*/*/native_EXECUTION.json'))];measured=[read(p) for p in pins];elapsed=time.monotonic()-start
                receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=elapsed*place['allocated_CPUs'],
                  native_update_intents=len(list(out.glob('tasks/*/*/START.json'))),native_child_execution_receipts=len(pins),successful_native_updates=sum(r['returncode']==0 for r in measured),
                  measured_child_CPU_seconds=sum(r.get('child_CPU_seconds',0) for r in measured),peak_child_RSS_KiB=max((r.get('peak_RSS_KiB',0) for r in measured),default=0),
                  resource_receipts=pins,new_SCC_calls=0,new_MACE_calls=0,new_DFT_calls=0)
                write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('prepare');q.add_argument('--preparation',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('execute');q.add_argument('--config',type=Path,required=True);q.add_argument('--config-sha',required=True)
    args=p.parse_args()
    if args.action=='prepare':prepare(args.preparation,args.output)
    else:execute(verify({'path':str(args.config.resolve()),'sha256':args.config_sha}))
