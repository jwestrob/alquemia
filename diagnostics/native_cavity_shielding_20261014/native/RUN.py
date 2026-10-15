#!/usr/bin/env python3
"""Twelve declared native cavity updates, with scout gates before continuation."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import socket
import subprocess
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE=ROOT/'diagnostics/native_cavity_shielding_20261014/native'
OLD=ROOT/'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/manifest.json'
TERMS=('Gelec','Ghb','Gsasa','Gshift');BLOCKS=('real_real','real_cap','cap_cap')
PROTOCOL='native_capped_union_electrostatic_scout_v1'
TOL=1e-9;CEILING=1e-6;HA=627.509474
SCOUT=('local','union','cap_removed','translation','rotation','permutation')
IDS=('sample0_Ca','sample0_La','sample3_Ca','sample3_La')

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
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
    s=importlib.util.spec_from_file_location('_qualified_native',verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def prepare(preparation,build,output):
    prep=read(preparation);br=read(build);bm=read(br['manifest']);old=read(OLD)
    require(prep['status']=='prepared' and prep['prepared_endpoints']==4,'Four actual prepared endpoints required')
    require(br['status']=='built' and not br['engine_rebuilt'],'Qualified wrapper build required')
    require(bm['implementation']['REPLAY.f90']['sha256']==pin(HERE/'CAVITY.f90')['sha256'],'Build wrapper differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('RUN.py','RUN.sbatch'):shutil.copyfile(HERE/name,out/name)
    config={'protocol':PROTOCOL,'implementation':pin(out/'RUN.py'),'wrapper':pin(out/'RUN.sbatch'),
      'plan':pin(HERE.parent/'PLAN.md'),'preparation':pin(preparation),'build':pin(build),'binary':br['binary'],
      'native_manifest':pin(OLD),'native_adapter':old['implementation']['NATIVE.py'],'native_source':bm['implementation']['REPLAY.f90'],
      'numerical_allowance_Eh':TOL,'informative_charge_bound_Eh':CEILING,'HA_TO_KCAL':HA,
      'scout_id':IDS[0],'scout_order':SCOUT,'continuation_order':IDS[1:],'maximum_cavity_updates':12,
      'native_threads':1,'reserved_CPUs':4,'new_SCC_calls':0,'new_MACE_calls':0,'new_DFT_calls':0}
    write(out/'CONFIG.json',config);print(json.dumps(pin(out/'CONFIG.json')))

def allocation():
    require(bool(os.environ.get('SLURM_JOB_ID')),'Actual Slurm allocation required')
    command=['scontrol','show','job','-o',os.environ['SLURM_JOB_ID']]
    text=subprocess.check_output(command,text=True)
    def field(k):
        m=re.search(r'(?:^|\s)'+re.escape(k)+r'=([^ ]+)',text);require(m is not None,'Missing scheduler '+k);return m.group(1)
    cpus=int(field('NumCPUs'));require(field('JobState')=='RUNNING' and int(field('NumNodes'))==1 and cpus>=4,'One actual shared node with at least4 CPUs required')
    tres=dict(x.split('=',1) for x in field('AllocTRES').split(',') if '=' in x)
    match=re.fullmatch(r'([0-9.]+)([KMGT]?)',tres.get('mem',''));require(match is not None,'Actual allocated RAM unavailable')
    memory=int(float(match[1])*{'':1,'K':1024,'M':1024**2,'G':1024**3,'T':1024**4}[match[2]])
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        require(os.environ.get(k)=='1','Serial native/BLAS thread policy differs: '+k)
    return {'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'scheduler':text,'allocated_CPUs':cpus,
      'allocated_memory_bytes':memory,'numerical_threads':1,'affinity':sorted(os.sched_getaffinity(0))}

def variant(v):
    for key in ('xyz','charges','masks','mapping'):verify(v[key])
    lines=verify(v['xyz']).read_text().splitlines();n=int(lines[0]);require(n==v['atom_count'] and len(lines)==n+2,'Cavity atom count differs')
    atoms=[(a[0],tuple(map(float,a[1:]))) for a in (line.split() for line in lines[2:])]
    q=[float(x) for x in verify(v['charges']).read_text().split()]
    masks=verify(v['masks']).read_text().splitlines();require(int(masks[0])==n and len(masks)==n+1,'Mask count differs')
    flags=[tuple(map(int,x.split())) for x in masks[1:]]
    require(len(q)==n and all(len(x)==2 and set(x)<=set((0,1)) for x in flags),'Charge/mask schema differs')
    require(all(len(x[1])==3 and all(map(math.isfinite,x[1])) for x in atoms) and all(map(math.isfinite,q)),'Nonfinite cavity input')
    require([i for i,x in enumerate(flags) if x[0]]==v['cap_indices'] and [i for i,x in enumerate(flags) if x[1]]==v['uncertain_charge_indices'],'Mask identities differ')
    require(all(uncertain or charge==0 for charge,(_,uncertain) in zip(q,flags)),'Added charge is not exact zero')
    return atoms,q,flags

def preflight(c,n,place):
    require(c['protocol']==PROTOCOL and c['numerical_allowance_Eh']==TOL and c['informative_charge_bound_Eh']==CEILING,'Frozen protocol differs')
    for key in ('implementation','wrapper','plan','native_source'):verify(c[key])
    prep=read(c['preparation']);read(prep['manifest']);br=read(c['build']);bm=read(br['manifest']);nm=read(c['native_manifest'])
    require(br['status']=='built' and br['binary']==c['binary'] and bm['implementation']['REPLAY.f90']==c['native_source'],'Built caller differs')
    require(bm['build']==nm['build'] and bm['helper']==c['native_adapter'],'Native build/helper differs')
    verify(c['binary'])
    for p in bm['build']['pins']+br['linked_libraries']:verify(p)
    require({x['path']:x['sha256'] for x in br['linked_libraries']}=={x['path']:x['sha256'] for x in nm['build']['linked_libraries']},'Native runtime libraries differ')
    require(prep['status']=='prepared' and prep['prepared_endpoints']==4 and prep['terminal'],'Prepared source mapping unavailable')
    rows={r['id']:r for r in prep['rows']};require(set(rows)==set(IDS),'Finite source inventory differs')
    records={};sizes=[]
    for key,row in rows.items():
        require(row['status']=='prepared' and row['core_atom_count']==192 and row['cap_count']==14,'Original capped state differs')
        for p in ('state_signature','source_preparation','protonated_source'):verify(row[p])
        old=read(row['archived_native_result']);require(old['status']=='available' and old['manifest']==c['native_manifest'],'Archived native endpoint differs')
        archived=n.task_preflight(nm,old['input'])
        require(all(x['passed'] for x in n.checks(nm,old['water'],archived)) and all(x['informative'] for x in old['vacuum_precision']),'Archived native qualification failed')
        require(row['saved_solvent']['component_v2_acceptance']['accepted'],'Original component-v2 gate unavailable')
        for medium in ('water','vacuum'):
            for p in ('execution','xyz','charges','cavity','born_matrix_with_HB','result'):verify(old[medium][p])
        data={name:variant(v) for name,v in row['variants'].items()}
        local,union,removed=[data[k] for k in ('local','union','cap_removed')]
        require(row['variants']['local']['xyz']['sha256']==old['input']['xyz']['sha256'],'Local geometry does not match archived state')
        require(local[1]==[float(q) for q in verify(old['vacuum']['charges']).read_text().split()],'Saved core charges changed')
        require(len(local[0])==192 and sum(x[0] for x in local[2])==14 and all(x[1] for x in local[2]),'Core charge/cap masks differ')
        require(all(union[k][:192]==local[k] for k in range(3)),'Union does not preserve complete local core')
        require(all(x==(0,0) for x in union[2][192:]),'Environment mask differs')
        keep=[i for i,x in enumerate(union[2]) if not x[0]]
        require(all(removed[k]==[union[k][i] for i in keep] for k in range(3)),'Cap removal changes a real center/charge')
        require(len(removed[0])==len(union[0])-14,'Cap removal count differs')
        if key==IDS[0]:
            require(set(data)==set(SCOUT),'Six prescribed scout variants required')
            for name in ('translation','rotation','permutation'):
                a,q,f=data[name]
                if name=='permutation':require((a,q,f)==tuple(list(reversed(x)) for x in union),'Permutation changes aligned identity/charges/masks')
                else:
                    expected=[(z,tuple(x+t for x,t in zip(p,(11,-7,3)))) if name=='translation' else (z,(-p[1],p[0],p[2])) for z,p in union[0]]
                    require(a==expected and q==union[1] and f==union[2],'Rigid control differs from frozen transformation')
        else:require(set(data)=={'local','union','cap_removed'},'Unexpected continuation variant')
        sizes.extend(v['atom_count'] for v in row['variants'].values())
        records[key]={'source':row,'archived':old}
    maxn=max(sizes);pairs=maxn*(maxn-1)//2
    arrays={'brdr_dsdrt_bornMat':56*maxn*maxn,'nnlists':4*maxn*maxn,
      'ppind_nnlistr_ddpair':52*pairs,'linear_arrays_and_grid_allowance':1024*maxn+32768}
    dense=sum(arrays.values());required=math.ceil(1.25*dense)+1024**3
    require(required<=place['allocated_memory_bytes'],'Native allocation estimate exceeds granted memory; no dispatch')
    return records,{'status':'passed','largest_atom_count':maxn,'source_array_estimate_bytes':arrays,
      'estimated_native_bytes':dense,'required_with_25percent_and_1GiB_headroom':required,'allocated_memory_bytes':place['allocated_memory_bytes'],
      'native_neighbor_list_parallel':False,'redundant_dense_wrapper_copy':False,'records':{k:records[k]['source']['archived_native_result'] for k in IDS}}

def fields(text,label,count):
    rows=re.findall(r'^'+re.escape(label)+r'\s+(.+)$',text,re.M);require(len(rows)==1,'Missing native output '+label)
    a=list(map(float,rows[0].split()));require(len(a)==count and all(map(math.isfinite,a)),'Invalid native '+label);return a

def evaluate(c,n,row,name,out,env):
    v=row['variants'][name];d=out/'tasks'/row['id']/name;d.mkdir(parents=True)
    write(d/'START.json',{'source':row['id'],'variant':name,'input':v,'binary':c['binary'],'cavity_updates':1})
    receipt={'status':'unavailable','source':row['id'],'variant':name,'input':v,'checks':[]}
    try:
        execution=n.command([str(verify(c['binary'])),str(verify(v['xyz'])),str(verify(v['charges'])),str(verify(v['masks']))],d,'native',dict(env,XTBHOME=str(d),XTBPATH=str(d)),read(c['native_manifest'])['build']['timer']['path'])
        text=verify(execution['stdout']).read_text();require('internal GFN2-xTB/ALPB' in text,'Internal factory not recorded')
        receipt.update(components_Eh=dict(zip(TERMS,fields(text,'COMPONENTS',4))),repeat_components_Eh=dict(zip(TERMS,fields(text,'REPEAT_COMPONENTS',4))),
          cap_zero_components_Eh=dict(zip(TERMS,fields(text,'CAP_ZERO_COMPONENTS',4))),electrostatic_blocks_Eh=dict(zip(BLOCKS,fields(text,'ELECTROSTATIC_BLOCKS',3))),
          shape_blocks_Eh=dict(zip(BLOCKS,fields(text,'SHAPE_BLOCKS',3))),shape_Eh=dict(zip(('full','cap_zero'),fields(text,'SHAPE',2))),
          rounding_bounds_Eh=dict(zip(('Gelec','Ghb','Gshape'),fields(text,'ROUNDING_BOUNDS',3))),
          cap_zero_rounding_bounds_Eh=dict(zip(('Gelec','Ghb','Gshape'),fields(text,'CAP_ZERO_ROUNDING_BOUNDS',3))),
          charge_sums_e=dict(zip(('full','real','cap'),fields(text,'CHARGE_SUMS',3))),
          cavity_global=dict(zip(('aDet','keps','alpbet','dielectric'),fields(text,'CAVITY_GLOBAL',4))),
          counts=dict(zip(('atoms','uncertain','caps'),map(int,fields(text,'COUNTS',3)))),execution=pin(d/'native_EXECUTION.json'),cavity=pin(d/'cavity.tsv'))
        checks=receipt['checks']
        def check(label,residual,tolerance=TOL):checks.append({'name':label,'residual_Eh':residual,'tolerance_Eh':tolerance,'passed':abs(residual)<=tolerance})
        check('native_block_closure',math.fsum(receipt['electrostatic_blocks_Eh'].values())-receipt['components_Eh']['Gelec'])
        check('cap_zero_real_block',receipt['electrostatic_blocks_Eh']['real_real']-receipt['cap_zero_components_Eh']['Gelec'])
        check('shape_block_closure',math.fsum(receipt['shape_blocks_Eh'].values())-receipt['shape_Eh']['full'])
        check('cap_zero_shape_block',receipt['shape_blocks_Eh']['real_real']-receipt['shape_Eh']['cap_zero'])
        for term in TERMS:check('repeat_'+term,receipt['repeat_components_Eh'][term]-receipt['components_Eh'][term])
        for label in ('rounding_bounds_Eh','cap_zero_rounding_bounds_Eh'):
            for term,bound in receipt[label].items():checks.append({'name':label+'_'+term,'bound_Eh':bound,'ceiling_Eh':CEILING,'passed':0<=bound<CEILING})
        expected={'atoms':v['atom_count'],'uncertain':len(v['uncertain_charge_indices']),'caps':len(v['cap_indices'])}
        require(receipt['counts']==expected,'Native mask inventory differs')
        require(all(x['passed'] for x in checks),'Native repeat/block/precision qualification failed')
        receipt['status']='available'
    except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
    finally:write(d/'RESULT.json',receipt)
    return receipt|{'result':pin(d/'RESULT.json')}

def scout_gate(results,record,place):
    checks=[];local=results['local'];union=results['union'];removed=results['cap_removed'];old=record['archived']['vacuum']
    for term in TERMS:
        residual=local['components_Eh'][term]-old['components_Eh'][term]
        checks.append({'name':'archived_local_'+term,'residual_Eh':residual,'tolerance_Eh':TOL,
          'recorded_charge_bound_Eh':old['rounding_bounds_Eh'].get(term,0),'same_rounded_input':True,'passed':abs(residual)<=TOL})
    differences={}
    for name in ('translation','rotation','permutation'):
        v=results[name];residual=v['components_Eh']['Gelec']-union['components_Eh']['Gelec']
        checks.append({'name':name+'_Gelec','residual_Eh':residual,'tolerance_Eh':TOL,'passed':abs(residual)<=TOL})
        differences[name]={term:v['components_Eh'][term]-union['components_Eh'][term] for term in TERMS}
    require(abs(union['charge_sums_e']['real']-removed['charge_sums_e']['full'])<=1e-12,'Cap control changes real charge sum')
    for label,bound in difference_bounds(results,old).items():
        checks.append({'name':label,'bound_Eh':bound,'ceiling_Eh':CEILING,'passed':bound<CEILING})
    peak=max(read(v['execution'])['peak_RSS_KiB'] for v in results.values())
    memory_ok=peak*1024*1.25+512*1024**2<=place['allocated_memory_bytes']
    return {'status':'passed' if all(x['passed'] for x in checks) and memory_ok else 'failed','checks':checks,
      'all_component_invariance_residuals_Eh':differences,'measured_maximum_RSS_KiB':peak,
      'continuation_memory_gate':{'passed':memory_ok,'allocated_memory_bytes':place['allocated_memory_bytes'],'headroom':'25percent of measured peak plus512MiB'}}

def difference_bounds(results,old):
    union=results['union'];removed=results['cap_removed']
    return {'union_minus_local':union['rounding_bounds_Eh']['Gelec']+old['rounding_bounds_Eh']['Gelec'],
      'cap_zero_union_minus_cap_removed':union['cap_zero_rounding_bounds_Eh']['Gelec']+removed['rounding_bounds_Eh']['Gelec']}

def old_shape(old):
    execution=read(old['execution']);text=verify(execution['stdout']).read_text();a,k,b,_=fields(text,'CAVITY_GLOBAL',4)
    q=math.fsum(map(float,verify(old['charges']).read_text().split()))
    return 0.5*k*b/a*q*q

def collect(c,records,results,gate):
    rows=[]
    for key in IDS:
        row={'id':key,'status':'unavailable','reason':'scout gate blocked or native update incomplete','updates':results.get(key,{}),'classification':None}
        if key in records:
            row['archived_local']=records[key]['source']['archived_native_result']
            union=results.get(key,{}).get('union');removed=results.get(key,{}).get('cap_removed')
            if gate and gate['status']=='passed' and union and removed and all(v['status']=='available' for v in (union,removed)):
                old=records[key]['archived']['vacuum'];transfer=union['components_Eh']['Gelec']-old['components_Eh']['Gelec']
                shape=union['shape_Eh']['full']-old_shape(old);cap_volume=union['cap_zero_components_Eh']['Gelec']-removed['components_Eh']['Gelec']
                cap_shape=union['shape_Eh']['cap_zero']-removed['shape_Eh']['full']
                bounds=difference_bounds(results[key],old);informative=all(b<CEILING for b in bounds.values())
                row.update(status='available' if informative else 'unavailable',reason=None if informative else 'Propagated difference charge bound uninformative',transfer_Eh={'Gelec':transfer,'Gshape':shape,'Gelec_without_shape':transfer-shape},
                  cap_volume_Eh={'Gelec':cap_volume,'Gshape':cap_shape,'Gelec_without_shape':cap_volume-cap_shape},
                  conservative_transfer_charge_bound_Eh=bounds['union_minus_local'],
                  conservative_cap_volume_charge_bound_Eh=bounds['cap_zero_union_minus_cap_removed'])
        rows.append(row)
    paired={};paired_bounds={};by={r['id']:r for r in rows}
    for sample in (0,3):
        ca,la=[by[f'sample{sample}_{z}'] for z in ('Ca','La')]
        if ca['status']==la['status']=='available':
            paired[str(sample)]={k:{term:ca[k][term]-la[k][term] for term in ('Gelec','Gshape','Gelec_without_shape')} for k in ('transfer_Eh','cap_volume_Eh')}
            paired_bounds[str(sample)]={k:ca[k]+la[k] for k in ('conservative_transfer_charge_bound_Eh','conservative_cap_volume_charge_bound_Eh')}
        else:paired[str(sample)]=None;paired_bounds[str(sample)]=None
    contrast=None
    if all(paired.values()):contrast={k:{term:paired['3'][k][term]-paired['0'][k][term] for term in ('Gelec','Gshape','Gelec_without_shape')} for k in ('transfer_Eh','cap_volume_Eh')}
    contrast_bounds={k:paired_bounds['0'][k]+paired_bounds['3'][k] for k in paired_bounds['0']} if all(paired_bounds.values()) else None
    return {'protocol':PROTOCOL,'terminal':True,'status':'complete' if all(r['status']=='available' for r in rows) else 'incomplete',
      'config_native_source':c['native_source'],'preparation':c['preparation'],'build':c['build'],'rows':rows,'scout':gate,
      'Ca_minus_La':paired,'sample3_minus_sample0':contrast,'hartree_to_kcal_mol':HA,
      'Ca_minus_La_charge_bounds_Eh':paired_bounds,'sample3_minus_sample0_charge_bounds_Eh':contrast_bounds,
      'paired_precision_informative':{k:all(v<CEILING for v in b.values()) if b else False for k,b in paired_bounds.items()},
      'contrast_precision_informative':all(v<CEILING for v in contrast_bounds.values()) if contrast_bounds else False,
      'classification':None,'interpretation':'Electrostatic capped-union sensitivity only; no addition of full-protein SASA/HB/shift to the composite; exact-zero environment is not a chemical state'}

def execute(path):
    place=allocation();c=read(path);require(pin(__file__)==c['implementation'],'Runner snapshot differs')
    out=Path(path).parent/'execution';out.mkdir(exist_ok=False);start=time.monotonic();records={};results={};gate=None
    receipt={'status':'failed','config':pin(path),'placement':place};n=load(c['native_adapter']);previous={}
    with (Path(path).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for sig in (signal.SIGTERM,signal.SIGHUP,signal.SIGINT):
            previous[sig]=signal.signal(sig,lambda *_:n.STOP.set())
        try:
            records,pf=preflight(c,n,place);write(out/'PREFLIGHT.json',pf)
            env=dict(os.environ);key=IDS[0];results[key]={}
            for name in SCOUT:results[key][name]=evaluate(c,n,records[key]['source'],name,out,env)
            gate=scout_gate(results[key],records[key],place);write(out/'SCOUT_QUALIFICATION.json',gate)
            require(gate['status']=='passed','Scout qualification failed; continuation blocked')
            for key in IDS[1:]:
                results[key]={}
                for name in ('union','cap_removed'):results[key][name]=evaluate(c,n,records[key]['source'],name,out,env)
                require(all(v<CEILING for v in difference_bounds(results[key],records[key]['archived']['vacuum']).values()),'Propagated difference charge bound uninformative; further updates blocked')
            receipt['status']='complete'
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
        finally:
            for sig,handler in previous.items():signal.signal(sig,handler)
            try:write(out/'FINAL_COLLECTION.json',{'config':pin(path),**collect(c,records,results,gate)})
            except BaseException as e:receipt.update(status='failed',collection_failure={'reason':str(e),'traceback':traceback.format_exc()});raise
            finally:
                receipts=[pin(p) for p in sorted(out.glob('tasks/*/*/native_EXECUTION.json'))]
                measured=[read(p) for p in receipts];elapsed=time.monotonic()-start
                receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=elapsed*place['allocated_CPUs'],
                  native_update_intents=len(list(out.glob('tasks/*/*/START.json'))),native_child_execution_receipts=len(measured),successful_native_updates=sum(r['returncode']==0 for r in measured),
                  measured_child_CPU_seconds=sum(r.get('child_CPU_seconds',0) for r in measured),peak_child_RSS_KiB=max((r.get('peak_RSS_KiB',0) for r in measured),default=0),
                  resource_receipts=receipts,new_SCC_calls=0,new_MACE_calls=0,new_DFT_calls=0)
                write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('prepare');q.add_argument('--preparation',type=Path,required=True);q.add_argument('--build',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('execute');q.add_argument('--config',type=Path,required=True);q.add_argument('--config-sha',required=True)
    args=p.parse_args()
    if args.action=='prepare':prepare(args.preparation,args.build,args.output)
    else:execute(verify({'path':str(args.config.resolve()),'sha256':args.config_sha}))
