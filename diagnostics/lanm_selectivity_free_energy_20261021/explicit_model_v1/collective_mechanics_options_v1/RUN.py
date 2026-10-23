#!/usr/bin/env python3
"""Exact16 reference/model trimer execution; old scientific helpers remain frozen."""
import argparse
import concurrent.futures
import fcntl
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shutil
import socket
import subprocess
import threading
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('trimer_preparation_helper',HERE/'PREPARE.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
pin=p.pin;verify=p.verify;read=p.read;require=p.require;module=p.module
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')


def save(path,value):
    path=Path(path);tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)


def default_config():
    prep=pin(p.W/'explicit_model_v1/collective_mechanics_options_v1/preparation_v1/RESULT.json')
    m=read(read(prep)['manifest']); prior=read(m['sources']['model_result'])
    return {'protocol':'consumed_trimer_execution_v1','implementation':pin(__file__),
        'execution_specification':pin(HERE/'EXECUTION_SPEC.md'),'preparation_source':pin(HERE/'PREPARE.py'),
        'preparation':prep,'manifest':read(prep)['manifest'],'declaration':m['declaration'],
        'common_adapter':m['sources']['common_adapter'],'common_panel':m['sources']['common_panel'],
        'model_adapter':m['sources']['model_adapter'],'native_helper':prior['config']['native_helper'],
        'induction_helper':prior['config']['induction_helper'],'orca':m['orca'],'analyze':m['model_analyze'],
        'reference_wrapper':pin(p.D/'reference_chemistry_v1/common_basis_decomposition_v1/RUN.sbatch'),
        'model_wrapper':pin(p.D/'explicit_model_v1/saved_redistribution_counterpart_v1/RUN.sbatch'),
        'maximum_reference_programs':16,'maximum_model_programs':16,'reference_ranks':4,'model_workers':4}


def context(cfg):
    require(cfg['protocol']=='consumed_trimer_execution_v1' and cfg['implementation']==pin(__file__) and
            cfg['preparation_source']==pin(HERE/'PREPARE.py') and cfg['execution_specification']==pin(HERE/'EXECUTION_SPEC.md'),'Execution source differs')
    require(cfg['maximum_reference_programs']==cfg['maximum_model_programs']==16 and cfg['reference_ranks']==cfg['model_workers']==4,'Finite scope differs')
    for x in cfg.values():
        if isinstance(x,dict) and set(x)=={'path','sha256'}:verify(x)
    prep=read(cfg['preparation']);m=read(cfg['manifest'])
    require(prep['status']=='prepared_unreleased' and prep['manifest']==cfg['manifest'] and prep['new_native_calls']==0 and
            m['declaration']==cfg['declaration'] and m['implementation']==cfg['preparation_source'],'Preparation differs')
    require(len(m['reference_cells'])==len(m['model_rows'])==16 and len(m['cycles'])==8,'Manifest denominator differs')
    require(m['sources']['common_adapter']==cfg['common_adapter'] and m['sources']['common_panel']==cfg['common_panel'] and
            m['sources']['model_adapter']==cfg['model_adapter'] and m['orca']==cfg['orca'] and m['model_analyze']==cfg['analyze'],'Helper/executable authority differs')
    common=module(cfg['common_adapter'],'trimer_execution_common');parent=read(cfg['common_panel'])
    for name,path in common.HELPERS.items():require(parent[name]==pin(path),'Native reference helper differs')
    for assets in parent['assets'].values():
        for x in assets.values():verify(x)
    for c in m['reference_cells']:require(verify(c['template_input']).read_text()==common.deck(parent,c,4,3000),'Prepared reference input differs')
    for row in m['model_rows'].values():
        for x in row['input'].values():verify(x)
    return m,common,parent,module(cfg['model_adapter'],'trimer_execution_model'),module(cfg['induction_helper'],'trimer_execution_induction')


def release(cfg,m,r):
    require(r['status']=='approved_consumed_trimer_energy_moment_diagnostic' and r['implementation']==cfg['implementation'] and
            read(r['configuration'])==cfg and r['manifest']==cfg['manifest'] and r['preparation']==cfg['preparation'] and r['declaration']==cfg['declaration'],'Root release differs')
    require(r['reference_cell_ids']==[c['cell_id'] for c in m['reference_cells']] and r['model_cell_ids']==list(m['model_rows']),'Released IDs differ')
    require(r['maximum_reference_programs']==r['maximum_model_programs']==16 and r['reference_ranks']==r['model_workers']==4 and
            r['fit_authorized'] is False and r['force_or_FD_authorized'] is False and r['finite_field_authorized'] is False,'Released scientific scope differs')
    for key in ('preparation_review','execution_review'):
        review=read(r[key]);require(review['status'].startswith('passed'),'Independent review unavailable')
    require(Path(r['reference_execution_root']).is_absolute() and Path(r['model_execution_root']).is_absolute() and
            r['reference_execution_root']!=r['model_execution_root'],'Distinct absolute output roots required')


def runtime(kind):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm allocation required')
    for k in THREADS:os.environ[k]='1'
    os.environ.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE')
    ranks=int(os.environ.get('SLURM_NTASKS','1'));cpus=int(os.environ.get('SLURM_CPUS_PER_TASK','0'))
    require((ranks,cpus)==((4,1) if kind=='reference' else (1,4)),'Actual CPU allocation differs')


def reference_claim(root,out,cell_id):
    require(out.resolve()==root.resolve()/cell_id,'Reference attempt path differs')
    root.mkdir(parents=True,exist_ok=True)
    with (root/'.admission.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        require(not (root/'STOP.json').exists(),'Reference global stop set')
        out.mkdir(exist_ok=False)


def reference_stop(root,value):
    with (root/'.admission.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if not (root/'STOP.json').exists():save(root/'STOP.json',value)


def admit_properties(folder,records):
    expected=pin(Path(folder)/'endpoint.property.txt')
    require([x for x in records if Path(x['path']).name=='endpoint.property.txt']==[expected],
            'Exact native property file missing/duplicated/changed')
    for x in records:
        require(Path(x['path']).parent.resolve()==Path(folder).resolve(),'Property file is outside native directory')
        verify(x)
    return expected


def run_reference(cfg,release_path,cid,out):
    m,common,parent,_,_=context(cfg);r=read(pin(release_path));release(cfg,m,r);runtime('reference')
    cells=[c for c in m['reference_cells'] if c['cell_id']==cid];require(len(cells)==1,'Unknown reference cell');c=cells[0]
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or 4*int(os.environ.get('SLURM_MEM_PER_CPU','0'))
    require(mem>0,'Actual granted memory required');maxcore=math.floor(.75*mem/4);require(maxcore>0,'Invalid MaxCore')
    root=Path(r['reference_execution_root']);reference_claim(root,out,cid)
    (out/'endpoint.inp').write_text(common.deck(parent,c,4,maxcore))
    command=[cfg['orca']['path'],'endpoint.inp','--bind-to none'];env=os.environ.copy()
    env['PATH']=str(Path(cfg['orca']['path']).parent)+os.pathsep+env.get('PATH','')
    start={'cell':c,'cell_id':cid,'adapter':cfg['implementation'],'config':r['configuration'],'manifest':cfg['manifest'],
           'release':pin(release_path),'input':pin(out/'endpoint.inp'),'command':command,'ranks':4,'memory_MiB':mem,
           'maxcore_MiB':maxcore,'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'threads':{k:env[k] for k in THREADS},
           'maximum_cell_native_programs':1}
    save(out/'START.json',start);t=time.monotonic();code=None;parsed=None;error=None;properties=[]
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:
            proc=subprocess.run(command,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False);code=proc.returncode
        parsed=common.parse((out/'endpoint.out').read_text(errors='replace'),c,code)
        require(parsed['status']=='complete' and all(parsed['checks'].values()),'Reference native admission failed')
        properties=[pin(f) for f in sorted(out.glob('*.property*')) if f.is_file()]
        admit_properties(out,properties)
    except Exception as exc:
        error=repr(exc);reference_stop(root,{'cell_id':cid,'error':error,'start':pin(out/'START.json')})
    elapsed=time.monotonic()-t
    result={'status':'complete' if error is None else 'failed','adapter':cfg['implementation'],'config':r['configuration'],
            'manifest':cfg['manifest'],'release':pin(release_path),'start':pin(out/'START.json'),'cell':c,'cell_id':cid,
            'input':start['input'],'output':pin(out/'endpoint.out') if (out/'endpoint.out').exists() else None,
            'stderr':pin(out/'endpoint.err') if (out/'endpoint.err').exists() else None,'returncode':code,'parsed':parsed,
            'native_program_invocations':1,'elapsed_seconds':elapsed,'allocated_elapsed_CPU_seconds_proxy':4*elapsed,
            'ranks':4,'job_id':start['job_id'],'error':error,
            'property_files':properties}
    save(out/'RESULT.json',result);return result


class Counter:
    def __init__(self):self.lock=threading.RLock();self.stop=threading.Event();self.count=0
    def reserve(self):
        with self.lock:
            require(not self.stop.is_set() and self.count<16,'Global stop or16-program ceiling')
            self.count+=1


def run_model(cfg,release_path,out):
    m,_,_,model,induction=context(cfg);r=read(pin(release_path));release(cfg,m,r);runtime('model')
    require(out.resolve()==Path(r['model_execution_root']).resolve(),'Model attempt path differs');out.mkdir(parents=True,exist_ok=False)
    native=module(cfg['native_helper'],'trimer_execution_native');count=Counter();tick=time.monotonic();cpu=time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN)
    result={'status':'running','adapter':cfg['implementation'],'config':r['configuration'],'manifest':cfg['manifest'],
            'release':pin(release_path),'job_id':os.environ['SLURM_JOB_ID'],'rows':{k:{'id':k,'status':'not_attempted'} for k in m['model_rows']},
            'native_program_invocations':0,'maximum_native_programs':16,'global_stop':False,'error':None}
    def persist():
        with count.lock:
            result.update(native_program_invocations=count.count,global_stop=count.stop.is_set());save(out/'RESULT.json',result)
    def worker(cid):
        if count.stop.is_set():return
        row=m['model_rows'][cid];folder=out/'cells'/cid;record={'id':cid,'status':'started','input':row['input'],'receipts':[]}
        try:
            folder.mkdir(parents=True,exist_ok=False)
            for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]:shutil.copyfile(verify(row['input'][k]),folder/name)
            record['runtime_input']={k:pin(folder/name) for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]}
            with count.lock:
                if count.stop.is_set():return
                count.reserve();save(out/'CALL_COUNT.json',{'native_program_invocations':count.count,'maximum':16});save(folder/'START.json',dict(record,job_id=result['job_id'],adapter=cfg['implementation'],manifest=cfg['manifest'],release=pin(release_path)))
            native.invoke([cfg['analyze']['path'],'system.xyz','EM'],folder,folder/'analyze.out',record['receipts'])
            parsed=model.parse_moments(folder/'analyze.out',model.induced_file(folder/'analyze.out'),row,induction.induction)
            record.update(status='complete',parsed_moments=parsed,start=pin(folder/'START.json'))
        except Exception as exc:
            count.stop.set();record.update(status='failed',error=repr(exc))
        with count.lock:
            result['rows'][cid]=record;save(folder/'RESULT.json',record)
            if record['status']=='failed':result['error']=record['error']
            persist()
    persist()
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(worker,m['model_rows']))
    ch=resource.getrusage(resource.RUSAGE_CHILDREN)
    result.update(status='complete' if all(x['status']=='complete' for x in result['rows'].values()) else 'failed',
        cost={'wall_seconds':time.monotonic()-tick,'parent_CPU_seconds':time.process_time()-cpu,
              'child_user_CPU_seconds':ch.ru_utime-ch0.ru_utime,'child_system_CPU_seconds':ch.ru_stime-ch0.ru_stime,
              'workers':4,'native_threads_per_worker':1})
    persist();return result


def admit_reference(cfg,manifest,result_path):
    m,common,parent,_,_=context(cfg);require(m==manifest,'Manifest argument differs');path=Path(result_path);x=json.loads(path.read_text())
    r=read(x['release']);release(cfg,m,r);c=next(c for c in m['reference_cells'] if c['cell_id']==x['cell_id'])
    require(x['status']=='complete' and x['error'] is None and x['returncode']==0 and x['native_program_invocations']==1 and
            x['adapter']==cfg['implementation'] and x['config']==r['configuration'] and x['manifest']==cfg['manifest'] and x['cell']==c,'Reference result not admitted')
    require(path.resolve()==Path(r['reference_execution_root']).resolve()/c['cell_id']/'RESULT.json','Reference result path differs')
    s=read(x['start']);require(x['start']==pin(path.parent/'START.json') and s['input']==x['input'] and s['cell']==c and s['release']==x['release'] and s['manifest']==cfg['manifest'] and s['config']==x['config'] and s['adapter']==cfg['implementation'],'Reference START differs')
    require(s['command']==[cfg['orca']['path'],'endpoint.inp','--bind-to none'] and s['ranks']==4 and s['maximum_cell_native_programs']==1 and
            s['maxcore_MiB']==math.floor(.75*s['memory_MiB']/4) and all(s['threads'][k]=='1' for k in THREADS),'Reference invocation differs')
    require(verify(x['input']).read_text()==common.deck(parent,c,4,s['maxcore_MiB']),'Executed reference deck differs');verify(x['stderr'])
    parsed=common.parse(verify(x['output']).read_text(),c,x['returncode'])
    require(parsed==x['parsed'] and parsed['status']=='complete' and all(parsed['checks'].values()),'Reference raw replay differs')
    admit_properties(path.parent,x['property_files'])
    return dict(x,result=pin(path),new_native_calls=0)


def admit_model(cfg,manifest,result):
    m,_,_,model,induction=context(cfg);require(m==manifest,'Manifest argument differs');r=read(result['release']);release(cfg,m,r)
    require(result['status']=='complete' and result['error'] is None and not result['global_stop'] and result['native_program_invocations']==16 and
            result['maximum_native_programs']==16 and result['adapter']==cfg['implementation'] and result['config']==r['configuration'] and
            result['manifest']==cfg['manifest'] and set(result['rows'])==set(m['model_rows']),'Model result incomplete/different')
    rows={}
    for cid,x in result['rows'].items():
        row=m['model_rows'][cid];folder=Path(r['model_execution_root'])/'cells'/cid
        require(x['status']=='complete' and x['input']==row['input'] and len(x['receipts'])==1,'Model row unavailable')
        start=read(x['start']);require(x['start']==pin(folder/'START.json') and start['input']==row['input'] and start['runtime_input']==x['runtime_input'] and start['release']==result['release'] and start['manifest']==cfg['manifest'] and start['adapter']==cfg['implementation'],'Model START differs')
        for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]:
            require(x['runtime_input'][k]==pin(folder/name) and verify(x['runtime_input'][k]).read_bytes()==verify(row['input'][k]).read_bytes(),'Model runtime input differs')
        rec=x['receipts'][0];require(rec['argv']==[cfg['analyze']['path'],'system.xyz','EM'] and rec['cwd']==str(folder) and rec['returncode']==0 and all(v=='1' for v in rec['threads'].values()) and rec['output']==pin(folder/'analyze.out'),'Model native receipt differs')
        parsed=model.parse_moments(folder/'analyze.out',model.induced_file(folder/'analyze.out'),row,induction.induction)
        require(parsed==x['parsed_moments'],'Model raw replay differs');rows[cid]=dict(x,new_native_calls=0)
    return rows


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['write-config','reference','model']);ap.add_argument('--config',required=True);ap.add_argument('--release');ap.add_argument('--cell');ap.add_argument('--out',type=Path)
    a=ap.parse_args()
    if a.mode=='write-config':p.write(Path(a.config),default_config())
    else:
        require(a.release is not None and a.out is not None,'Root release and output required');cfg=json.loads(Path(a.config).read_text())
        result=run_reference(cfg,a.release,a.cell,a.out.resolve()) if a.mode=='reference' else run_model(cfg,a.release,a.out.resolve())
        print(json.dumps({'status':result['status'],'native_program_invocations':result['native_program_invocations']}))
        raise SystemExit(0 if result['status']=='complete' else 1)
