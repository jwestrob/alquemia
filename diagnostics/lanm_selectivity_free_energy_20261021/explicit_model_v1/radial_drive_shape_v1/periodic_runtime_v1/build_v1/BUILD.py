"""Candidate-independent archive/wrapper build; never invokes a molecular program."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import socket
import subprocess
import time

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
MAXIMA={'configure':1,'archive':1,'compile':8,'link':3}

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(p):
    p=Path(p).resolve()
    with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
    return {'path':str(p),'sha256':h}
def verify(p):
    require(isinstance(p,dict) and set(p)=={'path','sha256'} and pin(p['path'])==p,'Missing or changed pin')
    return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def load(p,name):
    spec=importlib.util.spec_from_file_location(name,verify(p));x=importlib.util.module_from_spec(spec);spec.loader.exec_module(x);return x
def save(path,x):
    tmp=path.with_suffix('.json.tmp')
    with tmp.open('w') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(path)

def commands(h,out,cmake):
    """Pure command construction; the only parent recipe delta is archive workers=12."""
    out=Path(out).resolve();archive=out/'archive';objects=out/'objects';mods=objects/'mod'
    require(h['archive_build_argv'][-2:]==['--parallel','1'],'Frozen handoff parallelism')
    require(h['new_architecture_flags']==[] and len(h['compile_order'])==8,'Frozen compiler scope')
    configure=[str(cmake) if x=='cmake' else str(archive) if x=='<fresh-build-dir>' else x for x in h['archive_configuration_argv']]
    archive_command=[str(cmake) if x=='cmake' else str(archive) if x=='<fresh-build-dir>' else x for x in h['archive_build_argv']]
    archive_command[-1]='12'
    result=[{'label':'configure','kind':'configure','argv':configure},
            {'label':'archive','kind':'archive','argv':archive_command}]
    sources=h['source_inventory']['isolated_objects']
    for name in h['compile_order']:
        argv=['/usr/bin/f95',*h['wrapper_flags'],'-I'+str(archive/'mod'),'-I'+str(mods),'-J'+str(mods),
              '-c',sources[name]['path'],'-o',str(objects/(name+'.o'))]
        result.append({'label':'compile_'+name.replace('.','_'),'kind':'compile','argv':argv,'output':str(objects/(name+'.o'))})
    for role,sequence in h['link_objects_before_archive'].items():
        tail=[x.replace('<new-radial-build>',str(archive)) for x in h['link_tail']]
        argv=['/usr/bin/f95',*h['wrapper_flags'],*[str(objects/x) for x in sequence],*tail,'-o',str(out/role)]
        result.append({'label':'link_'+role,'kind':'link','role':role,'argv':argv,'output':str(out/role)})
    require(len(result)==13 and len({x['label'] for x in result})==13,'Exactly thirteen build commands')
    return result

def execute_plan(plan,result,out,invoke,persist):
    """One attempt per fixed command. A failed invocation immediately unwinds."""
    for item in plan:
        kind=item['kind'];require(kind in MAXIMA,'Command kind')
        result['attempts'][kind]+=1;require(result['attempts'][kind]<=MAXIMA[kind],'Finite command ceiling')
        result['phase']=item['label'];persist()
        archive=kind=='archive'
        invoke(result,out,item['label'],item['argv'],out,1800 if archive else 180,3*1024**3,3600 if archive else 180)
        if 'output' in item:require(Path(item['output']).is_file(),'Compiler/link did not produce its required output')

def source_context(config,release,out):
    cfg=json.loads(config.read_text());rel=json.loads(release.read_text())
    require(cfg['status']=='bound_periodic_v2_build_configuration','Unbound build')
    require(cfg['implementation']==pin(__file__),'Exact build implementation')
    require(rel['status']=='released_one_periodic_v2_build' and rel['configuration']==pin(config),'Exact build release')
    require(rel['implementation']==pin(__file__) and Path(rel['output']).resolve()==out.resolve(),'Released implementation/output')
    require(rel['maximum_archive_builds']==1 and rel['maximum_wrapper_compiles']==8 and rel['maximum_wrapper_links']==3,'Finite release')
    require(rel['maximum_native_programs']==0,'Build-only release')
    x=read(cfg['inputs']);require(x['status']=='prepared_candidate_independent_v2_build_inputs','Prepared build inputs')
    require(x['implementation']==pin(__file__) and cfg['candidate'] is None and cfg['composition'] is None,'Candidate-independent implementation')
    for p in x['sources'].values():verify(p)
    declaration=read(x['declaration']);h=read(x['handoff']);ready=read(x['runtime_ready'])
    require(declaration['parents']['BUILD_HANDOFF.json']==x['handoff'],'Declared handoff')
    require(declaration['prospective_build']['archive_workers']==cfg['archive_workers']==12,'Declared twelve archive workers')
    require(cfg['runtime_source_review']==x['runtime_source_review'] and cfg['runtime_root_admission']==x['runtime_root_admission'],'Exact runtime authorities')
    review=read(cfg['runtime_source_review']);root=read(cfg['runtime_root_admission'])
    require(review['status']=='passed_periodic_v2_runtime_source_review' and review['ready']==x['runtime_ready'],'Independent runtime source review')
    require(review['native_sources']==h['runtime_native_manifest']==ready['native_sources'],'Reviewed generated source')
    require(root['status']=='admitted_periodic_v2_runtime_source' and root['ready']==x['runtime_ready'] and root['independent_review']==cfg['runtime_source_review'],'Root runtime source admission')
    build_review=read(rel['build_source_review'])
    require(build_review['status']=='passed_periodic_v2_build_source_review' and build_review['implementation']==pin(__file__) and build_review['inputs']==cfg['inputs'],'Build adapter review')
    for p in h['source_inventory']['radial_archive_files']+list(h['source_inventory']['isolated_objects'].values())+h['compiler_and_link_dependencies']:verify(p)
    radial=read(h['radial_source_manifest']);native=read(h['runtime_native_manifest'])
    require(radial['source_outputs']==h['source_inventory']['radial_archive_files'] and native['files']==h['source_inventory']['isolated_objects'],'Complete inventory equality')
    require(len(radial['source_outputs'])==615,'Exact archive file count')
    for p in x['tools'].values():verify(p)
    require(pin('/usr/bin/f95')==h['compiler_and_link_dependencies'][0],'Compiler symlink target')
    baseline=h['old_baseline_reuse']
    for p in baseline.values():
        if isinstance(p,dict):verify(p)
    oldroot=read(baseline['root_admission']);oldbuild=read(baseline['build_result'])
    require(oldroot['result']==baseline['build_result'] and oldbuild['executables']['gradient_old']==baseline['executable'],'Admitted old executable reuse')
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread per process')
    require(not any(os.environ.get(k) for k in ('LD_PRELOAD','LD_LIBRARY_PATH')),'No loader injection')
    require(not any(os.environ.get(k) for k in ('FFLAGS','LDFLAGS','CMAKE_GENERATOR','CMAKE_TOOLCHAIN_FILE')),'No implicit build recipe override')
    return cfg,x,h

def run(config,release,out):
    out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'started','implementation':pin(__file__),'configuration':pin(config),'release':pin(release),
       'commands':{},'attempts':dict.fromkeys(MAXIMA,0),'native_programs':0,'automatic_retry':False,'candidate':None,'error':None}
    def persist():save(out/'RESULT.json',result)
    persist()
    try:
        cfg,x,h=source_context(config,release,out)
        helper=load(x['sources']['COMMON.py'],'build_common');allocation=load(x['sources']['allocation12'],'build_allocation')
        job=os.environ.get('SLURM_JOB_ID','');require(job.isdecimal() and job!='0','Actual Slurm job')
        sargv=[x['tools']['scontrol']['path'],'show','job',job,'-o']
        observed=subprocess.run(sargv,check=True,capture_output=True,text=True,timeout=15)
        require(not observed.stderr.strip(),'Scheduler readback stderr')
        alloc=allocation.parse_allocation(observed.stdout,os.environ,socket.gethostname(),os.getuid())
        require(len(os.sched_getaffinity(0))>=12,'Twelve allocated CPUs available')
        (out/'SCHEDULER_PREFLIGHT.txt').write_text(observed.stdout)
        result.update(inputs=cfg['inputs'],native_sources=h['runtime_native_manifest'],radial_source_manifest=h['radial_source_manifest'],
            allocation={**alloc,'argv':sargv,'raw':pin(out/'SCHEDULER_PREFLIGHT.txt')},phase='source_and_node_preflight',
            compiler=pin('/usr/bin/f95'),toolchain=x['tools'],compiler_and_link_dependencies=h['compiler_and_link_dependencies'])
        persist()
        (out/'objects/mod').mkdir(parents=True)
        plan=commands(h,out,x['tools']['cmake']['path']);result['planned_commands']=plan;persist()
        execute_plan(plan,result,out,helper.invoke,persist)
        executables={p['role']:pin(p['output']) for p in plan if p['kind']=='link'}
        executables['gradient_old']=h['old_baseline_reuse']['executable']
        loader=load(x['sources']['loader'],'build_native_loader');runtime={}
        for role,p in executables.items():
            label='loader_'+role
            helper.invoke(result,out,label,[x['tools']['ldd']['path'],verify(p)],out,30,3*1024**3,30)
            rec=result['commands'][label];require(not verify(rec['stderr']).read_text().strip(),'Loader stderr')
            runtime[role]=loader.native_loader(verify(rec['stdout']).read_text(),x['native_runtime'])
        require(result['attempts']==MAXIMA,'Completed exact command denominator')
        result.update(status='periodic_v2_build_complete',phase='complete',executables=executables,
           radial_library=pin(out/'archive/libtinker.a'),radial_polar_module=pin(out/'archive/mod/polar.mod'),
           old_baseline=h['old_baseline_reuse'],native_runtime=runtime,
           objects={Path(p['output']).name:pin(p['output']) for p in plan if p['kind']=='compile'},
           modules=[pin(p) for folder in (out/'archive/mod',out/'objects/mod') for p in sorted(folder.glob('*.mod'))],
           cache=pin(out/'archive/CMakeCache.txt'))
    except BaseException as exc:result.update(status='periodic_v2_build_failed',error=repr(exc))
    usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,
       'child_CPU_seconds':usage.ru_utime+usage.ru_stime,'peak_child_RSS_KiB':usage.ru_maxrss,
       'archive_build_workers':12,'wrapper_build_workers':1,'numerical_threads_per_process':1,
       'scheduler_accounting_separate':True,'peak_RSS_is_not_aggregate_job_peak':True}
    persist();return 0 if result['status']=='periodic_v2_build_complete' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();raise SystemExit(run(a.config.resolve(),a.release.resolve(),a.out.resolve()))
