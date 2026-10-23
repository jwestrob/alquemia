"""Finite Stage2 artifact/runtime joins; no scientific computation at import."""
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import signal
import socket
import subprocess
import sys
import time

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
def require(ok,message):
    if not ok: raise ValueError(message)
def pin(p):
    p=Path(p).resolve()
    with p.open('rb') as f: digest=hashlib.file_digest(f,'sha256').hexdigest()
    return {'path':str(p),'sha256':digest}
def verify(p):
    require(isinstance(p,dict) and set(p)=={'path','sha256'} and pin(p['path'])==p,'missing/changed pin '+str(p))
    return Path(p['path'])
def read(p): return json.loads(verify(p).read_text())
def load(p,name):
    spec=importlib.util.spec_from_file_location(name,verify(p));module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def save(p,x):
    p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp')
    def scalar(y):
        if type(y).__module__=='numpy' and hasattr(y,'item'):return y.item()
        raise TypeError(type(y).__name__)
    with tmp.open('w') as f:json.dump(x,f,indent=2,allow_nan=False,default=scalar);f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(p)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def environment(blas=1):
    out=dict(os.environ)
    for k in THREADS:out[k]='1'
    out['OPENBLAS_NUM_THREADS']=str(blas);out['PYTHONDONTWRITEBYTECODE']='1'
    for k in ('LD_LIBRARY_PATH','LD_PRELOAD'):out.pop(k,None)
    return out
def set_limits(memory,cpu):
    for kind,value in ((resource.RLIMIT_AS,memory),(resource.RLIMIT_CPU,cpu)):
        old=resource.getrlimit(kind);require(all(x==resource.RLIM_INFINITY or value<=x for x in old),'cannot relax inherited resource cap')
        resource.setrlimit(kind,(value,value))
def invoke(result,out,label,argv,cwd,wall,memory,cpu,blas=1):
    require(label not in result['commands'],'one attempt per command')
    rec={'status':'intent_durable_before_launch','argv':list(map(str,argv)),'cwd':str(cwd),'wall_seconds_limit':wall,'AS_bytes':memory,'CPU_seconds_limit':cpu,'PID':None,'started_utc':now()}
    result['commands'][label]=rec;save(out/'RESULT.json',result)
    tick=time.monotonic();u0=resource.getrusage(resource.RUSAGE_CHILDREN);proc=None
    try:
        with (out/(label+'.out')).open('x') as stdout,(out/(label+'.err')).open('x') as stderr:
            proc=subprocess.Popen(rec['argv'],cwd=cwd,env=environment(blas),stdout=stdout,stderr=stderr,start_new_session=True,preexec_fn=lambda:set_limits(memory,cpu))
            rec.update(status='launched',PID=proc.pid);save(out/'RESULT.json',result)
            try:code=proc.wait(timeout=wall)
            except subprocess.TimeoutExpired:
                rec['wall_timeout']=True;os.killpg(proc.pid,signal.SIGKILL);code=proc.wait()
            rec.update(status='complete' if code==0 else 'failed',returncode=code)
    except BaseException as exc:
        rec.update(status='failed',error=repr(exc))
        if proc is not None and proc.poll() is None:os.killpg(proc.pid,signal.SIGKILL);rec['returncode']=proc.wait()
        raise
    finally:
        u1=resource.getrusage(resource.RUSAGE_CHILDREN)
        rec['cost']={'wall_seconds':time.monotonic()-tick,'child_user_CPU_seconds':u1.ru_utime-u0.ru_utime,'child_system_CPU_seconds':u1.ru_stime-u0.ru_stime,'aggregate_child_peak_RSS_KiB':u1.ru_maxrss,'per_phase_peak_not_inferred':True}
        rec['stdout']=pin(out/(label+'.out')) if (out/(label+'.out')).exists() else None
        rec['stderr']=pin(out/(label+'.err')) if (out/(label+'.err')).exists() else None
        save(out/'RESULT.json',result)
    require(rec['status']=='complete',label+' failed; no retry')
    return rec

def allocation(manifest,cpus):
    helper=load(manifest['helpers']['slurm1' if cpus==1 else 'slurm12'],'stage2_allocation')
    job=os.environ.get('SLURM_JOB_ID','');require(job.isdecimal() and job!='0','actual Slurm job required')
    argv=['/usr/local/bin/scontrol','show','job',job,'-o'];raw=subprocess.run(argv,check=True,capture_output=True,text=True,timeout=15)
    require(not raw.stderr.strip(),'scontrol stderr')
    rec=helper.parse_allocation(raw.stdout,os.environ,socket.gethostname(),os.getuid());rec['argv']=argv
    require(len(os.sched_getaffinity(0))>=cpus,'CPU affinity');return rec

def source_context(config,release,out,implementation,kind,cpus,blas=1):
    cfg=json.loads(Path(config).read_text());rel=json.loads(Path(release).read_text())
    require(cfg['status']=='bound_stage2_'+kind+'_configuration','unbound configuration')
    require(cfg['implementation']==pin(implementation) and rel['implementation']==pin(implementation),'implementation pin')
    require(rel['status']=='released_stage2_'+kind and rel['configuration']==pin(config),'exact root release')
    require(Path(rel['output']).resolve()==Path(out).resolve(),'exact released output')
    manifest=read(cfg['manifest'])
    require(manifest['status']=='prepared_stage2_execution_inputs_no_numerical_calls','manifest status')
    for item in manifest['execution_sources'].values():verify(item)
    for item in manifest['helpers'].values():verify(item)
    review=read(rel['source_review'])
    require(review['status']=='passed_stage2_execution_source_review' and review['manifest']==cfg['manifest'],'source review binding')
    require(review['execution_sources']==manifest['execution_sources'],'reviewed execution sources')
    require(review['implementations'][kind]==cfg['implementation'],'reviewed exact phase implementation')
    stage1=read(manifest['stage1_actual_admission'])
    require(stage1['status']=='admitted_finite_periodic_operator_parity' and stage1['sampled_operator_parity_passed'] is True,'actual Stage1 parity')
    require(stage1['candidate']==manifest['candidate'] and stage1['static_input']==manifest['central_static_input'],'exact old center')
    require(read(manifest['native_source_review'])['status']=='root_reviewed_stage2_native_source_design','native source review')
    env=read(manifest['environment']);require(Path(sys.prefix).resolve()==Path(env['prefix']).resolve(),'exact isolated interpreter')
    for item in env['runtime_files']:verify(item)
    require(all(os.environ.get(k)==(str(blas) if k=='OPENBLAS_NUM_THREADS' else '1') for k in THREADS),'fixed thread environment')
    alloc=allocation(manifest,cpus)
    if kind!='build':
        build=read(cfg['build']);require(build['status']=='stage2_build_complete' and build['manifest']==cfg['manifest'],'exact successful build')
        require(build['source_manifest']==manifest['native_sources'],'exact native source generation')
        if cfg['state_id'] not in ('central','BYPASS_OLD_ABSENT','BYPASS_CLONE_ABSENT','BYPASS_CLONE_ZERO'):
            scout=read(cfg['central_scout_review'])
            require(scout['status']=='admitted_stage2_central_and_bypass_scout' and scout['build']==cfg['build'] and scout['manifest']==cfg['manifest'],'central/bypass scout gate')
            require(scout['all_three_bypass_tokens_equal'] is True and scout['central_dual_energy_passed'] is True,'scout actual gates')
    return cfg,manifest,rel,alloc

def native_build(cfg,manifest):
    b=read(cfg['build']);require(b['status']=='stage2_build_complete' and b['manifest']==cfg['manifest'],'bound build')
    for item in b['executables'].values():verify(item)
    return b

def initialization(cfg,manifest):
    x=read(cfg['initializer']);row=manifest['rows'][cfg['state_id']]
    require(row['role']=='ACTIVE' and cfg['state_id']!='central','only changed active initialization')
    require(x['status']=='stage2_initialization_complete' and x['state_id']==cfg['state_id'],'matching initializer')
    require(x['inputs']==row['inputs'] and x['build']==cfg['build'] and x['manifest']==cfg['manifest'],'initializer inputs/build')
    require(x['counts']=={'native_programs':1,'direct_initializations':1,'mutual_applications':0,'induced_solves':0,'gradient_calls':0,'energy_only_calls':0},'initializer denominator')
    require(x['PME_checks']['site_checks']['passed'] and x['PME_checks']['summary_comparison']['passed'],'PME admission')
    verify(x['state_output']);verify(x['pme_export']);return x

def certificate(cfg,manifest):
    c=read(cfg['certificate']);row=manifest['rows'][cfg['state_id']]
    require(c['status']=='stage2_state_nominal_certificate_passed' and c['state_id']==cfg['state_id'],'state-specific certificate')
    require(c['inputs']==row['inputs'] and c['manifest']==cfg['manifest'],'certificate literal input/manifest')
    require(c['candidate']==manifest['candidate'] and c['e_app'] is None and c['application_error_bound_established'] is False,'nominal scope')
    require(c['sufficient_numerical_test_passed'] is True,'positive sufficient test')
    verify(c['expected_state']);verify(c['pme_export'])
    if cfg['state_id']=='central':
        actual=read(manifest['stage1_actual_admission'])
        require(c['root_admission']==manifest['stage1_actual_admission'] and c['expected_state']==actual['expected_central_state'] and c['pme_export']==actual['initialized_PME'],'central actual parity reuse')
        require(c['factor_root_admission']==actual['static_root_admission'],'central factor reuse')
    else:
        require(c['build']==cfg['build'],'certificate native build')
        i=read(c['initializer']);a=read(c['assembly']);f=read(c['factor'])
        require(i['state_id']==a['state_id']==f['state_id']==cfg['state_id'],'certificate phase identity')
        require(i['inputs']==c['inputs'] and i['state_output']==c['expected_state'] and i['pme_export']==c['pme_export'],'certificate raw state')
        require(a['status']=='stage2_assembly_complete' and f['status']=='passed_nominal_static_sufficient_numerical_test','certificate phases complete')
        require(f['assembly']==c['assembly'] and f['sufficient_numerical_test_passed'] is True,'factor exact assembly')
        require(c['nominal']==f['operator'],'exact nominal shift')
        admission=read(cfg['certificate_admission'])
        require(admission['status']=='accepted_stage2_exact_state_certificate' and admission['certificate']==cfg['certificate'],'root certificate admission')
    return c
