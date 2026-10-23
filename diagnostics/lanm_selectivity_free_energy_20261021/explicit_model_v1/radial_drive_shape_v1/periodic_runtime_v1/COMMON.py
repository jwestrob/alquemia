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

# All inherited persistence, invocation, allocation and cap helpers above are exact.
DECLARATION=json.loads((Path(__file__).parent/'ROOT_DECLARATION.json').read_text())
IDENTITY_PIN=DECLARATION['parents']['identity_identity_validator']
if 'IDENTITY' in sys.modules:
    require(pin(sys.modules['IDENTITY'].__file__)==IDENTITY_PIN,'Existing IDENTITY import differs')
    I=sys.modules['IDENTITY']
else:
    I=load(IDENTITY_PIN,'IDENTITY');sys.modules['IDENTITY']=I
import POLICY as policy

def identity_state(m):
    require(m['identity_sources']=={k:DECLARATION['parents'][name] for k,name in
        [('IDENTITY.py','identity_identity_validator'),('STATE.py','identity_parser')]},'Frozen identity sources')
    require(pin(I.__file__)==m['identity_sources']['IDENTITY.py'],'Exact imported IDENTITY')
    return load(m['identity_sources']['STATE.py'],'periodic_v2_saved_state')

def source_context(config,release,out,implementation,kind,cpus,blas=1):
    cfg=json.loads(Path(config).read_text());rel=json.loads(Path(release).read_text())
    require(cfg['status']=='bound_periodic_v2_'+kind+'_configuration','Unbound V2 configuration')
    require(cfg['implementation']==pin(implementation) and rel['implementation']==pin(implementation),'Implementation pin')
    require(rel['status']=='released_periodic_v2_'+kind and rel['configuration']==pin(config),'Exact root release')
    require(Path(rel['output']).resolve()==Path(out).resolve(),'Exact released output')
    m=read(cfg['manifest']);require(m['status']=='bound_periodic_v2_runtime_inputs','Bound manifest')
    policy.model(m);identity_state(m)
    for p in m['execution_sources'].values():verify(p)
    for p in m['helpers'].values():verify(p)
    review=read(rel['source_review'])
    require(review['status']=='passed_periodic_v2_runtime_source_review','Runtime source review')
    require(review['execution_sources']==m['execution_sources'] and review['native_sources']==m['native_sources'],'Reviewed source bundle')
    require(review['implementations'][kind]==cfg['implementation'],'Reviewed phase implementation')
    require(m['native_source_admission']==DECLARATION['parents']['radial_port_source_review'],'Radial source authority')
    require(m['identity_source_admission']==DECLARATION['parents']['identity_root_review'],'V2 source authority')
    verify(m['native_source_admission']);verify(m['identity_source_admission'])
    env=read(m['environment']);require(Path(sys.prefix).resolve()==Path(env['prefix']).resolve(),'Exact interpreter')
    for p in env['runtime_files']:verify(p)
    require(all(os.environ.get(k)==(str(blas) if k=='OPENBLAS_NUM_THREADS' else '1') for k in THREADS),'Thread environment')
    alloc=allocation(m,cpus)
    build=native_build(cfg,m)
    sid=cfg['state_id'];row=m['rows'][sid]
    operation='certificate' if kind=='certificate' else cfg['operation']
    needs=policy.phase_needs(operation,sid,row['role'])
    if needs['central_parity']:
        p=read(cfg['central_parity_admission']);policy.parity_scope(p,m,cfg['build'])
        actual=read(p['result'])
        require(actual['status']=='periodic_v2_parity_complete' and actual['campaign_identity']==p['campaign_identity'],'Actual parity result')
        require(actual['identity']==p['central_identity'] and actual['certificate_admission']==p['certificate_admission'],'Actual central parity identity')
        replay=read(actual['replay']);require(replay['status']=='passed_finite_complete_periodic_operator_parity','Actual parity replay')
    if needs['central_scout']:
        policy.scout_scope(read(cfg['central_scout_review']),m,cfg['build'],cfg['central_parity_admission'])
    return cfg,m,rel,alloc

def native_build(cfg,m):
    b=read(cfg['build'])
    require(b['status']=='periodic_v2_build_complete' and b['native_sources']==m['native_sources'],'Exact V2 build/source')
    admission=read(cfg['build_admission'])
    require(admission['status']=='admitted_periodic_v2_build' and admission['result']==cfg['build'],'Explicit build admission')
    require(set(b['executables'])=={'initialize','parity','gradient_clone','gradient_old'},'Four executable roles')
    for p in list(b['executables'].values())+[b['radial_library'],b['radial_polar_module']]:verify(p)
    # The old baseline must carry its own matching module/library/gate/parser.
    old=b['old_baseline']
    require(old['executable']==b['executables']['gradient_old'],'Old baseline executable lineage')
    require(old['ticket_schema']=='STAGE2_BOUND_INPUT_V1' and old['state_schema']=='PERIODIC_PARITY_STATE_V1','Old baseline grammar')
    for k in ('library','polar_module','gate','reporter','parser','root_admission'):verify(old[k])
    return b

def initialization(cfg,m):
    x=read(cfg['initializer']);row=m['rows'][cfg['state_id']]
    require(row['role']=='ACTIVE','Active initializer required')
    require(x['status']=='stage2_initialization_complete' and x['state_id']==cfg['state_id'],'Matching initializer')
    require(x['inputs']==row['inputs'] and x['build']==cfg['build'] and x['manifest']==cfg['manifest'],'Initializer inputs/build')
    require(x['counts']=={'native_programs':1,'direct_initializations':1,'mutual_applications':0,'induced_solves':0,'gradient_calls':0,'energy_only_calls':0},'Initializer denominator')
    require(x['PME_checks']['site_checks']['passed'] and x['PME_checks']['summary_comparison']['passed'],'Actual PME readback')
    expected=policy.full_identity(m,cfg['build'],cfg['state_id'],x['state_output'],x['pme_export'])
    I.validate_identity(x['identity'],expected);verify(x['state_output']);verify(x['pme_export'])
    return x

def certificate(cfg,m):
    c=read(cfg['certificate']);sid=cfg['state_id'];row=m['rows'][sid]
    require(row['role']=='ACTIVE','Active certificate role')
    require(c['manifest']==cfg['manifest'] and c['inputs']==row['inputs'] and c['build']==cfg['build'],'Certificate scope')
    expected=policy.full_identity(m,cfg['build'],sid,c['expected_state'],c['pme_export'])
    pins={k:c[k] for k in ('initializer','assembly','factor')}
    pins.update(certificate=cfg['certificate'],root_admission=cfg['certificate_admission'])
    I.validate_bound_chain(expected,pins)
    i=read(c['initializer'])
    require(i['implementation']==m['execution_sources']['NATIVE.py'] and i['native_attempts']==1,'Exact initializer producer/attempt')
    require(i['commands']['native']['returncode']==0,'Successful initializer native call')
    return c
