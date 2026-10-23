"""One finite Slurm link, explicit periodic-operator application and saved replay."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import socket
import subprocess
import sys
import time
from decimal import Decimal

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
         'NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','PYTHONDONTWRITEBYTECODE')
RESOURCES=dict(nodes=1,tasks=1,cpus_per_task=1,memory_grant_MiB=4000,
               address_space_bytes=3*1024**3,link_wall_seconds=180,preflight_wall_seconds=60,
               native_wall_seconds=3600,replay_wall_seconds=1800,numerical_threads=1)


def require(condition,message):
    if not condition:raise ValueError(message)


def pin(path):
    path=Path(path).resolve()
    with path.open('rb') as stream:h=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(path),sha256=h)


def verify(item):
    require(isinstance(item,dict) and set(item)=={'path','sha256'} and pin(item['path'])==item,'Missing/changed pin')
    return Path(item['path'])


def read(item):return json.loads(verify(item).read_text())


def save(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    temporary.replace(path)
    descriptor=os.open(path.parent,os.O_RDONLY)
    try:os.fsync(descriptor)
    finally:os.close(descriptor)


def parse_allocation(raw,environ,hostname,uid):
    job=environ.get('SLURM_JOB_ID','');require(job.isdecimal() and job!='0','Actual Slurm job ID')
    fields=dict(word.split('=',1) for word in raw.split() if '=' in word)
    expected={'JobId':job,'JobState':'RUNNING','Account':'banfield','QOS':'standard','Partition':'standard',
              'NumNodes':'1','NumCPUs':'1','NumTasks':'1','CPUs/Task':'1','Requeue':'0','OverSubscribe':'NO'}
    for key,value in expected.items():require(fields.get(key)==value,'Allocation '+key)
    require(fields['BatchHost'].split('.')[0]==hostname.split('.')[0],'Allocated host')
    require(fields['UserId'].endswith('('+str(uid)+')'),'Allocated user')
    require(environ.get('SLURM_CPUS_PER_TASK')=='1' and environ.get('SLURM_NTASKS')=='1','Task environment')
    require(environ.get('SLURM_JOB_NUM_NODES',environ.get('SLURM_NNODES'))=='1','Single node')
    tres=dict(x.split('=',1) for x in fields['AllocTRES'].split(','))
    require(tres['cpu']=='1' and tres['node']=='1','Allocated TRES')
    def memory(text):
        units={'K':1024,'M':1024**2,'G':1024**3,'T':1024**4}
        require(text[-1] in units,'Memory units');value=Decimal(text[:-1])*units[text[-1]]
        require(value==int(value),'Integral memory');return int(value)
    grant=memory(tres['mem'])
    require(grant==4000*1024**2 and memory(fields['MinMemoryCPU'])==4000*1024**2,'Default memory grant')
    require(RESOURCES['address_space_bytes']<grant,'Address-space cap below grant')
    return dict(job_id=job,hostname=hostname,mode='slurm',account='banfield',partition='standard',qos='standard',
                allocated_CPUs=1,nodes=1,tasks=1,memory_grant_bytes=grant,default_memory_MiB_per_CPU=4000,
                raw_scontrol=raw,requeue=False,scheduler_allocated_CPU_seconds=None)


def child_environment():
    env=dict(os.environ)
    removed={key:env.pop(key) for key in ('LD_PRELOAD','LD_LIBRARY_PATH') if key in env}
    for key in THREADS:env[key]='1'
    return env,removed


def limits(cpu_seconds):
    def apply():
        for kind,cap in ((resource.RLIMIT_AS,RESOURCES['address_space_bytes']),(resource.RLIMIT_CPU,cpu_seconds)):
            soft,hard=resource.getrlimit(kind)
            require(all(v==resource.RLIM_INFINITY or cap<=v for v in (soft,hard)),'Cannot relax inherited limit')
            resource.setrlimit(kind,(cap,cap))
    return apply


def native_loader(stdout,expected):
    require('not found' not in stdout,'Native runtime missing')
    resolved={}
    for line in stdout.splitlines():
        fields=line.strip().split()
        if not fields or fields[0].startswith('linux-vdso'):continue
        if '=>' in fields:
            require(len(fields)>=3 and fields[2].startswith('/'),'Unexpected loader record')
            resolved[fields[0]]=pin(fields[2])
        else:
            require(fields[0].startswith('/'),'Unexpected loader path')
            resolved['loader']=pin(fields[0])
    require(resolved==expected,'Native loader/runtime differs from bound parent environment')
    return resolved


def context(config,release,out):
    cfg=json.loads(config.read_text());root=json.loads(release.read_text())
    require(cfg['status']=='bound_stage1_slurm_parity_configuration','Unbound configuration')
    require(cfg['implementation']==pin(__file__) and cfg['resources']==RESOURCES,'Exact runner/resources')
    require(root['status']=='released_one_slurm_periodic_operator_parity','Root release status')
    for key,expected in (('configuration',pin(config)),('implementation',pin(__file__)),('stage1_ready',cfg['stage1_ready'])):
        require(root[key]==expected,'Release '+key)
    require(root['resources']==RESOURCES and root['maximum_native_invocations']==1,'Root finite limits')
    require(Path(root['output']).resolve()==out.resolve(),'Released output')
    ready=read(cfg['stage1_ready']);source=read(ready['manifest']);prep=read(ready['preparation'])
    require(ready['status']=='source_ready_stage1_complete_periodic_operator_parity','Stage1 readiness')
    require(ready['normal_guards_unchanged'] and prep['counts']==dict(native_programs=1,complete_operator=30,
            separate_reciprocal=30,total_reciprocal=60,direct_initializations=1,induced_solves=0,energies=0,gradients=0),
            'Stage1 finite scope')
    require(ready['execution']==dict(native_programs=1,complete_operator_applications=30,
            separate_reciprocal_applications=30,total_reciprocal_applications=60,direct_initializations=1,
            induced_solves=0,energies=0,gradients=0,factorizations=0),'Declared action counts')
    for item in source['sources'].values():verify(item)
    for key in ('implementation','setup','replay','link_recipe','numerical_contract'):verify(ready[key])
    review=read(cfg['stage1_source_review'])
    require(review['status']=='passed_stage1_complete_periodic_operator_source_review','Source review status')
    for key,expected in (('ready',cfg['stage1_ready']),('implementation',ready['implementation']),
                         ('replay',ready['replay']),('preparation',ready['preparation'])):
        require(review[key]==expected,'Source review '+key)
    policy=read(root['adapter_source_review'])
    require(policy['status']=='passed_stage1_slurm_adapter_source_review' and policy['implementation']==pin(__file__),
            'Adapter review status/source')
    require(policy['configuration']==pin(config),'Adapter reviewed configuration')
    static=read(ready['static_input']);admission=read(prep['static_root_admission'])
    require(ready['candidate']==prep['candidate']==static['candidate'],'Exact old candidate')
    require(ready['static_input']==prep['static_input'],'Exact job20053 input')
    require(root['candidate']==ready['candidate'] and root['static_input']==ready['static_input'],'Released identity')
    require(admission==read(cfg['static_root_admission']) and cfg['static_root_admission']==prep['static_root_admission'],
            'Exact static admission authority')
    require(admission['status']=='admitted_nominal_saved_La_static_numerical_test' and admission['job_id']==20053
            and admission['sufficient_test_passed'] is True and admission['input_manifest']==ready['static_input'],
            'Successful exact static20053 admission')
    env=read(cfg['environment'])
    require(env['status']=='passed_exact_isolated_environment_audit','Admitted replay environment')
    require(Path(sys.prefix).resolve()==Path(env['prefix']).resolve(),'Exact Python environment')
    for item in env['runtime_files']:verify(item)
    recipe=read(ready['link_recipe'])
    require(recipe['library']==ready['library']==prep['library'],'Immutable archive')
    for item in [recipe['library'],*recipe['source'],*recipe['module_files'],*recipe['compiler_and_link_dependencies']]:verify(item)
    for item in cfg['native_runtime'].values():verify(item)
    verify(cfg['preflight']);verify(cfg['wrapper']);verify(cfg['loader_command'])
    require(cfg['scheduler_command']=='/usr/local/bin/scontrol','Read-only scheduler command')
    require(all(os.environ.get(key)=='1' for key in THREADS),'One-thread environment')
    return cfg,ready,prep,recipe,env


def run(config,release,out):
    tick,cpu=time.monotonic(),time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    out.mkdir(parents=True,exist_ok=False)
    result=dict(status='started',phase='context',implementation=pin(__file__),configuration=pin(config),release=pin(release),
                started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),hostname=socket.gethostname(),PID=os.getpid(),
                link_attempts=0,links_launched=0,native_attempts=0,native_invocations=0,replay_attempts=0,
                induced_solves=0,energies=0,gradients=0,factorizations=0,automatic_retry=False,commands={},error=None,
                e_app=None,uniform_application_bound=False,energy_force_qualified=False,physical_accuracy_qualified=False)
    def persist():save(out/'RESULT.json',result)
    def invoke(name,argv,cwd,wall_seconds,env):
        row=dict(status='intent_durable_before_launch',argv=list(map(str,argv)),cwd=str(cwd),PID=None,
                 wall_limit_seconds=wall_seconds,address_space_limit_bytes=RESOURCES['address_space_bytes'],
                 CPU_limit_seconds=wall_seconds,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        result['commands'][name]=row;persist()
        start=time.monotonic();prior=resource.getrusage(resource.RUSAGE_CHILDREN);process=None
        try:
            with (out/(name+'.out')).open('x') as stdout,(out/(name+'.err')).open('x') as stderr:
                process=subprocess.Popen(argv,cwd=cwd,env=env,stdout=stdout,stderr=stderr,
                                         start_new_session=True,preexec_fn=limits(wall_seconds))
                row.update(status='launched',PID=process.pid)
                if name=='native':result['native_invocations']+=1
                if name=='link':result['links_launched']+=1
                persist()
                try:code=process.wait(timeout=wall_seconds)
                except subprocess.TimeoutExpired:
                    row['wall_timeout']=True;os.killpg(process.pid,signal.SIGKILL);code=process.wait()
            row.update(status='complete' if code==0 else 'failed',returncode=code)
        except BaseException as error:
            row.update(status='failed',error=repr(error))
            if process is not None and process.poll() is None:
                os.killpg(process.pid,signal.SIGKILL);row['returncode']=process.wait()
            raise
        finally:
            after=resource.getrusage(resource.RUSAGE_CHILDREN)
            row['cost']=dict(wall_seconds=time.monotonic()-start,
                             child_user_CPU_seconds=after.ru_utime-prior.ru_utime,
                             child_system_CPU_seconds=after.ru_stime-prior.ru_stime,
                             aggregate_child_peak_RSS_so_far_KiB=after.ru_maxrss,
                             per_phase_peak_not_inferred=True)
            row['stdout']=pin(out/(name+'.out')) if (out/(name+'.out')).exists() else None
            row['stderr']=pin(out/(name+'.err')) if (out/(name+'.err')).exists() else None
            persist()
        require(row['status']=='complete',name+' failed; no retry')
        return row
    def stopped(signum,frame):raise RuntimeError('Received signal '+str(signum))
    signal.signal(signal.SIGTERM,stopped);signal.signal(signal.SIGINT,stopped)
    try:
        soft,hard=resource.getrlimit(resource.RLIMIT_AS);cap=RESOURCES['address_space_bytes']
        require(all(v==resource.RLIM_INFINITY or cap<=v for v in (soft,hard)),'Controller inherited memory limit')
        resource.setrlimit(resource.RLIMIT_AS,(cap,cap))
        persist();cfg,ready,prep,recipe,environment=context(config,release,out)
        result.update(stage1_ready=cfg['stage1_ready'],candidate=ready['candidate'],static_input=ready['static_input'],
                      static_root_admission=prep['static_root_admission'],preparation=ready['preparation'],
                      stage1_source_review=cfg['stage1_source_review'],library=recipe['library'],resources=RESOURCES)
        env,removed=child_environment();result['removed_child_loader_overrides']=removed
        job=os.environ.get('SLURM_JOB_ID','');require(job.isdecimal() and job!='0','Actual Slurm job ID')
        command=[cfg['scheduler_command'],'show','job',job,'-o']
        receipt=subprocess.run(command,capture_output=True,text=True,timeout=15,check=True)
        require(not receipt.stderr.strip(),'Scheduler stderr')
        allocation=parse_allocation(receipt.stdout,os.environ,socket.gethostname(),os.getuid())
        allocation['argv']=command;allocation['affinity_CPUs']=len(os.sched_getaffinity(0))
        require(allocation['affinity_CPUs']>=1,'CPU affinity');result['allocation']=allocation;persist()
        result['phase']='node_environment'
        invoke('preflight',[sys.executable,cfg['preflight']['path'],cfg['environment']['path'],str(out/'NODE_PREFLIGHT.json')],
               out,RESOURCES['preflight_wall_seconds'],env)
        node=json.loads((out/'NODE_PREFLIGHT.json').read_text())
        require(node['status']=='passed_stage1_node_runtime_preflight','Node preflight')
        result['node_preflight']=pin(out/'NODE_PREFLIGHT.json')
        result['phase']='single_link';result['link_attempts']+=1;require(result['link_attempts']==1,'One link attempt')
        executable=Path(recipe['output']);executable.parent.mkdir(parents=True,exist_ok=False)
        require(recipe['argv'][-2:]==['-o',str(executable)],'Exact output binding')
        invoke('link',recipe['argv'],executable.parent,RESOURCES['link_wall_seconds'],env)
        result['executable']=pin(executable)
        invoke('loader',[cfg['loader_command']['path'],str(executable)],out,RESOURCES['preflight_wall_seconds'],env)
        result['actual_native_runtime']=native_loader((out/'loader.out').read_text(),cfg['native_runtime'])
        folder=out/'native';folder.mkdir()
        for name,item in prep['raw_inputs'].items():
            shutil.copyfile(verify(item),folder/name)
            require(pin(folder/name)['sha256']==item['sha256'],'Exact input copy '+name)
        result['raw_inputs']={name:pin(folder/name) for name in prep['raw_inputs']}
        result['phase']='explicit_operator';result['native_attempts']+=1
        require(result['native_attempts']==1 and result['native_invocations']==0,'One native attempt')
        verify(result['executable'])
        invoke('native',[str(executable),'system.xyz'],folder,RESOURCES['native_wall_seconds'],env)
        for name,item in prep['raw_inputs'].items():require(pin(folder/name)['sha256']==item['sha256'],'Post-run input '+name)
        result['raw_outputs']={name:pin(folder/name) for name in
                              ('parity_state.dat','parity_apply.dat','parity_receipt.dat','pme_norm_export_001.dat')}
        result['phase']='saved_replay';result['replay_attempts']+=1
        require(result['replay_attempts']==1,'One saved replay')
        invoke('replay',[sys.executable,ready['replay']['path'],ready['preparation']['path'],str(folder),str(out/'replay')],
               out,RESOURCES['replay_wall_seconds'],env)
        comparison=json.loads((out/'replay/RESULT.json').read_text());result['replay_result']=pin(out/'replay/RESULT.json')
        require(comparison['status']=='passed_finite_complete_periodic_operator_parity','Saved parity status')
        require(comparison['implementation']==ready['replay'] and comparison['preparation']==ready['preparation'],
                'Saved parity source/preparation')
        require(comparison['candidate']==ready['candidate'] and comparison['static_input']==ready['static_input'],
                'Saved parity identity')
        require(comparison['counts']==prep['counts'] and result['native_invocations']==1,'Finite action denominator')
        result.update(status='stage1_periodic_operator_parity_complete',phase='terminal',sampled_operator_parity_passed=True)
    except BaseException as error:
        result.update(status='stage1_periodic_operator_parity_failed',error=repr(error),sampled_operator_parity_passed=False)
        if (out/'replay/RESULT.json').exists():result['replay_result']=pin(out/'replay/RESULT.json')
    finally:
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        result['cost']=dict(wall_seconds=time.monotonic()-tick,controller_CPU_seconds=time.process_time()-cpu,
                            child_user_CPU_seconds=after.ru_utime-before.ru_utime,child_system_CPU_seconds=after.ru_stime-before.ru_stime,
                            aggregate_child_peak_RSS_KiB=after.ru_maxrss,controller_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                            numerical_threads=1,concurrent_numerical_workers=1,scheduler_accounting_separate=True)
        persist()
    print(result['status']);return 0 if result['status']=='stage1_periodic_operator_parity_complete' else 1


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--release',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();raise SystemExit(run(args.config.resolve(),args.release.resolve(),args.out.resolve()))
