"""Synthetic allocation and pre-call failure controls, no scheduler/native/link."""
import ast
import copy
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import time
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('stage1_slurm',HERE/'RUN.py')
M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)


def run(out):
    tick,cpu=time.monotonic(),time.process_time();checks=[]
    def check(label,ok):
        checks.append(dict(label=label,passed=bool(ok)))
        if not ok:raise ValueError(label)
    def reject(label,call):
        try:call()
        except (ValueError,KeyError):check(label,True)
        else:check(label,False)
    raw=('JobId=123 JobState=RUNNING Account=banfield QOS=standard_expand Partition=standard '
         'NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 Requeue=0 OverSubscribe=NO '
         'BatchHost=node-64-768g-4 UserId=jwestrob(601) AllocTRES=cpu=1,mem=4000M,node=1,billing=1 MinMemoryCPU=4000M')
    env=dict(SLURM_JOB_ID='123',SLURM_CPUS_PER_TASK='1',SLURM_NTASKS='1',SLURM_JOB_NUM_NODES='1')
    actual=M.parse_allocation(raw,env,'node-64-768g-4',601)
    check('Exact1CPU/defaultmemory allocation fixture',actual['memory_grant_bytes']==4000*1024**2 and actual['allocated_CPUs']==1)
    for before,after in [('JobId=123','JobId=124'),('RUNNING','PENDING'),('Account=banfield','Account=other'),
                         ('QOS=standard_expand','QOS=standard'),('Partition=standard','Partition=other'),
                         ('NumNodes=1','NumNodes=2'),('NumCPUs=1','NumCPUs=2'),('NumTasks=1','NumTasks=2'),
                         ('CPUs/Task=1','CPUs/Task=12'),('Requeue=0','Requeue=1'),('OverSubscribe=NO','OverSubscribe=YES'),
                         ('mem=4000M','mem=3000M'),('MinMemoryCPU=4000M','MinMemoryCPU=8000M')]:
        reject('Reject '+before,lambda before=before,after=after:M.parse_allocation(raw.replace(before,after),env,'node-64-768g-4',601))
    reject('Reject wrong node',lambda:M.parse_allocation(raw,env,'biotite',601))
    reject('Reject wrong UID',lambda:M.parse_allocation(raw,env,'node-64-768g-4',602))
    reject('Reject missing real job environment',lambda:M.parse_allocation(raw,{},'node-64-768g-4',601))
    with patch.dict(os.environ,{'LD_PRELOAD':'manufactured.so','LD_LIBRARY_PATH':'/manufactured','OMP_NUM_THREADS':'8'}):
        child,removed=M.child_environment()
        check('Owned child loader overrides removed',not {'LD_PRELOAD','LD_LIBRARY_PATH'}&set(child) and len(removed)==2)
        check('One thread in each child',all(child[k]=='1' for k in M.THREADS))
        check('Parent loader environment unchanged',os.environ['LD_PRELOAD']=='manufactured.so')
    fake={'liba.so':{'path':'/fake/a','sha256':'a'},'loader':{'path':'/fake/ld','sha256':'b'}}
    with patch.object(M,'pin',side_effect=lambda path:fake['loader' if path=='/fake/ld' else 'liba.so']):
        check('Exact native loader fixture accepted',M.native_loader('linux-vdso.so.1 (x)\nliba.so => /fake/a (x)\n/fake/ld (x)\n',fake)==fake)
        reject('Missing loader dependency rejected',lambda:M.native_loader('liba.so => not found\n',fake))
        reject('Unexpected loader library rejected',lambda:M.native_loader('other.so => /fake/a (x)\n/fake/ld (x)\n',fake))
    # Exercise terminal failure before allocation/application; all subprocesses forbidden.
    with tempfile.TemporaryDirectory() as temporary:
        root=Path(temporary);config=root/'config';release=root/'release'
        config.write_text('{}');release.write_text('{}')
        with patch.object(M,'context',side_effect=ValueError('manufactured prerelease failure')), \
             patch.object(M.subprocess,'run',side_effect=AssertionError('unexpected scheduler')), \
             patch.object(M.subprocess,'Popen',side_effect=AssertionError('unexpected subprocess')), \
             patch.object(M.resource,'getrlimit',return_value=(M.resource.RLIM_INFINITY,M.resource.RLIM_INFINITY)), \
             patch.object(M.resource,'setrlimit'),patch.object(M.signal,'signal'):
            code=M.run(config,release,root/'out')
        result=json.loads((root/'out/RESULT.json').read_text())
        check('Pre-call failure retained as terminal',code==1 and result['status']=='stage1_periodic_operator_parity_failed')
        check('Pre-call failure has zero invocation counters',all(result[k]==0 for k in ('link_attempts','links_launched','native_attempts','native_invocations','replay_attempts')))
    tree=ast.parse((HERE/'RUN.py').read_text())
    run_node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run')
    invocations=[n for n in ast.walk(run_node) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='invoke']
    tags=[n.args[0].value for n in invocations]
    check('Exactly one link/native/replay call site',tags==['preflight','link','loader','native','replay'])
    invoke=next(n for n in run_node.body if isinstance(n,ast.FunctionDef) and n.name=='invoke')
    check('Intent persisted before Popen',ast.unparse(invoke).index('persist()')<ast.unparse(invoke).index('subprocess.Popen'))
    check('No scientific optimizer/factor APIs',not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
          and n.func.attr in ('factorize','minimize','solve','eigh') for n in ast.walk(tree)))
    script=(HERE/'RUN.sbatch').read_text()
    check('Default scheduler memory/time',all(flag not in script for flag in ('--mem','--time','--exclusive','--oversubscribe')))
    check('Exact1CPU standard_expand wrapper',all(flag in script for flag in ('--cpus-per-task=1','--nodes=1','--ntasks=1','--qos=standard_expand','--no-requeue')))
    check('No multithread BLAS in wrapper','OPENBLAS_NUM_THREADS=1' in script and 'OMP_NUM_THREADS=1' in script)
    for name in ('RUN.py','PREFLIGHT.py','TEST_POLICY.py'):
        ast.parse((HERE/name).read_text());check('Syntax '+name,True)
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    result=dict(status='passed_stage1_slurm_policy_controls',implementation=M.pin(HERE/'RUN.py'),test_source=M.pin(__file__),
                wrapper=M.pin(HERE/'RUN.sbatch'),checks=checks,check_count=len(checks),native_calls=0,builds=0,
                scheduler_calls=0,submissions=0,operator_applications=0,actual_box_reads=0,
                cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status'],len(checks))

if __name__=='__main__':
    import sys
    run(sys.argv[1])
