"""Finite allocation-policy fixtures; no scheduler or numerical function called."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('slurm_factor_source',HERE/'RUN.py');M=importlib.util.module_from_spec(s);s.loader.exec_module(M)

def run(out):
    tick,cpu=time.perf_counter(),time.process_time();checks=[]
    def check(name,x):
        assert x,name;checks.append(name)
    def reject(name,fn):
        try:fn()
        except (AssertionError,KeyError):checks.append(name);return
        raise AssertionError(name)
    raw='JobId=123 JobState=RUNNING Account=banfield QOS=standard_expand Partition=standard NumNodes=1 NumCPUs=12 NumTasks=1 CPUs/Task=12 Requeue=0 OverSubscribe=NO BatchHost=node-64-768g-4 UserId=jwestrob(601) AllocTRES=cpu=12,mem=48000M,node=1,billing=12 MinMemoryCPU=4000M'
    env=dict(SLURM_JOB_ID='123',SLURM_CPUS_PER_TASK='12',SLURM_NTASKS='1',SLURM_JOB_NUM_NODES='1')
    r=M.parse_allocation(raw,env,'node-64-768g-4',601)
    check('exact real-profile fixture accepted',r['memory_grant_bytes']==48000*1024**2 and r['allocated_CPUs']==12 and r['scheduler_allocated_CPU_seconds'] is None)
    for before,after in [('JobId=123','JobId=124'),('RUNNING','PENDING'),('Account=banfield','Account=other'),('QOS=standard_expand','QOS=standard'),('Partition=standard','Partition=other'),('NumNodes=1','NumNodes=2'),('NumCPUs=12','NumCPUs=24'),('NumTasks=1','NumTasks=2'),('CPUs/Task=12','CPUs/Task=6'),('Requeue=0','Requeue=1'),('OverSubscribe=NO','OverSubscribe=YES'),('mem=48000M','mem=16000M'),('MinMemoryCPU=4000M','MinMemoryCPU=8000M')]:
        reject('reject '+before,lambda before=before,after=after:M.parse_allocation(raw.replace(before,after),env,'node-64-768g-4',601))
    reject('reject login host',lambda:M.parse_allocation(raw,env,'biotite',601))
    reject('reject wrong user',lambda:M.parse_allocation(raw,env,'node-64-768g-4',602))
    reject('reject absent allocation environment',lambda:M.parse_allocation(raw,{},'node-64-768g-4',601))
    parent=HERE.parent/'writable_indices_recovery_v1'
    def method(p,name):return next(x for x in ast.walk(ast.parse(p.read_text())) if isinstance(x,ast.FunctionDef) and x.name==name)
    for name in ('nominal_parameters','nominal_provenance'):
        check('unchanged '+name,ast.dump(method(HERE/'RUN.py',name))==ast.dump(method(parent/'RUN.py',name)))
    new=copy.deepcopy(method(HERE/'RUN.py','worker'))
    removed=0
    for call in ast.walk(new):
        if isinstance(call,ast.Call):
            old=len(call.keywords);call.keywords=[k for k in call.keywords if not(k.arg=='allocation' and ast.unparse(k.value)=='ACTUAL_ALLOCATION')];removed+=old-len(call.keywords)
    check('worker only allocation-metadata delta',removed==1 and ast.dump(new)==ast.dump(method(parent/'RUN.py','worker')))
    cfg=json.loads((HERE/'CONFIG.json').read_text());old=json.loads((parent/'CONFIG_UNRELEASED.json').read_text())
    extra={'scheduler_command','parent_ready','parent_implementation','writable_index_controls','scheduling_declaration'}
    check('all scientific config fields unchanged',{k:v for k,v in cfg.items() if k not in extra|{'implementation','resources'}}=={k:v for k,v in old.items() if k not in {'implementation','resources'}})
    for k in ('address_space_bytes','wall_seconds','CPU_seconds','openblas_threads','openmp_threads'):check('unchanged cap '+k,cfg['resources'][k]==old['resources'][k])
    script=(HERE/'RUN.sbatch').read_text()
    check('default scheduler memory/time and shared node',all(x not in script for x in ('--mem','--time','--exclusive','--oversubscribe')))
    check('fixed scheduler profile',all(x in script for x in ('#SBATCH --nodes=1','#SBATCH --ntasks=1','#SBATCH --cpus-per-task=12','#SBATCH --qos=standard_expand','#SBATCH --partition=standard','#SBATCH --account=banfield','#SBATCH --no-requeue')))
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    result=dict(status='passed_slurm_static_factor_source_policy_controls',implementation=M.pin(HERE/'RUN.py'),configuration=M.pin(HERE/'CONFIG.json'),wrapper=M.pin(HERE/'RUN.sbatch'),test_source=M.pin(__file__),checks=checks,checks_passed=len(checks),scheduler_calls=0,submissions=0,full_matrix_reads=0,numeric_factorizations=0,solve_calls=0,native_calls=0,cost=dict(wall_seconds=time.perf_counter()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1))
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(result=M.pin(out/'RESULT.json'),checks=len(checks),cost=result['cost'])))

if __name__=='__main__':
    import sys
    run(sys.argv[1])
