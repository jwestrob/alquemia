"""Tiny source/command controls. No compiler, application, scheduler or loader call."""
import ast
import copy
import importlib.util
import json
import os
from pathlib import Path
import resource
import tempfile
import time
import BUILD as B

def run(inputs,out):
    out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'failed','implementation':B.pin(__file__),'inputs':B.pin(inputs),'checks':[],
      'compiler_calls':0,'native_calls':0,'scheduler_calls':0,'error':None}
    def check(name,ok):B.require(ok,name);result['checks'].append(name)
    def rejects(name,fn):
        try:fn()
        except (ValueError,AssertionError):check(name,True)
        else:raise AssertionError(name)
    try:
        x=json.loads(inputs.read_text());h=B.read(x['handoff'])
        for p in x['sources'].values():B.verify(p)
        check('all fixed source/helper pins',True)
        for p in h['source_inventory']['radial_archive_files']+list(h['source_inventory']['isolated_objects'].values()):B.verify(p)
        check('615 archive files and eight isolated source pins',len(h['source_inventory']['radial_archive_files'])==615 and len(h['source_inventory']['isolated_objects'])==8)
        for p in x['sources'].values():
            if p['path'].endswith('.py'):ast.parse(Path(p['path']).read_text())
        check('source AST parse',True)
        plan=B.commands(h,out/'manufactured',x['tools']['cmake']['path'])
        check('13 unique fixed commands',len(plan)==13 and len({p['label'] for p in plan})==13)
        check('12 archive workers only',plan[1]['argv'][-3:]==['tinker','--parallel','12'])
        check('tinker-only archive target',plan[1]['argv'][3:5]==['--target','tinker'])
        check('OpenMM disabled',plan[0]['argv'][-1]=='-DTINKER_OPENMM=OFF')
        check('eight sequential compile slots',sum(p['kind']=='compile' for p in plan)==8)
        check('three link slots',sum(p['kind']=='link' for p in plan)==3)
        check('no application or old-baseline rebuild',all(p['argv'][0] in (x['tools']['cmake']['path'],'/usr/bin/f95') for p in plan) and all('gradient_old' not in ' '.join(p['argv']) for p in plan))
        oldmod=str(Path(h['old_baseline_reuse']['polar_module']['path']).parent)
        check('old module directory excluded',all(oldmod not in ' '.join(p['argv']) for p in plan))
        links={p['role']:p['argv'] for p in plan if p['kind']=='link'}
        for role in ('initialize','gradient_clone'):
            check(role+' complete induce object precedes archive',next(i for i,v in enumerate(links[role]) if v.endswith('/induce_qualification.f.o'))<next(i for i,v in enumerate(links[role]) if v.endswith('/libtinker.a')))
        check('parity excludes qualification override',not any(v.endswith('/induce_qualification.f.o') for v in links['parity']))
        for role,a in links.items():check(role+' shared V2 writer linked',sum(v.endswith('/STATE_WRITER_V2.f90.o') for v in a)==1)
        check('no new architecture flags',all(not any(v.startswith(('-march','-mtune')) for v in p['argv']) for p in plan))
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);calls=[];persisted=[]
            def fake(result,out,label,argv,cwd,wall,memory,cpu):
                calls.append(label);result['commands'][label]={'status':'manufactured_only'}
                if '-o' in argv:
                    output=Path(argv[argv.index('-o')+1]);output.parent.mkdir(parents=True,exist_ok=True);output.write_text('manufactured, not an object or executable\n')
            r={'attempts':dict.fromkeys(B.MAXIMA,0),'commands':{}}
            B.execute_plan(B.commands(h,tmp/'success',x['tools']['cmake']['path']),r,tmp,fake,lambda:persisted.append(copy.deepcopy(r)))
            check('manufactured full command denominator',r['attempts']==B.MAXIMA and len(calls)==13)
            check('attempt persisted before each command',len(persisted)==13 and all(sum(v['attempts'].values())==i+1 for i,v in enumerate(persisted)))
            calls.clear();r={'attempts':dict.fromkeys(B.MAXIMA,0),'commands':{}}
            def fail(result,out,label,*args):
                calls.append(label)
                if label=='archive':raise ValueError('manufactured archive failure')
            rejects('archive failure stops before wrappers',lambda:B.execute_plan(B.commands(h,tmp/'failure',x['tools']['cmake']['path']),r,tmp,fail,lambda:None))
            check('no retry after archive failure',calls==['configure','archive'] and r['attempts']['compile']==r['attempts']['link']==0)
            r={'attempts':dict.fromkeys(B.MAXIMA,0),'commands':{}}
            rejects('missing object is terminal',lambda:B.execute_plan(B.commands(h,tmp/'missing',x['tools']['cmake']['path']),r,tmp,lambda *a:None,lambda:None))
            check('missing object stops first compile',r['attempts']=={'configure':1,'archive':1,'compile':1,'link':0})
            r={'attempts':{'configure':1,'archive':0,'compile':0,'link':0},'commands':{}}
            rejects('attempt ceiling cannot extend',lambda:B.execute_plan(B.commands(h,tmp/'cap',x['tools']['cmake']['path']),r,tmp,lambda *a:None,lambda:None))
        allocation=B.load(x['sources']['allocation12'],'mock_allocation')
        env={'SLURM_JOB_ID':'123','SLURM_CPUS_PER_TASK':'12','SLURM_NTASKS':'1','SLURM_JOB_NUM_NODES':'1'}
        raw='JobId=123 JobState=RUNNING Account=banfield QOS=standard_expand Partition=standard NumNodes=1 NumCPUs=12 NumTasks=1 CPUs/Task=12 Requeue=0 OverSubscribe=NO BatchHost=mock UserId=mock(1234) AllocTRES=cpu=12,mem=48000M,node=1 MinMemoryCPU=4000M'
        rec=allocation.parse_allocation(raw,env,'mock',1234)
        check('actual helper accepts manufactured default allocation',rec['allocated_CPUs']==12 and rec['memory_grant_bytes']==48000*1024**2)
        rejects('wrong allocation rejected',lambda:allocation.parse_allocation(raw.replace('NumCPUs=12','NumCPUs=1'),env,'mock',1234))
        rejects('memory override rejected',lambda:allocation.parse_allocation(raw.replace('48000M','50000M'),env,'mock',1234))
        wrapper=(Path(__file__).parent/'RUN.sbatch').read_text()
        check('no memory/time override',not any(x in wrapper for x in ('#SBATCH --mem','#SBATCH --time')))
        check('single12CPU no-requeue scheduler contract',all(x in wrapper for x in ('--partition=standard','--qos=standard_expand','--cpus-per-task=12','--ntasks=1','--no-requeue')))
        result['status']='passed_periodic_v2_build_source_controls'
    except BaseException as exc:result['error']=repr(exc)
    result['check_count']=len(result['checks']);result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    B.save(out/'RESULT.json',result);return result

if __name__=='__main__':
    import sys
    result=run(Path(sys.argv[1]),Path(sys.argv[2]));print(json.dumps({k:result[k] for k in ('status','check_count','error','cost')}));raise SystemExit(0 if result['error'] is None else 1)
