"""Tiny artificial reviewer controls; no actual execution/result is read."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time


def main(output):
    tick,cpu=time.monotonic(),time.process_time();out=Path(output);out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).with_name('REVIEW.py')
    spec=importlib.util.spec_from_file_location('saved_review',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    checks=[]
    def check(name,ok):
        checks.append(dict(label=name,passed=bool(ok)))
        if not ok:raise AssertionError(name)
    def rejects(name,call):
        try:call()
        except (ValueError,KeyError):check(name,True)
        else:check(name,False)
    try:
        tree=ast.parse(source.read_text())
        imports=[n.names[0].name for n in ast.walk(tree) if isinstance(n,ast.Import)]
        check('No scientific/native/producer import',not set(imports)&{'numpy','scipy','subprocess','ctypes','RUN','REPLAY'})
        event=dict(jobs=['20135'],records={'20135':dict(job_id='20135',state='COMPLETED',elapsed_seconds=104,
                   allocated_cpus=1,end='end',start='start',exit_code='0:0',nodes='node-example')})
        account=m.scheduler(event)
        check('Single allocation cost',account['allocated_CPU_seconds']==104)
        r=dict(returncode=0,stderr='',job_id='20135',allocated_CPU_seconds=104,
               stdout='20135|COMPLETED|0:0|104|1|01:40.251||node-example\n20135.batch|COMPLETED|0:0|104|1|01:40.251|947952K|node-example\n')
        cost=m.accounting(r,account)
        check('Batch not double counted',cost['allocated_CPU_seconds']==104 and cost['batch_cost_added'] is False)
        check('Measured CPU kept distinct',cost['measured_total_CPU_seconds']==100.251 and cost['batch_MaxRSS_KiB']==947952)
        for name,key,value in [('Wrong job','jobs',['20136']),('Nonterminal','records',{'20135':dict(event['records']['20135'],state='RUNNING')}),
                               ('Wrong allocation','records',{'20135':dict(event['records']['20135'],allocated_cpus=2)})]:
            mutated=copy.deepcopy(event);mutated[key]=value;rejects(name,lambda:m.scheduler(mutated))
        failed=copy.deepcopy(event);failed['records']['20135']['state']='FAILED';failed['records']['20135']['exit_code']='1:0'
        check('Terminal failure retained',m.scheduler(failed)['state']=='FAILED')
        for name,key,value in [('Accounting failed','returncode',1),('Double counted cost','allocated_CPU_seconds',208),
                               ('Extra accounting row','stdout',r['stdout']+'20136|COMPLETED|0:0|1|1|00:01||node-example\n'),
                               ('Different elapsed','stdout',r['stdout'].replace('|104|','|105|'))]:
            mutated=dict(r);mutated[key]=value;rejects(name,lambda:m.accounting(mutated,account))
        def row(label):return dict(label=label,max_abs_error=1e-15,scale=1.,bound=1e-11,passed=True)
        replay=dict(checks=[row(f'pdamp {i}') for i in range(1,11959)]+[row(x) for x in
                    ('native compact h','native compact H','normalized H versus job20053')],applications=[])
        labels=[f'artificial_{i}' for i in range(10)]
        for mode in (1,2,3):
            for probe in range(1,11):
                own=[row(f'{mode}/{probe}/{channel}/{kind}') for channel in ('real','imag') for kind in
                     ('reciprocal','complete','real+self+H')]
                replay['checks'].extend(own);replay['applications'].append(dict(mode=mode,probe=probe,label=labels[probe-1],checks=own))
        passed=m.parity_checks(replay,labels)
        check('Exact manufactured finite denominator',passed['check_count']==12141 and passed['applications']==30)
        def altered(name,modify):
            changed=copy.deepcopy(replay);modify(changed);rejects(name,lambda:m.parity_checks(changed,labels))
        altered('Missing check',lambda x:x['checks'].pop())
        altered('Duplicate check',lambda x:x['checks'].append(copy.deepcopy(x['checks'][0])))
        altered('Relaxed bound rejected',lambda x:x['checks'][0].update(bound=1e-10))
        altered('Nonfinite error rejected',lambda x:x['checks'][0].update(max_abs_error=float('nan')))
        altered('Failed gate rejected',lambda x:x['checks'][0].update(max_abs_error=2e-11))
        altered('False pass rejected',lambda x:x['checks'][0].update(passed=False))
        altered('Channel/probe identity rejected',lambda x:x['applications'][0].update(probe=2))
        altered('Duplicated application metadata must join',lambda x:x['applications'][0].update(checks=[]))
        status='passed_stage1_saved_reviewer_controls';error=None
    except Exception as exc:status='failed_stage1_saved_reviewer_controls';error=repr(exc)
    result=dict(status=status,error=error,implementation=m.pin(__file__),reviewer=m.pin(source),checks=checks,
                check_count=len(checks),actual_execution_results_read=0,new_native_calls=0,new_builds=0,new_FFT_or_sparse_replays=0,
                cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result))
    return 0 if error is None else 1


if __name__=='__main__':sys.exit(main(sys.argv[1]))
