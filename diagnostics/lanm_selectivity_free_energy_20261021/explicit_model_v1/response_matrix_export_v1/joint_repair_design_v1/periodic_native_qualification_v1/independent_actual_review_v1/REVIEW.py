"""Terminal saved-artifact audit; no producer import, native or numerical replay."""
import hashlib
import json
import math
from pathlib import Path
import time
from decimal import Decimal


def pin(path):
    path=Path(path).resolve()
    with path.open('rb') as stream:sha=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(path),sha256=sha)


def verify(item):
    if not isinstance(item,dict) or set(item)!= {'path','sha256'} or pin(item['path'])!=item:
        raise ValueError('Missing/changed pin')
    return Path(item['path'])


def read(item):return json.loads(verify(item).read_text())


def scheduler(event):
    if event['jobs']!=['20135']:raise ValueError('Exact released job20135')
    row=event['records']['20135'];state=row['state'].split()[0]
    if row['job_id']!='20135' or state not in ('COMPLETED','FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED','PREEMPTED','NODE_FAIL'):
        raise ValueError('Actual terminal scheduler record required')
    elapsed=int(row['elapsed_seconds']);cpus=int(row['allocated_cpus'])
    if cpus!=1 or elapsed<0 or not row['end']:raise ValueError('Scheduler resource/end fields')
    return dict(job_id='20135',state=state,exit_code=row['exit_code'],elapsed_seconds=elapsed,allocated_CPUs=cpus,
                allocated_CPU_seconds=elapsed*cpus,nodes=row['nodes'],start=row['start'],end=row['end'])


def accounting(receipt,account):
    if receipt['returncode']!=0 or receipt['stderr'].strip() or receipt['job_id']!='20135':
        raise ValueError('Successful saved scheduler accounting required')
    rows=[line.split('|') for line in receipt['stdout'].splitlines() if line]
    if [row[0] for row in rows]!=['20135','20135.batch'] or any(len(row)!=8 for row in rows):
        raise ValueError('Exact job and batch accounting rows')
    job,batch=rows
    if (job[1],job[2],int(job[3]),int(job[4]),job[7])!=(account['state'],account['exit_code'],
            account['elapsed_seconds'],account['allocated_CPUs'],account['nodes']):
        raise ValueError('Completion and scheduler accounting disagree')
    if receipt['allocated_CPU_seconds']!=int(job[3])*int(job[4]):raise ValueError('Allocation cost mismatch')
    def seconds(value):
        days,clock=value.split('-',1) if '-' in value else ('0',value)
        fields=list(map(Decimal,clock.split(':')))
        if len(fields) not in (2,3) or any(not x.is_finite() or x<0 for x in fields):raise ValueError('CPU time format')
        total=Decimal(days)*86400
        for i,x in enumerate(reversed(fields)):total+=x*(60**i)
        return float(total)
    rss=batch[6]
    if not rss.endswith('K') or not rss[:-1].isdecimal():raise ValueError('Saved batch RSS format')
    return dict(allocated_CPU_seconds=int(job[3])*int(job[4]),measured_total_CPU_seconds=seconds(job[5]),
                batch_MaxRSS_KiB=int(rss[:-1]),raw_job=job,raw_batch=batch,
                batch_cost_added=False,memory_scope='Slurm batch MaxRSS; no independent whole-job instantaneous peak inferred.')


def parity_checks(replay,labels):
    expected={f'pdamp {i}' for i in range(1,11959)}|{'native compact h','native compact H','normalized H versus job20053'}
    for mode in (1,2,3):
        for probe in range(1,11):
            for channel in ('real','imag'):
                for kind in ('reciprocal','complete','real+self+H'):
                    expected.add(f'{mode}/{probe}/{channel}/{kind}')
    checks=replay['checks'];observed={x['label'] for x in checks}
    if observed!=expected or len(checks)!=len(expected):raise ValueError('Exact persisted parity-check denominator')
    by_label={x['label']:x for x in checks}
    for row in checks:
        if not all(isinstance(row[k],(int,float)) and math.isfinite(row[k]) for k in ('max_abs_error','scale','bound')):
            raise ValueError('Finite recorded parity arithmetic')
        if row['scale']<1 or row['max_abs_error']<0 or row['bound']!=1e-11*row['scale']:
            raise ValueError('Unchanged declared parity bound')
        if row['passed'] is not True or row['max_abs_error']>row['bound']:raise ValueError('Recorded parity gate failed')
    applications=replay['applications']
    if len(applications)!=30 or len(labels)!=10:raise ValueError('Thirty applications')
    for index,row in enumerate(applications):
        mode,probe=index//10+1,index%10+1
        if (row['mode'],row['probe'],row['label'])!=(mode,probe,labels[probe-1]):raise ValueError('Probe/channel order')
        names=[f'{mode}/{probe}/{channel}/{kind}' for channel in ('real','imag')
               for kind in ('reciprocal','complete','real+self+H')]
        if row['checks']!=[by_label[x] for x in names]:raise ValueError('Application/global arithmetic joins')
    return dict(check_count=len(checks),applications=30,maximum_scaled_error=max(x['max_abs_error']/x['scale'] for x in checks),
                gate=1e-11,full_FFT_or_sparse_replay_repeated=False)


def review(config_path,checks):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config_path).read_text())
    def check(label,condition):
        checks.append(dict(label=label,passed=bool(condition)))
        if not condition:raise ValueError(label)
    check('Bound terminal reviewer configuration',cfg['status']=='bound_stage1_terminal_review_configuration')
    check('Exact reviewer source',cfg['implementation']==pin(__file__))
    event=read(cfg['completion']);account=scheduler(event)
    cost_receipt=read(cfg['scheduler_accounting']);account['saved_accounting']=accounting(cost_receipt,account)
    result=read(cfg['result']);root=read(cfg['root_release']);runtime=read(cfg['runtime_configuration'])
    ready=read(cfg['stage1_ready']);prep=read(ready['preparation']);submission=read(cfg['submission'])
    source_review=read(runtime['stage1_source_review']);adapter_review=read(root['adapter_source_review'])
    check('Exact submission release/job',submission['job_id']=='20135' and submission['returncode']==0 and submission['release']==cfg['root_release'])
    check('Submitted argv authority',submission['argv'][-4:]==[runtime['wrapper']['path'],cfg['runtime_configuration']['path'],cfg['root_release']['path'],root['output']])
    check('Exact released configuration and implementation',root['configuration']==cfg['runtime_configuration'] and
          root['implementation']==runtime['implementation']==result['implementation'] and result['configuration']==cfg['runtime_configuration'])
    check('Exact actual release pin',result['release']==cfg['root_release'])
    check('Actual result location',Path(cfg['result']['path']).parent.resolve()==Path(root['output']).resolve())
    check('Root finite declaration',root['status']=='released_one_slurm_periodic_operator_parity' and root['maximum_native_invocations']==1
          and root['maximum_links']==1 and root['maximum_saved_replays']==1)
    check('Root/source READY join',root['stage1_ready']==runtime['stage1_ready']==cfg['stage1_ready'])
    check('Unchanged allocation policy',root['resources']==runtime['resources'])
    check('Stage1 independent source review',source_review['status']=='passed_stage1_complete_periodic_operator_source_review'
          and source_review['ready']==cfg['stage1_ready'] and source_review['implementation']==ready['implementation']
          and source_review['replay']==ready['replay'] and source_review['preparation']==ready['preparation'])
    check('Launcher source review',adapter_review['status']=='passed_stage1_slurm_adapter_source_review'
          and adapter_review['implementation']==runtime['implementation'] and adapter_review['configuration']==cfg['runtime_configuration'])
    source=read(ready['manifest'])
    for item in source['sources'].values():verify(item)
    verify(runtime['implementation']);verify(runtime['preflight']);verify(runtime['wrapper'])
    static=read(ready['static_input']);admission=read(prep['static_root_admission'])
    check('Exact old e41f and static20053 identity',root['candidate']==ready['candidate']==prep['candidate']==static['candidate'] and
          root['static_input']==ready['static_input']==prep['static_input']==admission['input_manifest'])
    check('Existing nominal static admission reused',admission['status']=='admitted_nominal_saved_La_static_numerical_test'
          and admission['job_id']==20053 and admission['sufficient_test_passed'] is True)
    check('Static admission source join',runtime['static_root_admission']==prep['static_root_admission'])
    check('No extra science counters',all(result.get(k)==0 for k in ('induced_solves','energies','gradients','factorizations')))
    for k in ('link_attempts','links_launched','native_attempts','native_invocations','replay_attempts'):
        check('Finite count '+k,type(result[k]) is int and result[k] in (0,1))
    check('Attempt/invocation ordering',result['links_launched']<=result['link_attempts'] and result['native_invocations']<=result['native_attempts'])
    commands=result['commands'];order=['preflight','link','loader','native','replay']
    check('Sequential command prefix',list(commands)==order[:len(commands)])
    for name,row in commands.items():
        check('Retained command status '+name,row['status'] in ('intent_durable_before_launch','launched','complete','failed'))
        if row.get('stdout'):verify(row['stdout'])
        if row.get('stderr'):verify(row['stderr'])
        if row.get('PID') is not None:check('Launched PID '+name,type(row['PID']) is int and row['PID']>0)
        check('Declared memory cap '+name,row['address_space_limit_bytes']==3*1024**3)
    check('Native launch count from receipt',result['native_invocations']==int(commands.get('native',{}).get('PID') is not None))
    check('Link launch count from receipt',result['links_launched']==int(commands.get('link',{}).get('PID') is not None))
    if 'allocation' in result:
        allocation=result['allocation']
        check('Actual node/job/resources',allocation['job_id']=='20135' and allocation['allocated_CPUs']==1 and
              allocation['memory_grant_bytes']==4000*1024**2 and allocation['qos']=='standard_expand' and
              allocation['hostname'].split('.')[0] in account['nodes'].split(','))
    for category in ('raw_inputs','raw_outputs'):
        for item in result.get(category,{}).values():verify(item)
    if result.get('executable'):verify(result['executable'])
    complete=result['status']=='stage1_periodic_operator_parity_complete'
    check('Producer terminal or scheduler-proven incomplete attempt',complete or result['status'] in
          ('stage1_periodic_operator_parity_failed','started'))
    scientific={};replay_result=None
    if result.get('replay_result'):
        replay_result=read(result['replay_result'])
    if complete:
        check('Successful scheduler termination',account['state']=='COMPLETED' and account['exit_code']=='0:0')
        check('Complete producer has no error',result['error'] is None and result['sampled_operator_parity_passed'] is True)
        check('All single launches completed',all(result[k]==1 for k in
              ('link_attempts','links_launched','native_attempts','native_invocations','replay_attempts')) and
              list(commands)==order and all(row['status']=='complete' and row['returncode']==0 for row in commands.values()))
        recipe=read(ready['link_recipe'])
        check('Unchanged linked archive',recipe['library']==ready['library']==result['library'])
        verify(recipe['library'])
        check('Exact compile/link command',commands['link']['argv']==recipe['argv'])
        check('Exact native command',commands['native']['argv']==[recipe['output'],'system.xyz'])
        check('Native executable identity',result['executable']==pin(recipe['output']))
        check('Exact saved replay command',commands['replay']['argv'][1:]==[ready['replay']['path'],ready['preparation']['path'],
              str(Path(root['output'])/'native'),str(Path(root['output'])/'replay')])
        check('Node native runtime exact',result['actual_native_runtime']==runtime['native_runtime'])
        for item in result['actual_native_runtime'].values():verify(item)
        node=read(result['node_preflight']);environment=read(runtime['environment'])
        check('Node Python/thread/memory preflight',node['status']=='passed_stage1_node_runtime_preflight'
              and node['prefix']==environment['prefix'] and node['actual_threads']=={'openblas':1,'openmp':1}
              and node['address_space_limit']==[3*1024**3]*2)
        check('Node numerical versions',all(node['versions'][k]==environment['versions'][k] for k in ('python','numpy','scipy')))
        check('Exact input set',set(result['raw_inputs'])==set(prep['raw_inputs']))
        for name,item in result['raw_inputs'].items():
            check('Exact input '+name,item['sha256']==prep['raw_inputs'][name]['sha256'])
        check('Successful retained saved replay',replay_result is not None and
              replay_result['status']=='passed_finite_complete_periodic_operator_parity')
        check('Replay immutable lineage',replay_result['implementation']==ready['replay'] and replay_result['preparation']==ready['preparation']
              and replay_result['candidate']==ready['candidate'] and replay_result['static_input']==ready['static_input'])
        check('Same initialized PME state as admitted static matrix',replay_result['actual_snapshot']['sha256']==static['pme_export']['sha256']
              and replay_result['actual_snapshot']==result['raw_outputs']['pme_norm_export_001.dat'])
        check('Exact application counters',replay_result['counts']==prep['counts'])
        for name,item in replay_result['raw'].items():check('Replay raw authority '+name,item==result['raw_outputs'][name])
        expected_receipt='\n'.join(['PERIODIC_PARITY_RECEIPT_V1',
            'COUNTS full=30 separate_reciprocal=30 total_reciprocal=60 direct_initializations=1',
            'COUNTS induced_solves=0 energies=0 gradients=0 factorizations=0',
            'NORMAL_GUARDS_UNCHANGED TRUE','STATE_RESTORED TRUE','END',''])
        check('Native terminal counts/guards/state receipt',Path(result['raw_outputs']['parity_receipt.dat']['path']).read_text()==expected_receipt)
        scientific=parity_checks(replay_result,prep['labels'])
        check('Finite ordinary-double scope only',replay_result['e_app'] is None and not replay_result['uniform_application_bound']
              and not replay_result['energy_force_qualified'] and not replay_result['physical_accuracy_qualified'])
    elif replay_result is not None:
        scientific=dict(retained_replay_status=replay_result['status'],retained_error=replay_result.get('error'),
                        completed_check_records=len(replay_result.get('checks',[])),admitted=False)
    return dict(status='passed_stage1_periodic_operator_actual_review' if complete else 'reviewed_failed_stage1_periodic_operator_attempt',
                reviewer_implementation=pin(__file__),review_configuration=pin(config_path),result=cfg['result'],
                implementation=runtime['implementation'],configuration=cfg['runtime_configuration'],root_release=cfg['root_release'],
                completion=cfg['completion'],scheduler_accounting=cfg['scheduler_accounting'],submission=cfg['submission'],stage1_ready=cfg['stage1_ready'],
                candidate=ready['candidate'],static_input=ready['static_input'],static_root_admission=prep['static_root_admission'],
                operator_parity_admitted=complete,producer_status=result['status'],producer_error=result.get('error'),
                native_invocations=result['native_invocations'],link_attempts=result['link_attempts'],commands=commands,
                checks=checks,check_count=len(checks),scientific=scientific,scheduler=account,
                producer_cost=result.get('cost'),e_app=None,energy_force_qualified=False,physical_accuracy_qualified=False,
                limitations=['Existing reviewed numerical replay reused; no repeated FFT/sparse action or factor audit.',
                             'Finite sampled operator correspondence only; common direct+h solve, energies/forces and dynamics unqualified.',
                             'Author-authored execution audit with prior independent source review; no independent QM claim.'],
                new_native_calls=0,new_builds=0,new_sparse_or_FFT_replays=0,new_factorizations=0,
                cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))


if __name__=='__main__':
    import sys
    config,out=sys.argv[1:];out=Path(out);out.mkdir(parents=True,exist_ok=False);checks=[]
    try:result=review(config,checks)
    except Exception as error:
        result=dict(status='failed_stage1_artifact_review',error=repr(error),reviewer_implementation=pin(__file__),
                    configuration=pin(config),checks=checks,new_native_calls=0,new_sparse_or_FFT_replays=0)
    (out/'REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'REPORT.txt').write_text(result['status']+'\n'+json.dumps({k:result.get(k) for k in
        ('producer_status','producer_error','operator_parity_admitted','native_invocations','scientific','scheduler','limitations')},indent=2)+'\n')
    print(result['status']);sys.exit(1 if result['status']=='failed_stage1_artifact_review' else 0)
