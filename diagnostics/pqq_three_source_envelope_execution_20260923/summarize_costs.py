"""Actual integration receipts/accounting; nested stage times are not added twice."""
import argparse
from datetime import datetime
import json
from pathlib import Path
import statistics
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from affordable_common import read_json,record,verify,write_new


def seconds(value):
    days=0
    if '-' in value:d,value=value.split('-',1);days=int(d)
    return days*86400+sum(float(v)*60**i for i,v in enumerate(reversed(value.split(':'))))


def summarize(plan,accounting,output):
    root=Path(plan).parent;job=read_json(root/'SUBMISSION.json')['job_id']
    table={v[0]:v for line in Path(accounting).read_text().splitlines() if (v:=line.split('|'))}
    row=table[job];steps={k:v for k,v in table.items() if k.startswith(job+'.')}
    scalar=[];native=[];search=[];workers=[]
    for path in root.rglob('endpoint.out.execution.json'):
        r=read_json(path);a,b=r.get('started_at_utc'),r.get('finished_at_utc')
        wall=(datetime.fromisoformat(b.replace('Z','+00:00'))-datetime.fromisoformat(a.replace('Z','+00:00'))).total_seconds() if a and b else None
        scalar.append({'receipt':record(path),'returncode':r['returncode'],'normal_termination':r['normal_termination'],
            'SCF_converged':r['scf_converged'],'wall_seconds':wall,'parallelism':r['parallelism']})
    for path in root.rglob('mace_result.json'):
        r=read_json(path);native.append({'receipt':record(path),'stage':path.relative_to(root).parts[0],
            'status':r['status'],'wall_seconds':r['wall_seconds'],'model_call_started':r['model_call_started']})
    for path in root.glob('*/gpu_summary_*.json'):workers.append({'receipt':record(path),**read_json(path)})
    for path in (root/'proposals/proposals').glob('*/result.json'):
        r=read_json(path);search.append({'receipt':record(path),'task_id':r['task_id'],'status':r['status'],
            'reason':r.get('reason'),'wall_seconds':r['wall_seconds'],'optimizer':r['optimizer'],'boundary_flag':r.get('boundary_flag')})
    request=read_json(verify(read_json(plan)['request']));preps=[read_json(verify(g['preparation'])) for g in request['groups']]
    result={'plan':record(plan),'submission':record(root/'SUBMISSION.json'),'accounting':record(accounting),
        'job_id':job,'state':row[1],'exit_code':row[2],'allocated_CPUs':int(row[3]),
        'allocation_wall_seconds':int(row[4]),'allocated_core_seconds':int(row[5]),'requested_GPU_seconds':int(row[4]),
        'sum_distinct_step_reported_CPU_seconds':sum(seconds(v[6]) for v in steps.values()),
        'step_peak_RSS':{k:v[9] for k,v in steps.items()},'requested_memory':row[7],
        'separate_preparation_function_wall_seconds':sum(p['preparation_wall_seconds'] for p in preps),
        'historical_protonation_reused':True,'preparation_timing_excludes_imports_request_validation_final_dry_run':True,
        'GFN2_attempts':len(scalar),'GFN2_normal':sum(r['normal_termination'] for r in scalar),
        'GFN2_summed_attempt_wall_seconds':sum(r['wall_seconds'] or 0 for r in scalar),
        'GFN2_median_attempt_wall_seconds':statistics.median(r['wall_seconds'] for r in scalar if r['wall_seconds'] is not None) if scalar else None,
        'GFN2_receipts':scalar,'native_MACE_requests':len(native),'native_MACE_complete':sum(r['status']=='complete' for r in native),
        'native_MACE':native,'warm_workers':workers,'searches':search,
        'executor_receipts':[record(p) for p in root.glob('execution_*.json')],
        'stage_timings_nested_not_additional_allocations':True,'local_setup_tests_reporting_additional_unmetered':True,
        'new_DFT':0,'archive_energy_force_reuse':False,'scientific_retry':False}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('GFN2_receipts','native_MACE','warm_workers','searches')}


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for name in ('plan','accounting','output'):ap.add_argument('--'+name,required=True)
    print(json.dumps(summarize(**vars(ap.parse_args())),indent=2))
