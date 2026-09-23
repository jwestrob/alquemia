"""Actual matched-job accounting; stage timings overlap and are not additive costs."""
import argparse
import datetime as dt
import json
from pathlib import Path
import statistics
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from affordable_common import read_json,record,write_new


def seconds(value):
    days=0
    if '-' in value:
        d,value=value.split('-',1);days=int(d)
    return days*86400+sum(float(v)*60**i for i,v in enumerate(reversed(value.split(':'))))


def duration(receipt):
    a,b=receipt.get('started_at_utc'),receipt.get('finished_at_utc')
    return (dt.datetime.fromisoformat(b.replace('Z','+00:00'))-dt.datetime.fromisoformat(a.replace('Z','+00:00'))).total_seconds() if a and b else None


def summarize(plan,accounting):
    p=read_json(plan);root=Path(plan).parent;submission=read_json(root/'SUBMISSION.json')
    job=submission['stdout'].strip().split(';')[0]
    table={v[0]:v for line in Path(accounting).read_text().splitlines() if (v:=line.split('|'))}
    row=table[job];steps={k:v for k,v in table.items() if k.startswith(job+'.')}
    prep=root/'source_preparation/preparation.json';prepared=read_json(prep) if prep.exists() else None
    source_rows=(read_json(prepared['source_preparation']['path'])['cases']
                 if prepared and p['arm']=='union-candidate' else prepared['cases'] if prepared else [])
    scalar=[]
    for path in root.rglob('endpoint.out.execution.json'):
        receipt=read_json(path);scalar.append({'receipt':record(path),'case_id':receipt.get('task_id'),
            'returncode':receipt['returncode'],'normal_termination':receipt['normal_termination'],
            'SCF_converged':receipt.get('scf_converged'),'wall_seconds':duration(receipt)})
    native=[];workers=[];search=[];static_timings=[]
    if p['arm']=='released-static':
        for path in (root/'prepared_score/scoring/mace').glob('*/result.json'):
            value=read_json(path);native.append({'receipt':record(path),'status':value['status'],
                'task_id':value['task_id'],'wall_seconds':value['wall_seconds'],
                'initialization_and_input_seconds':value.get('model_load_seconds'),
                'evaluation_seconds':value.get('evaluation_seconds')})
        for path in (root/'prepared_score/scoring/timing').glob('*.json'):
            static_timings.append({'receipt':record(path),**read_json(path)})
    else:
        for path in root.rglob('mace_result.json'):
            value=read_json(path);native.append({'receipt':record(path),'status':value['status'],
                'stage':path.relative_to(root).parts[0],'model_call_started':value['model_call_started'],
                'wall_seconds':value['wall_seconds']})
        for path in root.glob('*/gpu_summary_*.json'):
            workers.append({'receipt':record(path),**read_json(path)})
        for path in (root/'proposals/proposals').glob('*/result.json'):
            value=read_json(path);search.append({'receipt':record(path),'task_id':value['task_id'],
                'status':value['status'],'reason':value.get('reason'),'wall_seconds':value['wall_seconds'],
                'optimizer':value['optimizer'],'boundary_flag':value.get('boundary_flag')})
    times=[v['wall_seconds'] for v in scalar if v['wall_seconds'] is not None]
    solver_events=[]
    for path in root.rglob('budget_events.jsonl'):
        events=[json.loads(line) for line in path.read_text().splitlines()]
        solver_events.extend({'events_file':record(path),**v} for v in events if 'allocated_core_seconds' in v)
    return {'plan':record(plan),'submission':record(root/'SUBMISSION.json'),'accounting':record(accounting),
        'job_id':job,'arm':p['arm'],'state':row[1],'exit_code':row[2],'allocated_CPUs':int(row[3]),
        'allocation_wall_seconds':int(row[4]),'allocated_core_seconds':int(row[5]),
        'requested_GPU_seconds':int(row[4]),'sum_distinct_step_reported_CPU_seconds':sum(seconds(v[6]) for v in steps.values()),
        'step_reported_CPU_seconds':{k:seconds(v[6]) for k,v in steps.items()},
        'step_peak_RSS':{k:v[9] for k,v in steps.items()},'requested_memory':row[7],
        'fresh_preparation_seconds':prepared['fresh_preparation_seconds'] if prepared else None,
        'per_source_preparation':[{'case_id':c['case_id'],'status':c['status'],'source_preparation_seconds':c.get('source_preparation_seconds'),
                                  'normalization_seconds':c.get('normalization_seconds'),'protonation_seconds':c.get('protonation_seconds')}
                                 for c in source_rows],
        'GFN2_attempts':len(scalar),'GFN2_normal':sum(x['normal_termination'] for x in scalar),
        'GFN2_summed_attempt_wall_seconds':sum(times),'GFN2_median_attempt_wall_seconds':statistics.median(times) if times else None,
        'GFN2_receipts':scalar,'solver_execution_events':solver_events,'native_MACE_requests':len(native),'native_MACE_complete':sum(x['status'] in ('computed','complete') for x in native),
        'native_MACE':native,'warm_workers':workers,'searches':search,'static_per_source_score_timings':static_timings,
        'stage_timings_are_nested_or_overlapping_not_additional_allocations':True,
        'same_host_two_user_paths_not_isolated_MPI_or_accommodation_speedup':True,
        'new_DFT':0,'local_setup_tests_and_reporting_additional_unmetered':True}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plans',nargs=2,type=Path,required=True);ap.add_argument('--accounting',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    rows=[summarize(p,a.accounting) for p in a.plans]
    result={'arms':rows,'allocated_core_seconds':sum(r['allocated_core_seconds'] for r in rows),
        'requested_GPU_seconds':sum(r['requested_GPU_seconds'] for r in rows),'new_GFN2_attempts':sum(r['GFN2_attempts'] for r in rows),
        'new_MACE_requests':sum(r['native_MACE_requests'] for r in rows),'new_DFT':0,'archive_energy_force_reuse':False}
    write_new(a.output,result);print(json.dumps({k:v for k,v in result.items() if k!='arms'}))

if __name__=='__main__':main()
