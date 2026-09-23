"""Read actual rank job/step accounting and every attempted task receipt."""
import argparse
import datetime as dt
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new


def seconds(text):
    days=0
    if '-' in text:day,text=text.split('-',1);days=int(day)
    fields=[float(x) for x in text.split(':')]
    return days*86400+sum(x*60**i for i,x in enumerate(reversed(fields)))


def main():
    p=argparse.ArgumentParser();p.add_argument('--design',type=Path,required=True);p.add_argument('--accounting',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    design=read_json(a.design);root=a.design.parent
    accounting={v[0]:v for line in a.accounting.read_text().splitlines() if (v:=line.split('|'))}
    ranks=[]
    for pin in design['manifests']:
        mp=verify(pin);m=read_json(mp);n=m['ranks'];submitted=read_json(root/f'SUBMISSION_rank_{n}.json');job=submitted['stdout'].strip().split(';')[0];rec=accounting[job]
        steps={key:row for key,row in accounting.items() if key.startswith(job+'.')}
        receipts=[read_json(Path(t['output_path']+'.execution.json')) for t in m['tasks'] if Path(t['output_path']+'.execution.json').exists()]
        wall=[]
        for r in receipts:
            parse=lambda x:dt.datetime.fromisoformat(x.replace('Z','+00:00'))
            wall.append((parse(r['finished_at_utc'])-parse(r['started_at_utc'])).total_seconds())
        events=[json.loads(line) for line in (mp.parent/'budget_events.jsonl').read_text().splitlines()] if (mp.parent/'budget_events.jsonl').exists() else []
        done=next((e for e in reversed(events) if 'allocated_core_seconds' in e),None)
        ranks.append({'ranks':n,'job_id':job,'job_state':rec[1],'exit_code':rec[2],'allocated_CPUs':int(rec[3]),'allocation_wall_seconds':int(rec[4]),
          'allocated_core_seconds':int(rec[5]),'top_level_reported_CPU_seconds':seconds(rec[6]),'sum_distinct_step_reported_CPU_seconds':sum(seconds(row[6]) for row in steps.values()),
          'step_CPU_seconds':{key:seconds(row[6]) for key,row in steps.items()},'step_peak_RSS_KiB':{key:row[9] for key,row in steps.items()},'requested_memory':rec[7],
          'attempt_receipts':len(receipts),'process_returns':sum(r['returncode'] is not None for r in receipts),'normal_termination':sum(r['normal_termination'] for r in receipts),
          'all_attempt_summed_wall_seconds':sum(wall),'all_attempt_task_rank_seconds':sum(wall)*n,'executor_accounting':done,
          'collection':record(mp.parent/'collection.json') if (mp.parent/'collection.json').exists() else None})
    result={'design':record(a.design),'accounting':record(a.accounting),'ranks':ranks,'total_allocated_core_seconds':sum(x['allocated_core_seconds'] for x in ranks),
       'total_process_returns':sum(x['process_returns'] for x in ranks),'total_distinct_step_CPU_seconds':sum(x['sum_distinct_step_reported_CPU_seconds'] for x in ranks),
       'new_GPU_MACE_DFT_calls':0,'requested_GPU_seconds':0,'local_preparation_testing_reporting_CPU_unmetered_additional':True,
       'whole_allocation_includes_staging_collection_and_final_rank8_comparison':True,'performance_is_one_matched_batch_per_rank_not_universal_scaling':True}
    if result['total_process_returns']>24:raise ValueError('finite call scope exceeded')
    write_new(a.output,result)


if __name__=='__main__':main()
