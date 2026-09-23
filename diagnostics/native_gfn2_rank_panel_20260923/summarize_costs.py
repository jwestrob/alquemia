"""Summarize actual full-job/step accounting and fixed-panel execution receipts."""
import argparse
import datetime as dt
import json
from pathlib import Path
import statistics
import sys

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,write_new


def seconds(value):
    days=0
    if '-' in value:
        day,value=value.split('-',1);days=int(day)
    return days*86400+sum(float(x)*60**i for i,x in enumerate(reversed(value.split(':'))))


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--accounting',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    m=read_json(a.manifest);root=a.manifest.parent
    job=read_json(root/'SUBMISSION.json')['stdout'].strip().split(';')[0]
    rows={v[0]:v for line in a.accounting.read_text().splitlines() if (v:=line.split('|'))}
    r=rows[job];steps={k:v for k,v in rows.items() if k.startswith(job+'.')}
    receipts=[];times=[]
    for t in m['tasks']:
        rp=Path(t['output_path']+'.execution.json')
        if not rp.exists():continue
        receipt=read_json(rp);receipts.append(receipt)
        parse=lambda s:dt.datetime.fromisoformat(s.replace('Z','+00:00'))
        if receipt['finished_at_utc']:
            times.append((parse(receipt['finished_at_utc'])-parse(receipt['started_at_utc'])).total_seconds())
    events=[json.loads(line) for line in (root/'budget_events.jsonl').read_text().splitlines()]
    execution=next((x for x in reversed(events) if 'allocated_core_seconds' in x),None)
    result={'manifest':record(a.manifest),'accounting':record(a.accounting),'job_id':job,
        'state':r[1],'exit_code':r[2],'allocated_CPUs':int(r[3]),'allocation_wall_seconds':int(r[4]),
        'allocated_core_seconds':int(r[5]),'top_level_reported_CPU_seconds':seconds(r[6]),
        'sum_distinct_step_reported_CPU_seconds':sum(seconds(v[6]) for v in steps.values()),
        'step_reported_CPU_seconds':{k:seconds(v[6]) for k,v in steps.items()},
        'step_peak_RSS_KiB':{k:v[9] for k,v in steps.items()},'requested_memory':r[7],
        'attempt_receipts':len(receipts),'process_returns':sum(x['returncode'] is not None for x in receipts),
        'normal_terminations':sum(x['normal_termination'] for x in receipts),
        'summed_attempt_wall_seconds':sum(times),'median_attempt_wall_seconds':statistics.median(times) if times else None,
        'executor_accounting':execution,'GPU_MACE_DFT_calls':0,'requested_GPU_seconds':0,
        'whole_allocation_includes_collection':True,'matched_hardware_speed_claim':False,
        'local_preparation_tests_reporting_are_additional_unmetered_CPU':True}
    if result['attempt_receipts']>336:raise ValueError('finite cell count exceeded')
    write_new(a.output,result)


if __name__=='__main__':main()
