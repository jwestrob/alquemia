"""Read terminal pilot receipts; no molecular calls or source modification."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,write_new


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',required=True,type=Path);ap.add_argument('--accounting',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    m=read_json(a.run/'manifest.json');execution=read_json(a.run/'EXECUTION.json')
    receipts=[(p,read_json(p)) for p in a.run.rglob('execution.json')]
    native=[(p,read_json(p)) for p in a.run.rglob('mace_result.json')]
    points=[read_json(p) for p in (a.run/'searches').glob('*/result.json')]
    rows=a.accounting.read_text().splitlines();job=next(line.split('|') for line in rows if line.split('|')[0]==execution['job_id'])
    gpu_path=a.run/('gpu_summary_'+execution['job_id']+'.json')
    gpu=read_json(gpu_path) if gpu_path.exists() else None
    # Accounting columns: JobIDRaw,State,ExitCode,AllocCPUS,ElapsedRaw,CPUTimeRAW,ReqMem,AllocTRES,MaxRSS.
    costs={'job_id':job[0],'state':job[1],'exit_code':job[2],'allocated_CPUs':int(job[3]),'elapsed_seconds':int(job[4]),
        'allocated_core_seconds':int(job[5]),'requested_GPU_seconds':int(job[4]),'accounting':record(a.accounting),
        'requested_host_memory':job[6],
        'scheduler_step_peak_RSS_KiB':{line.split('|')[0]:line.split('|')[8] for line in rows if len(line.split('|'))>8 and line.split('|')[8]},
        'execution':record(a.run/'EXECUTION.json'),'manifest':record(a.run/'manifest.json'),
        'standalone_receipts':len(receipts),'standalone_model_calls_started':sum(r['model_call_started'] for _,r in receipts),
        'standalone_successful_returns':sum(r['returncode']==0 for _,r in receipts),
        'new_MACE_calls_started':sum(r['model_call_started'] for _,r in native),'new_MACE_complete':sum(r['status']=='complete' for _,r in native),
        'search_records':len(points),'search_status':dict(Counter(p['status'] for p in points)),
        'search_termination':dict(Counter(p['termination'] for p in points)),
        'distinct_composite_points':sum(len(p['points']) for p in points),
        'baseline_standalone_attempts':sum('/baseline/' in str(p) for p,_ in receipts),
        'search_standalone_attempts':sum('/searches/' in str(p) for p,_ in receipts),
        'cross_standalone_attempts':sum('/cross/' in str(p) for p,_ in receipts),
        'maximum_distinct_evaluations_used':max((len(p['points']) for p in points),default=0),
        'standalone_summed_call_wall_seconds':sum(r['wall_seconds'] for _,r in receipts),
        'new_MACE_summed_call_wall_seconds':sum(r['wall_seconds'] for _,r in native),
        'GPU_worker_summary':record(gpu_path) if gpu else None,
        'GPU_model_load_seconds':gpu['model_load_seconds'] if gpu else None,
        'GPU_worker_span_including_idle_seconds':gpu['wall_seconds'] if gpu else None,
        'peak_cuda_allocated_bytes':gpu['peak_cuda_allocated_bytes'] if gpu else None,
        'peak_GPU_worker_host_RSS_KiB':gpu['peak_host_RSS_KiB'] if gpu else None,
        'historical_baseline_reused_cells':56,'new_DFT_calls':0,'production_changed':False,
        'whole_allocation_includes_in_job_staging_startup_collection':True,
        'local_preflight_preparation_tests_reporting_CPU_unmetered_additional':True,
        'no_matched_hardware_speed_claim':True}
    if costs['standalone_receipts']>1000 or costs['new_MACE_calls_started']>480 or costs['maximum_distinct_evaluations_used']>40:raise RuntimeError('frozen scientific scope exceeded')
    write_new(a.output,costs);print(json.dumps(costs,indent=2))

if __name__=='__main__':main()
