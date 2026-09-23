"""Actual uniform32 allocation and molecular receipt cost; no model calls."""
import argparse,json,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new

p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--job',required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();receipt_rows=[];runtime=0.;counts={}
for stage in (1,2):
    m=read_json(a.run/f'stage{stage}/manifest.json');rows=[]
    for t in m['tasks']:
        output=Path(t['output_path']);receipt=Path(str(output)+'.execution.json')
        r={'stage':stage,'task_id':t['task_id'],'receipt':record(receipt) if receipt.exists() else None,'ORCA_seconds':None}
        if receipt.exists():
            rc=read_json(receipt);r.update(returncode=rc['returncode'],normal_termination=rc['normal_termination'],scf_converged=rc['scf_converged'],ranks=rc['parallelism']['nprocs'])
        if output.exists():
            hits=re.findall(r'TOTAL RUN TIME: (\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec',output.read_text())
            if len(hits)==1:
                d,h,minute,s,ms=map(int,hits[0]);r['ORCA_seconds']=d*86400+h*3600+minute*60+s+ms/1000;runtime+=r['ORCA_seconds']
        rows.append(r)
    counts[str(stage)]={'declared_fresh':len(m['tasks']),'receipts':sum(r['receipt'] is not None for r in rows),'reused':len(m['reused']),'missing':len(m['missing'])};receipt_rows.extend(rows)
raw=subprocess.check_output(['sacct','-j',a.job,'--parsable2','--noheader','--format=JobID,State,ElapsedRaw,AllocCPUS,ReqMem,AllocTRES,MaxRSS'],text=True)
row=next(x.split('|') for x in raw.splitlines() if x.split('|')[0]==a.job)
if row[1] not in ('COMPLETED','FAILED','CANCELLED','TIMEOUT'):raise RuntimeError('allocation not terminal')
result={'job_id':a.job,'state':row[1],'wall_seconds':int(row[2]),'allocated_cpus':int(row[3]),'allocated_core_seconds':int(row[2])*int(row[3]),'requested_memory':row[4],'actual_allocation':row[5],'GPU_seconds':0,'stage_counts':counts,'summed_ORCA_seconds':runtime,'per_call_receipts':receipt_rows,'sacct':raw,'new_MACE_DFT_optimizations':0,'historical_reuse_cost_included':False,'local_preparation_test_report_metered':False}
write_new(a.output,result);print(json.dumps({k:v for k,v in result.items() if k!='per_call_receipts'},indent=2))
