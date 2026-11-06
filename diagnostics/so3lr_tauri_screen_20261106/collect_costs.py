"""Account for explicit owned IDs from scheduler records, including failures."""
import os,argparse,subprocess
from common import *
ap=argparse.ArgumentParser();ap.add_argument('--jobs',nargs='+',required=True);args=ap.parse_args()
assert all(j.isdecimal() for j in args.jobs)
fields=['JobIDRaw','JobName','State','Start','End','ElapsedRaw','AllocCPUS','CPUTimeRAW','TotalCPU','AllocTRES','MaxRSS']
cmd=['sacct','-n','-P','-j',','.join(args.jobs),'--format='+','.join(fields)]
proc=subprocess.run(cmd,env=dict(os.environ,TZ='UTC'),capture_output=True,text=True,check=True)
(W/'SCHEDULER_ACCOUNTING.tsv').write_text('\t'.join(fields)+'\n'+'\n'.join('\t'.join(line.split('|')[:len(fields)]) for line in proc.stdout.splitlines())+'\n')
def seconds(s):
 day=0
 if '-' in s:ds,s=s.split('-',1);day=int(ds)
 vals=[float(v) for v in s.split(':')]
 return day*86400+sum(v*60**i for i,v in enumerate(reversed(vals)))
rows=[]
for line in proc.stdout.splitlines():
 r=dict(zip(fields,line.split('|')))
 if r['JobIDRaw'] not in args.jobs:continue
 assert r['State'].split()[0].rstrip('+') in ('COMPLETED','FAILED','CANCELLED','OUT_OF_MEMORY','TIMEOUT','NODE_FAIL')
 r.update(allocated_CPU_seconds=int(r['ElapsedRaw'])*int(r['AllocCPUS']),measured_CPU_seconds=seconds(r['TotalCPU']))
 gpu=next((int(t.split('=')[-1]) for t in r['AllocTRES'].split(',') if t.startswith('gres/gpu=')),0)
 r['allocated_GPU_seconds']=int(r['ElapsedRaw'])*gpu;rows.append(r)
assert len(rows)==len(args.jobs),(len(rows),args.jobs)
save(RESULTS_DIR/'COMPUTE_COST.json',dict(jobs=rows,number_jobs=len(rows),
 allocated_CPU_seconds=sum(r['allocated_CPU_seconds'] for r in rows),
 measured_CPU_seconds=sum(r['measured_CPU_seconds'] for r in rows),
 allocated_GPU_seconds=sum(r['allocated_GPU_seconds'] for r in rows),
 raw_accounting=pin(W/'SCHEDULER_ACCOUNTING.tsv'),timezone='UTC',
 includes='All explicit owned scientific/audit jobs, output failures, parser checks and pending GPU cancellation.',
 excludes='Unmetered source/metadata preparation, downloads and new-environment setup. Other-session jobs are excluded.',
 monetary_cost='No billing rate supplied; report resource use without invented currency cost.'))
print('accounted',len(rows),'jobs')
