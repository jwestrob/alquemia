"""Read exact outer allocation costs, including preserved empty failed stages."""
from pathlib import Path
from collections import Counter
import argparse,json,subprocess,sys,re
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,write_new,verify

def main(workspace,output):
    b=Path(workspace);o=Path(output);o.mkdir(parents=True,exist_ok=False)
    receipts=sorted(b.glob('run_v2/SUBMISSION_*.json'))+sorted(b.glob('run_v[12]/shard_*/SUBMISSION_*.json'))+sorted(b.glob('run_v1/shard_*/origins/SUBMISSION_*.json'))
    jobs={read_json(p)['job_id']:p for p in receipts}
    command=['sacct','-j',','.join(jobs),'-X','--parsable2','-n','--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,AllocTRES,NodeList,Start,End']
    text=subprocess.check_output(command,text=True);(o/'SACCT.txt').write_text(text);rows=[]
    for line in text.splitlines():
        f=line.split('|');job,state,elapsed,cpus,tres,node,start,end=f[:8]
        if job not in jobs:continue
        if state in ('RUNNING','PENDING','CONFIGURING','COMPLETING'):raise RuntimeError('allocation not terminal: '+job)
        allocation=dict(x.split('=',1)for x in tres.split(',')if '='in x)
        submission=read_json(jobs[job]);wrapper=verify(submission['wrapper']).read_text()
        requested=re.findall(r'^#SBATCH --gres=gpu:(\d+)$',wrapper,re.M)
        ngpu=int(requested[0])if requested else 0
        if any(x.startswith('--gres')for x in submission['command']):raise RuntimeError('unaccounted command GPU override')
        rows.append({'job_id':job,'state':state,'elapsed_seconds':int(elapsed),'allocated_CPUs':int(cpus),
            'allocated_CPU_seconds':int(elapsed)*int(cpus),'requested_GPU_seconds':int(elapsed)*ngpu,
            'requested_GPUs_from_pinned_wrapper':ngpu,'allocation_TRES':allocation,'node':node,'start':start,'end':end,'submission':record(jobs[job])})
    if {r['job_id']for r in rows}!=set(jobs):raise RuntimeError('accounting row missing')
    result={'jobs':rows,'job_count':len(rows),'allocated_CPU_seconds':sum(r['allocated_CPU_seconds']for r in rows),
        'requested_GPU_seconds':sum(r['requested_GPU_seconds']for r in rows),'summed_elapsed_not_makespan_seconds':sum(r['elapsed_seconds']for r in rows),
        'states':dict(Counter(r['state']for r in rows)),'GPU_accounting':'requested GPU-seconds from pinned SBATCH wrapper; this cluster omits GPU from AllocTRES, not utilization','original_empty_stages_included':True,'nested_runner_cost_added':False,
        'sacct':record(o/'SACCT.txt'),'command':command,'implementation':record(__file__)}
    write_new(o/'COSTS.json',result);print(json.dumps({k:v for k,v in result.items()if k!='jobs'},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workspace',required=True);p.add_argument('--output',required=True);main(**vars(p.parse_args()))
