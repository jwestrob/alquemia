"""Summarize the five actual finite pilot allocations and saved stage timers."""
import csv,json,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,write_new
BASE=ROOT/'workspaces/motion_envelope_20260923'
STAGES=[('origin_MACE','pilot34_origins_v1',1211294,1),('native_search','pilot34_searches_v1',1211342,1),
        ('origin_scalar','pilot34_scalar_origins_v1',1211343,0),('cross_MACE','pilot34_pool_v1',1211421,1),
        ('candidate_scalar','pilot34_pool_v1/solvent/shard_0',1211422,0)]

def main(output):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    command=['sacct','-X','-j',','.join(str(s[2])for s in STAGES),'--parsable2','--format=JobID,JobName,Partition,State,ExitCode,ElapsedRaw,AllocCPUS,ReqMem,AllocTRES,NodeList,Submit,Start,End']
    text=subprocess.run(command,check=True,text=True,capture_output=True).stdout
    rows={int(r['JobID']):r for r in csv.DictReader(text.splitlines(),delimiter='|')}
    if set(rows)!={s[2]for s in STAGES} or any(r['State']in('RUNNING','PENDING','COMPLETING')for r in rows.values()):raise InvalidArtifact('all five allocations must be terminal')
    (out/'SACCT_FINAL.txt').write_text(text);stages=[]
    for label,relative,job,gpus in STAGES:
        root=BASE/relative;row=rows[job];timers={};pins=[]
        for pattern in ('EXECUTION_*.json','MACE_EXECUTION_*.json','gpu_summary_*.json','*execution_summary*.json','execution*.json','mace/execution_*.json'):
            for p in root.glob(pattern):
                d=read_json(p);pins.append(record(p));timers[p.name]={k:v for k,v in d.items()if k in('wall_seconds','elapsed_seconds','model_load_seconds','evaluation_wall_seconds','new_MACE_calls','failed_requests','error','job_id','device')}
        if label in ('origin_scalar','candidate_scalar'):
            log=ROOT/'diagnostics/motion_envelope_20260923'/((('origin_cpu_'if label=='origin_scalar'else'candidate_cpu_')+str(job)+'.out'))
            textlog=log.read_text();pins.append(record(log))
            for match in re.finditer(r'^\{',textlog,re.M):
                try:d,_=json.JSONDecoder().raw_decode(textlog[match.start():])
                except ValueError:continue
                if 'elapsed_seconds'in d:timers['native_runner']={k:d[k]for k in ('elapsed_seconds','allocated_core_seconds')if k in d}
        if label=='native_search':
            proposals=[read_json(p)for p in root.glob('proposals/*/result.json')];evaluations=[read_json(p)for p in root.glob('proposals/*/evaluations/*/mace_result.json')]
            sums={'search_count':len(proposals),'available_proposals':sum(r['status']=='proposal_available'for r in proposals),
                  'boundary_proposals':sum(r.get('boundary_flag',False)for r in proposals),'sum_per_search_wall_seconds':sum(r['wall_seconds']for r in proposals),
                  'maximum_per_search_wall_seconds':max(r['wall_seconds']for r in proposals),'MACE_evaluation_count':len(evaluations),
                  'sum_MACE_evaluation_wall_seconds':sum(r['wall_seconds']for r in evaluations)}
            exe=read_json(root/f'EXECUTION_{job}.json');gpu=read_json(root/f'gpu_summary_{job}.json')
            sums['executor_minus_summed_search_wall_seconds']=exe['wall_seconds']-sums['sum_per_search_wall_seconds']
            sums['executor_minus_recorded_GPU_worker_interval_seconds']=exe['wall_seconds']-gpu['wall_seconds'];timers['search_timer_aggregation']=sums
        stages.append({'stage':label,'job_id':job,'scheduler':row,'allocated_core_seconds':int(row['ElapsedRaw'])*int(row['AllocCPUS']),
                       'requested_GPU_seconds':int(row['ElapsedRaw'])*gpus,'outer_elapsed_seconds':int(row['ElapsedRaw']),'timer_receipts':pins,'timers':timers})
    result={'stages':stages,'allocated_core_seconds':sum(s['allocated_core_seconds']for s in stages),'requested_GPU_seconds':sum(s['requested_GPU_seconds']for s in stages),
            'outer_elapsed_seconds_summed_not_makespan':sum(s['outer_elapsed_seconds']for s in stages),'scheduler_command':command,'scheduler_receipt':record(out/'SACCT_FINAL.txt'),
            'implementation':record(__file__),'additional_local_preparation_wall_seconds':73.1082,
            'limits':['Nested runner allocations are not added to outer allocation cost.','Local preparation/preflight/analysis is additional and lacks an allocated-core receipt.',
                      'Model load is measured once per batch; other initialization/collection residuals are not assigned an unmeasured cause.',
                      'These34-source batch timings do not measure an isolated three-source request. GPU-seconds mean requested allocation, not active utilization.']}
    write_new(out/'COSTS_v1.json',result);print(json.dumps({k:v for k,v in result.items()if k in('allocated_core_seconds','requested_GPU_seconds','outer_elapsed_seconds_summed_not_makespan')},indent=2))
if __name__=='__main__':main(sys.argv[1])
