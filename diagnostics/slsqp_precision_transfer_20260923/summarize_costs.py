"""Read actual terminal receipts; never executes a scientific calculation."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify,write_new


def main(submissions,output,sacct_output):
    sub=read_json(submissions);ids=[r[stage]['job_id'] for r in sub['shards'] for stage in ('GPU','CPU')]
    cmd=['sacct','-j',','.join(ids),'--noheader','--parsable2','--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,ReqMem,AllocTRES,MaxRSS']
    proc=subprocess.run(cmd,text=True,capture_output=True,check=True)
    Path(sacct_output).write_text(proc.stdout);rows={r[0]:r for line in proc.stdout.splitlines() if (r:=line.split('|'))[0] in ids}
    if set(rows)!=set(ids) or any(r[1] in ('RUNNING','PENDING','COMPLETING','CONFIGURING') for r in rows.values()):
        raise InvalidArtifact('all eight terminal allocations required')
    jobs=[];summaries=[];searches=[];low=[];events=[]
    for shard in sub['shards']:
        mp=verify(shard['GPU']['manifest']);m=read_json(mp);base=mp.parent.parent
        for stage in ('GPU','CPU'):
            jid=shard[stage]['job_id'];row=rows[jid];elapsed=int(row[2]);cpus=int(row[3]);gpus=1 if stage=='GPU' else 0
            jobs.append({'job_id':jid,'stage':stage,'state':row[1],'elapsed_seconds':elapsed,'allocated_CPUs':cpus,
                'allocated_core_seconds':elapsed*cpus,'requested_GPU_seconds':elapsed*gpus,'requested_GPUs':gpus,
                'requested_memory':row[4],'allocated_TRES':row[5]})
        for stage,folder in (('search',base/'proposals'),('cross',base/'pool')):
            for path in folder.glob('gpu_summary_*.json'):summaries.append({**read_json(path),'stage':stage,'artifact':record(path)})
        for t in m['tasks']:
            path=mp.parent/'proposals'/t['task_id']/'result.json';r=read_json(path) if path.exists() else None
            oldsrc=next(s for s in m['sources'] if s['case_id']==t['case_id']);old=read_json(verify(oldsrc['old_receipts'][t['metal']]))
            searches.append({'case_id':t['case_id'],'metal':t['metal'],'receipt':record(path) if r else None,
                'status':r['status'] if r else 'not_run','wall_seconds':r.get('wall_seconds') if r else None,
                'optimizer':r.get('optimizer') if r else None,'boundary_flag':r.get('boundary_flag') if r else None,
                'reused_request_count':sum(q.get('reused',False) for q in r['requests']) if r else None,
                'scientific_origin_reused_count':sum(q.get('scientific_origin_reused',False) for q in r['requests']) if r else None,
                'old_status':old['status'],'old_wall_seconds':old['wall_seconds'],'old_optimizer':old['optimizer']})
        lowroot=base/'pool/solvent/shard_0';lm=read_json(lowroot/'manifest.json')
        for t in lm['tasks']:
            path=Path(t['output_path']+'.execution.json');r=read_json(path) if path.exists() else None
            low.append({'task_id':t['task_id'],'receipt':record(path) if r else None,'normal_termination':r['normal_termination'] if r else None,
                'scf_converged':r['scf_converged'] if r else None,'returncode':r['returncode'] if r else None,
                'ranks':r['parallelism']['nprocs'] if r else None})
        event=lowroot/'budget_events.jsonl'
        if event.exists():events.append({'artifact':record(event),'events':[json.loads(line) for line in event.read_text().splitlines()]})
    totals={'allocated_core_seconds':sum(j['allocated_core_seconds'] for j in jobs),
        'requested_GPU_seconds':sum(j['requested_GPU_seconds'] for j in jobs),
        'new_MACE_search_calls':sum(r['new_MACE_calls'] for r in summaries if r['stage']=='search'),
        'new_MACE_cross_calls':sum(r['new_MACE_calls'] for r in summaries if r['stage']=='cross'),
        'new_GFN2_attempt_receipts':sum(r['receipt'] is not None for r in low),
        'normal_GFN2_terminations':sum(r['normal_termination'] is True for r in low),
        'new_search_count':len(searches),'available_searches':sum(r['status']=='proposal_available' for r in searches),
        'optimizer_function_evaluations':sum(r['optimizer']['function_evaluations'] for r in searches if r['optimizer']),
        'search_wall_seconds':sum(r['wall_seconds'] for r in searches if r['wall_seconds'] is not None),
        'old_optimizer_function_evaluations_same202_sources':sum(r['old_optimizer']['function_evaluations'] for r in searches),
        'old_search_wall_seconds_same202_sources':sum(r['old_wall_seconds'] for r in searches),
        'new_q0_calls':0,'new_DFT_calls':0,'six_reused_pools_cost_excluded':True}
    result={'submissions':record(submissions),'sacct':record(sacct_output),'jobs':jobs,'GPU_summaries':summaries,
        'searches':searches,'GFN2_receipts':low,'low_execution_events':events,'totals':totals,'implementation':record(__file__),
        'cost_scope':'Eight new allocations, including in-allocation validation/loading/collection. Shared-node concurrency changes elapsed latency, not summed allocation. Local preparation/tests/report time and historical reuse cost are separate and unmetered here.'}
    write_new(output,result);print(json.dumps(totals,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('submissions','output','sacct_output'):p.add_argument('--'+key.replace('_','-'),required=True,type=Path)
    main(**vars(p.parse_args()))
