#!/usr/bin/env python3
"""Eight analytic programs using the frozen joint evaluator's pure parser."""
import argparse
import concurrent.futures
import copy
import importlib.util
import json
import os
from pathlib import Path
import resource
import threading
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('radial_binding',HERE/'BIND.py');B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B)

def parse(ctx,rec,row):return ctx[2].Evaluator._parse(ctx[9],rec,row['runtime'])

def inputs(ctx,rec,row):
    ec=ctx[3];B.require(set(rec['program_inputs'])=={'analyze','testgrad'} and len(rec['receipts'])==2 and rec['render']==row['render'],'Exact two native program inputs')
    for program,receipt in zip(('analyze','testgrad'),rec['receipts']):
        inp=rec['program_inputs'][program];folder=B.verify(inp['xyz']).parent
        for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
            B.require(B.verify(inp[k])==folder/n and inp[k]['sha256']==row['files'][k]['sha256'],'Literal runtime input')
        args=['system.xyz','EM'] if program=='analyze' else ['system.xyz','Y','N']
        B.require(receipt['argv']==[ec['executables'][program]['path']]+args and receipt['cwd']==str(folder) and receipt['output']==B.pin(folder/(program+'.out')) and receipt['returncode']==0,'Exact native argv/output')
        B.require(set(receipt['threads'])==set(B.THREADS) and all(receipt['threads'][k]=='1' for k in B.THREADS),'Single-thread native receipt')
    B.require(B.verify(rec['program_inputs']['analyze']['xyz']).parent!=B.verify(rec['program_inputs']['testgrad']['xyz']).parent,'Separate program folders')

def admit_result(bp,result):
    b,ctx,r=B.release(bp,result['release'],'model')
    B.require(result['implementation']==B.pin(__file__) and result['bound']==bp and result['candidate_result']==b['candidate_result'] and result['selected']==b['selected'],'Actual model candidate/source')
    B.require(set(result['rows'])==set(B.IDS) and 0<=result['native_program_invocations']<=8 and result['new_FD_energies']==0,'Actual eight-program/four-cell scope')
    rows={}
    for gid,rec in result['rows'].items():
        try:
            B.require(rec['status']=='complete','Native cell incomplete');inputs(ctx,rec,b['rows'][gid]);p=parse(ctx,rec,b['rows'][gid])
            B.require(p['model']==rec['model'] and p['energy']==rec['energy'] and p['local_response_readback']==rec['local_response_readback'],'Raw E/F and local response replay')
            B.require(p['model']['native']['finite_difference'] is None,'No reserved Cartesian finite differences')
            rows[gid]=dict(status='complete',model=p['model'],receipts=rec['receipts'],native_readback=p['local_response_readback'])
        except Exception as exc:rows[gid]=dict(status='unavailable',model=None,reason=repr(exc))
    complete=sum(v['status']=='complete' for v in rows.values())
    B.require(complete<4 or result['status']=='water_coordination_radial_model_complete' and result['error'] is None and result['global_stop'] is False and result['native_program_invocations']==8,'Contradictory terminal model result')
    return dict(status='saved_model_admission_complete' if complete==4 else 'saved_model_admission_partial',rows=rows,complete=complete,anchor_model=b['anchor_model'],new_native_calls=0)

def run(bp,rp,out):
    tick,cpu=time.monotonic(),time.process_time();children=resource.getrusage(resource.RUSAGE_CHILDREN)
    b,ctx,r=B.release(bp,rp,'model');ec=ctx[3];native=ctx[9].h['native_helper']
    B.require(os.environ.get('SLURM_JOB_ID','').isdigit() and os.environ.get('SLURM_NTASKS','1')=='1' and os.environ.get('SLURM_CPUS_PER_TASK')=='4' and all(os.environ.get(k)=='1' for k in B.THREADS),'Exactly four single-thread Slurm workers')
    out=Path(out).resolve();B.require(str(out)==r['output_directory'] and not out.exists(),'Fresh exact released model output directory');out.mkdir(parents=True)
    result=dict(status='running',implementation=B.pin(__file__),bound=bp,release=rp,candidate_result=b['candidate_result'],selected=b['selected'],parameters=b['parameters'],
                rows={g:dict(status='not_attempted',receipts=[]) for g in B.IDS},native_program_invocations=0,new_FD_energies=0,error=None,job_id=os.environ['SLURM_JOB_ID'])
    lock=threading.RLock();stop=threading.Event()
    def persist():B.save(out/'RESULT.json',result)
    def worker(gid):
        rec=dict(status='not_attempted',cell_id=gid,receipts=[],program_inputs={},render=b['rows'][gid]['render']);claimed=False
        if stop.is_set():return
        try:
            row=b['rows'][gid]
            for program in ('analyze','testgrad'):
                if stop.is_set():return
                folder=out/gid/program;folder.mkdir(parents=True,exist_ok=False)
                for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):(folder/n).write_bytes(B.verify(row['files'][k]).read_bytes())
                inp={k:B.pin(folder/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
                B.require(all(inp[k]['sha256']==row['files'][k]['sha256'] for k in inp),'Copied model input differs')
                with lock:
                    if stop.is_set():return
                    B.require(result['native_program_invocations']<8,'Eight-program ceiling');claimed=True
                    rec.update(status='started',active_program=program);rec['program_inputs'][program]=inp
                    B.save(folder/'START.json',dict(program=program,input=inp,executable=ec['executables'][program],bound=bp,job_id=result['job_id']))
                    result['native_program_invocations']+=1;result['rows'][gid]=copy.deepcopy(rec);persist()
                argv=[ec['executables'][program]['path'],'system.xyz']+(['EM'] if program=='analyze' else ['Y','N'])
                native.invoke(argv,folder,folder/(program+'.out'),rec['receipts'])
                rec.update(parse(ctx,rec,row))
                with lock:result['rows'][gid]=copy.deepcopy(rec);persist()
            rec['status']='complete'
        except Exception as exc:
            stop.set();rec.update(status='failed',error=repr(exc))
            with lock:result['error']=repr(exc);result['rows'][gid]=copy.deepcopy(rec);persist()
        finally:
            if claimed:
                with lock:result['rows'][gid]=copy.deepcopy(rec);persist()
    try:
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(worker,B.IDS))
        B.require(not stop.is_set() and all(r['status']=='complete' for r in result['rows'].values()) and result['native_program_invocations']==8,'Complete eight-program panel')
        result['status']='water_coordination_radial_model_complete'
    except Exception as exc:stop.set();result.update(status='native_or_admission_failure',error=repr(exc))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN)
        result.update(global_stop=stop.is_set(),elapsed_seconds=time.monotonic()-tick,parent_CPU_seconds=time.process_time()-cpu,aggregate_child_user_CPU_seconds=ch.ru_utime-children.ru_utime,aggregate_child_system_CPU_seconds=ch.ru_stime-children.ru_stime,per_program_child_CPU_may_overlap_not_summed=True);persist()
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bound',required=True);p.add_argument('--release',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    r=run(B.pin(a.bound),B.pin(a.release),a.out);print(json.dumps({'status':r['status'],'error':r['error']}));raise SystemExit(0 if r['status']=='water_coordination_radial_model_complete' else 1)
