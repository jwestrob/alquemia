#!/usr/bin/env python3
"""Frozen paired-pool code with one explicit certified-search admission bridge."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import time
import traceback
import types

PROTOCOL='lbt_hybrid_response_one_certified_slot_paired_pool_v1'
RESPONSE_SHA='05b1c412feebfc72048405c8723d426168872325ce32d956732f5c898c028981'
RECOVERY_SHA='ba02f7ac712db0e198a5f77b68518ff7838c14179488b2c61d70d3e35b633bbf'
STATES=(('7CCO','La'),('7CCO','Dy'),('7CCN','La'),('7CCN','Dy'))


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(path):
    p=Path(path).resolve()
    with p.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
    return {'path':str(p),'sha256':sha}


def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'],'Changed pin: '+p['path'])
    return Path(p['path'])


def read(p):return json.loads(verify(p).read_text())


def write(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
    return pin(path)


def same(a,b):
    verify(a);verify(b)
    return a['sha256']==b['sha256']


def module(name,p):
    spec=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m);return m


def prepare(cfg,out,result):
    """Original complete-parent preparation; no scientific model or gate change."""
    require(cfg['protocol']==PROTOCOL and same(cfg['implementation'],pin(__file__)), 'Pool source/protocol differs')
    verify(cfg['specification'])
    require(cfg['response_implementation']['sha256']==RESPONSE_SHA
            and cfg['recovery_implementation']['sha256']==RECOVERY_SHA,'Pinned protocol dependency differs')
    response=module('certified_pool_original_response',cfg['response_implementation'])
    recovery=module('certified_pool_explicit_recovery',cfg['recovery_implementation'])
    original=read(cfg['original_response_config'])
    require(same(original['implementation'],cfg['response_implementation']) and cfg['region']=='baseline'
            and cfg['states']==[list(s) for s in STATES] and len(cfg['run_results'])==4
            and cfg['recovered_result'] is not None
            and same(cfg['recovered_result'],cfg['run_results'][2]), 'Exact four-state recovery binding required')
    for key in original:
        if key not in ('protocol','implementation','specification','run_results'):
            require(recovery.normalized(cfg[key])==recovery.normalized(original[key]),
                    'Original response scientific config differs: '+key)
    r=read(cfg['recovered_result']);rcfg=read(r['config'])
    require(r['protocol']==recovery.PROTOCOL and same(r['implementation'],cfg['recovery_implementation'])
            and same(rcfg['implementation'],cfg['recovery_implementation'])
            and same(rcfg['original_response_config'],cfg['original_response_config']),
            'Recovered slot implementation/original source differs')
    d=out/'original_response_replay';d.mkdir()
    loaded=response.prepare(original,d,result)
    loaded['_bridge']={'original':response,'recovery':recovery,'original_config':original,
                       'config':copy.deepcopy(cfg),'out':out}
    result.update(original_response_implementation=cfg['response_implementation'],
        recovery_implementation=cfg['recovery_implementation'],recovered_result=cfg['recovered_result'],
        preserved_failed_result=rcfg['original_failed_result'],
        saved_admission_policy={'ordinary_slots':[list(s) for s in STATES if s!=('7CCN','La')],
            'certified_slot':['7CCN','La'],'unknown_protocols_allowed':False})
    return loaded


def saved_search(record,cfg,loaded,state,stages):
    """Exactly one new protocol; all other slots retain frozen original admission."""
    bridge=loaded['_bridge'];state=tuple(state)
    require(state in STATES and stages==('run',) and same(record,cfg['run_results'][STATES.index(state)])
            and same(cfg['implementation'],bridge['config']['implementation']), 'Ordered saved search slot differs')
    r=read(record)
    if state==('7CCN','La'):
        require(same(record,cfg['recovered_result']) and r['protocol']==bridge['recovery'].PROTOCOL,
                'Only the explicitly certified new result can fill the recovered slot')
        return bridge['recovery'].admit_recovered_search(record,bridge['original_config'],loaded,
                                                       bridge['out']/'recovered_search_admission')
    require(r['protocol']==bridge['original'].PROTOCOL,'Unknown protocol in ordinary response slot')
    return bridge['original'].saved_search(record,bridge['original_config'],loaded,state,stages)


def paired_pool(cfg,loaded,out,result):
    original=loaded['_bridge']['original']
    require(pin(original.__file__)['sha256']==RESPONSE_SHA,'Frozen pool source changed')
    # Only this function object's global name lookup is adapted. No module
    # globals or historical source are modified; its executable code is exact.
    private=dict(original.paired_pool.__globals__)
    private['saved_search']=saved_search
    adapted=types.FunctionType(original.paired_pool.__code__,private,
                               original.paired_pool.__name__,original.paired_pool.__defaults__,original.paired_pool.__closure__)
    adapted(cfg,loaded,out,result)
    result['paired_pool_code_unchanged']=True
    result['only_private_global_substitution']='saved_search'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=('preflight','pool'),required=True)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    tick=time.monotonic();before=resource.getrusage(resource.RUSAGE_SELF)
    result={'protocol':PROTOCOL,'stage':a.stage,'status':'execution_failure','implementation':pin(__file__),
        'config':pin(a.config),'slurm_job_id':os.environ['SLURM_JOB_ID'],'region':'baseline',
        'declared_search_states':[list(s) for s in STATES],'declared_search_denominator':4,
        'maximum_new_cross_composites':4,'new_geometries':0,'optimizer_steps':0,'new_DFT_calls':0,
        'new_pool_available':False}
    try:
        cfg=json.loads(a.config.read_text())
        require(cfg['protocol']==PROTOCOL and same(cfg['implementation'],pin(__file__)),'Pool source/protocol differs')
        if cfg.get('recovered_result') is None or len(cfg.get('run_results',[]))!=4 or any(p is None for p in cfg['run_results']):
            result.update(status='unbound',missing_required_evidence='four bound search results including the certified recovery')
        else:
            loaded=prepare(cfg,out,result)
            if a.stage=='preflight':
                slots=[{'source_id':s,'metal':m,'status':'unavailable','result':r}
                       for (s,m),r in zip(STATES,cfg['run_results'])]
                result['terminal_slots']=slots
                for state,slot in zip(STATES,slots):
                    try:
                        _,s,_,point=saved_search(slot['result'],cfg,loaded,state,('run',))
                        slot.update(status='admitted',point=point['provenance'],termination=s['termination'],converged=s['converged'])
                    except Exception as exc:slot.update(failure=repr(exc),traceback=traceback.format_exc())
                require(loaded['ev'].cost['composites_attempted']==0 and loaded['ev'].native is None,
                        'Saved four-slot preflight made a molecular call')
                admitted=all(s['status']=='admitted' for s in slots)
                result.update(status='complete' if admitted else 'unavailable_terminal_evidence',all_four_admitted=admitted)
            else:
                require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles'),'Allocated GPU required')
                paired_pool(cfg,loaded,out,result)
        verify(result['config'])
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc(),new_pool_available=False)
    ev=result.pop('_evaluator',None)
    if ev is not None:
        result.update(actual_cost=ev.cost,requests=ev.requests,points=ev.points)
        if ev.native is not None:result.update(native_MACE_calls=ev.native.calls,native_MACE_device=ev.native.device,native_MACE_versions=ev.native.versions)
    after=resource.getrusage(resource.RUSAGE_SELF)
    result.update(elapsed_seconds=time.monotonic()-tick,parent_CPU_seconds=(after.ru_utime+after.ru_stime)-(before.ru_utime+before.ru_stime),
        allocation={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','CUDA_VISIBLE_DEVICES')})
    write(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')}))
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
