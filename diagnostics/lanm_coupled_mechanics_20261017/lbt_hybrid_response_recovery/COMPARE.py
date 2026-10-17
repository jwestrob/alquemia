#!/usr/bin/env python3
"""Saved comparison with one explicitly certified response-search protocol."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
import traceback
import types

PROTOCOL='lanm_hybrid_saved_comparison_certified_response_v1'
ORIGINAL_SHA='cc617b98eeedef6e42786ca763278b5ab477349c5969633873dda034dcca407d'
POOL_SHA='91487e1390e81f12aa6651d0269a874adbe555dd05981f81a3b37a2ba6ae8d6c'
POOL_PROTOCOL='lbt_hybrid_response_one_certified_slot_paired_pool_v1'


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


def compare(cfg,out,result):
    require(cfg['protocol']==PROTOCOL and cfg['implementation']['sha256']==pin(__file__)['sha256']
            and cfg['original_comparison_implementation']['sha256']==ORIGINAL_SHA
            and cfg['response_implementation']['sha256']==POOL_SHA,'Explicit comparison bridge differs')
    spec=importlib.util.spec_from_file_location('certified_response_frozen_comparison',verify(cfg['original_comparison_implementation']))
    original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
    verify(cfg['original_comparison_specification'])
    # The full original comparison/replay/selection code is retained. Its own
    # private declaration names this new source/protocol and the new pool reader.
    private=dict(original.compare.__globals__)
    private.update(PROTOCOL=PROTOCOL,__file__=__file__,FROZEN=dict(original.FROZEN,response_implementation=POOL_SHA))
    adapted=types.FunctionType(original.compare.__code__,private,original.compare.__name__,
                               original.compare.__defaults__,original.compare.__closure__)
    adapted(cfg,out,result)
    result['original_comparison_implementation']=cfg['original_comparison_implementation']
    result['comparison_math_and_replay_code_unchanged']=True
    result['private_declaration_changes']=['implementation identity','protocol identity','explicit response-pool implementation']
    if result['status']=='complete':
        pool=read(cfg['response_pool_result']);pcfg=read(pool['config'])
        require(pool['protocol']==POOL_PROTOCOL and len(pool['terminal_slots'])==4
                and all(s['status']=='admitted' for s in pool['terminal_slots']), 'Complete new-protocol four-slot pool required')
        recovered=read(pcfg['recovered_result']);rcfg=read(recovered['config'])
        old=read(rcfg['original_failed_result'])
        require(old['status']=='unavailable' and old['slurm_job_id']=='13019'
                and recovered['protocol']=='lbt_response_single_embedded_solver_recovery_v1'
                and recovered['status']=='complete' and recovered['continuation_execution_and_physical_admission'],
                'Explicit recovered/failed lineage differs')
        lineage={'recovered_result':pcfg['recovered_result'],'recovery_implementation':pcfg['recovery_implementation'],
            'original_failed_result':rcfg['original_failed_result'],'original_failed_SEARCH':old['search'],
            'original_failed_request':recovered['original_failed_request'],'old_failure_remains_unavailable':True,
            'solver_exception_certificate':recovered['solver_exception_certificate'],
            'certification_result':rcfg['scc_result'],'no_global_solver_fallback':True}
        result['explicit_response_recovery']=lineage
        result['parent_costs']['response_recovery_work_separate']={
            'original_failed_result':rcfg['original_failed_result'],'original_failed_cost':old['actual_cost'],
            'certification_result':rcfg['scc_result'],'certification_cost':recovered['certification_cost'],
            'recovery_result':pcfg['recovered_result'],'completion_component_cost':recovered['completion_actual_cost'],
            'continuation_component_cost':recovered['continuation_actual_cost'],
            'completion_and_continuation_partition_the_already_listed_recovery_allocation':True,
            'do_not_add_the_partition_to_response_search_cost_twice':True}
    original.report(out,result)
    report=out/'REPORT.md'
    if result['status']=='complete':
        with report.open('a') as f:
            f.write('\nThis comparison uses an explicit certified recovery protocol for one 7CCN/La search. '
                    'Its original failed search and four-call SCC diagnostic remain separately recorded. '
                    'No global solver fallback or changed numerical comparison is implied.\n')
    result['report']=pin(report)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Saved numerical comparison requires Slurm')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    tick,cpu=time.monotonic(),time.process_time()
    result={'protocol':PROTOCOL,'status':'execution_failure','implementation':pin(__file__),'config':pin(a.config),
        'slurm_job_id':os.environ['SLURM_JOB_ID'],'new_molecular_calls':0,'new_geometries':0,'optimizer_steps':0,
        'declared_sources':['7CCO','7CCN'],'declared_stages':['Module2_whole','hybrid_same_Module2_selection','hybrid_parent_pool','hybrid_response_pool'],
        'declared_response_states':[['7CCO','La'],['7CCO','Dy'],['7CCN','La'],['7CCN','Dy']]}
    try:
        cfg=json.loads(a.config.read_text());compare(cfg,out,result);verify(result['config'])
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick,process_CPU_seconds=time.process_time()-cpu)
    with (out/'RESULT.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')}))
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
