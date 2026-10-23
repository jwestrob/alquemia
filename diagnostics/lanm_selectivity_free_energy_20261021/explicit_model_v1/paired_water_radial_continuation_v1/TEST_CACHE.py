#!/usr/bin/env python3
"""Read-only inherited-cache and finite-budget controls; native invoke forbidden."""
import argparse
import concurrent.futures
import copy
import importlib.util
import json
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('continuation_evaluator_controls',HERE/'EVALUATE.py');e=importlib.util.module_from_spec(s);s.loader.exec_module(e)


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False);cfg=json.loads((HERE/'CONFIG_EVALUATOR.json').read_text());prior=e.read(cfg['prior_result']);checks=[]
    def check(ok,name):
        assert ok,name;checks.append({'case':name,'passed':True})
    ev=e.Evaluator(cfg,out/'cache');actual_invokes=[]
    def forbidden(*a,**k):actual_invokes.append(True);raise AssertionError('Native invoke forbidden by saved control')
    ev.native.invoke=forbidden
    current=prior['selected']['parameters'];m=ev.evaluate(current)
    check(ev.snapshot()['native_program_invocations']==0,'Exact selected start uses zero calls')
    check(ev.snapshot()['inherited_native_cells_indexed']==672,'All672 complete prior native cells indexed')
    for gid in e.CELLS:
        check(all(m[gid][k]==prior['selected']['model'][gid][k]for k in ('energy','forces','native')),'Exact selected E/F '+gid)
        p=m[gid]['provenance'];check(p['native_authority']==prior['selected']['model'][gid]['provenance']['native_authority'] and p['inherited_native'] and p['native_implementation']==cfg['prior_evaluator'],'Prior authority retained '+gid)
    oldpoint=prior['vectors'][prior['iterations'][-1]['derivative_keys'][0]]['parameters']
    hybrid=oldpoint[:2]+current[2:]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4)as pool:models=list(pool.map(ev.evaluate,[current,hybrid,current,hybrid]))
    check(models[0]==models[2] and models[1]==models[3],'Repeated/concurrent saved requests exactly replay')
    for gid in e.CELLS:
        if gid.endswith('_Dy'):check(models[1][gid]['provenance']['native_authority']==m[gid]['provenance']['native_authority'],'Opposite-metal exact reuse '+gid)
    baseline=ev.evaluate(list(e.INITIAL))
    check(all(v['provenance']['reused_baseline'] and v['provenance']['native_implementation']is None and v['provenance']['request_cache_implementation']==cfg['implementation']for v in baseline.values()),'Old17578 baseline authority not mislabeled as new native source')
    check(ev.snapshot()['native_program_invocations']==0 and not actual_invokes,'All real saved controls remain zero native calls')
    # Budget guards below use explicitly synthetic counters and a forbidden native stub.
    novel=current[:];novel[0]+=.000001234567
    original_keys=copy.deepcopy(ev._ledger['vector_keys']);ev._ledger['vector_keys']=['synthetic']*385
    try:ev.evaluate(novel)
    except ValueError as exc:check('Parameter-vector ceiling' in str(exc),'385-total-vector guard rejects a new request')
    else:raise AssertionError('Vector guard absent')
    ev._ledger['vector_keys']=original_keys;ev._ledger['native_program_invocations']=4096;ev._gate=lambda:None
    os.environ['SLURM_JOB_ID']='saved_control_no_scheduler_job'
    try:ev.evaluate(novel)
    except ValueError as exc:check('Native call ceiling' in str(exc),'4096-program guard rejects before invoke')
    else:raise AssertionError('Native call guard absent')
    check(ev.stop_event.is_set() and ev._ledger['stopped'],'Program-cap failure sets durable global stop')
    check(set(ev._ledger['requests'][-1]['cells'])==set(e.CELLS),'Failed request retains full eight-cell denominator')
    check(not actual_invokes,'Forbidden native stub was never reached')
    for p in (cfg['prior_result'],prior['config'],prior['evaluator']['ledger']):e.verify(p)
    check(True,'Original result/config/ledger pins unchanged')
    inherited=HERE.parent/'paired_water_radial_repair_v1'
    previous_controls=Path(cfg['training_result']['path']).parents[1]/'evaluator_validation_v1/RESULT.json'
    old=e.read(e.pin(previous_controls));check(old['status']=='passed_saved_only_evaluator_controls' and old['passed_checks']==88 and old['new_native_calls']==0,'Prior88 controls retained as inherited evidence')
    result={'status':'passed','implementation':e.pin(__file__),'evaluator':e.pin(HERE/'EVALUATE.py'),'config':e.pin(HERE/'CONFIG_EVALUATOR.json'),
        'prior_result':cfg['prior_result'],'prior_controls':e.pin(previous_controls),'records':checks,'new_native_calls':0,
        'synthetic_counter_controls_separate_from_actual_cost':True,'maximum_local_workers':4,
        'elapsed_seconds':time.monotonic()-start,'reserved_validation_opened':False}
    e.save(out/'RESULT.json',result);print(json.dumps({'status':'passed','checks':len(checks),'result':e.pin(out/'RESULT.json')}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args();main(a.out)
