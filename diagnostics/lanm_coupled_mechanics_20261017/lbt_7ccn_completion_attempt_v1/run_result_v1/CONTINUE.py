#!/usr/bin/env python3
"""One exact step92-to172/292 continuation, with no automatic follow-up."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import socket
import sys
import time
import traceback

PROTOCOL='lbt_7ccn_completion_attempt_v1'
STATE=('7CCN','La')
LIMITS={'accepted_steps':172,'evaluations':292}
ADDITIONAL={'accepted_steps':80,'evaluations':200}
FIXED={
    'previous_implementation':'0c14a03e0c3d13884739d8be33a6c2f2c39f734c5374a36de0f4d511d7573b8e',
    'previous_config':'7127a7fcba909e71a0ec3675a89159eecf627e23ce61def59e8aebc6935a14f7',
    'previous_result':'c445cf4b622fa2dbc7d7c47456e4b20b7c8db1a0b7718365fab2310c3c90a8c6',
    'checkpoint':'317c2b24eff75a2af4a6fb7066f1754d1df1b9a284405a8aad2ea3b20d5c74cd',
    'point':'eb96df292085e05f3eed1788e506c54babab646c372698a41443705945fa2a1b',
    'reference_result':'5908cb2d3b3fcb07b8ac764e3ddfb78fbc56beb29196790388cf79213af4ad72',
    'optimizer':'c5092f2480fbf393d6ca0d67733ebf2131580c922a62a67925dea8945eb5acb9',
}
CANONICAL_CHECKPOINT='981ec084b57fa86599125619b9ecce7613c883597b23e38d3521be1eaa5c0da2'


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve()
    with p.open('rb') as h: sha=hashlib.file_digest(h,'sha256').hexdigest()
    return {'path':str(p),'sha256':sha}


def verify(record):
    require(pin(record['path'])['sha256']==record['sha256'],'Changed pin: '+record['path'])
    return Path(record['path'])


def read(record): return json.loads(verify(record).read_text())


def write(path,value):
    with Path(path).open('x') as h: json.dump(value,h,indent=2,allow_nan=False); h.write('\n')
    return pin(path)


def module(name,record):
    spec=importlib.util.spec_from_file_location(name,verify(record))
    obj=importlib.util.module_from_spec(spec); spec.loader.exec_module(obj); return obj


def binding(cfg):
    return {k:(v['sha256'] if k=='implementation' else v) for k,v in cfg.items()
            if k not in ('preflight_result','run_decision')}


def prepare(cfg,out,result):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    require(cfg['protocol']==PROTOCOL and cfg['state']==list(STATE) and cfg['cumulative_limits']==LIMITS
            and cfg['additional_limits']==ADDITIONAL and pin(__file__)['sha256']==cfg['implementation']['sha256'],
            'Declared single completion scope differs')
    for k,sha in FIXED.items(): require(cfg[k]['sha256']==sha,'Frozen predecessor differs: '+k)
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}: verify(v)
    old=module('completion_previous_continuation',cfg['previous_implementation']); pcfg=read(cfg['previous_config'])
    require(cfg['python']==pcfg['python'] and cfg['python_prefix']==pcfg['python_prefix']
            and cfg['qualified_host']==pcfg['qualified_host'] and cfg['plan']==pcfg['plan'],'Original runtime/source binding differs')
    # This unchanged helper admits the original80/80 fixture and seeds its entire cache, without evaluating.
    inherited={}
    cart,loaded,_,ev,initial,oldterminal,oldcp,oldevidence=old.prepare(pcfg,out,inherited)
    previous=read(cfg['previous_result']); row=previous['rows'][0]; search=read(row['search']); cp=read(cfg['checkpoint'])
    pre=read(previous['preflight_result']); preconfig=read(pre['config']); decision=read(previous['run_decision'])
    require(previous['status']=='complete' and previous['stage']=='run' and previous['slurm_job_id']=='15254'
            and previous['implementation']==cfg['previous_implementation'] and previous['config']==cfg['previous_config']
            and previous['scientific_binding']==old.binding(pcfg)==old.binding(preconfig)
            and pre['status']=='preflight_passed' and pre['stage']=='preflight'
            and previous['preflight_result']==pcfg['preflight_result'] and previous['run_decision']==pcfg['run_decision']
            and pre['preflight_evidence']==previous['preflight_evidence']==oldevidence
            and decision['status']=='approved_single_continuation_scout' and decision['scientific_binding']==old.binding(pcfg)
            and decision['preflight_result']==previous['preflight_result'], 'Actual predecessor execution chain differs')
    require(previous['declared_states']==[list(STATE)] and len(previous['rows'])==1
            and (row['source_id'],row['metal'])==STATE and row['status']=='complete' and row['physical_state_admitted']
            and not row['converged'] and row['termination']=='accepted_step_limit'
            and row['checkpoint']==search['checkpoint']==cfg['checkpoint']
            and row['terminal_point']==search['terminal_point']==cp['point']==cfg['point']
            and row['additional_accepted_steps']==row['additional_distinct_evaluations']==12
            and search['accepted_steps']==cp['accepted_steps']==search['distinct_evaluations']==cp['distinct_evaluations']==92
            and cp['limits']=={'accepted_steps':92,'evaluations':120} and cp['initial_point']==oldcp['initial_point']
            and cp['failure'] is None and hashlib.sha256(json.dumps(cp,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()==CANONICAL_CHECKPOINT
            and search['original_source']==oldevidence['original_source'] and not search['geometry_source_rebased'],
            'Exact92/92 checkpoint and original source/start required')
    require(previous['new_composite_evaluations']==previous['new_MACE_calls']==previous['new_ALPB_calls']==12
            and all(previous[k]==0 for k in ('new_DFT_calls','new_GFN_calls','new_MM_calls','opposite_metal_calls','other_searches')),
            'Predecessor call scope differs')
    original=read(pcfg['original_result'])['rows'][0]
    require(len(original['requests'])==80 and len(previous['requests'])==12
            and cp['evaluation_requests']==original['requests']+previous['requests']
            and all(r['status']=='complete' for r in cp['evaluation_requests']), 'Cumulative92 request ledger differs')
    records=[]
    for record in original['points']+previous['points']:
        if record not in records: records.append(record)
    require(len(records)==113 and all(p in records for p in cp['evaluated_points'])
            and all(r['point'] in records for r in cp['evaluation_requests']), 'Complete historical point union differs')
    cache={}
    for record in records:
        p,_,_=ev.base.admit_point(ev.world,ev.source,ev.metal,record,True)
        cache[p['cache_key']]=record
    require(len(cache)==109 and set(ev.engine.cache)<=set(cache),'Prior coordinate cache coverage differs')
    # Keep immutable archived authorities, excluding fresh aliases produced by the inherited saved loader.
    ev.engine.cache=cache; ev.engine.records=records
    terminal=ev.read_point(cfg['point']); core=module('completion_budget_optimizer',cfg['optimizer'])
    import numpy as np
    require(core.SETTINGS==cp['settings']==oldcp['settings']==loaded['core'].SETTINGS
            and terminal['physical_state']==oldterminal['physical_state'],'Original numerical settings/state differ')
    for point in (initial,oldterminal,terminal):
        key=ev.engine.key(point['positions_A']); require(key in cache,'Saved endpoint missing from merged cache')
        restored=ev.read_point(cache[key])
        require(np.array_equal(restored['positions_A'],point['positions_A'])
                and np.array_equal(restored['gradient_kcal_mol_A'],point['gradient_kcal_mol_A'])
                and restored['native_components']==point['native_components'],'Exact endpoint cache replay differs')
        ev.gradients(point['provenance'])
    require(ev.difference(terminal,oldterminal)==row['work_from_checkpoint'],'Saved incremental work differs')
    guard=loaded['guards'][STATE].check(terminal['positions_A'])
    require(guard==row['final_geometry']==search['final_geometry'] and guard['pass'] and guard['state_compatible'],
            'Original-source terminal guard replay differs')
    ref=read(cfg['reference_result']); science=read(cfg['science_decision'])
    require(ref['status']=='complete' and ref['stage']=='collect' and ref['slurm_job_id']=='15456' and ref['full_two_available']
            and ref['input_evidence']['new_point']==cfg['point'] and ref['comparison']['status']=='complete'
            and science['protocol']==PROTOCOL and science['status']=='declared_single_finite_model_continuation'
            and science['previous_result']==cfg['previous_result'] and science['checkpoint']==cfg['checkpoint']
            and science['point']==cfg['point'] and science['incremental_reference_comparison']==cfg['reference_result']
            and science['source_state']==list(STATE) and science['additional_limits']==ADDITIONAL and science['cumulative_limits']==LIMITS,
            'Actual reference rationale/declaration differs')
    require(ev.engine.native is None and not ev.evaluation_requests and all(v==0 for v in ev.engine.cost.values()),
            'Saved preparation must make zero calls and initialize no model')
    evidence=dict(status='passed',previous_result=cfg['previous_result'],checkpoint=cfg['checkpoint'],initial_point=cp['initial_point'],
        terminal_point=cfg['point'],original_source=oldevidence['original_source'],source_guard_rebased=False,
        settings=core.SETTINGS,previous_limits=cp['limits'],cumulative_limits=LIMITS,additional_limits=ADDITIONAL,
        prior_accepted_steps=92,prior_distinct_evaluations=92,history_pairs=len(cp['history']),trust_A=cp['trust_A'],terminal_guard=guard,
        imported_point_records=113,imported_unique_keys=109,imported_request_records=92,
        cache_rows=[{'point':r,'cache_key':read(r)['cache_key'],'cache_authority':cache[read(r)['cache_key']]} for r in records],
        recipe=terminal['physical_state'],versions=oldevidence['versions'],geometry_helper=oldevidence['geometry_helper'],
        original_preflight_result=previous['preflight_result'],qualification_result=oldevidence['qualification_result'],
        no_model_initialized=True,new_molecular_calls=0)
    result.update(preflight_evidence=evidence,rationale={'incremental_reference':cfg['reference_result'],
                  'comparison':ref['comparison'],'accuracy_gate':False,'force_accuracy_qualified':False})
    return loaded,core,ev,initial,terminal,cp,evidence


def run_gate(cfg,evidence):
    require(cfg.get('preflight_result') and cfg.get('run_decision'),'Actual preflight and root run decision required')
    p=read(cfg['preflight_result']); pcfg=read(p['config']); d=read(cfg['run_decision'])
    require(p['status']=='preflight_passed' and p['stage']=='preflight' and p['protocol']==PROTOCOL
            and p['implementation']['sha256']==cfg['implementation']['sha256'] and binding(pcfg)==binding(cfg)==p['scientific_binding']
            and p['preflight_evidence']==evidence
            and p['new_composite_evaluations']==p['new_MACE_calls']==p['new_ALPB_calls']==0,'Exact actual zero-call preflight required')
    require(d['protocol']==PROTOCOL and d['status']=='approved_single_completion_attempt'
            and d['scientific_binding']==binding(cfg) and d['preflight_result']==cfg['preflight_result']
            and d['checkpoint']==cfg['checkpoint'] and d['state']==list(STATE) and d['cumulative_limits']==LIMITS,
            'Root finite completion decision differs')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==12
            and os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles'),'12CPUs and allocated GPU required')


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--stage',choices=('preflight','run'),required=True)
    parser.add_argument('--config',required=True); parser.add_argument('--output',required=True); args=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=False); tick=time.monotonic(); ev=None
    result=dict(protocol=PROTOCOL,stage=args.stage,status='unavailable',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),declared_states=[list(STATE)],
        rows=[dict(source_id='7CCN',metal='La',status='unavailable',reason='not_attempted')],new_composite_evaluations=0,
        new_MACE_calls=0,new_ALPB_calls=0,new_DFT_calls=0,new_GFN_calls=0,new_MM_calls=0,opposite_metal_calls=0,other_searches=0,
        physical_state_qualified=False,force_accuracy_qualified=False)
    try:
        cfg=read(result['config']); result['scientific_binding']=binding(cfg)
        loaded,core,ev,initial,terminal,cp,evidence=prepare(cfg,out,result)
        if args.stage=='preflight':
            require(cfg.get('preflight_result') is None and cfg.get('run_decision') is None,'No future success bindings in preflight')
            result.update(status='preflight_passed',rows=[dict(source_id='7CCN',metal='La',status='saved_inputs_admitted')])
        else:
            run_gate(cfg,evidence); result.update(preflight_result=cfg['preflight_result'],run_decision=cfg['run_decision'])
            write(out/'START.json',result)
            search=core.search(ev,loaded['guards'][STATE],initial,out/'search',LIMITS,loaded['world']['lbfgs'],resume=cp)
            search.update(source_id='7CCN',metal='La',protocol_id=PROTOCOL,previous_result=cfg['previous_result'],
                original_source=evidence['original_source'],geometry_source_rebased=False,continuation_checkpoint=cfg['checkpoint'])
            sp=write(out/'SEARCH.json',search); final=ev.read_point(search['terminal_point'])
            result['rows']=[dict(source_id='7CCN',metal='La',status=search['status'],physical_state_admitted=search['physical_state_admitted'],
                termination=search['termination'],converged=search['converged'],search=sp,terminal_point=search['terminal_point'],
                checkpoint=search['checkpoint'],additional_accepted_steps=search['accepted_steps']-92,
                additional_distinct_evaluations=search['distinct_evaluations']-92,cumulative_accepted_steps=search['accepted_steps'],
                cumulative_distinct_evaluations=search['distinct_evaluations'],work_from_checkpoint=ev.difference(final,terminal),
                work_from_original_start=search['work_from_start'],final_full_force=search['final_full_force'],final_geometry=search['final_geometry'])]
            require(92<=search['accepted_steps']<=172 and search['distinct_evaluations']==92+ev.new_evaluations<=292,
                    'Declared additional/cumulative allowance exceeded')
            result['status']=search['status']
        verify(result['config'])
    except Exception as exc: result.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
    finally:
        if ev is not None:
            cost=ev.engine.cost
            result.update(actual_costs=cost,points=ev.engine.records,requests=ev.engine.requests,model_initialized=ev.engine.native is not None,
                new_composite_evaluations=cost['new_composite_evaluations'],new_MACE_calls=cost['MACE_calls_attempted'],new_ALPB_calls=cost['ALPB_calls_attempted'])
            if any(result[k]>200 for k in ('new_composite_evaluations','new_MACE_calls','new_ALPB_calls')):
                result.update(status='unavailable',failure='Declared200-evaluation ceiling exceeded')
    result.update(elapsed_seconds=time.monotonic()-tick,allocation={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','CUDA_VISIBLE_DEVICES')},
        scientific_accuracy_claim=False,no_automatic_followup=True,other_state_Dy_force_status='unchanged unavailable in reference cross-scores')
    write(out/'RESULT.json',result); return 0 if result['status'] in ('complete','preflight_passed') else 1


if __name__=='__main__': raise SystemExit(main())
