#!/usr/bin/env python3
"""One exact step172 continuation; completed saved audit gates lazy native startup."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import socket
import time
import traceback

PROTOCOL='lbt_7ccn_generous_completion_v1'
STATE=('7CCN','La')
LIMITS={'accepted_steps':1172,'evaluations':2172}
ADDITIONAL={'accepted_steps':1000,'evaluations':2000}
FIXED={
    'previous_implementation':'551d964badc41141cfc53a85d704135ec74c21f185e12b6cc50d578dc47910d2',
    'previous_config':'8a78864e5baa4836c1fef00240023cdc941493f41a9dad2e2a47aef31ea9a997',
    'previous_result':'9159ea43cf76692fb9ad40acbf199c1a6b146acdf8dc21ff68a045350a764a10',
    'previous_search':'4d59031becf63f38381a33733bbda62fef4ebe48e8ba0cceef8266bd627304ca',
    'checkpoint':'14cd220d924040caf046060c0b6fe3869941bff4bfeb9b210cc1a1e74a5eeacc',
    'point':'67690a33f7233a3e1b17b1abb687d41717fefd5903e3c04199a5060b1d6b7d81',
    'saved_audit':'228ef78599c1ecf1bdd67120c32d5b83cc9f861cf22b47e4bef07f7436c1f1fb',
    'science_decision':'e8c4bda02f6714f452251b2c4edd83dbbdebfdbafa81b84cfdef75cf5574cc80',
    'optimizer':'6635c0f062ec9f41cd5bc84e12606af2490bd34a5ed216843d1735ea7b017db2',
}
CANONICAL_CHECKPOINT='2d9a27465af93f883099a3e9d0562af1ffe98f41ad8d590587cb10dce3a06535'


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


def prepare(cfg,out,result):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==12
            and os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles'),
            '12CPUs and allocated GPU required')
    require(cfg['protocol']==PROTOCOL and cfg['state']==list(STATE) and cfg['cumulative_limits']==LIMITS
            and cfg['additional_limits']==ADDITIONAL and pin(__file__)['sha256']==cfg['implementation']['sha256'],
            'Declared single completion scope differs')
    for k,sha in FIXED.items(): require(cfg[k]['sha256']==sha,'Frozen predecessor differs: '+k)
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}: verify(v)
    old=module('generous_previous_completion',cfg['previous_implementation']); pcfg=read(cfg['previous_config'])
    require(cfg['python']==pcfg['python'] and cfg['python_prefix']==pcfg['python_prefix']
            and cfg['qualified_host']==pcfg['qualified_host'] and cfg['plan']==pcfg['plan'],
            'Original runtime/source binding differs')
    # The inherited loader creates one lazy evaluator and admits the full step92 chain.
    inherited={}
    loaded,_,ev,initial,oldterminal,oldcp,oldevidence=old.prepare(pcfg,out,inherited)
    old.run_gate(pcfg,oldevidence)
    previous=read(cfg['previous_result']); require(len(previous['rows'])==1,'Single predecessor row required')
    row=previous['rows'][0]; search=read(cfg['previous_search']); cp=read(cfg['checkpoint'])
    require(previous['status']=='complete' and previous['stage']=='run' and previous['slurm_job_id']=='15458'
            and previous['implementation']==cfg['previous_implementation'] and previous['config']==cfg['previous_config']
            and previous['scientific_binding']==old.binding(pcfg) and previous['preflight_evidence']==oldevidence
            and previous['preflight_result']==pcfg['preflight_result'] and previous['run_decision']==pcfg['run_decision'],
            'Actual predecessor execution chain differs')
    require(previous['declared_states']==[list(STATE)] and (row['source_id'],row['metal'])==STATE
            and row['status']=='complete' and row['physical_state_admitted'] and not row['converged']
            and row['termination']=='accepted_step_limit' and row['search']==cfg['previous_search']
            and row['checkpoint']==search['checkpoint']==cfg['checkpoint']
            and row['terminal_point']==search['terminal_point']==cp['point']==cfg['point']
            and row['additional_accepted_steps']==row['additional_distinct_evaluations']==80
            and search['accepted_steps']==cp['accepted_steps']==search['distinct_evaluations']==cp['distinct_evaluations']==172
            and cp['limits']=={'accepted_steps':172,'evaluations':292} and cp['initial_point']==oldcp['initial_point']
            and cp['failure'] is None and hashlib.sha256(json.dumps(cp,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()==CANONICAL_CHECKPOINT
            and search['original_source']==oldevidence['original_source'] and not search['geometry_source_rebased'],
            'Exact172/172 checkpoint and original source/start required')
    require(previous['new_composite_evaluations']==previous['new_MACE_calls']==previous['new_ALPB_calls']==80
            and all(previous[k]==0 for k in ('new_DFT_calls','new_GFN_calls','new_MM_calls','opposite_metal_calls','other_searches')),
            'Predecessor call scope differs')
    require(len(oldcp['evaluation_requests'])==92 and len(previous['requests'])==80
            and cp['evaluation_requests']==oldcp['evaluation_requests']+previous['requests']
            and all(r['status']=='complete' for r in cp['evaluation_requests']), 'Cumulative172 request ledger differs')
    records=list(previous['points'])
    require(len(records)==193 and len({(r['path'],r['sha256']) for r in records})==193
            and all(r in records for r in ev.engine.records)
            and len(cp['evaluated_points'])==172 and all(r in records for r in cp['evaluated_points'])
            and all(r['point'] in records for r in cp['evaluation_requests']), 'Complete historical point union differs')
    cache={}
    for record in records:
        p,_,_=ev.base.admit_point(ev.world,ev.source,ev.metal,record,True)
        cache[p['cache_key']]=record
    require(len(cache)==189 and set(ev.engine.cache)<=set(cache),'Prior coordinate cache coverage differs')
    ev.engine.cache=cache; ev.engine.records=records
    terminal=ev.read_point(cfg['point']); core=module('generous_budget_optimizer',cfg['optimizer'])
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
    audit=read(cfg['saved_audit']); ai=audit['input_evidence']; science=read(cfg['science_decision'])
    verify(audit['implementation']); verify(audit['config'])
    require(audit['protocol']=='lbt_7ccn_trace_audit_v1' and audit['status']=='complete'
            and audit['slurm_job_id']=='15489' and audit['replay_passed'] and audit['scientific_force_criteria_unchanged']
            and audit['declared_points']==audit['admitted_points']==len(audit['points'])==81
            and audit['declared_transitions']==audit['admitted_transitions']==len(audit['transitions'])==80
            and len(audit['replay_checks'])==2828 and all(c['pass'] for c in audit['replay_checks'])
            and not audit['model_initialized'] and not audit['force_accuracy_qualified']
            and all(audit[k]==0 for k in ('new_molecular_calls','model_loads','optimizer_steps','geometry_trials'))
            and ai['run']==cfg['previous_result'] and ai['initial_checkpoint']==pcfg['checkpoint']
            and ai['final_checkpoint']==cfg['checkpoint'] and ai['initial_point']==cp['initial_point']
            and ai['source']==oldevidence['original_source'] and ai['geometry_helper']==oldevidence['geometry_helper']
            and ai['recipe']==terminal['physical_state'] and ai['settings']==core.SETTINGS and ai['optimizer']==pcfg['optimizer'],
            'Completed saved audit15489 and exact chain required')
    require(science['protocol']==PROTOCOL and science['status']=='declared_single_generous_model_completion'
            and science['state']==list(STATE) and science['additional_limits']==ADDITIONAL and science['cumulative_limits']==LIMITS
            and science['previous_result']==cfg['previous_result'] and science['previous_search']==cfg['previous_search']
            and science['checkpoint']==cfg['checkpoint'] and science['point']==cfg['point']
            and science['previous_implementation']==cfg['previous_implementation'] and science['previous_configuration']==cfg['previous_config']
            and science['previous_optimizer']==pcfg['optimizer'] and science['saved_audit']==cfg['saved_audit']
            and science['canonical_checkpoint_sha256']==CANONICAL_CHECKPOINT,'Root single-completion declaration differs')
    require(ev.engine.native is None and not ev.evaluation_requests and all(v==0 for v in ev.engine.cost.values()),
            'Saved admission must make zero calls and initialize no model')
    evidence=dict(status='passed',previous_result=cfg['previous_result'],checkpoint=cfg['checkpoint'],initial_point=cp['initial_point'],
        terminal_point=cfg['point'],original_source=oldevidence['original_source'],source_guard_rebased=False,
        settings=core.SETTINGS,previous_limits=cp['limits'],cumulative_limits=LIMITS,additional_limits=ADDITIONAL,
        prior_accepted_steps=172,prior_distinct_evaluations=172,history_pairs=len(cp['history']),trust_A=cp['trust_A'],terminal_guard=guard,
        imported_point_records=193,imported_unique_keys=189,imported_request_records=172,
        recipe=terminal['physical_state'],versions=oldevidence['versions'],geometry_helper=oldevidence['geometry_helper'],
        saved_audit=cfg['saved_audit'],saved_audit_replay_checks=2828,previous_preflight=previous['preflight_result'],
        qualification_result=oldevidence['qualification_result'],no_model_initialized=True,new_molecular_calls=0)
    result.update(saved_admission=evidence,saved_audit=cfg['saved_audit'],science_decision=cfg['science_decision'],
        inherited_reference_rationale=inherited['rationale'],accuracy_gate=False)
    return loaded,core,ev,initial,terminal,cp,evidence


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True); parser.add_argument('--output',required=True); args=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=False); tick=time.monotonic(); ev=None
    result=dict(protocol=PROTOCOL,stage='run',status='unavailable',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),declared_states=[list(STATE)],
        rows=[dict(source_id='7CCN',metal='La',status='unavailable',reason='not_attempted')],new_composite_evaluations=0,
        new_MACE_calls=0,new_ALPB_calls=0,new_DFT_calls=0,new_GFN_calls=0,new_MM_calls=0,opposite_metal_calls=0,other_searches=0,
        model_initialized=False,physical_state_qualified=False,force_accuracy_qualified=False)
    try:
        cfg=read(result['config']); result['scientific_binding']=cfg
        loaded,core,ev,initial,terminal,cp,evidence=prepare(cfg,out,result)
        write(out/'START.json',result)
        search=core.search(ev,loaded['guards'][STATE],initial,out/'search',LIMITS,loaded['world']['lbfgs'],resume=cp)
        search.update(source_id='7CCN',metal='La',protocol_id=PROTOCOL,previous_result=cfg['previous_result'],
            original_source=evidence['original_source'],geometry_source_rebased=False,continuation_checkpoint=cfg['checkpoint'])
        sp=write(out/'SEARCH.json',search); final=ev.read_point(search['terminal_point'])
        result['rows']=[dict(source_id='7CCN',metal='La',status=search['status'],physical_state_admitted=search['physical_state_admitted'],
            termination=search['termination'],converged=search['converged'],search=sp,terminal_point=search['terminal_point'],
            checkpoint=search['checkpoint'],additional_accepted_steps=search['accepted_steps']-172,
            additional_distinct_evaluations=search['distinct_evaluations']-172,cumulative_accepted_steps=search['accepted_steps'],
            cumulative_distinct_evaluations=search['distinct_evaluations'],work_from_checkpoint=ev.difference(final,terminal),
            work_from_original_start=search['work_from_start'],final_full_force=search['final_full_force'],final_geometry=search['final_geometry'])]
        require(172<=search['accepted_steps']<=1172 and search['distinct_evaluations']==172+ev.new_evaluations<=2172,
                'Declared additional/cumulative allowance exceeded')
        result['status']=search['status']; verify(result['config'])
    except Exception as exc: result.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
    finally:
        if ev is not None:
            cost=ev.engine.cost
            result.update(actual_costs=cost,points=ev.engine.records,requests=ev.engine.requests,model_initialized=ev.engine.native is not None,
                new_composite_evaluations=cost['new_composite_evaluations'],new_MACE_calls=cost['MACE_calls_attempted'],new_ALPB_calls=cost['ALPB_calls_attempted'])
            if any(result[k]>2000 for k in ('new_composite_evaluations','new_MACE_calls','new_ALPB_calls')):
                result.update(status='unavailable',failure='Declared2000-evaluation ceiling exceeded')
    result.update(elapsed_seconds=time.monotonic()-tick,allocation={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','CUDA_VISIBLE_DEVICES')},
        scientific_accuracy_claim=False,no_automatic_followup=True,other_state_Dy_force_status='unchanged unavailable in reference cross-scores')
    write(out/'RESULT.json',result); return 0 if result['status']=='complete' else 1


if __name__=='__main__': raise SystemExit(main())
