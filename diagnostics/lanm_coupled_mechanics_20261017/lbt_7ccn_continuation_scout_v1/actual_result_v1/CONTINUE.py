#!/usr/bin/env python3
"""One exact 7CCN-La checkpoint continuation; root owns execution decisions."""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import socket
import sys
import time
import traceback

PROTOCOL = 'lbt_7ccn_La_checkpoint_continuation_scout_v1'
STATE = ('7CCN','La')
LIMITS = {'accepted_steps':92,'evaluations':120}
FIXED = {
    'original_implementation':'f2419e88d3c43d12f9b476bb2b5175308657f723c70a3e52eb08398f064e8f4b',
    'original_optimizer':'d93821da485f8d078e42bf7031018931632858ea720af617f26da2dab77617ee',
    'optimizer':'3badac29106412f8b29e79e874678c330af620fe5e6d4cc548dd0c3b3f0abfbc',
    'original_result':'fa72a4df6c7504aec669cda631366f64cc2fc82d32ba6dd30eee56fcb77fac0b',
    'checkpoint':'4e3a612380f93e9ba7f374636ed883a659c8b0ade6b35e2f7099f5d7643e9ecf',
    'point':'cf59ef7add95948eb46e3b61f86014e473e0f04035d7a06c15885163e4042022',
    'reference_result':'fff95c2f235947acb24e91991915b320731992f48c2b3d6c9ca3912a6e927f38',
    'direction_result':'cc102d11faeb3d7ee7afccd9f772b72a99d2866448fa42e46223f3bbe53bc987',
}
CHECKPOINT_CANONICAL = '2be9f7cd5ca1da78aebec332268a30a26a350b659af4f7d567ad00cb22847436'


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
    obj=importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def binding(cfg):
    return {k:(v['sha256'] if k=='implementation' else v) for k,v in cfg.items()
            if k not in ('preflight_result','run_decision')}


def prepare(cfg,out,result):
    require(cfg['protocol']==PROTOCOL and cfg['state']==list(STATE) and cfg['cumulative_limits']==LIMITS
            and cfg['additional_limits']=={'accepted_steps':12,'evaluations':40},'Declared single-scout scope differs')
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Continuation source differs')
    for k,sha in FIXED.items(): require(cfg[k]['sha256']==sha,'Frozen predecessor differs: '+k)
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}: verify(v)
    require(pin(sys.executable)['sha256']==cfg['python']['sha256']
            and str(Path(sys.prefix).resolve())==cfg['python_prefix'],'Original Python environment required')
    require(socket.gethostname().split('.')[0]==cfg['qualified_host']=='node-224-2t-8gpu-1',
            'Qualified historical-asset host required')
    original=read(cfg['original_result']); ocfg=read(original['config']); row=original['rows'][0]
    search=read(row['search']); cp=read(cfg['checkpoint']); p=read(cfg['point'])
    require(original['status']=='complete' and original['stage']=='run' and original['slurm_job_id']=='12968'
            and original['implementation']==ocfg['implementation']==cfg['original_implementation']
            and original['optimizer']==ocfg['optimizer']==cfg['original_optimizer']
            and len(original['rows'])==1 and (row['source_id'],row['metal'])==STATE
            and row['status']=='complete' and row['physical_state_admitted'],'Original completed search required')
    require(search['status']=='complete' and search['physical_state_admitted'] and not search['converged']
            and search['termination']=='accepted_step_limit' and search['checkpoint']==cfg['checkpoint']
            and search['terminal_point']==row['terminal_point']==cp['point']==cfg['point']
            and search['start_point']==cp['initial_point'] and search['accepted_steps']==cp['accepted_steps']==80
            and search['distinct_evaluations']==cp['distinct_evaluations']==80 and cp['failure'] is None
            and search['limits']==cp['limits']=={'accepted_steps':80,'evaluations':200}
            and hashlib.sha256(json.dumps(cp,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()==CHECKPOINT_CANONICAL,
            'Exact80/80 terminal checkpoint required')
    reference=read(cfg['reference_result']); audit=read(cfg['direction_result'])
    require(reference['status']=='complete' and reference['slurm_job_id']=='15173' and reference['full_four_available']
            and audit['status']=='complete' and audit['slurm_job_id']=='15176'
            and audit['inputs']['collector']==cfg['reference_result'] and audit['inputs']['checkpoint']==cfg['checkpoint']
            and audit['inputs']['point']==cfg['point'] and audit['molecular_calls']==audit['optimizer_steps']==0,
            'Actual scientific rationale records required')
    result['rationale']={'reference_comparison':reference['comparison'],'conditional_direction':audit['direction']['derivatives'],
        'accuracy_gate':False,'force_accuracy_qualified':False,'caveats':reference['caveats']}
    cart=module('single_continuation_original_cartesian',cfg['original_implementation'])
    loaded=cart.load_inputs(ocfg,out)
    core=module('single_continuation_budget_optimizer',cfg['optimizer'])
    import numpy as np
    require(core.SETTINGS==loaded['core'].SETTINGS==search['settings']==cp['settings']
            and loaded['geometry_helper']==original['geometry_helper'],'Unchanged optimizer/direction helper required')
    qualification=cart.parent_gate(ocfg,'qualification_result','qualify')
    require(qualification['qualification_passed'] and qualification['declared_directions']==8
            and qualification['declared_cells']==32
            and all(r['status']=='passed' for r in qualification['qualification_rows']), 'Original full derivative gate required')
    scout=cart.parent_gate(ocfg,'scout_result','scout')
    require(scout['continuation_execution_and_physical_admission'],'Original scout admission differs')
    require(loaded['terminals'][STATE]['point_pin']==cp['initial_point']
            and loaded['terminals'][STATE]['source']['prepared']['physical_atoms']==search['original_source']
            and search['geometry_source_rebased'] is False,'Original start/source guard cannot be rebased')
    versions={k:importlib.metadata.version(k) for k in ('torch','mace-torch','ase')}
    require(versions==loaded['world']['mcfg']['expected_versions']==p['recipe']['MACE_versions'],'Original native versions differ')
    require(len(row['points'])==99 and len(row['requests'])==80
            and cp['evaluation_requests']==row['requests']
            and all(q['status']=='complete' and q['point'] in row['points'] for q in row['requests'])
            and all(q in row['points'] for q in cp['evaluated_points']),'Original cache/request denominator differs')
    ev=cart.Evaluator(loaded,STATE,out/'evaluations',additional_points=row['points'])
    initial=ev.read_point(cp['initial_point']); terminal=ev.read_point(cp['point'])
    require(terminal['physical_state']==p['recipe'] and p['source_id']=='7CCN' and p['metal']=='La'
            and p['recipe']['physical_charge']==-2 and p['recipe']['physical_multiplicity']==1
            and p['recipe']['target_index']==255,'Exact whole physical state differs')
    cache_rows=[]
    for record in row['points']:
        old=read(record); key=old['cache_key']
        require(key in ev.engine.cache and ev.engine.cache[key] is not None,'Previously admitted cache entry unavailable')
        cache_rows.append({'point':record,'cache_key':key,'cache_authority':ev.engine.cache[key]})
    for point in (initial,terminal):
        key=ev.engine.key(point['positions_A'])
        require(key in ev.engine.cache and ev.engine.cache[key] is not None,'Initial/terminal cache identity missing')
        restored=ev.read_point(ev.engine.cache[key])
        require(np.array_equal(restored['positions_A'],point['positions_A'])
                and np.array_equal(restored['gradient_kcal_mol_A'],point['gradient_kcal_mol_A'])
                and restored['native_components']==point['native_components'],'Exact endpoint cache replay differs')
        ev.gradients(point['provenance'])
    guard=loaded['guards'][STATE]; initial_guard=guard.check(initial['positions_A']); terminal_guard=guard.check(terminal['positions_A'])
    require(all(a['pass'] and a['state_compatible'] for a in (initial_guard,terminal_guard)), 'Original-source physical admission failed')
    require(ev.engine.native is None and not ev.evaluation_requests and all(v==0 for v in ev.engine.cost.values()),
            'Saved preflight must not initialize a model or call any component')
    evidence=dict(status='passed',initial_point=cp['initial_point'],terminal_point=cp['point'],checkpoint=cfg['checkpoint'],
        original_source=search['original_source'],source_guard_rebased=False,settings=core.SETTINGS,
        cumulative_limits=LIMITS,prior_accepted_steps=80,prior_distinct_evaluations=80,
        history_pairs=len(cp['history']),trust_A=cp['trust_A'],initial_guard=initial_guard,terminal_guard=terminal_guard,
        cache_rows=cache_rows,imported_point_records=len(row['points']),imported_unique_keys=len({r['cache_key'] for r in cache_rows}),
        recipe=p['recipe'],versions=versions,geometry_helper=loaded['geometry_helper'],
        qualification_result=ocfg['qualification_result'],original_scout_result=ocfg['scout_result'],
        no_model_initialized=True,new_molecular_calls=0)
    result['preflight_evidence']=evidence
    return cart,loaded,core,ev,initial,terminal,cp,evidence


def run_gate(cfg,evidence):
    require(cfg.get('preflight_result') and cfg.get('run_decision'),'Actual zero-call preflight and root review required')
    pre=read(cfg['preflight_result']); pcfg=read(pre['config']); decision=read(cfg['run_decision'])
    require(pre['status']=='preflight_passed' and pre['stage']=='preflight' and pre['protocol']==PROTOCOL
            and binding(pcfg)==binding(cfg) and pre['scientific_binding']==binding(cfg)
            and pre['preflight_evidence']==evidence and pre['new_composite_evaluations']==pre['new_MACE_calls']==pre['new_ALPB_calls']==0,
            'Actual exact zero-call preflight required')
    require(decision['protocol']==PROTOCOL and decision['status']=='approved_single_continuation_scout'
            and decision['scientific_binding']==binding(cfg) and decision['preflight_result']==cfg['preflight_result']
            and decision['checkpoint']==cfg['checkpoint'] and decision['state']==list(STATE)
            and decision['cumulative_limits']==LIMITS,'Root single-scout execution review differs')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==12
            and os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles'), '12CPUs and allocated GPU required')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('preflight','run'),required=True)
    parser.add_argument('--config',required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=1,'Allocated CPU required')
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();ev=None
    result=dict(protocol=PROTOCOL,stage=args.stage,status='unavailable',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),declared_states=[list(STATE)],
        rows=[{'source_id':'7CCN','metal':'La','status':'unavailable','reason':'not_attempted'}],
        new_composite_evaluations=0,new_MACE_calls=0,new_ALPB_calls=0,new_DFT_calls=0,new_GFN_calls=0,new_MM_calls=0,
        opposite_metal_calls=0,other_searches=0,physical_state_qualified=False,force_accuracy_qualified=False)
    try:
        cfg=read(result['config']);result['scientific_binding']=binding(cfg)
        cart,loaded,core,ev,initial,terminal,cp,evidence=prepare(cfg,out,result)
        if args.stage=='preflight':
            require(cfg.get('preflight_result') is None and cfg.get('run_decision') is None,'Preflight must be unbound to future success')
            result.update(status='preflight_passed',rows=[dict(source_id='7CCN',metal='La',status='saved_inputs_admitted')])
        else:
            run_gate(cfg,evidence);result.update(preflight_result=cfg['preflight_result'],run_decision=cfg['run_decision'])
            write(out/'START.json',result)
            search=core.search(ev,loaded['guards'][STATE],initial,out/'search',LIMITS,loaded['world']['lbfgs'],resume=cp)
            search.update(source_id='7CCN',metal='La',protocol_id=PROTOCOL,original_search=read(cfg['original_result'])['rows'][0]['search'],
                original_source=evidence['original_source'],geometry_source_rebased=False,continuation_checkpoint=cfg['checkpoint'])
            sp=write(out/'SEARCH.json',search);final=ev.read_point(search['terminal_point'])
            result['rows']=[dict(source_id='7CCN',metal='La',status=search['status'],physical_state_admitted=search['physical_state_admitted'],
                termination=search['termination'],converged=search['converged'],search=sp,terminal_point=search['terminal_point'],
                checkpoint=search['checkpoint'],additional_accepted_steps=search['accepted_steps']-80,
                additional_distinct_evaluations=search['distinct_evaluations']-80,
                cumulative_accepted_steps=search['accepted_steps'],cumulative_distinct_evaluations=search['distinct_evaluations'],
                work_from_checkpoint=ev.difference(final,terminal),work_from_original_start=search['work_from_start'],
                final_full_force=search['final_full_force'],final_geometry=search['final_geometry'])]
            require(80<=search['accepted_steps']<=92 and search['distinct_evaluations']==80+ev.new_evaluations<=120,
                    'Declared cumulative/additional allowance exceeded')
            result['status']=search['status']
        verify(result['config'])
    except Exception as exc:result.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
    finally:
        if ev is not None:
            cost=ev.engine.cost
            result.update(actual_costs=cost,points=ev.engine.records,requests=ev.engine.requests,
                new_composite_evaluations=cost['new_composite_evaluations'],new_MACE_calls=cost['MACE_calls_attempted'],
                new_ALPB_calls=cost['ALPB_calls_attempted'],model_initialized=ev.engine.native is not None)
            if any(result[k]>40 for k in ('new_composite_evaluations','new_MACE_calls','new_ALPB_calls')):
                result.update(status='unavailable',failure='Declared40-call ceiling exceeded')
    result.update(elapsed_seconds=time.monotonic()-tick,allocation={k:os.environ.get(k) for k in
                  ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','CUDA_VISIBLE_DEVICES')},
                  scientific_accuracy_claim=False,other_state_Dy_force_status='unchanged unavailable in reference cross-scores')
    write(out/'RESULT.json',result)
    return 0 if result['status'] in ('complete','preflight_passed') else 1


if __name__=='__main__':raise SystemExit(main())
