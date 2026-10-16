#!/usr/bin/env python3
"""Paired prior/control/candidate menus for the frozen LBT3 common-basis test.

CLI: LBT3_BASIS_POOL.py run --config CONFIG --output FRESH
Config pins: implementation, plan, range_pool_result, basis_scout_result,
basis_remaining_result. Reuse the complete range pool and exact own endpoints.
At most eight cross-state composite energies, in two sequential batches using
the unchanged four-cell pool executor. No new gradients or optimization.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback

RANGE_POOL_SHA = '2abe9640b29f2335e89c05724666a2152a24e3fa9ca364b85438f25bfb3dec25'
PROTOCOL = 'lbt3_common_basis_terminal_pool_v1'
SEARCH_PROTOCOL = 'lbt3_common_basis_response_v1'
SOURCES, METALS, ARMS = ('7CCO', '7CCN'), ('La', 'Dy'), ('control', 'candidate')
EXTRA = [9, 38, 32]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def read(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: ' + record['path'])
    return json.loads(Path(record['path']).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False); handle.write('\n')


def module(name, record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed helper')
    spec = importlib.util.spec_from_file_location(name, record['path'])
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def empty_cost():
    return {'new_composite_evaluations': 0, 'MACE_calls_attempted': 0, 'GFN2_calls_attempted': 0,
            'model_load_seconds': 0., 'MACE_call_seconds': 0., 'GFN2_call_seconds': 0.}


def inputs(config_path):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    cfg = json.loads(config_path.read_text())
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Basis-pool adapter differs')
    prior = read(cfg['range_pool_result']); prior_cfg = read(prior['config'])
    require(prior['status'] == 'complete' and prior['stage'] == 'run'
            and prior['implementation'] == prior_cfg['implementation']
            and prior['implementation']['sha256'] == RANGE_POOL_SHA, 'Complete prior range pool required')
    range_pool = module('lbt3_basis_pool_frozen_range_pool', prior['implementation'])
    (_, pool, world, old, candidates, aliases, cells, old_slots, range_slots,
     range_module, base, range_parents) = range_pool.inputs(Path(prior['config']['path']))
    range_admitted, missing = range_pool.add_terminals(pool, world, candidates, aliases, cells,
                                                      old_slots, range_slots, range_module, base)
    require(all(s['status'] == 'admitted' for s in range_admitted)
            and prior['searches'] == range_admitted and prior['aliases'] == aliases
            and prior['candidates'] == [{k: v for k, v in c.items() if k != 'x'} for c in candidates]
            and prior['required_cross_cells'] == [{'candidate_id': i, 'metal': m} for i, m in missing],
            'Prior range candidate/admission replay differs')
    require(len(prior['cross_rows']) == len(missing) and
            {(r['candidate_id'], r['metal']) for r in prior['cross_rows']} == set(missing), 'Prior range matrix incomplete')
    for row in prior['cross_rows']:
        require(row['status'] == 'accepted' and read(row['cell']) == {k: v for k, v in row.items() if k != 'cell'}, 'Prior range cell differs')
        cells[(row['candidate_id'], row['metal'])] = row
    require(len(cells) == 2*len(candidates) and range_pool.summarize(pool, old, candidates, aliases, cells, range_admitted) == prior['summary'],
            'Complete prior range matrix/report replay differs')
    pool.verify(cfg['plan'])
    scout, remaining = read(cfg['basis_scout_result']), read(cfg['basis_remaining_result'])
    require(scout['stage'] == 'scout' and remaining['stage'] == 'remaining'
            and remaining['basis_scout_result'] == cfg['basis_scout_result'], 'Matched basis stage lineage differs')
    for parent in (scout, remaining):
        parent_cfg = read(parent['config']); pool.verify(parent['implementation'])
        require(parent['status'] in ('complete', 'unavailable') and parent['protocol_id'] == SEARCH_PROTOCOL
                and parent['preflight']['status'] == 'passed' and parent['preflight']['declared_arm_states'] == 8
                and parent['preflight']['molecular_calls'] == 0 and parent['preflight']['model_loaded'] is False
                and read(parent['preflight_result']) == parent['preflight'], 'Basis execution/preflight failure cannot admit arms')
        require((parent['status'] == 'complete') == all(r['status'] == 'complete' and r['candidate_admitted'] for r in parent['rows']),
                'Basis parent/row availability differs')
        require(parent['plan'] == parent_cfg['plan'] == cfg['plan'] and parent['implementation'] == parent_cfg['implementation']
                and parent['range_scout_result'] == parent_cfg['range_scout_result'] == prior['range_scout_result']
                and parent['range_remaining_result'] == parent_cfg['range_remaining_result'] == prior['range_remaining_result']
                and parent['profile_result'] == range_parents[0]['profile_result'] and parent['fd_result'] == range_parents[0]['fd_result']
                and parent['frozen_scout_implementation'] == range_parents[0]['frozen_scout_implementation']
                and parent['frozen_range_implementation'] == range_parents[0]['implementation'], 'Basis/source/recipe lineage differs')
    require(scout['implementation']['sha256'] == remaining['implementation']['sha256']
            and scout['settings'] == remaining['settings'] == range_parents[0]['settings']
            and scout['fallback_rule'] == remaining['fallback_rule'], 'Matched arm search settings differ')
    global np
    import numpy as np
    require(pin(sys.executable)['sha256'] == world['mcfg']['python']['sha256'], 'Native Python differs')
    slots = [dict(row, parent=cfg[key]) for key, parent in (('basis_scout_result', scout), ('basis_remaining_result', remaining)) for row in parent['rows']]
    require(len(slots) == 8 and {(r['source_id'], r['metal'], r['arm']) for r in slots}
            == {(s, m, a) for s in SOURCES for m in METALS for a in ARMS}, 'Eight-arm denominator differs')
    return cfg, pool, world, prior, candidates, aliases, cells, range_admitted, slots, range_module, base, (scout, remaining)


def add_terminals(pool, world, candidates, aliases, cells, prior_slots, slots, range_module, base, fallback_rule):
    from scaffold_accommodation import Geometry, SETTINGS as GEOMETRY_SETTINGS
    old_by = {(s['source_id'], s['metal']): s for s in prior_slots}
    by_coord = {c['coordinate_sha256']: c for c in candidates}; admitted = []; transforms = {}
    old_active = world['profile']['active_candidate_basis']['mode_indices_selection_order']
    for row in sorted(slots, key=lambda r: (SOURCES.index(r['source_id']), METALS.index(r['metal']), ARMS.index(r['arm']))):
        source_id, metal, arm = row['source_id'], row['metal'], row['arm']; source = world['sources'][source_id]
        previous = old_by[(source_id, metal)]; old_search = read(previous['search'])
        item = {'source_id': source_id, 'metal': metal, 'arm': arm, 'parent': row['parent'], 'status': 'unavailable'}
        if row.get('search') is None:
            item['failure'] = 'No arm search receipt'; admitted.append(item); continue
        search = read(row['search']); item['search'] = row['search']
        require((search['source_id'], search['metal'], search['arm']) == (source_id, metal, arm), 'Arm slot identity differs')
        if row['status'] != 'complete' or not search['candidate_admitted']:
            item['failure'] = 'Arm unavailable; archived points are not replacement terminal candidates'; admitted.append(item); continue
        active = old_active if arm == 'control' else old_active+EXTRA
        require(search['protocol_id'] == SEARCH_PROTOCOL and search['status'] == 'finite_proposal_available'
                and search['candidate'] in row['points'] and search['candidate'] == search['trace'][-1]['point']
                and search['continuation_start'] == search['trace'][0]['point'] == previous['point']
                and search['previous_search'] == previous['search'] and search['source_rebased'] is False
                and search['coordinate_offset_q64'] == old_search['full_q']
                and search['authoritative_origin'] == old_search['authoritative_origin']
                and search['active_indices'] == active and search['active_names'] == [source['kin'].modes[k]['name'] for k in active]
                and search['initial_history_length'] == search['initial_new_evaluations'] == 0
                and search['initial_offset_v_A'] == [0.]*len(active) and search['fallback_rule'] == fallback_rule
                and len(search['trace']) == search['accepted_steps']+1 and search['origin_guard']['pass']
                and search['continuation_start_guard']['pass'] and search['final_geometry']['pass']
                and all(t['geometry']['pass'] for t in search['trace']), 'Matched start/reference/arm admission differs')
        point = read(search['candidate']); x, components = pool.point_components(world, source, metal, point)
        q, q0, v, transform = [np.asarray(search[k], dtype=float) for k in
                              ('full_q', 'coordinate_offset_q64', 'offset_v_A', 'coordinate_scaling_T')]
        require(q.shape == q0.shape == (64,) and v.shape == (len(active),) and transform.shape == (len(active),len(active))
                and all(np.isfinite(a).all() for a in (q, q0, v, transform))
                and np.array_equal(q[[k for k in range(64) if k not in active]], q0[[k for k in range(64) if k not in active]])
                and np.allclose(q[active], q0[active]+transform@v, rtol=0, atol=1e-12), 'Arm absolute/offset coordinates differ')
        replay, _, _, _ = source['kin'].evaluate(q); _, _, jac0, _ = source['kin'].evaluate(np.zeros(64))
        matrix = jac0[active].reshape(len(active), -1).T
        require(np.allclose(replay, x, rtol=0, atol=1e-12)
                and np.allclose((matrix@transform).T@(matrix@transform), np.eye(len(active)), rtol=0, atol=1e-10), 'Arm kinematics/source whitening differs')
        key_arm = (source_id, arm)
        if key_arm in transforms:
            require(np.array_equal(transforms[key_arm], transform), 'Metals do not share the same source/arm transform')
        transforms[key_arm] = transform
        guard = range_module.range_guard(base, source, active, Geometry, GEOMETRY_SETTINGS, metal); check = guard.check(q,x,x)
        require(check['pass'] and check['source_baseline'] == source['prepared']['physical_atoms']
                and check['source_displacement_bound_A'] == .40, 'Arm original-source guard failed')
        key = hashlib.sha256(np.where(x == 0.,0.,x).astype('<f8').tobytes()).hexdigest()
        if key not in by_coord:
            candidate = {'candidate_id': 'geometry_%02d' % len(candidates), 'coordinate_sha256': key,
                'source_id': source_id, 'physical_source': point['physical_source'], 'aliases': [], 'x': x}
            candidates.append(candidate); by_coord[key] = candidate
        candidate = by_coord[key]; require(np.array_equal(candidate['x'],x), 'Coordinate hash collision')
        alias = source_id+':'+arm+'_terminal_'+metal; candidate['aliases'].append(alias); aliases[alias] = candidate['candidate_id']
        cell_key = (candidate['candidate_id'],metal)
        if cell_key in cells:
            require(cells[cell_key]['components'] == components, 'Conflicting exact-state cached energies')
        else:
            cells[cell_key] = {'candidate_id': candidate['candidate_id'], 'metal': metal, 'status': 'accepted', 'components': components,
                'reused': True, 'evidence': [{'terminal_point': search['candidate'], 'search': row['search'], 'parent': row['parent']}]}
        work = pool.difference(components,world['origin_components'][(source_id,metal)])
        increment = pool.difference(components,read(previous['point'])['components'])
        require(work == search['work_from_origin'] and increment == search['work_from_continuation_start']
                and increment['composite_kcal_mol'] <= 1e-8, 'Arm work/reference differs')
        item.update(status='admitted', candidate_id=candidate['candidate_id'], point=search['candidate'], previous_point=previous['point'],
            physical_check=check, own_work_from_original=work, own_incremental_work=increment,
            termination=search['termination'], accepted_steps=search['accepted_steps'], accepted_fallback_steps=search['accepted_fallback_steps'],
            final_force_diagnostic=search['final_force_diagnostic'], boundary_flags=search['preserved_boundary_flags'])
        admitted.append(item)
    # A saved trial may score an already admitted geometry without becoming a
    # candidate itself. Include both arms' shared old/new physical-point caches.
    seen = set(); blocked = {}
    for row in slots:
        source = world['sources'][row['source_id']]; metal = row['metal']
        for record in row.get('points', []):
            identity = (record['path'],record['sha256'])
            if identity in seen:
                continue
            seen.add(identity); point = read(record)
            require(point['source_id'] == row['source_id'] and point['metal'] == metal
                    and point['recipe'] == pool.recipe(world,source,metal), 'Saved trial recipe/state differs')
            x = pool.check_coordinates(world,source,metal,point['xyz'],point['physical_source'],point['positions'])
            require(point['cache_key'] == pool.digest({'recipe':point['recipe'],'positions_A':x.tolist()}), 'Saved trial cache identity differs')
            coordinate_key = hashlib.sha256(np.where(x == 0.,0.,x).astype('<f8').tobytes()).hexdigest()
            if coordinate_key not in by_coord:
                continue
            candidate = by_coord[coordinate_key]; require(np.array_equal(candidate['x'],x), 'Saved trial coordinate hash collision')
            cell_key = (candidate['candidate_id'],metal)
            if point['status'] != 'complete':
                blocked[cell_key] = {'candidate_id':cell_key[0],'metal':metal,'point':record,
                    'reason':'Earlier identical-state molecular attempt failed; no automatic retry'}
                continue
            _,components = pool.point_components(world,source,metal,point)
            if cell_key in cells:
                require(cells[cell_key]['components'] == components, 'Conflicting complete trial/candidate energies')
            else:
                cells[cell_key] = {'candidate_id':cell_key[0],'metal':metal,'status':'accepted','components':components,
                    'reused':True,'evidence':[{'exact_saved_trial':record,'parent':row['parent'],'trial_added_as_candidate':False}]}
    blocked = {k:v for k,v in blocked.items() if k not in cells}
    missing = [(c['candidate_id'],m) for c in candidates for m in METALS
               if (c['candidate_id'],m) not in cells and (c['candidate_id'],m) not in blocked]
    require(len(candidates) <= 18 and len(missing) <= 8, 'Declared expanded pool/cross-call ceiling exceeded')
    return admitted,missing,list(blocked.values())


def compact(menu):
    return {k:v for k,v in menu.items() if k in ('status','candidate_ids','selected','exact_ties','missing_cells','reason',
        'independently_selected_Dy_minus_La','La_selected_geometry_Dy_minus_La','Dy_selected_geometry_Dy_minus_La')}


def summarize(pool, prior, candidates, aliases, cells, slots):
    local = {s:{'prior': prior['summary']['source_local_menus'][s]['expanded']} for s in SOURCES}
    global_menus = {'prior': prior['summary']['expanded_global_pool']}; deltas = {}; spread = {}
    for arm in ARMS:
        arm_slots = [s for s in slots if s['arm'] == arm]
        for source in SOURCES:
            names = [source+':'+arm+'_terminal_'+m for m in METALS]
            local[source][arm] = pool.pool_summary(local[source]['prior']['candidate_ids']+[aliases[n] for n in names],cells) if all(n in aliases for n in names) else {
                'status':'unavailable','reason':'A declared source/arm search is unavailable'}
        names = [s+':'+arm+'_terminal_'+m for s in SOURCES for m in METALS]
        global_menus[arm] = pool.pool_summary(global_menus['prior']['candidate_ids']+[aliases[n] for n in names],cells) if all(s['status']=='admitted' for s in arm_slots) else {
            'status':'unavailable','reason':'A complete arm requires its four declared searches'}
    for name in ('prior',*ARMS):
        pair = [local[s][name] for s in SOURCES]
        if all(m['status']=='available' for m in pair):
            values = [m['independently_selected_Dy_minus_La']['composite_kcal_mol'] for m in pair]
            spread[name] = {'7CCN_minus_7CCO_kcal_mol':values[1]-values[0],'absolute_spread_kcal_mol':abs(values[1]-values[0])}
        else:
            spread[name] = {'status':'unavailable'}
    def change(a,b):
        if a['status']!='available' or b['status']!='available':
            return {'status':'unavailable'}
        return {'status':'available','contrast_change_kcal_mol':b['independently_selected_Dy_minus_La']['composite_kcal_mol']-a['independently_selected_Dy_minus_La']['composite_kcal_mol'],
            'selection_changed':{m:a['selected'][m]!=b['selected'][m] for m in METALS},
            'own_metal_work':{m:pool.difference(cells[(b['selected'][m],m)]['components'],cells[(a['selected'][m],m)]['components']) for m in METALS}}
    for before,after in (('prior','control'),('prior','candidate'),('control','candidate')):
        name = before+'_to_'+after
        deltas[name] = {'global':change(global_menus[before],global_menus[after]),
                       'source_local':{s:change(local[s][before],local[s][after]) for s in SOURCES}}
        if all('absolute_spread_kcal_mol' in spread[k] for k in (before,after)):
            deltas[name]['absolute_source_spread_change_kcal_mol'] = spread[after]['absolute_spread_kcal_mol']-spread[before]['absolute_spread_kcal_mol']
    complete = all(s['status']=='admitted' for s in slots) and all(global_menus[a]['status']=='available' for a in ARMS)
    return {'paired_arm_comparison_available':complete, 'source_local_menus':{s:{a:compact(m) for a,m in menus.items()} for s,menus in local.items()},
        'global_menus':{a:compact(m) for a,m in global_menus.items()},'source_spread':spread,'changes':deltas,
        'contrast_convention':'Dy minus La; changes are conditional composite energy response, not affinity',
        'coverage_rule':'Each arm requires four admitted searches and complete paired cells; paired arm comparison requires all eight.',
        'global_menu_caveat':'One shared global menu removes source choice by construction; source-local menus diagnose remaining source sensitivity.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['run'])
    for key in ('config','output'):
        parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args()
    cfg,pool,world,prior,candidates,aliases,cells,prior_slots,slots,range_module,base,parents=inputs(args.config.resolve())
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    result={'protocol_id':PROTOCOL,'stage':'run','status':'unavailable','implementation':pin(__file__),'config':pin(args.config),
        'plan':cfg['plan'],'range_pool_result':cfg['range_pool_result'],'basis_scout_result':cfg['basis_scout_result'],
        'basis_remaining_result':cfg['basis_remaining_result'],'pool_helper':prior['pool_helper'],
        'slurm_job_id':os.environ['SLURM_JOB_ID'],'prior_candidate_count':len(candidates),'prior_cell_count':len(cells),
        'declared_arm_states':8,'new_gradients_requested':0,'optimizer_steps':0,'actual_cost':empty_cost(),'cross_rows':[],'cross_batches':[],
        'interpretation_limits':'Consumed peptide test, conditional inventory/termini; no equilibrium, affinity, classifier or natural LanM validation.'}
    write(out/'START.json',result)
    try:
        admitted,missing,blocked=add_terminals(pool,world,candidates,aliases,cells,prior_slots,slots,range_module,base,parents[0]['fallback_rule'])
        result.update(searches=admitted,aliases=aliases,candidates=[{k:v for k,v in c.items() if k!='x'} for c in candidates],
            required_cross_cells=[{'candidate_id':i,'metal':m} for i,m in missing],unavailable_cached_cells=blocked)
        write(out/'POOL.json',{k:result[k] for k in ('searches','aliases','candidates','required_cross_cells')});result['pool']=pin(out/'POOL.json')
        for offset in range(0,len(missing),4):
            directory=out/('batch_%d'%(offset//4));directory.mkdir()
            batch={'actual_cost':empty_cost(),'cross_rows':[]}
            rows=pool.execute(world,candidates,missing[offset:offset+4],directory,batch)
            for row in rows:
                require(read(row['cell'])=={k:v for k,v in row.items() if k!='cell'},'New cross-cell receipt differs')
                if row['status']=='accepted':
                    cells[(row['candidate_id'],row['metal'])]=row
            result['cross_rows'].extend(rows)
            result['cross_batches'].append({k:v for k,v in batch.items() if k!='cross_rows'})
            for key,value in batch['actual_cost'].items():
                result['actual_cost'][key]+=value
        result['summary']=summarize(pool,prior,candidates,aliases,cells,admitted)
        result['paired_matrix']=[{'candidate_id':c['candidate_id'],'aliases':c['aliases'],
            'cells':{m:cells.get((c['candidate_id'],m),{'status':'unavailable'}) for m in METALS}} for c in candidates]
        result['status']='complete' if result['summary']['paired_arm_comparison_available'] else 'unavailable'
        result['arm_search_costs']=[p['actual_cost'] for p in parents]
        result['incremental_call_counts_including_searches']={k:result['actual_cost'][k]+sum(p['actual_cost'][k] for p in parents)
            for k in ('new_composite_evaluations','MACE_calls_attempted','GFN2_calls_attempted')}
        require(result['actual_cost']['MACE_calls_attempted']<=8 and result['actual_cost']['GFN2_calls_attempted']<=16
                and result['incremental_call_counts_including_searches']['MACE_calls_attempted']<=200
                and result['incremental_call_counts_including_searches']['GFN2_calls_attempted']<=400,'Declared call ceiling exceeded')
        write(out/'SUMMARY.json',dict(result['summary'],search_slots=[{k:r[k] for k in ('source_id','metal','arm','status')} for r in admitted],
            incremental_call_counts=result['incremental_call_counts_including_searches'],range_pool_result=cfg['range_pool_result'],
            basis_scout_result=cfg['basis_scout_result'],basis_remaining_result=cfg['basis_remaining_result']))
        result['summary_artifact']=pin(out/'SUMMARY.json');pool.verify(result['config'])
    except Exception as exc:
        result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-started
    result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ['SLURM_CPUS_PER_TASK'])
    result['application_allocated_GPU_seconds']=result['elapsed_seconds']
    result['peak_host_RSS_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['allocation']={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','SLURM_MEM_PER_CPU','CUDA_VISIBLE_DEVICES')}
    write(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')}),flush=True)
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':
    raise SystemExit(main())
