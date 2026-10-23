#!/usr/bin/env python3
"""One conditional saved SOCP; no native evaluator, candidate or fit API."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

import clarabel
import numpy as np
from scipy import sparse

REMOVED = {'block/old_EF': (0, 414), 'block/pair_mu': (414, 517),
           'block/C_mu': (517, 530), 'block/chi': (530, 535)}
RUNTIME = ('terminal_campaign', 'full_comparison', 'full_actual_review',
           'full_actual_reviewer', 'source_review', 'root_release')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def read(p):
    require(isinstance(p, dict) and pin(p['path']) == p, 'Changed/unbound pin')
    return json.loads(Path(p['path']).read_text())


def load_module(p, name):
    require(pin(p['path']) == p, 'Changed implementation')
    spec = importlib.util.spec_from_file_location(name, p['path'])
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def fingerprint(cfg):
    x = {k: v for k, v in cfg.items() if k not in RUNTIME}
    return hashlib.sha256(json.dumps(x, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def context(cfg):
    require(cfg['implementation'] == pin(__file__), 'Exact diagnostic source')
    base = load_module(cfg['comparison_implementation'], 'frozen_full_conic_comparison')
    require(all(os.environ.get(k) == '1' for k in base.THREADS), 'One numerical thread')
    read(cfg['declaration']) if cfg['declaration']['path'].endswith('.json') else require(
        pin(cfg['declaration']['path']) == cfg['declaration'], 'Exact declaration')
    # Every following authority is mandatory. Null pins prevent early execution.
    review = read(cfg['full_actual_review'])
    require(pin(cfg['full_actual_reviewer']['path']) == cfg['full_actual_reviewer'],
            'Exact full-comparison actual reviewer')
    require(review['status'] == 'passed_saved_conic_actual_review'
            and review['implementation'] == cfg['full_actual_reviewer']
            and review['result'] == cfg['full_comparison']
            and review['source_review'] == cfg['full_source_review']
            and review['certificate_checks_passed'] is True
            and review['solver_success'] is True and review['solver_status'] == 'Solved',
            'Successful independently reviewed full-conic baseline required')
    full = read(cfg['full_comparison'])
    require(cfg['full_comparison']['path'] == cfg['expected_full_result_path'],
            'Only the frozen full-comparison result')
    require(full['status'] == 'completed_saved_conic_comparison'
            and full['error'] is None and full['solver_calls'] == 1
            and full['solver_success'] is True and full['solver_status'] == 'Solved'
            and full['implementation'] == cfg['comparison_implementation'],
            'A failed/reduced baseline is not an interpretable counterfactual')
    full_cfg = read(full['configuration'])
    require(all(full_cfg[k] == cfg[k] for k in ('inputs', 'environment', 'settings')),
            'Exact inherited scientific inputs/environment/settings')
    require(full_cfg['source_review'] == cfg['full_source_review']
            and full_cfg['terminal_campaign'] == cfg['terminal_campaign']
            and full['terminal_campaign'] == cfg['terminal_campaign'],
            'Exact baseline source review and terminal campaign')
    # Reuse the admitted environment, terminal, source and root-release gates.
    inputs, fit, settings = base.context(full_cfg)
    require(full['origin_key'] == inputs['origin_key']
            and full['derivative_keys'] == inputs['derivative_keys'],
            'Same origin24 and exact24 probes')
    require(full['effective_settings'] == base.attributes(settings),
            'Actual baseline settings unchanged')
    own_review = read(cfg['source_review']); release = read(cfg['root_release'])
    require(own_review['status'] == 'passed_joint_ceiling_removal_source_review'
            and own_review['implementation'] == cfg['implementation']
            and own_review['configuration_fingerprint'] == fingerprint(cfg),
            'Exact own source review')
    require(release['status'] == 'released_one_joint_ceiling_removal_diagnostic'
            and release['configuration_fingerprint'] == fingerprint(cfg)
            and all(release[k] == cfg[k] for k in RUNTIME if k != 'root_release')
            and release['maximum_solver_calls'] == 1 and release['new_native_calls'] == 0,
            'One exact conditional diagnostic release')
    return base, full, inputs, settings


def problem_arrays(saved):
    p = copy.deepcopy(saved)
    for key in ('center', 'lower', 'upper', 'width', 'z', 'bounds_lower',
                'bounds_upper', 'spans', 'Jacobian_original531', 'Jacobian_components'):
        p[key] = np.asarray(p[key], dtype=float)
    for b in p['blocks'].values():
        b['r'] = np.asarray(b['r'], dtype=float); b['J'] = np.asarray(b['J'], dtype=float)
    for b in p['forces'].values():
        for key in ('r', 'J', 'weights'): b[key] = np.asarray(b[key], dtype=float)
    for key in ('r', 'J'): p['energy'][key] = np.asarray(p['energy'][key], dtype=float)
    return p


def remove_joint_ceilings(full_model):
    require(full_model['A'].shape == (1817, 13), 'Original1817x13 form')
    keep = np.ones(1817, dtype=bool); removed = []; layout = []; cones = []; cursor = 0
    for row in full_model['layout']:
        count = row['stop']-row['start']
        if row['name'] in REMOVED:
            require(row['kind'] == 'SOC' and (row['start'], row['stop']) == REMOVED[row['name']],
                    'Exact removed cone/rows')
            keep[row['start']:row['stop']] = False; removed.append(copy.deepcopy(row)); continue
        entry = dict(row, start=cursor, stop=cursor+count); layout.append(entry); cursor += count
        if row['kind'] == 'SOC': cones.append(clarabel.SecondOrderConeT(count))
    require({r['name'] for r in removed} == set(REMOVED) and (~keep).sum() == 535,
            'Only all four specified ceilings removed jointly')
    require(len(cones) == 17 and sum(r['stop']-r['start'] for r in layout if r['kind']=='SOC') == 426,
            'P/B2 plus forceRMS15 cones')
    require(sum(r['kind']=='nonnegative' for r in layout) == 856 and cursor == 1282,
            'All original linear sides and1282 retained rows')
    cones.append(clarabel.NonnegativeConeT(856))
    return {'P': full_model['P'].copy(), 'c': full_model['c'].copy(),
            'A': full_model['A'][keep, :].tocsc(), 'b': full_model['b'][keep].copy(),
            'cones': cones, 'layout': layout}, removed, np.flatnonzero(keep).tolist()


def correspondence(p, model, x, core):
    values, _, names = core.constraints(p, np.r_[x[:12], x[-1]**2])
    expected = {name: value for name, value in zip(names, values) if name not in REMOVED}
    slack = model['b']-model['A']@x; observed = {}
    for row in model['layout']:
        s = slack[row['start']:row['stop']]
        if row['kind']=='SOC': observed[row['name']] = float(s[0]**2-s[1:]@s[1:])
        elif not row['name'].startswith('bound/'): observed[row['name']] = float(s[0])
    require(set(expected)==set(observed) and len(expected)==847, 'All847 retained source inequalities')
    a=np.asarray([observed[k] for k in expected]);b=np.asarray(list(expected.values()))
    core.near(a,b,'Reduced cone/source correspondence')
    return {'retained_source_inequalities':847,'maximum_arithmetic_difference':float(np.max(abs(a-b))),
            'minimum_retained_source_margin':float(np.min(b)),'nonnegative_t':bool(x[-1]>=0)}


def reports(p, x, status, core):
    require(np.shape(x)==(13,) and np.all(np.isfinite(x)) and x[-1]>=0,
            'Finite nonnegative epigraph; no clipping')
    answer=[]
    for f in core.FRACTIONS:
        slot=core.certify_fraction(p,np.r_[x[:12],x[-1]**2],f,status=='Solved')
        cert=slot['certificate'];margins=cert['residual_constraints']
        retained={k:v for k,v in margins.items() if k not in REMOVED}
        block_rows={}
        for name,b in p['blocks'].items():
            origin_norm=math.fsum(float(v)*float(v) for v in b['r'])
            predicted=cert['blocks'][name]['predicted_squared_norm']
            # CORE's original531 vector gives quarter-weighted old block norms.
            factor=1. if name in ('P_energy','B_energy') else 4.
            block_rows[name]={'origin_score':factor*origin_norm,
                'predicted_score':factor*predicted,'predicted_change':factor*(predicted-origin_norm),
                'ratio_to_frozen_normalizer':cert['blocks'][name]['ratio'],
                'removed_ceiling':('block/'+name in REMOVED),
                'original_ceiling_pass':cert['blocks'][name]['ratio']<=1 if 'block/'+name in REMOVED else None}
        retained_ok=all(v>=0 for v in retained.values())
        if status!='Solved': retained_ok=all(v>0 for v in retained.values())
        answer.append({'fraction':f,'source_certificate_under_original_policy':slot,
            'all_six_blocks':block_rows,'retained_source_margins':retained,
            'removed_source_margins':{k:margins[k] for k in REMOVED},
            'retained_predicted_feasible':cert['in_box_and_trust'] and cert['all25_predicted'] and retained_ok,
            'all25_original_gates':cert['original25'],'native_evaluated':False,
            'selection_authorized':False})
    return answer


def run(config_path):
    wall=time.monotonic();cpu=time.process_time();cfg=json.loads(Path(config_path).read_text())
    base,full,inputs,settings=context(cfg)
    out=Path(cfg['output_directory']);out.mkdir(parents=True,exist_ok=False)
    result={'status':'started','implementation':pin(__file__),'configuration':pin(config_path),
        'configuration_fingerprint':fingerprint(cfg),'terminal_campaign':cfg['terminal_campaign'],
        'full_comparison':cfg['full_comparison'],'full_actual_review':cfg['full_actual_review'],
        'solver_calls':0,'new_native_calls':0,'new_candidate_evaluations':0,
        'selected_candidate':None,'current_fit_policy_changed':False,'error':None}
    base.save(out/'RESULT.json',result)
    try:
        core=base.module(inputs['core'])
        saved=read(full['linearized_problem']);p=problem_arrays(saved)
        full_model=base.assemble(p)
        require(pin(full['conic_data']['path'])==full['conic_data'],'Exact full sparse authority')
        with np.load(full['conic_data']['path'],allow_pickle=False) as raw:
            for key,array in [('A_data',full_model['A'].data),('A_indices',full_model['A'].indices),
                              ('A_indptr',full_model['A'].indptr),('b',full_model['b']),('c',full_model['c'])]:
                require(np.array_equal(raw[key],array),'Same full-comparison array '+key)
            require(tuple(raw['A_shape'])==full_model['A'].shape,'Same full shape')
        require(read(full['row_layout'])==full_model['layout'],'Exact full row ordering')
        model,removed,retained_indices=remove_joint_ceilings(full_model)
        result.update(origin_key=full['origin_key'],derivative_keys=full['derivative_keys'],
            linearized_problem=full['linearized_problem'],removed_cones=removed,
            retained_original_row_indices=retained_indices,effective_settings=base.attributes(settings),
            fixed_denominators=full['fixed_denominators'],removed_fixed_ceilings=full['fixed_ceilings'])
        result['zero_step_correspondence']=correspondence(p,model,np.r_[np.zeros(12),math.sqrt(p['q_current'])],core)
        fx=np.asarray(full['solution']['x'],dtype=float)
        result['full_constraint_baseline']={'solver_status':full['solver_status'],'returned_q':full['returned_q'],
            'fractions':reports(p,fx,full['solver_status'],core),'info':full['info'],
            'independent_linear_residuals':full['independent_linear_residuals']}
        np.savez(out/'CONIC_DATA.npz',A_data=model['A'].data,A_indices=model['A'].indices,
            A_indptr=model['A'].indptr,A_shape=model['A'].shape,b=model['b'],c=model['c'])
        base.save(out/'ROW_LAYOUT.json',model['layout'])
        result['conic_data']=pin(out/'CONIC_DATA.npz');result['row_layout']=pin(out/'ROW_LAYOUT.json')
        with open(cfg['one_attempt_marker'],'x') as stream:
            json.dump({'configuration':pin(config_path),'full_comparison':cfg['full_comparison'],
                       'pid':os.getpid()},stream)
        result['attempt_reserved']=True;base.save(out/'RESULT.json',result)
        solver=clarabel.DefaultSolver(model['P'],model['c'],model['A'],model['b'],model['cones'],settings)
        require(base.attributes(solver.get_settings())==base.attributes(settings),'Exact inherited backend settings')
        solver.print_to_buffer();result['solver_calls']=1;base.save(out/'RESULT.json',result)
        solution=solver.solve()  # Sole additional call; no other solve or retry.
        info=solver.get_info();(out/'SOLVER.log').write_text(solver.get_print_buffer())
        result['solution']=base.plain(base.attributes(solution))
        info_raw=base.attributes(info);info_raw['linsolver']=base.attributes(info.linsolver)
        result['info']=base.plain(info_raw);result['solver_status']=str(solution.status)
        result['solver_success']=result['solver_status']=='Solved'
        result['reduced_accuracy_status']=result['solver_status'].startswith('Almost')
        base.save(out/'RESULT.json',result)
        result['independent_linear_residuals']=base.plain(base.diagnostics(model,solution))
        x=np.asarray(solution.x,dtype=float)
        if x.shape==(13,) and np.all(np.isfinite(x)) and x[-1]>=0:
            result['returned_q']=float(x[-1]**2)
            result['q_difference_relaxed_minus_full']=result['returned_q']-full['returned_q']
            result['retained_correspondence']=correspondence(p,model,x,core)
            result['fractions']=reports(p,x,result['solver_status'],core)
        else:
            result['fractions']=[];result['return_rejection']='Malformed/nonfinite/negative t; no repair'
        result['status']='completed_joint_ceiling_removal_diagnostic'
        result['interpretation']='Conditional affine tradeoff only. All removed block changes and retained25 gates reported. No candidate, native validation, policy change, or nonlinear/global capacity claim.'
    except Exception as exc:
        result.update(status='failed_joint_ceiling_removal_diagnostic',error=repr(exc))
    finally:
        result['cost']={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,
                        'threads':1,'scheduler_job_id':None}
        base.save(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'solver_calls':result['solver_calls'],
                      'result':pin(out/'RESULT.json'),'error':result['error']}))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True,type=Path)
    args=parser.parse_args();raise SystemExit(bool(run(args.config)['error']))
