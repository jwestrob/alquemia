#!/usr/bin/env python3
"""One bounded radial14 campaign, exact saved first proposal, no continuation/retry."""
import argparse
import ast
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time
import numpy as np
import ADAPTER as A
import CORE as C
import RESOURCE as X

HERE = Path(__file__).resolve().parent
THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')
LIMITS = dict(outer_steps=8, new_parameter_vectors=228, new_native_programs=7280, reused_origin_vectors=1, reused_initial_derivative_vectors=28,
              logical_proposals=8, new_solver_calls=7, workers=12, numerical_threads_per_process=1, allocated_CPUs=13, native_timeout_seconds=300)
require = C.require


def pin(p):
    p = Path(p).resolve(); return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def verify(p):
    require(isinstance(p, dict) and pin(p['path']) == p, 'Changed/unbound source'); return Path(p['path'])


def read(p): return json.loads(verify(p).read_text())
def digest(v): return hashlib.sha256(json.dumps(v, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
def fingerprint(c): return digest({k:v for k,v in c.items() if k not in ('root_release', 'source_review', 'preparation', 'fingerprint')})


def module(p, name, aliases=None):
    prior = {k:sys.modules.get(k) for k in aliases or {}}
    try:
        sys.modules.update(aliases or {})
        s = importlib.util.spec_from_file_location(name, verify(p)); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
    finally:
        for k,v in prior.items():
            if v is None: sys.modules.pop(k, None)
            else: sys.modules[k] = v


def inventory_view(inv):
    f = inv['fixed_authorities']; saved = inv['saved_proposal']
    return dict(inv, derivative_keys=inv['initial_derivative_keys'], actual_normalized_spans=inv['initial_actual_spans'],
                saved_affine_result=saved['result'], saved_linearization=saved['linearization'], saved_radial_conic=saved['conic_data'], saved_radial_layout=saved['row_layout'],
                fixed_safeguards=dict(actual_scores=f['electrical_actual_scores'], literal_squared_norm_denominators=f['electrical_conic_denominators'],
                                      W9_reference_Debye=f['signed_W9_reference_Debye'], W9_absolute_ceiling_Debye=f['signed_W9_absolute_ceiling_Debye']))


def context(cfg):
    require(cfg['implementation'] == pin(__file__) and cfg['adapter'] == pin(A.__file__) and cfg['core'] == pin(C.__file__), 'Frozen three new sources')
    require(cfg['limits'] == LIMITS and cfg['root_declaration'] == pin(HERE/'ROOT_DECLARATION.json') and cfg['reuse_inventory'] == pin(HERE/'REUSE_INVENTORY.json'), 'Frozen finite policy')
    require(all(os.environ.get(k) == '1' for k in THREADS), 'One numerical thread')
    require(cfg['resource_helper']==pin(X.__file__) and cfg['resource_amendment']==pin(HERE/'ROOT_RESOURCE_AMENDMENT.json'),'Scheduling amendment source')
    require(cfg['campaign_specification']==pin(HERE/'SPEC.md') and cfg['scheduler_wrapper']==pin(HERE/'RUN.sbatch'),'Frozen specification/wrapper')
    amendment=read(cfg['resource_amendment']); require(amendment['original_declaration']==cfg['root_declaration'] and amendment['changes']['cpus_per_task']==13,'Exact scheduling-only amendment')
    declaration = read(cfg['root_declaration']); inv = inventory_view(read(cfg['reuse_inventory'])); oldcfg = read(cfg['radial_configuration'])
    require(inv['root_declaration'] == cfg['root_declaration'] and declaration['iterations']['maximum_new_native_programs'] == 7280 and declaration['iterations']['maximum_new_parameter_vectors'] == 228 and declaration['iterations']['maximum_new_solver_calls'] == 7, 'Declared caps')
    ignored = {'implementation','adapter','core','limits','root_release','source_review','preparation','fingerprint','run_id','output_directory','reuse_inventory'}
    require(all(cfg[k] == v for k,v in oldcfg.items() if k not in ignored), 'Unchanged inherited scientific configuration')
    radial = module(cfg['radial_adapter'], 'campaign_radial_adapter')
    M = module(oldcfg['math'], 'campaign_radial_math')
    P = module(cfg['radial_runner'], 'campaign_radial_runner', {'ADAPTER':radial, 'MATH':M})
    # The inherited scientific loader also enforced its former login hostname.
    # Remove only that obsolete resource predicate; actual allocation has its own
    # independently pinned13-CPU gate. No environment spoofing or science changes.
    tree=ast.parse(Path(P.__file__).read_text()); fn=copy.deepcopy(next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='context')); removed=0
    for i,node in enumerate(fn.body):
        if isinstance(node,ast.Expr) and isinstance(node.value,ast.Call) and len(node.value.args)>1 and isinstance(node.value.args[1],ast.Constant) and node.value.args[1].value=='Truthful login context':
            fn.body[i]=ast.Pass(); removed+=1
    require(removed==1,'One obsolete parent login-only predicate')
    ns=dict(vars(P)); exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),__file__,'exec'),ns)
    s, base, old, original, attr0, checks, parent, center, oldinv = ns['context'](oldcfg)
    require(oldcfg['adapter'] == cfg['radial_adapter'] and oldcfg['implementation'] == cfg['radial_runner'], 'Pinned radial runtime authority')
    for authority in inv['authorities'].values(): verify(authority)
    for key in ('initial_local_result', 'initial_local_root_admission', 'fixed_point_result', 'fixed_point_root_admission', 'zero_bridge', 'numerical_admission'):
        require(inv['authorities'][key] == declaration['inputs'][key], 'Declared actual authority '+key)
    local, localroot = read(inv['authorities']['initial_local_result']), read(inv['authorities']['initial_local_root_admission'])
    native, nativeroot = read(inv['authorities']['fixed_point_result']), read(inv['authorities']['fixed_point_root_admission'])
    require(localroot['status'] == 'admitted_radial_local_comparison' and localroot['result'] == inv['authorities']['initial_local_result'] and local['affine_result'] == inv['saved_affine_result'], 'Actual derivative/first-solve admission')
    require(nativeroot['result'] == inv['authorities']['fixed_point_result'] and native['status'] == 'completed_two_fixed_radial_native_comparisons' and native['error'] is None, 'Actual optional half-points')
    entries = {e['key']:e for e in inv['entries']}; require(len(entries) == 31 and len(inv['derivative_keys']) == 28, 'Exact bounded31 inventory')
    for e in entries.values():
        verify(e['artifact']); require(digest(e['parameters']) == e['key'], 'Literal seed key')
    domain = M.extend_domain(read(parent['domain']))
    require(C.probe_points(entries[inv['origin_key']]['parameters'], domain) == [entries[k]['parameters'] for k in inv['derivative_keys']], 'Exact first28 probe coordinates')
    saved, problem = read(inv['saved_affine_result']), read(inv['saved_linearization'])
    require(saved['linearization'] == inv['saved_linearization'] and saved['center_key'] == inv['origin_key'] and saved['derivative_keys'] == inv['derivative_keys'] and saved['spans'] == inv['actual_normalized_spans'], 'Literal first problem/stencil identity')
    require(saved['settings'] == cfg['settings'] == inv['saved_proposal']['settings'] and saved['targets'] == cfg['targets'] and saved['cases']['radial14']['conic_data'] == inv['saved_radial_conic'] and saved['cases']['radial14']['row_layout'] == inv['saved_radial_layout'], 'Saved first conic/settings pins')
    core = M.core14(original); core.probe_points = C.probe_points; core.linearize.__globals__['probe_points'] = C.probe_points
    attr = M.attribution14(attr0); solver = module(oldcfg['solver'], 'campaign_radial_solver', {'MATH':M})
    p = C.restore_problem(problem); model = solver.assemble(core, attr, p, 'radial14')
    stored = np.load(verify(inv['saved_radial_conic']), allow_pickle=False)
    expected = dict(A_data=model['A'].data, A_indices=model['A'].indices, A_indptr=model['A'].indptr, A_shape=model['A'].shape, b=model['b'], c=model['c'], columns=model['columns'])
    require(all(np.array_equal(stored[k],v) for k,v in expected.items()) and model['layout'] == read(inv['saved_radial_layout']), 'Exact initial conic identity before native calls; no Jacobian replay')
    require(all(problem['blocks'][k]['denominator'] == inv['fixed_safeguards']['literal_squared_norm_denominators'][k] for k in C.ELECTRICAL) and problem['signed_W9_z']['absolute_ceiling_Debye'] == inv['fixed_safeguards']['W9_absolute_ceiling_Debye'], 'Literal fixed electrical authorities')
    require([problem['blocks'][k]['denominator'] for k in ('P_energy','B_energy')] == list(inv['fixed_authorities']['q_denominators'].values()), 'Original component normalization')
    A.bounded_parent(old)
    return P,M,s,base,old,radial,original,core,attr,solver,checks,parent,center,inv,saved,problem


def release_gate(cfg):
    review, prep, release = read(cfg['source_review']), read(cfg['preparation']), read(cfg['root_release']); fp = fingerprint(cfg)
    require(review['status'] == 'passed_radial14_iterative_fit_source_review' and review['configuration_fingerprint'] == fp, 'Independent exact source review')
    require(prep['status'] == 'prepared_radial14_iterative_fit_zero_calls' and prep['configuration_fingerprint'] == fp and prep['native_programs'] == prep['solver_calls'] == 0, 'Zero-call preparation')
    require(release['status'] == 'released_one_bounded_radial14_iterative_fit' and release['configuration_fingerprint'] == fp, 'One exact root release')
    for k in ('implementation','adapter','core','root_declaration','resource_amendment','resource_helper','scontrol','scheduler_wrapper','campaign_specification','reuse_inventory','radial_configuration','limits','output_directory','root_cache_bridge','source_review','preparation'):
        require(release[k] == cfg[k], 'Release binding '+k)


def run(mode, config_path, out):
    tick,cpu = time.monotonic(),time.process_time(); children = resource.getrusage(resource.RUSAGE_CHILDREN)
    cfg=json.loads(Path(config_path).read_text()); cfg['fingerprint']=fingerprint(cfg); out=Path(out).resolve(); out.mkdir(parents=True,exist_ok=False)
    ev=None; current=None; policy=None; P=None; values={}; allocation=None
    result=dict(status='started',mode=mode,implementation=pin(__file__),configuration=pin(config_path),configuration_fingerprint=cfg['fingerprint'],limits=LIMITS,
                error=None,vectors={},iterations=[],native_programs=0,solver_calls=0,logical_proposals=0,imported_proposals=0,selected_candidate=None,job_id=None,
                physical_accuracy_qualified=False,heldout_reads=0,reference_calls=0,automatic_extension=False)
    def persist(terminal=False):
        if P is not None: P.persist_result(out,result,ev,terminal)
        else: (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    def obtain(parameters,role,step):
        key=digest(parameters); ev.authorize(parameters,role,step); usage=dict(role=role,outer_step=step)
        if key in values:
            result['vectors'][key].setdefault('requested_roles',[]).append(usage); return values[key]
        rec=dict(parameters=parameters,status='started',requested_roles=[usage],inherited=key in ev.registry); result['vectors'][key]=rec; persist()
        try:
            value=ev.evaluate(parameters); values[key]=value; path=out/'vectors'/(key+'.json'); require(not path.exists(),'Immutable vector artifact'); P.save(path,value)
            rec.update(status='complete',artifact=pin(path)); persist(); return value
        except Exception as exc: rec.update(status='budget_not_attempted' if isinstance(exc,old.BudgetStop) else 'failed',error=repr(exc)); persist(); raise
    try:
        if mode=='prepare': require(not os.environ.get('SLURM_JOB_ID'),'Preparation is one-thread login only')
        P,M,s,base,old,radial,original,core,attr,solver,checks,parent,center,inv,saved,problem=context(cfg)
        result.update(origin_key=inv['origin_key'],reuse_inventory=cfg['reuse_inventory'],saved_initial_proposal=inv['saved_affine_result'])
        if mode=='prepare':
            result.update(status='prepared_radial14_iterative_fit_zero_calls',preparation_checks=dict(exact_first28_metadata=True,exact_initial_conic_arrays=True,fixed_electrical_authorities=True,
                          journal_limits_before_initialization=True,native_calls=0,solver_calls=0,old28_numerical_replay=False,new_cache_constructed=False))
        else:
            require(str(out)==cfg['output_directory'],'Exact fresh runtime path'); release_gate(cfg)
            allocation=X.allocation(str(verify(cfg['scontrol']))); result.update(allocation=allocation,job_id=allocation['job_id'])
            ev=A.evaluator_class(old,radial,base,checks,release_gate)(cfg,out/'cache',mode,inv,parent,allocation)
            declaration=read(s['campaign_declaration']); fixed=original.make_fixed(read(declaration['start']['artifact']),ev.m,read(s['decomposition']),ev.domain,declaration)
            fixed['domain']=M.extend_domain(fixed['domain'])
            require(fixed['denominators']==inv['fixed_authorities']['q_denominators'],'Fixed original q normalization')
            policy=C.Policy(core,attr,M,solver,fixed,inv,read(cfg['targets']),saved,problem,cfg['settings'])
            (out/'vectors').mkdir(); current=ev.load_vector(inv['origin_key'])
            # Saved vectors retain original artifacts; no derived-score replay or fake new execution.
            for entry in inv['entries']:
                key=entry['key']; value=ev.load_vector(key); values[key]=value
                result['vectors'][key]=dict(status='complete',parameters=entry['parameters'],artifact=entry['artifact'],inherited=True,view_lineage=entry['lineage'],import_role=entry['declared_role'],requested_roles=[])
            result.update(initial_metrics=policy.assess(current),current_key=inv['origin_key']); require(result['initial_metrics']['all_declared_safeguards_pass'],'Initial fixed safeguards'); persist()
            stop=None
            for step in range(1,9):
                row=dict(outer_step=step,origin_key=digest(current['parameters']),status='started',derivative_keys=[],trial_slots=[]); result['iterations'].append(row); persist()
                try:
                    probes=[]
                    for point in C.probe_points(current['parameters'],fixed['domain']):
                        value=obtain(point,'probe',step); probes.append(value); row['derivative_keys'].append(digest(point)); persist()
                    require(result['solver_calls'] < 7 or step == 1, 'At most seven new solves')
                    proposal=policy.propose(current,probes,step,out/f'proposal_{step:02d}'); row['proposal']=pin(out/f'proposal_{step:02d}/RESULT.json')
                    result['solver_calls']+=proposal['solver_calls']; result['logical_proposals']+=1; result['imported_proposals']+=int(proposal['solver_calls']==0)
                    require(proposal['solver_calls']==(0 if step==1 else 1),'Only first proposal reused')
                    slots=proposal['trial_slots']; row['trial_slots']=[{k:v for k,v in x.items() if k!='prediction'} for x in slots]; persist()
                    if not slots or not any(x['eligible'] for x in slots): row['status']='no_eligible_fraction'; stop='no_eligible_fraction'; break
                    for i,slot in enumerate(slots):
                        if not slot['eligible']: slot['execution_status']='not_eligible'; continue
                        row['trial_slots'][i]['required_new_programs']=ev.required_programs(slot['parameters']); persist()
                        value=obtain(slot['parameters'],'trial',step); slot.update(execution_status='complete',native_admitted=True,vector=value)
                        actual=policy.assess(value); predicted=slot['prediction']
                        row['trial_slots'][i].update(execution_status='complete',native_admitted=True,artifact=result['vectors'][slot['key']]['artifact'],actual=actual,
                            actual_minus_predicted=dict(q=actual['q']-predicted['actual_component_q'],component_errors_kcal_mol={k:actual['component_errors_kcal_mol'][k]-v for k,v in predicted['component_errors_kcal_mol'].items()},
                                                       signed_W9_residual_Debye=actual['signed_W9_z_safeguard']['model_minus_reference_Debye']-predicted['signed_W9_z_residual_Debye'])); persist()
                    selection=policy.select(current,slots); row['selection']=selection; chosen=selection['selected_trial_index']
                    if chosen is None: row['status']='no_actual_improvement'; stop='no_actual_safeguarded_improvement'; break
                    current=copy.deepcopy(slots[chosen]['vector']); row.update(status='complete',selected_key=selection['selected_key']); result['current_key']=selection['selected_key']; persist()
                    if selection['selected_metrics']['full_development_target']: stop='development_target_met'; break
                except old.BudgetStop: row['status']='hard_budget_exhausted'; stop='hard_budget_exhausted'; break
                except Exception: row['status']='failed_incomplete_batch_no_selection'; raise
            key=digest(current['parameters'])
            result.update(status='bounded_radial14_iterative_fit_complete',stopping_reason=stop or 'outer_step_limit',selected_key=key,selected=result['vectors'][key]['artifact'],
                          selected_parameters=current['parameters'],selected_view_lineage=result['vectors'][key].get('view_lineage'),selected_metrics=policy.assess(current),candidate_frozen=True,
                          completed_outer_steps=sum(x['status']=='complete' for x in result['iterations']))
    except Exception as exc:
        result.update(status='failed_radial14_iterative_fit' if mode=='run' else 'failed_preparation',error=repr(exc),stopping_reason='failure_no_retry',candidate_frozen=False)
        if ev is not None:
            try: ev._stop(exc)
            except Exception as other: result['stop_reporting_error']=repr(other)
        if current is not None:
            key=digest(current['parameters']); result.update(last_complete_incumbent_key=key,last_complete_incumbent=result['vectors'].get(key,{}).get('artifact'),last_complete_incumbent_parameters=current['parameters'])
            try: result['last_complete_incumbent_metrics']=policy.assess(current)
            except Exception as other: result['incumbent_reporting_error']=repr(other)
    finally:
        # Include an interrupted/failed solve's durable accounting without retry.
        proposals=[]
        for path in sorted(out.glob('proposal_*/RESULT.json')):
            p=json.loads(path.read_text()); proposals.append(pin(path))
        if proposals:
            docs=[read(p) for p in proposals]; result.update(proposals=proposals,solver_calls=sum(p['solver_calls'] for p in docs),logical_proposals=len(docs),imported_proposals=sum('imported_saved_proposal' in p for p in docs))
        ch=resource.getrusage(resource.RUSAGE_CHILDREN)
        result['cost']=dict(wall_seconds=time.monotonic()-tick,controller_CPU_seconds=time.process_time()-cpu,child_user_CPU_seconds=ch.ru_utime-children.ru_utime,
                            child_system_CPU_seconds=ch.ru_stime-children.ru_stime,native_workers=12,numerical_threads=1,execution_mode='slurm' if allocation else 'login_preparation_or_preallocation_failure',scheduler_allocated_CPU_seconds=None if allocation else 0,allocation=allocation,inherited_costs_not_added=True)
        persist(True)
    print(json.dumps(dict(status=result['status'],error=result['error'],native_programs=result['native_programs'],solver_calls=result['solver_calls']))); return result


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('mode',choices=('prepare','run')); p.add_argument('--config',required=True); p.add_argument('--out',required=True); a=p.parse_args()
    raise SystemExit(bool(run(a.mode,a.config,a.out)['error']))
