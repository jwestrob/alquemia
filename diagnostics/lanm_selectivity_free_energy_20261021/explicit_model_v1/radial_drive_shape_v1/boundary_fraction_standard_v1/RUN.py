#!/usr/bin/env python3
"""Prepare/carry out exactly one saved1/16 witness; no solve or fit interface."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
import types
import numpy as np
import ADAPTER as A
HERE = Path(__file__).resolve().parent
THREADS = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
LIMITS = dict(new_parameter_vectors=1,new_native_programs=56,new_cells=46,new_solver_calls=0,new_derivative_probes=0,new_reference_calls=0,workers=11,numerical_threads_per_process=1,allocated_CPUs=12,native_timeout_seconds=300)
FRACTION = .0625

def require(ok,message):
    if not ok: raise ValueError(message)
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verify(p):
    require(isinstance(p,dict) and pin(p['path'])==p,'Exact bound file');return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def fingerprint(c):return digest({k:v for k,v in c.items() if k not in ('root_release','source_review','preparation','fingerprint')})
def module(p,name,aliases=None):
    prior={k:sys.modules.get(k) for k in aliases or {}}
    try:
        sys.modules.update(aliases or {});s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
    finally:
        for k,v in prior.items():
            if v is None:sys.modules.pop(k,None)
            else:sys.modules[k]=v


def scheduling_bridge(cfg):
    """Read metadata only; retain the original scalar certificate and fingerprint."""
    require(cfg['scheduling_declaration']==pin(HERE/'ROOT_DECLARATION.json'), 'Exact scheduling declaration')
    declaration=read(cfg['scheduling_declaration'])
    require(declaration['status']=='declared_source_only_scheduling_recovery', 'Declared recovery')
    submission=read(declaration['original_submission'])
    require(submission['job_id']=='20734' and submission['configuration']==cfg['original_configuration'], 'Original unmodified submission')
    old=read(cfg['original_configuration']);bridge=read(cfg['scheduling_bridge'])
    require(bridge['status']=='metadata_only_boundary_scheduling_bridge' and bridge['original_configuration']==cfg['original_configuration']
            and bridge['scheduling_declaration']==cfg['scheduling_declaration'], 'Exact bridge authority')
    changed={'implementation','adapter','resource_helper','scheduler_wrapper','limits','run_id','output_directory','source_review','root_release','fingerprint'}
    added={'scheduling_declaration','scheduling_bridge','original_configuration'}
    require(set(cfg)==set(old)|added and all(cfg[k]==v for k,v in old.items() if k not in changed), 'Only declared scheduling configuration differences')
    oldlimits=dict(LIMITS,workers=12,allocated_CPUs=13)
    require(old['limits']==oldlimits and cfg['limits']==LIMITS and bridge['old_limits']==oldlimits and bridge['new_limits']==LIMITS, 'Only concurrency limits changed')
    require(cfg['resource_helper']==pin(HERE/'RESOURCE.py'), 'Versioned actual allocation helper')
    for k in ('implementation','adapter','resource_helper','scheduler_wrapper'):
        require(bridge[k]==cfg[k], 'Bridge source '+k);verify(cfg[k])
    prep=read(cfg['preparation']);review=read(old['source_review']);witness=read(prep['witness'])
    require(fingerprint(old)==bridge['original_configuration_fingerprint']==prep['configuration_fingerprint']==review['configuration_fingerprint'], 'Preserved original preparation fingerprint')
    require(prep['status']=='prepared_boundary_fraction_witness_zero_calls' and prep['implementation']==old['implementation']
            and prep['native_programs']==prep['solver_calls']==0 and prep['eligible'], 'Original zero-call eligible preparation')
    require(review['status']=='passed_boundary_fraction_witness_source_review' and review['implementation']==old['implementation']
            and review['adapter']==old['adapter'] and review['preparation']==cfg['preparation'] and review['witness']==prep['witness'], 'Original reviewed scalar evidence')
    require(bridge['preparation']==cfg['preparation'] and bridge['witness']==prep['witness']==read(old['root_release'])['witness']
            and witness['eligible'] and witness['fraction']==FRACTION and witness['key']==prep['witness_key']==bridge['witness_key'], 'Unchanged exact prepared point')
    return bridge


def metadata(cfg,probe_pins=False):
    require(cfg['implementation']==pin(__file__) and cfg['adapter']==pin(A.__file__),'Exact new sources')
    scheduling_bridge(cfg)
    original=HERE.parent/'boundary_fraction_witness_v1'
    require(cfg['root_declaration']==pin(original/'ROOT_DECLARATION.json') and cfg['witness_specification']==pin(original/'SPEC.txt') and cfg['scheduler_wrapper']==pin(HERE/'RUN.sbatch'),'Frozen original science and new wrapper')
    require(cfg['limits']==LIMITS and cfg['threads']==list(THREADS) and all(os.environ.get(k)=='1' for k in THREADS),'Exact finite limits/one numerical thread')
    dec=read(cfg['root_declaration']);require(dec['status']=='declared_source_and_scalar_preparation_only' and dec['policy']['fraction']==FRACTION,'One exact prospective fraction')
    require(dec['policy']['limits']==dict(new_solver_calls=0,new_derivative_probes=0,new_vectors=1,new_cells=46,new_native_programs=56,new_reference_calls=0,native_program_timeout_seconds=300),'Root finite scope')
    require(cfg['witness_inputs']==dec['inputs']['INPUTS.json'],'Root exact saved evidence')
    inp=read(cfg['witness_inputs']);fitcfg=read(cfg['fit_configuration'])
    require(cfg['fit_configuration']==inp['terminal_configuration'],'Exact old runtime configuration')
    changed={'implementation','adapter','limits','run_id','output_directory','root_release','source_review','preparation','fingerprint','root_declaration','scheduler_wrapper','resource_helper'}
    require(all(cfg[k]==v for k,v in fitcfg.items() if k not in changed),'Inherited scientific/native/safety configuration unchanged')
    root=read(inp['root_actual_admission']);producer=read(inp['terminal_producer']);origin=read(inp['origin']);proposal=read(inp['step8_proposal']);linear=read(inp['linearization'])
    require(root['status']=='admitted_bounded_radial14_development_fit' and root['producer']==inp['terminal_producer'] and root['comparison']==inp['terminal_comparison'] and root['selected_artifact']==inp['origin'],'Admitted terminal authority')
    require(producer['status']=='bounded_radial14_iterative_fit_complete' and producer['error'] is None and producer['configuration']==cfg['fit_configuration'] and producer['selected']==inp['origin'] and producer['proposals'][-1]==inp['step8_proposal'],'Exact terminal/last proposal')
    require(producer['stopping_reason']=='no_actual_safeguarded_improvement' and producer['completed_outer_steps']==7,'Preserved original finite stop')
    require(proposal['status']=='completed_proposal' and proposal['error'] is None and proposal['outer_step']==8 and proposal['solver_status']=='AlmostSolved' and proposal['solver_success'] is False,'Exact eighth AlmostSolved return')
    require(proposal['linearization']==inp['linearization'] and proposal['conic_data']==inp['conic_data'] and proposal['row_layout']==inp['row_layout'],'Exact saved affine authority')
    require(proposal['origin_key']==inp['origin_key']==producer['selected_key']==digest(origin['parameters']) and origin['status']=='complete' and linear['center']==origin['parameters'],'Same exact terminal/linearization origin')
    rows=inp['ordered_derivative_artifacts'];require(len(rows)==28 and [r['key'] for r in rows]==proposal['derivative_keys']==producer['iterations'][-1]['derivative_keys'],'Ordered last28 provenance')
    if probe_pins:
        for row in rows:verify({k:row[k] for k in ('path','sha256')})
    require(linear['spans']==producer['iterations'][-1].get('spans',linear['spans']) and len(linear['spans'])==14,'Literal14 spans retained')
    require(root['all25_physical_pass'] and root['electrical_safeguards_pass'],'Admitted incumbent safeguards')
    return inp,fitcfg,producer,origin,proposal,linear


def context(cfg,fitcfg,origin,linear):
    # Use pure loaders/assessment and qualified native methods only. No old context,
    # first-stencil replay, evaluator construction, solver assembly or constructor.
    C=module(fitcfg['core'],'boundary_campaign_core');X=module(cfg['resource_helper'],'boundary_slurm_resource')
    CA=module(fitcfg['adapter'],'boundary_campaign_adapter')
    D=module(fitcfg['implementation'],'boundary_campaign_driver',{'ADAPTER':CA,'CORE':C,'RESOURCE':X})
    radial_cfg=read(fitcfg['radial_configuration']);sources=read(radial_cfg['inputs'])['sources']
    for p in sources.values():verify(p)
    radial=module(fitcfg['radial_adapter'],'boundary_radial_adapter');M=module(fitcfg['math'],'boundary_radial_math')
    P=module(fitcfg['radial_runner'],'boundary_radial_persistence',{'ADAPTER':radial,'MATH':M})
    parent=read(cfg['parent_configuration']);base=module(cfg['parent_evaluator'],'boundary_base_evaluator');old=module(cfg['original_adapter'],'boundary_journal_adapter')
    require(parent['implementation']==cfg['parent_evaluator'] and pin(old.J.__file__)==sources['journal'],'Qualified evaluator/journal')
    original=module(sources['core'],'boundary_original_core');attr0=module(sources['attribution_math'],'boundary_attribution')
    sys.path.insert(0,str(Path(sources['qualification_checks']['path']).parent));checks=module(sources['qualification_checks'],'boundary_native_readbacks')
    solver=module(cfg['solver'],'boundary_scalar_certificate',{'MATH':M})
    env=read(cfg['environment']);require(str(Path(sys.prefix).resolve())==env['prefix'],'Pinned numerical environment')
    import importlib.metadata
    for k,v in env['versions'].items():require(importlib.metadata.version(k)==v,'Pinned package '+k)
    for p in env['files']:verify(p)
    build=read(cfg['build']);require(build['status']=='complete' and build['manifest']==cfg['native_source_manifest'] and build['executables']==cfg['executables'],'Qualified radial build')
    for p in cfg['executables'].values():verify(p)
    P.bridge_gate(radial_cfg)  # Reuses source/admission records, never native qualification.
    m=read(parent['manifest']);domain=read(parent['domain']);decl=read(sources['campaign_declaration'])
    fixed=original.make_fixed(read(decl['start']['artifact']),m,read(sources['decomposition']),domain,decl);fixed['domain']=M.extend_domain(fixed['domain'])
    core=M.core14(original);attr=M.attribution14(attr0);inv=D.inventory_view(read(cfg['reuse_inventory']))
    scalar=types.SimpleNamespace(certificate=lambda c,a,p,x,f,ok:fixed_fraction_certificate(c,a,p,x,f,ok,solver.certificate))
    policy=C.Policy(core,attr,M,scalar,fixed,inv,read(cfg['targets']),{'settings':cfg['settings']},None,cfg['settings'])
    p=C.restore_problem(linear)
    require(fixed['denominators']==inv['fixed_authorities']['q_denominators'],'Original q denominators')
    require(all(p['blocks'][k]['denominator']==inv['fixed_safeguards']['literal_squared_norm_denominators'][k] for k in C.ELECTRICAL),'Fixed feca conic ceilings')
    require(p['signed_W9_z']['absolute_ceiling_Debye']==inv['fixed_safeguards']['W9_absolute_ceiling_Debye'] and p['signed_W9_z']['reference_Dy_minus_La_z_Debye']==inv['fixed_safeguards']['W9_reference_Debye'],'Fixed signed W9 authority')
    require(len(m['rows'])==46 and len(m['static'])==48 and set(origin['cells'])==set(m['rows']) and all(v['status']=='complete' for v in origin['cells'].values()),'Complete admitted incumbent46+48')
    return C,X,P,base,old,radial,checks,parent,policy,p


def fixed_fraction_certificate(core,attr,p,point,fraction,success,qualified):
    require(fraction==FRACTION, 'Only declared1/16 external fraction')
    original=core.certify_fraction
    require(original.__globals__['FRACTIONS']==(1.,.5,.25,.125), 'Exact inherited fraction policy')
    ns=dict(original.__globals__,FRACTIONS=(FRACTION,))
    raw=types.FunctionType(original.__code__,ns,original.__name__,original.__defaults__,original.__closure__)
    facade=types.SimpleNamespace(**vars(core));calls=[]
    def routed(problem,returned,f,solver_success):
        require(problem is p and solver_success is success, 'Same exact certificate problem/status')
        if not calls:
            require(f==FRACTION, 'Raw fraction is exactly1/16');calls.append(f)
            return raw(problem,returned,f,solver_success)
        require(calls==[FRACTION] and f==1., 'Only internal full check of fractioned auxiliary witness')
        calls.append(f);return original(problem,returned,f,solver_success)
    facade.certify_fraction=routed
    answer=qualified(facade,attr,p,point,fraction,success)
    require(calls==[FRACTION,1.], 'Exactly raw and auxiliary scalar checks')
    return answer


def certify(C,policy,p,proposal):
    x=np.asarray(proposal['solution']['x'],dtype=float)
    require(x.shape==(15,) and np.isfinite(x).all() and x[-1]>=0,'Finite saved normalized14+t return')
    point=np.r_[x[:14],x[-1]**2]
    slot=policy.fraction(p,point,FRACTION,False)
    require(slot['fraction']==FRACTION and slot['solver_success'] is False and digest(slot['parameters'])==slot['key'],'Exact fraction and provenance')
    return slot


def release_gate(cfg):
    fp=fingerprint(cfg);review=read(cfg['source_review']);prep=read(cfg['preparation']);release=read(cfg['root_release'])
    require(review['status']=='passed_boundary_fraction_witness_source_review' and review['configuration_fingerprint']==fp and review['implementation']==cfg['implementation'] and review['adapter']==cfg['adapter'],'Exact independent source review')
    scheduling_bridge(cfg)
    require(review['scheduling_bridge']==cfg['scheduling_bridge'] and review['scheduling_declaration']==cfg['scheduling_declaration'], 'Reviewed scheduling bridge')
    require(prep['status']=='prepared_boundary_fraction_witness_zero_calls' and prep['native_programs']==prep['solver_calls']==0,'Original eligible exact scalar preparation')
    require(release['original_job_20734_cancelled_unstarted'] is True, 'Root admits old job cancellation before replacement');verify(release['original_job_resolution'])
    witness=read(prep['witness']);require(witness['eligible'] and witness['fraction']==FRACTION and witness['key']==prep['witness_key'],'Exact admitted one-point witness')
    require(release['status']=='released_one_fixed_boundary_fraction_witness' and release['configuration_fingerprint']==fp and release['witness']==prep['witness'],'One exact root release')
    for k in ('implementation','adapter','root_declaration','witness_inputs','fit_configuration','core','resource_helper','scheduler_wrapper','witness_specification','limits','output_directory','source_review','preparation','scheduling_declaration','scheduling_bridge','original_configuration'):
        require(release[k]==cfg[k],'Release binding '+k)
    return witness


def margins(metrics):
    physical={}
    for k,v in metrics['original25'].items():
        physical[k]=({'absolute_energy_margin_kcal_mol':1-abs(v['error_kcal_mol'])} if k.startswith('energy/') else
                     {'RMS_margin_kcal_mol_A':2-v['rigid_RMS_kcal_mol_A'],'maximum_margin_kcal_mol_A':5-v['maximum_body_component_error_kcal_mol_A']})
    return dict(physical=physical,electrical={k:v['ceiling']-v['score'] for k,v in metrics['fixed_electrical_safeguards'].items()},
                signed_W9_margin_Debye=metrics['signed_W9_z_safeguard']['absolute_ceiling_Debye']-abs(metrics['signed_W9_z_safeguard']['model_minus_reference_Debye']))


def run(mode,config_path,out):
    require(mode=='run','Scheduling recovery forbids scalar preparation replay')
    tick,cpu=time.monotonic(),time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN)
    cfg=json.loads(Path(config_path).read_text());cfg['fingerprint']=fingerprint(cfg);out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    ev=None;P=None;allocation=None
    result=dict(status='started',mode=mode,implementation=pin(__file__),configuration=pin(config_path),configuration_fingerprint=cfg['fingerprint'],limits=LIMITS,
                error=None,native_programs=0,solver_calls=0,new_derivative_probes=0,reference_calls=0,heldout_reads=0,selected_candidate=None,automatic_promotion=False,automatic_extension=False,job_id=None)
    def persist(terminal=False):
        if P is not None:P.persist_result(out,result,ev,terminal)
        else:save(out/'RESULT.tmp',result);os.replace(out/'RESULT.tmp',out/'RESULT.json')
        if terminal and ev is not None:
            # P.persist_result has completed the final authoritative materialization.
            # Pin afterward, then serialize without another snapshot/ledger mutation.
            if 'terminal_snapshot_error' not in result:
                result['journal_reconciliation']=dict(journal=pin(ev.cache_dir/'EVENTS.jsonl'),ledger=pin(ev.ledger_path),
                    final_authoritative_snapshot=True,new_vectors=len(ev.ledger['vector_keys']),new_native_programs=ev.ledger['native_program_invocations'],
                    returned_programs=ev.ledger['native_returncodes_observed'],imported_cells=0)
            else:
                result.pop('journal_reconciliation',None);result['complete_native_admission']=False
                result['witness_scientific_passed']=False;result['status']='failed_boundary_fraction_witness'
            save(out/'RESULT.tmp',result);os.replace(out/'RESULT.tmp',out/'RESULT.json')
    try:
        if mode=='prepare':require(not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID'),'One-thread login preparation')
        inp,fitcfg,producer,origin,proposal,linear=metadata(cfg,mode=='prepare')
        C,X,P,base,old,radial,checks,parent,policy,p=context(cfg,fitcfg,origin,linear)
        result.update(origin=inp['origin'],step8_proposal=inp['step8_proposal'],original_fit_stop=producer['stopping_reason'],solver_status_retained=proposal['solver_status'],solver_success=False,
                      scalar_authority=inp['linearization'],fixed_safeguards=policy.safety,incumbent_metrics=producer['selected_metrics'],incumbent_margins=margins(producer['selected_metrics']))
        require(p['raw_center_q']==producer['selected_metrics']['q'],'Literal incumbent q authority')
        if mode=='prepare':
            slot=certify(C,policy,p,proposal);save(out/'WITNESS.json',slot);result.update(witness=pin(out/'WITNESS.json'),witness_key=slot['key'],witness_parameters=slot['parameters'],eligible=slot['eligible'],
                status='prepared_boundary_fraction_witness_zero_calls' if slot['eligible'] else 'scalar_rejected_boundary_witness_unrun',old28_replayed=False,new_cache_constructed=False)
        else:
            require(str(out)==cfg['output_directory'],'Exact fresh execution directory');slot=release_gate(cfg)
            prep=read(cfg['preparation']);result.update(witness=prep['witness'],witness_key=slot['key'],witness_parameters=slot['parameters'],eligible=True)
            allocation=X.allocation(str(verify(cfg['scontrol'])));result.update(allocation=allocation,job_id=allocation['job_id']);persist()
            ev=A.evaluator_class(old,radial,base,checks,release_gate)(cfg,out/'cache',slot,parent,allocation);persist()
            answer=ev.evaluate(slot['parameters']);require(answer['status']=='complete' and len(answer['cells'])==46 and all(v['status']=='complete' for v in answer['cells'].values()),'Complete native-admitted cells')
            save(out/'VECTOR.json',answer);actual=policy.assess(answer);current=producer['selected_metrics'];tol=1e-10*max(1.,abs(current['q']));pred=slot['prediction']
            passed=actual['all_declared_safeguards_pass'] and actual['q']<current['q']-tol
            result.update(status='completed_boundary_fraction_native_test',vector=pin(out/'VECTOR.json'),actual_metrics=actual,actual_margins=margins(actual),actual_q_improvement=current['q']-actual['q'],
                          strict_q_improvement_tolerance=tol,witness_scientific_passed=passed,complete_native_admission=True,
                          actual_minus_predicted=dict(q=actual['q']-pred['actual_component_q'],component_errors_kcal_mol={k:actual['component_errors_kcal_mol'][k]-v for k,v in pred['component_errors_kcal_mol'].items()},
                                                      signed_W9_residual_Debye=actual['signed_W9_z_safeguard']['model_minus_reference_Debye']-pred['signed_W9_z_residual_Debye']))
            require(ev.ledger['native_program_invocations']==56 and ev.ledger['native_returncodes_observed']==56 and len(ev.ledger['vector_keys'])==1,'Exact one-vector56 completed attempts')
    except Exception as exc:
        result.update(status='failed_boundary_fraction_witness' if mode=='run' else 'failed_preparation',error=repr(exc))
        if ev is not None:
            try:ev._stop(exc)
            except Exception as other:result['stop_reporting_error']=repr(other)
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN)
        result['cost']=dict(wall_seconds=time.monotonic()-tick,controller_CPU_seconds=time.process_time()-cpu,child_user_CPU_seconds=ch.ru_utime-child.ru_utime,child_system_CPU_seconds=ch.ru_stime-child.ru_stime,
                            execution_mode='slurm' if allocation else 'login_preparation_or_preallocation_failure',allocation=allocation,scheduler_allocated_CPU_seconds=None if allocation else 0,
                            parent_scheduler_accounting_pending=bool(allocation),native_workers=11 if allocation else 0,numerical_threads=1,inherited_costs_not_added=True,per_receipt_child_CPU_not_summed=True)
        persist(True)
    print(json.dumps(dict(status=result['status'],error=result['error'],native_programs=result['native_programs'],solver_calls=0)));return result

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('mode',choices=('run',));a.add_argument('--config',required=True);a.add_argument('--out',required=True);args=a.parse_args()
    raise SystemExit(bool(run(args.mode,args.config,args.out)['error']))
