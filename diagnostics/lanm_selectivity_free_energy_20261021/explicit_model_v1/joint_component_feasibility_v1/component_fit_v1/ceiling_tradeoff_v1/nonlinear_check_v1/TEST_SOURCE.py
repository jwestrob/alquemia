#!/usr/bin/env python3
"""Focused candidate binding and physical-unit report controls; zero native calls."""
import ast
import copy
import json
import math
from pathlib import Path
import sys
import time

import RUN as R

HERE=Path(__file__).resolve().parent


def main(out):
    wall,cpu=time.monotonic(),time.process_time();out.mkdir(parents=True,exist_ok=False);checks=[]
    result={'status':'started','implementation':R.pin(HERE/'RUN.py'),'test_implementation':R.pin(__file__),
            'native_calls':0,'optimizer_calls':0,'checks':checks,'error':None}
    def check(name,yes):
        if not yes:raise AssertionError(name)
        checks.append(name)
    def reject(name,call):
        try:call()
        except (ValueError,KeyError):checks.append(name);return
        raise AssertionError('Unexpected acceptance '+name)
    try:
        cfg=json.loads((HERE/'CONFIG_PREPARE.json').read_text());actual=R.read(cfg['tradeoff_result']);policy=R.read(cfg['followup_policy'])
        w=R.policy_candidate(actual,policy)
        check('exact_predeclared_first_full_fraction',w['fraction']==1. and w['constructed_source_certificate']['parameters']==cfg['candidate']['parameters'])
        bad=copy.deepcopy(actual);bad['witnesses'][0]['strict_retained_feasible']=False
        reject('no_fallback_half_after_full_rejection',lambda:R.policy_candidate(bad,policy))
        bad=copy.deepcopy(actual);bad['witnesses'][0]['actual_component_rho']=bad['baseline_actual_component_rho']-.049
        reject('subthreshold_full_does_not_launch',lambda:R.policy_candidate(bad,policy))
        bad=copy.deepcopy(actual);bad['witnesses'][0]['constructed_source_certificate']['key']='changed'
        reject('different_parameter_key_rejected',lambda:R.policy_candidate(bad,policy))
        bad=copy.deepcopy(policy);bad['selection_rule']['order']=[.5,1.,.25,.125]
        reject('order_cannot_change',lambda:R.policy_candidate(actual,bad))
        old=R.verify(cfg['original_pilot_implementation']).read_text();new=(HERE/'RUN.py').read_text()
        def fun(src):return {x.name:ast.dump(x,include_attributes=False) for x in ast.walk(ast.parse(src)) if isinstance(x,ast.FunctionDef)}
        a,b=fun(old),fun(new)
        for name in ('allowed','_stop','_invoke','_cell','evaluate','snapshot','one'):
            check('unchanged_execution_or_cycle_'+name,a[name]==b[name])
        check('actual_review_gate_preserves_failed_solver_status',"actual['solver_success']==source['solver_success']" in new and "actual['certificate_checks_passed'] is True" in new)
        check('all_expected_actualreview_fields_present',all(k in new for k in ['prospective_policy','prospective_trigger_fraction','prospective_trigger_witness']))
        check('one_vector_56_attempts12_workers',R.LIMITS=={'new_parameter_vectors':1,'native_program_attempts':56,'workers':12,'aggregate_login_thread_budget':16})
        parent=R.read(cfg['parent_configuration']);origin=R.read(cfg['origin'])
        report=R.component_report(cfg,parent,origin,origin)
        check('saved_origin_cycle_identity',report['candidate']==report['origin'] and report['actual_q']==1.)
        check('old_ceilings_descriptive_only',not report['old_ceilings_required_for_reporting'] and not report['model_admitted'])
        details=origin['score']['response']['per_observable'];pair=[]
        for parent_id in sorted({v['parent'] for v in details['pair'].values()}):
            values=[x for row in details['pair'].values() if row['parent']==parent_id for x in row['error_Debye']]
            pair.append(math.fsum(x*x for x in values)/len(values))
        expected={'equal_parent_pair_cartesian_RMS_Debye':math.sqrt(math.fsum(pair)/4),
                  'four_parent_C_cartesian_RMS_Debye':math.sqrt(math.fsum(x*x for row in details['C'].values() for x in row['error_Debye'])/12),
                  'four_state_susceptibility_RMS_a0_cubed':math.sqrt(math.fsum(row['error_a0_cubed']**2 for row in details['chi'].values())/4)}
        for k,v in expected.items():check('literal_residual_physical_unit_'+k,abs(report['electrical_residual_RMS']['origin'][k]-v)<=1e-10*max(1,abs(v)))
        predicted=report['predicted_blocks'];order=('P8','B8','P9','B9')
        err=report['predicted_component_errors_kcal_mol']
        check('predicted_component_scalar_closure',all(abs(sum(err[k]**2 for k in order if k.startswith(ch))/2-predicted[ch+'_energy']['predicted_squared_norm'])<=1e-10*max(1,predicted[ch+'_energy']['predicted_squared_norm']) for ch in ('P','B')))
        check('candidate_pin_changes_fingerprint',R.fingerprint(cfg)!=R.fingerprint(dict(cfg,candidate={'changed':True})))
        check('runtime_review_bound_by_release_not_preparation_fingerprint',R.fingerprint(cfg)==R.fingerprint(dict(cfg,tradeoff_actual_review=None)))
        result['status']='passed_ceiling_tradeoff_native_source_controls'
    except Exception as exc:result.update(status='failed_ceiling_tradeoff_native_source_controls',error=repr(exc))
    result['check_count']=len(checks);result['cost']={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'result':R.pin(out/'RESULT.json'),'checks':len(checks),'error':result['error']}));return result


if __name__=='__main__':raise SystemExit(bool(main(Path(sys.argv[1]))['error']))
