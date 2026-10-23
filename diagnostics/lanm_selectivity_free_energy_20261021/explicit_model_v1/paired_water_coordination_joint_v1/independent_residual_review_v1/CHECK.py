#!/usr/bin/env python3
"""Independent SVD of the immutable selected-origin saved finite differences."""
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[5]
REL=Path('lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_coordination_joint_v1')
W=ROOT/'workspaces'/REL

def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    start,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    paths={'fit':W/'run_v1/execution/RESULT.json','raw_review':W/'selected_review_v1/run_v1/execution/REVIEW.json','summary':W/'residual_review_v1/actual_v1/RESULT.json'}
    expected={'fit':'5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832','raw_review':'fd7ea77167c833c3990d7624ff136042ce03f3ccf9dc4169ba7cff2aa51bb6a2','summary':'fe4ad9fa4748f307c4a2bf8825e505635fb2cd78d11ff7518f7a6ca92e70f68d'}
    pins={k:pin(p) for k,p in paths.items()}
    assert all(pins[k]['sha256']==expected[k] for k in pins)
    fit,raw,summary=(json.loads(paths[k].read_text()) for k in ('fit','raw_review','summary'))
    checks=[]
    def check(label,test):
        assert test,label
        checks.append(label)
    check('Exact terminal source-result joins',raw['result']==summary['parent_result']==pins['fit'])
    check('Successful raw review preserves scientific failure',raw['status']=='selected_ten_raw_review_passed' and raw['complete']==10 and raw['all_selected_raw_replays_exact'] and raw['exact_selected_score_replayed'] and raw['all_development_controls_passed']==fit['development_checks_passed']==False)
    entry=fit['iterations'][-1];records=fit['vectors'];origin=records[entry['origin_key']];selected=fit['selected']
    check('Last exact Jacobian origin is selected',origin['parameters']==selected['parameters'] and origin['score']==selected['score'] and summary['last_Jacobian']['origin_equals_selected'])
    lo=np.array([x['lower'] for x in fit['family']['parameters']]);hi=np.array([x['upper'] for x in fit['family']['parameters']]);width=hi-lo
    z=(np.array(selected['parameters'])-lo)/width
    check('Every selected parameter is strictly interior',bool(np.all((z>0)&(z<1))))
    columns=[];spans=[]
    for j in range(6):
        minus,plus=(records[k] for k in entry['derivative_keys'][2*j:2*j+2])
        a,b=np.array(minus['normalized_parameters']),np.array(plus['normalized_parameters'])
        mask=np.arange(6)!=j
        check('Derivative coordinate isolation '+str(j),bool(np.allclose(a[mask],z[mask],atol=1e-14,rtol=0) and np.allclose(b[mask],z[mask],atol=1e-14,rtol=0)))
        span=b[j]-a[j];spans.append(float(span))
        check('Derivative actual span '+str(j),span==entry['derivative_actual_widths'][j] and span>0)
        columns.append((np.array(plus['score']['weighted_residuals'])-np.array(minus['score']['weighted_residuals']))/span)
    J=np.column_stack(columns);r=np.array(origin['score']['weighted_residuals'])
    check('Exact stored 413 by6 Jacobian reconstruction',J.shape==(413,6) and np.array_equal(J,np.array(entry['Jacobian_normalized'])))
    U,s,Vt=np.linalg.svd(J,full_matrices=False);keep=s>s[0]*1e-6
    projected_change=-U[:,keep]@(U[:,keep].T@r)
    delta=-Vt[keep,:].T@((U[:,keep].T@r)/s[keep])
    rp=r+projected_change
    sj=summary['last_Jacobian'];origin_loss=float(r@r);projected_loss=float(rp@rp)
    check('Independent SVD rank/singular values',int(keep.sum())==sj['rank']==6 and np.allclose(s,sj['singular_values'],rtol=1e-12,atol=1e-12))
    check('Independent unconstrained direction agrees',np.allclose(delta,sj['unrestricted_normalized_direction'],rtol=1e-7,atol=1e-12))
    check('Independent projected loss agrees',abs(projected_loss-sj['unrestricted_linear_projection']['weighted_loss'])<1e-12)
    check('Small unconstrained direction is not clipped',np.max(np.abs(delta))<.025 and np.all((z+delta>0)&(z+delta<1)) and not sj['unrestricted_exceeds_trust'] and not sj['unrestricted_exceeds_parameter_box'])
    check('Stop occurred before all finite caps',fit['stopping_reason']=='no_actual_improvement_at_complete_trial_batch' and fit['completed_outer_steps']==22<24 and fit['new_parameter_vectors']==338<384 and fit['new_native_program_invocations']==4120<5120)
    trials=[]
    for slot,key in zip(entry['trial_slots'],entry['trial_keys']):
        value=records[key];dz=(np.array(value['parameters'])-np.array(origin['parameters']))/width
        pred=r+J@dz;loss=value['score']['objective']
        trials.append({'id':slot['id'],'predicted_loss':float(pred@pred),'actual_loss':loss,'actual_change_from_selected':loss-origin_loss,'development_passed':value['score']['development_checks_passed']})
    check('All four final trials complete and fail to improve',len(trials)==4 and all(t['actual_change_from_selected']>=0 and not t['development_passed'] for t in trials))
    energy_scale=np.sqrt(.25*.5/4)
    numeric={'origin_loss':origin_loss,'linear_projected_loss':projected_loss,'linear_loss_reduction':origin_loss-projected_loss,'relative_linear_loss_reduction':(origin_loss-projected_loss)/origin_loss,
             'retained_residual_norm_fraction':float(np.linalg.norm(rp)/np.linalg.norm(r)),
             'normalized_parameter_positions':z.tolist(),'closest_normalized_bound_distance':float(np.min(np.minimum(z,1-z))),
             'actual_derivative_spans':spans,'singular_values':s.tolist(),'condition':float(s[0]/s[-1]),'rank':int(keep.sum()),
             'unconstrained_normalized_direction':delta.tolist(),'maximum_absolute_direction':float(np.max(np.abs(delta))),'trust_box_halfwidth':.025,
             'J_transpose_residual_norm':float(np.linalg.norm(J.T@r)),'gradient_norm_2JTr':float(np.linalg.norm(2*J.T@r)),
             'actual_W9_exchange_error':float(r[3]/energy_scale),'linear_projected_W9_exchange_error':float(rp[3]/energy_scale),
             'W9_exchange_weighted_loss_contribution':float(r[3]**2),'W9_exchange_fraction_of_total_loss':float(r[3]**2/origin_loss),
             'actual_final_trials':trials,'selected_control_counts':{'energy_passed':summary['selected']['energy_passed_count'],'energy_total':10,'force_passed':summary['selected']['force_passed_count'],'force_total':15}}
    result={'status':'independent_saved_interpretation_passed','implementation':pin(__file__),'specification':pin(Path(__file__).with_name('SPEC.md')),'inputs':pins,'checks_passed':len(checks),'checks':checks,'numerical_findings':numeric,
            'conclusion':'No evidence that this terminal failure is a simple iteration-cap, parameter-bound or trust-clipping limitation. The saved selected-point finite-difference linearization has almost all weighted residual outside its local six-column span; all evaluated final trials fail to improve.',
            'limits':['Finite-difference Jacobian at normalized probe halfwidth0.005 is not an exact infinitesimal derivative.','Engineering residuals and paired/individual observations are correlated; this projection is not a statistical confidence interval.','First-order local incompatibility is not a proof that no distant or nonlinear parameter set can pass.','Microscopic noise and finite-difference error can affect the last1e-9 loss changes; none explains the remaining4.338 kcal/mol W9 exchange error.','No unique missing physical interaction follows from this parameter-space projection.'],
            'new_native_calls':0,'objective_or_raw_replay_repeated':False,'candidate_changed':False,'fresh_validation_outcomes_read':False,'numerical_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    out=W/'independent_residual_review_v1';out.mkdir(parents=True,exist_ok=True);dest=out/'RESULT.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'result':pin(dest),'checks_passed':len(checks),'findings':numeric,'wall_seconds':result['wall_seconds'],'CPU_seconds':result['CPU_seconds']}))

if __name__=='__main__':main()
