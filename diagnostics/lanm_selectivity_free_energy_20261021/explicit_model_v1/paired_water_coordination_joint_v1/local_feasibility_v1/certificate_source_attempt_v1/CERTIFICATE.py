#!/usr/bin/env python3
"""Persist a sparse dual certificate from saved LP output; no solver rerun."""
import hashlib,json,os,time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[5]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_coordination_joint_v1'
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def main():
    start,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    path=W/'local_feasibility_v1/RESULT.json';result=json.loads(path.read_text());fit_path=Path(result['inputs']['fit']['path']);assert pin(fit_path)==result['inputs']['fit'];fit=json.loads(fit_path.read_text())
    assert result['classification']=='necessary_local_linear_system_infeasible_within_trust_box'
    entry=fit['iterations'][-1];r=np.array(fit['selected']['score']['weighted_residuals']);J=np.array(entry['Jacobian_normalized'])
    assert fit['vectors'][entry['origin_key']]['parameters']==fit['selected']['parameters']
    indices={'balanced_exchange/water_central':1,'balanced_exchange/aquo9_central':3}
    records=result['LP']['active_linear_constraints'];assert len(records)==2 and all(v['sign']==-1 and v['quantity'] in indices for v in records)
    # For these negative residual inequalities: -E_J delta -t <= E_0.
    A=np.array([np.r_[-np.sqrt(32)*J[indices[v['quantity']]],-1.] for v in records]);b=np.array([np.sqrt(32)*r[indices[v['quantity']]] for v in records]);y=np.array([v['dual_marginal'] for v in records])
    x=np.r_[result['LP']['predicted_point']['normalized_delta'],result['LP']['threshold_multiplier_t']];c=np.r_[np.zeros(6),1.]
    lb=np.r_[result['normalized_delta_lower'],0.];ub=np.r_[result['normalized_delta_upper'],np.inf]
    g=c-A.T@y;lower=np.zeros(7);upper=np.zeros(7)
    for i in range(7):
        if abs(x[i]-lb[i])<1e-12:lower[i]=g[i]
        elif np.isfinite(ub[i]) and abs(x[i]-ub[i])<1e-12:upper[i]=g[i]
    stationarity=c-A.T@y-lower-upper;dual=float(b@y+lb@lower+ub[:6]@upper[:6]);slack=b-A@x
    residuals={'stationarity_max_abs':float(np.max(np.abs(stationarity))),'dual_sign_violation':max(0.,float(np.max(y)),float(np.max(-lower)),float(np.max(upper))),
               'complementarity_max_abs':max(float(np.max(np.abs(slack*y))),float(np.max(np.abs(lower*(x-lb)))),float(np.max(np.abs(upper[:6]*(ub[:6]-x[:6])))))),
               'primal_minus_dual_gap':float(x[-1]-dual)}
    assert all(abs(v)<1e-7 for v in residuals.values()) and dual>1+1e-7
    out={'status':'saved_sparse_dual_certificate_passed','implementation':pin(__file__),'LP_result':pin(path),'fit':pin(fit_path),'inequality_convention':'A*x<=b, dual marginals y<=0; lower-bound marginals>=0 and upper-bound marginals<=0',
         'active_linear_constraints':records,'active_A':A.tolist(),'active_b':b.tolist(),'active_inequality_duals':y.tolist(),'lower_bound_marginals':lower.tolist(),'upper_bound_marginals':upper.tolist(),
         'normalized_delta_lower':lb[:6].tolist(),'normalized_delta_upper':ub[:6].tolist(),'dual_objective_lower_bound_on_any_necessary_LP_t':dual,'residuals':residuals,
         'scope':'Local selected-point linear model and trust/original bounds only. This certificate uses only two energy inequalities; all other necessary inequality marginals are zero.',
         'new_solver_calls':0,'new_native_calls':0,'numerical_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    p=W/'local_feasibility_v1/SPARSE_DUAL_CERTIFICATE.json';assert not p.exists();p.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'certificate':pin(p),'dual_objective':dual,'residuals':residuals,'wall_seconds':out['wall_seconds'],'CPU_seconds':out['CPU_seconds']}))
if __name__=='__main__':main()
