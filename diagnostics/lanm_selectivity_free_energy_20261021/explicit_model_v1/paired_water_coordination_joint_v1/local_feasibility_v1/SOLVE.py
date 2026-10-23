#!/usr/bin/env python3
"""One saved-Jacobian local gate-feasibility diagnostic; no molecular imports."""
import hashlib
import json
import os
from pathlib import Path
import time
import warnings
import numpy as np
import scipy
from scipy.optimize import linprog, minimize

ROOT=Path(__file__).resolve().parents[5]
REL=Path('lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_coordination_joint_v1')
W=ROOT/'workspaces'/REL
GEOMETRIES=('water_minus','water_central','water_plus','aquo8_central','aquo9_central')
SIZES=(9,9,9,51,57)
NUMERICAL_TOL=1e-7

def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def run():
    wall,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    paths={'fit':W/'run_v1/execution/RESULT.json','raw_review':W/'selected_review_v1/run_v1/execution/REVIEW.json','independent_residual_review':W/'independent_residual_review_v1/RESULT.json'}
    expected={'fit':'5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832','raw_review':'fd7ea77167c833c3990d7624ff136042ce03f3ccf9dc4169ba7cff2aa51bb6a2','independent_residual_review':'c8c0535469e8e49c013307b3239762d7a4d3f300d84dcf370332dc553864acdc'}
    pins={k:pin(p) for k,p in paths.items()}
    assert all(pins[k]['sha256']==expected[k] for k in pins)
    fit,review,previous=(json.loads(paths[k].read_text()) for k in ('fit','raw_review','independent_residual_review'))
    proof=[]
    def check(name,value):
        assert value,name
        proof.append(name)
    check('Exact successful raw-review join',review['result']==pins['fit'] and review['status']=='selected_ten_raw_review_passed' and review['all_selected_raw_replays_exact'] and review['exact_selected_score_replayed'])
    check('Prior independent local review join',previous['inputs']['fit']==pins['fit'] and previous['status']=='independent_saved_interpretation_passed')
    selected=fit['selected'];entry=fit['iterations'][-1];records=fit['vectors'];origin=records[entry['origin_key']];score=selected['score']
    check('Selected point is exact final Jacobian origin',selected['parameters']==origin['parameters'] and selected['score']==origin['score'])
    lo=np.array([p['lower'] for p in fit['family']['parameters']]);hi=np.array([p['upper'] for p in fit['family']['parameters']]);width=hi-lo
    p0=np.array(selected['parameters']);z=(p0-lo)/width
    lower=np.maximum(-.025,-z);upper=np.minimum(.025,1-z)
    check('Valid original and trust-box intersection',np.all(width>0) and np.all(lower<=0) and np.all(upper>=0))
    r=np.array(score['weighted_residuals']);columns=[]
    for j in range(6):
        minus,plus=(records[k] for k in entry['derivative_keys'][2*j:2*j+2])
        a,b=np.array(minus['normalized_parameters']),np.array(plus['normalized_parameters']);span=b[j]-a[j];mask=np.arange(6)!=j
        check('Isolated derivative and recorded span '+str(j),span==entry['derivative_actual_widths'][j] and span>0 and np.allclose(a[mask],z[mask],atol=1e-14,rtol=0) and np.allclose(b[mask],z[mask],atol=1e-14,rtol=0))
        columns.append((np.array(plus['score']['weighted_residuals'])-np.array(minus['score']['weighted_residuals']))/span)
    J=np.column_stack(columns)
    check('Exact413 by6 Jacobian and finite residual',J.shape==(413,6) and np.array_equal(J,np.array(entry['Jacobian_normalized'])) and r.shape==(413,) and np.isfinite(J).all() and np.isfinite(r).all())
    E=np.zeros((10,413));E[:8,:8]=np.eye(8)*np.sqrt(32);E[8]=E[6]-E[4];E[9]=E[7]-E[5]
    enames=['balanced_exchange/'+q for q in ('water_minus','water_central','water_plus','aquo9_central')]+['individual_radial_work/'+m+'_'+s for m in ('La','Dy') for s in ('minus','plus')]+['paired_radial_work/'+s for s in ('minus','plus')]
    eb,EA=E@r,E@J
    check('Exact ten raw energy rows and paired signs',all(abs(eb[i]-score['energy_controls'][group][key]['error_kcal_mol'])<1e-12 for i,name in enumerate(enames) for group,key in [name.split('/')]))
    groups=[];scales=[];component_labels=[];start=8
    for kind in ('La','Dy','paired'):
        for q,size in zip(GEOMETRIES,SIZES):
            nw=(size-3)//6
            scale=np.array([np.sqrt(480)]*3+[np.sqrt(960*nw)]*(size-3))
            segment=slice(start,start+size);raw=scale*r[segment];record=score['per_fixture'][kind][q]
            check('Raw body component reconstruction '+kind+'/'+q,np.allclose(raw,record['raw_body_components_kcal_mol_A'],atol=2e-12,rtol=0))
            check('Independent RMS/max reconstruction '+kind+'/'+q,abs(np.sqrt(80)*np.linalg.norm(r[segment])-record['rigid_RMS_kcal_mol_A'])<1e-12 and abs(np.max(np.abs(raw))-record['maximum_body_component_error_kcal_mol_A'])<2e-12)
            groups.append({'name':kind+'/'+q,'slice':segment,'scale':scale});scales.extend(scale.tolist())
            component_labels.extend([kind+'/'+q+'/component_'+str(i) for i in range(size)]);start+=size
    scales=np.array(scales)
    check('Full unchanged25-gate denominator',start==413 and len(groups)==15 and len(eb)==10 and len(scales)==405)
    cb=scales*r[8:];CA=scales[:,None]*J[8:,:]
    offset=np.r_[eb,cb/5];matrix=np.vstack((EA,CA/5));labels=enames+component_labels
    # Minimize common threshold multiplier; +/- rows avoid an abs kink.
    Aub=np.vstack((np.c_[matrix,-np.ones(415)],np.c_[-matrix,-np.ones(415)]));bub=np.r_[-offset,offset]
    objective=np.r_[np.zeros(6),1.];bounds=list(zip(lower,upper))+[(0,None)]
    solver_start=time.monotonic()
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter('always')
        lp=linprog(objective,A_ub=Aub,b_ub=bub,bounds=bounds,method='highs-ds',options={'threads':1,'maxiter':10000,'time_limit':30.,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    solver_wall=time.monotonic()-solver_start
    def assess(delta):
        v=r+J@delta;energy=E@v;energy_checks={name:{'predicted_error_kcal_mol':float(e),'passed':bool(abs(e)<=1)} for name,e in zip(enames,energy)}
        force={}
        for g in groups:
            segment=v[g['slice']];raw=g['scale']*segment;rms=np.sqrt(80)*np.linalg.norm(segment);maximum=np.max(np.abs(raw))
            force[g['name']]={'predicted_RMS_kcal_mol_A':float(rms),'predicted_maximum_body_component_kcal_mol_A':float(maximum),'passed':bool(rms<=2 and maximum<=5),'predicted_body_components_kcal_mol_A':raw.tolist()}
        ratios=[abs(x) for x in energy]+[v for g in force.values() for v in (g['predicted_RMS_kcal_mol_A']/2,g['predicted_maximum_body_component_kcal_mol_A']/5)]
        return {'normalized_delta':delta.tolist(),'predicted_parameters':(p0+width*delta).tolist(),'within_original_and_trust_bounds':bool(np.all(delta>=lower-1e-12) and np.all(delta<=upper+1e-12)),'energy':energy_checks,'force':force,'energy_passed':sum(x['passed'] for x in energy_checks.values()),'force_passed':sum(x['passed'] for x in force.values()),'energy_total':10,'force_total':15,'all25_predicted_gates_pass':all(x['passed'] for x in energy_checks.values()) and all(x['passed'] for x in force.values()),'maximum_normalized_full_gate_ratio':float(max(ratios)),'predicted_weighted_loss':float(v@v),'native_evaluation_performed':False}
    lpr={'status':int(lp.status),'success':bool(lp.success),'message':str(lp.message),'iterations':int(lp.nit),'solver_wall_seconds':solver_wall,'warnings':[str(x.message) for x in captured],'method':'highs-ds','options':{'threads':1,'maxiter':10000,'time_limit':30,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9},'optimality_residuals':None,'predicted_point':None}
    certified=False
    if lp.x is not None:
        x=np.array(lp.x);lpr.update(threshold_multiplier_t=float(x[-1]),normalized_slack_over_original_threshold=float(x[-1]-1),predicted_point=assess(x[:6]))
        if lp.success:
            dual=np.array(lp.ineqlin.marginals);lm=np.array(lp.lower.marginals);um=np.array(lp.upper.marginals)
            lower7=np.r_[lower,0.];upper6=upper
            stationarity=objective-Aub.T@dual-lm-um
            dual_objective=float(bub@dual+lower7@lm+upper6@um[:6])
            primal_violation=max(0.,float(np.max(Aub@x-bub)),float(np.max(lower7-x)),float(np.max(x[:6]-upper6)))
            dual_sign_violation=max(0.,float(np.max(dual)),float(np.max(-lm)),float(np.max(um)))
            gap=float(lp.fun-dual_objective)
            comp=max(float(np.max(np.abs(dual*(bub-Aub@x)))),float(np.max(np.abs(lm*(x-lower7)))),float(np.max(np.abs(um[:6]*(upper6-x[:6])))))
            residuals={'primal_max_violation':primal_violation,'dual_sign_max_violation':dual_sign_violation,'stationarity_max_abs':float(np.max(np.abs(stationarity))),'complementarity_max_abs':comp,'primal_objective':float(lp.fun),'dual_objective':dual_objective,'primal_minus_dual_gap':gap}
            certified=all(abs(v)<=NUMERICAL_TOL for k,v in residuals.items() if k not in ('primal_objective','dual_objective'))
            lpr.update(optimality_residuals=residuals,numerically_certified_optimal=certified)
            active=np.where(np.abs(bub-Aub@x)<1e-7)[0]
            lpr['active_linear_constraints']=[{'quantity':labels[i%415],'sign':1 if i<415 else -1,'dual_marginal':float(dual[i])} for i in active]
    norm={'status':'not_run','reason':None};classification='unresolved_solver_status'
    if lp.success and lp.x is not None:
        t=float(lp.x[-1]);point=lpr['predicted_point']
        if certified and t>1+NUMERICAL_TOL:
            classification='necessary_local_linear_system_infeasible_within_trust_box';norm['reason']='Certified necessary LP optimum exceeds1; adding RMS constraints cannot restore local feasibility.'
        elif point['all25_predicted_gates_pass']:
            classification='predicted_local_point_satisfies_all25_gates';norm['reason']='LP point already passes all exact projected energy/RMS/component gates.'
        elif t<=1+NUMERICAL_TOL and any(x['predicted_RMS_kcal_mol_A']>2 for x in point['force'].values()):
            # One convex-norm minimax attempt, not a new molecular fit.
            q0=np.r_[lp.x[:6],max(t,point['maximum_normalized_full_gate_ratio'])];state={'x':q0.copy(),'callbacks':0};ns=time.monotonic()
            def fun(x):return float(x[-1])
            def jac(x):return objective
            def constraints(x):
                v=r+J@x[:6]
                return np.r_[bub-Aub@x,[x[-1]-np.sqrt(20)*np.linalg.norm(v[g['slice']]) for g in groups]]
            def cj(x):
                v=r+J@x[:6];rows=[]
                for g in groups:
                    sl=g['slice'];length=np.linalg.norm(v[sl]);gradient=np.zeros(6) if length==0 else -np.sqrt(20)*(J[sl].T@v[sl])/length
                    rows.append(np.r_[gradient,1.])
                return np.vstack((-Aub,np.array(rows)))
            def callback(x):
                state['x']=x.copy();state['callbacks']+=1
                if time.monotonic()-ns>30:raise TimeoutError('Declared30s SLSQP callback limit')
            try:
                n=minimize(fun,q0,jac=jac,method='SLSQP',bounds=bounds,constraints={'type':'ineq','fun':constraints,'jac':cj},callback=callback,options={'maxiter':300,'ftol':1e-10,'disp':False})
                final=n.x;norm={'status':int(n.status),'success':bool(n.success),'message':str(n.message),'iterations':int(n.nit),'objective_evaluations':int(n.nfev),'threshold_multiplier_t':float(n.fun),'normalized_slack_over_original_threshold':float(n.fun-1),'max_constraint_violation':float(max(0.,-np.min(constraints(final)))),'multipliers':np.asarray(n.get('multipliers',[])).tolist()}
            except Exception as exc:
                final=state['x'];norm={'status':'exception','success':False,'error':repr(exc),'iterations_completed_callback':state['callbacks'],'threshold_multiplier_t':float(final[-1]),'max_constraint_violation':float(max(0.,-np.min(constraints(final))))}
            norm.update(solver_wall_seconds=time.monotonic()-ns,max_iterations=300,ftol=1e-10,callback_wall_limit_seconds=30,predicted_point=assess(final[:6]))
            classification='predicted_local_point_satisfies_all25_gates' if norm['predicted_point']['all25_predicted_gates_pass'] else 'norm_constrained_attempt_did_not_supply_all25_passing_point'
        else:
            classification='necessary_LP_boundary_or_uncertified_status_unresolved';norm['reason']='No certified infeasibility and no exact passing point; preserve numerical boundary.'
    else:norm['reason']='Necessary LP not successfully solved; no feasibility claim or retry.'
    check('Parent scientific artifacts unchanged',all(pin(paths[k])==pins[k] for k in pins))
    result={'status':'saved_local_feasibility_diagnostic_complete','classification':classification,'implementation':pin(__file__),'declaration':pin(Path(__file__).with_name('DECLARATION.md')),'inputs':pins,'source_proof_checks':proof,'source_proof_count':len(proof),'residual_count':413,'energy_gate_count':10,'force_gate_count':15,'raw_component_count':405,'linear_inequality_count':830,'selected_parameters':p0.tolist(),'original_parameter_lower':lo.tolist(),'original_parameter_upper':hi.tolist(),'normalized_delta_lower':lower.tolist(),'normalized_delta_upper':upper.tolist(),'LP':lpr,'conditional_norm_solve':norm,'scipy_version':scipy.__version__,'numpy_version':np.__version__,'numerical_certification_tolerance':NUMERICAL_TOL,'new_native_calls':0,'new_reference_outcomes_read':False,'candidate_changed':False,'physical_family_infeasibility_claimed':False,'auto_evaluation_released':False,'numerical_threads':1,'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu}
    out=W/'local_feasibility_v1';out.mkdir(parents=True,exist_ok=True);dest=out/'RESULT.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'result':pin(dest),'classification':classification,'proof_checks':len(proof),'LP':{k:v for k,v in lpr.items() if k not in ('predicted_point','active_linear_constraints')},'norm_status':norm['status'],'wall_seconds':result['wall_seconds'],'CPU_seconds':result['CPU_seconds']}))

if __name__=='__main__':run()
