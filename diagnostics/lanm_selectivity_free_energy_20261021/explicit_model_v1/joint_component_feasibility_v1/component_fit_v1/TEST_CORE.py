"""Bounded saved/mock controls: no real optimizer or native evaluator is called."""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
from unittest.mock import patch

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    if os.environ.get(name)!='1':raise RuntimeError('One numerical thread required: '+name)
import numpy as np
from scipy.optimize import OptimizeResult

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('component_core',HERE/'CORE.py')
C=importlib.util.module_from_spec(spec);spec.loader.exec_module(C)
checks=[]
def check(name,ok,detail=None):
    if not bool(ok):raise AssertionError(name)
    checks.append({'name':name,'passed':True,'detail':detail})
def reject(name,fn):
    try:fn()
    except (ValueError,KeyError):check(name,True);return
    raise AssertionError('Expected rejection: '+name)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):
    assert pin(p['path'])==p
    return json.loads(Path(p['path']).read_text())

def main(out):
    wall=time.monotonic();cpu=time.process_time()
    inputs=json.loads((HERE/'CORE_CONTROL_INPUTS.json').read_text())
    data={k:read(v) for k,v in inputs.items()}
    fixed=C.make_fixed(data['origin'],data['manifest'],data['decomposition'],data['domain'],data['declaration'])
    origin=data['origin'];a=C.assess(origin,fixed)
    check('saved_origin25_and_fixed4',a['physical_constraints_pass'] and len(a['original25'])==25 and abs(a['q']-1)<1e-12)
    check('original_initial_milestone_separate_from_quarter',set(a['original_initial_component_RMS'])=={'P','B'} and not a['substantive_recovery_milestone'] and not a['q_at_or_below_quarter'])
    iteration=data['fit']['iterations'][-1]
    probes=[read(data['fit']['vectors'][key]['artifact']) for key in iteration['derivative_keys']]
    p=C.linearize(origin,probes,fixed)
    check('exact24_probe_order_actual_spans',p['spans']==iteration['derivative_spans'] and [C.token(v['parameters']) for v in probes]==iteration['derivative_keys'])
    jerr=float(np.max(abs(p['Jacobian_original531']-np.asarray(iteration['Jacobian_normalized']))))
    cerr=float(np.max(abs(p['Jacobian_components']-np.asarray(data['saved_diagnostic']['origins'][-1]['Jacobian_components']))))
    check('saved531_and_component_J_replay',jerr<1e-10 and cerr<1e-10,{'531_max_abs':jerr,'component_max_abs':cerr})
    reject('reordered_probe_rejected',lambda:C.linearize(origin,[probes[1],probes[0]]+probes[2:],fixed))
    pilot={'status':'complete','parameters':data['pilot']['candidate']['parameters'],**data['pilot']['answer']}
    pilot_metrics=C.assess(pilot,fixed)
    check('admitted_pilot_exact25_fixedceilings_components',pilot_metrics['physical_constraints_pass'] and pilot_metrics['q']<a['q'])
    old=data['failed_proposal'];returned=np.r_[old['normalized_direction'],1-old['common_fractional_squared_norm_reduction']]
    mock_result=OptimizeResult(x=returned,success=False,status=8,message=old['solver']['message'],nit=77,nfev=329,njev=73)
    with patch.object(C,'minimize',return_value=mock_result) as mocked:
        proposal=C.propose(origin,probes,fixed)
        check('one_mocked_SLSQP_no_retry',mocked.call_count==1 and mocked.call_args.kwargs['options']==C.OPTIONS)
    check('failed_solver_status_preserved',proposal['solver']['success'] is False and proposal['solver']['status']==8)
    half=proposal['trial_slots'][1]
    check('known_strict_half_witness_identity',half['eligible'] and half['key']==data['half_review']['half_parameter_token'] and half['parameters']==data['half_review']['physical_half_parameters'])
    check('known_full_nearfeasible_failure_rejected',not proposal['trial_slots'][0]['eligible'])
    check('four_fixed_fractions_and_strict_residuals',[s['fraction'] for s in proposal['trial_slots']]==list(C.FRACTIONS) and all(not s['eligible'] or s['certificate']['all_residual_margins_strictly_positive'] for s in proposal['trial_slots']))
    x=np.r_[returned[:12]*.5,(p['q_current']+returned[-1])*.5]
    v,j,names=C.constraints(p,x)
    scalar=half['certificate']['residual_constraints'];err=max(abs(float(y)-scalar[n]) for n,y in zip(names,v))
    check('independent_scalar851_certificate',len(v)==851 and err<1e-10,err)
    # Central differences evaluate only the already-frozen affine/quadratic math.
    numerical=np.empty_like(j);eps=1e-6
    for k in range(13):
        delta=np.zeros(13);delta[k]=eps
        numerical[:,k]=(C.constraints(p,x+delta)[0]-C.constraints(p,x-delta)[0])/(2*eps)
    jdiff=float(np.max(abs(j-numerical)))
    check('all851_analytic_constraint_derivatives',jdiff<1e-7,jdiff)
    bad=returned.copy();bad[0]=2
    check('out_of_box_fraction_rejected',not C.certify_fraction(p,bad,1,True)['eligible'])
    reject('undeclared_fraction_rejected',lambda:C.certify_fraction(p,returned,.75,False))
    check('no_q_improvement_rejected',not C.certify_fraction(p,np.r_[np.zeros(12),p['q_current']],1,True)['eligible'])
    with patch.object(C,'minimize',side_effect=RuntimeError('synthetic solver error')) as mocked:
        failed=C.propose(origin,probes,fixed)
        check('solver_exception_preserved_no_retry',mocked.call_count==1 and failed['trial_slots']==[] and failed['status']=='proposal_solver_failure_without_qualified_return')
    # Selection-only mocks exercise policy; these are not molecular candidates.
    current={'parameters':[0.]*12,'metric':{'q':1.,'physical_constraints_pass':True,'key':C.token([0.]*12)}}
    def slots(qs,passed=None):
        result=[]
        for i,(fraction,q) in enumerate(zip(C.FRACTIONS,qs)):
            values=[float(i+1)]*12;key=C.token(values)
            vector={'parameters':values,'metric':{'q':q,'physical_constraints_pass':True if passed is None else passed[i],'key':key}}
            result.append({'fraction':fraction,'parameters':values,'key':key,'eligible':True,'execution_status':'complete','native_admitted':True,'vector':vector})
        return result
    with patch.object(C,'assess',side_effect=lambda v,f:v['metric']):
        choice=C.select(current,slots([.9,.8,.7,.6]),fixed)
        check('actual_minimax_selection',choice['selected_trial_index']==3)
        check('current_wins_ties',C.select(current,slots([1.,1.,1.,1.]),fixed)['selected_trial_index'] is None)
        check('earlier_fraction_wins_ties',C.select(current,slots([.8,.8,.8,.8]),fixed)['selected_trial_index']==0)
        check('physical_gate_failure_filtered',C.select(current,slots([.5,.8,.9,.95],[False,True,True,True]),fixed)['selected_trial_index']==1)
        incomplete=slots([.5,.8,.9,.95]);incomplete[-1]['execution_status']='failed'
        reject('no_selection_from_incomplete_batch',lambda:C.select(current,incomplete,fixed))
        wrong=slots([.5,.8,.9,.95]);wrong[0]['vector']['parameters']=[9.]*12
        reject('certified_vector_identity_required',lambda:C.select(current,wrong,fixed))
        reject('all_four_slots_required',lambda:C.select(current,slots([.5,.8,.9,.95])[:3],fixed))
    json.dumps([a,pilot_metrics,proposal,failed,choice],allow_nan=False)
    check('public_outputs_JSON_finite',True)
    result={'status':'passed_component_core_saved_mock_controls','implementation':pin(HERE/'CORE.py'),'controls_source':pin(__file__),
            'inputs':pin(HERE/'CORE_CONTROL_INPUTS.json'),'checks':checks,'passed':len(checks),'actual_optimizer_calls':0,'mock_optimizer_calls':2,
            'new_native_calls':0,'new_reference_calls':0,'new_parameter_evaluations':0,
            'cost':{'wall_seconds':time.monotonic()-wall,'process_CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'status':result['status'],'passed':len(checks),'cost':result['cost']}))

if __name__=='__main__':main(sys.argv[1])
