"""Focused two-point policy/reporting controls; no evaluator cache/native/solver."""
import argparse,copy,json,math,subprocess,threading,time,types
from pathlib import Path
from unittest.mock import patch
import RUN as R

def main(out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False);checks=[]
    def check(name,x):R.require(x,name);checks.append(name)
    def rejected(fn):
        try:fn()
        except (ValueError,KeyError):return True
        return False
    def forbid(*a,**kw):raise AssertionError('Native/solver forbidden')
    import clarabel
    cfg=json.loads((R.HERE/'CONFIG_TEMPLATE.json').read_text())
    with patch.object(subprocess,'run',forbid),patch.object(subprocess,'Popen',forbid),patch.object(clarabel,'DefaultSolver',forbid):
        declaration,parent,points=R.metadata(cfg)
        check('two_literal_half_points_join_saved_affine',len(points)==2 and [x['case'] for x in points]==list(R.CASES))
        check('zero_extension_matched_control',points[0]['parameters'][12:]==[0.,0.])
        check('solver_status_retained_not_success',all(x['solver_status']=='AlmostSolved' for x in points))
        changed=copy.deepcopy(cfg);changed['limits']['new_parameter_vectors']=3
        check('scope_expansion_rejected',rejected(lambda:R.metadata(changed)))
        changed=copy.deepcopy(cfg);changed['inherited_actual_spans'][7]=.01
        check('scientific_drift_rejected',rejected(lambda:R.metadata(changed)))
        base=R.module(cfg['parent_evaluator'],'two_native_test_base');old=R.module(cfg['original_adapter'],'two_native_test_old');E=R.A.factory(R.P.A)(old,base,types.SimpleNamespace(),forbid)
        ev=E.__new__(E);ev.domain=R.read(R.read(cfg['parent_configuration'])['domain']);ev.allowed={x['key']:{'parameters':x['parameters']} for x in points};ev.ledger={'vector_keys':[],'native_program_invocations':0};ev.required_programs=lambda v:56
        third=points[0]['parameters'][:];third[0]+=1e-6
        check('undeclared_third_vector_rejected_before_dispatch',rejected(lambda:ev.evaluate(third)))
        ev.ledger['vector_keys']=['one','two'];check('two_vector_cap_before_dispatch',rejected(lambda:ev.evaluate(points[0]['parameters'])))
        ev.ledger.update(vector_keys=[],native_program_invocations=57);check('whole56_vector_budget_before_dispatch',rejected(lambda:ev.evaluate(points[0]['parameters'])))
        ev.ledger['native_program_invocations']=112;ev.slots=threading.Semaphore(1);ev.lock=threading.RLock();ev.stop_event=threading.Event();ev.mode='run'
        check('actual_inherited_invoke_boundary112_cap',rejected(lambda:ev._invoke('analyze',out,{})))
        original={'original_block_scores':{'pair_mu':1.,'C_mu':1.,'chi':1.}}
        metrics={'original_block_scores':{'pair_mu':1.1,'C_mu':.8,'chi':.9},'all25_original_physical_gates':False}
        center={'response':{'C_dipoles_Debye':{'aquo9_central_Dy':[0.,0.,2.5],'aquo9_central_La':[0.,0.,2.]}}}
        answer={'response':{'C_dipoles_Debye':{'aquo9_central_Dy':[0.,0.,2.6],'aquo9_central_La':[0.,0.,2.]}}}
        prediction={'signed_W9_z_absolute_ceiling_Debye':math.fsum([2.5,-2.,-.2]),'signed_W9_z_residual_Debye':.25}
        helper=types.SimpleNamespace(summarize=lambda *a:{'metrics':copy.deepcopy(metrics),'W9_C_Dy_minus_La_Debye':{'reference':[0.,0.,.2]}});core=types.SimpleNamespace(assess=lambda *a:original)
        report=R.summarize(answer,prediction,core,{},None,center,helper)
        check('all_three_norms_reported_without_new_tolerance',set(report['electrical_norm_nonincrease'])=={'pair_mu','C_mu','chi'} and not report['electrical_norm_nonincrease']['pair_mu']['nonincrease'])
        check('signed_residual_failure_retained',not report['signed_W9_z_safeguard']['passed'])
        check('physical_and_electrical_failure_not_selected_or_admitted',not report['all_declared_safeguards_pass'] and report['native_candidate_admission'] is False)
    result={'status':'passed_two_fixed_radial_native_policy_controls','implementation':R.pin(__file__),'runner':R.pin(R.__file__),'adapter':R.pin(R.A.__file__),'configuration':R.pin(R.HERE/'CONFIG_TEMPLATE.json'),'checks':checks,'checks_passed':len(checks),'native_programs':0,'solver_calls':0,'actual_vectors_evaluated':0,'old46_replay':False,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    R.save(out/'RESULT.json',result);print(json.dumps(result))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);main(p.parse_args().out)
