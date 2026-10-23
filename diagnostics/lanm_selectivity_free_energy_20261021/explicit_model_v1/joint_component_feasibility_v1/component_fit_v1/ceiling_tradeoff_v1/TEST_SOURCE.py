#!/usr/bin/env python3
"""Narrow manufactured controls; never instantiate the real solver."""
import ast
import copy
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import numpy as np
import BASELINE as B
import TRADEOFF as T

HERE = Path(__file__).resolve().parent


def main(out):
    wall, cpu = time.monotonic(), time.process_time()
    out.mkdir(parents=True, exist_ok=False); checks=[]
    result={'status':'started','implementation':T.pin(HERE/'TRADEOFF.py'),
            'baseline_implementation':B.pin(HERE/'BASELINE.py'),'test_implementation':B.pin(__file__),
            'checks':checks,'solver_calls':0,'native_calls':0,'scientific_vectors_read':0,'error':None}
    def check(name, yes):
        if not yes: raise AssertionError(name)
        checks.append(name)
    def rejects(name, fn):
        try: fn()
        except (ValueError, TypeError, FileExistsError): checks.append(name); return
        raise AssertionError('Unexpected pass '+name)
    try:
        parent=B.module(B.pin(HERE.parent/'conic_saved_comparison_v1/TEST_SOURCE.py'),'manufactured_authority')
        old=B.module(B.pin(HERE.parent/'conic_joint_ceiling_removal_v1/DIAGNOSE.py'),'row_removal')
        def forbidden(*args, **kwargs): raise AssertionError('Real solver construction forbidden')
        T.clarabel.DefaultSolver=forbidden
        p=parent.artificial();full={'solution':{'x':[0.]*13},'solver_status':'AlmostSolved','solver_success':False,'returned_q':0.}
        cert=B.certify(p,full,parent.core)
        check('raw_epigraph_failure_retained',cert['raw_certificate']['certificate']['minimum_epigraph_margin']<0)
        check('only_q_changed',cert['raw_certificate']['parameters']==cert['certificate']['parameters'] and not cert['physical_vector_changed'])
        check('all851_and25_strict',cert['certificate']['eligible'] and len(cert['certificate']['certificate']['residual_constraints'])==851)
        check('declared_cushion_exact',cert['auxiliary_q']==cert['actual_component_rho']+1e-8)
        rejects('nonfinite_ratio_rejected',lambda:B.auxiliary_q([float('nan'),1.]))
        rejects('negative_ratio_rejected',lambda:B.auxiliary_q([-.1,1.]))
        bad=copy.deepcopy(p);bad['energy']['r'][0]=2.
        rejects('baseline_physical_failure_not_repaired',lambda:B.certify(bad,full,parent.core))
        bad=copy.deepcopy(p);bad['blocks']['old_EF']['r']*=100
        rejects('baseline_ceiling_failure_not_repaired',lambda:B.certify(bad,full,parent.core))
        raw=np.r_[np.zeros(12),cert['auxiliary_q']]
        w=T.witness(p,raw,1.,False,parent.core,old.REMOVED,cert['actual_component_rho'])
        check('cushion_does_not_create_gain',w['gain_vs_baseline_actual_rho']==0 and not w['feasible_gain_exceeds_numerical_comparison'])
        check('almost_status_can_have_strict_witness',w['strict_retained_feasible'] and not w['constructed_source_certificate']['certificate']['solver_success'])
        w=T.witness(p,np.zeros(13),1.,False,parent.core,old.REMOVED,cert['actual_component_rho'])
        check('raw_epigraph_failure_retained_but_auxiliary_witness_valid',w['raw_minimum_epigraph_margin']<0 and w['strict_retained_feasible'] and not w['physical_fraction_changed'])
        w=T.witness(p,np.zeros(13),.5,False,parent.core,old.REMOVED,cert['actual_component_rho'])
        check('fraction_reexpression_does_not_change_physical_point',w['raw_source_certificate']['parameters']==w['constructed_source_certificate']['parameters'])
        bad=copy.deepcopy(p);bad['energy']['r'][0]=1.;bad['energy']['J'][0]=0
        check('solved_status_cannot_waive_strict_margin',not T.witness(bad,raw,1.,True,parent.core,old.REMOVED,1.)['strict_retained_feasible'])
        bad=copy.deepcopy(p);bad['blocks']['old_EF']['r']*=100
        w=T.witness(bad,raw,1.,False,parent.core,old.REMOVED,1.)
        check('removed_ceiling_visible_and_not_physical_gate',w['strict_retained_feasible'] and w['removed_ceiling_ratios']['old_EF']>1)
        check('witness_never_native_admission',not w['candidate_admitted'] and not w['native_evaluated'])
        bad=copy.deepcopy(p);bad['energy']['r'][0]=2.
        check('retained_energy_violation_rejected',not T.witness(bad,raw,1.,False,parent.core,old.REMOVED,1.)['strict_retained_feasible'])
        bad=copy.deepcopy(p);next(iter(bad['forces'].values()))['r'][0]=10.
        check('retained_force_violation_rejected',not T.witness(bad,raw,1.,False,parent.core,old.REMOVED,1.)['strict_retained_feasible'])
        badx=raw.copy();badx[0]=.03
        check('trust_violation_rejected',not T.witness(p,badx,1.,False,parent.core,old.REMOVED,1.)['strict_retained_feasible'])
        full_model=parent.C.assemble(p); model,removed,keep=old.remove_joint_ceilings(full_model)
        check('literal535_row_removal_only',keep==list(range(535,1817)) and np.array_equal(model['A'].toarray(),full_model['A'].toarray()[535:]))
        check('17SOC_856_linear',len(model['cones'])==18 and model['A'].shape==(1282,13))
        tree=ast.parse((HERE/'TRADEOFF.py').read_text())
        check('one_solve_site',sum(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='solve' for n in ast.walk(tree))==1)
        check('source_review_and_release_do_not_change_fingerprint',T.fingerprint({'x':1,'source_review':None,'root_release':None})==T.fingerprint({'x':1,'source_review':'r','root_release':'s'}))
        check('baseline_pin_changes_fingerprint',T.fingerprint({'baseline_certificate':1})!=T.fingerprint({'baseline_certificate':2}))
        settings=SimpleNamespace(dummy=1); mocks=[]
        class MockSolver:
            def __init__(self,*args):mocks.append('construct')
            def get_settings(self):return settings
            def print_to_buffer(self):pass
            def solve(self):
                mocks.append('solve');x=np.r_[np.zeros(12),.5]
                return SimpleNamespace(x=x,s=model['b']-model['A']@x,z=np.zeros(1282),status='AlmostSolved')
            def get_info(self):return SimpleNamespace(linsolver=SimpleNamespace(name='mock'),status='AlmostSolved')
            def get_print_buffer(self):return 'manufactured mock only\n'
        T.clarabel.DefaultSolver=MockSolver
        T.context=lambda cfg:({},parent.C,old,parent.core,p,full_model,
            {'origin_key':'manufactured','derivative_keys':[], 'linearized_problem':{'manufactured':True}},settings,{'baseline':cert})
        cfg={'output_directory':str(out/'mock1'),'one_attempt_marker':str(out/'ONCE.json'),
             'inputs':{'manufactured':True},'baseline_certificate':{'manufactured':True}}
        # The reporting authority is a mock pin; no scientific input is loaded.
        original_context=T.context
        T.context=lambda cfg:({'full_result':{'manufactured':True}},*original_context(cfg)[1:])
        cp=out/'MOCK_CONFIG.json';cp.write_text(json.dumps(cfg));a=T.run(cp)
        check('mock_preserves_AlmostSolved',a['status']=='completed_saved_ceiling_tradeoff' and a['solver_status']=='AlmostSolved' and not a['solver_success'])
        check('mock_one_attempt_only',mocks==['construct','solve'] and a['solver_calls']==1)
        cfg['output_directory']=str(out/'mock2');cp=out/'MOCK_CONFIG2.json';cp.write_text(json.dumps(cfg));a=T.run(cp)
        check('persistent_marker_rejects_second_attempt',a['status']=='failed_saved_ceiling_tradeoff' and a['solver_calls']==0 and mocks==['construct','solve'])
        result['mock_solver_calls']=1;result['status']='passed_ceiling_tradeoff_manufactured_controls'
    except Exception as exc:
        result.update(status='failed_ceiling_tradeoff_manufactured_controls',error=repr(exc))
    result['check_count']=len(checks);result['cost']={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'result':B.pin(out/'RESULT.json'),'checks':len(checks),'error':result['error']}))
    return result


if __name__=='__main__':raise SystemExit(bool(main(Path(sys.argv[1]))['error']))
