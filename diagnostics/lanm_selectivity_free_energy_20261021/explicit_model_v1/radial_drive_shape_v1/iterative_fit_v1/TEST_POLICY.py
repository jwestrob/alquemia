"""Manufactured changed-path controls. No molecular vector or solver execution."""
import concurrent.futures
import copy
import json
import os
import socket
import subprocess
import threading
import time
import types
from pathlib import Path
from unittest.mock import patch
import numpy as np
import clarabel
import DRIVER as R
import CORE as C
import ADAPTER as A


def run(out):
    out=Path(out).resolve(); out.mkdir(parents=True,exist_ok=False); tick,cpu=time.monotonic(),time.process_time(); checks=[]
    def check(name,ok): R.require(ok,name); checks.append(name)
    def forbid(*a,**kw): raise AssertionError('Real subprocess/solver forbidden')
    cfg0=json.loads((R.HERE/'CONFIG_DRAFT.json').read_text()); parentcfg=R.read(cfg0['radial_configuration'])
    old=R.module(parentcfg['original_adapter'],'test_campaign_old'); real=R.module(parentcfg['parent_evaluator'],'test_campaign_base'); radial=R.module(parentcfg['adapter'],'test_radial_base')
    domain={'parameters':[dict(lower=-1.,upper=1.) for _ in range(14)]}
    with patch.object(subprocess,'run',forbid),patch.object(subprocess,'Popen',forbid),patch.object(clarabel,'DefaultSolver',forbid):
        v=[.2]*12+[0.,0.]; points=C.probe_points(v,domain)
        check('literal_zero_radial_probe_plus_minus',points[24][12]==-.01 and points[25][12]==.01 and points[26][13]==-.01 and points[27][13]==.01)
        v[12:]=[.43,-.998]; points=C.probe_points(v,domain)
        check('nonzero14_probe_order_and_unchanged_coordinates',len(points)==28 and all(all(p[k]==v[k] for k in range(14) if k!=i//2) for i,p in enumerate(points)))
        check('nonzero_radial_shift_and_clipped_actual_span',points[24][12]<v[12]<points[25][12] and points[26][13]==-1. and 0<(points[27][13]-points[26][13])/2<.01)
        policy=object.__new__(C.Policy)
        safety=dict(actual_scores={k:1. for k in C.ELECTRICAL},literal_squared_norm_denominators={k:.25 for k in C.ELECTRICAL},W9_reference_Debye=-.1,W9_absolute_ceiling_Debye=.2)
        policy.safety=safety; policy.fixed={}; policy.targets={}
        policy.core=types.SimpleNamespace(assess=lambda v,f:copy.deepcopy(v['metrics']),token=R.digest)
        def vector(q=1.,score=1.,signed=.2):
            params=[q]+[0.]*13
            return dict(parameters=params,response={'C_dipoles_Debye':{'aquo9_central_Dy':[0.,0.,signed-.1],'aquo9_central_La':[0.,0.,0.]}},
                        metrics=dict(key=R.digest(params),q=q,original_block_scores={k:score for k in C.ELECTRICAL},all25_original_physical_gates=True,physical_constraints_pass=True,full_development_target=False))
        current=vector(); electric_bad=policy.assess(vector(.5,1.01))
        check('physical25_meaning_not_overwritten_by_electrical_failure',electric_bad['physical_constraints_pass'] and not electric_bad['all_declared_safeguards_pass'])
        check('fixed_signed_guard_rejects_actual_deterioration',not policy.assess(vector(.5,signed=.2001))['all_declared_safeguards_pass'])
        def slots(vs):
            return [dict(fraction=f,parameters=v['parameters'],key=R.digest(v['parameters']),eligible=True,execution_status='complete',native_admitted=True,vector=v) for f,v in zip(C.FRACTIONS,vs)]
        result=policy.select(current,slots([vector(),vector(1.-5e-11),vector(.9),vector(.9)]))
        check('incumbent_epsilon_and_earlier_fraction_tie_policy',result['selected_trial_index']==2 and result['strict_q_improvement_tolerance']==1e-10)
        check('electrical_failure_cannot_win_q_selection',policy.select(current,slots([vector(.1,1.01),vector(.9),vector(),vector()]))['selected_trial_index']==1)
        incomplete=slots([vector(.8)]*4); incomplete[-1]['execution_status']='started'
        try: policy.select(current,incomplete)
        except ValueError: pass
        else: raise AssertionError('Incomplete batch accepted')
        check('incomplete_batch_never_selected',True)
        p={'blocks':{k:{'denominator':9.} for k in C.ELECTRICAL},'signed_W9_z':{'absolute_ceiling_Debye':99.,'reference_Dy_minus_La_z_Debye':-.1}}
        policy.M=types.SimpleNamespace(problem=lambda *a:copy.deepcopy(p)); policy.attr=None
        fixed=policy.problem(None,None)
        check('later_current_norms_and_signed_ceiling_never_reset',all(fixed['blocks'][k]['denominator']==.25 for k in C.ELECTRICAL) and fixed['signed_W9_z']['absolute_ceiling_Debye']==.2)
        response=dict(retained_exact_sign_feasible=True,retained_strict_positive=False,all25_predicted=True,signed_W9_safeguard_pass=True,actual_component_q=.8,parameters=[0.]*14,key=R.digest([0.]*14),physical_fraction_changed=False)
        policy.S=types.SimpleNamespace(certificate=lambda *a:copy.deepcopy(response))
        cert=policy.fraction({'raw_center_q':1.},None,.5,False)
        check('exact_zero_margin_witness_allowed_without_solver_success',cert['eligible'] and cert['solver_success'] is False)
        response['retained_exact_sign_feasible']=False
        check('negative_margin_not_relaxed_to_one_e_minus8',not policy.fraction({'raw_center_q':1.},None,.5,False)['eligible'])

        # Real public evaluator -> real journal request -> real invocation intent;
        # controlled native helper raises before any external process can start.
        files=out/'manufactured'; files.mkdir()
        for n in ('x','p','o'): (files/n).write_text('manufactured '+n+'\n')
        pp,xyz=R.pin(files/'p'),R.pin(files/'x')
        rows={f'{m}_{i:02d}':dict(metal=m,needs_gradient=i<5,needs_export=False,fixed_POLPAIR=0.,source_input=dict(xyz=xyz,parameters=pp)) for m in ('La','Dy') for i in range(23)}
        calls=[]
        class ControlledStop(RuntimeError): pass
        def invoke(*a): calls.append(a[0]); raise ControlledStop('Stop before subprocess')
        eff=lambda v,m,d:list(v[:6] if m=='La' else v[6:12])
        def render(v,row,d,resolver,exports=None):
            x=eff(v,row['metal'],d); return 'MOCK '+json.dumps(x)+'\n',dict(physical_values_a_b_c_d_R_epsilon=x)
        class FakeBase:
            def __init__(self,cfg,path,root_release=None):
                self.config=cfg; self.root_release=root_release; self.cache_dir=Path(path); self.cache_dir.mkdir(); self.ledger_path=self.cache_dir/'LEDGER.json'
                self.m=dict(rows=rows,units={}); self.domain={}; self.family={}; self.initial_vector=[0.]*12
                self.lock=threading.RLock(); self.evaluation_lock=threading.RLock(); self.slots=threading.Semaphore(1); self.stop_event=threading.Event(); self.cell_locks={}; self.cache={}; self._gated=False
                self.h=dict(resolver=None,native_helper=types.SimpleNamespace(invoke=invoke)); self.ledger=dict(native_program_invocations=0,vector_keys=[],requests=[],stopped=False)
            def render(self,v,cid,exports=None): return render(v,rows[cid],{},None,exports)
            cache_identity=real.Evaluator.cache_identity
            def snapshot(self): return dict(manufactured=True,native_program_invocations=self.ledger['native_program_invocations'])
        base=types.SimpleNamespace(**vars(real)); base.Evaluator=FakeBase; base.render=render; base.effective=eff; base.checked=lambda v,d:list(v)
        entries=[dict(key=R.digest([i/1000.]*14),parameters=[i/1000.]*14) for i in range(31)]
        inv=dict(entries=entries,origin_key=entries[0]['key'])
        parent={k:pp for k in real.HELPERS}; parent.update(implementation=pp,source_manifest=pp,executables=dict(analyze=pp,testgrad=pp))
        cfg=dict(fingerprint='manufactured',adapter=R.pin(A.__file__),implementation=R.pin(R.__file__),parent_evaluator=parentcfg['parent_evaluator'],run_id='manufactured',
                 native_source_manifest=pp,executables=parent['executables'],build=pp,qualification_review=pp,limits=copy.deepcopy(R.LIMITS),threads=list(R.THREADS),hostname=socket.gethostname(),parent_configuration=pp,reuse_inventory=pp)
        gates=[]; E=A.evaluator_class(old,radial,base,types.SimpleNamespace(),lambda c:gates.append(True))
        pool=concurrent.futures.ThreadPoolExecutor
        class OneWorker(pool):
            def __init__(self,*a,**kw): super().__init__(max_workers=1)
        env=dict(SLURM_JOB_ID='123',SLURM_NTASKS='1',SLURM_JOB_NUM_NODES='1',SLURM_CPUS_PER_TASK='13',SLURM_RESTART_COUNT='0')
        text='JobId=123 JobState=RUNNING Account=banfield Partition=standard QOS=standard_expand NumNodes=1 NumTasks=1 NumCPUs=13 CPUs/Task=13 Requeue=0'
        allocation=R.X.parse(text,env)
        check('exact13CPU_scheduler_readback',allocation['allocated_CPUs']==13 and allocation['job_id']=='123')
        try: R.X.parse(text.replace('NumCPUs=13','NumCPUs=12'),env)
        except ValueError: pass
        else: raise AssertionError('Wrong grant accepted')
        check('wrong_scheduler_grant_rejected',True)
        with patch.object(old.concurrent.futures,'ThreadPoolExecutor',OneWorker),patch.dict(os.environ,env):
            ev=E(cfg,out/'cache','run',inv,parent,allocation)
            check('limits_committed_before_first_request',ev.persistence.state.ledger['limits']==R.LIMITS and ev.ledger['limits'] is not cfg['limits'])
            v=[.123]*14; ev.authorize(v,'probe',2)
            try: ev.evaluate(v)
            except ControlledStop: pass
            else: raise AssertionError('Boundary stop required')
            events=[json.loads(line)['event'] for line in (ev.cache_dir/'EVENTS.jsonl').read_text().splitlines()]
            check('public_evaluate_reaches_one_durable_attempt',sum(e['kind']=='request_start' for e in events)==1 and sum(e['kind']=='attempt_intent' for e in events)==1 and len(calls)==1)
            ev.snapshot()
            check('failed_attempt_preserves46slots_and_journal_identity',len(ev.ledger['requests'][0]['cells'])==46 and ev.ledger['requests'][0]['status']=='failed' and ev.persistence.state.ledger==ev.ledger)
            ev.stop_event.clear(); ev.ledger['native_program_invocations']=7280
            try: ev._gate()
            except ValueError: pass
            else: raise AssertionError('7280 cap bypass')
            check('invocation_gate_rejects_7280_without_second_helper_call',len(calls)==1)
            ev2=E(cfg,out/'import_failure_cache','run',inv,parent,allocation)
            def bad_import(*a): raise ValueError('manufactured inherited identity mismatch')
            ev2._import_saved=bad_import; ev2.authorize(entries[0]['parameters'],'probe',1)
            try: ev2.evaluate(entries[0]['parameters'])
            except ValueError as exc: check('saved_metal_identity_fails_before_request_or_invocation','inherited identity' in str(exc) and not ev2.ledger['requests'] and ev2.ledger['native_program_invocations']==0 and len(calls)==1)
            else: raise AssertionError('Unverified metal import accepted')
        # Failure writer is the admitted v2 function, with no journal reconciliation.
        M=R.module(parentcfg['math'],'test_persist_math'); P=R.module(cfg0['radial_runner'],'test_persistence',{'ADAPTER':radial,'MATH':M})
        def bad_snapshot(): raise ValueError('manufactured journal mismatch')
        ev.snapshot=bad_snapshot
        apparent=dict(mode='run',status='bounded_radial14_iterative_fit_complete',error=None,evaluator={'stale':True})
        P.persist_result(out,apparent,ev,True)
        check('terminal_snapshot_failure_cannot_claim_success',apparent['error'] is not None and apparent['evaluator_snapshot_available'] is False and 'evaluator' not in apparent and apparent['status'].startswith('failed'))
    result=dict(status='passed_radial14_changed_policy_controls',implementation=R.pin(__file__),sources={k:R.pin(R.HERE/(k+'.py')) for k in ('CORE','ADAPTER','DRIVER','RESOURCE')},checks=checks,count=len(checks),native_calls=0,solver_calls=0,scientific_vectors_read=0,
                mock_invocation_boundary_attempts=1,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(dict(status=result['status'],checks=len(checks))))


if __name__=='__main__':
    import sys
    try: run(sys.argv[1])
    except Exception as exc:
        out=Path(sys.argv[1]); out.mkdir(parents=True,exist_ok=True); (out/'FAILURE.json').write_text(json.dumps(dict(error=repr(exc),native_calls=0,solver_calls=0),indent=2)+'\n'); raise
