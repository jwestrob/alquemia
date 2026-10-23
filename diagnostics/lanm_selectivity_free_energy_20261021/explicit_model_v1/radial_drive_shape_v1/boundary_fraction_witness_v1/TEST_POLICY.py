#!/usr/bin/env python3
"""Tiny public request/invocation and scalar-wiring controls; no scientific calls."""
import argparse,ast,concurrent.futures,copy,json,os,subprocess,threading,time,types
from pathlib import Path
from unittest.mock import patch
import RUN as R
import ADAPTER as A

def main(out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    cfg0=json.loads((R.HERE/'CONFIG_PREPARE.json').read_text());old=R.module(cfg0['original_adapter'],'test_original_adapter');real=R.module(cfg0['parent_evaluator'],'test_original_evaluator');radial=R.module(cfg0['radial_adapter'],'test_radial_adapter')
    checks=[]
    def check(name,ok):R.require(ok,name);checks.append(name)
    def rejected(name,fn):
        try:fn()
        except (ValueError,RuntimeError):check(name,True)
        else:raise AssertionError(name)
    def forbid(*a,**k):raise AssertionError('Actual subprocess/solver forbidden')
    class ControlledStop(RuntimeError):pass
    calls=[]
    def invoke(*args):calls.append(args[0]);raise ControlledStop('Manufactured pre-subprocess stop')
    f=out/'manufactured';f.mkdir();(f/'x').write_text('fake xyz\n');(f/'p').write_text('fake prm\n');pp=R.pin(f/'p');xyz=R.pin(f/'x')
    rows={f'{m}_{i:02d}':dict(metal=m,needs_gradient=i<5,needs_export=False,fixed_POLPAIR=0.,source_input=dict(xyz=xyz,parameters=pp)) for m in ('La','Dy') for i in range(23)}
    eff=lambda v,m,d:list(v[:6] if m=='La' else v[6:12])
    def render(v,row,d,h,exports=None):return 'FAKE\n',dict(physical_values_a_b_c_d_R_epsilon=eff(v,row['metal'],d))
    class FakeBase:
        def __init__(self,cfg,path,root_release=None):
            self.config=cfg;self.root_release=root_release;self.cache_dir=Path(path);self.cache_dir.mkdir();self.ledger_path=self.cache_dir/'LEDGER.json'
            self.m={'rows':rows,'units':{}};self.domain={};self.family={};self.initial_vector=[0.]*12
            self.lock=threading.RLock();self.evaluation_lock=threading.RLock();self.slots=threading.Semaphore(1);self.stop_event=threading.Event();self.cell_locks={};self.cache={};self._gated=False
            self.h={'resolver':None,'native_helper':types.SimpleNamespace(invoke=invoke)}
            self.ledger={'native_program_invocations':0,'vector_keys':[],'requests':[],'stopped':False}
        def render(self,v,cid,exports=None):return render(v,rows[cid],{},None,exports)
        cache_identity=real.Evaluator.cache_identity
        def snapshot(self):return dict(native_program_invocations=self.ledger['native_program_invocations'])
    base=types.SimpleNamespace(**vars(real));base.Evaluator=FakeBase;base.render=render;base.effective=eff;base.checked=lambda v,d:list(v)
    parent={k:pp for k in real.HELPERS};parent.update(implementation=pp,source_manifest=pp,executables=dict(analyze=pp,testgrad=pp))
    point=dict(parameters=[.2]*12+[.01,-.01]);point['key']=R.digest(point['parameters']);allocation=dict(mode='slurm',job_id='123',allocated_CPUs=13)
    cfg=dict(fingerprint='mock',adapter=R.pin(A.__file__),implementation=R.pin(R.__file__),parent_evaluator=cfg0['parent_evaluator'],run_id='mock',native_source_manifest=pp,executables=parent['executables'],build=pp,qualification_review=pp,limits=copy.deepcopy(R.LIMITS),threads=list(R.THREADS),parent_configuration=pp,reuse_inventory=pp)
    gates=[];E=A.evaluator_class(old,radial,base,types.SimpleNamespace(),lambda c:gates.append(True))
    import clarabel
    pool=concurrent.futures.ThreadPoolExecutor
    class Serial(pool):
        def __init__(self,*a,**kw):super().__init__(max_workers=1)
    env={k:'1' for k in R.THREADS};env.update(SLURM_JOB_ID='123',SLURM_CPUS_PER_TASK='13',SLURM_NTASKS='1')
    with patch.dict(os.environ,env),patch.object(subprocess,'run',forbid),patch.object(subprocess,'Popen',forbid),patch.object(clarabel,'DefaultSolver',forbid),patch.object(old.concurrent.futures,'ThreadPoolExecutor',Serial):
        ev=E(cfg,out/'cache',point,parent,allocation)
        first=json.loads((ev.cache_dir/'EVENTS.jsonl').read_text().splitlines()[0])['event']
        check('detached1_56_limits_precede_initialize',first['data']['ledger']['limits']==R.LIMITS and ev.ledger['limits'] is not cfg['limits'])
        check('fresh_no_import_exact46_56',not ev.cache and not ev.imported_cell_keys and ev.required_programs(point['parameters'])==56)
        bad=point['parameters'].copy();bad[-1]+=.00001
        rejected('near_point_rejected_before_request',lambda:ev.evaluate(bad));check('no_request_for_wrong_point',not ev.ledger['requests'])
        try:ev.evaluate(point['parameters'])
        except ControlledStop:pass
        else:raise AssertionError('Expected mock stop')
        events=[json.loads(x)['event'] for x in (ev.cache_dir/'EVENTS.jsonl').read_text().splitlines()]
        check('one_durable_request_and_attempt',sum(e['kind']=='request_start' for e in events)==1 and sum(e['kind']=='attempt_intent' for e in events)==1 and len(calls)==1)
        check('56_program_request_preserved',ev.ledger['requests'][0]['required_new_programs']==56 and len(ev.ledger['requests'][0]['cells'])==46)
        check('first_failure_stops_queue',ev.stop_event.is_set() and ev.ledger['native_program_invocations']==1 and ev.ledger['native_returncodes_observed']==0)
        check('failed_attempt_status_preserved',ev.program_attempts[0]['status']=='failed_native_attempt' and ev.ledger['requests'][0]['status']=='failed')
        snap=ev.snapshot();check('journal_exact_after_failure',ev.persistence.state.ledger==ev.ledger and snap['max_parameter_vectors']==1 and snap['max_native_programs']==56)
        rejected('failed_vector_cannot_retry',lambda:ev.evaluate(point['parameters']))
        cell=next((ev.cache_dir/'cells').glob('*/RESULT.json'));check('cell_actual_scheduler_metadata',json.loads(cell.read_text())['job_id']=='123')
        ev2=E(cfg,out/'cap_cache',point,parent,allocation);ev2.ledger['native_program_invocations']=56
        rejected('invoke56_ceiling_before_helper',lambda:ev2._invoke('analyze',f,{}));check('no57th_attempt_or_helper',not ev2.program_attempts and len(calls)==1)
        class Policy:
            def fraction(self,p,pt,fraction,success):
                self.received=(pt.tolist(),fraction,success);v=[.1]*14
                return dict(parameters=v,key=R.digest(v),fraction=fraction,solver_success=success,eligible=True)
        policy=Policy();R.certify(None,policy,{},dict(solution={'x':[.2]*14+[.3]}))
        check('fraction_and_q_square_not_t',policy.received==([.2]*14+[.09],.0625,False))
        rejected('nonfinite_saved_solution_rejected',lambda:R.certify(None,policy,{},dict(solution={'x':[float('nan')]*15})))
        rejected('negative_t_rejected',lambda:R.certify(None,policy,{},dict(solution={'x':[.2]*14+[-.3]})))
        bounded=A.bounded_parent(old);P=bounded.evaluator_class(types.SimpleNamespace(**vars(base)))
        check('inherited_native_and_scientific_methods',E._invoke.__code__.co_code==P._invoke.__code__.co_code and E._parse.__code__.co_code==radial.evaluator_class(bounded,base,types.SimpleNamespace(),lambda c:None)._parse.__code__.co_code)
        X=R.module(cfg0['resource_helper'],'test_resource')
        text='JobId=123 JobState=RUNNING Account=banfield Partition=standard QOS=standard_expand NumNodes=1 NumTasks=1 NumCPUs=13 CPUs/Task=13 Requeue=0'
        env2=dict(env,SLURM_JOB_NUM_NODES='1');check('real_allocation_parser13',X.parse(text,env2)['allocated_CPUs']==13)
        rejected('12cpu_allocation_rejected',lambda:X.parse(text.replace('NumCPUs=13','NumCPUs=12'),env2))
    result=dict(status='passed_boundary_fraction_manufactured_controls',implementation=R.pin(__file__),run=R.pin(R.__file__),adapter=R.pin(A.__file__),checks=checks,checks_passed=len(checks),native_programs=0,solver_calls=0,scientific_vectors_read=0,mock_invocation_attempts=1,threads=1,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu))
    R.save(out/'RESULT.json',result);print(json.dumps(dict(status=result['status'],checks=len(checks))))
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--out',required=True);main(a.parse_args().out)
