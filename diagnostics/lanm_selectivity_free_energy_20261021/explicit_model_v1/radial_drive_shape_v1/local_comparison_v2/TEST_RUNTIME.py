#!/usr/bin/env python3
"""Manufactured public-evaluate boundary and terminal persistence controls; no science."""
import argparse,copy,concurrent.futures,json,os,socket,subprocess,threading,time,types
from pathlib import Path
from unittest.mock import patch
import RUN as R
import ADAPTER as A

def main(out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    cfg0=json.loads((R.HERE.parent/'local_comparison_v1/CONFIG_TEMPLATE.json').read_text())
    old=R.module(cfg0['original_adapter'],'runtime_test_old');real=R.module(cfg0['parent_evaluator'],'runtime_test_base')
    checks=[]
    def check(name,ok):
        R.require(ok,name);checks.append(name)
    def forbid(*a,**kw):raise AssertionError('Actual subprocess/solver forbidden')
    class ControlledStop(RuntimeError):pass
    mock_calls=[]
    def invoke(*a):
        mock_calls.append(a[0]);raise ControlledStop('Controlled stop before subprocess')
    eff=lambda v,m,d:list(v[:6] if m=='La' else v[6:12])
    def render(v,row,domain,resolver,exports=None):
        x=eff(v,row['metal'],domain);return 'MOCK '+json.dumps(x)+'\n',{'physical_values_a_b_c_d_R_epsilon':x}
    files=out/'manufactured';files.mkdir();(files/'x').write_text('manufactured coordinates\n');(files/'p').write_text('manufactured parameters\n');(files/'o').write_text('manufactured historical output\n')
    p=R.pin(files/'p');xyz=R.pin(files/'x');output=R.pin(files/'o')
    rows={f'{m}_{i:02d}':{'metal':m,'needs_gradient':i<5,'needs_export':False,'fixed_POLPAIR':0.,'source_input':{'xyz':xyz,'parameters':p}} for m in ('La','Dy') for i in range(23)}
    parent={k:p for k in real.HELPERS};parent.update(implementation=p,source_manifest=p,executables={'analyze':p,'testgrad':p})
    center={'parameters':[.2]*12,'cells':{}}
    for cid,row in rows.items():
        key,meta=render(center['parameters'],row,{},None);kp=files/(cid+'.key');kp.write_text(key)
        q={'cell_id':cid,'effective_values':meta['physical_values_a_b_c_d_R_epsilon'],'xyz':xyz,'parameters':p,'key_sha256':R.hashlib.sha256(key.encode()).hexdigest(),'implementation':p,'source_manifest':p,'executables':parent['executables'],'helpers':{k:p for k in real.HELPERS}}
        center['cells'][cid]={'status':'complete','cache_identity':q,'cache_key':R.digest(q),'request_provenance':{'render':meta},'program_inputs':{'analyze':{'xyz':xyz,'parameters':p,'key':R.pin(kp)}},'receipts':[{'output':output,'returncode':0}]}
    class FakeBase:
        def __init__(self,cfg,path,root_release=None):
            self.config=cfg;self.root_release=root_release;self.cache_dir=Path(path);self.cache_dir.mkdir();self.ledger_path=self.cache_dir/'LEDGER.json'
            self.m={'rows':rows,'units':{}};self.domain={};self.family={};self.initial_vector=[0.]*12
            self.lock=threading.RLock();self.evaluation_lock=threading.RLock();self.slots=threading.Semaphore(1);self.stop_event=threading.Event();self.cell_locks={};self.cache={};self._gated=False
            self.h={'resolver':None,'native_helper':types.SimpleNamespace(invoke=invoke)}
            self.ledger={'native_program_invocations':0,'vector_keys':[],'requests':[],'stopped':False}
        def render(self,v,cid,exports=None):return render(v,rows[cid],{},None,exports)
        cache_identity=real.Evaluator.cache_identity
        def snapshot(self):return {'manufactured':True,'native_program_invocations':self.ledger['native_program_invocations']}
    base=types.SimpleNamespace(**vars(real));base.Evaluator=FakeBase;base.render=render;base.effective=eff;base.checked=lambda v,d:list(v)
    probes=[]
    for j in (12,13):
        for sign in (-1,1):
            v=center['parameters']+[0.,0.];v[j]=sign*.01;probes.append({'key':R.digest(v),'parameters':v,'ordinal':len(probes)+1})
    cfg={'fingerprint':'manufactured','adapter':R.pin(A.__file__),'implementation':R.pin(R.__file__),'parent_evaluator':cfg0['parent_evaluator'],'run_id':'manufactured_boundary','native_source_manifest':p,'executables':parent['executables'],'build':p,'limits':copy.deepcopy(R.LIMITS),'center':p,'threads':list(R.THREADS),'hostname':socket.gethostname(),'parent_configuration':p,'reuse_inventory':p}
    gate=[];E=A.evaluator_class(old,base,types.SimpleNamespace(),lambda c:gate.append(True))
    # Public wrapper/constructor and real journal are retained; only native helper/scientific data are manufactured.
    import clarabel
    original_pool=concurrent.futures.ThreadPoolExecutor
    class SingleWorker(original_pool):
        def __init__(self,*a,**kw):super().__init__(max_workers=1)
    with patch.object(subprocess,'run',forbid),patch.object(subprocess,'Popen',forbid),patch.object(clarabel,'DefaultSolver',forbid),patch.object(old.concurrent.futures,'ThreadPoolExecutor',SingleWorker):
        ev=E(cfg,out/'cache','run',{'new_probes':probes},center,parent)
        first=json.loads((ev.cache_dir/'EVENTS.jsonl').read_text().splitlines()[0])
        check('limits_present_in_committed_initialize',first['event']['data']['ledger']['limits']==R.LIMITS)
        check('limits_are_detached_from_configuration',ev.ledger['limits'] is not cfg['limits'])
        try:ev.evaluate(probes[0]['parameters'])
        except ControlledStop:pass
        else:raise AssertionError('Controlled invocation stop required')
        events=[json.loads(x)['event'] for x in (ev.cache_dir/'EVENTS.jsonl').read_text().splitlines()]
        check('public_request_start_committed',sum(x['kind']=='request_start' for x in events)==1)
        check('one_durable_attempt_before_mock_stop',sum(x['kind']=='attempt_intent' for x in events)==1 and len(mock_calls)==1)
        check('one_failed_attempt_no_returncode',ev.ledger['native_program_invocations']==1 and ev.ledger['native_returncodes_observed']==0 and ev.program_attempts[0]['status']=='failed_native_attempt')
        check('failed_request_keeps_all46_cells',len(ev.ledger['requests'])==1 and len(ev.ledger['requests'][0]['cells'])==46 and ev.ledger['requests'][0]['status']=='failed')
        ev.snapshot();check('real_journal_matches_live_after_controlled_stop',ev.persistence.state.ledger==ev.ledger)
        check('literal_limit_counts_unchanged',ev.ledger['limits']==R.LIMITS and gate==[True])
        # A failed dry request append reproduces the original uncommitted/live mismatch.
        ev2=E(cfg,out/'mismatch_cache','run',{'new_probes':probes},center,parent);append=ev2.persistence.append
        def fail_request(kind,data):
            if kind=='request_start':raise KeyError('manufactured_preappend')
            return append(kind,data)
        ev2.persistence.append=fail_request
        try:ev2.evaluate(probes[0]['parameters'])
        except KeyError as exc:primary=repr(exc);ev2._stop(exc)
        else:raise AssertionError('Injected journal rejection required')
        before=(ev2.cache_dir/'EVENTS.jsonl').read_bytes();res={'mode':'run','status':'failed_radial_local_comparison','error':primary,'solver_calls':0}
        R.persist_result(out,res,ev2,terminal=True)
        check('primary_failure_preserved',res['error']==primary and 'Journal/live ledger mismatch' in res['terminal_snapshot_error'])
        check('unavailable_snapshot_is_not_admitted',res['evaluator_snapshot_available'] is False and 'evaluator' not in res and res['uncertified_live_accounting']['journal_reconciled'] is False)
        check('journal_failure_not_repaired',(ev2.cache_dir/'EVENTS.jsonl').read_bytes()==before and len(before.splitlines())==2)
        try:R.persist_result(out,res,ev2)
        except old.J.JournalError:pass
        else:raise AssertionError('Nonterminal journal mismatch must raise')
        check('nonterminal_snapshot_still_fails_closed',True)
        success={'mode':'run','status':'completed_radial_local_comparison','error':None,'evaluator':{'stale':True}}
        R.persist_result(out,success,ev2,terminal=True)
        check('apparent_success_becomes_terminal_failure',success['status']=='failed_radial_local_comparison' and success['error']==success['terminal_snapshot_error'] and 'evaluator' not in success)
        check('no_second_mock_invocation_after_journal_failure',len(mock_calls)==1)
    result={'status':'passed_radial_runtime_recovery_controls','implementation':R.pin(__file__),'run':R.pin(R.__file__),'adapter':R.pin(A.__file__),'original_adapter':cfg0['original_adapter'],'original_evaluator':cfg0['parent_evaluator'],'checks':checks,'checks_passed':len(checks),'mock_invocation_boundary_attempts':1,'native_programs':0,'solver_calls':0,'scientific_vectors_read':0,'test_worker_threads':1,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu}}
    R.save(out/'RESULT.json',result);print(json.dumps(result))
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--out',required=True);main(a.parse_args().out)
