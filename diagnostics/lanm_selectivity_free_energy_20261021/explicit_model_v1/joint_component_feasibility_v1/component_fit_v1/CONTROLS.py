#!/usr/bin/env python3
"""Small mock execution-policy checks. No native/solver invocation."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import threading
import time
import types
from unittest.mock import patch
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('policy_control',HERE/'ADAPTER.py');A=importlib.util.module_from_spec(s);s.loader.exec_module(A)
checks=[]
def check(name,ok):A.require(ok,name);checks.append(name)
def rejects(name,f):
    try:f()
    except (ValueError,A.BudgetStop,A.NotAttempted):checks.append(name);return
    raise AssertionError(name)
class Base:
    def snapshot(self):return {}
class SerialPool:
    def __init__(self,**kw):self.pool=ORIGINAL_POOL(max_workers=1)
    def __enter__(self):return self.pool.__enter__()
    def __exit__(self,*args):return self.pool.__exit__(*args)
ORIGINAL_POOL=A.concurrent.futures.ThreadPoolExecutor
C=A.evaluator_class(types.SimpleNamespace(Evaluator=Base))
def instance(folder,fail=None):
    ev=object.__new__(C);ev.cache_dir=Path(folder);ev.cache_dir.mkdir();ev.lock=threading.RLock();ev.evaluation_lock=threading.Lock();ev.slots=threading.BoundedSemaphore(12);ev.stop_event=threading.Event()
    ev.ledger={'native_program_invocations':0,'native_returncodes_observed':0,'vector_keys':[],'requests':[],'stopped':False};ev.ledger_path=ev.cache_dir/'LEDGER.json';A.save(ev.ledger_path,ev.ledger)
    ev.program_attempts=[];ev.cell_attempts={};ev.imported_cell_keys=set();ev.started_cells=set();ev.allowed={};ev.complete_vectors={};ev.imported_vector_keys=set();ev.cache={}
    ev.batch={'adapter':A.pin(HERE/'ADAPTER.py')};ev.batch_fingerprint='mock';ev.executor_identity={'mode':'mock','job_id':None};ev.config={'executables':{'analyze':{'path':'MOCK_analyze'},'testgrad':{'path':'MOCK_testgrad'}}}
    ev.m={'rows':{f'c{i:02}':{'needs_gradient':i<10} for i in range(46)}};ev.checked_values=lambda x:list(x);ev._gate=lambda: A.require(not ev.stop_event.is_set(),'Stopped mock gate')
    ev.required_programs=lambda values:56
    ev.cache_identity=lambda values,cid:(cid,{},'',{})
    def invoke(argv,folder,output,receipts):
        if fail=='native':raise ValueError('mock native failure')
        receipts.append({'argv':argv,'returncode':0,'cwd':str(folder)})
    ev.h={'native_helper':types.SimpleNamespace(invoke=invoke)}
    def cell(values,cid):
        if ev.stop_event.is_set():raise A.NotAttempted('mock untouched')
        if fail=='precall':raise ValueError('mock input failure')
        rec={'cell_id':cid,'cache_key':cid,'receipts':[],'status':'started'}
        for prog in ['analyze']+(['testgrad'] if ev.m['rows'][cid]['needs_gradient'] else []):ev._invoke(prog,ev.cache_dir/cid,rec)
        rec['status']='complete';return rec
    ev._cell=cell;ev._assemble=lambda cells:{'cells':cells,'model':{},'response':{},'score':{'objective':0}}
    return ev

def main():
    tick,cpu=time.monotonic(),time.process_time()
    with tempfile.TemporaryDirectory() as td,patch.object(A.concurrent.futures,'ThreadPoolExecutor',SerialPool):
        p=Path(td);v=[0.]*12
        e=instance(p/'success');rejects('unauthorized_vector_rejected',lambda:e.evaluate(v));e.authorize(v,'trial',1);r=e.evaluate(v)
        check('full46_admitted',len(r['cells'])==46);check('full56_actual_boundaries',e.ledger['native_program_invocations']==56);check('full56_returncodes',e.ledger['native_returncodes_observed']==56)
        check('one_new_vector',len(e.ledger['vector_keys'])==1);e.authorize(v,'probe',2);r2=e.evaluate(v);check('same_vector_reuses_without_programs',e.ledger['native_program_invocations']==56 and r==r2)
        for name in ('precall','native'):
            e=instance(p/name,name);e.authorize(v,'trial',1);rejects(name+'_stops_campaign',lambda:e.evaluate(v));check(name+'_stopped',e.stop_event.is_set());check(name+'_correct_actual_attempts',e.ledger['native_program_invocations']==(name=='native'))
            check(name+'_all46_statuses_preserved',len(e.ledger['requests'][0]['cells'])==46)
            check(name+'_45_queued_untouched',sum(x['status']=='not_attempted' for x in e.ledger['requests'][0]['cells'].values())==45)
            rejects(name+'_no_retry',lambda:e.evaluate(v))
        e=instance(p/'cap');e.authorize(v,'trial',1);e.ledger['native_program_invocations']=6495;rejects('whole_vector_native_budget_stop',lambda:e.evaluate(v));check('budget_no_partial_launch',not e.ledger['vector_keys'] and not e.program_attempts)
        e.ledger['native_program_invocations']=0;e.ledger['vector_keys']=[str(i) for i in range(200)];rejects('newvector200_cap',lambda:e.evaluate(v));check('vectorcap_no_calls',not e.program_attempts)
        e.complete_vectors[A.digest(v)]={'status':'complete','parameters':v,'score':{'stored':True}};check('inherited_reuse_at_cap',e.evaluate(v)['score']=={'stored':True})
        e=instance(p/'gate');e.mode='prepare';e._gate=types.MethodType(C._gate,e);rejects('preparation_native_guard',e._gate)
        check('inherited_scientific_methods_untouched',all(name not in C.__dict__ for name in ('render','cache_identity','_parse','_assemble')))
    out=Path(__import__('sys').argv[1]);out.mkdir(parents=True,exist_ok=False)
    result={'status':'passed_component_execution_controls','implementation':A.pin(__file__),'adapter':A.pin(HERE/'ADAPTER.py'),'checks':checks,'passed':len(checks),'native_programs':0,'real_proposal_calls':0,'mock_invocations_are_not_native':True,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    A.save(out/'RESULT.json',result);print(json.dumps({'result':A.pin(out/'RESULT.json'),'passed':len(checks)}))
if __name__=='__main__':main()
