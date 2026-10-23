"""One-thread artificial persistence checks; native helper is a pure mock."""
import ast
import copy
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import sys
import threading
import time

HERE=Path(__file__).resolve().parent
OUT=HERE/'controls_v1'

def load(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
J=load(HERE/'JOURNAL.py','test_journal')
A=load(HERE/'ADAPTER.py','test_future_adapter')

def ledger(cap=4):return {'config_sha256':'synthetic','vector_keys':[],'native_program_invocations':0,'native_returncodes_observed':0,'requests':[],'stopped':False,'limits':{'new_parameter_vectors':2,'new_native_programs':cap}}
def request():return {'values':[0.0],'vector_key':'v1','role':{'role':'trial','outer_step':1},'status':'started','required_new_programs':2,'cells':{'c1':{'status':'not_attempted'},'c2':{'status':'not_attempted'}}}
def attempt(number):return {'number':number,'cache_key':'token'+str(number),'cell_id':'c'+str(number),'program':'analyze','folder':'synthetic'+str(number),'status':'invoke_boundary_entered'}
def end(number,failed=False):return {**attempt(number),'status':'failed_native_attempt' if failed else 'returned_success','returncode_observed':True}
def start_event(j):j.append('request_start',{'index':0,'request':request()})
def begin(j,number):j.append('attempt_intent',{'attempt':attempt(number),'request_index':0})
def finish(j,number,failed=False):j.append('attempt_return',{'attempt':end(number,failed),'receipt_delta':[{'returncode':1 if failed else 0}]})


def main():
 OUT.mkdir(exist_ok=False);checks=[];tick,cpu=time.monotonic(),time.process_time()
 def ck(name,value):
  checks.append({'name':name,'passed':bool(value)})
  if not value:raise AssertionError(name)
 def reject(name,fn):
  try:fn()
  except J.JournalError:ck(name,True)
  else:ck(name,False)
 def journal(name,cap=4):return J.Journal(OUT/(name+'.jsonl'),ledger(cap),threading.RLock(),{'test_only':True})
 j=journal('complete');start_event(j);begin(j,1);begin(j,2);finish(j,2);finish(j,1)
 for n in (1,2):j.append('request_cell',{'index':0,'cell_id':'c'+str(n),'state':{'status':'complete','cache_key':'token'+str(n),'immutable_import':False}})
 j.append('request_update',{'index':0,'updates':{'status':'complete','objective':.1}})
 expected=ledger();expected.update(vector_keys=['v1'],native_program_invocations=2,native_returncodes_observed=2,native_attempts=[end(1),end(2)])
 req=request();req.update(status='complete',objective=.1)
 for n in (1,2):req['cells']['c'+str(n)]={'status':'complete','cache_key':'token'+str(n),'immutable_import':False}
 expected['requests']=[req]
 ck('lossless_incremental_state',j.state.ledger==expected)
 j.materialize(OUT/'LEDGER.json',expected)
 ck('legacy_ledger_bytes_semantics',json.loads((OUT/'LEDGER.json').read_text())==expected)
 replay=J.replay(j.path);ck('lossless_full_replay',replay['ledger']==expected and not replay['unresolved_attempt_numbers'] and replay['sequence']==9)
 before=j.path.read_bytes();reject('duplicate_intent_no_append',lambda:begin(j,1));ck('duplicate_preserves_journal',j.path.read_bytes()==before)
 reject('duplicate_return_no_append',lambda:finish(j,1))
 bad=copy.deepcopy(expected);bad['native_program_invocations']=3
 reject('materialize_detects_live_mismatch',lambda:j.materialize(OUT/'BAD_LEDGER.json',bad))
 j.close()
 # Intent-only crash remains charged. No constructor/replay supports continuation.
 pending=journal('pending');start_event(pending);begin(pending,1);pending.close()
 r=J.replay(pending.path);ck('unresolved_intent_charged_no_retry',r['budget_charged_attempts']==1 and r['unresolved_attempt_numbers']==[1] and r['automatic_retry_permitted'] is False)
 with pending.path.open('ab') as f:f.write(b'{"sequence":4,"event":')
 r=J.replay(pending.path);ck('partial_suffix_retained',r['partial_suffix_bytes']>0 and r['committed_prefix_only'] and r['budget_charged_attempts']==1)
 ck('partial_file_not_rewritten',pending.path.read_bytes().endswith(b'{"sequence":4,"event":'))
 corrupted=OUT/'corrupted.jsonl';corrupted.write_bytes(before.replace(b'"synthetic"',b'"changedxx"',1))
 reject('committed_corruption_rejected',lambda:J.replay(corrupted))
 # Stop permits return evidence for an already active attempt, never a new launch.
 stopped=journal('stopped');start_event(stopped);begin(stopped,1);stopped.append('stop',{'error':'synthetic failure'})
 reject('stop_rejects_new_intent',lambda:begin(stopped,2));finish(stopped,1,True)
 stopped.append('request_cell',{'index':0,'cell_id':'c1','state':{'status':'failed_partial','native_attempts':1,'error':'synthetic'}})
 stopped.append('request_cell',{'index':0,'cell_id':'c2','state':{'status':'not_attempted','native_attempts':0,'error':'global stop'}})
 stopped.append('request_update',{'index':0,'updates':{'status':'failed','error':'synthetic'}})
 stopped.materialize(OUT/'FAILED_LEDGER.json');stopped.close()
 sr=J.replay(stopped.path);ck('failure_and_untouched_cell_preserved',sr['ledger']['stopped'] and sr['ledger']['failure']=='synthetic failure' and sr['ledger']['requests'][0]['cells']['c2']['status']=='not_attempted' and sr['ledger']['native_program_invocations']==1)
 capped=journal('capped',1);start_event(capped);begin(capped,1);reject('hard_attempt_cap',lambda:begin(capped,2));capped.close()
 # Exhaustive valid interleavings: sequence numbers assigned when starts enter
 # the same lock; logical worker B may finish before A. No OS worker is added.
 orders=[]
 for order in itertools.permutations(('b1','e1','b2','e2')):
  if order.index('b1')<order.index('e1') and order.index('b2')<order.index('e2') and order.index('b1')<order.index('b2'):orders.append(order)
 for i,order in enumerate(orders):
  z=journal('interleaved_'+str(i));start_event(z)
  for ev in order:(begin if ev[0]=='b' else finish)(z,int(ev[1]))
  z.close();rr=J.replay(z.path)
  ck('interleaving_'+''.join(order),rr['ledger']['native_program_invocations']==2 and rr['ledger']['native_returncodes_observed']==2 and not rr['unresolved_attempt_numbers'])
 # Actual future _invoke function with pure receipt mock, no subprocess.
 class FakeBase:
  class Evaluator:pass
 C=A.evaluator_class(FakeBase);ev=object.__new__(C)
 ev.lock=threading.RLock();ev.slots=threading.BoundedSemaphore(1);ev.stop_event=threading.Event();ev._gate=lambda:None
 ev.ledger=ledger();ev.program_attempts=[];ev.cell_attempts={};ev.executor_identity={'mode':'mock_only'}
 ev.config={'executables':{'analyze':{'path':'never_executed'}}};ev.ledger_path=OUT/'MOCK_LEDGER.json'
 ev.persistence=J.Journal(OUT/'mock.jsonl',ev.ledger,ev.lock,{'test_only':True})
 req=request();ev.ledger['vector_keys'].append('v1');ev.ledger['requests'].append(req);start_event(ev.persistence)
 called=[]
 class NativeMock:
  @staticmethod
  def invoke(argv,folder,stdout,receipts):
   rr=J.replay(ev.persistence.path)
   assert rr['budget_charged_attempts']==1 and rr['unresolved_attempt_numbers']==[1]
   called.append('durable_intent_precedes_mock_call');receipts.append({'returncode':0,'mock':True})
 ev.h={'native_helper':NativeMock}
 record={'cache_key':'token1','cell_id':'c1','receipts':[],'program_inputs':{'analyze':{'test_only':True}}}
 ev._invoke('analyze',OUT/'never_executed',record)
 ck('adapter_intent_before_helper_and_return_replay',called==['durable_intent_precedes_mock_call'] and ev.ledger==J.replay(ev.persistence.path)['ledger'])
 ck('return_executor_receipt_retained',J.replay(ev.persistence.path)['ledger']['native_returncodes_observed']==1 and record['receipts'][0]['executor']==ev.executor_identity)
 ev.persistence.materialize(ev.ledger_path,ev.ledger);ev.persistence.close()
 # Source isolation and exact scientific/native scheduling structure.
 manifest=json.loads((HERE/'MANIFEST.json').read_text())
 for k in ('original_adapter','original_driver','future_adapter','future_driver','journal'):
  ck('source_pin_'+k,J.pin(manifest[k]['path'])==manifest[k])
 old=Path(manifest['original_adapter']['path']).read_text();new=(HERE/'ADAPTER.py').read_text()
 def defs(s):return {n.name:ast.dump(n,include_attributes=False) for n in ast.walk(ast.parse(s)) if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
 od,nd=defs(old),defs(new)
 for name in ('authorize','_gate','_cell','required_programs','cache_identity' if 'cache_identity' in od else 'cell_digest'):
  ck('unchanged_AST_'+name,od[name]==nd[name])
 ck('driver_byte_identical',Path(manifest['original_driver']['path']).read_bytes()==(HERE/'DRIVER.py').read_bytes())
 def withs(s):
  fn=next(n for n in ast.walk(ast.parse(s)) if isinstance(n,ast.FunctionDef) and n.name=='_invoke')
  return [ast.dump(x.context_expr,include_attributes=False) for n in ast.walk(fn) if isinstance(n,ast.With) for x in n.items]
 ck('native_slots_and_lock_structure_unchanged',withs(old)==withs(new))
 ck('one_original_native_helper_call',new.count("self.h['native_helper'].invoke(")==old.count("self.h['native_helper'].invoke(")==1)
 result={'status':'passed_artificial_persistence_controls','implementation':J.pin(Path(__file__)),'manifest':J.pin(HERE/'MANIFEST.json'),
  'checks':checks,'checks_passed':len(checks),'wall_seconds':time.monotonic()-tick,'cpu_seconds':time.process_time()-cpu,
  'numerical_threads':1,'OS_worker_threads_created':0,'logical_writer_interleavings':len(orders),
  'native_calls':0,'scientific_evaluations':0,'mock_helper_calls':len(called),
  'limits':'Interleavings tested sequentially under one thread; this is not scheduler stress testing or a production crash/fsync guarantee on the shared filesystem.'}
 with (OUT/'RESULT.json').open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps({'status':result['status'],'checks':len(checks),'cpu_seconds':result['cpu_seconds'],'result':J.pin(OUT/'RESULT.json')}))
if __name__=='__main__':main()
