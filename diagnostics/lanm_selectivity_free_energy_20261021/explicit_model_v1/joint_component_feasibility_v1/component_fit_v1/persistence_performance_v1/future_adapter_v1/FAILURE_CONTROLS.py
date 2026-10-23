"""Focused future _invoke failure mocks; one thread, zero native calls."""
import importlib.util
import json
from pathlib import Path
import threading
import time

HERE=Path(__file__).resolve().parent;OUT=HERE/'failure_controls_v1'
def load(p,n):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
A=load(HERE/'ADAPTER.py','failure_adapter');J=A.J
class Base:
 class Evaluator:pass
C=A.evaluator_class(Base)

def make(name):
 e=object.__new__(C);d=OUT/name;d.mkdir()
 e.lock=threading.RLock();e.slots=threading.BoundedSemaphore(1);e.stop_event=threading.Event();e._gate=lambda:None
 e.ledger={'vector_keys':['v'],'native_program_invocations':0,'native_returncodes_observed':0,'requests':[{'status':'started','vector_key':'v','cells':{'c':{'status':'not_attempted'}}}],'stopped':False,'limits':{'new_parameter_vectors':2,'new_native_programs':3}}
 initial={**e.ledger,'vector_keys':[],'requests':[]}
 e.program_attempts=[];e.cell_attempts={};e.executor_identity={'mode':'mock'};e.config={'executables':{'analyze':{'path':'not_executed'}}};e.ledger_path=d/'LEDGER.json'
 e.persistence=J.Journal(d/'EVENTS.jsonl',initial,e.lock,{'mock':True})
 e.persistence.append('request_start',{'index':0,'request':e.ledger['requests'][0]})
 return e,d,{'cache_key':'token','cell_id':'c','receipts':[],'program_inputs':{'analyze':{'mock':True}}}

def main():
 OUT.mkdir(exist_ok=False);checks=[];start=time.process_time()
 def ck(n,b):
  checks.append({'name':n,'passed':bool(b)})
  if not b:raise AssertionError(n)
 e,d,r=make('native_failure');calls=[]
 class Failure:
  @staticmethod
  def invoke(argv,folder,stdout,receipts):
   calls.append('mock');receipts.append({'returncode':7,'mock':True});raise RuntimeError('mock native failure')
 e.h={'native_helper':Failure}
 try:e._invoke('analyze',d,r)
 except RuntimeError as ex:ck('original_mock_failure_propagates',str(ex)=='mock native failure')
 else:raise AssertionError('failure swallowed')
 rr=J.replay(e.persistence.path)
 ck('stop_and_return_evidence_preserved',e.stop_event.is_set() and rr['ledger']['stopped'] and rr['ledger']['native_attempts'][0]['status']=='failed_native_attempt' and rr['ledger']['native_returncodes_observed']==1)
 e.persistence.materialize(e.ledger_path,expected=e.ledger)
 ck('failed_live_ledger_replays_exactly',rr['ledger']==e.ledger)
 before=e.persistence.path.read_bytes()
 try:e._invoke('analyze',d,r)
 except A.NotAttempted:pass
 else:raise AssertionError('stopped queued call launched')
 ck('queued_after_stop_untouched',len(calls)==1 and e.persistence.path.read_bytes()==before)
 e.persistence.close()
 e,d,r=make('write_failure');calls=[]
 class Never:
  @staticmethod
  def invoke(*args):calls.append('forbidden')
 e.h={'native_helper':Never};original=J.os.write
 def denied(*args):raise OSError('synthetic persistence failure')
 try:
  J.os.write=denied
  try:e._invoke('analyze',d,r)
  except OSError as ex:ck('intent_write_failure_propagates',str(ex)=='synthetic persistence failure')
  else:raise AssertionError('persistence failure swallowed')
 finally:J.os.write=original
 ck('failed_intent_prevents_helper',calls==[] and e.persistence.failed)
 replay=J.replay(e.persistence.path)
 ck('committed_prefix_kept_on_failed_intent',replay['budget_charged_attempts']==0 and replay['sequence']==2)
 try:e.persistence.append('stop',{'error':'must not append after I/O failure'})
 except J.JournalError:ck('journal_closed_to_further_events_after_IO_failure',True)
 else:raise AssertionError('continued after persistence failure')
 e.persistence.close()
 result={'status':'passed_focused_adapter_failure_controls','implementation':J.pin(Path(__file__)),'adapter':J.pin(HERE/'ADAPTER.py'),'journal':J.pin(HERE/'JOURNAL.py'),
 'checks':checks,'checks_passed':len(checks),'cpu_seconds':time.process_time()-start,'threads':1,'OS_worker_threads_created':0,'native_calls':0,'scientific_evaluations':0,'mock_failure_calls':1}
 with (OUT/'RESULT.json').open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps({'status':result['status'],'checks':len(checks),'cpu_seconds':result['cpu_seconds'],'result':J.pin(OUT/'RESULT.json')}))
if __name__=='__main__':main()
