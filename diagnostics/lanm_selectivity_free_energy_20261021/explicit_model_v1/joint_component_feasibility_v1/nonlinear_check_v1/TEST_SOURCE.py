"""Focused source/guard/concurrency mocks; native helper is never called."""
import copy
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import threading
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('fixed_half_adapter',HERE/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def main(out):
    tick,cpu=time.monotonic(),time.process_time();checks={};cfg=json.loads((HERE/'CONFIG_PREPARE.json').read_text())
    base,parent,cert,origin=m.context(cfg);cls=m.evaluator_class(base)
    for name in ('render','cache_identity','_parse','_assemble','_saved_inputs'):
        checks['unchanged_'+name]=getattr(cls,name) is getattr(base.Evaluator,name)
    def rejects(fn,text):
        try:fn()
        except Exception as exc:return text in str(exc)
        return False
    with tempfile.TemporaryDirectory() as temp:
        p=Path(temp);ev=cls(cfg,p/'fresh','prepare')
        checks['exact_allowlisted_vector']=ev.allowed(cfg['candidate']['parameters'])==cert['physical_half_parameters']
        changed=list(cfg['candidate']['parameters']);changed[0]+=1e-5
        checks['unlisted_vector_rejected']=rejects(lambda:ev.allowed(changed),'Only exact')
        changed[0]=100
        checks['out_of_box_rejected']=rejects(lambda:ev.allowed(changed),'Parameter outside')
        checks['prepare_native_gate_closed']=rejects(ev._gate,'Preparation forbids')
        checks['old_cache_rejected']=rejects(lambda:cls(cfg,p/'fresh','prepare'),'Fresh batch cache')
        wrong=copy.deepcopy(cfg);wrong['executor']['hostname']='wrong-host'
        checks['wrong_launcher_rejected']=rejects(lambda:m.context(wrong),'truthful login')
        previous=os.environ.get('SLURM_JOB_ID');os.environ['SLURM_JOB_ID']='123'
        try:checks['fake_slurm_rejected']=rejects(lambda:m.context(cfg),'cannot impersonate')
        finally:
            if previous is None:os.environ.pop('SLURM_JOB_ID')
            else:os.environ['SLURM_JOB_ID']=previous
        altered=copy.deepcopy(parent);qual=m.read(parent['qualification_review']);qual['force_checks_passed']=False
        qp=p/'bad_qualification.json';qp.write_text(json.dumps(qual));altered['qualification_review']=m.pin(qp)
        ev.config=altered;ev.mode='run';ev.root_release={'path':'never_read','sha256':'none'}
        checks['changed_native_qualification_rejected']=rejects(ev._gate,'All four native')
        ev.config=parent;ev.mode='prepare'
        calls=[]
        def native_mock(argv,folder,output,receipts):
            calls.append(str(folder));receipts.append({'returncode':0,'argv':argv,'cwd':str(folder),'mock_only':True})
        ev.h['native_helper']=SimpleNamespace(invoke=native_mock);ev._gate=lambda:None
        ev.ledger['native_program_invocations']=56
        rec={'cell_id':'x','receipts':[]}
        checks['attempt56_ceiling_before_helper']=rejects(lambda:ev._invoke('analyze',p,rec),'Finite56') and not calls
        ev.ledger['native_program_invocations']=0
        def one(i):
            folder=p/str(i);folder.mkdir();r={'cell_id':str(i),'receipts':[]};ev._invoke('analyze',folder,r);return r
        with ThreadPoolExecutor(max_workers=12) as pool:rows=list(pool.map(one,range(12)))
        checks['concurrent12_attempt_publication']=len(calls)==12 and ev.ledger['native_program_invocations']==12 and ev.ledger['native_returncodes_observed']==12 and all(r['receipts'][0]['executor']['job_id'] is None for r in rows)
        checks['truthful_host_pid']=all(r['receipts'][0]['executor']['hostname']==m.socket.gethostname() and r['receipts'][0]['executor']['controller_pid']==os.getpid() for r in rows)
        def fail(argv,folder,output,receipts):
            receipts.append({'returncode':1,'mock_only':True});raise RuntimeError('synthetic first failure')
        ev.h['native_helper']=SimpleNamespace(invoke=fail)
        checks['first_failure_sets_stop']=rejects(lambda:ev._invoke('analyze',p,{'cell_id':'failure','receipts':[]}),'synthetic first failure') and ev.stop_event.is_set()
        before=ev.ledger['native_program_invocations']
        checks['queued_attempt_preserved']=rejects(lambda:ev._invoke('analyze',p,{'cell_id':'queued','receipts':[]}),'Global stop before') and ev.ledger['native_program_invocations']==before
        checks['queued_cell_never_prepared']=rejects(lambda:ev._cell(cfg['candidate']['parameters'],next(iter(ev.m['rows']))),'Untouched cell')
        ev._stop(RuntimeError('later failure'))
        checks['first_cause_preserved']='synthetic first failure' in ev.ledger['failure']
        checks['shared_attempts_saved_consistently']=json.loads(ev.ledger_path.read_text())['native_attempts']==ev.program_attempts
    m.require(all(checks.values()),str(checks))
    result={'status':'passed_fixed_half_source_controls','implementation':m.pin(HERE/'RUN.py'),'configuration':m.pin(HERE/'CONFIG_PREPARE.json'),'tests':m.pin(__file__),'checks':checks,'count':len(checks),'mock_native_attempts':13,'real_native_calls':0,'solver_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
    with Path(out).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();main(a.out)
