#!/usr/bin/env python3
"""Synthetic process/result and delivery controls; never queue a real wake."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import time
from types import SimpleNamespace
import WATCH as w


def main(output):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    tick,cpu=time.monotonic(),time.process_time();checks=[]
    def check(ok,name):
        w.require(ok,name);checks.append(name)
    def reject(fn,name):
        try:fn()
        except ValueError:checks.append(name);return
        raise ValueError('Accepted corrupt control: '+name)
    def forbidden(*a,**kw):raise AssertionError('Real subprocess forbidden')
    subprocess.Popen=forbidden
    result={'status':'failed_login_completion_watcher_controls','error':None,'source':w.pin(w.__file__),
            'implementation':w.pin(__file__),'native_calls':0,'real_queue_attempts':0}
    try:
        proc=out/'synthetic_proc';p=proc/'77';p.mkdir(parents=True)
        fields=['0']*50;fields[0]='S';fields[19]='1234'
        (p/'stat').write_text('77 (name with ) embedded) '+' '.join(fields))
        (p/'cmdline').write_bytes(b'python\0RUN.py\0')
        c={'controller':{'PID':77,'start_ticks':1234},'result_path':str(out/'synthetic_RESULT.json'),
           'runtime_source':'S','runtime_configuration':'C','runtime_release':'R','thread':w.THREAD,
           'codex':{'path':'/synthetic/codex'}};rc={'manifest':'M'}
        check(w.process(77,proc)['start_ticks']==1234,'PID_start_ticks_with_parentheses')
        check(not w.observe(c,rc,proc)['terminal'],'live_missing_result_waits')
        row={'implementation':'S','configuration':'C','release':'R','manifest':'M','phase':'remaining','status':'started'}
        Path(c['result_path']).write_text(json.dumps(row))
        check(not w.observe(c,rc,proc)['terminal'],'live_started_waits')
        bad=copy.deepcopy(c);bad['controller']['start_ticks']=999
        check(w.observe(bad,rc,proc)['reason']=='controller_not_live_without_terminal_result','PID_reuse_detected')
        for status in w.TERMINAL:
            row['status']=status;Path(c['result_path']).write_text(json.dumps(row))
            check(w.observe(c,rc,proc)['reason']=='saved_terminal_result','saved_terminal_'+status)
        row['configuration']='changed';Path(c['result_path']).write_text(json.dumps(row))
        reject(lambda:w.observe(c,rc,proc),'changed_result_identity')
        row['configuration']='C';row['status']='started';Path(c['result_path']).write_text(json.dumps(row))
        fields[0]='Z';(p/'stat').write_text('77 (name) '+' '.join(fields))
        check(w.observe(c,rc,proc)['terminal'],'zombie_is_not_running')
        (p/'stat').unlink()
        check(w.observe(c,rc,proc)['terminal'],'gone_controller_without_terminal_result')
        delivered=[];folder=out/'mock_delivery';folder.mkdir()
        def success(argv,**kw):
            check((folder/'QUEUE_INTENT.json').exists(),'intent_before_queue')
            delivered.append(argv);return SimpleNamespace(returncode=0,stdout='synthetic queued',stderr='')
        obs={'phase_status':'reference_precision_all14_complete','reason':'saved_terminal_result'}
        check(w.deliver(c,folder,obs,success)['status']=='completion_wake_queued','one_successful_wake')
        check(w.deliver(c,folder,obs,forbidden)['status']=='prior_queue_intent_present_no_duplicate_attempt','no_duplicate_delivery')
        check(len(delivered)==1 and delivered[0][:4]==['/synthetic/codex','queue','--thread',w.THREAD],'exact_thread_single_attempt')
        folder2=out/'mock_failed_delivery';folder2.mkdir()
        def failed(argv,**kw):return SimpleNamespace(returncode=7,stdout='',stderr='synthetic refusal')
        check(w.deliver(c,folder2,obs,failed)['status']=='queue_attempt_failed','queue_failure_preserved')
        check(w.deliver(c,folder2,obs,forbidden)['status']=='prior_queue_intent_present_no_duplicate_attempt','failed_delivery_no_retry')
        result['status']='passed_login_completion_watcher_controls'
    except Exception as exc:result['error']=repr(exc)
    result.update(checks=checks,check_count=len(checks),cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1})
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'error':result['error']}))
    return 0 if result['error'] is None else 1


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);args=ap.parse_args();raise SystemExit(main(args.out))
