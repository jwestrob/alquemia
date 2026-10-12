#!/usr/bin/env python3
"""Thin login-host client supervisor; records termination and queues root once."""
import datetime,json,pathlib,subprocess,sys
launch_path=pathlib.Path(sys.argv[1]);launch=json.loads(launch_path.read_text())
base=launch_path.parent
with (base/'launcher.log').open('ab') as log:
    p=subprocess.Popen(launch['command'],cwd=launch['cwd'],env=None,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
    (base/'client_pid.json').write_text(json.dumps({'pid':p.pid})+'\n')
    code=p.wait()
receipt={'returncode':code,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'task_output':launch['task_output'],'session_uuid':launch['session_uuid'],'status':'enqueue_pending','science_accepted':False}
path=base/'completion.json';path.write_text(json.dumps(receipt,indent=2)+'\n')
message=('[AUTOMATED CLAUDE TASK COMPLETION; not a human message]\n'
         f"Task: orca-parameter-runtime-20261012; process exit {code}.\nReceipt: {path}\n"
         f"Read actual task outputs in {launch['task_output']}. Process exit is not scientific acceptance. "
         'Continue the user-authorized parameter/corrected-discrimination work under the aggregate16-thread login limit. Preserve production and original failures.')
try:
    p=subprocess.run([launch['codex'],'queue','--thread',launch['thread'],'--message',message],capture_output=True,text=True,timeout=45)
    receipt.update(status='queued' if p.returncode==0 else 'enqueue_failed',queue_returncode=p.returncode,queue_stdout=p.stdout,queue_stderr=p.stderr)
except Exception as exc:
    receipt.update(status='enqueue_uncertain',queue_error=repr(exc))
path.write_text(json.dumps(receipt,indent=2)+'\n')
