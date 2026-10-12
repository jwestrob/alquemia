#!/usr/bin/env python3
"""Run one authorized login-host Claude client; queue its actual exit to root."""
import datetime
import json
import os
import pathlib
import subprocess
import sys

launch_path = pathlib.Path(sys.argv[1]).resolve()
launch = json.loads(launch_path.read_text())
base = launch_path.parent
if os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Claude clients belong on the login host')
env = os.environ.copy()
env['CLAUDE_CONFIG_DIR'] = launch['claude_config_dir']
with (base / 'launcher.log').open('ab') as log:
    child = subprocess.Popen(launch['command'], cwd=launch['cwd'], env=env,
                             stdin=subprocess.DEVNULL, stdout=log,
                             stderr=subprocess.STDOUT)
    (base / 'client_pid.json').write_text(json.dumps({'pid': child.pid}) + '\n')
    code = child.wait()
receipt = {
    'task': launch['task'], 'returncode': code,
    'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'task_output': launch['task_output'], 'session_uuid': launch['session_uuid'],
    'status': 'enqueue_pending', 'task_accepted': False,
}
path = base / 'completion.json'
path.write_text(json.dumps(receipt, indent=2) + '\n')
message = (
    '[AUTOMATED CLAUDE TASK COMPLETION; not a human message]\n'
    f"Task: {launch['task']}; process exit {code}.\nReceipt: {path}\n"
    f"Read actual task outputs in {launch['task_output']}. "
    'Process exit is not scientific acceptance. All new scientific computation '
    'uses Slurm. Continue the already-authorized staged work, preserving '
    'production, completed results, and current ownership.'
)
try:
    queued = subprocess.run(
        [launch['codex'], 'queue', '--thread', launch['thread'], '--message', message],
        capture_output=True, text=True, timeout=45)
    receipt.update(status='queued' if queued.returncode == 0 else 'enqueue_failed',
                   queue_returncode=queued.returncode,
                   queue_stdout=queued.stdout, queue_stderr=queued.stderr)
except Exception as exc:
    receipt.update(status='enqueue_uncertain', queue_error=repr(exc))
path.write_text(json.dumps(receipt, indent=2) + '\n')
