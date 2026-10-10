#!/usr/bin/env python3
"""Login-host orchestration only: one tiny Slurm probe, then one Codex wake."""
import argparse
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import time

TERMINAL = {'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL',
            'BOOT_FAIL','DEADLINE','PREEMPTED','REVOKED'}


def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write(path, value):
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2)+'\n'); tmp.replace(path)


def run(cmd):
    # Timeout is for login-host scheduler/queue RPCs, never a scientific process.
    try:
        # Do not inherit an unrelated allocation's sbatch overrides.
        env = {k:v for k,v in os.environ.items() if not k.startswith('SBATCH_')}
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=45, env=env)
        return {'returncode':p.returncode, 'stdout':p.stdout, 'stderr':p.stderr}
    except (OSError, subprocess.TimeoutExpired) as e:
        return {'returncode':None, 'stdout':'', 'stderr':str(e)}


def record_for(job, query):
    if query['returncode'] != 0:
        return None
    for line in query['stdout'].splitlines():
        fields = line.split('|')
        if len(fields) >= 6 and fields[0] == job:
            return dict(zip(('job_id','state','exit_code','elapsed_seconds','allocated_cpus','nodes'), fields[:6]))
    return None


def successful(record, marker_present):
    return bool(record and record['state'].split()[0].rstrip('+') == 'COMPLETED'
                and record['exit_code'] == '0:0' and marker_present)


def submission_outcome(query):
    match = re.fullmatch(r'([0-9]+)(?:;[^\s]+)?',query['stdout'].strip())
    if match:
        return 'watching', match[1]
    if query['returncode'] is not None and query['returncode'] != 0 and 'Batch job submission failed' in query['stderr']:
        return 'waiting', None
    return 'submission_uncertain', None


def notify(state, path, message):
    # Persist intent BEFORE queueing. An uncertain result is not retried blindly.
    state.update(status='enqueue_pending', wake_message=message)
    write(path,state)
    result = run([state['codex'],'queue','--thread',state['thread'],'--message',message])
    state['notification'] = result; state['notification_time'] = stamp()
    state['status'] = 'queued' if result['returncode'] == 0 else 'enqueue_uncertain'
    write(path,state)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--thread',required=True);p.add_argument('--codex',required=True)
    p.add_argument('--interval',type=int,default=1800)
    p.add_argument('--instructions',type=Path,required=True)
    a = p.parse_args()
    if a.interval < 1800:p.error('Do not submit probes more often than every 30 minutes')
    if not a.instructions.is_file():p.error('Readable on-ready instructions required')
    if not os.access(a.codex,os.X_OK):p.error('Executable Codex path required')
    if not re.fullmatch(r'[A-Za-z0-9_-]+',a.directory.name):p.error('Use a simple directory name')
    directory=a.directory.resolve();directory.mkdir(parents=True,exist_ok=True)
    path=directory/'state.json'
    with (directory/'monitor.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        identity={'thread':a.thread,'codex':a.codex,'interval_seconds':a.interval,
                  'instructions':str(a.instructions.resolve())}
        state=json.loads(path.read_text()) if path.exists() else dict(identity,attempts=0,status='waiting',next_attempt_epoch=0)
        if any(state.get(k)!=v for k,v in identity.items()):raise SystemExit('Existing monitor identity differs')
        if state['status'] in ('queued','enqueue_pending','enqueue_uncertain'):
            raise SystemExit('Notification already sent or uncertain; inspect before rearming')
        state.update(pid=os.getpid(),monitor_started_utc=stamp());write(path,state)
        while True:
            if state['status'] in ('submitting','submission_uncertain'):
                # Includes a crash between scheduler acceptance and storing its ID.
                state['status']='submission_uncertain';write(path,state)
                notify(state,path,'[AUTOMATED SLURM READINESS MONITOR ALERT]\n'
                       f'Submission outcome uncertain; inspect {path} and reconcile the unique probe job name before any retry. No readiness or scientific launch is claimed.')
                return
            if state.get('active_job'):
                job=state['active_job']
                q=run(['sacct','-n','-P','-j',job,'--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList'])
                rec=record_for(job,q); state.update(last_checked_utc=stamp(),last_query=q,last_job_record=rec)
                if rec and rec['state'].split()[0].rstrip('+') in TERMINAL:
                    output=Path(state['output_pattern'].replace('%j',job))
                    try:
                        marker_present=state['marker'] in output.read_text().splitlines()
                    except OSError as exc:
                        state['last_output_error']=str(exc)
                        marker_present=False
                    write(Path(state['attempt_directory'])/'terminal.json',dict(record=rec,marker_present=marker_present,observed_utc=stamp()))
                    if successful(rec,marker_present):
                        state['ready_utc']=stamp()
                        notify(state,path,'[AUTOMATED SLURM READINESS EVENT; not a human message]\n'
                               f'Tiny test-partition probe {job} COMPLETED with exit0:0 and its unique output marker.\nReceipt: {path}\n'
                               f'Read and follow {a.instructions.resolve()}. Jacob authorized automatic resumption of owned LanM work once scheduling works. '
                               'Inspect current ownership/results; start or resume the prepared Claude tasks and submit scientifically ready owned jobs without another permission request. '
                               'Claude stays on the login host; all scientific computation goes through Slurm. Do not repeat completed chemistry or touch PQQ production.')
                        return
                    # A completed job with missing/late output is rechecked, not duplicated.
                    if rec['state'].split()[0].rstrip('+')=='COMPLETED' and rec['exit_code']=='0:0':
                        state['status']='awaiting_success_output';write(path,state);time.sleep(30);continue
                    state.update(active_job=None,status='waiting')
                write(path,state)
                if state.get('active_job'):time.sleep(30);continue
            delay=state['next_attempt_epoch']-time.time()
            if delay>0:time.sleep(min(delay,60));continue
            state['attempts']+=1
            attempt=directory/f"attempt_{state['attempts']:05d}";attempt.mkdir(exist_ok=True)
            token=f"NIKASHA_SLURM_READY_{directory.name}_{state['attempts']:05d}"
            name=f"nikready-{directory.name[-10:]}-{state['attempts']}"
            script=attempt/'probe.sbatch'
            script.write_text('#!/bin/bash\n#SBATCH --partition=test\n'
                              f'#SBATCH --job-name={name}\nset -euo pipefail\nhostname\ndate -Is\n'
                              f"printf '%s\\n' '{token}'\n")
            state.update(status='submitting',attempt_directory=str(attempt),job_name=name,marker=token,
                         output_pattern=str(attempt/'probe_%j.out'),next_attempt_epoch=time.time()+a.interval,
                         last_submission_utc=stamp())
            write(path,state)
            q=run(['sbatch','--parsable','--partition=test','--nodes=1','--ntasks=1','--output',state['output_pattern'],'--error',str(attempt/'probe_%j.err'),str(script)])
            write(attempt/'submission.json',q)
            status,job=submission_outcome(q)
            state.update(active_job=job,status=status,last_submission_error=q['stderr'])
            write(path,state)


if __name__=='__main__':main()
