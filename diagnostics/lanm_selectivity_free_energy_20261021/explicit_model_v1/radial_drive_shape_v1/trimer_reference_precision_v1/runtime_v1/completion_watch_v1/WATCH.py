#!/usr/bin/env python3
"""Observe one login phase; queue at most one completion wake, never compute."""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import time

THREAD = '01a0a63a-ee36-7483-b726-1f7f5b7f75cc'
CODEX = Path('/home/jwestrob/.local/bin/codex')
RUNTIME = Path(__file__).resolve().parents[1]/'RUN.py'
TERMINAL = {'reference_precision_all14_complete', 'reference_precision_phase_failed'}


def require(ok, text):
    if not ok:
        raise ValueError(text)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(p):
    require(pin(p['path']) == p, 'Changed pinned file: '+p['path'])
    return json.loads(Path(p['path']).read_text())


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, value):
    p = Path(path); temporary = p.with_name('.'+p.name+'.tmp')
    with temporary.open('w') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(temporary, p)


def process(pid, proc=Path('/proc')):
    p = proc/str(pid)
    try:
        fields = (p/'stat').read_text().rsplit(')', 1)[1].split()
        return {'PID': pid, 'start_ticks': int(fields[19]), 'state': fields[0],
                'argv': (p/'cmdline').read_bytes().decode().strip('\0').split('\0')}
    except (FileNotFoundError, ProcessLookupError):
        return None


def context(c):
    require(not os.environ.get('SLURM_JOB_ID'), 'Login watcher must not claim a Slurm allocation')
    require(c['status'] == 'bound_reference_precision_login_completion_watch' and
            c['implementation'] == pin(__file__) and c['runtime_source'] == pin(RUNTIME) and
            c['thread'] == THREAD and c['poll_seconds'] == 60 and c['codex'] == pin(CODEX), 'Watcher/source contract differs')
    rc, rel = read(c['runtime_configuration']), read(c['runtime_release'])
    launch = read(c['launch_record'])
    require(launch['configuration'] == c['runtime_configuration'] and launch['release'] == c['runtime_release'] and
            launch['controller_PID'] == c['controller']['PID'] and int(launch['process_start_ticks']) == c['controller']['start_ticks'] and
            launch['host'] == c['controller']['host'] and launch['boot_id'] == c['controller']['boot_id'] and
            launch['result_path'] == c['result_path'] and launch['scheduler_job_id'] is None, 'Saved launch identity differs')
    require(rc['status'] == 'bound_reference_precision_login_configuration' and rc['phase'] == 'remaining' and
            rc['implementation'] == c['runtime_source'] and
            rel['status'] == 'approved_reference_precision_login_phase' and rel['phase'] == 'remaining' and
            rel['configuration'] == c['runtime_configuration'] and rel['implementation'] == c['runtime_source'] and
            rel['manifest'] == rc['manifest'] and rel['maximum_phase_attempts'] == 12 and
            rel['scheduler_job_id'] is None, 'Unchanged released remaining-phase binding required')
    require(c['result_path'] == str(Path(rc['output_root'])/'RESULT.json') and
            c['controller']['host'] == socket.gethostname() and
            c['controller']['boot_id'] == Path('/proc/sys/kernel/random/boot_id').read_text().strip(), 'Result/host/boot differs')
    return rc


def observe(c, rc, proc=Path('/proc')):
    snapshot = process(c['controller']['PID'], proc)
    same = snapshot is not None and snapshot['start_ticks'] == c['controller']['start_ticks'] and snapshot['state'] not in ('Z', 'X')
    p = Path(c['result_path']); result = json.loads(p.read_text()) if p.exists() else None
    if result is not None:
        require(result['implementation'] == c['runtime_source'] and
                result['configuration'] == c['runtime_configuration'] and result['release'] == c['runtime_release'] and
                result['manifest'] == rc['manifest'] and result['phase'] == 'remaining', 'Observed result identity differs')
    status = result['status'] if result else 'result_not_yet_written'
    ended = status in TERMINAL or not same
    return {'checked_utc': utc(), 'controller_observation': snapshot, 'same_live_process': same,
            'phase_status': status, 'terminal': ended,
            'reason': 'saved_terminal_result' if status in TERMINAL else
                      ('controller_not_live_without_terminal_result' if not same else 'running'),
            'result': pin(p) if ended and p.exists() else None,
            'scientific_admission_performed': False, 'scheduler_job_id': None}


def deliver(c, folder, observation, runner=subprocess.run):
    intent = folder/'QUEUE_INTENT.json'
    if intent.exists():
        return {'status': 'prior_queue_intent_present_no_duplicate_attempt'}
    message = ('RP REFERENCE PRECISION LOGIN PHASE ENDED: '+observation['phase_status']+
               '. Reason: '+observation['reason']+'. Read '+c['result_path']+' and '+str(folder/'OBSERVATION.json')+
               '. This is a completion wake for the released remaining12 only; it does not admit chemistry or release new calls.')
    argv = [c['codex']['path'], 'queue', '--thread', c['thread'], '--message', message]
    # At-most-once attempt: crash after intent never causes an automatic requeue.
    with intent.open('x') as f:
        json.dump({'utc': utc(), 'argv': argv, 'observation': observation}, f, indent=2)
        f.write('\n'); f.flush(); os.fsync(f.fileno())
    directory_fd = os.open(folder, os.O_RDONLY|os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    tick = time.monotonic()
    try:
        p = runner(argv, capture_output=True, text=True, check=False, timeout=60)
        result = {'status': 'completion_wake_queued' if p.returncode == 0 else 'queue_attempt_failed',
                  'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}
    except Exception as exc:
        result = {'status': 'queue_attempt_failed', 'error': repr(exc)}
    result.update(utc=utc(), argv=argv, wall_seconds=time.monotonic()-tick, automatic_retry=False)
    save(folder/'DELIVERY.json', result)
    return result


def watch(config):
    cp = pin(config); c = read(cp); rc = context(c)
    folder = Path(c['receipts']); folder.mkdir(parents=True, exist_ok=True)
    with (folder/'WATCH.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        require(not (folder/'WATCH.json').exists() or read(pin(folder/'WATCH.json'))['configuration'] == cp, 'Receipt directory belongs to a different watcher')
        save(folder/'WATCH.json', {'configuration': cp, 'PID': os.getpid(), 'host': socket.gethostname(),
                                 'started_utc': utc(), 'scheduler_job_id': None})
        if (folder/'QUEUE_INTENT.json').exists():
            return {'status': 'prior_queue_intent_present_no_duplicate_attempt'}
        while True:
            try:
                observation = observe(c, rc)
            except Exception as exc:
                observation = {'checked_utc': utc(), 'phase_status': 'observation_failed', 'reason': repr(exc),
                               'terminal': True, 'scientific_admission_performed': False, 'scheduler_job_id': None}
            save(folder/'OBSERVATION.json', observation)
            if observation['terminal']:
                return deliver(c, folder, observation)
            time.sleep(c['poll_seconds'])


def bind(launch_path, output, receipts):
    lp = pin(launch_path); launch = read(lp); rc = read(launch['configuration'])
    c = {'status': 'bound_reference_precision_login_completion_watch', 'implementation': pin(__file__),
         'runtime_source': pin(RUNTIME), 'launch_record': lp,
         'runtime_configuration': launch['configuration'], 'runtime_release': launch['release'],
         'result_path': str(Path(rc['output_root'])/'RESULT.json'), 'thread': THREAD, 'codex': pin(CODEX), 'poll_seconds': 60,
         'controller': {'PID': launch['controller_PID'], 'start_ticks': int(launch['process_start_ticks']), 'host': launch['host'],
                        'boot_id': launch['boot_id']},
         'receipts': str(Path(receipts).resolve())}
    context(c)
    require(not Path(receipts).exists(), 'Fresh watcher receipt directory required at binding')
    with Path(output).open('x') as f:
        json.dump(c, f, indent=2); f.write('\n')
    return {'status': 'watcher_configuration_bound_not_launched', 'configuration': pin(output)}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='operation', required=True)
    b = sub.add_parser('bind'); b.add_argument('--launch', required=True)
    b.add_argument('--output', required=True); b.add_argument('--receipts', required=True)
    w = sub.add_parser('watch'); w.add_argument('--config', required=True)
    args = ap.parse_args()
    result = (watch(args.config) if args.operation == 'watch' else
              bind(args.launch, args.output, args.receipts))
    print(json.dumps(result), flush=True)
    raise SystemExit(1 if result['status'] == 'queue_attempt_failed' else 0)
