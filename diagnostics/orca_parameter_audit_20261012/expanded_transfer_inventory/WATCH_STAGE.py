"""Completion wake for one owned corrected-reference stage; never launches science."""
import argparse
import json
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--worker-pid', required=True, type=int)
    parser.add_argument('--worker-start', required=True)
    parser.add_argument('--directory', required=True)
    parser.add_argument('--stage', required=True, choices=['canonical', 'transfer'])
    parser.add_argument('--thread', required=True)
    parser.add_argument('--codex', required=True)
    args = parser.parse_args()
    base = Path(args.directory).resolve()
    receipt = base / (args.stage + '_completion_event.json')
    if receipt.exists():
        raise RuntimeError('Completion receipt exists; inspect before rearming.')
    result = base / ('CANONICAL_COLLECTION.json' if args.stage == 'canonical' else 'FINAL_COLLECTION.json')
    while True:
        try:
            data = json.loads(result.read_text())
            state = 'collection_available'
            break
        except (FileNotFoundError, json.JSONDecodeError):
            pass
        try:
            stat = Path(f'/proc/{args.worker_pid}/stat').read_text().rsplit(')', 1)[1].split()
            alive = stat[0] != 'Z' and stat[19] == args.worker_start
        except FileNotFoundError:
            alive = False
        if not alive:
            state = 'worker_ended_without_collection'
            data = None
            break
        time.sleep(5)
    event = {
        'event': 'corrected-expanded-' + args.stage + '-v1',
        'state': state, 'result': str(result), 'worker_pid': args.worker_pid,
        'worker_start_ticks': args.worker_start, 'created_epoch': time.time(),
        'thread': args.thread, 'scientific_jobs_modified': False,
        'queue_attempts': [],
    }
    if data:
        for key in ('calibration', 'completed_inputs', 'declared_inputs', 'available_pools', 'declared_pools'):
            event[key] = data.get(key)
    message = (
        '[AUTOMATED OWNED COMPLETION EVENT; not a human message]\n'
        f'Event: {event["event"]}\nReceipt: {receipt}\nState: {state}\nResult: {result}\n'
        'Read actual corrected reference results and failures. Continue the already-authorized '
        'staged evaluation only if its frozen canonical gate passes. Preserve production, '
        'historical failures, and the aggregate 16-thread login limit.'
    )

    def save():
        temporary = receipt.with_suffix('.tmp')
        temporary.write_text(json.dumps(event, indent=2) + '\n')
        temporary.replace(receipt)

    save()
    for _ in range(3):
        try:
            run = subprocess.run(
                [args.codex, 'queue', '--thread', args.thread, '--message', message],
                text=True, capture_output=True, timeout=30,
            )
            attempt = {'returncode': run.returncode, 'stdout': run.stdout, 'stderr': run.stderr}
        except Exception as exc:
            attempt = {'returncode': None, 'error': str(exc)}
        event['queue_attempts'].append(attempt)
        save()
        if attempt['returncode'] == 0:
            break
        time.sleep(5)
    print(json.dumps(event), flush=True)


if __name__ == '__main__':
    main()
