#!/usr/bin/env python3
"""Run the four pinned upstream dynamics cases once in an actual GPU allocation."""
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import resource
import socket
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ASSERTIONS = {'NVE-Verlet-ArBox': 40, 'NPT-Monte-Iso': 40,
              'NPT-Bussi-Iso': 40, 'NPT-VRESPA': 7}

def pin(path):
    path = Path(path).resolve(strict=True)
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(record):
    if pin(record['path']) != record:
        raise ValueError('Changed artifact: ' + record['path'])
    return Path(record['path'])


def parse(xml, expected):
    doc = ET.fromstring(xml)
    cases = doc.findall('.//TestCase')
    if doc.tag != 'Catch2TestRun' or doc.get('filters') != ','.join('"' + name + '"' for name in expected):
        raise ValueError('Wrong native suite/filter')
    if len(cases) != len(expected) or {c.get('name') for c in cases} != set(expected):
        raise ValueError('Incomplete or unexpected native cases')
    rows = []
    for case in cases:
        overall = case.findall('OverallResult')
        if len(overall) != 1 or overall[0].get('success') != 'true':
            raise ValueError('Native case failed: ' + str(case.get('name')))
        expressions = case.findall('.//Expression')
        if len(expressions) != ASSERTIONS[case.get('name')] or any(e.get('success') != 'true' for e in expressions):
            raise ValueError('Incomplete/failed native case assertions: ' + str(case.get('name')))
        rows.append({'case': case.get('name'), 'passed': True,
                     'successful_assertions': len(expressions),
                     'scope': 'initialization_only' if case.get('name') == 'NPT-VRESPA' else 'trajectory_regression'})
    totals = doc.findall('OverallResults')
    counts = doc.findall('OverallResultsCases')
    if len(totals) != 1 or len(counts) != 1:
        raise ValueError('Missing native completion totals')
    for item in totals + counts:
        if int(item.get('failures', '-1')) != 0 or int(item.get('expectedFailures', '-1')) != 0:
            raise ValueError('Native assertion failure')
    assertions = int(totals[0].get('successes', '0'))
    if assertions <= 0 or int(counts[0].get('successes', '0')) != len(expected):
        raise ValueError('Empty or incomplete native execution')
    if doc.findall('.//Exception') or doc.findall('.//FatalErrorCondition'):
        raise ValueError('Native exception')
    return {'cases': rows, 'successful_assertions': assertions, 'failed_assertions': 0}


def main():
    config_path = Path(sys.argv[1]).resolve(strict=True)
    config = json.loads(config_path.read_text())
    if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('CUDA_VISIBLE_DEVICES'):
        raise RuntimeError('Actual Slurm GPU allocation required')
    for name in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        if os.environ.get(name) != '1':
            raise RuntimeError('Unbounded native threads: ' + name)
    for record in [config['specification'], config['binary'], config['runtime'], config['runner']] + config['sources']:
        verify(record)
    runtime = json.loads(verify(config['runtime']).read_text())
    for record in runtime['files']:
        verify(record['staged'])
    out = Path(config['output'])
    out.mkdir(parents=True, exist_ok=False)
    xml_path = out / 'native.xml'
    console = out / 'native.log'
    argv = [config['binary']['path'], ','.join(config['test_cases']), '--success', '--reporter', 'xml', '--out', str(xml_path)]
    result = {'status': 'running', 'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
              'cuda_visible_devices': os.environ['CUDA_VISIBLE_DEVICES'],
              'configuration': pin(config_path), 'argv': argv, 'cwd': str(out),
              'physical_model_qualification': False, 'native_invocations': 0}
    (out / 'START.json').write_text(json.dumps(result, indent=2) + '\n')
    start = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        result['native_invocations'] = 1
        with console.open('w') as handle:
            proc = subprocess.run(argv, cwd=out, stdin=subprocess.DEVNULL, stdout=handle, stderr=subprocess.STDOUT)
        result['returncode'] = proc.returncode
        if proc.returncode != 0:
            raise RuntimeError('Native process exit ' + str(proc.returncode))
        result['tests'] = parse(xml_path.read_text(), config['test_cases'])
        result['status'] = 'passed'
    except Exception as exc:
        result['status'] = 'failed'
        result['error'] = repr(exc)
    finally:
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        result['elapsed_seconds'] = time.monotonic() - start
        result['native_user_cpu_seconds'] = after.ru_utime - before.ru_utime
        result['native_system_cpu_seconds'] = after.ru_stime - before.ru_stime
        result['finished_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
        result['outputs'] = [pin(p) for p in (xml_path, console) if p.exists()]
        (out / 'RESULT.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return 0 if result['status'] == 'passed' else 1


if __name__ == '__main__':
    sys.exit(main())
