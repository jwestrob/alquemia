"""Run only the eight real-artifact collector checks; require their Slurm receipt."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import time
import xml.etree.ElementTree as ET

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
TEST_FILE = Path(__file__).resolve().parents[3] / 'tests/test_metal_environment_lady_embedded.py'
PYTEST_PYTHON = '/home/jwestrob/.pyenv/versions/3.12.0/bin/python'
MODULES = ('metal_environment_lady_embedded.py', 'metal_environment_lady_embedded_prepare.py',
           'metal_environment_frozen_f.py', 'metal_environment_lanm_reference.py',
           'metal_environment_reference.py', 'metal_environment_interrupted_seed.py',
           'metal_environment_orbital_seed.py', 'ggr_sensitivity.py',
           'affordable_common.py', 'affordable_workflow.py', 'run_orca_task_manifest.py')
TESTS = {'test_real_accepted_moread_la_endpoint_collects', 'test_real_accepted_trah_dy_endpoint_collects',
         'test_relocated_moread_la_endpoint_collects', 'test_relocated_trah_dy_endpoint_collects',
         'test_silent_pmodel_fallback_rejected', 'test_missing_moread_file_evidence_rejected',
         'test_missing_trah_evidence_rejected', 'test_point_charge_count_mismatch_rejected'}


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def outcomes(path):
    cases = ET.parse(path).getroot().findall('.//testcase')
    names = [case.attrib['name'] for case in cases]
    return {'tests': len(cases), 'names': names, 'passed': len(cases) == 8 and set(names) == TESTS
            and all(not any(case.find(tag) is not None for tag in ('failure', 'error', 'skipped')) for case in cases)}


def run(directory):
    assert os.environ.get('SLURM_JOB_ID'), 'Collector qualification requires Slurm'
    directory.mkdir(parents=True, exist_ok=False)
    source_pins = {name: pin(ROOT / 'scripts' / name) for name in MODULES}
    test_pin = pin(TEST_FILE)
    command = [PYTEST_PYTHON, '-m', 'pytest', str(TEST_FILE), '-q', '-o', 'addopts=',
               '--junitxml', str(directory / 'pytest.xml'), '--basetemp', str(directory / 'tmp')]
    env = os.environ.copy()
    env.update({key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                    'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')})
    env['PYTEST_DISABLE_PLUGIN_AUTOLOAD'] = '1'
    started = time.time()
    with (directory / 'pytest.log').open('x') as log:
        completed = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    result = outcomes(directory / 'pytest.xml') if (directory / 'pytest.xml').exists() else {'passed': False}
    stable = test_pin == pin(TEST_FILE) and all(value == pin(value['path']) for value in source_pins.values())
    receipt = {'status': 'passed' if completed.returncode == 0 and result['passed'] and stable else 'failed',
               'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
               'allocated_cpus': os.environ.get('SLURM_CPUS_ON_NODE'), 'command': command,
               'started_epoch': started, 'finished_epoch': time.time(), 'returncode': completed.returncode,
               'outcomes': result, 'test_file': test_pin, 'source_pins': source_pins,
               'gate': pin(__file__), 'log': pin(directory / 'pytest.log'), 'new_molecular_calls': 0,
               'junit': pin(directory / 'pytest.xml') if (directory / 'pytest.xml').exists() else None}
    with (directory / 'RECEIPT.json').open('x') as handle:
        json.dump(receipt, handle, indent=2)
    if receipt['status'] != 'passed':
        raise RuntimeError('Collector tests did not qualify; read preserved receipt/log')


def check(receipt_path, manifest_path):
    assert os.environ.get('SLURM_JOB_ID'), 'Molecular test gate requires Slurm'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'passed' and receipt['returncode'] == 0
    assert receipt['test_file'] == pin(TEST_FILE) and receipt['gate'] == pin(__file__)
    for key in ('log', 'junit'):
        assert receipt[key] == pin(receipt[key]['path']), 'Changed test evidence'
    assert receipt['outcomes'] == outcomes(receipt['junit']['path']) and receipt['outcomes']['passed']
    manifest = json.loads(manifest_path.read_text())
    for name in MODULES:
        tested = receipt['source_pins'][name]
        assert tested == pin(tested['path'])
        assert manifest['implementation'][name]['sha256'] == tested['sha256'], 'Untested implementation: ' + name
        assert manifest['implementation'][name] == pin(manifest['implementation'][name]['path'])
    job = str(receipt['slurm_job_id'])
    accounting = subprocess.check_output(['sacct', '-n', '-P', '-j', job,
                                         '--format=JobIDRaw,State,ExitCode'], text=True)
    rows = [row.split('|') for row in accounting.splitlines() if row.split('|')[0] == job]
    assert len(rows) == 1 and rows[0][1:3] == ['COMPLETED', '0:0'], 'Test Slurm job did not complete successfully'
    with (manifest_path.parent / ('COLLECTOR_TEST_GATE_' + os.environ['SLURM_JOB_ID'] + '.json')).open('x') as handle:
        json.dump({'receipt': pin(receipt_path), 'manifest': pin(manifest_path),
                   'qualification_accounting': accounting, 'checked_epoch': time.time()}, handle, indent=2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='op', required=True)
    one = sub.add_parser('run'); one.add_argument('--directory', type=Path, required=True)
    two = sub.add_parser('check'); two.add_argument('--receipt', type=Path, required=True); two.add_argument('--manifest', type=Path, required=True)
    args = parser.parse_args()
    if args.op == 'run':
        run(args.directory.resolve())
    else:
        check(args.receipt.resolve(), args.manifest.resolve())
