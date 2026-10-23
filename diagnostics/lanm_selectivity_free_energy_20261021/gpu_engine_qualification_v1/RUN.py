#!/usr/bin/env python3
"""Bounded GPU parity qualification; immutable inputs, no native fallback."""
import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import socket
import subprocess
import time
import xml.etree.ElementTree as ET

PROTOCOL = 'gpu_engine_qualification_v1'
SOURCE_COMMIT = '44bcd7c898e68f6827c112df79b8361e443730d2'
CONTROL_IDS = ('acetate_bi_central_La', 'acetate_bi_central_Dy',
               'water_central_La', 'water_central_Dy')
SUITE_NAMES = {'PolPair-Ewald', 'PolPair-NonEwald',
               'PolPair-Ewald-A', 'PolPair-NonEwald-A'}
GATES = {'energy_absolute_kcal_mol': 1e-4,
         'force_max_component_kcal_mol_A': 1e-3,
         'finite_difference_max_component_kcal_mol_A': 1e-3}
NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
FAILURE = re.compile(r'not converged|failed to converge|undefined|fatal|'
                     r'error termination|uncaught exception|segmentation fault|'
                     r'cuda error|\b(?:nan|inf|infinity)\b', re.I)
THREAD_KEYS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
               'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def pin(path):
    p = Path(path).resolve(strict=True)
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(isinstance(record, dict) and set(record) == {'path', 'sha256'},
            'Unbound or malformed artifact pin')
    require(pin(record['path']) == record, f'Artifact changed: {record}')
    return Path(record['path'])


def load(record):
    return json.loads(verify(record).read_text())


def save(path, data):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def number(value):
    x = float(value.replace('D', 'E').replace('d', 'e'))
    require(math.isfinite(x), 'Nonfinite native number')
    return x


def parse_gradient(text, atoms):
    require(not FAILURE.search(text), 'Native failure/nonfinite marker')
    energies = re.findall(r'^\s*Total Potential Energy\s*:\s*(' + NUMBER +
                          r')\s+Kcal/mole\s*$', text, re.M)
    require(len(energies) == 1 and text.count('Total Potential Energy') == 1,
            'Expected exactly one finite native total energy in Kcal/mole')
    arrays = {'Anlyt': [], 'Numer': []}
    norms = {'Anlyt': [], 'Numer': []}
    pattern = (r'\s*(Anlyt|Numer)\s+(\d+)\s+(' + NUMBER + r')\s+(' +
               NUMBER + r')\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s*')
    for line in text.splitlines():
        match = re.fullmatch(pattern, line)
        if match:
            kind, index = match[1], int(match[2])
            require(index == len(arrays[kind]) + 1, 'Duplicate/nonsequential gradient atom')
            arrays[kind].append([number(match[k]) for k in (3, 4, 5)])
            norms[kind].append(number(match[6]))
    for kind in arrays:
        require(len(arrays[kind]) == atoms, f'Incomplete {kind} gradient array')
        summary = re.findall(kind + r'\s+Total Gradient Norm Value\s*:?\s*(' + NUMBER + ')', text)
        require(len(summary) == 1, f'Missing/duplicate {kind} final norm summary')
        require(number(summary[0]) >= 0, 'Negative gradient norm')
    fd = max(abs(a-b) for aa, bb in zip(arrays['Anlyt'], arrays['Numer'])
             for a, b in zip(aa, bb))
    return {'total_energy_kcal_mol': number(energies[0]),
            'gradient_kcal_mol_A': arrays['Anlyt'],
            'numerical_gradient_kcal_mol_A': arrays['Numer'],
            'force_kcal_mol_A': [[-v for v in xyz] for xyz in arrays['Anlyt']],
            'finite_difference': {
                'full_separation_A': 1e-5, 'one_sided_displacement_A': 5e-6,
                'max_absolute_component_error_kcal_mol_A': fd,
                'passed': fd <= GATES['finite_difference_max_component_kcal_mol_A']}}


def compare(native, reference, atoms):
    require(len(native['force_kcal_mol_A']) == atoms, 'Native atom count')
    require(len(reference['force_kcal_mol_A']) == atoms, 'Reference atom count')
    differences = []
    for actual, ref in zip(native['force_kcal_mol_A'], reference['force_kcal_mol_A']):
        require(len(actual) == len(ref) == 3, 'Reference Cartesian dimensions')
        require(all(math.isfinite(v) for v in ref + actual), 'Nonfinite force component')
        differences.extend(a-b for a, b in zip(actual, ref))
    delta = native['total_energy_kcal_mol'] - reference['total_energy_kcal_mol']
    require(math.isfinite(delta), 'Nonfinite energy difference')
    maximum = max(abs(v) for v in differences)
    return {'energy_GPU_minus_CPU_kcal_mol': delta,
            'energy_passed': abs(delta) <= GATES['energy_absolute_kcal_mol'],
            'force_max_absolute_component_error_kcal_mol_A': maximum,
            'force_component_RMS_error_kcal_mol_A': math.sqrt(sum(v*v for v in differences)/(3*atoms)),
            'force_RMS_denominator': 3*atoms,
            'force_passed': maximum <= GATES['force_max_component_kcal_mol_A'],
            'finite_difference_passed': native['finite_difference']['passed'],
            'passed': abs(delta) <= GATES['energy_absolute_kcal_mol'] and
                      maximum <= GATES['force_max_component_kcal_mol_A'] and
                      native['finite_difference']['passed']}


def parse_suite(xml_text, console_text):
    require(not FAILURE.search(console_text), 'Native suite failure marker')
    doc = ET.fromstring(xml_text)
    require(doc.tag == 'Catch2TestRun', 'Unexpected suite XML root')
    require(doc.get('filters') == '[polpair]', 'Unexpected suite filter')
    cases = doc.findall('.//TestCase')
    require(len(cases) == 4 and {c.get('name') for c in cases} == SUITE_NAMES,
            'Suite must execute exactly the four POLPAIR cases')
    for case in cases:
        results = case.findall('OverallResult')
        require(len(results) == 1 and results[0].get('success') == 'true',
                f'Failed native case: {case.get("name")}')
    totals = doc.findall('OverallResults')
    counts = doc.findall('OverallResultsCases')
    require(len(totals) == len(counts) == 1, 'Missing suite totals')
    for item in (totals[0], counts[0]):
        require(int(item.get('failures', '-1')) == 0 and
                int(item.get('expectedFailures', '-1')) == 0, 'Native suite failures')
    require(int(totals[0].get('successes', '0')) > 0, 'No native assertions')
    require(int(counts[0].get('successes', '0')) == 4, 'Incomplete case count')
    require(not doc.findall('.//Exception') and not doc.findall('.//FatalErrorCondition'),
            'Native suite exception')
    for output in doc.findall('.//StdOut') + doc.findall('.//StdErr'):
        require(not FAILURE.search(output.text or ''), 'Native failure in captured suite output')
    return {'passed': True, 'cases': sorted(SUITE_NAMES),
            'successful_assertions': int(totals[0].get('successes')),
            'test_case_count': 4, 'failed_assertions': 0,
            'native_tolerances_changed': False}


def inspect_inputs(config, bound):
    require(config['protocol'] == PROTOCOL, 'Wrong protocol')
    require(config['gates'] == GATES, 'Changed gates')
    require(config['control_ids'] == list(CONTROL_IDS), 'Changed control panel/order')
    require(config['timeouts_seconds'] == {'suite': 900, 'control': 300}, 'Changed timeout')
    require(verify(config['adapter']).resolve() == Path(__file__).resolve(), 'Wrong adapter')
    verify(config['specification'])
    controls = load(config['cpu_controls_result'])
    require(controls['status'] == 'controls_passed' and controls['complete_rows'] == 4,
            'CPU controls not admitted')
    rows = {r['id']: r for r in controls['rows'] if r['status'] == 'complete'}
    require(set(rows) == set(CONTROL_IDS), 'CPU completed panel mismatch')
    for row_id in CONTROL_IDS:
        row = rows[row_id]
        require(row['input']['id'] == row_id, 'CPU row identity mismatch')
        for label in ('xyz', 'key', 'parameters'):
            require(row['input'][label] == config['control_input_pins'][row_id][label],
                    'Control input provenance mismatch')
            verify(row['input'][label])
        require(row['input']['derivative_test'] and row['native']['finite_difference']['passed'],
                'Unqualified CPU gradient reference')
        text = verify(row['native']['gradient_output']).read_text()
        replay = parse_gradient(text, row['input']['atoms'])
        require(replay['total_energy_kcal_mol'] == row['native']['total_energy_kcal_mol'],
                'CPU energy fails exact native replay')
        require(replay['force_kcal_mol_A'] == row['native']['force_kcal_mol_A'],
                'CPU force fails exact native replay')
    source = config['source']
    require(source['commit'] == SOURCE_COMMIT, 'Wrong source commit')
    root = Path(source['path'])
    for relative, commit in [('.', SOURCE_COMMIT)] + list(source['submodules'].items()):
        directory = root / relative
        head = subprocess.check_output(['git', '-C', str(directory), 'rev-parse', 'HEAD'], text=True).strip()
        require(head == commit, f'Wrong source HEAD: {directory}')
        difference = subprocess.check_output(
            ['git', '-C', str(directory), 'diff', '--binary', '--no-ext-diff',
             '--ignore-submodules=all', 'HEAD', '--'])
        expected = source['working_tree_diff_sha256'].get(relative)
        if bound or expected is not None:
            require(hashlib.sha256(difference).hexdigest() == expected,
                    f'Unbound/changed tracked source patch: {directory}')
    for artifact in source['artifacts']:
        verify(artifact)
    if bound:
        require(config.get('binding_status') == 'root_bound', 'Executable/build pins are not bound')
        build = config['completed_build']
        for label in ('receipt', 'cmake_cache', 'compile_commands'):
            verify(build[label])
        for artifact in build['additional_provenance']:
            verify(artifact)
        receipt = load(build['receipt'])
        require(receipt.get('status') == 'build_passed', 'Build receipt not passed')
        require(receipt.get('source_commit') == SOURCE_COMMIT, 'Build source mismatch')
        require(receipt.get('source_working_tree_diff_sha256') ==
                source['working_tree_diff_sha256'], 'Build/source patch binding mismatch')
        cache = verify(build['cmake_cache']).read_text()
        for setting in ('GPU_LANG:STRING=CUDA', 'HOST:BOOL=OFF',
                        'PREC:STRING=double', 'COMPUTE_CAPABILITY:STRING=90'):
            require(setting in cache.splitlines(), f'Wrong build setting: {setting}')
        for label in ('all_tests', 'tinker9'):
            executable = verify(config['executables'][label])
            require(os.access(executable, os.X_OK), 'Nonexecutable native binary')
            require(receipt.get('executables', {}).get(label) == config['executables'][label],
                    'Build receipt/executable mismatch')
    return rows


def invoke(argv, cwd, output, timeout, result):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.monotonic()
    receipt = {'argv': [str(x) for x in argv], 'cwd': str(cwd),
               'started_utc': stamp(), 'timeout_seconds': timeout}
    result['native_invocations'].append(receipt)
    try:
        with open(output, 'w') as handle:
            process = subprocess.run(argv, cwd=cwd, stdin=subprocess.DEVNULL,
                                     stdout=handle, stderr=subprocess.STDOUT,
                                     timeout=timeout, check=False)
        receipt['returncode'] = process.returncode
        require(process.returncode == 0, f'Native return code {process.returncode}')
    finally:
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt.update({'ended_utc': stamp(), 'elapsed_seconds': time.monotonic()-start,
                        'child_user_CPU_seconds': after.ru_utime-before.ru_utime,
                        'child_system_CPU_seconds': after.ru_stime-before.ru_stime})
        if Path(output).exists():
            receipt['output'] = pin(output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--saved-preflight-only', action='store_true')
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    result = {'protocol': PROTOCOL, 'status': 'started', 'started_utc': stamp(),
              'config': pin(args.config), 'adapter': pin(__file__),
              'hostname': socket.gethostname(), 'native_invocations': [],
              'native_suite': {'status': 'not_run'},
              'rows': [{'id': x, 'status': 'not_run'} for x in CONTROL_IDS],
              'physical_accuracy_qualified': False, 'production_promoted': False,
              'gpu_execution_performed': False,
              'allocated_CPUs': None, 'allocated_GPUs': None,
              'allocation': {k: os.environ.get(k) for k in
                  ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_JOB_PARTITION',
                   'CUDA_VISIBLE_DEVICES', 'SLURM_JOB_GPUS')},
              'thread_environment': {k: os.environ.get(k) for k in THREAD_KEYS}}
    save(output/'START.json', result)
    try:
        config = load(result['config'])
        controls = inspect_inputs(config, bound=not args.saved_preflight_only)
        if args.saved_preflight_only:
            result.update(status='saved_preflight_passed', gpu_execution_performed=False)
            return 0
        require(os.environ.get('SLURM_JOB_ID'), 'Actual Slurm allocation required')
        require(os.environ.get('SLURM_CPUS_PER_TASK') == '4', 'Require exactly four allocated CPUs')
        require(os.environ.get('SLURM_JOB_PARTITION') == 'gpu_h200', 'Require H200 partition')
        visible = os.environ.get('CUDA_VISIBLE_DEVICES', '').split(',')
        require(len(visible) == 1 and visible[0] not in ('', '-1', 'NoDevFiles'),
                'Require exactly one allocated visible GPU')
        require(all(os.environ.get(k) == '1' for k in THREAD_KEYS), 'Require thread caps of one')
        result['allocated_CPUs'], result['allocated_GPUs'] = 4, 1
        result['gpu_execution_performed'] = False
        inventory = subprocess.run(['nvidia-smi', '--query-gpu=index,name,uuid',
                                    '--format=csv,noheader'], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True, timeout=30, check=True)
        (output/'GPU_INVENTORY.txt').write_text(inventory.stdout + inventory.stderr)
        require('H200' in inventory.stdout, 'H200 absent from host inventory')
        result['gpu_inventory'] = pin(output/'GPU_INVENTORY.txt')
        suite = output/'native_polpair'
        suite.mkdir()
        result['native_suite'] = {'status': 'running'}
        xml = suite/'report.xml'
        result['gpu_execution_performed'] = True
        invoke([str(verify(config['executables']['all_tests'])), '[polpair]',
                '--reporter', 'xml', '--out', str(xml)], suite, suite/'console.out', 900, result)
        suite_result = parse_suite(xml.read_text(), (suite/'console.out').read_text())
        result['native_suite'] = dict(suite_result, status='passed',
                                      xml=pin(xml), console=pin(suite/'console.out'))
        save(output/'RESULT.json', result)
        for entry in result['rows']:
            control = controls[entry['id']]
            directory = output/entry['id']
            directory.mkdir()
            entry['status'] = 'running'
            entry['input'] = control['input']
            for label, name in (('xyz', 'system.xyz'), ('key', 'system.key'),
                                ('parameters', 'amoeba09.prm')):
                artifact = control['input'][label]
                shutil.copyfile(verify(artifact), directory/name)
                require(pin(directory/name)['sha256'] == artifact['sha256'], 'Copied input mismatch')
            invoke([str(verify(config['executables']['tinker9'])), 'testgrad',
                    'system.xyz', 'Y', 'Y', '0.00001'], directory,
                   directory/'gradient.out', 300, result)
            entry['native'] = parse_gradient((directory/'gradient.out').read_text(),
                                              control['input']['atoms'])
            entry['native']['gradient_output'] = pin(directory/'gradient.out')
            for label, name in (('xyz', 'system.xyz'), ('key', 'system.key'),
                                ('parameters', 'amoeba09.prm')):
                require(pin(directory/name)['sha256'] == control['input'][label]['sha256'],
                        'Native execution changed copied input')
            entry['comparison'] = compare(entry['native'], control['native'], control['input']['atoms'])
            entry['status'] = 'passed' if entry['comparison']['passed'] else 'failed_gate'
            save(output/'RESULT.json', result)
            require(entry['comparison']['passed'], f'Failed GPU qualification gate: {entry["id"]}')
        inspect_inputs(config, bound=True)
        verify(result['config'])
        result.update(status='gpu_parity_passed', gpu_execution_performed=True,
                      source_and_inputs_unchanged_after_execution=True)
        return 0
    except Exception as exc:
        result.update(status='failed', error_type=type(exc).__name__, error=str(exc))
        for row in result['rows']:
            if row['status'] == 'running':
                row['status'] = 'failed_execution_or_parse'
        if result['native_suite']['status'] == 'running':
            result['native_suite']['status'] = 'failed_execution_or_parse'
        return 1
    finally:
        result.update(ended_utc=stamp(), elapsed_seconds=time.monotonic()-started,
                      native_program_invocation_count=len(result['native_invocations']),
                      complete_controls=sum(row['status'] == 'passed' for row in result['rows']))
        if result['allocated_CPUs'] is not None:
            result['allocated_CPU_seconds_adapter_interval'] = 4*result['elapsed_seconds']
            result['allocated_GPU_seconds_adapter_interval'] = result['elapsed_seconds']
        save(output/'RESULT.json', result)
        print(json.dumps({'status': result['status'], 'result': str(output/'RESULT.json'),
                          'native_program_invocations': len(result['native_invocations'])}))


if __name__ == '__main__':
    raise SystemExit(main())
