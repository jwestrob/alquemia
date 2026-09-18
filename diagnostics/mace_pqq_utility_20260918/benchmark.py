#!/usr/bin/env python3
"""Matched prepared-PQQ timing, delegating science to the existing frozen runners."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import socket
import statistics
import subprocess
import sys
import time


def read(path):
    return json.loads(Path(path).read_text())


def record(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(ref):
    if record(ref['path']) != ref:
        raise ValueError('artifact changed: ' + ref['path'])
    return Path(ref['path'])


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, sort_keys=True)
        f.write('\n')


def classify(value, ca, la):
    return 'Ca' if value <= ca else 'La' if value >= la else 'inconclusive'


def dft_score(energies, baseline):
    delta = energies['Ca']-energies['La']
    raw = delta*627.509474
    gauge = baseline['aquo_reporting_gauge']['delta_E_aquo_hartree']
    b = baseline['calibration']
    return {'score': (delta-gauge)*627.509474, 'raw_R_kcal_mol': raw,
            'class': classify(raw,b['U_max_Ca_kcal_mol'],b['L_min_La_kcal_mol']),
            'classification_scale': 'released_primary_R_bands'}


def prepare(root, output, frozen_mace):
    root, out = Path(root).resolve(), Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    impl = out/'implementation'
    impl.mkdir()
    for name in ('run_orca_task_manifest.py', 'render_orca_runtime_input.py', 'result_protocol.py'):
        shutil.copyfile(root/'scripts'/name, impl/name)
    shutil.copyfile(__file__, impl/'benchmark.py')
    frozen = Path(frozen_mace).resolve()
    base = root/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json'
    work = root/'workspaces/mace_omol_20260917'
    reference = work/'masked_calibration_v1/reference.json'
    prep = work/'intact_panel_prepared_v2/preparation_manifest.json'
    baseline, masked = read(base), read(reference)
    rows = {r['case_id']: r for r in read(prep)['rows']}
    scores = {r['case_id']: r for r in masked['scores']}
    cases = []
    for d in baseline['scores']:
        case = d['panel_id']; m = scores[case]; p = rows[case]
        assert d['class'] == m['expected_class'] == p['expected_class']
        assert p['status'] == 'prepared' and m['evaluation_role'] == 'calibration'
        endpoints = {}
        for metal in ('Ca', 'La'):
            artifacts = d['artifacts'][metal]
            old_receipt = read(Path(artifacts['output']['path'] + '.execution.json'))
            inp = old_receipt['artifacts']['template_input']
            endpoints[metal] = {'input': inp, 'xyz': artifacts['xyz']}
            verify(inp); verify(artifacts['xyz'])
        cases.append({'case_id': case, 'expected_class': d['class'],
                      'sequence_group': m['sequence_accession_group'],
                      'evidence_role': 'consumed_calibration', 'mace_preparation': p['preparation'],
                      'dft_endpoints': endpoints, 'archived_DFT_S': d['S_aquo_gauge_kcal_mol'],
                      'archived_MACE_R': m['two_call_model_kcal']})
    assert len(cases) == 25 and len({c['case_id'] for c in cases}) == 25
    manifest = {'schema': 'alquemia.pqq_utility_timing.v1', 'cases': cases,
                'agreement': record(root/'diagnostics/mace_pqq_utility_20260918/GOAL.md'),
                'baseline': record(base), 'calibration': record(reference),
                'factorization': masked['factorization'],
                'development': record(work/'charge_ablation_development_v2/collection_job_1200828.json'),
                'frozen_MACE': {p.name: record(p) for p in frozen.glob('*.py')},
                'implementation': {p.name: record(p) for p in impl.glob('*.py')},
                'driver_python': record(sys.executable),
                'orca': record('/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca'),
                'expected_host': 'node-128-512g-8gpu-1',
                'MACE_CPUs': 16, 'DFT_CPUs': 32, 'DFT_MPI_ranks_per_endpoint': 16,
                'score_reproduction_tolerance': 0.01, 'speed_target': 1.5,
                'fresh_endpoint_counts': {'MACE': 50, 'DFT': 50},
                'new_calibration': False, 'baseline_changed': False}
    write(out/'manifest.json', manifest)
    return {'manifest': record(out/'manifest.json'), 'cases': len(cases)}


def run_dft(m, c, out):
    from run_orca_task_manifest import run_manifest
    started = time.monotonic(); tasks = []
    for metal, e in c['dft_endpoints'].items():
        inp, xyz = verify(e['input']), verify(e['xyz'])
        shutil.copyfile(inp, out/inp.name); shutil.copyfile(xyz, out/xyz.name)
        tasks.append({'task_id': metal, 'input': record(out/inp.name),
                      'xyz': record(out/xyz.name), 'output_path': str(out/(metal+'.out'))})
    manifest = out/'manifest.json'
    write(manifest, {'tasks': tasks, 'protocol_id': read(verify(m['baseline']))['protocol_id'],
                    'source_case': c, 'execution_policy': {
                        'task_runner': m['implementation']['run_orca_task_manifest.py'],
                        'runtime_renderer': m['implementation']['render_orca_runtime_input.py']}})
    prepared = time.monotonic()
    results = run_manifest(manifest, orca_path=verify(m['orca']), workers=2,
                           nprocs=m['DFT_MPI_ranks_per_endpoint'],
                           expected_runner_sha256=m['implementation']['run_orca_task_manifest.py']['sha256'])
    if any(r['status'] != 'complete' for r in results):
        raise ValueError('fresh converged endpoints required: '+str(results))
    calculated = time.monotonic(); energies = {}
    for metal in ('Ca', 'La'):
        path = out/(metal+'.out'); receipt = read(str(path)+'.execution.json')
        assert receipt['normal_termination'] and receipt['scf_converged']
        verify(receipt['artifacts']['output'])
        values = re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+0-9.]+)', path.read_text())
        if len(values) != 1:
            raise ValueError('expected one final single-point energy')
        energies[metal] = float(values[0])
    baseline = read(verify(m['baseline']))
    result = {'status': 'complete', 'energy_hartree': energies, **dft_score(energies,baseline),
              'score_units': 'kcal/mol', 'preparation_seconds': prepared-started,
              'execution_seconds': calculated-prepared,
              'report_seconds': time.monotonic()-calculated}
    write(out/'result.json', result)
    return result


def run_mace(m, c, out):
    # Imports come from the immutable, previously qualified implementation.
    from mace_omol_prepared import prepare as prepare_case, report
    from mace_hybrid import execute
    from mace_file_checks import cached_file_checks

    @cached_file_checks
    def operation():
        started = time.monotonic(); case_dir = out/'score'
        prepare_case(verify(c['mace_preparation']), verify(m['development']),
                     verify(m['agreement']), case_dir, factorization=verify(m['factorization']))
        prepared = time.monotonic()
        execute(case_dir/'manifest.json', 'native')
        calculated = time.monotonic()
        report(case_dir/'manifest.json', out/'report', calibration=verify(m['calibration']))
        r = read(out/'report/result.json')
        if r['status'] != 'complete' or not r['numerical_gate_pass']:
            raise ValueError('incomplete or numerically invalid MACE result')
        result = {'status': 'complete', 'score': r['R_mask_model_kcal'],
                  'class': r['calibrated_class'], 'score_units': r['score_unit'],
                  'preparation_seconds': prepared-started,
                  'execution_seconds': calculated-prepared,
                  'report_seconds': time.monotonic()-calculated,
                  'scientific_report': record(out/'report/result.json')}
        write(out/'result.json', result)
        return result
    return operation()


def execute(manifest, method):
    m = read(manifest); root = Path(manifest).resolve().parent
    if not os.environ.get('SLURM_JOB_ID') or socket.gethostname().split('.')[0] != m['expected_host']:
        raise ValueError('matched allocated host required')
    if int(os.environ['SLURM_CPUS_ON_NODE']) != m[method+'_CPUs']:
        raise ValueError('unexpected CPU allocation')
    for ref in [m['agreement'], m['driver_python'], *m['implementation'].values(), *m['frozen_MACE'].values()]:
        verify(ref)
    if record(sys.executable) != m['driver_python']:
        raise ValueError('recorded driver required')
    # Avoid imports from the mutable worktree or paused experimental prototype.
    sys.path.insert(0, str(verify(m['frozen_MACE']['mace_omol_prepared.py']).parent)
                    if method == 'MACE' else str(root/'implementation'))
    destination = root/method
    destination.mkdir(exist_ok=False)
    failures = 0
    for c in m['cases']:
        out = destination/c['case_id']; out.mkdir()
        start = time.monotonic(); cpu = time.process_time()
        receipt = {'case_id': c['case_id'], 'method': method, 'manifest': record(manifest),
                   'job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
                   'allocated_cpus': int(os.environ['SLURM_CPUS_ON_NODE']),
                   'allocated_memory_MiB': os.environ.get('SLURM_MEM_PER_NODE'),
                   'started_UTC': dt.datetime.now(dt.timezone.utc).isoformat()}
        try:
            receipt['result'] = (run_mace if method == 'MACE' else run_dft)(m,c,out)
            receipt['status'] = 'complete'
        except Exception as e:
            import traceback
            traceback.print_exc()
            receipt.update(status='failed', error=repr(e)); failures += 1
        receipt.update(wall_seconds=time.monotonic()-start,
                       controller_CPU_seconds=time.process_time()-cpu)
        write(out/'timing.json', receipt)
        print(json.dumps({'case':c['case_id'],'method':method,'status':receipt['status'],
                          'wall_seconds':receipt['wall_seconds']}), flush=True)
        # Repeated preparation/backend errors should not generate 25 identical attempts.
        if failures >= 2:
            break
    return {'status': 'failed' if failures else 'complete', 'failures': failures}


def collect(manifest, output):
    m = read(manifest); root = Path(manifest).resolve().parent; rows = []
    for c in m['cases']:
        row = {'case_id':c['case_id'],'expected_class':c['expected_class']}
        for method in ('DFT','MACE'):
            path = root/method/c['case_id']/'timing.json'
            r = read(path) if path.exists() else {'status':'unavailable'}
            if r['status'] == 'complete':
                old = c['archived_DFT_S' if method=='DFT' else 'archived_MACE_R']
                delta = r['result']['score']-old
                r.update(score_difference=delta, reproduced=abs(delta)<=m['score_reproduction_tolerance'],
                         literal_class_correct=r['result']['class']==c['expected_class'])
            row[method] = r
        rows.append(row)
    pairs = [r for r in rows if all(r[k]['status']=='complete' for k in ('DFT','MACE'))]
    result = {'manifest':record(manifest),'status':'complete' if len(pairs)==25 else 'incomplete',
              'complete_pairs':len(pairs),'rows':rows,'baseline_changed':False}
    if pairs:
        totals = {k:sum(r[k]['wall_seconds'] for r in pairs) for k in ('DFT','MACE')}
        result['paired_total_seconds'] = totals
        result['total_speedup'] = totals['DFT']/totals['MACE']
        result['median_pair_speedup'] = statistics.median(r['DFT']['wall_seconds']/r['MACE']['wall_seconds'] for r in pairs)
    write(output,result)
    return {k:v for k,v in result.items() if k!='rows'}


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare');a.add_argument('--root',required=True);a.add_argument('--output',required=True);a.add_argument('--frozen-mace',required=True)
    a=sub.add_parser('execute');a.add_argument('--manifest',required=True);a.add_argument('--method',choices=('DFT','MACE'),required=True)
    a=sub.add_parser('collect');a.add_argument('--manifest',required=True);a.add_argument('--output',required=True)
    args=vars(p.parse_args());command=args.pop('command')
    result=globals()[command](**args)
    print(json.dumps(result,indent=2))
    if command=='execute' and result['status']!='complete':sys.exit(1)
