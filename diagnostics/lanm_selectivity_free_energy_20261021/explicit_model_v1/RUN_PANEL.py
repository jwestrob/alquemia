#!/usr/bin/env python3
"""Finite native Tinker panel. Serial CPU only; preserve every requested slot."""
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

PROTOCOL = 'explicit_model_v1_native_panel'
THREAD_KEYS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
               'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'

def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):
    assert pin(p['path']) == p, f'Artifact changed: {p}'
    return Path(p['path'])

def read(p):
    return json.loads(verify(p).read_text())

def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')

def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def parse(analyze_path, gradient_path, row):
    a, g = Path(analyze_path).read_text(), Path(gradient_path).read_text()
    for text in (a, g):
        assert 'Tinker' in text and 'Version 26.2' in text
        assert not re.search(r'(?i)not converged|undefined|fatal|error termination|\bNaN\b|\bInfinity\b', text)
    energies = []
    for text in (a, g):
        m = re.findall(r'Total Potential Energy\s*:\s*(' + NUMBER + ')', text)
        assert len(m) == 1
        energies.append(float(m[0]))
    assert energies[0] == energies[1], 'Energy differs between native analytic programs'
    interm = re.findall(r'Intermolecular Energy\s*:\s*(' + NUMBER + ')', a)
    assert len(interm) <= 1
    components = {}
    tail = a.split('Energy Component Breakdown :', 1)[1]
    for line in tail.splitlines():
        m = re.fullmatch(r'\s*([A-Za-z][A-Za-z /()\-]+?)\s+(' + NUMBER + r')\s+(\d+)\s*', line)
        if m:
            components[m[1].strip()] = {'energy_kcal_mol': float(m[2]), 'interactions': int(m[3])}
    assert components
    tail = g.split('Cartesian Gradient Breakdown over Individual Atoms :', 1)[1]
    arrays = {'Anlyt': [], 'Numer': []}
    for line in tail.splitlines():
        m = re.fullmatch(r'\s*(Anlyt|Numer)\s+(\d+)\s+(' + NUMBER + r')\s+(' + NUMBER +
                         r')\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s*', line)
        if m:
            assert int(m[2]) == len(arrays[m[1]]) + 1
            arrays[m[1]].append([float(m[k]) for k in (3, 4, 5)])
    assert len(arrays['Anlyt']) == row['atoms']
    assert len(arrays['Numer']) == (row['atoms'] if row['derivative_test'] else 0)
    values = energies + [v['energy_kcal_mol'] for v in components.values()]
    values += [v for arr in arrays.values() for xyz in arr for v in xyz]
    assert all(math.isfinite(v) for v in values)
    fd = None
    if row['derivative_test']:
        error = max(abs(x-y) for ax, nx in zip(arrays['Anlyt'], arrays['Numer']) for x,y in zip(ax,nx))
        fd = {'step_A': 1e-5, 'max_absolute_error_kcal_mol_A': error,
              'absolute_tolerance_kcal_mol_A': 1e-3, 'passed': error <= 1e-3,
              'representation': 'native printed eight-decimal gradient components'}
    return {'total_energy_kcal_mol': energies[0],
            'native_intermolecular_energy_kcal_mol': float(interm[0]) if interm else None,
            'components': components,
            'printed_component_sum_remainder_kcal_mol': energies[0] - sum(v['energy_kcal_mol'] for v in components.values()),
            'gradient_kcal_mol_A': arrays['Anlyt'],
            'force_kcal_mol_A': [[-v for v in xyz] for xyz in arrays['Anlyt']],
            'finite_difference': fd, 'analyze_output': pin(analyze_path), 'gradient_output': pin(gradient_path),
            'polarization': {'method': 'mutual induced dipoles', 'declared_polar_eps': 1e-8,
                             'no_native_nonconvergence_warning': True,
                             'printed_iteration_residual': None}}

def invoke(argv, cwd, output, receipts):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.monotonic()
    receipt = {'argv': [str(v) for v in argv], 'cwd': str(cwd), 'started_utc': stamp(),
               'threads': {k: os.environ[k] for k in THREAD_KEYS}, 'timeout_seconds': 300}
    receipts.append(receipt)
    try:
        with open(output, 'w') as handle:
            proc = subprocess.run(argv, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT,
                                  stdin=subprocess.DEVNULL, timeout=300, check=False)
        receipt['returncode'] = proc.returncode
        assert proc.returncode == 0, f'Native program exited {proc.returncode}'
    finally:
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt.update(ended_utc=stamp(), elapsed_seconds=time.monotonic()-start,
                       child_user_cpu_seconds=after.ru_utime-before.ru_utime,
                       child_system_cpu_seconds=after.ru_stime-before.ru_stime,
                       output=pin(output))

def inputs(config):
    assert config['protocol'] == PROTOCOL
    assert config['implementation'] == pin(__file__)
    verify(config['specification'])
    manifest = read(config['manifest'])
    for key in ('implementation', 'panel', 'published_La_SI', 'Dy_author_parameters'):
        verify(manifest[key])
    for p in config['executables'].values():
        verify(p)
    assert manifest['ALPB'] is False and manifest['implicit_solvent'] is False
    assert manifest['transfer_qualified'] is False
    rows = manifest['rows']
    assert len(rows) == 32 and len({r['id'] for r in rows}) == 32
    assert sum(r['kind'] != 'free_fragment' for r in rows) == 28
    assert sum(r['derivative_test'] for r in rows) == 4
    for row in rows:
        for key in ('xyz', 'key', 'parameters'):
            verify(row[key])
        verify(row['source_geometry']['xyz'])
        assert len(row['types']) == row['atoms'] == len(row['bonds_zero_based'])
    assert all(os.environ.get(k) == '1' for k in THREAD_KEYS)
    return manifest

def comparison(rows, results):
    by_id = {r['id']: r for r in results}
    contrasts = []
    for row in rows:
        if row['metal'] != 'La':
            continue
        other = row['id'][:-2] + 'Dy'
        pair = [by_id[row['id']], by_id[other]]
        c = {'La_id': row['id'], 'Dy_id': other, 'status': 'unavailable'}
        if all(p['status'] == 'complete' for p in pair):
            c.update(status='complete', Dy_minus_La_energy_kcal_mol=pair[1]['native']['total_energy_kcal_mol']-pair[0]['native']['total_energy_kcal_mol'])
        contrasts.append(c)
    interactions = []
    for row in rows:
        if row['fragment_id'] is None:
            continue
        c = {'complex': row['id'], 'fragment': row['fragment_id'], 'status': 'unavailable'}
        a,b = by_id[row['id']], by_id[row['fragment_id']]
        if a['status'] == b['status'] == 'complete':
            c.update(status='complete', interaction_energy_kcal_mol=a['native']['total_energy_kcal_mol']-b['native']['total_energy_kcal_mol'],
                     isolated_spherical_ion_self_energy_kcal_mol=0.0)
        interactions.append(c)
    return {'same_geometry_metal_contrasts': contrasts, 'dimer_fragment_interactions': interactions,
            'aquo_binding_free_energy': None, 'affinity_or_selectivity_prediction': False}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False)
    config = json.loads(Path(args.config).read_text())
    start = time.monotonic()
    result = {'protocol': PROTOCOL, 'status': 'unavailable', 'started_utc': stamp(),
              'config': pin(args.config), 'implementation': pin(__file__), 'host': socket.gethostname(),
              'slurm_job_id': os.environ.get('SLURM_JOB_ID'), 'native_program_invocations': 0,
              'rows': [], 'error': None}
    try:
        manifest = inputs(config)
        stage = config['stage']; assert stage in ('controls', 'panel')
        if stage == 'panel':
            assert os.environ.get('SLURM_JOB_ID'), 'Full fixed panel requires root-owned Slurm submission'
            parent = read(config['controls'])
            assert parent['status'] == 'controls_passed' and parent['manifest'] == config['manifest']
            assert parent['implementation'] == pin(__file__)
            old = {r['id']: r for r in parent['rows'] if r['status'] == 'complete'}
            assert len(old) == 4
        else:
            old = {}
        result.update(stage=stage, manifest=config['manifest'], rows=[
            {'id': r['id'], 'status': 'not_requested' if stage == 'controls' and not r['derivative_test'] else 'unavailable',
             'input': r, 'native': None, 'receipts': []} for r in manifest['rows']])
        stopped = False
        for row, target in zip(manifest['rows'], result['rows']):
            if target['status'] == 'not_requested' or stopped:
                continue
            try:
                if row['id'] in old:
                    previous = old[row['id']]
                    assert previous['input'] == row
                    a = verify(previous['native']['analyze_output'])
                    g = verify(previous['native']['gradient_output'])
                    target.update(native=parse(a,g,row), reused_from=config['controls'])
                elif stage == 'controls' and row['id'] == 'water_central_La':
                    saved = config['saved_La_water']
                    assert saved['xyz'] == row['xyz'] and saved['key'] == row['key'] and saved['parameters'] == row['parameters']
                    target.update(native=parse(verify(saved['analyze_output']),verify(saved['gradient_output']),row),
                                  reused_from='original completed one-thread smoke; exact output pins')
                else:
                    cwd = out / row['id']; cwd.mkdir()
                    for key, name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
                        shutil.copyfile(verify(row[key]), cwd/name)
                        assert pin(cwd/name)['sha256'] == row[key]['sha256']
                    for name, argv in (
                        ('analyze', [config['executables']['analyze']['path'], 'system.xyz', 'E']),
                        ('gradient', [config['executables']['testgrad']['path'], 'system.xyz', 'Y'] +
                         (['Y','0.00001'] if row['derivative_test'] else ['N']))):
                        result['native_program_invocations'] += 1
                        invoke(argv, cwd, cwd/(name+'.out'), target['receipts'])
                    target['native'] = parse(cwd/'analyze.out',cwd/'gradient.out',row)
                target['status'] = 'complete'
                fd = target['native']['finite_difference']
                if fd is not None and not fd['passed']:
                    target['status'] = 'failed_derivative_gate'; stopped = True
            except Exception as exc:
                target.update(status='unavailable', error=repr(exc)); stopped = True
            save(out/'RESULT.json', result)
        selected = [r for r in result['rows'] if r['status'] != 'not_requested']
        success = all(r['status']=='complete' for r in selected)
        result['status'] = ('controls_passed' if stage=='controls' else 'panel_complete') if success else 'incomplete_or_failed'
        result['comparison'] = comparison(manifest['rows'], result['rows'])
        result['complete_rows'] = sum(r['status']=='complete' for r in result['rows'])
        result['transfer_qualified'] = False
        result['caveats'] = ['Dy common-water/donor transfer is unqualified.',
                             'Fixed clusters are not explicit bulk-water sampling or free energies.',
                             'No protein, ALPB, GPU, parameter fitting or optimization.',
                             'FD checks use printed gradients; polarization residuals are not printed.']
    except Exception as exc:
        result['error'] = repr(exc)
    finally:
        result.update(ended_utc=stamp(), elapsed_seconds=time.monotonic()-start)
        save(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] in ('controls_passed','panel_complete') else 1

if __name__ == '__main__':
    raise SystemExit(main())
