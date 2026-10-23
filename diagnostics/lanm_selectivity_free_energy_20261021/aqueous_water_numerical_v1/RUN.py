#!/usr/bin/env python3
"""Two fixed aqueous-namespace controls; no dynamics or automatic recovery."""
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import socket
import sys
import time

PROTOCOL = 'aqueous_water_numerical_v1'
GATES = {'energy_kcal_mol': 1e-4, 'gradient_component_kcal_mol_A': 1e-3,
         'component_energy_kcal_mol': 1e-4}


def require(value, message):
    if not value:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve(strict=True)
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path']) == record, 'Changed artifact: ' + str(record))
    return Path(record['path'])


def save(path, value):
    tmp = Path(path).with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def xyz(path):
    lines = Path(path).read_text().splitlines()
    require(int(lines[0].split()[0]) == 25 and len(lines) == 26, 'Expected aquo8, no box')
    rows = []
    for index, line in enumerate(lines[1:], 1):
        p = line.split()
        require(int(p[0]) == index, 'Nonsequential XYZ atoms')
        r = {'element': p[1], 'position': list(map(float, p[2:5])),
             'type': int(p[5]), 'bonds': list(map(int, p[6:]))}
        require(all(math.isfinite(x) for x in r['position']), 'Nonfinite coordinates')
        rows.append(r)
    return rows


def urey_delta(rows):
    """Exact BIO-minus-SI correction at fixed coordinates, kcal/mol and Angstrom."""
    k, old_r, new_r = -7.60, 1.5337, 1.5537
    gradient = [[0.0] * 3 for _ in rows]
    energy = 0.0
    pairs = []
    for row in rows:
        if row['type'] != 402:
            continue
        require(len(row['bonds']) == 2, 'Water oxygen topology')
        i, j = (v-1 for v in row['bonds'])
        require(rows[i]['type'] == rows[j]['type'] == 403, 'Water hydrogen topology')
        require(i != j and 0 <= i < 25 and 0 <= j < 25, 'Invalid hydrogen indices')
        d = [a-b for a, b in zip(rows[i]['position'], rows[j]['position'])]
        r = math.sqrt(sum(v*v for v in d))
        require(r > 0, 'Coincident hydrogens')
        energy += k*((r-new_r)**2-(r-old_r)**2)
        for a in range(3):
            g = 2*k*(old_r-new_r)*d[a]/r
            gradient[i][a] += g
            gradient[j][a] -= g
        pairs.append([i, j])
    require(len(pairs) == 8 and len({i for pair in pairs for i in pair}) == 16,
            'Eight distinct water pairs required')
    return {'energy_kcal_mol': energy, 'gradient_kcal_mol_A': gradient,
            'hydrogen_pairs_zero_based': pairs, 'force_constant': k,
            'SI_equilibrium_A': old_r, 'BIO_equilibrium_A': new_r,
            'ureyunit': 1.0, 'cubic': 0.0, 'quartic': 0.0}


def metadata(config):
    require(not sys.flags.optimize, 'Native helper requires assertion checks')
    require(config['protocol'] == PROTOCOL and config['gates'] == GATES, 'Changed protocol/gates')
    require(verify(config['implementation']).resolve() == Path(__file__).resolve(), 'Wrong source')
    for item in config['authorities']:
        verify(item)
    verify(config['declaration'])
    require(config['maximum_native_programs'] == 4, 'Changed finite program inventory')
    require(config['urey_constants'] == {'unit': 1.0, 'cubic': 0.0, 'quartic': 0.0},
            'Unqualified Urey constants')
    rows = config['rows']
    require(len(rows) == 2 and [r['metal'] for r in rows] == ['La', 'Dy'], 'Changed cells')
    for row in rows:
        for field in ('xyz', 'key', 'parameters', 'reference_xyz', 'reference_key', 'reference_parameters'):
            verify(row[field])
        for field in ('analyze_output', 'gradient_output'):
            verify(row['reference_native'][field])
        a, b = xyz(row['xyz']['path']), xyz(row['reference_xyz']['path'])
        require([r['element'] for r in a] == [r['element'] for r in b], 'Changed elements/order')
        require([r['position'] for r in a] == [r['position'] for r in b], 'Changed geometry')
        require([r['bonds'] for r in a] == [r['bonds'] for r in b], 'Changed connectivity')
        mapping = {36: 402, 37: 403, 500: 900, 501: 901}
        require([r['type'] for r in a] == [mapping[r['type']] for r in b], 'Wrong namespace transfer')
        require(sum(r['type'] == (900 if row['metal'] == 'La' else 901) for r in a) == 1,
                'Wrong metal')
        for prm, pattern in ((row['parameters'], r'^ureybrad\s+74\s+73\s+74\s+(-?\S+)\s+(\S+)'),
                             (row['reference_parameters'], r'^ureybrad\s+35\s+34\s+35\s+(-?\S+)\s+(\S+)')):
            matches = re.findall(pattern, Path(prm['path']).read_text(), re.M | re.I)
            expected = (-7.6, 1.5537 if prm == row['parameters'] else 1.5337)
            require(len(matches) == 1 and tuple(map(float, matches[0])) == expected, 'Urey record differs')
        key = Path(row['key']['path']).read_text()
        p = re.findall(r'^parameters\s+(\S+)\s*$', key, re.M | re.I)
        require(p == [row['parameter_filename']] and Path(p[0]).name == p[0], 'Parameter filename')
        require(not re.search(r'^\s*(?:a-axis|b-axis|c-axis|ewald|tinfoil|solvate)\b', key, re.M | re.I),
                'Unexpected boundary or solvent')
        require(re.findall(r'^polar-eps\s+(\S+)\s*$', key, re.M | re.I) == ['1.0e-8'], 'Changed polarization tolerance')
        require(re.findall(r'^openmp-threads\s+(\S+)\s*$', key, re.M | re.I) == ['1'], 'Changed native thread count')
        for field in ('key', 'parameters', 'reference_key', 'reference_parameters'):
            body = Path(row[field]['path']).read_text()
            for name, expected in [('ureyunit', 1.0), ('urey-cubic', 0.0), ('urey-quartic', 0.0)]:
                values = re.findall(r'^\s*'+name+r'\s+(\S+)', body, re.M | re.I)
                require(len(values) <= 1 and all(float(v) == expected for v in values), 'Changed Urey global coefficient')
        urey_delta(a)
    for record in config['executables'].values():
        verify(record)
    helper_path = verify(config['native_helper'])
    spec = importlib.util.spec_from_file_location('aqueous_native_helper', helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


def comparison(actual, reference, correction):
    de = actual['total_energy_kcal_mol'] - reference['total_energy_kcal_mol']
    errors = [a-b-c for ar, br, cr in zip(actual['gradient_kcal_mol_A'],
              reference['gradient_kcal_mol_A'], correction['gradient_kcal_mol_A'])
              for a, b, c in zip(ar, br, cr)]
    require(len(errors) == 75 and all(math.isfinite(x) for x in errors), 'Gradient dimensions/nonfinite')
    ac, rc = actual['components'], reference['components']
    require(set(ac) == set(rc), 'Changed component coverage')
    require('Urey-Bradley' in ac, 'Missing Urey component')
    component_errors = {}
    for name in ac:
        require(ac[name]['interactions'] == rc[name]['interactions'], 'Changed interaction count: ' + name)
        component_errors[name] = (ac[name]['energy_kcal_mol']-rc[name]['energy_kcal_mol']
                                  -(correction['energy_kcal_mol'] if name == 'Urey-Bradley' else 0))
    energy_error = de-correction['energy_kcal_mol']
    largest = max(abs(v) for v in errors)
    passed = abs(energy_error) <= GATES['energy_kcal_mol'] and largest <= GATES['gradient_component_kcal_mol_A']
    passed = passed and all(abs(v) <= GATES['component_energy_kcal_mol'] for v in component_errors.values())
    return {'passed': passed, 'uncorrected_BIO_minus_SI_energy_kcal_mol': de,
            'known_urey_difference': correction, 'energy_remainder_kcal_mol': energy_error,
            'max_gradient_remainder_kcal_mol_A': largest,
            'gradient_remainder_RMS_kcal_mol_A': math.sqrt(sum(v*v for v in errors)/75),
            'component_energy_remainders_kcal_mol': component_errors}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    result = {'status': 'unavailable', 'protocol': PROTOCOL, 'config': pin(args.config),
              'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'host': socket.gethostname(),
              'slurm_job': os.environ.get('SLURM_JOB_ID'), 'native_programs': 0, 'rows': [],
              'bulk_physics_qualified': False, 'protein_qualified': False, 'error': None}
    try:
        helper = metadata(cfg)
        if args.preflight_only:
            result['status'] = 'metadata_preflight_passed'
        else:
            require(os.environ.get('SLURM_JOB_ID'), 'Native execution requires Slurm')
            require(all(os.environ.get(k) == '1' for k in helper.THREAD_KEYS), 'Thread allowance')
            release = json.loads(verify(cfg['root_release']).read_text())
            require(release['implementation'] == cfg['implementation'] and release['maximum_native_programs'] == 4,
                    'Unbound release')
            base_cfg = {k: v for k, v in cfg.items() if k != 'root_release'}
            require(release['config_sha256_without_release'] == hashlib.sha256(
                    json.dumps(base_cfg, sort_keys=True).encode()).hexdigest(), 'Release does not bind actual configuration')
            for row in cfg['rows']:
                target = {'id': row['id'], 'metal': row['metal'], 'status': 'unavailable', 'receipts': []}
                result['rows'].append(target)
                folder = out / row['id']
                folder.mkdir()
                for field, name in [('xyz', 'system.xyz'), ('key', 'system.key'), ('parameters', row['parameter_filename'])]:
                    shutil.copyfile(verify(row[field]), folder / name)
                reference = helper.parse(verify(row['reference_native']['analyze_output']),
                                         verify(row['reference_native']['gradient_output']),
                                         {'atoms': 25, 'derivative_test': row['reference_native']['derivative_test']})
                for name, argv in [('analyze', ['system.xyz', 'E']), ('testgrad', ['system.xyz', 'Y', 'N'])]:
                    result['native_programs'] += 1
                    save(out/'RESULT.json', result)
                    helper.invoke([cfg['executables'][name]['path']] + argv, folder,
                                  folder/(name+'.out'), target['receipts'])
                actual = helper.parse(folder/'analyze.out', folder/'testgrad.out', {'atoms': 25, 'derivative_test': False})
                target.update(native=actual, reference=reference,
                              comparison=comparison(actual, reference, urey_delta(xyz(row['xyz']['path']))))
                target['status'] = 'passed' if target['comparison']['passed'] else 'failed_mapping_gate'
                save(out/'RESULT.json', result)
                if target['status'] != 'passed':
                    break
            metadata(cfg)
            result['status'] = ('mapping_controls_passed' if len(result['rows']) == 2 and
                                all(r['status'] == 'passed' for r in result['rows']) else 'mapping_controls_failed')
    except Exception as exc:
        result['error'] = repr(exc)
    finally:
        result['elapsed_seconds'] = time.monotonic()-start
        result['ended_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
        save(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] in ('mapping_controls_passed', 'metadata_preflight_passed') else 1


if __name__ == '__main__':
    raise SystemExit(main())
