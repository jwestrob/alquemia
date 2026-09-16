"""Contained exact-density coupling and permanent-field response diagnostics."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import fcntl
import math
import os
from pathlib import Path
import re
import shutil

import numpy as np

from affordable_common import (BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, digest,
                               energy, paired, read_json, record, verify, write_new)
from affordable_solver import run_command
from affordable_workflow import dry_run, execute
from affordable_environment import mbis_charges
from global_electrostatic import CASES, TABI_COULOMB_KCAL_A, collect_endpoints, validate_vacuum_output

PROTOCOL = 'native_r2scan3c_permanent_field_density_diagnostic_v1'


def prepare(states_manifest, endpoint_manifest, agreement, output):
    states_manifest, endpoint_manifest, output = map(Path, (states_manifest, endpoint_manifest, output))
    if output.exists():
        raise InvalidArtifact('existing preparation is immutable')
    sm, qm = read_json(states_manifest), read_json(endpoint_manifest)
    if sm['endpoint_manifest'] != record(endpoint_manifest):
        raise InvalidArtifact('state and endpoint manifests differ')
    verified = collect_endpoints(endpoint_manifest)
    if verified['status'] != 'endpoints_complete':
        raise InvalidArtifact('saved quantum endpoints are not verified')
    output = output.resolve(); output.mkdir(parents=True)
    tasks, potentials = [], []
    for old in qm['tasks']:
        key = old['task_id']; state_record = sm['states'][key]
        state = read_json(verify(state_record))
        quality = read_json(verify(state['charge_quality']['receipt']))
        d = output / key; d.mkdir()
        shutil.copyfile(verify(old['xyz']), d / 'core.xyz')
        env = state['environment_atoms']
        pc = d / 'environment.pc'
        pc.write_text(str(len(env)) + '\n' + ''.join(
            ' '.join(format(v, '.17g') for v in (a['charge_e'], *a['xyz_A'])) + '\n' for a in env))
        inp = d / 'endpoint.inp'
        inp.write_text('! r2SCAN-3c NoAutostart DefGrid3 MBIS\n'
                       '%method\n MBIS_LARGEPRINT true\n DoEQ false\nend\n'
                       '%pointcharges "environment.pc"\n'
                       f"* xyzfile {old['charge']} 1 core.xyz\n")
        tasks.append({'task_id': key, 'case': old['case'], 'metal': old['metal'],
                      'charge': old['charge'], 'multiplicity': 1, 'input': record(inp),
                      'xyz': record(d / 'core.xyz'), 'output_path': str(d / 'endpoint.out'),
                      'pointcharges': record(pc), 'source_state': state_record,
                      'vacuum_task': old})
        p = d / 'vacuum_potential'; p.mkdir()
        for field, filename in [('gbw', 'endpoint.runtime.gbw'), ('density', 'endpoint.runtime.densities')]:
            shutil.copyfile(verify(quality[field]), p / filename)
        points = p / 'points_bohr.xyz'
        coords = np.array([a['xyz_A'] for a in env]) / BOHR_TO_A
        points.write_text(str(len(coords)) + '\n' + ''.join(' '.join(f'{v:.12f}' for v in row) + '\n' for row in coords))
        potentials.append({'task_id': key, 'state': state_record, 'points': record(points),
                           'gbw': record(p / 'endpoint.runtime.gbw'),
                           'density': record(p / 'endpoint.runtime.densities'),
                           'source_quality': state['charge_quality']['receipt'],
                           'utility': quality['utility'], 'directory': str(p)})
    for case in CASES:
        pair = {t['metal']: t for t in tasks if t['case'] == case}
        paired(pair['La']['xyz']['path'], pair['Ca']['xyz']['path'], pair['La']['charge'], pair['Ca']['charge'])
        if pair['La']['pointcharges']['sha256'] != pair['Ca']['pointcharges']['sha256']:
            raise InvalidArtifact('paired permanent charges differ')
    manifest = {'schema_version': 'alquemia.density_embedding.v1', 'protocol_id': PROTOCOL,
                'tasks': tasks, 'orca': qm['orca'], 'execution_policy': qm['execution_policy'],
                'agreement': record(agreement), 'source_states': record(states_manifest),
                'source_quantum': record(endpoint_manifest), 'implementation': record(__file__),
                'energy_expression': 'E_embedded - E_vacuum - sum(Q * phi_vacuum_density)',
                'solvation': None, 'absolute_reference': None, 'classification': None,
                'compute_budget': None, 'wall_time_limit': None}
    write_new(output / 'manifest.json', manifest)
    write_new(output / 'potential_manifest.json', {'protocol_id': PROTOCOL, 'tasks': potentials,
              'agreement': record(agreement), 'quantum_manifest': record(output / 'manifest.json'),
              'implementation': record(__file__), 'units': 'points in Bohr, potential in atomic units'})
    return {'status': 'prepared', 'quantum_tasks': len(tasks), 'potential_tasks': len(potentials)}


def parse_potential(path, points):
    values = np.loadtxt(path, skiprows=1)
    if values.shape != (len(points), 4) or not np.isfinite(values).all():
        raise InvalidArtifact('invalid potential output dimensions/values')
    if np.allclose(values[:, :3], points, atol=1e-5, rtol=0):
        return values[:, 3]
    if np.allclose(values[:, 1:], points, atol=1e-5, rtol=0):
        return values[:, 0]
    raise InvalidArtifact('potential atom order/coordinates differ')


def potential_task(task):
    for key in ('state', 'points', 'gbw', 'density', 'source_quality', 'utility'):
        verify(task[key])
    directory = Path(task['directory'])
    if (directory / 'execution.json').exists():
        raise InvalidArtifact('potential attempt exists; collect it explicitly')
    result = run_command([str(verify(task['utility'])), str(verify(task['gbw'])),
                          'endpoint.runtime.scfp', str(verify(task['points'])), str(directory / 'potential.out')],
                         directory, directory / 'utility.log', directory / 'resources.txt')
    result.update(task=task, implementation=record(__file__))
    for key in ('gbw', 'density', 'points'):
        verify(task[key])
    if (directory / 'potential.out').exists():
        result['potential'] = record(directory / 'potential.out')
    write_new(directory / 'execution.json', result)
    return {'task_id': task['task_id'], 'returncode': result['returncode']}


def execute_potentials(manifest):
    if not os.environ.get('SLURM_JOB_ID'):
        raise InvalidArtifact('potential evaluation requires an allocation')
    m = read_json(manifest); verify(m['agreement']); verify(m['implementation'])
    with Path(manifest).with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        workers = min(int(os.environ['SLURM_CPUS_ON_NODE']), len(m['tasks']))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            return list(pool.map(potential_task, m['tasks']))


def collect_potentials(manifest):
    m = read_json(manifest); rows = []
    for task in m['tasks']:
        row = {'task_id': task['task_id'], 'status': 'unavailable'}
        try:
            for key in ('state', 'points', 'gbw', 'density', 'source_quality', 'utility'):
                verify(task[key])
            path = Path(task['directory']) / 'execution.json'; receipt = read_json(path)
            if receipt['task'] != task or receipt['returncode'] != 0 or not receipt['slurm_job_id']:
                raise InvalidArtifact('invalid utility execution receipt')
            state = read_json(verify(task['state']))
            points = np.loadtxt(verify(task['points']), skiprows=1)
            expected = np.array([a['xyz_A'] for a in state['environment_atoms']]) / BOHR_TO_A
            if points.shape != expected.shape or not np.allclose(points, expected, atol=1e-11, rtol=0):
                raise InvalidArtifact('potential probes differ from actual environment positions')
            phi = parse_potential(verify(receipt['potential']), points)
            core = state['core_atoms']; env = state['environment_atoms']
            coords = np.array([a['xyz_A'] for a in core]) / BOHR_TO_A
            distances = np.linalg.norm(points[:, None, :] - coords[None, :, :], axis=2)
            if np.min(distances) * BOHR_TO_A < 1.0:
                raise InvalidArtifact('charge/nucleus overlap violates frozen boundary policy')
            phi_mbis = np.sum(np.array([a['charge_e'] for a in core])[None, :] / distances, axis=1)
            weights = np.array([a['charge_e'] for a in env])
            contributions = weights * phi * HA_TO_KCAL
            approximations = weights * phi_mbis * HA_TO_KCAL
            by_residue = {}
            for atom, value, approx in zip(env, contributions, approximations):
                key = atom['id'].rsplit('/', 1)[0]
                bucket = by_residue.setdefault(key, {'density_kcal_mol': 0., 'mbis_kcal_mol': 0.})
                bucket['density_kcal_mol'] += float(value); bucket['mbis_kcal_mol'] += float(approx)
            actual = float(np.sum(contributions)); approx = float(np.sum(approximations))
            row.update(status='computed', execution_receipt=record(path),
                       direct_density_kcal_mol=actual, direct_mbis_atomic_units_kcal_mol=approx,
                       density_minus_mbis_kcal_mol=actual-approx,
                       direct_mbis_tabi_convention_kcal_mol=approx*TABI_COULOMB_KCAL_A/(BOHR_TO_A*HA_TO_KCAL),
                       minimum_core_environment_distance_A=float(np.min(distances)*BOHR_TO_A),
                       per_residue=by_residue, point_count=len(points),
                       potential_RMS_error_au=float(np.sqrt(np.mean((phi-phi_mbis)**2))))
        except (ValueError, OSError, KeyError) as exc:
            row.update(status='unavailable', reason=str(exc))
        rows.append(row)
    pairs = {}
    for case in CASES:
        pair = {r['task_id'].rsplit('_', 1)[1]: r for r in rows if r['task_id'].startswith(case+'_')}
        if all(pair[metal]['status'] == 'computed' for metal in ('Ca', 'La')):
            pairs[case] = {key: pair['Ca'][key]-pair['La'][key] for key in
                           ('direct_density_kcal_mol', 'direct_mbis_atomic_units_kcal_mol', 'density_minus_mbis_kcal_mol')}
    partition = {key: pairs[CASES[1]][key]-pairs[CASES[0]][key] for key in pairs.get(CASES[0], {})} if len(pairs)==2 else None
    return {'protocol_id': PROTOCOL, 'manifest': record(manifest), 'rows': rows, 'paired': pairs,
            'partition_differences': partition, 'status': 'complete' if len(pairs)==2 else 'incomplete',
            'claim': 'Exact-density permanent coupling diagnostic, not a new global score or affinity prediction.'}


def execute_embedded(manifest):
    m = read_json(manifest)
    for task in m['tasks']:
        verify(task['pointcharges']); verify(task['source_state'])
    result = execute(manifest)
    for task in m['tasks']:
        verify(task['pointcharges'])
    return result


def collect_embedded(manifest, potentials):
    from run_orca_task_manifest import load_manifest_tasks, _completed_attempt_is_valid
    m, normalized = load_manifest_tasks(Path(manifest)); pc = read_json(potentials)
    pm = read_json(verify(pc['manifest']))
    if pm['quantum_manifest'] != record(manifest):
        raise InvalidArtifact('density diagnostic belongs to another endpoint manifest')
    source = read_json(verify(m['source_quantum']))
    vacuum = {r['task_id']: r for r in collect_endpoints(verify(m['source_quantum']))['rows']}
    coupling = {r['task_id']: r for r in pc['rows']}
    rows = []
    for task, normal in zip(m['tasks'], normalized):
        row = {'task_id': task['task_id'], 'status': 'unavailable', 'energy_hartree': None}
        output = Path(task['output_path']); receipt = Path(str(output)+'.execution.json')
        try:
            verify(task['pointcharges']); verify(task['source_state']); verify(task['xyz'])
            expected = ('! r2SCAN-3c NoAutostart DefGrid3 MBIS\n'
                        '%method\n MBIS_LARGEPRINT true\n DoEQ false\nend\n'
                        '%pointcharges "environment.pc"\n'
                        f"* xyzfile {task['charge']} 1 core.xyz\n")
            if verify(task['input']).read_text() != expected:
                raise InvalidArtifact('embedded Hamiltonian declaration changed')
            if not _completed_attempt_is_valid(receipt, output, manifest_sha256=digest(manifest), task=normal,
                    runner_identity=m['execution_policy']['task_runner'],
                    runtime_renderer_identity=m['execution_policy']['runtime_renderer']):
                raise InvalidArtifact('invalid/incomplete embedded execution receipt')
            # Reuse only the native method/ECP/electronic-state checks of the
            # matching vacuum task. New input/PC provenance is checked above.
            electronic = validate_vacuum_output(task['vacuum_task'], output)
            electronic['permanent_external_field'] = True
            text = output.read_text()
            if 'environment.pc' not in text or 'POINT CHARGE' not in text.upper():
                raise InvalidArtifact('missing native external-charge evidence')
            v = vacuum[task['task_id']]; c = coupling[task['task_id']]
            if v['status'] != 'vacuum_endpoint_and_mbis_available' or c['status'] != 'computed':
                raise InvalidArtifact('matching vacuum density/energy unavailable')
            embedded = energy(output)
            response = (embedded-v['energy_hartree'])*HA_TO_KCAL-c['direct_density_kcal_mol']
            row.update(status='computed', energy_hartree=embedded, electronic_state=electronic,
                       output=record(output), execution_receipt=record(receipt),
                       charges=mbis_charges(output, verify(task['xyz']), task['charge']),
                       vacuum_energy_hartree=v['energy_hartree'],
                       direct_density_kcal_mol=c['direct_density_kcal_mol'],
                       density_minus_mbis_kcal_mol=c['density_minus_mbis_kcal_mol'],
                       response_energy_kcal_mol=response,
                       variational_diagnostic='consistent' if response<=0.05 else 'requires_energy_accounting_review')
        except (ValueError, OSError, KeyError) as exc:
            row.update(reason=str(exc))
        rows.append(row)
    pairs = {}
    for case in CASES:
        p = {r['task_id'].rsplit('_', 1)[1]: r for r in rows if r['task_id'].startswith(case+'_')}
        if all(p[metal]['status']=='computed' for metal in ('Ca','La')):
            pairs[case] = {'R_embedded_kcal_mol': (p['Ca']['energy_hartree']-p['La']['energy_hartree'])*HA_TO_KCAL,
                           'R_vacuum_kcal_mol': (p['Ca']['vacuum_energy_hartree']-p['La']['vacuum_energy_hartree'])*HA_TO_KCAL}
            for key in ('direct_density_kcal_mol','density_minus_mbis_kcal_mol','response_energy_kcal_mol'):
                pairs[case][key] = p['Ca'][key]-p['La'][key]
    partition = {key:pairs[CASES[1]][key]-pairs[CASES[0]][key] for key in pairs.get(CASES[0],{})} if len(pairs)==2 else None
    return {'protocol_id': PROTOCOL, 'manifest':record(manifest), 'potentials':record(potentials),
            'rows':rows, 'paired':pairs, 'partition_differences':partition,
            'status':'complete' if len(pairs)==2 else 'incomplete', 'global_score':None,
            'claim':'Permanent-field response diagnostic only; no solvent, calibration or affinity decision.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='operation', required=True)
    p = sub.add_parser('prepare')
    for name in ('states-manifest', 'endpoint-manifest', 'agreement', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    for name in ('execute-potentials', 'collect-potentials', 'execute-embedded', 'dry-run'):
        p = sub.add_parser(name); p.add_argument('--manifest', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('collect-embedded')
    for name in ('manifest','potentials','output'):
        p.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    if args.operation == 'prepare':
        result = prepare(args.states_manifest, args.endpoint_manifest, args.agreement, args.output)
    elif args.operation == 'collect-embedded':
        result = collect_embedded(args.manifest,args.potentials); write_new(args.output,result)
    else:
        fn = {'execute-potentials': execute_potentials, 'collect-potentials': collect_potentials,
              'execute-embedded': execute_embedded, 'dry-run': dry_run}[args.operation]
        result = fn(args.manifest); write_new(args.output, result)
    print(result.get('status', 'recorded') if isinstance(result, dict) else 'recorded')


if __name__ == '__main__':
    main()
