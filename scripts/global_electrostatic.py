"""Approved frozen vacuum-QM / intact-protein electrostatic experiment.

Uses the existing ORCA task runner. No CPCM or isolated-core PB subtraction
belongs to this model. Scientific endpoints and failed attempts are immutable.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import math
import os
from pathlib import Path
import re
import shutil
import time

import numpy as np

from affordable_common import (InvalidArtifact, HA_TO_KCAL, read_json, write_new,
                               record, verify, digest, xyz, paired, cache_key)
from affordable_environment import mbis_charges, physical_boundary_key
from affordable_workflow import dry_run, execute

PROTOCOL = 'vacuum_r2scan3c_mbis_global_tabi_electrostatic_v1'
TABI_COULOMB_KCAL_A = 1389.3875744 / 4.184  # pinned TABI constants.h; old baseline unchanged
CASES = ('1h4i_qm33', '1h4i_qm36')
RADII = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'Ca': 1.8, 'La': 1.8}


def input_state(path, *, cpcm=False):
    text = Path(path).read_text()
    lines = [line.strip() for line in text.splitlines() if line.strip().startswith('!')]
    expected = {'r2scan-3c', 'noautostart', 'defgrid3', 'mbis'}
    if cpcm:
        expected.add('cpcm(water)')
    if len(lines) != 1 or set(lines[0].lower().split()[1:]) != expected:
        raise InvalidArtifact('unexpected endpoint Hamiltonian')
    if re.search(r'%pointcharges|%basis|%pal|\bOpt\b|\bNumGrad\b|\bEnGrad\b', text, re.I):
        raise InvalidArtifact('unsupported endpoint addition')
    match = re.search(r'^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+(\S+)\s*$', text, re.M | re.I)
    if not match or match[2] != '1':
        raise InvalidArtifact('missing singlet coordinate declaration')
    return {'charge': int(match[1]), 'multiplicity': 1, 'xyz_filename': match[3]}


def prepare_partition(source_pilot, environment, output, agreement):
    started = time.monotonic()
    source_pilot, environment, output = (Path(x).resolve() for x in (source_pilot, environment, output))
    if output.exists():
        raise InvalidArtifact('refusing an existing preparation directory')
    m = read_json(source_pilot)
    selected = [t for t in m['tasks'] if t['case'] in CASES]
    if {(t['case'], t['metal']) for t in selected} != {(c, metal) for c in CASES for metal in ('La', 'Ca')}:
        raise InvalidArtifact('partition source task set differs')
    # Verify all scientific inputs before creating a new preparation.
    for t in selected:
        state = input_state(verify(t['input']), cpcm=True)
        verify(t['xyz'])
        if state['charge'] != t['charge']:
            raise InvalidArtifact('source input charge mismatch')
        s = read_json(environment / t['task_id'] / 'skeleton.json')
        verify(s['source']); verify(s['boundary_mapping'])
        if s['source_endpoint_xyz'] != t['xyz'] or s['core_total_charge_e'] != t['charge']:
            raise InvalidArtifact('environment skeleton differs from source endpoint')
    output.mkdir(parents=True)
    tasks = []
    for old in selected:
        d = output / old['case'] / old['metal']; d.mkdir(parents=True)
        xp = d / 'core.xyz'; xp.write_bytes(verify(old['xyz']).read_bytes())
        ip = d / 'endpoint.inp'
        ip.write_text('! r2SCAN-3c NoAutostart DefGrid3 MBIS\n%method\n MBIS_LARGEPRINT true\nend\n'
                      f"* xyzfile {old['charge']} 1 core.xyz\n")
        tasks.append({'task_id': old['task_id'], 'case': old['case'], 'metal': old['metal'],
                      'charge': old['charge'], 'multiplicity': 1,
                      'input': record(ip), 'xyz': record(xp), 'output_path': str(d / 'endpoint.out'),
                      'source_input': old['input'], 'source_xyz': old['xyz'],
                      'baseline_output': record(old['output_path']),
                      'environment_skeleton': record(environment / old['task_id'] / 'skeleton.json')})
    cases = []
    for name in CASES:
        old = next(c for c in m['cases'] if c['case'] == name)
        ts = {t['metal']: t for t in tasks if t['case'] == name}
        invariant = paired(ts['La']['xyz']['path'], ts['Ca']['xyz']['path'], ts['La']['charge'], ts['Ca']['charge'])
        cases.append({**old, 'paired_invariants': invariant, 'use': 'consumed_method_development'})
    implementation = output / 'implementation'; implementation.mkdir()
    deps = [Path(__file__), *Path(__file__).parent.glob('affordable_*.py'),
            Path(__file__).with_name('run_orca_task_manifest.py'), Path(__file__).with_name('render_orca_runtime_input.py')]
    pins = {}
    for p in deps:
        target = implementation / p.name; shutil.copyfile(p, target)
        pins[p.name] = {'source': record(p), 'snapshot': record(target)}
    write_new(implementation / 'inventory.json', pins)
    result = {'schema_version': 'alquemia.global_electrostatic_endpoints.v1', 'protocol_id': PROTOCOL,
              'status': 'approved_prepared', 'stage': 'partition', 'tasks': tasks, 'cases': cases,
              'agreement': record(agreement), 'source_pilot': record(source_pilot), 'orca': m['orca'],
              'execution_policy': {name: record(Path(__file__).with_name(filename)) for name, filename in
                                   [('task_runner', 'run_orca_task_manifest.py'), ('runtime_renderer', 'render_orca_runtime_input.py')]},
              'implementation_snapshot': record(implementation / 'inventory.json'),
              'model': {'quantum': 'native_ORCA_6.1.1_r2SCAN-3c_vacuum_DefGrid3', 'charge_model': 'MBIS',
                        'solvation': 'global_TABI_reaction_energy_only', 'solute_dielectric': 1., 'solvent_dielectric': 78.54,
                        'salt_molar': 0., 'temperature_K': 298.15, 'probe_radius_A': 1.4, 'radii_A': RADII,
                        'energy_expression': 'E_QM_vac + Coulomb(core,environment) + G_RF_full',
                        'reference': None, 'calibrated_threshold': None},
              'cost_tracking': {'compute_budget': None, 'wall_time_limit': None,
                                'preparation_wall_seconds': time.monotonic() - started}}
    write_new(output / 'manifest.json', result)
    return result


def validate_vacuum_output(task, output):
    state = input_state(verify(task['input']))
    if (state['charge'], state['multiplicity']) != (task['charge'], task['multiplicity']):
        raise InvalidArtifact('input/manifest state mismatch')
    text = Path(output).read_text()
    if not re.search(r'Program Version\s+6\.1\.1\b', text):
        raise InvalidArtifact('unexpected ORCA version')
    # ORCA always credits the SMD authors in its banner, including vacuum runs.
    # Reject actual calculation section headings, not those attribution lines.
    if re.search(r'^\s*(?:CPCM SOLVATION MODEL|SMD SOLVATION(?: MODEL)?|COSMO SOLVATION(?: MODEL)?|ORCA NUMERICAL GRADIENT(?: CALCULATION)?)\s*$', text, re.I | re.M):
        raise InvalidArtifact('non-vacuum or unsupported endpoint')
    ecp = re.findall(r'Type\s+(\w+)\s+ECP\s+(\S+)\s+\(replacing\s+(\d+)\s+core electrons', text)
    expected_ecp = [('La', 'Def2-ECP', '46')] if task['metal'] == 'La' else []
    if ecp != expected_ecp:
        raise InvalidArtifact('native ECP convention mismatch')
    numbers = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'S': 16, 'Ca': 20, 'La': 57}
    atoms = xyz(verify(task['xyz']))
    electrons = sum(numbers[a[0]] for a in atoms) - task['charge'] - (46 if task['metal'] == 'La' else 0)
    for pattern, expected in [(r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)', task['charge']),
                              (r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)', 1),
                              (r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)', electrons)]:
        hits = re.findall(pattern, text)
        if len(hits) != 1 or int(hits[0]) != expected:
            raise InvalidArtifact('actual electronic state mismatch')
    if electrons % 2 or not re.search(r'DFTD4', text) or not re.search(r'gCP correction\s+[-+0-9.]', text):
        raise InvalidArtifact('parity/native composite correction mismatch')
    return {'orca_version': '6.1.1', 'explicit_electrons': electrons, 'ecp_core_electrons': 46 if ecp else 0,
            'vacuum': True, 'native_D4_gCP': True}


def collect_endpoints(manifest):
    from affordable_common import energy
    from run_orca_task_manifest import load_manifest_tasks, _completed_attempt_is_valid
    m, normal = load_manifest_tasks(Path(manifest)); rows = []
    for t, n in zip(m['tasks'], normal):
        row = {'task_id': t['task_id'], 'case': t['case'], 'metal': t['metal'],
               'status': 'not_run', 'energy_hartree': None, 'charges': None}
        op = Path(t['output_path']); ep = Path(str(op) + '.execution.json')
        if op.exists():
            try:
                if not _completed_attempt_is_valid(ep, op, manifest_sha256=digest(manifest), task=n,
                        runner_identity=m['execution_policy']['task_runner'],
                        runtime_renderer_identity=m['execution_policy']['runtime_renderer']):
                    raise InvalidArtifact('absent/invalid execution receipt or partial attempt')
                row['electronic_state'] = validate_vacuum_output(t, op)
                row['energy_hartree'] = energy(op)
                row['charges'] = mbis_charges(op, verify(t['xyz']), t['charge'])
                row['status'] = 'vacuum_endpoint_and_mbis_available'
            except (ValueError, OSError) as exc:
                row['status'] = 'unscorable'; row['reason'] = str(exc)
            row['output'] = record(op)
            row['execution_receipt'] = record(ep) if ep.exists() else None
        rows.append(row)
    return {'protocol_id': PROTOCOL, 'manifest': record(manifest), 'rows': rows,
            'status': 'endpoints_complete' if all(r['status'] == 'vacuum_endpoint_and_mbis_available' for r in rows) else 'incomplete',
            'global_score_status': 'awaiting_ESP_and_surface_solver_checks', 'global_score': None}


def endpoint_components(energy_hartree, direct_kcal_mol, reaction_kj_mol):
    if any(v is None for v in (energy_hartree, direct_kcal_mol, reaction_kj_mol)):
        return {'status': 'unavailable', 'total_kcal_mol': None}
    if not all(math.isfinite(v) for v in (energy_hartree, direct_kcal_mol, reaction_kj_mol)):
        raise InvalidArtifact('nonfinite energy component')
    return {'status': 'computed_uncalibrated_descriptor', 'vacuum_energy_hartree': energy_hartree,
            'direct_coulomb_kcal_mol': direct_kcal_mol, 'reaction_field_kJ_mol': reaction_kj_mol,
            'reaction_field_kcal_mol': reaction_kj_mol / 4.184,
            'total_kcal_mol': energy_hartree * HA_TO_KCAL + direct_kcal_mol + reaction_kj_mol / 4.184}


def paired_components(ca, la):
    if any(x.get('status') != 'computed_uncalibrated_descriptor' for x in (ca, la)):
        return {'status': 'unavailable', 'R_global_kcal_mol': None, 'S_global_kcal_mol': None}
    # Subtract in native units first, then convert exactly once.
    quantum = (ca['vacuum_energy_hartree'] - la['vacuum_energy_hartree']) * HA_TO_KCAL
    direct = ca['direct_coulomb_kcal_mol'] - la['direct_coulomb_kcal_mol']
    reaction = (ca['reaction_field_kJ_mol'] - la['reaction_field_kJ_mol']) / 4.184
    return {'status': 'computed_uncalibrated_descriptor', 'R_quantum_kcal_mol': quantum,
            'R_direct_kcal_mol': direct, 'R_reaction_kcal_mol': reaction,
            'R_global_kcal_mol': quantum + direct + reaction, 'S_global_kcal_mol': None,
            'decision': 'uncalibrated_protocol', 'reference_status': 'unavailable'}


def run_esp(manifest, output):
    """Actual ORCA electrostatic-potential checks, with immutable task receipts."""
    from affordable_solver import esp_check
    if not os.environ.get('SLURM_JOB_ID'):
        raise InvalidArtifact('ESP execution requires a compute allocation')
    m = read_json(manifest); c = collect_endpoints(manifest)
    if c['status'] != 'endpoints_complete':
        raise InvalidArtifact('ESP requires four successfully collected vacuum endpoints')
    output = Path(output).resolve(); output.mkdir(parents=True, exist_ok=False)
    write_new(output / 'task_manifest.json', {'protocol_id': PROTOCOL, 'endpoint_manifest': record(manifest),
              'tasks': [t['task_id'] for t in m['tasks']], 'implementation': record(__file__),
              'esp_implementation': record(Path(__file__).with_name('affordable_solver.py')),
              'rule': 'relative_RMS<=0.10 OR absolute_RMS<=0.005_au',
              'compute_budget': None, 'wall_time_limit': None})
    cpus = int(os.environ.get('SLURM_CPUS_ON_NODE', '1')); start = time.monotonic()
    def one(t):
        try:
            q, rec = esp_check(t, verify(m['orca']).parent, RADII, output / t['task_id'])
            return {'task_id': t['task_id'], 'status': q['status'], 'quality': rec}
        except Exception as exc:
            return {'task_id': t['task_id'], 'status': 'unavailable', 'reason': str(exc)}
    with ThreadPoolExecutor(max_workers=min(cpus, len(m['tasks']))) as pool:
        rows = list(pool.map(one, m['tasks']))
    elapsed = time.monotonic() - start
    result = {'protocol_id': PROTOCOL, 'endpoint_manifest': record(manifest), 'rows': rows,
              'status': 'passed' if all(r['status'] == 'passed' for r in rows) else 'failed_or_unavailable',
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'allocated_cpus': cpus,
              'wall_seconds': elapsed, 'allocated_core_seconds': cpus * elapsed}
    write_new(output / 'esp_result.json', result)
    return result


def prepare_states(manifest, esp_result, output):
    """Attach only the new vacuum MBIS charges to the archived physical system."""
    from affordable_state import validate_skeleton_pair
    m = read_json(manifest); ep = read_json(esp_result)
    collection = collect_endpoints(manifest)
    if ep['endpoint_manifest'] != record(manifest) or collection['status'] != 'endpoints_complete':
        raise InvalidArtifact('ESP/endpoint manifest mismatch or incomplete endpoints')
    if ep['status'] != 'passed':
        raise InvalidArtifact('ESP checks failed; frozen-charge model unavailable')
    output = Path(output).resolve()
    if output.exists():
        raise InvalidArtifact('refusing existing state preparation')
    states = {}
    for t in m['tasks']:
        qrec = next(r['quality'] for r in ep['rows'] if r['task_id'] == t['task_id'])
        q = read_json(verify(qrec))
        for key in ('actual_quantum_potential', 'points', 'gbw', 'density', 'utility', 'execution_receipt'):
            verify(q[key])
        if q['status'] != 'passed' or q['charges']['source_output'] != record(t['output_path']) or q['charges']['source_xyz'] != t['xyz']:
            raise InvalidArtifact('charge provenance or quality differs')
        s = read_json(verify(t['environment_skeleton']))
        atoms = xyz(verify(t['xyz']))
        if len(s['core_atoms']) != len(atoms) or len(atoms) != len(q['charges']['charge_e']):
            raise InvalidArtifact('MBIS/source core atom count differs')
        for a, row, charge in zip(s['core_atoms'], atoms, q['charges']['charge_e']):
            if a['element'] != row[0] or a['xyz_A'] != list(row[1:]):
                raise InvalidArtifact('core mapping/coordinates changed')
            a['charge_e'] = charge
        s.update(protocol_id=PROTOCOL, endpoint_input=t['input'], endpoint_output=record(t['output_path']),
                 source_endpoint_xyz=t['xyz'], environment_skeleton=t['environment_skeleton'],
                 charge_quality={'status': 'passed', 'receipt': qrec},
                 settings={k: v for k, v in m['model'].items() if k in
                           ('solute_dielectric', 'solvent_dielectric', 'salt_molar', 'temperature_K', 'probe_radius_A')})
        s['settings']['radii_policy'] = 'Bondi_CHNOS_common_1p8A_metal_v1'
        s['settings']['direct_coulomb_constant_kcal_A_per_mol_e2'] = TABI_COULOMB_KCAL_A
        s['physical_boundary_hash'] = physical_boundary_key(s['physical_atoms'])
        states[t['task_id']] = s
    paired_checks = {name: validate_skeleton_pair(states[name + '_La'], states[name + '_Ca']) for name in CASES}
    if len({s['physical_boundary_hash'] for s in states.values()}) != 1:
        raise InvalidArtifact('numerical partitions do not share the physical protein surface')
    output.mkdir(parents=True)
    pins = {}
    for name, s in states.items():
        p = output / f'{name}.json'; write_new(p, s); pins[name] = record(p)
    result = {'protocol_id': PROTOCOL, 'endpoint_manifest': record(manifest), 'esp_result': record(esp_result),
              'states': pins, 'paired_checks': paired_checks, 'status': 'prepared_for_surface_solver'}
    write_new(output / 'states_manifest.json', result)
    return result


def prepare_surface_campaign(states_manifest, software_manifest, numerics, output):
    from affordable_tabi import prepare
    sm = read_json(states_manifest); output = Path(output).resolve()
    if sm['status'] != 'prepared_for_surface_solver' or set(sm['states']) != {c + '_' + m for c in CASES for m in ('La', 'Ca')}:
        raise InvalidArtifact('unexpected state set or status')
    if output.exists():
        raise InvalidArtifact('refusing existing surface campaign')
    verify(record(numerics)); output.mkdir(parents=True)
    schedule = []
    for name, rec in sm['states'].items():
        for level in ('primary', 'refined', 'tree'):
            schedule.append((name, level, rec, level, 'identity', 'total', 'paired_numerical_checks'))
    for metal in ('La', 'Ca'):
        name = '1h4i_qm33_' + metal
        for transform in ('translated', 'rotated'):
            schedule.append((name, transform, sm['states'][name], 'primary', transform, 'total', 'rigid_transform_check'))
    name = '1h4i_qm33_La'
    schedule.append((name, 'repeat', sm['states'][name], 'primary', 'identity', 'total', 'independent_repeat'))
    for metal in ('La', 'Ca'):
        name = '1h4i_qm33_' + metal
        isolated = copy.deepcopy(read_json(verify(sm['states'][name])))
        isolated.update(physical_atoms=[dict(a, charge_e=0.) for a in isolated['core_atoms']],
                        environment_atoms=[], expected_environment_charge_e=0., cavity_role='isolated_reduction',
                        base_state=sm['states'][name], numerical_control='isolated_core_environment_reduction')
        isolated['physical_boundary_hash'] = physical_boundary_key(isolated['physical_atoms'])
        ip = output / 'isolated_states' / f'{name}.json'; write_new(ip, isolated)
        schedule.append((name, 'isolated', record(ip), 'primary', 'identity', 'total', 'isolated_reduction'))
    for name, rec in sm['states'].items():
        schedule.append((name, 'core_only', rec, 'primary', 'identity', 'core', 'full_cavity_component_audit'))
    for case in CASES:
        name = case + '_La'
        schedule.append((name, 'environment_only', sm['states'][name], 'primary', 'identity', 'environment', 'full_cavity_component_audit'))
    if len(schedule) != 25:
        raise InvalidArtifact('unexpected approved solver count')
    tasks = []
    for name, label, state, level, transform, component, purpose in schedule:
        directory = output / name / label
        prepare(verify(state), directory, software_manifest, level=level, transform=transform,
                component=component, numerics=numerics)
        tasks.append({'task_id': name + '/' + label, 'state_key': name, 'level': level, 'transform': transform,
                      'component': component, 'purpose': purpose, 'manifest': record(directory / 'tabi_manifest.json')})
    result = {'protocol_id': PROTOCOL, 'status': 'prepared_not_executed', 'states_manifest': record(states_manifest),
              'endpoint_manifest': sm['endpoint_manifest'], 'numerics': record(numerics), 'software': record(software_manifest),
              'implementation': record(__file__), 'tasks': tasks,
              'thresholds': {'numerical_kcal_mol': .5, 'partition_kcal_mol': 2., 'reduction_kcal_mol': .01},
              'execution_policy': {'compute_budget': None, 'wall_time_limit': None,
                                   'initial_tasks': ['1h4i_qm33_La/primary', '1h4i_qm33_Ca/primary'],
                                   'initial_reason': 'measure native solver memory/cost before concurrent admission of remaining fixed tasks'}}
    write_new(output / 'campaign_manifest.json', result)
    return result


def execute_surfaces(manifest, selected=None):
    from affordable_tabi import execute as execute_tabi
    from concurrent.futures import as_completed
    import fcntl
    if not os.environ.get('SLURM_JOB_ID'):
        raise InvalidArtifact('surface execution requires a compute allocation')
    mp = Path(manifest).resolve(); m = read_json(mp)
    if m['protocol_id'] != PROTOCOL:
        raise InvalidArtifact('surface campaign protocol differs')
    for key in ('states_manifest', 'endpoint_manifest', 'numerics', 'software'):
        verify(m[key])
    tasks = m['tasks'] if selected is None else [t for t in m['tasks'] if t['task_id'] in selected]
    if not tasks or (selected is not None and {t['task_id'] for t in tasks} != set(selected)):
        raise InvalidArtifact('unknown/empty surface task subset')
    for t in tasks:
        verify(t['manifest'])
    cpus = int(os.environ.get('SLURM_CPUS_ON_NODE', '1'))
    lock = (mp.parent / 'campaign.lock').open('a+')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    job = os.environ['SLURM_JOB_ID']; started = time.monotonic()
    admission = {'campaign': record(mp), 'selected_tasks': [t['task_id'] for t in tasks],
                 'allocated_cpus': cpus, 'workers': min(cpus, len(tasks)), 'slurm_job_id': job,
                 'compute_budget': None, 'wall_time_limit': None, 'implementation': record(__file__)}
    write_new(mp.parent / f'admission_{job}.json', admission)
    rows = []
    def one(task):
        try:
            result = execute_tabi(verify(task['manifest']), reuse=True)
            return {'task_id': task['task_id'], **result}
        except Exception as exc:
            return {'task_id': task['task_id'], 'status': 'unavailable', 'reason': str(exc)}
    with ThreadPoolExecutor(max_workers=min(cpus, len(tasks))) as pool:
        futures = [pool.submit(one, t) for t in tasks]
        for f in as_completed(futures):
            row = f.result(); rows.append(row)
            write_new(mp.parent / f'completion_{job}' / (row['task_id'].replace('/', '__') + '.json'), row)
            print(row['task_id'], row['status'], flush=True)
    elapsed = time.monotonic() - started
    result = {'campaign': record(mp), 'slurm_job_id': job, 'rows': rows,
              'wall_seconds': elapsed, 'allocated_cpus': cpus, 'allocated_core_seconds': cpus * elapsed,
              'status': 'selected_tasks_finished'}
    write_new(mp.parent / f'execution_{job}.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='op', required=True)
    q = sub.add_parser('prepare-partition')
    for flag in ('source-pilot', 'environment', 'output', 'agreement'):
        q.add_argument('--' + flag, type=Path, required=True)
    for name in ('dry-run', 'execute', 'collect'):
        q = sub.add_parser(name); q.add_argument('--manifest', type=Path, required=True)
        q.add_argument('--output', type=Path)
    q = sub.add_parser('esp'); q.add_argument('--manifest', type=Path, required=True); q.add_argument('--output', type=Path, required=True)
    q = sub.add_parser('prepare-states')
    for name in ('manifest', 'esp-result', 'output'):
        q.add_argument('--' + name, type=Path, required=True)
    q = sub.add_parser('prepare-surfaces')
    for name in ('states-manifest', 'software-manifest', 'numerics', 'output'):
        q.add_argument('--' + name, type=Path, required=True)
    q = sub.add_parser('execute-surfaces'); q.add_argument('--manifest', type=Path, required=True)
    q.add_argument('--tasks', nargs='+')
    a = p.parse_args()
    if a.op == 'prepare-partition':
        r = prepare_partition(a.source_pilot, a.environment, a.output, a.agreement)
    elif a.op == 'esp':
        r = run_esp(a.manifest, a.output)
    elif a.op == 'prepare-states':
        r = prepare_states(a.manifest, a.esp_result, a.output)
    elif a.op == 'prepare-surfaces':
        r = prepare_surface_campaign(a.states_manifest, a.software_manifest, a.numerics, a.output)
    elif a.op == 'execute-surfaces':
        r = execute_surfaces(a.manifest, a.tasks)
    else:
        r = {'dry-run': dry_run, 'execute': execute, 'collect': collect_endpoints}[a.op](a.manifest)
        if a.output:
            write_new(a.output, r)
    print(r.get('status', 'prepared'))


if __name__ == '__main__':
    main()
