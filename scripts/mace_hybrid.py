"""Opt-in, manifested MACE-POLAR capability pilot; never changes baseline scores."""
from __future__ import annotations

import argparse
import contextlib
import copy
import fcntl
import importlib.metadata
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import traceback

import numpy as np

from affordable_common import (HA_TO_KCAL, InvalidArtifact, cache_key, digest,
                               energy, read_json, record, verify, write_new, xyz)

PROTOCOL = 'mace_polar_1m_vacuum_r2scan3c_subtractive_pilot_v1'
EV_TO_KCAL = 23.06054783061903  # eV * exact SI elementary charge * exact Avogadro / 4184
TASK_IDS = [f'1h4i_qm{size}_{metal}' for size in (33, 36) for metal in ('La', 'Ca')]
Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'Na': 11, 'S': 16, 'Ca': 20, 'La': 57}


def rotation():
    axis = np.array([1., 2., 3.]); axis /= np.linalg.norm(axis)
    x, y, z = axis
    cross = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    angle = np.deg2rad(37.)
    return np.eye(3) * np.cos(angle) + (1 - np.cos(angle)) * np.outer(axis, axis) + np.sin(angle) * cross


def write_xyz(path, atoms):
    with Path(path).open('x') as f:
        f.write(f'{len(atoms)}\nPinned real geometry; identity recorded in the MACE manifest\n')
        for a in atoms:
            f.write(a[0] + ' ' + ' '.join(format(v, '.17g') for v in a[1:]) + '\n')


def check_atoms(atoms, charge, *, background_calcium_indices=()):
    if any(a[0] not in Z for a in atoms):
        raise InvalidArtifact('unsupported element')
    background = set(background_calcium_indices)
    if (len(background) != len(background_calcium_indices)
            or any(not isinstance(i, int) or not 0 <= i < len(atoms) or atoms[i][0] != 'Ca' for i in background)):
        raise InvalidArtifact('invalid explicitly mapped background calcium')
    if sum(a[0] in ('Ca', 'La') for i,a in enumerate(atoms) if i not in background) != 1:
        raise InvalidArtifact('exactly one selected metal required')
    electrons = sum(Z[a[0]] for a in atoms) - charge
    if electrons % 2:
        raise InvalidArtifact('singlet electron parity mismatch')
    if not np.isfinite(np.array([a[1:] for a in atoms])).all():
        raise InvalidArtifact('nonfinite coordinates')
    return {'atoms': len(atoms), 'all_electron_count': electrons, 'multiplicity': 1}


def prepare(root, software, output, agreement):
    from global_electrostatic import collect_endpoints
    root, output = Path(root).resolve(), Path(output).resolve()
    sm = read_json(software)
    for key in ('checkpoint', 'python', 'requirements', 'backend_source_inventory'):
        verify(sm[key])
    source = root / 'workspaces/global_electrostatic_20260916'
    em = source / 'partition_tasks_v1/manifest.json'
    collection = collect_endpoints(em)
    if collection['status'] != 'endpoints_complete':
        raise InvalidArtifact('archived vacuum endpoints are not valid')
    state_manifest = read_json(source / 'states_v1/states_manifest.json')
    states = {name: read_json(verify(state_manifest['states'][name])) for name in TASK_IDS}
    physical = None; mapping = {}; sources = []; coordinate_roundoff = 0.
    for name, s in states.items():
        for key in ('source', 'endpoint_input', 'endpoint_output', 'source_endpoint_xyz', 'boundary_mapping'):
            verify(s[key])
        current = [(a['id'], a['element'], a['xyz_A']) for a in s['physical_atoms']]
        if physical is None:
            physical = current
        else:
            if [(a[0], a[1]) for a in current] != [(a[0], a[1]) for a in physical]:
                raise InvalidArtifact('physical atom identity changes across endpoints/partitions')
            error = float(np.max(np.abs(np.array([a[2] for a in current])-np.array([a[2] for a in physical]))))
            coordinate_roundoff = max(coordinate_roundoff, error)
            if error > 1e-12:
                raise InvalidArtifact('physical coordinates change across endpoints/partitions')
        b = read_json(verify(s['boundary_mapping']))
        if abs(b['ff_chain_charge_e'] + 8) > 1e-6 or s['explicit_waters']:
            raise InvalidArtifact('unapproved chain charge/water inventory')
        expected = -8 - 3 + (3 if s['metal'] == 'La' else 2)
        if abs(s['core_total_charge_e'] + s['expected_environment_charge_e'] - expected) > 1e-6:
            raise InvalidArtifact('formal full-system charge closure failed')
        if s['assembly'] != 'deposited_catalytic_chain_A':
            raise InvalidArtifact('assembly mismatch')
        lookup = {a['id']: i for i, a in enumerate(s['physical_atoms'])}
        if len(lookup) != len(s['physical_atoms']):
            raise InvalidArtifact('duplicate physical source IDs')
        mapping[name] = [{'core_index': i, 'id': a['id'], 'kind': a.get('kind'),
                          'physical_index': lookup.get(a['id'])}
                         for i, a in enumerate(s['core_atoms'])]
        for a, m in zip(s['core_atoms'], mapping[name]):
            if a.get('kind') == 'cap':
                if m['physical_index'] is not None:
                    raise InvalidArtifact('cap present in full protein')
            elif m['physical_index'] is None or a['xyz_A'] != s['physical_atoms'][m['physical_index']]['xyz_A']:
                raise InvalidArtifact('core source-coordinate mapping mismatch')
        sources.append(state_manifest['states'][name])
    output.mkdir(parents=True, exist_ok=False)
    write_new(output / 'archived_endpoints.json', collection)
    write_new(output / 'atom_mappings.json', mapping)
    impl = output / 'implementation'; impl.mkdir()
    pins = {}
    for name in ('mace_hybrid.py', 'affordable_common.py', 'mace_realspace_compat.py', 'mace_blocked.py', 'mace_local_memory.py'):
        p = Path(__file__).with_name(name); q = impl / name
        shutil.copyfile(p, q); pins[name] = record(q)
    tasks = []
    for name in TASK_IDS:
        s = states[name]; atoms = xyz(verify(s['source_endpoint_xyz']))
        p = output / f'{name}.xyz'; write_xyz(p, atoms)
        tasks.append({'task_id': name, 'system': name.rsplit('_', 1)[0], 'kind': 'core',
                      'variant': 'primary', 'metal': s['metal'], 'charge': s['core_total_charge_e'],
                      'spin_multiplicity': 1, 'xyz': record(p), 'state': check_atoms(atoms, s['core_total_charge_e']),
                      'rotation_matrix': np.eye(3).tolist()})
    center = np.array(next(a[2] for a in physical if a[0] == 'metal'))
    for variant in ('primary', 'repeat', 'translate', 'rotate'):
        for metal in ('La', 'Ca'):
            matrix = rotation() if variant == 'rotate' else np.eye(3)
            shift = np.array([10., -7., 3.]) if variant == 'translate' else np.zeros(3)
            atoms = []
            for identity, elem, coords in physical:
                pos = (np.array(coords) - center) @ matrix.T + center + shift if variant in ('rotate', 'translate') else np.array(coords)
                atoms.append((metal if elem == 'M' else elem, *pos.tolist()))
            charge = -8 if metal == 'La' else -9
            name = f'1h4i_full_{metal}_{variant}'; p = output / f'{name}.xyz'; write_xyz(p, atoms)
            tasks.append({'task_id': name, 'system': '1h4i_full', 'kind': 'full', 'variant': variant,
                          'metal': metal, 'charge': charge, 'spin_multiplicity': 1,
                          'xyz': record(p), 'state': check_atoms(atoms, charge),
                          'rotation_matrix': matrix.tolist()})
    model = {'checkpoint': sm['checkpoint'], 'model_type': 'PolarMACE', 'dtype': 'float64',
             'pbc_handling': 'realspace', 'periodic': False, 'external_field': [0., 0., 0.],
             'solvent': None, 'assembly': 'deposited_catalytic_chain_A', 'microstate': states[TASK_IDS[0]]['microstate'],
             'explicit_waters': [], 'spin_input_convention': 'multiplicity_1_closed_shell',
             'backend_source_inventory': sm['backend_source_inventory']}
    from mace_realspace_compat import ADAPTER_ID
    model['execution_adapter'] = ADAPTER_ID
    for t in tasks:
        t['cache_key'] = cache_key({'task': t, 'model': model, 'software': record(software), 'implementation': pins})
    result = {'schema_version': 'alquemia.mace_pilot.v1', 'protocol_id': PROTOCOL,
              'status': 'approved_prepared', 'agreement': record(agreement), 'software': record(software),
              'source_states': sources, 'source_endpoint_manifest': record(em), 'archived_endpoints': record(output/'archived_endpoints.json'),
              'atom_mappings': record(output/'atom_mappings.json'), 'implementation': pins,
              'archived_physical_coordinate_roundoff_max_A': coordinate_roundoff,
              'model': model, 'tasks': tasks, 'reference': None, 'calibrated_decision': None,
              'evidence_use': 'consumed_method_development',
              'tolerances': {'energy_kcal_mol': .01, 'force_max_eV_A': .001, 'total_charge_e': 1e-5, 'partition_kcal_mol': 2.},
              'run_inventory': {'distinct_mace_energy_force_calls': 12, 'new_DFT_endpoints': 0},
              'compute_budget': None, 'wall_time_limit': None}
    write_new(output/'manifest.json', result)
    return {'status': 'prepared', 'manifest': record(output/'manifest.json'), 'tasks': len(tasks)}


def dry_run(manifest):
    if read_json(manifest).get('schema_version') == 'alquemia.mace_omol.v1':
        from mace_omol import validate
        return validate(manifest)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_canonical.v1', 'alquemia.mace_canonical_gb.v1'):
        from mace_canonical_run import validate
        return validate(manifest)
    if read_json(manifest).get('schema_version') == 'alquemia.mace_mechanics_minimum_short.v1':
        from mace_mechanics_minimum import validate
        return validate(manifest)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_mechanics_core.v1', 'alquemia.mace_mechanics_gb.v1', 'alquemia.mace_mechanics_short.v1'):
        from mace_mechanics_run import validate
        return validate(manifest)
    if read_json(manifest).get('schema_version') == 'alquemia.mace_short_engine.v1':
        from mace_short_engine import validate
        return validate(manifest)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_curvature.v1', 'alquemia.mace_curvature_gb.v1'):
        from mace_curvature import validate
        return validate(manifest)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_local_correction.v1', 'alquemia.mace_local_gb.v1'):
        from mace_local_correction import validate
        return validate(manifest)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_global_benchmark.v1', 'alquemia.mace_global_gb.v1'):
        from mace_global_benchmark import validate
        return validate(manifest)
    if read_json(manifest).get('schema_version') == 'alquemia.mace_gb.v1':
        from mace_gb import validate
        return validate(manifest)
    m = read_json(manifest); sm = read_json(verify(m['software']))
    for ref in [m['agreement'], m['archived_endpoints'], m['atom_mappings'], *m['source_states'],
                *m['implementation'].values(), sm['checkpoint'], sm['python'], sm['requirements'], sm['backend_source_inventory']]:
        verify(ref)
    inventory = read_json(verify(sm['backend_source_inventory']))
    for ref in inventory['files']:
        verify(ref)
    for key in ('memory_agreement', 'kernel_validation_reference'):
        if m.get(key):
            verify(m[key])
    if m.get('schema_version') == 'alquemia.mace_hydrogen.v1':
        from mace_hydrogen import validate
        validate(m)
    elif m.get('schema_version') == 'alquemia.mace_response_trace.v1':
        from mace_response_trace import validate
        validate(m)
    elif m.get('schema_version') == 'alquemia.mace_analytic.v1':
        from mace_analytic_pilot import validate
        validate(m)
    elif m.get('schema_version') == 'alquemia.mace_rotation.v1':
        from mace_rotation import validate
        validate(m)
    elif len(m['tasks']) != 12 or len({t['task_id'] for t in m['tasks']}) != 12:
        raise InvalidArtifact('approved task set changed')
    for t in m['tasks']:
        check_atoms(xyz(verify(t['xyz'])), t['charge'])
        base = {k: v for k, v in t.items() if k != 'cache_key'}
        expected = cache_key({'task': base, 'model': m['model'], 'software': m['software'], 'implementation': m['implementation']})
        if t['cache_key'] != expected:
            raise InvalidArtifact('task scientific cache key mismatch')
    return {'status': 'pass', 'tasks': len(m['tasks']), 'new_DFT_endpoints': 0, 'manifest': record(manifest)}


def repair_interface(manifest, output, pair_tile=None, memory_agreement=None, reference_collection=None, edge_tile=None, node_tile=None):
    """Version only the execution adapter; reuse every exact pinned input byte."""
    from mace_realspace_compat import ADAPTER_ID
    dry_run(manifest)
    m = copy.deepcopy(read_json(manifest))
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    impl = out/'implementation'; impl.mkdir()
    pins = {}
    for name in ('mace_hybrid.py', 'affordable_common.py', 'mace_realspace_compat.py', 'mace_blocked.py', 'mace_local_memory.py'):
        destination = impl/name
        shutil.copyfile(Path(__file__).with_name(name), destination)
        pins[name] = record(destination)
    m['implementation'] = pins
    m['model']['execution_adapter'] = ADAPTER_ID
    m['technical_revision_of'] = record(manifest)
    if pair_tile is not None:
        if pair_tile < 1 or not memory_agreement or not reference_collection:
            raise InvalidArtifact('blocked revision requires positive tile, agreement and actual reference collection')
        m['model']['pair_kernel'] = {'kernel_id': 'graph044_blocked_analytic_backward_v1', 'tile_sites': pair_tile}
        if edge_tile is not None:
            if edge_tile < 1: raise InvalidArtifact('positive edge tile required')
            m['model']['pair_kernel']['edge_tile'] = edge_tile
            m['model']['pair_kernel']['core_edge_tile'] = min(edge_tile,128)
        if node_tile is not None:
            if node_tile < 1 or edge_tile is None: raise InvalidArtifact('node tile requires edge checkpointing')
            m['model']['pair_kernel']['node_tile'] = node_tile
            m['model']['pair_kernel']['core_node_tile'] = min(node_tile,17)
        m['memory_agreement'] = record(memory_agreement)
        m['kernel_validation_reference'] = record(reference_collection)
        m['kernel_validation_tolerances'] = {'energy_eV': 1e-6, 'force_eV_A': 1e-6, 'density': 1e-8}
    for task in m['tasks']:
        task.pop('cache_key')
        task['cache_key'] = cache_key({'task': task, 'model': m['model'],
                                      'software': m['software'], 'implementation': pins})
    write_new(out/'manifest.json', m)
    return dry_run(out/'manifest.json')


def worker(manifest, task_id, output, memory_mode):
    if read_json(manifest).get('schema_version') == 'alquemia.mace_omol.v1':
        from mace_omol import worker as omol_worker
        return omol_worker(manifest, task_id, output, memory_mode)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_gb.v1', 'alquemia.mace_global_gb.v1', 'alquemia.mace_local_gb.v1', 'alquemia.mace_curvature_gb.v1', 'alquemia.mace_mechanics_gb.v1', 'alquemia.mace_canonical_gb.v1'):
        from mace_gb import worker as gb_worker
        return gb_worker(manifest, task_id, output, memory_mode)
    import torch
    from ase import Atoms
    from mace.calculators import mace_polar
    m = read_json(manifest); t = next(t for t in m['tasks'] if t['task_id'] == task_id)
    output = Path(output); result = {'task_id': task_id, 'cache_key': t['cache_key'], 'memory_mode': memory_mode,
                                    'manifest': record(manifest), 'status': 'unavailable'}
    start = time.monotonic()
    try:
        if not os.environ.get('SLURM_JOB_ID') or not torch.cuda.is_available():
            raise InvalidArtifact('allocated CUDA device required')
        torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK']))
        torch.set_default_dtype(torch.float64)
        torch.cuda.reset_peak_memory_stats()
        props = torch.cuda.get_device_properties(0)
        result['device'] = {'name': props.name, 'total_memory_bytes': props.total_memory,
                            'CUDA_VISIBLE_DEVICES': os.environ.get('CUDA_VISIBLE_DEVICES')}
        print(json.dumps({'event': 'model_load_start', 'task': task_id, 'device': result['device']}), flush=True)
        calc = mace_polar(model=str(verify(m['model']['checkpoint'])), device='cuda',
                          default_dtype='float64', pbc_handling='realspace')
        if m['model'].get('execution_adapter'):
            from mace_realspace_compat import ADAPTER_ID, configure
            if m['model']['execution_adapter'] != ADAPTER_ID:
                raise InvalidArtifact('unsupported execution adapter')
            result['execution_adapter'] = configure(calc)
        if m['model'].get('pair_kernel'):
            settings = m['model']['pair_kernel']
            edge_tile = settings.get('core_edge_tile', settings.get('edge_tile')) if t['kind'] == 'core' else settings.get('edge_tile')
            node_tile = settings.get('core_node_tile', settings.get('node_tile')) if t['kind'] == 'core' else settings.get('node_tile')
            if settings['kernel_id'] == 'graph044_analytic_multipoles_blocked_v1':
                from mace_analytic import configure as configure_analytic
                result['pair_kernel'] = configure_analytic(calc, settings['tile_atoms'], edge_tile, node_tile)
            else:
                from mace_blocked import KERNEL_ID, configure as configure_blocked
                if settings['kernel_id'] != KERNEL_ID:
                    raise InvalidArtifact('unsupported blocked pair kernel')
                result['pair_kernel'] = configure_blocked(calc, settings['tile_sites'], edge_tile, node_tile)
        result['model_load_seconds'] = time.monotonic() - start
        charge_observer = None
        if m['model'].get('response_trace'):
            from mace_response_trace import ChargeTrace, TRACE_ID
            if m['model']['response_trace'] != TRACE_ID:
                raise InvalidArtifact('unsupported response observer')
            charge_observer = ChargeTrace(calc.models[0])
        rows = xyz(verify(t['xyz']))
        atoms = Atoms([a[0] for a in rows], positions=[a[1:] for a in rows], pbc=False)
        atoms.info.update(charge=t['charge'], spin=t['spin_multiplicity'], external_field=[0., 0., 0.])
        atoms.calc = calc
        print(json.dumps({'event': 'evaluation_start', 'task': task_id, 'atoms': len(atoms), 'charge': t['charge'],
                          'memory_mode': memory_mode}), flush=True)
        torch.cuda.synchronize(); evaluation_start = time.monotonic()
        cm = torch.autograd.graph.save_on_cpu(pin_memory=True) if memory_mode == 'host_offload' else contextlib.nullcontext()
        with cm:
            if m['model'].get('short_adapter'):
                from mace_short_engine import ADAPTER, evaluate
                if m['model']['short_adapter'] != ADAPTER or t.get('energy_component') != 'interaction_energy':
                    raise InvalidArtifact('unsupported short component adapter')
                value, forces = evaluate(calc, atoms)
            else:
                value = float(atoms.get_potential_energy())
                forces = np.asarray(atoms.get_forces(), dtype=np.float64)
        torch.cuda.synchronize()
        result['evaluation_seconds'] = time.monotonic() - evaluation_start
        if not math.isfinite(value) or forces.shape != (len(atoms), 3) or not np.isfinite(forces).all():
            raise InvalidArtifact('invalid MACE energy/forces')
        fp = output / 'forces_eV_A.npy'; np.save(fp, forces)
        if m['model'].get('short_adapter'):
            result.update(status='computed', energy_eV=value, forces=record(fp),
                          energy_component='interaction_energy', density_coefficients=None,
                          charge_check=None, charge_response_status='not_part_of_this_component',
                          energy_definition='trained_local_interaction_energy_only',
                          force_definition='negative_Cartesian_gradient_of_interaction_energy',
                          input_state_check=check_atoms(rows,t['charge']))
        else:
            density = np.asarray(calc.results['density_coefficients'])
            if density.shape != (len(atoms), 4) or not np.isfinite(density).all():
                raise InvalidArtifact(f'unexpected density shape {density.shape}')
            dp = output / 'density_coefficients.npy'; np.save(dp, density)
            total_charge = float(density[:, 0].sum())
            result.update(status='computed', energy_eV=value, forces=record(fp), density_coefficients=record(dp),
                          density_total_charge_e=total_charge, total_charge_error_e=total_charge-t['charge'],
                          charge_check=abs(total_charge-t['charge']) <= m['tolerances']['total_charge_e'])
            if charge_observer is not None:
                result['charge_trace'] = charge_observer.finish(output, t['charge'], density)
        if m.get('schema_version') == 'alquemia.mace_rotation.v1':
            from mace_rotation import probe
            result['rotation_probe'] = probe(calc, m, t, output)
        if m.get('schema_version') in ('alquemia.mace_rotation.v1', 'alquemia.mace_analytic.v1', 'alquemia.mace_response_trace.v1', 'alquemia.mace_hydrogen.v1', 'alquemia.mace_global_benchmark.v1', 'alquemia.mace_local_correction.v1', 'alquemia.mace_curvature.v1', 'alquemia.mace_mechanics_core.v1', 'alquemia.mace_canonical.v1'):
            result['energy_components_eV'] = {key: float(calc.results[key]) for key in
                ('interaction_energy', 'electrostatic_energy', 'electron_energy')}
    except Exception as exc:
        result.update(status='failed', reason=str(exc), exception_type=type(exc).__name__)
        traceback.print_exc()
    finally:
        result.update(wall_seconds=time.monotonic()-start, peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_available() else None,
                      peak_cuda_reserved_bytes=torch.cuda.max_memory_reserved() if torch.cuda.is_available() else None,
                      slurm_job_id=os.environ.get('SLURM_JOB_ID'), allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),
                      allocated_host_mem_MiB=os.environ.get('SLURM_MEM_PER_NODE'),
                      versions={name: importlib.metadata.version(name) for name in ('torch', 'mace-torch', 'ase', 'e3nn')})
        write_new(output/'result.json', result)
    print(json.dumps({k: result.get(k) for k in ('task_id', 'status', 'reason', 'evaluation_seconds', 'peak_host_RSS_KiB', 'peak_cuda_allocated_bytes')}), flush=True)
    return result


def accepted_attempt(attempt, task, manifest):
    try:
        r = read_json(attempt/'receipt.json'); result = read_json(verify(r['result']))
        if (r['returncode'] != 0 or r['manifest'] != record(manifest) or r['task_id'] != task['task_id']
                or result['status'] != 'computed' or result['cache_key'] != task['cache_key']):
            return None
        omol_components=('MACE_OMOL_total_vacuum_energy','MACE_OMOL_charge_feature_ablated_descriptor')
        if not (task.get('energy_only') is True and task.get('energy_component') in omol_components):
            verify(result['forces'])
        if task.get('energy_component') in omol_components:
            from mace_omol import accepted_state
            if not accepted_state(result, task):
                return None
        elif task.get('energy_component') == 'interaction_energy':
            if (result.get('energy_component') != 'interaction_energy' or result.get('density_coefficients') is not None
                    or not result.get('input_state_check') or not math.isfinite(result['energy_eV'])):
                return None
        elif not result.get('charge_check'):
            return None
        else:
            verify(result['density_coefficients'])
        if result.get('charge_trace'):
            trace = read_json(verify(result['charge_trace']))
            verify(trace['arrays'])
        return result
    except (OSError, ValueError, KeyError, TypeError):
        return None


def execute(manifest, memory_mode, selected=None):
    if not os.environ.get('SLURM_JOB_ID'):
        raise InvalidArtifact('execute requires an allocation')
    dry_run(manifest)
    mp = Path(manifest).resolve(); m = read_json(mp); sm = read_json(verify(m['software']))
    root = mp.parent/'execution'; root.mkdir(exist_ok=True)
    with (root/'execute.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        selected = set(selected or [t['task_id'] for t in m['tasks']])
        if selected - {t['task_id'] for t in m['tasks']}:
            raise InvalidArtifact('unknown task selection')
        passed_ablation_gates=set()
        for t in m['tasks']:
            if t['task_id'] not in selected:
                continue
            if m.get('stage')=='ablation_development' and t['kind']=='full':
                from mace_omol_ablation_run import execution_gate
                gate='core' if t['ablation_group']=='numerical' else 'numerical'
                if gate not in passed_ablation_gates:
                    if not execution_gate(mp,gate):
                        raise InvalidArtifact('descriptor execution blocked by incomplete/failed '+gate+' gate')
                    passed_ablation_gates.add(gate)
            if t['kind'] == 'full' and m.get('schema_version') == 'alquemia.mace_gb.v1':
                from mace_gb import core_gate
                if core_gate(mp)['status'] != 'pass':
                    raise InvalidArtifact('full GB execution requires native/custom core agreement')
            if t['kind'] == 'full' and m['model'].get('pair_kernel'):
                if m.get('schema_version') in ('alquemia.mace_global_benchmark.v1', 'alquemia.mace_short_engine.v1', 'alquemia.mace_mechanics_short.v1', 'alquemia.mace_mechanics_minimum_short.v1'):
                    from mace_global_benchmark import numerical_parent_gate
                    validation = numerical_parent_gate(m)
                elif m.get('schema_version') == 'alquemia.mace_hydrogen.v1':
                    from mace_hydrogen import validate
                    validation = validate(m)
                elif m.get('schema_version') == 'alquemia.mace_response_trace.v1':
                    from mace_response_trace import core_gate
                    validation = core_gate(mp)
                elif m.get('schema_version') == 'alquemia.mace_analytic.v1':
                    from mace_analytic_pilot import core_gate
                    validation = core_gate(mp)
                else:
                    validation = compare_cores(mp)
                if validation['status'] != 'pass':
                    raise InvalidArtifact('full-system execution requires its declared core validation gate')
            td = root/t['task_id']; td.mkdir(exist_ok=True)
            previous = sorted(td.glob('attempt_*'))
            if any(accepted_attempt(a, t, mp) is not None for a in previous):
                print(json.dumps({'task_id': t['task_id'], 'status': 'verified_cache_reuse'}), flush=True)
                continue
            attempt = td / f'attempt_{len(previous)+1:03d}'; attempt.mkdir()
            command = [str(verify(sm['python'])), str(verify(m['implementation']['mace_hybrid.py'])),
                       'worker', '--manifest', str(mp), '--task-id', t['task_id'], '--output', str(attempt),
                       '--memory-mode', memory_mode]
            write_new(attempt/'started.json', {'manifest': record(mp), 'task_id': t['task_id'], 'command': command,
                      'memory_mode': memory_mode, 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'started_unix': time.time()})
            start = time.monotonic()
            with (attempt/'worker.log').open('x') as log:
                p = subprocess.run(['/usr/bin/time', '-v', '-o', str(attempt/'resources.txt'), *command],
                                   stdout=log, stderr=subprocess.STDOUT, check=False)
            receipt = {'manifest': record(mp), 'task_id': t['task_id'], 'returncode': p.returncode,
                       'wall_seconds': time.monotonic()-start, 'memory_mode': memory_mode,
                       'allocated_cpus': int(os.environ['SLURM_CPUS_PER_TASK']), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
                       'log': record(attempt/'worker.log'), 'resource_usage': record(attempt/'resources.txt'),
                       'result': record(attempt/'result.json') if (attempt/'result.json').exists() else None}
            receipt['allocated_core_seconds'] = receipt['allocated_cpus'] * receipt['wall_seconds']
            write_new(attempt/'receipt.json', receipt)
            valid = accepted_attempt(attempt, t, mp)
            print(json.dumps({'task_id': t['task_id'], 'status': 'computed' if valid else 'failed',
                              'wall_seconds': receipt['wall_seconds']}), flush=True)
            if valid is None:
                return {'status': 'partial_failure', 'failed_task': t['task_id'], 'attempt': str(attempt)}
    return {'status': 'selected_tasks_complete'}


def compare_cores(manifest):
    m = read_json(manifest)
    reference = read_json(verify(m['kernel_validation_reference']))
    reference_m = read_json(verify(reference['manifest']))
    if reference_m['model']['checkpoint'] != m['model']['checkpoint']:
        raise InvalidArtifact('core-validation checkpoint mismatch')
    current = collect(manifest)
    tolerances = m['kernel_validation_tolerances']
    rows = []
    for name in TASK_IDS:
        task = next(t for t in m['tasks'] if t['task_id'] == name)
        old_task = next(t for t in reference_m['tasks'] if t['task_id'] == name)
        for field in ('xyz', 'charge', 'spin_multiplicity'):
            if task[field] != old_task[field]:
                raise InvalidArtifact('core-validation input mismatch')
        old, new = reference['rows'][name], current['rows'][name]
        if old['status'] != 'computed' or new['status'] != 'computed':
            rows.append({'task_id': name, 'status': 'unavailable', 'pass': False})
            continue
        de = new['energy_eV'] - old['energy_eV']
        df = float(np.max(np.abs(np.load(verify(new['forces']))-np.load(verify(old['forces'])))))
        dd = float(np.max(np.abs(np.load(verify(new['density_coefficients']))-np.load(verify(old['density_coefficients'])))))
        rows.append({'task_id': name, 'energy_delta_eV': de, 'force_max_delta_eV_A': df,
                     'density_max_delta': dd, 'pass': abs(de) <= tolerances['energy_eV'] and
                     df <= tolerances['force_eV_A'] and dd <= tolerances['density']})
    return {'status': 'pass' if all(r['pass'] for r in rows) else 'failed',
            'manifest': record(manifest), 'reference': m['kernel_validation_reference'],
            'tolerances': tolerances, 'rows': rows}


def collect(manifest):
    if read_json(manifest).get('schema_version') == 'alquemia.mace_omol.v1':
        from mace_omol import collect as collect_omol
        return collect_omol(manifest)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_canonical.v1', 'alquemia.mace_canonical_gb.v1'):
        from mace_canonical_run import collect as collect_canonical
        return collect_canonical(manifest)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_mechanics_core.v1', 'alquemia.mace_mechanics_gb.v1', 'alquemia.mace_mechanics_short.v1', 'alquemia.mace_mechanics_minimum_short.v1'):
        from mace_mechanics_run import collect_mechanics
        return collect_mechanics(manifest)
    if read_json(manifest).get('schema_version') == 'alquemia.mace_short_engine.v1':
        from mace_short_engine import collect_short
        return collect_short(manifest)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_curvature.v1', 'alquemia.mace_curvature_gb.v1'):
        from mace_curvature import collect_curvature
        return collect_curvature(manifest)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_local_correction.v1', 'alquemia.mace_local_gb.v1'):
        from mace_local_correction import collect_local
        return collect_local(manifest)
    if read_json(manifest).get('schema_version') in ('alquemia.mace_global_benchmark.v1', 'alquemia.mace_global_gb.v1'):
        from mace_global_benchmark import collect_global
        return collect_global(manifest)
    if read_json(manifest).get('schema_version') == 'alquemia.mace_gb.v1':
        from mace_gb import collect_gb
        return collect_gb(manifest)
    dry_run(manifest); mp = Path(manifest).resolve(); m = read_json(mp)
    if m.get('schema_version') == 'alquemia.mace_hydrogen.v1':
        from mace_hydrogen import collect_hydrogen
        return collect_hydrogen(mp)
    if m.get('schema_version') == 'alquemia.mace_rotation.v1':
        from mace_rotation import collect_rotation
        return collect_rotation(mp)
    rows = {}; attempts = []
    for t in m['tasks']:
        valid = []
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            result = accepted_attempt(a, t, mp)
            attempts.append({'task_id': t['task_id'], 'path': str(a), 'accepted': result is not None,
                             'receipt': record(a/'receipt.json') if (a/'receipt.json').exists() else None})
            if result is not None:
                valid.append(result)
        rows[t['task_id']] = valid[-1] if valid else {'status': 'unavailable', 'energy_eV': None}
    complete = all(r['status'] == 'computed' for r in rows.values())
    result = {'protocol_id': m['protocol_id'], 'manifest': record(mp), 'collection_implementation': record(__file__),
              'status': 'complete' if complete else 'incomplete',
              'rows': rows, 'attempts': attempts, 'reference': None, 'S_kcal_mol': None,
              'decision': 'uncalibrated_vacuum_diagnostic', 'direct': None, 'hybrid': None, 'checks': None,
              'core_partition': None}
    if all(rows[name]['status'] == 'computed' for name in TASK_IDS):
        archive = {r['task_id']: r for r in read_json(verify(m['archived_endpoints']))['rows']}
        components = {}
        for size in (33, 36):
            ca, la = f'1h4i_qm{size}_Ca', f'1h4i_qm{size}_La'
            dft = (archive[ca]['energy_hartree']-archive[la]['energy_hartree'])*HA_TO_KCAL
            mace = (rows[ca]['energy_eV']-rows[la]['energy_eV'])*EV_TO_KCAL
            components[f'qm{size}'] = {'DFT_core_R_kcal_mol': dft, 'MACE_core_R_kcal_mol': mace,
                                      'DFT_minus_MACE_R_kcal_mol': dft-mace}
        delta = components['qm36']['DFT_minus_MACE_R_kcal_mol'] - components['qm33']['DFT_minus_MACE_R_kcal_mol']
        result['core_partition'] = {'components': components, 'hybrid_partition_shift_kcal_mol': delta,
                                    'partition_check': abs(delta) <= m['tolerances']['partition_kcal_mol'],
                                    'identity': 'common full-system MACE term cancels exactly in the partition difference',
                                    'full_system_evaluated': False}
    full = {metal: rows[f'1h4i_full_{metal}_primary'] for metal in ('La', 'Ca')}
    if all(r['status'] == 'computed' for r in full.values()):
        if result['core_partition'] is not None:
            result['core_partition']['full_system_evaluated'] = True
        direct_ev = full['Ca']['energy_eV'] - full['La']['energy_eV']
        result['direct'] = {'R_eV': direct_ev, 'R_kcal_mol': direct_ev * EV_TO_KCAL}
        archive = {r['task_id']: r for r in read_json(verify(m['archived_endpoints']))['rows']}
        hybrid = {}
        for size in (33, 36):
            names = {metal: f'1h4i_qm{size}_{metal}' for metal in ('La', 'Ca')}
            if not all(rows[n]['status'] == 'computed' for n in names.values()):
                continue
            core_r_ev = rows[names['Ca']]['energy_eV'] - rows[names['La']]['energy_eV']
            dft_r_ha = archive[names['Ca']]['energy_hartree'] - archive[names['La']]['energy_hartree']
            hybrid[f'qm{size}'] = {'DFT_core_R_hartree': dft_r_ha, 'DFT_core_R_kcal_mol': dft_r_ha * HA_TO_KCAL,
                                  'MACE_core_R_eV': core_r_ev, 'MACE_core_R_kcal_mol': core_r_ev * EV_TO_KCAL,
                                  'MACE_environment_contribution_kcal_mol': (direct_ev-core_r_ev)*EV_TO_KCAL,
                                  'hybrid_R_kcal_mol': dft_r_ha*HA_TO_KCAL + (direct_ev-core_r_ev)*EV_TO_KCAL}
        result['hybrid'] = hybrid
        if len(hybrid) == 2:
            result['partition_shift_kcal_mol'] = hybrid['qm36']['hybrid_R_kcal_mol'] - hybrid['qm33']['hybrid_R_kcal_mol']
            result['partition_check'] = abs(result['partition_shift_kcal_mol']) <= m['tolerances']['partition_kcal_mol']
        checks = []
        for variant in ('repeat', 'translate', 'rotate'):
            pair = {}
            for metal in ('La', 'Ca'):
                name = f'1h4i_full_{metal}_{variant}'; r = rows.get(name, {'status': 'unavailable'})
                if r['status'] != 'computed':
                    continue
                t = next(t for t in m['tasks'] if t['task_id'] == name)
                delta = (r['energy_eV']-full[metal]['energy_eV'])*EV_TO_KCAL
                forces = np.load(verify(r['forces'])) @ np.array(t['rotation_matrix'])
                force_delta = float(np.max(np.abs(forces-np.load(verify(full[metal]['forces'])))))
                checks.append({'variant': variant, 'metal': metal, 'energy_delta_kcal_mol': delta,
                               'force_max_delta_eV_A': force_delta,
                               'pass': abs(delta) <= m['tolerances']['energy_kcal_mol'] and force_delta <= m['tolerances']['force_max_eV_A']})
                pair[metal] = r['energy_eV']
            if len(pair) == 2:
                delta = ((pair['Ca']-pair['La'])-direct_ev)*EV_TO_KCAL
                checks.append({'variant': variant, 'paired_R_delta_kcal_mol': delta,
                               'pass': abs(delta) <= m['tolerances']['energy_kcal_mol']})
        result['checks'] = checks
    if m.get('schema_version') == 'alquemia.mace_analytic.v1':
        from mace_analytic_pilot import add_analysis
        add_analysis(result, m)
    if m.get('schema_version') == 'alquemia.mace_response_trace.v1':
        from mace_response_trace import add_analysis
        add_analysis(result, m)
    return result


def report(collection, output):
    c = read_json(collection)
    if c.get('protocol_id') == 'mace_frozen_monopole_obc2_v1':
        from mace_gb import report_gb
        return report_gb(collection, output)
    if c.get('protocol_id') == 'mace_polar_1m_realspace_rotation_attribution_v1':
        from mace_rotation import report_rotation
        return report_rotation(collection, output)
    m = read_json(verify(c['manifest']))
    lines = ['# MACE / DFT vacuum capability pilot', '', f"Status: **{c['status']}**.", '',
             f"Protocol: `{c['protocol_id']}`. Baseline/default unchanged.", '',
             'These are consumed 1H4I development cases. No solution-phase reference,',
             'calibrated class or predictive-accuracy claim is available.', '',
             '| Task | Status | MACE E (eV) | Evaluation (s) | Peak GPU allocated (GiB) |',
             '|---|---|---:|---:|---:|']
    for name, r in c['rows'].items():
        e = r.get('energy_eV'); seconds = r.get('evaluation_seconds'); memory = r.get('peak_cuda_allocated_bytes')
        fields = [name, r['status'], 'unavailable' if e is None else repr(e),
                  'unavailable' if seconds is None else f'{seconds:.6f}',
                  'unavailable' if memory is None else f'{memory/2**30:.6f}']
        lines.append('| ' + ' | '.join(fields) + ' |')
    lines.extend(['', f"Attempts recorded: {len(c['attempts'])}; successful tasks: "
                  f"{sum(r['status']=='computed' for r in c['rows'].values())}/{len(m['tasks'])}.",
                  'New DFT endpoints: 0. '+('No matched DFT core reference available for changed geometry.' if c.get('geometry_only_pilot') else 'Four archived converged vacuum endpoints reused.'), ''])
    if c.get('geometry_changes'):
        lines.extend(['## Hydrogen preparation effects', '', '```json',
                      json.dumps({'direct': c['direct'], 'endpoint_changes': c['geometry_changes'],
                                  'hybrid_unavailable_reason': c['hybrid_unavailable_reason']}, indent=2), '```', ''])
    if c.get('core_partition'):
        cp = c['core_partition']
        lines.extend(['## Core-only partition identity', '',
                      'The common full-protein MACE term cancels in the partition difference.',
                      f"Hybrid partition shift: {cp['hybrid_partition_shift_kcal_mol']:.9f} kcal/mol.",
                      f"Frozen 2 kcal/mol diagnostic: {'PASS' if cp['partition_check'] else 'FAIL'}.",
                      f"Full-system pair evaluated: {cp['full_system_evaluated']}.",
                      'This identity does not measure full-system runtime or validate a biological prediction.', ''])
    if c.get('hybrid'):
        lines.extend(['## Paired components', '',
                      '| Partition | DFT core R | MACE core R | MACE full-minus-core R | Hybrid R |',
                      '|---|---:|---:|---:|---:|'])
        for name, r in c['hybrid'].items():
            values = [r[k] for k in ('DFT_core_R_kcal_mol','MACE_core_R_kcal_mol',
                                     'MACE_environment_contribution_kcal_mol','hybrid_R_kcal_mol')]
            lines.append('| '+name+' | '+' | '.join(f'{v:.9f}' for v in values)+' |')
        lines.extend(['', 'All table components are kcal/mol; R = E(Ca) - E(La).',
                      f"Partition shift: {c.get('partition_shift_kcal_mol', 'unavailable')} kcal/mol; "
                      f"frozen {m['tolerances']['partition_kcal_mol']} kcal/mol check: {c.get('partition_check', 'unavailable')}.", ''])
    if c.get('checks'):
        lines.extend(['## Repeat and rigid-transform checks', '', '```json',
                      json.dumps(c['checks'], indent=2), '```', ''])
    if c.get('analytic_comparison'):
        lines.extend(['## Analytic-kernel core checks and method changes', '', '```json',
                      json.dumps(c['analytic_comparison'], indent=2), '```', ''])
    if c.get('trace_equivalence'):
        lines.extend(['## Read-only response trace equivalence', '', '```json',
                      json.dumps(c['trace_equivalence'], indent=2), '```', ''])
    lines.extend(['## Interpretation', '',
                  '- Numerical credibility: inspect charge, rigid-transform and partition checks separately.',
                  '- Scientific informativeness: this single consumed vacuum system cannot establish accuracy.',
                  '- Affordability: task timings exclude queue wait; whole-job accounting includes startup/failures.',
                  '- Structural response and solution-phase scoring remain unavailable.',
                  '- Retain the baseline. No automatic promotion follows this report.', ''])
    with Path(output).open('x') as f:
        f.write('\n'.join(lines))
    return {'status': 'report_written', 'report': record(output)}


def main():
    p = argparse.ArgumentParser(description=__doc__); sp = p.add_subparsers(dest='command', required=True)
    q = sp.add_parser('prepare')
    for key in ('root', 'software', 'output', 'agreement'): q.add_argument('--'+key, required=True)
    for command in ('dry-run', 'collect', 'compare-cores'):
        q = sp.add_parser(command); q.add_argument('--manifest', required=True); q.add_argument('--output')
    q = sp.add_parser('execute'); q.add_argument('--manifest', required=True)
    q.add_argument('--memory-mode', choices=('native', 'host_offload'), default='native'); q.add_argument('--task-id', action='append')
    q = sp.add_parser('worker')
    for key in ('manifest', 'task-id', 'output'): q.add_argument('--'+key, required=True)
    q.add_argument('--memory-mode', choices=('native', 'host_offload'), required=True)
    q = sp.add_parser('report'); q.add_argument('--collection', required=True); q.add_argument('--output', required=True)
    q = sp.add_parser('repair-interface'); q.add_argument('--manifest', required=True); q.add_argument('--output', required=True)
    q.add_argument('--pair-tile', type=int); q.add_argument('--memory-agreement'); q.add_argument('--reference-collection')
    q.add_argument('--edge-tile', type=int)
    q.add_argument('--node-tile', type=int)
    a = p.parse_args()
    if a.command == 'prepare': r = prepare(a.root, a.software, a.output, a.agreement)
    elif a.command == 'worker': r = worker(a.manifest, a.task_id, a.output, a.memory_mode)
    elif a.command == 'execute': r = execute(a.manifest, a.memory_mode, a.task_id)
    elif a.command == 'report': r = report(a.collection, a.output)
    elif a.command == 'repair-interface': r = repair_interface(a.manifest, a.output, a.pair_tile, a.memory_agreement, a.reference_collection, a.edge_tile, a.node_tile)
    else:
        if a.command == 'dry-run': r = dry_run(a.manifest)
        elif a.command == 'compare-cores': r = compare_cores(a.manifest)
        else: r = collect(a.manifest)
        if a.output: write_new(a.output, r)
    print(json.dumps({k: r[k] for k in ('status', 'tasks', 'manifest', 'failed_task') if k in r}, indent=2))
    if r.get('status') in ('failed', 'partial_failure'): sys.exit(1)


if __name__ == '__main__':
    main()
