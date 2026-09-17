"""Saved-output fixed-field DFT/MACE screen; no scientific executor or solvent score."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL
from mace_curvature import actual_collection
from density_embedding import collect_chargefit, collect_potentials
from global_electrostatic import collect_endpoints

PROTOCOL = 'native_DFT_CHELPG_MACE_short_vacuum_screen_v1'
TOL = {'charge_e': 1e-5, 'coupling_kcal_mol': .01, 'rigid_kcal_mol': .01,
       'partition_kcal_mol': 2., 'isolated_identity_kcal_mol': 1e-8}
IDS = {f'1h4i_qm{size}_{metal}' for size in (33, 36) for metal in ('La', 'Ca')}


def archive_equal(saved, current):
    """Preserve exact metadata; allow only floating summation roundoff on replay."""
    if isinstance(saved, dict) and isinstance(current, dict):
        return saved.keys() == current.keys() and all(archive_equal(saved[k], current[k]) for k in saved)
    if isinstance(saved, list) and isinstance(current, list):
        return len(saved) == len(current) and all(archive_equal(a, b) for a, b in zip(saved, current))
    if isinstance(saved, float) and isinstance(current, float):
        return math.isfinite(saved) and math.isfinite(current) and abs(saved-current) <= 1e-9
    return type(saved) is type(current) and saved == current


def direct_coupling(state, charges):
    core = np.array([a['xyz_A'] for a in state['core_atoms']])
    env = np.array([a['xyz_A'] for a in state['environment_atoms']]).reshape(-1, 3)
    if not len(env):
        return 0.
    distances = np.linalg.norm(env[:, None, :] - core[None, :, :], axis=2)
    if distances.min() < 1.:
        raise InvalidArtifact('core/environment charge overlap')
    weights = np.array([a['charge_e'] for a in state['environment_atoms']])
    return float(weights @ (np.sum(np.array(charges)[None, :] / distances, axis=1)) * BOHR_TO_A * HA_TO_KCAL)


def coordinate_check(atoms, expected, metal):
    elements = [metal if a['element'] == 'M' else a['element'] for a in expected]
    if [a[0] for a in atoms] != elements:
        raise InvalidArtifact('source physical/core atom ordering differs')
    if np.max(np.abs(np.array([a[1:] for a in atoms])-np.array([a['xyz_A'] for a in expected]))) > 1e-10:
        raise InvalidArtifact('source physical/core coordinates differ')


def contrast(ca, la):
    dft = (ca['DFT_vacuum_hartree']-la['DFT_vacuum_hartree'])*HA_TO_KCAL
    coupling = ca['CHELPG_coupling_kcal_mol']-la['CHELPG_coupling_kcal_mol']
    exact = ca['density_coupling_kcal_mol']-la['density_coupling_kcal_mol']
    short_full = (ca['short_full_eV']-la['short_full_eV'])*EV_TO_KCAL
    short_core = (ca['short_core_eV']-la['short_core_eV'])*EV_TO_KCAL
    return {'DFT_R_kcal_mol': dft, 'CHELPG_R_kcal_mol': coupling,
            'density_R_kcal_mol': exact, 'short_full_R_kcal_mol': short_full,
            'short_core_R_kcal_mol': short_core, 'short_correction_R_kcal_mol': short_full-short_core,
            'R_vacuum_screen_kcal_mol': dft+coupling+short_full-short_core,
            'density_diagnostic_R_kcal_mol': dft+exact+short_full-short_core,
            'solution_score': None, 'calibrated_class': None}


def compare(medium, large, charge_fit, agreement, output):
    started = time.perf_counter(); cpu = time.process_time()
    saved_fit = read_json(charge_fit); cm = verify(saved_fit['manifest'])
    current_fit = collect_chargefit(cm)
    if not archive_equal(saved_fit, current_fit) or current_fit['status'] != 'complete':
        raise InvalidArtifact('charge-fit archive differs from actual execution')
    charge_manifest = read_json(cm); saved_density = read_json(verify(charge_manifest['potential_result']))
    density = collect_potentials(verify(saved_density['manifest']))
    if not archive_equal(saved_density, density) or density['status'] != 'complete':
        raise InvalidArtifact('density archive differs from actual execution')
    fit_rows = {r['task_id']: r for r in current_fit['rows']}
    density_rows = {r['task_id']: r for r in density['rows']}
    if set(fit_rows) != IDS or set(density_rows) != IDS:
        raise InvalidArtifact('fixed four-state source inventory differs')
    states = {t['task_id']: read_json(verify(t['source_potential_task']['state'])) for t in charge_manifest['tasks']}
    state_refs = {t['task_id']: t['source_potential_task']['state'] for t in charge_manifest['tasks']}
    checks = []; models = {}; manifests = {}
    for label, path in (('medium', medium), ('large', large)):
        c, m = actual_collection(path); manifests[label] = m
        if {r['sha256'] for r in m['source_states']} != {r['sha256'] for r in state_refs.values()}:
            raise InvalidArtifact('MACE and charge-fit physical preparations differ')
        for ref in [*m['implementation'].values(), m['software'], m['model']['checkpoint'], m['atom_mappings']]:
            verify(ref)
        endpoints = collect_endpoints(verify(m['source_endpoint_manifest']))
        if endpoints['status'] != 'endpoints_complete' or not archive_equal(endpoints, read_json(verify(m['archived_endpoints']))):
            raise InvalidArtifact('DFT source lacks matching actual archived execution')
        dft = {r['task_id']: r for r in endpoints['rows']}
        tasks = {t['task_id']: t for t in m['tasks']}; rows = {}
        for name in sorted(IDS):
            s = states[name]; metal = s['metal']; core_task = tasks[name]
            full_id = '1h4i_full_'+metal+'_primary'; full_task = tasks[full_id]
            coordinate_check(xyz(verify(core_task['xyz'])), s['core_atoms'], metal)
            coordinate_check(xyz(verify(full_task['xyz'])), s['physical_atoms'], metal)
            if (core_task['charge'] != s['core_total_charge_e'] or
                    abs(full_task['charge']-s['core_total_charge_e']-s['expected_environment_charge_e']) > TOL['charge_e'] or
                    core_task['spin_multiplicity'] != 1 or full_task['spin_multiplicity'] != 1 or s['explicit_waters']):
                raise InvalidArtifact('charge, spin or water state changed')
            if {a['id'] for a in s['core_atoms']} & {a['id'] for a in s['environment_atoms']}:
                raise InvalidArtifact('duplicate QM/environment charge owner')
            q = fit_rows[name]['charge_e']; charge_error = sum(q)-s['core_total_charge_e']
            coupling = direct_coupling(s, q)
            fit_error = coupling-fit_rows[name]['direct_fit_kcal_mol']
            checks += [{'name': label+'_'+name+'_charge', 'error_e': charge_error,
                        'pass': abs(charge_error) <= TOL['charge_e']},
                       {'name': label+'_'+name+'_coupling_replay', 'error_kcal_mol': fit_error,
                        'pass': abs(fit_error) <= TOL['coupling_kcal_mol']}]
            rows[name] = {'DFT_vacuum_hartree': dft[name]['energy_hartree'],
                          'CHELPG_coupling_kcal_mol': coupling,
                          'density_coupling_kcal_mol': density_rows[name]['direct_density_kcal_mol'],
                          'short_core_eV': c['rows'][name]['energy_components_eV']['interaction_energy'],
                          'short_full_eV': c['rows'][full_id]['energy_components_eV']['interaction_energy'],
                          'source_state': state_refs[name], 'DFT_receipt': dft[name]['execution_receipt'],
                          'charge_receipt': fit_rows[name]['execution_receipt'],
                          'density_receipt': density_rows[name]['execution_receipt'],
                          'core_xyz': core_task['xyz'], 'full_xyz': full_task['xyz']}
            isolated = dict(s, environment_atoms=[])
            core_e = rows[name]['short_core_eV']
            identity = (core_e-core_e)*EV_TO_KCAL+direct_coupling(isolated, q)
            checks.append({'name': label+'_'+name+'_algebraic_isolated_identity', 'error_kcal_mol': identity,
                           'pass': abs(identity) <= TOL['isolated_identity_kcal_mol']})
        for size in (33, 36):
            la, ca = (states[f'1h4i_qm{size}_{metal}'] for metal in ('La', 'Ca'))
            if la['environment_atoms'] != ca['environment_atoms'] or la['assembly'] != ca['assembly']:
                raise InvalidArtifact('paired permanent environment changed')
        for task in m['tasks']:
            if task['variant'] != 'rotate':
                continue
            name = task['task_id']; primary = name.replace('_rotate', '_primary' if task['kind']=='full' else '')
            error = (c['rows'][name]['energy_components_eV']['interaction_energy']-
                     c['rows'][primary]['energy_components_eV']['interaction_energy'])*EV_TO_KCAL
            checks.append({'name': label+'_'+name+'_saved_short_rigid', 'error_kcal_mol': error,
                           'pass': abs(error) <= TOL['rigid_kcal_mol']})
        paired = {f'qm{size}': contrast(rows[f'1h4i_qm{size}_Ca'], rows[f'1h4i_qm{size}_La']) for size in (33, 36)}
        delta = {k: paired['qm36'][k]-paired['qm33'][k] for k in paired['qm33'] if k.endswith('_kcal_mol')}
        # Evaluate the equivalent full-term cancellation directly to expose sign errors.
        closure = delta['DFT_R_kcal_mol']+delta['CHELPG_R_kcal_mol']-delta['short_core_R_kcal_mol']
        checks.append({'name': label+'_partition_component_closure',
                       'error_kcal_mol': closure-delta['R_vacuum_screen_kcal_mol'],
                       'pass': abs(closure-delta['R_vacuum_screen_kcal_mol']) <= 1e-8})
        models[label] = {'source': record(path), 'checkpoint': m['model']['checkpoint'], 'rows': rows,
                         'paired': paired, 'partition': delta,
                         'partition_gate_pass': abs(delta['R_vacuum_screen_kcal_mol']) <= TOL['partition_kcal_mol']}
    for key in ('source_states', 'source_endpoint_manifest', 'archived_endpoints', 'atom_mappings'):
        if manifests['medium'][key] != manifests['large'][key]:
            raise InvalidArtifact('checkpoint comparison does not preserve physical source')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    implementation = out/Path(__file__).name; shutil.copyfile(__file__, implementation)
    result = {'status': 'complete', 'protocol_id': PROTOCOL, 'agreement': record(agreement),
              'implementation': record(implementation), 'charge_fit': record(charge_fit),
              'density_source': charge_manifest['potential_result'], 'tolerances': TOL,
              'models': models, 'checks': checks, 'accounting_checks_pass': all(c['pass'] for c in checks),
              'primary_checkpoint': 'medium', 'primary_partition_gate_pass': models['medium']['partition_gate_pass'],
              'solution_score': None, 'solvent_status': 'not_implemented_in_vacuum_screen',
              'calibrated_class': None, 'new_DFT_calls': 0, 'new_MACE_calls': 0, 'new_solver_calls': 0,
              'baseline_changed': False, 'prospectively_blind': False,
              'analysis_wall_seconds': time.perf_counter()-started, 'analysis_CPU_seconds': time.process_time()-cpu}
    write_new(out/'result.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('medium', 'large', 'charge-fit', 'agreement', 'output'):
        p.add_argument('--'+key, required=True)
    result = compare(**vars(p.parse_args()))
    print(json.dumps({'accounting_checks_pass': result['accounting_checks_pass'],
                      'models': {k: {'partition': v['partition'], 'gate': v['partition_gate_pass']} for k,v in result['models'].items()}}, indent=2))
