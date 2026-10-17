#!/usr/bin/env python3
"""One certified saved-component SCC recovery, then 37 untouched pool cells."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import time
import traceback

PROTOCOL = 'lbt_complete_pool_single_isolated_solver_recovery_v1'
POOL_SHA = '43f081396c56351ea5832b4187c27874ec708dbfbe9f6577911656f2b7ad9b02'
SCOUT_SHA = 'd1ecf6c53f0d5047a4273cf5a7dfde7c808c9002cbd329c6721ea220036fe460'
EVALUATOR_SHA = '686eba84fd3dfcf28724d2d4e7957528c987445d7137b18c88a07dfec055f661'
FAILED_SHA = '81edb534dcd3dd2447f31d410d20fa07d02dc59ec442c4639328799c9bdc6768'
PRIMARY_SHA = 'f40dbdf2cd3c45bb89deb42400fa9d45898d428c8c6a92d18cd234f755c47651'
TARGET = ('cartesian_03', 'baseline', 'La')
POLICY = {'target': list(TARGET), 'source_id': '7CCN', 'primary_case': 'failed_005',
          'complete_prior_cells': 66, 'saved_component_replacement_cells': 1,
          'maximum_new_composites': 37, 'full_cells_required': 104,
          'new_cells_use_original_solver': True, 'stop_at_first_new_backend_failure': True}
_POOL = None


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        sha = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    return pin(path)


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def normalized(value):
    """Verified byte identity permits immutable snapshots of the same inputs."""
    if isinstance(value, dict):
        if set(value) == {'path', 'sha256'}:
            verify(value)
            return {'sha256': value['sha256']}
        return {k: normalized(v) for k, v in value.items()}
    if isinstance(value, list):
        return [normalized(v) for v in value]
    return value


def binding(cfg):
    value = normalized({k: v for k, v in cfg.items() if k != 'recovery_preflight_result'})
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def slot_key(slot):
    return slot['candidate_id'], slot['region'], slot['metal']


def cell_from_point(slot, point, reused):
    return _POOL.cell_from_point(slot, point, reused)


def menu(*args, **kwargs):
    return _POOL.menu(*args, **kwargs)


def summarize(*args, **kwargs):
    return _POOL.summarize(*args, **kwargs)


def admit_saved_cell(ev, qualifier, slot, candidate):
    point = qualifier.restore(ev, slot['point'])
    require((point['source_id'], point['region'], point['metal']) ==
            (candidate['source_id'], slot['region'], slot['metal'])
            and ev.np.array_equal(point['positions_A'], candidate['_x'])
            and all(slot[k] == point[k] for k in
                    ('native_energies', 'component_energies_kcal_mol', 'total_energy_kcal_mol'))
            and slot['cache_key'] == point['provenance']['cache_key'], 'Saved pool cell differs')
    return point


def certify_scc(cfg, out, partial, failed_slot):
    """Rerun saved native admission and the three frozen agreement checks; no calls."""
    scout = module('recovery_saved_scc', cfg['scc_implementation'])
    require(cfg['scc_implementation']['sha256'] == SCOUT_SHA, 'SCC diagnostic source differs')
    result = read(cfg['scc_result'])
    scfg = read(result['config'])
    require(result['protocol'] == scout.PROTOCOL and result['stage'] == 'run'
            and result['status'] == 'complete' and result['qualification_passed']
            and result['implementation']['sha256'] == scfg['implementation']['sha256'] == SCOUT_SHA
            and result['slurm_job_id'] == '13009' and result['declared_calls'] == 4
            and result['actual_attempted_GFN_calls'] == result['completed_processes'] == 4
            and result['new_MACE_MM_ALPB_DFT_calls'] == 0
            and normalized(scfg['original_failed_pool']) == normalized(cfg['original_partial_result']),
            'Actual four-call SCC qualification unavailable')
    require(scfg['fixtures']['failed']['point'] == failed_slot['attempted_point'] == partial['requests'][-1]
            and scfg['fixtures']['failed']['point']['sha256'] == FAILED_SHA
            and scfg['fixtures']['preceding']['point'] == partial['requests'][-2], 'SCC fixtures not linked to exact partial pool')
    d = out / 'scc_saved_replay'
    d.mkdir()
    scratch = {}
    helpers = scout.prepare(scfg, d, scratch)
    backend, gradients, admission, _, np, fixtures = helpers
    require(len(result['cells']) == 4, 'Four SCC diagnostic cells required')
    replayed = {}
    for row, (case_id, fixture_name, damping) in zip(result['cells'], scout.CASES):
        require((row['case_id'], row['fixture'], row['broydamp']) == (case_id, fixture_name, damping)
                and row['status'] == 'accepted' and row['attempted'] is True
                and row['printed_initial_broydamp'] == damping, 'SCC case state differs')
        scout.verify_tree(row)
        require(read(row['receipt']) == {k: v for k, v in row.items() if k != 'receipt'}, 'SCC row receipt differs')
        fixture = fixtures[fixture_name]
        task, execution, native = read(row['task']), read(row['execution']), read(row['native_receipt'])
        staged = task['staged_files']
        directory = verify(row['task']).parent
        require(set(staged) == {'core.xyz', 'parameters', 'control.inp'}
                and staged['core.xyz']['sha256'] == fixture['task']['xyz']['sha256']
                and staged['parameters']['sha256'] == scfg['runtime']['parameters']['sha256']
                and verify(staged['control.inp']).read_text() == scout.controls_text(damping),
                'Changed SCC input bytes/control')
        command = [str(verify(scfg['runtime']['executable']))] + scout.ARGS + ['--input', 'control.inp']
        expected = copy.deepcopy(fixture['task'])
        expected.update(xyz=staged['core.xyz'], command=command, staged_files=staged,
                        original_task=fixture['native']['task'], solver_diagnostic_protocol=scout.PROTOCOL,
                        solver_control={'broydamp': damping, 'control': staged['control.inp'],
                                        'native_adaptive_behavior_unchanged': True})
        require(task == expected and execution['command'] == command, 'SCC changed settings beyond declared damping')
        values = re.findall(r'Broyden damping\s+([-+\d.Ee]+)', (directory / 'stdout.txt').read_text())
        require(values and float(values[-1]) == damping, 'Actual requested damping not replayed')
        admitted = admission.admit(directory, task, execution, backend=backend, gradients=gradients, np=np)
        require(admitted == row['native'] and all(native[k] == v for k, v in admitted.items())
                and native['status'] == 'accepted' and native['attempted'] is True
                and native['task'] == row['task'] and native['execution'] == row['execution']
                and native['solver_diagnostic_protocol'] == scout.PROTOCOL
                and native['solver_control'] == task['solver_control']
                and native['original_native_result'] == fixture['pins']['native_result']
                and native['original_failure_not_relabelled'] is True, 'Flat native certificate differs')
        scout.verify_tree(native)
        replayed[case_id] = admitted
    checks = []
    for name, left, right, fixture_name in (
            ('target_branch_consistency', 'failed_005', 'failed_010', 'failed'),
            ('preceding_control', 'preceding_005', None, 'preceding'),
            ('origin_control', 'origin_005', None, 'origin')):
        reference = fixtures[fixture_name]['reference'] if right is None else replayed[right]
        check = {'name': name, 'left': left,
                 'right': right or scfg['fixtures'][fixture_name]['native_result'], 'status': 'complete'}
        check.update(scout.compare(replayed[left], reference, fixtures[fixture_name]['mapping'], helpers))
        checks.append(check)
    require(checks == result['checks'] and all(c['pass'] for c in checks), 'Frozen SCC agreement checks do not replay')
    primary = result['cells'][0]['native_receipt']
    require(primary['sha256'] == PRIMARY_SHA and result['cells'][0]['case_id'] == POLICY['primary_case'],
            'Predeclared primary isolated result changed')
    evidence = {'result': cfg['scc_result'], 'implementation': cfg['scc_implementation'],
                'config': result['config'], 'specification': scfg['specification'],
                'primary': primary, 'secondary': result['cells'][1]['native_receipt'],
                'all_four_native_receipts': [r['native_receipt'] for r in result['cells']],
                'checks': checks, 'agreement_gates': scfg['agreement_gates'],
                'new_molecular_calls': 0, 'no_ground_state_or_global_solver_stability_claim': True}
    return evidence, fixtures, helpers


def assemble(cfg, base, ev, failed_slot, candidate, evidence, fixtures, helpers):
    """The V2 component arithmetic, using only existing raw artifacts."""
    np = ev.np
    failed = read(failed_slot['attempted_point'])
    require(failed['status'] == 'unavailable' and 'point' not in failed
            and failed['identity_sha256'] == ev.identity_hash
            and failed['state'] == ev.states[('7CCN', 'baseline', 'La')]
            and failed['cache_key'] == ev.key('7CCN', 'baseline', 'La', candidate['_x']), 'Failed-point state differs')
    require(set(failed['component_receipts']) == {'MACE', 'MM', 'MM_execution', 'GFN_embedded', 'GFN_isolated'},
            'Unexpected failed-point component inventory')
    base.verify_tree(failed['component_receipts'])
    mapping, task, positions = read(failed['mapping']), read(failed['task']), read(failed['positions'])
    require(np.array_equal(positions['positions_A'], candidate['_x'])
            and positions['source_ids'] == [a['id'] for a in mapping['physical_atoms']]
            and read(task['physical_source']) == mapping['physical_atoms']
            and task['physical_source']['sha256'] == candidate['physical_source']['sha256']
            and failed['source_task'] == ev.tasks[('7CCN', 'baseline', 'La')]
            and task['origin_task'] == failed['source_task'], 'Failed physical/task lineage differs')
    sy, xyz = ev.admit.xyz(task['xyz'])
    require(sy == ['La'] + [a['element'] for a in mapping['core'][1:]]
            and np.array_equal(xyz, [a['xyz_A'] for a in mapping['core']]), 'Failed current core XYZ differs')
    receipts = copy.deepcopy(failed['component_receipts'])
    ml, mm, embedded = [read(receipts[k]) for k in ('MACE', 'MM', 'GFN_embedded')]
    backend, gradients, admission, _, _, _ = helpers
    etask, execution = read(embedded['task']), read(embedded['execution'])
    admission_replay = admission.admit(verify(embedded['task']).parent, etask, execution,
                                      backend=backend, gradients=gradients, np=np)
    require(embedded['status'] == 'accepted' and all(embedded[k] == v for k, v in admission_replay.items())
            and etask['kind'] == 'embedded' and etask['origin_task'] == failed['source_task']
            and etask['xyz']['sha256'] == task['xyz']['sha256'], 'Old embedded native component unavailable')
    require(ml['task'] == failed['task'] and ml['xyz'] == task['xyz'] and ml['physical_source'] == task['physical_source']
            and ml['input_state_check']['charge'] == task['charge']
            and ml['input_state_check']['spin_multiplicity'] == task['multiplicity']
            and ml['input_state_check']['atoms'] == task['atom_count'] == 172
            and ml['force_convention'] == 'negative vacuum energy derivative', 'Old MACE state differs')
    mex = read(receipts['MM_execution'])
    require(mex['returncode'] == 0 and mex['owned_process_group_empty']
            and mm['status'] == 'complete' and mm['positions'] == failed['positions']
            and mm['system'] == ev.partitions[('7CCN', 'baseline')][0]['systems']['La']['system']
            and mm['implementation']['sha256'] == EVALUATOR_SHA
            and mm['configurations_attempted'] == 1 and mm['group_queries_attempted'] == 4,
            'Old MM ledger unavailable')
    q = read(task['global_solvent_charges'])['charges_e']
    solvent_key = ev.solvent_key('7CCN', 'La', candidate['_x'], q)
    require(solvent_key in ev.solvent_cache, 'Exact full-system solvent missing; no fallback call')
    sol = ev.solvent_cache[solvent_key]
    receipts['GFN_isolated'] = evidence['primary']
    receipts['ALPB'] = sol['record']
    isolated = read(evidence['primary'])
    zero = np.zeros((len(mapping['field']), 3))
    eg = ev.mapping.scatter(mapping, np.asarray(embedded['qm_gradient']['hartree_per_bohr']),
                            np.asarray(embedded['pointcharge_gradient']['hartree_per_bohr'])) * base.HA / base.BOHR
    ig = ev.mapping.scatter(mapping, np.asarray(isolated['qm_gradient']['hartree_per_bohr']), zero) * base.HA / base.BOHR
    gradients = {'MACE_vacuum': ev.mapping.scatter(mapping, -np.asarray(ml['forces_eV_A']) * base.EV, zero),
                 'GFN2_field_difference': eg - ig, 'MM_ledger': np.asarray(mm['gradient_kcal_mol_A']),
                 'fixed_FF_charge_ALPB': np.asarray(sol['gradient_hartree_per_bohr']) * base.HA / base.BOHR}
    native = {'MACE_eV': ml['energy_eV'], 'GFN_embedded_Eh': embedded['energy_hartree'],
              'GFN_isolated_Eh': isolated['energy_hartree'], 'MM_groups_kcal_mol': mm['component_energies_kcal_mol'],
              'ALPB_Eh': sol['energy_hartree']}
    components = {'MACE_vacuum': native['MACE_eV'] * base.EV,
                  'GFN2_field_difference': (native['GFN_embedded_Eh'] - native['GFN_isolated_Eh']) * base.HA,
                  'MM_ledger': mm['energy_kcal_mol'], 'fixed_FF_charge_ALPB': native['ALPB_Eh'] * base.HA}
    gradient = sum(gradients[k] for k in base.COMPONENTS)
    require(gradient.shape == (256, 3) and np.isfinite(gradient).all()
            and all(math.isfinite(v) for v in components.values()), 'Saved assembly is nonfinite')
    point = {'positions_A': candidate['_x'].tolist(), 'gradient_kcal_mol_A': gradient.tolist(),
             'total_energy_kcal_mol': math.fsum(components.values()), 'component_energies_kcal_mol': components,
             'component_gradients_kcal_mol_A': {k: gradients[k].tolist() for k in base.COMPONENTS},
             'native_energies': native, 'source_id': '7CCN', 'region': 'baseline', 'metal': 'La', 'state': failed['state']}
    certificate = {'protocol': PROTOCOL, 'binding_sha256': binding(cfg), 'policy': POLICY,
                   'Hamiltonian_identity_sha256': ev.identity_hash, 'Hamiltonian_cache_key': failed['cache_key'],
                   'same_Hamiltonian': True, 'identical_solver_protocol': False,
                   'original_failed_pool': cfg['original_partial_result'], 'original_failed_slot': failed_slot,
                   'original_failed_POINT': failed_slot['attempted_point'], 'old_failure_remains_unavailable': True,
                   'original_isolated_component': failed['component_receipts']['GFN_isolated'],
                   'saved_scc_qualification': evidence, 'primary_selection': 'predeclared failed_005, not energy selection',
                   'state': failed['state'], 'positions': failed['positions'], 'mapping': failed['mapping'], 'task': failed['task'],
                   'component_receipts': receipts, 'solvent_cache_key': solvent_key,
                   'assembly_implementation': cfg['implementation'], 'unchanged_Hamiltonian_implementation': ev.cfg['implementation'],
                   'new_molecular_calls_for_assembly': 0,
                   'solver_exception_scope': 'Only this exact isolated component; all remaining cells use unchanged V2 solver'}
    raw = {k: copy.deepcopy(failed[k]) for k in
           ('cache_key', 'identity_sha256', 'state', 'positions', 'mapping', 'source_task', 'task')}
    raw.update(status='complete', point=point, component_receipts=receipts,
               ALPB_reuse='exact_admitted_Module2_cell', assembly_protocol=PROTOCOL,
               original_failed_POINT=failed_slot['attempted_point'], new_molecular_calls=0, elapsed_seconds=0.)
    return certificate, raw


def prepare(cfg, out, result):
    """Saved-only parent API; certifies the archived exception before any restore."""
    global _POOL
    require(cfg['protocol'] == PROTOCOL and cfg['implementation']['sha256'] == pin(__file__)['sha256']
            and cfg['recovery_policy'] == POLICY, 'Recovery source/declaration differs')
    original_cfg = read(cfg['original_pool_config'])
    partial = read(cfg['original_partial_result'])
    require(original_cfg['implementation']['sha256'] == POOL_SHA
            and partial['implementation']['sha256'] == POOL_SHA
            and normalized(partial['config']) == normalized(cfg['original_pool_config'])
            and partial['status'] == 'unavailable_backend' and partial['stage'] == 'run'
            and partial['slurm_job_id'] == '12992' and len(partial['cells']) == 104
            and len(partial['requests']) == 59 and len(partial['points']) == 58,
            'Exact 12992 partial execution required')
    for k, value in original_cfg.items():
        if k not in ('implementation', 'specification', 'protocol'):
            require(normalized(cfg[k]) == normalized(value), 'Original scientific config changed: ' + k)
    _POOL = module('certified_recovery_original_pool', original_cfg['implementation'])
    require(partial['protocol'] == original_cfg['protocol'] == _POOL.PROTOCOL, 'Original pool protocol differs')
    verify(cfg['specification'])
    prior_dir = out / 'original_preparation_replay'
    prior_dir.mkdir()
    scratch = {}
    try:
        base, ev, candidates, cells, old_cells, whole = _POOL.prepare(original_cfg, prior_dir, scratch)
    finally:
        result.update(scratch)
    require(base.pin(base.__file__)['sha256'] == EVALUATOR_SHA and ev.identity_hash == partial['identity_sha256']
            and scratch['candidates'] == partial['candidates'] and len(candidates) == 26,
            'Original model/candidate identity differs')
    qualifier = module('certified_recovery_raw_restore', cfg['qualification_implementation'])
    expected = [slot_key(r) for r in scratch['cells']]
    actual = [slot_key(r) for r in partial['cells']]
    require(expected == actual and len(set(actual)) == 104, 'Original ordered full denominator differs')
    by_id = {c['candidate_id']: c for c in candidates}
    complete = [r for r in partial['cells'] if r['status'] == 'complete']
    failed = [r for r in partial['cells'] if r['status'] == 'unavailable_attempted']
    pending = [r for r in partial['cells'] if r['status'] == 'unavailable_not_attempted']
    require(len(complete) == 66 and len(failed) == 1 and len(pending) == 37
            and slot_key(failed[0]) == TARGET and failed[0]['source_id'] == '7CCN'
            and [r['point'] for r in complete if not r['reused']] == partial['points'] == partial['requests'][:-1],
            '66 accepted / one failed / 37 unattempted ledger differs')
    for slot in complete:
        admit_saved_cell(ev, qualifier, slot, by_id[slot['candidate_id']])
    cells = {slot_key(r): copy.deepcopy(r) for r in partial['cells']}
    require(not summarize(base, candidates, cells, old_cells, whole, result)
            and result['comparison'] == partial['comparison'], 'Original unavailable comparison does not replay')
    evidence, fixtures, helpers = certify_scc(cfg, out, partial, failed[0])
    certificate, raw = assemble(cfg, base, ev, failed[0], by_id[TARGET[0]], evidence, fixtures, helpers)
    previous = cfg.get('recovery_preflight_result')
    if previous is None:
        certificate_pin = write(out / 'SOLVER_EXCEPTION_CERTIFICATE.json', certificate)
        raw['solver_exception_certificate'] = certificate_pin
        raw['point']['solver_execution_provenance'] = certificate_pin
        point_pin = write(out / 'ASSEMBLED_POINT.json', raw)
    else:
        preflight = read(previous)
        pcfg = read(preflight['config'])
        require(preflight['protocol'] == PROTOCOL and preflight['status'] == 'complete'
                and preflight['stage'] == 'preflight' and preflight['implementation']['sha256'] == cfg['implementation']['sha256']
                and pcfg.get('recovery_preflight_result') is None and binding(pcfg) == binding(cfg)
                and preflight['actual_cost']['composites_attempted'] == 0,
                'Actual zero-call recovery preflight required')
        certificate_pin, point_pin = preflight['solver_exception_certificate'], preflight['assembled_point']
        require(read(certificate_pin) == certificate, 'Archived solver certificate no longer replays')
        raw['solver_exception_certificate'] = certificate_pin
        raw['point']['solver_execution_provenance'] = certificate_pin
        require(read(point_pin) == raw, 'Archived replacement raw components/assembly no longer replay')
    # Explicit certificate above authorizes exactly this immutable saved POINT.
    # Frozen restore performs every ordinary state/map/raw-gradient/cache gate.
    point = qualifier.restore(ev, point_pin)
    require(point['solver_execution_provenance'] == certificate_pin, 'Exception provenance was lost')
    slot = {k: failed[0][k] for k in ('candidate_id', 'source_id', 'region', 'metal', 'physical_source')}
    cell_from_point(slot, point, True)
    slot.update(solver_exception_certificate=certificate_pin, original_failed_slot=failed[0],
                reused_components_without_new_molecular_calls=True)
    cells[TARGET] = slot
    pending_keys = [list(slot_key(r)) for r in pending]
    require(sum(r['status'] == 'complete' for r in cells.values()) == 67
            and all(ev.key(r['source_id'], r['region'], r['metal'], by_id[r['candidate_id']]['_x']) not in ev.cache
                    for r in pending)
            and all(ev.cost[k] == 0 for k in ('composites_attempted', 'MACE_calls_attempted', 'GFN_calls_attempted',
                    'MM_configurations_attempted', 'MM_group_queries_attempted', 'ALPB_calls_attempted'))
            and ev.native is None and len(result['required_prior_solvent_cells']) == 52
            and all(r['status'] == 'complete' for r in result['required_prior_solvent_cells']),
            'Recovery preflight repeated work or lost the exact pending denominator')
    result.update(cells=[cells[k] for k in expected], compatible_hybrid_cells_reused=67, maximum_new_composites=37,
                  component_call_ceiling={'MACE': 37, 'GFN2': 74, 'MM_configurations': 37, 'MM_group_queries': 148, 'ALPB': 0},
                  original_partial_result=cfg['original_partial_result'], original_failure_preserved=failed[0],
                  original_preflight=result['preflight'], solver_exception_certificate=certificate_pin,
                  assembled_point=point_pin, remaining_original_slots=pending_keys,
                  recovery_binding_sha256=binding(cfg), actual_assembly_molecular_calls=0,
                  parent_costs=dict(result['parent_costs'], partial_pool=partial['actual_cost'],
                                    four_call_solver_diagnostic={'result': cfg['scc_result'], 'GFN_calls': 4}))
    result['preflight'] = write(out / 'PREFLIGHT.json', {'status': 'passed', 'new_molecular_calls': 0,
        'binding_sha256': result['recovery_binding_sha256'], 'original_complete_cells': 66,
        'saved_component_replacement_cells': 1, 'assembled_point': point_pin,
        'solver_exception_certificate': certificate_pin, 'required_full_cells': 104,
        'remaining_original_slots': pending_keys, 'component_call_ceiling': result['component_call_ceiling'],
        'all_52_prior_solvents_admitted': True, 'original_failure_preserved': failed[0]['attempted_point']})
    return base, ev, candidates, cells, old_cells, whole


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('preflight', 'run'), required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick = time.monotonic()
    usage0 = resource.getrusage(resource.RUSAGE_SELF)
    result = {'protocol': PROTOCOL, 'status': 'execution_failure', 'stage': args.stage,
              'implementation': pin(__file__), 'config': pin(args.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
              'new_DFT_calls': 0, 'optimizer_steps': 0, 'optimization_qualified': False}
    try:
        cfg = json.loads(args.config.read_text())
        require(args.stage != 'run' or cfg.get('recovery_preflight_result') is not None,
                'Run requires the actual successful saved-only preflight')
        base, ev, candidates, cells, old_cells, whole = prepare(cfg, out, result)
        if args.stage == 'preflight':
            result.update(status='complete', preflight_only=True)
        else:
            result['status'] = 'unavailable_backend'
            by_id = {c['candidate_id']: c for c in candidates}
            for key in result['remaining_original_slots']:
                slot = cells[tuple(key)]
                require(slot['status'] == 'unavailable_not_attempted', 'Only original unattempted cells may execute')
                before, count = copy.deepcopy(ev.cost), len(ev.requests)
                slot['status'] = 'unavailable_attempted'
                try:
                    value = ev.evaluate(slot['source_id'], slot['region'], slot['metal'], by_id[slot['candidate_id']]['_x'])
                    cell_from_point(slot, value, False)
                except Exception as exc:
                    slot.update(failure=repr(exc), traceback=traceback.format_exc())
                    break
                finally:
                    slot['cost_delta'] = {k: ev.cost[k] - before.get(k, 0) for k in ev.cost}
                    slot['attempted_point'] = ev.requests[-1] if len(ev.requests) > count else None
            require(ev.cost['composites_attempted'] <= 37 and ev.cost['MACE_calls_attempted'] <= 37
                    and ev.cost['GFN_calls_attempted'] <= 74 and ev.cost['MM_configurations_attempted'] <= 37
                    and ev.cost['MM_group_queries_attempted'] <= 148 and ev.cost['ALPB_calls_attempted'] == 0,
                    'Original solver continuation exceeded the declared component ceiling')
            if summarize(base, candidates, cells, old_cells, whole, result):
                require(len(cells) == 104 and all(r['status'] == 'complete' for r in cells.values()), 'Full 104 cells required')
                result['status'] = 'complete'
        require(binding(cfg) == result['recovery_binding_sha256'], 'Bound recovery inputs changed')
        verify(result['config'])
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    ev = result.pop('_evaluator', None)
    if ev is not None:
        result.update(actual_cost=ev.cost, requests=ev.requests, points=ev.points)
        if ev.native is not None:
            result.update(native_MACE_calls=ev.native.calls, native_MACE_device=ev.native.device,
                          native_MACE_versions=ev.native.versions,
                          peak_CUDA_allocated_bytes=ev.native.torch.cuda.max_memory_allocated())
    usage = resource.getrusage(resource.RUSAGE_SELF)
    result.update(elapsed_seconds=time.monotonic() - tick, peak_host_RSS_KiB=usage.ru_maxrss,
                  observed_parent_CPU_seconds=usage.ru_utime + usage.ru_stime - usage0.ru_utime - usage0.ru_stime)
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds'] * int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
