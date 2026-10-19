#!/usr/bin/env python3
"""Eight fixed LBT states: saved-component preflight and one whole-GFN scout.

This first frozen execution version implements only these two stages. The
declared development/challenge cells remain not_requested; no implicit retry.
Every stage requires Slurm before importing numerical or scientific helpers.
"""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import socket
import sys
import time
import traceback

PROTOCOL = 'lbt_whole_GFN_regional_vacuum_delta_v1'
IDS = tuple(s+'_'+t for s in ('7CCO', '7CCN') for t in ('before_La', 'after_La', 'before_Dy', 'after_Dy'))
SCOUT = '7CCO_before_La'
MISSING = '7CCO_after_Dy'
LIMITS = {'whole_GFN_calls': 8, 'regional_GFN_calls': 1, 'regional_MACE_calls': 1,
          'ALPB_calls': 0, 'MM_calls': 0, 'embedded_GFN_calls': 0, 'DFT_calls': 0,
          'optimizer_steps': 0, 'total_molecular_calls': 10, 'scout_wall_seconds': 600.,
          'development_per_metal_work_error_kcal_mol': 1.}
PROFILE_SHA = '635e111bac31e7c5755bf612c8e148ce73b401c110ff31887fd5e31b590fa85f'
PARAM_SHA = 'f6f05c6c5264cb2f75c54f8ac6838044541c1e493b16f1a4d4ee4fa64de5f545'
INPUT_SHA = 'd94e1cadcfa12e15013cac8e5fbc1679a10c4c566ab64e47ccae4dec09a4cc6a'
EV, HA, BOHR = 23.06054783061903, 627.509474, .529177210903
ARGS = ['core.xyz', '--gfn', '2', '--chrg', '-2', '--uhf', '0', '--etemp', '300',
        '--acc', '0.0001', '--iterations', '500', '--norestart', '--json', '--pop', '--grad']
BIND_KEYS = ('protocol', 'implementation', 'specification', 'design', 'design_review', 'inventory',
             'python', 'python_prefix', 'recovery_implementation', 'recovery_config',
             'whole_reference_helper', 'mapping_replay_helper', 'limits', 'reference_results',
             'directions', 'scope')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as stream:
        sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: '+record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return pin(path)


def verify_tree(value):
    if isinstance(value, dict):
        if set(value) == {'path', 'sha256'}:
            verify(value)
        else:
            for item in value.values():
                verify_tree(item)
    elif isinstance(value, list):
        for item in value:
            verify_tree(item)


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def binding(cfg):
    return {k: cfg[k]['sha256'] if k == 'implementation' else cfg[k] for k in BIND_KEYS}


def xyz_text(symbols, x):
    return str(len(symbols))+'\nFixed declared whole-base regional-delta fixture\n'+''.join(
        '%s %.17g %.17g %.17g\n' % (symbol, *position) for symbol, position in zip(symbols, x))


def shape(value, n, np):
    a = np.asarray(value, dtype=float)
    require(a.shape == (n, 3) and np.isfinite(a).all(), 'Finite full gradient/coordinate shape differs')
    return a


def native_replay(record, task_expected, h, exception=False):
    raw = read(record)
    require(raw['status'] == 'accepted' and raw['attempted'], 'Saved native GFN unavailable')
    task, execution = read(raw['task']), read(raw['execution'])
    verify_tree(raw)
    require(all(task[k] == v for k, v in task_expected.items()), 'Native GFN fixture/state differs')
    require(task['runtime'] == h['runtime'] and task['full_corrected_profile'] == h['pins']['GFN_profile']
            and task['parameter_loading_source'] == h['pins']['parameter_loading_source'], 'Native parameter/runtime lineage differs')
    command = [str(verify(h['runtime']['executable']))]+ARGS
    if exception:
        command += ['--input', 'control.inp']
        require(task['solver_control']['broydamp'] == .05, 'Certified isolated solver differs')
    else:
        require('solver_control' not in task, 'Unexpected solver variant')
    require(task['command'] == execution['command'] == command
            and task['staged_files']['parameters']['sha256'] == PARAM_SHA, 'Saved native command/parameters differ')
    got = h['admit'].admit(verify(raw['task']).parent, task, execution,
                          backend=h['backend'], gradients=h['gradient'], np=h['np'])
    # Original raw quantities remain authoritative. Derived cap sums are a new
    # candidate assembly, not a bitwise assertion about the old hybrid gradient.
    require(all(got[k] == raw[k] for k in ('energy_hartree', 'charges_e', 'scc_iterations', 'xtbout'))
            and got['qm_gradient']['artifact'] == raw['qm_gradient']['artifact']
            and got['qm_gradient']['hartree_per_bohr'] == raw['qm_gradient']['hartree_per_bohr'],
            'Raw native GFN replay differs')
    return raw


def prepare(cfg, out, result):
    require(cfg['protocol'] == PROTOCOL and cfg['limits'] == LIMITS and cfg['inventory']['sha256'] == INPUT_SHA
            and cfg['implementation']['sha256'] == pin(__file__)['sha256'], 'Frozen candidate declaration differs')
    verify_tree(cfg)
    require(pin(sys.executable)['sha256'] == cfg['python']['sha256']
            and str(Path(sys.prefix).resolve()) == cfg['python_prefix'], 'Bound Python environment differs')
    inv = read(cfg['inventory']); pins = inv['helpers']; hc = read(inv['parents']['qualified_hybrid_config'])
    require([r['cell_id'] for r in inv['cells']] == list(IDS)
            and all(hc[k] == v for k, v in pins.items()), 'Exact eight-state/helper inventory differs')
    runtime = pins['GFN_runtime']; profile = read(pins['GFN_profile'])
    require(pins['GFN_profile']['sha256'] == PROFILE_SHA and runtime['parameters']['sha256'] == PARAM_SHA,
            'Corrected full La/Dy parameter profile differs')
    import numpy as np
    sys.path.insert(0, str(verify(pins['common_helper']).parent))
    common = module('base_delta_common', pins['common_helper'])
    gradient = module('base_delta_gradient', pins['gradient_helper'])
    mapper = module('base_delta_mapper', pins['preparation_implementation'])
    mapper.np = np; mapper.cap_jacobians = gradient.cap_jacobians
    admit = module('base_delta_native_admission', pins['admission'])
    backend = module('base_delta_native_backend', pins['GFN_backend'])
    solvent = module('base_delta_solvent_state', pins['ALPB_helper'])
    whole = module('base_delta_whole_reference', cfg['whole_reference_helper'])
    replay = module('base_delta_map_replay', cfg['mapping_replay_helper'])
    require(common.HA_TO_KCAL == HA and common.BOHR_TO_A == BOHR
            and read(pins['GFN_parent_config'])['runtime'] == runtime
            and read(pins['GFN_parent_config'])['tolerances'] == admit.TOL, 'Original units/native gates differ')
    qual = read(pins['GFN_qualification']); qm = read(qual['manifest']); qc = read(qual['config'])
    require(qual['status'] == 'passed' and qual['slurm_job_id'] == '12488'
            and qual['accepted'] == qual['attempted'] == 8 and qm['profile'] == pins['GFN_profile']
            and all(r['passed'] for r in qual['parity']+qual['solvent_component_parity'])
            and qc['staging']['standalone_settings'] == {'accuracy': .0001, 'maxiter': 500,
                'electronic_temperature_K': 300, 'alpb_reference': 'gsolv'}, 'Qualified native lineage unavailable')
    for r in qual['rows']:
        if r['task']['backend'] == 'standalone':
            require(r['status'] == 'accepted' and all(r['task'][k] == v for k, v in runtime.items()),
                    'Qualified runtime differs')
    alpb_parent = read(pins['ALPB_evidence'])
    require(alpb_parent['status'] == 'complete' and alpb_parent['build']['binary'] == pins['ALPB_binary'],
            'Fixed-FF solvent lineage differs')
    verify_tree(alpb_parent['build'])
    mace_parent = read(pins['MACE_evidence']); mcfg = read(pins['MACE_config'])
    require(mace_parent['status'] == 'complete' and mace_parent['model_parameters_unchanged']
            and mace_parent['model'] == pins['model'] == mcfg['model']
            and mace_parent['versions'] == pins['expected_MACE_versions'], 'Native MACE recipe differs')
    verify_tree(mcfg['helper_pins'])
    h = dict(np=np, pins=pins, runtime=runtime, profile=profile, common=common, gradient=gradient,
             mapper=mapper, admit=admit, backend=backend, solvent=solvent, whole=whole, replay=replay)
    # Reuse the already-reviewed saved-only certificate verifier, not the pool
    # evaluator or recovery run path. This makes zero molecular calls.
    recovery = module('base_delta_saved_certificate', cfg['recovery_implementation'])
    rcfg = read(cfg['recovery_config']); partial = read(rcfg['original_partial_result'])
    cert_row = next(r for r in inv['cells'] if r['cell_id'] == '7CCN_after_La')
    cert = read(cert_row['solver_exception_certificate'])
    evidence, _, _ = recovery.certify_scc(rcfg, out, partial, cert['original_failed_slot'])
    require(evidence == cert['saved_scc_qualification'] and cert['same_Hamiltonian']
            and not cert['identical_solver_protocol'] and cert['old_failure_remains_unavailable']
            and evidence['primary'] == cert_row['local_isolated_GFN2'], 'Certified local solver exception differs')
    result['saved_solver_certificate_replayed'] = cert_row['solver_exception_certificate']
    declarations = {s: read(inv['parents'][s+'_inventory']) for s in ('7CCO', '7CCN')}
    supplement = read(inv['parents']['supplemental_after_Dy_ALPB_collection'])
    require(supplement['status'] == 'complete', 'Completed whole-model supplement unavailable')
    cells = {}; manifest = []
    for row in inv['cells']:
        cid = row['cell_id']; sid, metal = row['source_id'], row['metal']; mult = 1 if metal == 'La' else 6
        declared = next(r for r in declarations[sid]['cells'] if r['cell_id'] == cid[len(sid)+1:])
        declared_physical = declared['physical_geometry_parent'] if cid == MISSING else declared['physical']
        require(row['physical'] == declared_physical and row['positions'] == declared['positions']
                and row['metal'] == declared['metal'], 'Frozen reference physical fixture differs')
        if cid == MISSING:
            sr = next(r for r in supplement['rows'] if r['cell_id'] == 'after_Dy')
            require(declared['ALPB'] is None and sr['ALPB'] == row['ALPB']
                    and read(sr['physical_source']) == read(row['physical']),
                    'Exact supplemental after-Dy solvent differs')
        else:
            require(row['ALPB'] == declared['ALPB'], 'Frozen solvent receipt differs')
        atoms = read(row['physical']); ids = [a['id'] for a in atoms]
        x = shape([a['xyz_A'] for a in atoms], 256, np)
        xp = verify(row['positions'])
        saved_x = np.load(xp, allow_pickle=False) if xp.suffix == '.npy' else read(row['positions'])['positions_A']
        require(np.array_equal(x, saved_x) and ids[-1] == 'metal' and len(set(ids)) == 256,
                'Exact physical coordinate/order replay differs')
        symbols = [a['element'] for a in atoms[:-1]]+[metal]
        require(set(symbols[:-1]) <= {'H', 'C', 'N', 'O'} and atoms[-1]['element'] == 'M', 'Whole physical atom inventory differs')
        mapping = read(row['mapping'])
        require(mapping['physical_atoms'] == atoms and len(mapping['core']) == 172
                and len(mapping['caps']) == 2 and len(mapping['field']) == 84,
                'Archived full/regional map differs')
        stask = row.get('source_task', row.get('source_task_template_La'))
        partition = read(stask['partition'])
        regenerated = mapper.update_geometry(partition, x)
        mapping_check = replay.compare(mapping, regenerated)
        charges = read(stask['global_solvent_charges'])
        require(charges['physical_ids'] == ids and charges['total_e'] == -2
                and abs(math.fsum(charges['charges_e'])+2) <= admit.TOL['charge_closure_e'], 'Fixed FF charges differ')
        syq = [metal]+[a['element'] for a in mapping['core'][1:]]
        xq = shape([a['xyz_A'] for a in mapping['core']], 172, np)
        require(mapping['core'][0]['physical_index'] == 255 and stask['charge'] == -2 and stask['region'] == 'baseline',
                'Baseline target/charge differs')
        counts = {}
        for label, sy, n in (('whole', symbols, 256), ('region', syq, 172)):
            electrons = sum(sum(profile['element'][s]['refocc']) for s in sy)+2
            all_e = sum(admit.Z[s] for s in sy)+2
            require(electrons == int(electrons) and electrons > 0 and electrons % 2 == 0
                    and (all_e-mult+1) % 2 == 0, 'GFN/physical state parity differs')
            counts[label] = {'atoms': n, 'charge': -2, 'GFN_uhf': 0, 'GFN_effective_multiplicity': 1,
                             'GFN_valence_electrons': int(electrons), 'physical_all_electrons': all_e,
                             'MACE_physical_multiplicity': mult, 'physical_state_qualified': False}
        require(counts['whole']['GFN_valence_electrons'] == 750
                and counts['region']['GFN_valence_electrons'] == 538, 'Qualified full/regional valence inventory differs')
        cell = dict(row=row, atoms=atoms, ids=ids, x=x, symbols=symbols, mapping=mapping,
                    core_symbols=syq, core_x=xq, state=counts, charges=charges, physical_parent=row['physical'])
        cell['solvent'] = whole.admit_alpb(cell, {'ALPB_binary': pins['ALPB_binary']}, np, common)
        require(cell['solvent']['status'] == 'available', 'Declared solvent unavailable')
        local = None
        if cid != MISSING:
            point = read(row['local_POINT']); ml = read(row['local_MACE']); nt = read(row['task'])
            require(point['status'] == 'complete' and point['mapping'] == row['mapping']
                    and point['task'] == row['task'] and point['component_receipts']['MACE'] == row['local_MACE']
                    and point['component_receipts']['GFN_isolated'] == row['local_isolated_GFN2']
                    and point['point']['positions_A'] == x.tolist() and nt['physical_source']['sha256'] == row['physical']['sha256']
                    and nt['metal'] == metal and nt['charge'] == -2 and nt['multiplicity'] == mult, 'Cached local POINT/task differs')
            nsy, nx = admit.xyz(nt['xyz'])
            require(nsy == syq and np.array_equal(nx, xq), 'Exact regional XYZ differs')
            require(ml['task'] == row['task'] and ml['xyz'] == nt['xyz']
                    and ml['physical_source'] == nt['physical_source'], 'MACE task/coordinate provenance differs')
            solvent.state(ml['input_state_check'], 172, -2, mult, 0)
            require(ml['model'] == pins['model'] and ml['force_convention'] == 'negative vacuum energy derivative'
                    and math.isfinite(ml['energy_eV']), 'Regional MACE recipe/energy differs')
            gm = -shape(ml['forces_eV_A'], 172, np)*EV
            native = native_replay(row['local_isolated_GFN2'],
                {'source_id': sid, 'region': 'baseline', 'metal': metal, 'kind': 'isolated', 'atom_count': 172,
                 'charge': -2, 'uhf': 0, 'multiplicity': 1, 'physical_multiplicity': mult,
                 'valence_electrons': counts['region']['GFN_valence_electrons'], 'pointcharge_count': 0,
                 'field_charge_e': 0.}, h, cid == '7CCN_after_La')
            require(ml['energy_eV'] == point['point']['native_energies']['MACE_eV']
                    and native['energy_hartree'] == point['point']['native_energies']['GFN_isolated_Eh'], 'Cached local scalar replay differs')
            native_task = read(native['task'])
            require(native_task['xyz']['sha256'] == nt['xyz']['sha256']
                    and native_task['physical_source'] == nt['physical_source'], 'Regional native/MACE geometry differs')
            gi = shape(native['qm_gradient']['hartree_per_bohr'], 172, np)*HA/BOHR
            gm_full = mapper.scatter(mapping, gm, np.zeros((84, 3)))
            minus_gi_full = -mapper.scatter(mapping, gi, np.zeros((84, 3)))
            gd = gm_full+minus_gi_full
            shape(gd, 256, np)
            local = {'MACE_eV': ml['energy_eV'], 'GFN_hartree': native['energy_hartree'],
                     'gradient_difference_kcal_mol_A': gd.tolist(), 'MACE': row['local_MACE'],
                     'GFN': row['local_isolated_GFN2'], 'mapping': row['mapping'],
                     'component_gradients_kcal_mol_A': {'MACE_region': gm_full.tolist(),
                                                       'negative_GFN_region': minus_gi_full.tolist()},
                     'solver_exception_certificate': row['solver_exception_certificate'],
                     'old_hybrid_gradient_not_used': True}
        else:
            parent = read(row['map_parent_POINT'])
            require(parent['status'] == 'complete' and parent['mapping'] == row['mapping']
                    and parent['point']['positions_A'] == x.tolist() and parent['state']['metal'] == 'La'
                    and row['local_MACE'] is row['local_isolated_GFN2'] is row['local_POINT'] is None,
                    'Missing after-Dy parent/map differs')
        cell['local'] = local; cells[cid] = cell
        directory = out/cid; directory.mkdir()
        (directory/'whole.xyz').write_text(xyz_text(symbols, x))
        (directory/'region.xyz').write_text(xyz_text(syq, xq))
        manifest.append({'cell_id': cid, 'source_id': sid, 'geometry': row['geometry'], 'metal': metal,
            'physical': row['physical'], 'positions': row['positions'], 'physical_ids': ids,
            'whole_xyz': pin(directory/'whole.xyz'), 'region_xyz': pin(directory/'region.xyz'),
            'mapping': row['mapping'], 'mapping_replay': mapping_check, 'state': counts,
            'ALPB': cell['solvent'], 'local_components': local,
            'whole_GFN_status': 'not_requested', 'regional_status': 'saved_admitted' if local else 'not_requested',
            'new_call_allowlist': ['whole_GFN']+(['regional_GFN', 'regional_MACE'] if cid == MISSING else [])})
    for sid in ('7CCO', '7CCN'):
        for geometry in ('before', 'after'):
            a, b = cells[sid+'_'+geometry+'_La'], cells[sid+'_'+geometry+'_Dy']
            require(a['ids'] == b['ids'] and np.array_equal(a['x'], b['x'])
                    and a['charges'] == b['charges'] and a['mapping'] == b['mapping'], 'Paired metal geometry/map differs')
    require(sum(c['local'] is not None for c in cells.values()) == 7, 'Seven saved local pairs required')
    result['manifest'] = write(out/'MANIFEST.json', {'protocol': PROTOCOL, 'scientific_binding': binding(cfg),
        'limits': LIMITS, 'cells': manifest, 'saved_solver_certificate': result['saved_solver_certificate_replayed'],
        'native_admission': pins['admission'], 'new_molecular_calls': 0, 'model_loaded': False,
        'derivative_qualified': False, 'physical_state_qualified': False})
    result['prepared_states'] = [{'cell_id': r['cell_id'], 'state': r['state'], 'regional_status': r['regional_status'],
                                 'solvent_status': r['ALPB']['status']} for r in manifest]
    return cells, h


def scout(cfg, cells, h, out, result):
    pf = read(cfg['preflight_result'])
    require(pf['protocol'] == PROTOCOL and pf['stage'] == 'preflight' and pf['status'] == 'preflight_passed'
            and pf['scientific_binding'] == binding(cfg) and pf['new_molecular_calls'] == 0
            and not pf['model_loaded'], 'Actual zero-call preflight required')
    old = read(pf['manifest']); current = read(result['manifest'])
    def signature(m):
        return [{k: v for k, v in row.items() if k not in ('whole_xyz', 'region_xyz', 'mapping_replay')}
                | {'whole_xyz_sha256': row['whole_xyz']['sha256'], 'region_xyz_sha256': row['region_xyz']['sha256']}
                for row in m['cells']]
    require(signature(old) == signature(current), 'Actual preflight signatures differ')
    row = next(r for r in current['cells'] if r['cell_id'] == SCOUT)
    d = out/'whole_GFN_scout'; d.mkdir()
    (d/'core.xyz').write_bytes(verify(row['whole_xyz']).read_bytes())
    (d/'param_gfn2-xtb.txt').write_bytes(verify(h['runtime']['parameters']).read_bytes())
    staged = {'core.xyz': pin(d/'core.xyz'), 'parameters': pin(d/'param_gfn2-xtb.txt')}
    require(staged['parameters']['sha256'] == PARAM_SHA and staged['core.xyz']['sha256'] == row['whole_xyz']['sha256'],
            'Native scout staging differs')
    field = write(d/'FIELD.json', {'rows': [], 'external_field': False, 'physical_atom_count': 256})
    command = [str(verify(h['runtime']['executable']))]+ARGS
    task = {'task_id': SCOUT+'_whole_GFN', 'source_id': '7CCO', 'metal': 'La', 'region': 'whole',
            'kind': 'isolated', 'atom_count': 256, 'charge': -2, 'uhf': 0, 'multiplicity': 1,
            'physical_multiplicity': 1, 'valence_electrons': row['state']['whole']['GFN_valence_electrons'],
            'pointcharge_count': 0, 'field_charge_e': 0., 'field': field, 'xyz': staged['core.xyz'],
            'physical_source': row['physical'], 'physical_ids': row['physical_ids'], 'target_index': 255,
            'runtime': h['runtime'], 'full_corrected_profile': h['pins']['GFN_profile'],
            'parameter_loading_source': h['pins']['parameter_loading_source'], 'command': command,
            'staged_files': staged, 'protocol': PROTOCOL, 'scientific_binding': binding(cfg),
            'electronic_state': 'published spin-independent GFN2 valence model; effective UHF0; no DFT ECP',
            'physical_state_qualified': False, 'derivative_qualified': False}
    task_pin = write(d/'TASK.json', task)
    env = h['backend'].clean_env(staged['parameters']['path'])
    for k, v in h['runtime']['environment'].items():
        if v is None: env.pop(k, None)
        else: env[k] = v
    env['XTBPATH'] = str(d.resolve())
    native = {'status': 'unavailable_attempted', 'attempted': True, 'task': task_pin}
    result['new_molecular_calls'] = result['whole_GFN_calls_attempted'] = 1
    result['slots'][0]['status'] = 'unavailable_attempted'
    try:
        execution = h['backend'].call(command, d, env)
        native['execution'] = pin(d/'EXECUTION.json'); verify_tree(staged)
        native.update(h['admit'].admit(d, task, execution, backend=h['backend'], gradients=h['gradient'], np=h['np']), status='accepted')
        shape(native['qm_gradient']['hartree_per_bohr'], 256, h['np'])
        result['scout_cost_gate_passed'] = execution['wall_seconds'] <= LIMITS['scout_wall_seconds']
        result['scout_native_admitted'] = True
        result['scout_wall_seconds'] = execution['wall_seconds']
        result['slots'][0]['status'] = 'accepted' if result['scout_cost_gate_passed'] else 'accepted_excess_cost'
    except Exception as exc:
        native.update(failure=repr(exc), traceback=traceback.format_exc())
        if (d/'EXECUTION.json').is_file(): native['execution'] = pin(d/'EXECUTION.json')
    native['artifacts'] = {p.name: pin(p) for p in d.iterdir() if p.is_file()}
    result['scout_native_result'] = write(d/'RESULT.json', native)
    result['slots'][0]['whole_GFN'] = result['scout_native_result']
    require(native['status'] == 'accepted', 'Whole GFN scout native admission failed; no retry')
    require(result['scout_cost_gate_passed'], 'Whole GFN scout exceeded predeclared 600-second cost limit')
    result['status'] = 'scout_passed'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('preflight', 'run-scout'), required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required before scientific imports')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    tick = time.monotonic()
    result = {'protocol': PROTOCOL, 'stage': args.stage, 'status': 'execution_failure', 'implementation': pin(__file__),
        'config': pin(args.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
        'limits': LIMITS, 'declared_cells': list(IDS), 'slots': [{'cell_id': cid, 'status': 'not_requested'} for cid in IDS],
        'new_molecular_calls': 0, 'whole_GFN_calls_attempted': 0, 'regional_GFN_calls_attempted': 0,
        'regional_MACE_calls_attempted': 0, 'ALPB_calls': 0, 'MM_calls': 0, 'embedded_GFN_calls': 0,
        'DFT_calls': 0, 'optimizer_steps': 0, 'model_loaded': False,
        'development_implemented_in_this_version': False, 'challenge_implemented_in_this_version': False,
        'derivative_qualified': False, 'physical_state_qualified': False, 'accuracy_claim': False}
    try:
        cfg = read(result['config']); result['scientific_binding'] = binding(cfg)
        cells, helpers = prepare(cfg, out, result)
        if args.stage == 'preflight': result['status'] = 'preflight_passed'
        else: scout(cfg, cells, helpers, out, result)
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick, peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  python=sys.version, python_prefix=sys.prefix,
                  allocation={k: os.environ.get(k) for k in ('SLURM_CPUS_ON_NODE', 'SLURM_JOB_CPUS_PER_NODE', 'SLURM_JOB_GPUS')})
    write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json'), 'new_molecular_calls': result['new_molecular_calls']}))
    return 0 if result['status'] in ('preflight_passed', 'scout_passed') else 1


if __name__ == '__main__':
    sys.exit(main())
