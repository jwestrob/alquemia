#!/usr/bin/env python3
"""Saved-only, four-state local force diagnostic; no molecular calls."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
import traceback

IDS = ('7CCO_La', '7CCN_La', '7CCO_Dy', '7CCN_Dy')
HA = 627.509474
BOHR = 0.529177210903
EV = 23.06054783061903
GROUPS = ('all_atoms', 'target_metal', 'other_real_atoms', 'synthetic_caps')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    p = Path(record['path'])
    require(pin(p)['sha256'] == record['sha256'], 'Changed pin: ' + str(p))
    return p


def read(record):
    return json.loads(verify(record).read_text())


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def inputs(cfg):
    verify(cfg['plan']); verify(cfg['implementation'])
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Comparator source differs')
    parent = read(cfg['energy_comparison_result']); ecfg = read(parent['config'])
    verify(parent['implementation']); verify(ecfg['implementation'])
    require(parent['status'] == 'complete' and parent['slurm_job_id'] == '12936'
            and parent['implementation']['sha256'] == ecfg['implementation']['sha256'],
            'Completed energy comparison lineage differs')
    energy = module('local_force_energy_admission', ecfg['implementation'])
    replay = {}; energy.compare(ecfg, replay)
    require(replay['status'] == 'complete' and all(parent[k] == v for k, v in replay.items()),
            'Saved MACE/GFN/PBE0 energy admission does not replay')
    refcfg = read(ecfg['reference_config'])
    reference = module('local_force_reference_admission', ecfg['reference_adapter'])
    _, model, _, _, _, _, _ = reference.inputs(refcfg)
    require(model.HA_TO_KCAL == HA and model.BOHR_TO_A == BOHR, 'Native conversion constants differ')
    forces = read(cfg['MACE_force_result']); fcfg = read(forces['config'])
    verify(forces['implementation']); verify(fcfg['implementation']); verify(fcfg['plan'])
    require(forces['status'] == 'complete' and forces['slurm_job_id'] == '12937'
            and forces['implementation']['sha256'] == fcfg['implementation']['sha256']
            and forces['implementation']['sha256'] ==
            'cba062fa6a58db6dc19d52b628bfec9cb4b8df5c95bc661d91818334c077731e',
            'Exact recovered four-force implementation/result required')
    old = read(ecfg['MACE_result']); mcfg = read(old['config'])
    require(forces['energy_reference'] == fcfg['energy_reference'] == ecfg['MACE_result']
            and fcfg['preparation'] == refcfg['preparation_result'] == mcfg['preparation'],
            'Force/reference preparation differs')
    precedent = read(fcfg['force_precedent']); pcfg = read(precedent['config'])
    verify(pcfg['implementation'])
    require(precedent['status'] == 'complete' and precedent['slurm_job_id'] == '12500'
            and pcfg['energy_replay_tolerance_eV'] == fcfg['energy_replay_tolerance_eV']
            == forces['energy_replay_tolerance_eV'] == 1e-7, 'Native force parity rule differs')
    for key in ('model', 'helper_pins', 'expected_versions', 'python', 'preparation', 'native_precedent', 'native_reuse'):
        require(fcfg[key] == pcfg[key], 'Force precedent recipe differs: ' + key)
    for p in fcfg['helper_pins'].values():
        verify(p)
    verify(fcfg['model']['checkpoint'])
    require(forces['model'] == fcfg['model'] == old['model'] == mcfg['model']
            and forces['versions'] == old['versions'] == fcfg['expected_versions']
            and forces['analytic_forces_requested'] and forces['model_parameters_unchanged']
            and forces['declared_cells'] == forces['new_model_calls_attempted'] == 4,
            'Force model/runtime/call admission differs')
    require(fcfg['tasks'] == [t for t in mcfg['tasks'] if t['representation'] == 'matched_core'],
            'Exact original force task set/order differs')
    frows = {r['source_id'] + '_' + r['metal']: r for r in forces['rows']}
    oldrows = {r['source_id'] + '_' + r['metal']: r for r in old['rows'] if r['representation'] == 'matched_core'}
    require(len(forces['rows']) == 4 and set(frows) == set(oldrows) == set(IDS), 'Force denominator differs')
    for task in fcfg['tasks']:
        tid = task['source_id'] + '_' + task['metal']; row = frows[tid]; prior = oldrows[tid]
        require(all(row.get(k) == v for k, v in task.items()) and
                {k: v for k, v in task.items() if k != 'cache_key'} == refcfg['tasks'][tid], 'Force task/source differs')
        require(row['status'] == 'complete' and row['call_attempted'] and row['energy_replay_passed']
                and row['input_state_check'] == prior['input_state_check']
                and row['original_energy_eV'] == prior['energy_eV']
                and row['energy_replay_delta_eV'] == row['energy_eV'] - prior['energy_eV']
                and abs(row['energy_replay_delta_eV']) <= 1e-7
                and row['force_convention'] == 'negative derivative of the same native vacuum energy'
                and row['target_index'] == row['metal_index'] == 0, 'Force state/energy/sign admission differs')
        require(read(pin(verify(cfg['MACE_force_result']).parent / (row['task_id'] + '.json'))) == row,
                'Combined/individual force receipt differs')
    return parent, ecfg, refcfg, model, forces, frows


def atom_map(task, model):
    physical = read(task['physical_source']); mapping = read(task['source_cap_map'])
    verify(mapping['legacy_builder'])
    atoms = model.xyz(verify(task['xyz']))
    require(len(physical) == 256 and len(atoms) == 48 and mapping['shared_real_atoms_exactly_equal'],
            'Prepared physical inventory differs')
    require([r['core_index'] for r in mapping['mapping']] == list(range(48)), 'Core map order differs')
    groups = {'all_atoms': list(range(48)), 'target_metal': [], 'other_real_atoms': [], 'synthetic_caps': []}
    records = []
    real_ids = []
    for row, atom in zip(mapping['mapping'], atoms):
        i = row['core_index']; symbol = atom[0]; xyz = list(atom[1:])
        if row['kind'] == 'source':
            source = physical[row['physical_index']]
            require(source['id'] == row['physical_id'] and source['xyz_A'] == xyz
                    and (task['metal'] if source['element'] == 'M' else source['element']) == symbol,
                    'Mapped real atom identity/coordinate differs')
            if i == 0:
                require(source['id'] == 'metal' and row['physical_index'] == 255
                        and source['element'] == 'M', 'Target map differs')
                group = 'target_metal'
            else:
                require(source['element'] != 'M', 'Extra target in real map')
                group = 'other_real_atoms'
            real_ids.append(source['id'])
        else:
            require(row['kind'] == 'synthetic_cap' and row['whole_physical_atom'] is False
                    and row['xyz_A'] == xyz and row['definition']['element'] == symbol, 'Cap map differs')
            group = 'synthetic_caps'
        groups[group].append(i)
        records.append(dict(row, element=symbol, xyz_A=xyz, diagnostic_group=group))
    require(len(real_ids) == len(set(real_ids)) and all(groups[g] for g in GROUPS), 'Invalid group partition')
    require(sorted(groups['target_metal'] + groups['other_real_atoms'] + groups['synthetic_caps']) == list(range(48)),
            'Atom partition incomplete')
    return records, groups


def absent_gfn(record):
    saved = read(record); task = saved['task']; native = read(saved['standalone_json'])
    no_derivative = lambda name: not any(token in name.lower() for token in ('gradient', 'forces', 'force_e', 'engrad'))
    require('--grad' not in task['command'] and all(no_derivative(k) for k in saved)
            and all(no_derivative(k) for k in native) and all(no_derivative(k) for k in saved['artifacts']),
            'Unexpected GFN derivative payload requires explicit admission')
    files = sorted(p.name for p in Path(task['directory']).iterdir()
                   if any(token in p.name.lower() for token in ('grad', 'force')))
    require(files == [], 'Unregistered GFN derivative artifact requires explicit admission')
    return {'status': 'unavailable', 'reason': 'Exact archived local calls are energy-only; no saved gradient or force array',
            'result': record, 'task': task['task_record'], 'native_json': saved['standalone_json'],
            'command': task['command'], 'artifact_names': sorted(saved['artifacts']),
            'native_json_fields': list(native), 'derivative_files': files, 'new_calls': 0}


def compare(cfg, result):
    parent, ecfg, refcfg, model, forces, frows = inputs(cfg)
    np = model.np
    result['parent_admission_replayed'] = True
    result['MACE_energy_replay_deltas_eV'] = {tid: frows[tid]['energy_replay_delta_eV'] for tid in IDS}
    energy_rows = {r['task_id']: r for r in parent['matrix']}
    for output in result['matrix']:
        tid = output['task_id']; task = refcfg['tasks'][tid]; mr = frows[tid]
        try:
            maps, groups = atom_map(task, model)
            dr = read(ecfg['reference_results'][tid])['rows'][tid]
            require(dr['status'] == 'available' and dr['quantity'] == 'gradient_not_force'
                    and dr['engrad'] == energy_rows[tid]['PBE0']['engrad']
                    and dr['energy_hartree'] == energy_rows[tid]['PBE0']['energy_hartree'], 'Native gradient parent differs')
            native = model.read_engrad(verify(dr['engrad']))
            gradient = np.asarray(dr['gradient_kcal_mol_per_A'], dtype=float)
            fm = np.asarray(mr['forces_eV_A'], dtype=float)
            require(gradient.shape == fm.shape == (48, 3) and np.isfinite(gradient).all() and np.isfinite(fm).all(),
                    'Force/gradient shape or finiteness differs')
            require(np.array_equal(native['gradient_Ha_per_bohr'] * HA / BOHR, gradient),
                    'Native gradient unit conversion does not replay')
            require(mr['target_force_eV_A'] == mr['forces_eV_A'][0], 'Target force differs from full array')
            fd = -gradient / EV; difference = fm - fd
            stats = {}
            for name, indices in groups.items():
                stats[name] = {'atom_count': len(indices), 'indices': indices,
                    'MACE_mean_force_eV_A': fm[indices].mean(axis=0).tolist(),
                    'PBE0_mean_force_eV_A': fd[indices].mean(axis=0).tolist(),
                    'mean_MACE_minus_PBE0_eV_A': difference[indices].mean(axis=0).tolist(),
                    'RMS_atom_vector_difference_eV_A': float(np.sqrt(np.sum(difference[indices] ** 2) / len(indices)))}
            output.update(status='available', source_id=task['source_id'], metal=task['metal'],
                xyz=task['xyz'], physical_source=task['physical_source'], source_cap_map=task['source_cap_map'],
                physical_multiplicity=task['multiplicity'], charge=task['charge'], atom_map=maps,
                MACE_forces_eV_A=fm.tolist(), PBE0_gradient_Hartree_bohr=native['gradient_Ha_per_bohr'].tolist(),
                PBE0_gradient_kcal_mol_A=gradient.tolist(), PBE0_forces_eV_A=fd.tolist(),
                MACE_minus_PBE0_eV_A=difference.tolist(), groups=stats,
                target={'index': 0, 'MACE_force_eV_A': fm[0].tolist(), 'PBE0_force_eV_A': fd[0].tolist(),
                        'MACE_minus_PBE0_eV_A': difference[0].tolist()},
                PBE0_result=ecfg['reference_results'][tid], PBE0_engrad=dr['engrad'],
                PBE0_electronic_state=dr['electronic_state'], PBE0_SCF_evidence=dr['SCF_evidence'],
                GFN2=absent_gfn(energy_rows[tid]['GFN2_result']), force_consistency_qualified=False)
        except Exception as exc:
            output.update(status='unavailable', reason=repr(exc), traceback=traceback.format_exc())
    result.update(status='complete' if all(r['status'] == 'available' for r in result['matrix']) else 'unavailable',
                  available_MACE_PBE0_cells=sum(r['status'] == 'available' for r in result['matrix']),
                  available_GFN_force_cells=0,
                  source_execution_costs={'MACE_result': cfg['MACE_force_result'],
                      'MACE_new_calls_in_parent': forces['new_model_calls_attempted'],
                      'MACE_elapsed_seconds_in_parent': forces['elapsed_seconds'],
                      'PBE0_results': ecfg['reference_results']})


def report(result):
    lines = ['# Local matched-core force diagnostic', '', 'Status: ' + result['status'] + '.', '',
             'Fixed Cartesian target vectors in eV/Å; signed difference is MACE minus PBE0.', '',
             '| State | MACE target (x, y, z) | PBE0 target (x, y, z) | Difference (x, y, z) |',
             '|---|---|---|---|']
    vector = lambda values: ', '.join(f'{v:.6f}' for v in values)
    for row in result['matrix']:
        if row['status'] != 'available':
            lines.append(f"| {row['task_id']} | unavailable | unavailable | unavailable |")
        else:
            t = row['target']
            lines.append('| ' + row['task_id'] + ' | ' + ' | '.join(vector(t[k]) for k in
                         ('MACE_force_eV_A', 'PBE0_force_eV_A', 'MACE_minus_PBE0_eV_A')) + ' |')
    lines += ['', 'RMS atom-vector differences (eV/Å), with no acceptance threshold:', '',
              '| State | All atoms | Target | Other real atoms | Synthetic caps |', '|---|---:|---:|---:|---:|']
    for row in result['matrix']:
        vals = [f"{row['groups'][g]['RMS_atom_vector_difference_eV_A']:.6f}" for g in GROUPS] if row['status'] == 'available' else ['unavailable'] * 4
        lines.append('| ' + row['task_id'] + ' | ' + ' | '.join(vals) + ' |')
    lines += ['', 'All 48 atom vectors, maps, native gradients, group means and SCF evidence are retained in RESULT.json.', '',
              'GFN force arrays are unavailable from these exact energy-only core receipts. No new calls or substitute geometries are used.', '',
              'PBE0 force consistency and electronic stability remain unqualified. Its frozen-f effective singlet differs from physical Dy sextet metadata used by MACE; spin–orbit accuracy is unqualified. These raw cluster forces, including synthetic caps, are descriptive and do not establish whole-peptide physical forces, affinity or transferable preference accuracy.']
    return '\n'.join(lines) + '\n'


def run(args):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required for saved numerical comparison')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) == 1, 'One CPU required')
    cfg = json.loads(Path(args.config).read_text())
    output = Path(args.output).resolve(); output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); cpu = time.process_time()
    result = {'status': 'execution_failure', 'protocol_id': 'lbt_local_force_comparison_v1',
              'config': pin(args.config), 'implementation': pin(__file__), 'plan': cfg['plan'],
              'energy_comparison_result': cfg['energy_comparison_result'], 'MACE_force_result': cfg['MACE_force_result'],
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'declared_cells': 4,
              'matrix': [{'task_id': tid, 'status': 'unavailable'} for tid in IDS],
              'conversions': {'Hartree_to_kcal_mol': HA, 'bohr_to_angstrom': BOHR, 'eV_to_kcal_mol': EV,
                              'force_definition': 'negative gradient'},
              'new_molecular_calls': 0, 'new_geometries': 0, 'force_accuracy_threshold': None,
              'force_consistency_qualified': False, 'physical_Dy_state_qualified': False}
    try:
        compare(cfg, result)
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-started, process_CPU_seconds=time.process_time()-cpu)
    write(output / 'RESULT.json', result); (output / 'REPORT.md').write_text(report(result))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True); parser.add_argument('--output', required=True)
    raise SystemExit(run(parser.parse_args()))
