#!/usr/bin/env python3
"""Four saved same-coordinate whole/hybrid force comparisons; no evaluations."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
import traceback

PROTOCOL = 'lbt_four_selected_starts_same_coordinate_force_audit_v1'
STATES = (('7CCO', 'La'), ('7CCO', 'Dy'), ('7CCN', 'La'), ('7CCN', 'Dy'))
COMPONENTS = ('MACE_vacuum', 'GFN2_field_difference', 'MM_ledger', 'fixed_FF_charge_ALPB')
RESPONSE_SHA = '05b1c412feebfc72048405c8723d426168872325ce32d956732f5c898c028981'
PERCENTILES = (0, 25, 50, 75, 95, 100)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        return {'path': str(p), 'sha256': hashlib.file_digest(f, 'sha256').hexdigest()}


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


def same(a, b):
    verify(a)
    verify(b)
    return a['sha256'] == b['sha256']


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def inputs(cfg, out, result):
    require(cfg['protocol'] == PROTOCOL and same(cfg['implementation'], pin(__file__)), 'Audit source/protocol differs')
    verify(cfg['specification'])
    prior = read(cfg['response_preflight_result'])
    rcfg = read(prior['config'])
    require(prior['status'] == 'complete' and prior['stage'] == 'preflight' and prior['slurm_job_id'] == '13015'
            and prior['region'] == rcfg['region'] == 'baseline'
            and prior['implementation']['sha256'] == rcfg['implementation']['sha256'] == RESPONSE_SHA,
            'Exact admitted four-state response preflight required')
    verify(prior['implementation'])
    response = module('force_audit_response', rcfg['implementation'])
    d = out / 'hybrid_replay'
    d.mkdir()
    scratch = {}
    loaded = response.prepare(rcfg, d, scratch)
    ev = loaded['ev']
    require(scratch['starts'] == prior['starts'] and scratch['identity_sha256'] == prior['identity_sha256']
            and [(s['source_id'], s['metal']) for s in prior['starts']] == list(STATES), 'Selected start replay differs')
    hcfg = read(loaded['parent']['config'])
    whole = read(hcfg['module2_pool_result'])
    wcfg = read(whole['config'])
    require(whole['status'] == 'complete' and whole['declared_cells'] == 52
            and same(whole['implementation'], wcfg['implementation']), 'Complete whole-model parent required')
    cart = module('force_audit_whole_cartesian', wcfg['cartesian_implementation'])
    ccfg = read(wcfg['cartesian_config'])
    require(same(ccfg['implementation'], wcfg['cartesian_implementation']), 'Whole admission source differs')
    wd = out / 'whole_replay'
    wd.mkdir()
    old = cart.load_inputs(ccfg, wd)
    require(loaded['base'].EV == old['base'].EV and loaded['base'].HA == old['base'].HA
            and loaded['base'].BOHR == old['base'].BOHR, 'Gradient unit constants differ')
    require(ev.mcfg['model'] == old['world']['mcfg']['model'], 'Whole/hybrid MACE parameters differ')
    result.update(response_preflight=cfg['response_preflight_result'], hybrid_parent=rcfg['parent_pool_result'],
                  whole_parent=hcfg['module2_pool_result'], response_implementation=rcfg['implementation'],
                  whole_admission_implementation=wcfg['cartesian_implementation'],
                  selected_start_metadata=prior['starts'], Hamiltonian_identity_sha256=ev.identity_hash,
                  whole_terminal_metadata=whole['terminal_slots'], unchanged_MACE_model=ev.mcfg['model'],
                  conversion_constants={'eV_to_kcal_mol': old['base'].EV, 'Hartree_to_kcal_mol': old['base'].HA,
                                        'Bohr_to_A': old['base'].BOHR})
    by = {(r['candidate_id'], m): p for r in whole['paired_matrix'] for m, p in r['cells'].items()}
    return loaded, cart, old, by


def vector_summary(np, values):
    values = np.asarray(values, dtype=float)
    require(values.ndim == 2 and values.shape[1] == 3 and len(values) and np.isfinite(values).all(), 'Invalid force population')
    norms = np.linalg.norm(values, axis=1)
    return {'atoms': len(values), 'vector_sum': np.sum(values, axis=0).tolist(),
            'sum_atom_norms': float(np.sum(norms)), 'L2': float(np.linalg.norm(values)),
            'RMS_atom_vector': float(np.sqrt(np.mean(norms ** 2))), 'maximum_atom': float(norms.max()),
            'norm_percentiles_linear': {str(p): float(v) for p, v in zip(PERCENTILES, np.percentile(norms, PERCENTILES, method='linear'))}}


def labels(ev, point, source):
    raw = read(point['provenance']['point'])
    mapping = read(raw['mapping'])
    physical = mapping['physical_atoms']
    ids = [a['id'] for a in physical]
    require(ids == [a['id'] for a in source['atoms']] and len(ids) == 256 and ids[-1] == 'metal', 'Physical atom map differs')
    part = ev.partitions[(point['source_id'], 'baseline')][1]
    electronic = {a['physical_index'] for a in mapping['core'] if 'physical_index' in a}
    require(electronic == set(part['qm_source_indices']) | {255} and len(electronic) == 170
            and len(mapping['caps']) == 2, 'Baseline electronic/cap inventory differs')
    roles = {}
    for cap in mapping['caps']:
        for key, role, member in (('retained', 'retained_boundary_anchor', True), ('omitted', 'omitted_boundary_anchor', False)):
            i = cap[key + '_index']
            require(ids[i] == cap[key + '_source_id'] and (i in electronic) is member and i not in roles,
                    'Physical cap-anchor mapping differs')
            roles[i] = role
    omitted = {i for i, role in roles.items() if role == 'omitted_boundary_anchor'}
    field = {a['source_index'] for a in mapping['field']}
    require(len(field) == len(mapping['field'])
            and field == set(range(256)) - electronic - omitted
            and field | omitted == set(range(256)) - electronic,
            'Electronic/shifted-field/omitted-anchor physical coverage differs')
    rows = []
    for i, a in enumerate(physical):
        category = 'target_metal' if i == 255 else roles.get(i, 'electronic_real_other' if i in electronic else 'MM_other')
        rows.append({'index': i, 'id': a['id'], 'residue_id': a['id'].rsplit('/', 1)[0] if i != 255 else 'metal',
                     'resname': a.get('resname', 'metal'), 'element': point['metal'] if i == 255 else a['element'],
                     'source_element': a['element'], 'electronic_real': i in electronic,
                     'boundary_role': roles.get(i), 'category': category})
    return rows, raw['mapping']


def one_state(loaded, cart, old, by, slot, out):
    ev, np = loaded['ev'], loaded['ev'].np
    state = slot['source_id'], slot['metal']
    selected = next(s for s in slot['start_metadata'] if (s['source_id'], s['metal']) == state)
    cell = by[(selected['selected_candidate_id'], state[1])]
    require(cell['status'] == 'complete', 'Whole-model paired state unavailable')
    slot.update(candidate_id=selected['selected_candidate_id'], hybrid_point=selected['point'], whole_cell=cell)
    require('point' in cell, 'Whole cell has no saved full-gradient POINT')
    wp = read(cell['point'])
    require('gradient_kcal_mol_A' in wp and isinstance(read(wp['MACE']).get('forces_eV_A'), list),
            'Whole gradient is absent; no reconstruction from energy-only data')
    evaluator = cart.Evaluator(old, state, out / ('whole_' + '_'.join(state)))
    whole = evaluator.read_point(cell['point'])
    wg = evaluator.gradients(cell['point'])
    hybrid = loaded['qualifier'].restore(ev, selected['point']['point'])
    require(hybrid['provenance'] == selected['point'] and (hybrid['source_id'], hybrid['metal'], hybrid['region']) == (*state, 'baseline')
            and np.array_equal(whole['positions_A'], hybrid['positions_A'])
            and wp['native_energies'] == cell['native_energies'], 'Same-coordinate paired point differs')
    require(same(wp['ALPB'], hybrid['provenance']['component_receipts']['ALPB']), 'Same-coordinate raw ALPB receipt differs')
    source = old['world']['sources'][state[0]]
    atoms, mapping_pin = labels(ev, hybrid, source)
    zeros = np.zeros((256, 3))
    whole_components = {'MACE_vacuum': -wg['MACE_kcal_mol'], 'fixed_FF_charge_ALPB': -wg['fixed_ALPB_kcal_mol'],
                        'GFN2_field_difference': zeros, 'MM_ledger': zeros}
    hybrid_components = {k: -np.asarray(hybrid['component_gradients_kcal_mol_A'][k]) for k in COMPONENTS}
    old_force = -np.asarray(whole['gradient_kcal_mol_A'])
    new_force = -np.asarray(hybrid['gradient_kcal_mol_A'])
    require(np.array_equal(sum(whole_components[k] for k in COMPONENTS), old_force)
            and np.array_equal(sum(hybrid_components[k] for k in COMPONENTS), new_force), 'Raw component force closure differs')
    delta = new_force - old_force
    component_delta = {k: hybrid_components[k] - whole_components[k] for k in COMPONENTS}
    population = {'whole': old_force, 'hybrid': new_force, 'delta': delta}
    components = {'whole': whole_components, 'hybrid': hybrid_components, 'delta': component_delta}
    for a in atoms:
        i = a['index']
        a.update(forces={k: v[i].tolist() for k, v in population.items()},
                 norms={k: float(np.linalg.norm(v[i])) for k, v in population.items()},
                 components={scope: {k: v[i].tolist() for k, v in ps.items()} for scope, ps in components.items()})
    groups = []
    for group in ('category', 'residue_id'):
        for name in dict.fromkeys(a[group] for a in atoms):
            indices = [a['index'] for a in atoms if a[group] == name]
            groups.append({'grouping': group, 'name': name, 'physical_indices': indices,
                           'totals': {k: vector_summary(np, v[indices]) for k, v in population.items()},
                           'components': {scope: {k: vector_summary(np, v[indices]) for k, v in ps.items()} for scope, ps in components.items()}})
    top_atoms = sorted(atoms, key=lambda a: (-a['norms']['delta'], a['index']))[:20]
    residues = [g for g in groups if g['grouping'] == 'residue_id']
    top_residues = sorted(residues, key=lambda g: (-g['totals']['delta']['L2'], g['physical_indices'][0]))[:10]
    full = {'atom_rows': atoms, 'all_groups': groups, 'component_force_units': 'kcal/mol/A', 'mapping': mapping_pin,
            'component_attribution_is_not_residue_energy': True, 'boundary_flags_overlap_electronic_membership': True}
    slot.update(status='available', full_tables=write(out / ('FORCES_' + '_'.join(state) + '.json'), full),
                distributions={k: vector_summary(np, v) for k, v in population.items()},
                component_distributions={scope: {k: vector_summary(np, v) for k, v in ps.items()} for scope, ps in components.items()},
                component_sum_minus_total_delta=vector_summary(np, sum(component_delta.values()) - delta),
                same_raw_ALPB_receipt=wp['ALPB'], ALPB_conversion_order_difference=vector_summary(np, component_delta['fixed_FF_charge_ALPB']),
                whole_full_force=loaded['core'].force_metrics(-old_force), hybrid_full_force=loaded['core'].force_metrics(-new_force),
                top20_atoms=top_atoms, top10_residues=top_residues,
                hybrid_state=hybrid['state'], whole_state=whole['physical_state'])
    require(slot['hybrid_full_force'] == selected['initial_full_force'] and evaluator.new_evaluations == 0,
            'Saved hybrid force replay differs or whole evaluator made a call')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Saved force arithmetic requires Slurm')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick, cpu = time.monotonic(), time.process_time()
    result = {'protocol': PROTOCOL, 'status': 'execution_failure', 'implementation': pin(__file__), 'config': pin(a.config),
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'new_molecular_calls': 0, 'new_geometries': 0,
              'force_convention': 'F = -gradient; delta F = hybrid - whole; kcal/mol/A',
              'declared_states': [list(s) for s in STATES], 'rows': [{'source_id': s, 'metal': m, 'status': 'unavailable_not_admitted'} for s, m in STATES]}
    try:
        cfg = json.loads(a.config.read_text())
        loaded, cart, old, by = inputs(cfg, out, result)
        for slot in result['rows']:
            slot['start_metadata'] = result['selected_start_metadata']
            try:
                one_state(loaded, cart, old, by, slot, out)
            except Exception as exc:
                slot.update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
            finally:
                slot.pop('start_metadata', None)
        require(loaded['ev'].cost['composites_attempted'] == 0 and loaded['ev'].native is None, 'Hybrid replay made a molecular call')
        result.update(status='complete' if all(r['status'] == 'available' for r in result['rows']) else 'unavailable',
                      available_states=sum(r['status'] == 'available' for r in result['rows']))
        verify(result['config'])
    except Exception as exc:
        result.update(failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic() - tick, process_CPU_seconds=time.process_time() - cpu)
    lines = ['# Same-coordinate whole/hybrid force audit', '', 'Status: ' + result['status'] + '.', '',
             'Forces are kcal/mol/Å; delta = hybrid − whole. No coordinate or selection changed.', '',
             '| Source / metal | Whole max force | Hybrid max force | Delta RMS | Largest delta-force atom |',
             '| --- | ---: | ---: | ---: | --- |']
    for r in result['rows']:
        if r['status'] == 'available':
            lines.append('| %s / %s | %.4f | %.4f | %.4f | %s |' % (r['source_id'], r['metal'],
                r['whole_full_force']['maximum_atom_kcal_mol_A'], r['hybrid_full_force']['maximum_atom_kcal_mol_A'],
                r['distributions']['delta']['RMS_atom_vector'], r['top20_atoms'][0]['id']))
        else:
            lines.append('| %s / %s | unavailable | unavailable | unavailable | — |' % (r['source_id'], r['metal']))
    lines += ['', 'All atoms, residues, map-defined categories and component distributions are retained in the JSON tables.',
              'Force differences and vector sums are descriptive; they do not establish which model is accurate.',
              'Full convergence flags remain explicit in RESULT.json. No new forces, optimization or force-selected motion.', '']
    (out / 'REPORT.md').write_text('\n'.join(lines))
    result['report'] = pin(out / 'REPORT.md')
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
