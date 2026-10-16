#!/usr/bin/env python3
"""One-CPU saved-only diagnostics of four completed intact-LBT3 searches.

CLI: python SAVED_RESPONSE_DIAGNOSTIC.py --config CONFIG --output FRESH_DIRECTORY
Config pins: implementation, plan, scout_result, remaining_result. No engines,
optimization, geometry generation, repairs, or new scientific acceptance rules.
Top-12 atom/mode lists are display summaries; full bond/donor tables are retained.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import time

PINS = {}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(record):
    actual = pin(record['path'])
    require(actual['sha256'] == record['sha256'], 'Changed input: ' + record['path'])
    PINS[actual['path']] = actual
    return Path(actual['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def point_arrays(record, source, metal, ids, np):
    point = read(record)
    require(point['status'] == 'complete' and point['source_id'] == source and point['metal'] == metal, 'Point identity/status differs')
    recipe = point['recipe']
    require(recipe['physical_ids'] == ids and recipe['physical_charge'] == -2
            and recipe['physical_multiplicity'] == (1 if metal == 'La' else 6)
            and recipe['GFN2_effective_multiplicity'] == 1, 'Point physical/effective state differs')
    x = np.load(verify(point['positions']), allow_pickle=False)
    g = np.load(verify(point['gradient_kcal_mol_A']), allow_pickle=False)
    require(x.shape == g.shape == (256, 3) and np.isfinite(x).all() and np.isfinite(g).all(), 'Invalid saved arrays')
    physical = read(point['physical_source'])
    require([a['id'] for a in physical] == ids and np.array_equal(x, [a['xyz_A'] for a in physical]), 'Point physical mapping differs')
    lines = verify(point['xyz']).read_text().splitlines(); rows = [line.split() for line in lines[2:] if line.strip()]
    require(int(lines[0]) == len(rows) == 256 and all(len(r) == 4 for r in rows), 'Point XYZ layout differs')
    require([r[0] for r in rows] == recipe['symbols']
            == [metal if a['element'] == 'M' else a['element'] for a in physical]
            and np.array_equal(x, [[float(v) for v in r[1:]] for r in rows]), 'Point XYZ differs')
    return point, x, g


def one_search(owner, profile_rows, settings, np):
    search = read(owner['search']); source, metal = owner['source_id'], owner['metal']
    require(owner['status'] == 'complete' and search['status'] == 'finite_proposal_available'
            and search['candidate_admitted'] and search['final_geometry']['pass']
            and (search['source_id'], search['metal']) == (source, metal), 'Search unavailable or mismatched')
    profile = profile_rows[source]; prepared = read(profile['source_preparation']); modes = read(profile['modes'])
    require(prepared['source_id'] == modes['source_id'] == source and prepared['atoms'] == 256
            and prepared['charge'] == -2 and prepared['waters'] == prepared['spectators'] == prepared['caps'] == [], 'Source state differs')
    atoms = read(prepared['physical_atoms']); ids = [a['id'] for a in atoms]; by_id = {name: i for i, name in enumerate(ids)}
    require(len(ids) == len(by_id) == 256 and ids[-1] == 'metal' and ids == modes['physical_source_ids'], 'Source atom ordering differs')
    origin_pin = search['trace'][0]['point']; endpoint_pin = search['candidate']
    require(origin_pin in owner['points'] and endpoint_pin in owner['points']
            and search['trace'][-1]['point'] == endpoint_pin, 'Terminal/origin point is not a recorded accepted search point')
    origin, x0, g0 = point_arrays(origin_pin, source, metal, ids, np)
    endpoint, x1, g1 = point_arrays(endpoint_pin, source, metal, ids, np)
    require(origin['xyz'] == prepared['endpoints']['whole'][metal]['xyz'] == search['authoritative_origin']
            and origin['recipe'] == endpoint['recipe'] and np.array_equal(x0, [a['xyz_A'] for a in atoms]), 'Origin or endpoint recipe differs')
    delta = np.linalg.norm(x1-x0, axis=1); force = [-g0, -g1]
    norms = [np.linalg.norm(f, axis=1) for f in force]
    def label(i):
        return {'index': int(i), 'id': ids[i], 'residue': ids[i].rsplit('/', 1)[0],
                'resname': atoms[i].get('resname'), 'element': metal if i == 255 else atoms[i]['element']}
    bonds = read(prepared['covalent_bonds']); bond_rows = []; adjacency = [[] for _ in ids]
    for b in bonds:
        i, j = b['native_particle_indices']
        require([ids[i], ids[j]] == [b['atom_a_id'], b['atom_b_id']], 'Covalent bond atom mapping differs')
        lengths = [float(np.linalg.norm(x[j]-x[i])) for x in (x0, x1)]; equilibrium = b['forcefield_equilibrium_A']
        row = {'atom_a': label(i), 'atom_b': label(j), 'involves_H': 'H' in (atoms[i]['element'], atoms[j]['element']),
               'source_A': lengths[0], 'endpoint_A': lengths[1], 'change_A': lengths[1]-lengths[0],
               'forcefield_reference_A': equilibrium, 'source_minus_reference_A': lengths[0]-equilibrium,
               'endpoint_minus_reference_A': lengths[1]-equilibrium}
        bond_rows.append(row)
        for a, neighbor in ((i, j), (j, i)):
            adjacency[a].append({'neighbor': label(neighbor), **{k: v for k, v in row.items() if k not in ('atom_a', 'atom_b')},
                'relative_force_formula': '(F_neighbor-F_atom) dot unit(atom->neighbor); descriptor, not a bond-energy derivative',
                'source_relative_force_kcal_mol_A': float((force[0][neighbor]-force[0][a]) @ ((x0[neighbor]-x0[a])/lengths[0])),
                'endpoint_relative_force_kcal_mol_A': float((force[1][neighbor]-force[1][a]) @ ((x1[neighbor]-x1[a])/lengths[1]))})
    def force_row(i):
        return {**label(i), 'source_force_kcal_mol_A': force[0][i].tolist(), 'endpoint_force_kcal_mol_A': force[1][i].tolist(),
                'source_norm_kcal_mol_A': float(norms[0][i]), 'endpoint_norm_kcal_mol_A': float(norms[1][i]), 'bonded_neighbors': adjacency[i]}
    order = sorted(range(256), key=lambda i: (-float(delta[i]), i)); residues = {}
    for i in range(256):
        residues.setdefault(ids[i].rsplit('/', 1)[0], []).append(i)
    residue_rows = [{'residue': r, 'resname': atoms[ii[0]].get('resname'), 'max_displacement_A': float(delta[ii].max()),
                     'RMS_displacement_A': float(np.sqrt(np.mean(delta[ii]**2))), 'limiting_atom': label(max(ii, key=lambda i: delta[i]))}
                    for r, ii in residues.items()]
    contacts_pin = prepared['endpoints']['matched_core'][metal]['source_cap_map']; contacts = read(contacts_pin)['contacts']; donors = []
    for c in contacts:
        key = f"{c['chain']}/{c['resnum']}//{c['atom']}"; i = by_id[key]
        distances = [float(np.linalg.norm(x[i]-x[255])) for x in (x0, x1)]
        donors.append({**label(i), 'declared_fragment': c['fragment_id'], 'source_A': distances[0],
                       'endpoint_A': distances[1], 'change_A': distances[1]-distances[0]})
    diagnostics = {}
    for phase, saved, g in (('source', search['trace'][0]['force_diagnostic'], g0), ('endpoint', search['final_force_diagnostic'], g1)):
        loads = saved['mode_force_loads_kcal_mol_A']
        require(len(loads) == 64 and [r['name'] for r in loads] == modes['mode_order']
                and [r['mode_index'] for r in loads] == list(range(64))
                and {r['mode_index'] for r in loads if r['active']} == set(search['active_indices'])
                and np.isclose(np.linalg.norm(g), saved['full_gradient_L2_kcal_mol_A'], atol=1e-8, rtol=1e-12), 'Saved force diagnostic differs')
        omitted = sorted((r for r in loads if not r['active']), key=lambda r: (-abs(r['force_load']), r['mode_index']))
        diagnostics[phase] = {'saved_active_and_residual': {k: v for k, v in saved.items() if k != 'mode_force_loads_kcal_mol_A'},
            'outside_to_active_L2_ratio': saved['outside_active_span_L2_kcal_mol_A']/saved['active_span_L2_kcal_mol_A'] if saved['active_span_L2_kcal_mol_A'] else None,
            'all64_mode_loads': loads, 'top_omitted_modes': omitted[:12],
            'top_omitted_backbone_modes': [r for r in omitted if r['family'] in ('phi', 'psi')][:12]}
    trials = [{'accepted_step': row['accepted_step'], 'multiplier': t['multiplier'], 'reason': t.get('reason'),
               'molecular_evaluation_requested': t['molecular_evaluation_requested'], 'geometry': t['geometry']}
              for row in search['trace'] for t in row.get('trials', []) if not t['accepted']]
    return {'source_id': source, 'metal': metal, 'search': owner['search'], 'origin_point': origin_pin, 'endpoint_point': endpoint_pin,
        'preparation': profile['source_preparation'], 'modes': profile['modes'], 'declared_donor_map': contacts_pin,
        'termination': search['termination'], 'accepted_steps': search['accepted_steps'], 'work_from_origin': search['work_from_origin'],
        'saved_final_geometry': search['final_geometry'], 'rejected_trials': trials,
        'extent': {'bound_A': settings['all_atom_displacement_A'], 'maximum_A': float(delta.max()),
            'target_A': float(delta[255]), 'closest_to_extent_bound_atoms': [{**label(i), 'displacement_A': float(delta[i]),
                'headroom_A': settings['all_atom_displacement_A']-float(delta[i])} for i in order[:12]],
            'residues': sorted(residue_rows, key=lambda r: (-r['max_displacement_A'], r['residue']))},
        'donor_selection': 'Every atom in the pinned existing matched-core contact map, including both named carboxylate oxygens; no distance cutoff or bonding inference',
        'donor_distances': donors, 'force_projection': diagnostics, 'bonds': bond_rows,
        'top_full_force_atoms': {phase: [force_row(i) for i in sorted(range(256), key=lambda i: (-float(norms[k][i]), i))[:12]]
                                 for k, phase in enumerate(('source', 'endpoint'))},
        'top_endpoint_H_force_atoms': [force_row(i) for i in sorted((i for i, a in enumerate(atoms) if a['element'] == 'H'), key=lambda i: (-float(norms[1][i]), i))[:12]],
        'interpretation': 'Extent/contact rejection, omitted-mode loads and bond/H force descriptors are separate observations; no causal attribution or repair prescription.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) == 1, 'Exactly one allocated CPU required')
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[key] = '1'
    cfg_pin = pin(args.config); cfg = read(cfg_pin); verify(cfg['implementation']); verify(cfg['plan'])
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation differs')
    scout, remaining = read(cfg['scout_result']), read(cfg['remaining_result'])
    require(scout['status'] == remaining['status'] == 'complete' and scout['slurm_job_id'] == '12533'
            and remaining['slurm_job_id'] == '12538' and remaining['scout_result'] == cfg['scout_result']
            and scout['plan'] == remaining['plan'] == cfg['plan'] and scout['settings'] == remaining['settings']
            and scout['profile_result'] == remaining['profile_result'], 'Completed search lineage differs')
    for result in (scout, remaining):
        verify(result['config']); verify(result['implementation'])
    profile = read(scout['profile_result']); modes = read(profile['modes_result']); preparation = read(modes['preparation'])
    require(profile['status'] == modes['status'] == preparation['status'] == 'complete', 'Source lineage incomplete')
    profile_rows = {r['source_id']: r for r in profile['rows']}
    owners = [scout]
    for embedded in remaining['rows']:
        row = read(embedded['result'])
        require(row == {k: v for k, v in embedded.items() if k != 'result'}, 'Remaining row receipt differs')
        owners.append(row)
    require(len(owners) == 4 and {(r['source_id'], r['metal']) for r in owners} == {(s, m) for s in ('7CCO', '7CCN') for m in ('La', 'Dy')}, 'Exact four-search denominator required')
    import numpy as np
    start = time.monotonic(); out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    result = {'schema': 'lbt3_saved_response_diagnostic.v1', 'status': 'incomplete', 'config': cfg_pin,
              'implementation': pin(__file__), 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
              'new_molecular_calls': 0, 'geometry_generation': False, 'numpy_version': np.__version__, 'rows': []}
    result['rows'] = [one_search(owner, profile_rows, scout['settings'], np) for owner in owners]
    for record in list(PINS.values()):
        verify(record)
    result.update(status='complete', declared_searches=4, source_pins=list(PINS.values()), elapsed_seconds=time.monotonic()-start)
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds']; write(out/'RESULT.json', result)
    lines = ['# Saved intact-peptide response diagnostics', '', 'Four searches; zero new molecular calls. Forces are negative saved composite gradients.', '']
    for row in result['rows']:
        top = row['extent']['closest_to_extent_bound_atoms'][0]; force = row['top_full_force_atoms']['endpoint'][0]
        lines += [f"## {row['source_id']} / {row['metal']}", '', f"Recorded stop: {row['termination']}; accepted steps: {row['accepted_steps']}.",
                  f"Largest source displacement: {top['id']} ({top['displacement_A']:.6f} Å); target {row['extent']['target_A']:.6f} Å.",
                  f"Largest endpoint full force: {force['id']} ({force['endpoint_norm_kcal_mol_A']:.6f} kcal/mol/Å).", '']
    lines += ['Full atom/residue extent, declared donor distances, all64 mode loads, rejected guards, and source/endpoint bond details are in RESULT.json.',
              'Omitted loads are overlapping coordinate diagnostics. Force projections and bond reference deviations do not identify a causal defect or prescribe repair.']
    with (out/'REPORT.md').open('x') as stream:
        stream.write('\n'.join(lines)+'\n')
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))


if __name__ == '__main__':
    main()
