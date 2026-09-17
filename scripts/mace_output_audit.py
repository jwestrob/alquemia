"""Compare real saved full-protein MACE outputs without evaluating a model."""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from affordable_common import InvalidArtifact, read_json, verify, record, write_new
from mace_hybrid import EV_TO_KCAL, xyz


def load_pair(collection):
    c = read_json(collection)
    if c['status'] != 'complete':
        raise InvalidArtifact('complete collection required')
    m = read_json(verify(c['manifest']))
    data = {}
    atoms = None
    for metal in ('La', 'Ca'):
        task = next(t for t in m['tasks'] if t['kind'] == 'full' and
                    t['metal'] == metal and t['variant'] == 'primary')
        source_ref = next(ref for ref in m['source_states']
                          if read_json(verify(ref))['metal'] == metal)
        source = read_json(verify(source_ref))
        physical = source['physical_atoms']
        rows = xyz(verify(task['xyz']))
        positions = np.array([a['xyz_A'] for a in physical])
        expected = [metal if a['element'] == 'M' else a['element'] for a in physical]
        if ([a[0] for a in rows] != expected or
                not np.allclose(np.array([r[1:] for r in rows]), positions, atol=1e-10, rtol=0)):
            raise InvalidArtifact('source-to-physical coordinate/element mismatch')
        if atoms is not None and physical != atoms:
            raise InvalidArtifact('paired physical source atoms differ')
        atoms = physical
        row = c['rows'][task['task_id']]
        if row['status'] != 'computed' or not row['charge_check']:
            raise InvalidArtifact('computed charge-closed endpoint required')
        forces = np.load(verify(row['forces']))
        density = np.load(verify(row['density_coefficients']))
        if (forces.shape != (len(atoms), 3) or density.shape != (len(atoms), 4) or
                not np.isfinite(forces).all() or not np.isfinite(density).all()):
            raise InvalidArtifact('invalid saved force/density array')
        if abs(density[:, 0].sum() - task['charge']) > m['tolerances']['total_charge_e']:
            raise InvalidArtifact('saved density charge closure failed')
        data[metal] = {'row': row, 'forces': forces, 'charges': density[:, 0]}
    return c, atoms, data


def audit(medium_collection, large_collection, plan, output):
    collections = {'medium': medium_collection, 'large': large_collection}
    loaded = {k: load_pair(p) for k, p in collections.items()}
    atoms = loaded['medium'][1]
    if atoms != loaded['large'][1]:
        raise InvalidArtifact('checkpoint physical geometries differ')
    ids = [a['id'] for a in atoms]
    if len(set(ids)) != len(ids) or ids.count('metal') != 1:
        raise InvalidArtifact('nonunique physical mapping or missing metal')
    positions = np.array([a['xyz_A'] for a in atoms])
    distances = np.linalg.norm(positions - positions[ids.index('metal')], axis=1)
    nearest_dist, nearest_idx = cKDTree(positions).query(positions, k=2)
    if np.any(nearest_dist[:, 1] == 0):
        raise InvalidArtifact('coincident physical atoms')
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    result = {'status': 'complete', 'analysis': 'saved_full_output_audit_v1',
              'implementation': record(__file__), 'plan': record(plan),
              'collections': {k: record(p) for k, p in collections.items()},
              'new_model_calls': 0, 'atoms': len(atoms), 'models': {},
              'gradient_definition': 'grad(E_Ca-E_La) = F_La-F_Ca, eV/Angstrom',
              'hybrid_gradient': None, 'relaxation_correction': None,
              'interpretation': 'direct vacuum MACE only; exploratory consumed case; no calibrated class'}
    columns = {'atom_index': np.arange(len(atoms)), 'source_id': ids,
               'element': [a['element'] for a in atoms], 'metal_distance_A': distances}

    def atom_row(i, value):
        return {'atom_index': int(i), 'source_id': ids[i], 'value_eV_A': float(value),
                'metal_distance_A': float(distances[i]),
                'nearest_atom_id': ids[nearest_idx[i, 1]],
                'nearest_atom_distance_A': float(nearest_dist[i, 1])}

    for name, (collection, _, data) in loaded.items():
        la, ca = data['La'], data['Ca']
        gradient = la['forces'] - ca['forces']
        response = ca['charges'] - la['charges']
        norm = np.linalg.norm(gradient, axis=1)
        gp = out / f'{name}_direct_R_gradient_eV_A.npy'
        np.save(gp, gradient)
        components = {key: (ca['row']['energy_components_eV'][key] -
                            la['row']['energy_components_eV'][key]) * EV_TO_KCAL
                      for key in ca['row']['energy_components_eV']}
        raw = (ca['row']['energy_eV'] - la['row']['energy_eV']) * EV_TO_KCAL
        components['remaining_reference_term'] = raw - sum(components.values())
        shells = []
        for lower, upper in ((0, 5), (5, 10), (10, 20), (20, 40), (40, None)):
            mask = (distances >= lower) & (distances < (np.inf if upper is None else upper))
            shells.append({'lower_A': lower, 'upper_A': upper, 'atoms': int(mask.sum()),
                           'Ca_minus_La_charge_sum_e': float(response[mask].sum()),
                           'Ca_minus_La_charge_L1_e': float(np.abs(response[mask]).sum()),
                           'direct_R_gradient_squared_norm_eV2_A2': float(np.square(gradient[mask]).sum())})
        result['models'][name] = {'protocol_id': collection['protocol_id'], 'raw_R_kcal_mol': raw,
                                 'R_components_kcal_mol': components, 'shells': shells,
                                 'gradient': record(gp), 'charge_response_sum_e': float(response.sum()),
                                 'largest_direct_R_gradients': [atom_row(i, norm[i]) for i in np.argsort(-norm)[:10]]}
        for metal, state in data.items():
            columns[f'{name}_{metal}_charge_e'] = state['charges']
            columns[f'{name}_{metal}_force_norm_eV_A'] = np.linalg.norm(state['forces'], axis=1)
        for axis, j in zip('xyz', range(3)):
            columns[f'{name}_direct_R_gradient_{axis}_eV_A'] = gradient[:, j]
    for metal in ('La', 'Ca'):
        delta = loaded['large'][2][metal]['forces'] - loaded['medium'][2][metal]['forces']
        norm = np.linalg.norm(delta, axis=1)
        result[f'{metal}_largest_checkpoint_force_changes'] = [atom_row(i, norm[i]) for i in np.argsort(-norm)[:10]]
    result['large_minus_medium_R_components_kcal_mol'] = {
        k: result['models']['large']['R_components_kcal_mol'][k] -
        result['models']['medium']['R_components_kcal_mol'][k]
        for k in result['models']['large']['R_components_kcal_mol']}
    result['large_minus_medium_R_kcal_mol'] = result['models']['large']['raw_R_kcal_mol'] - result['models']['medium']['raw_R_kcal_mol']
    table = out / 'physical_atom_outputs.tsv'
    with table.open('x') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))
    result['all_atom_table'] = record(table)
    write_new(out / 'result.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('medium-collection', 'large-collection', 'plan', 'output'):
        p.add_argument('--' + key, required=True)
    a = p.parse_args()
    r = audit(a.medium_collection, a.large_collection, a.plan, a.output)
    print(f"{r['atoms']} physical atoms; no model calls; results in {a.output}")


if __name__ == '__main__':
    main()
