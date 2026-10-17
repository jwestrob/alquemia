#!/usr/bin/env python3
"""Saved Mex mechanics and deterministic two-metal segment selection."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import time
import traceback

PROTOCOL = 'mex_paired_saved_motion_segment_selection_v1'
METALS = ('La', 'Dy')
GROUPS = ('retained_bonded', 'MM_LJ', 'MM_Coulomb', 'QM_MM_LJ')
COMPONENTS = ('MACE_vacuum', 'GFN2_field_difference', 'MM_ledger', 'fixed_FF_charge_ALPB')


def require(test, message):
    if not test:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: '+record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')
    return pin(path)


def equal_pin(a, b):
    verify(a); verify(b)
    return a['sha256'] == b['sha256']


def close(a, b, tolerance, label):
    require(math.isfinite(a) and math.isfinite(b) and abs(a-b) <= tolerance, label)


def array(values, size=2049):
    a = np.asarray(values, dtype=float)
    require(a.shape == (size, 3) and np.isfinite(a).all(), 'Invalid full physical array')
    return a


def residue(atom_id):
    chain, number, insertion, _ = atom_id.split('/')
    return chain, int(number), insertion


def residue_label(key):
    return '/'.join((key[0], str(key[1]), key[2]))


def point(provenance, metal):
    p = read(provenance['point']); request = read(provenance['cartesian_request'])
    require(p['status'] == 'complete' and p['metal'] == metal and request['state']['metal'] == metal
        and request['state'] == provenance['state'] and request['identity_sha256'] == provenance['identity_sha256']
        and request['cache_key'] == provenance['cache_key'] and equal_pin(request['point'], provenance['point'])
        and request['native_energies'] == p['native_energies'], 'Point/state/request admission differs')
    physical = read(p['physical']); ids = [a['id'] for a in physical]
    require(len(ids) == len(set(ids)) == 2049 and ids == request['source_ids'], 'Physical inventory differs')
    x = array([a['xyz_A'] for a in physical]); require(np.array_equal(x, array(request['all_physical_coordinates_A'])), 'Raw request coordinates differ')
    mapping = read(p['mapping']); require(mapping['metal'] == metal and mapping['physical_atoms'] == physical, 'Mapping coordinates differ')
    gradient = read(p['assembled_gradient']); require(gradient['source_ids'] == ids, 'Gradient source order differs')
    g = array(gradient['gradient_kcal_mol_A']); cg = {k: array(gradient['component_gradients_kcal_mol_A'][k]) for k in COMPONENTS}
    require(np.max(np.abs(sum(cg.values())-g)) <= 1e-10, 'Total gradient closure differs')
    mm = read(p['MM_result']); trial = read(p['trial'])
    rows = [r for r in mm['rows'] if r['task_id'] == trial['task_id'] and r['metal'] == metal]
    require(mm['status'] == 'complete' and len(rows) == 1 and rows[0]['status'] == 'complete', 'Native MM row unavailable')
    mmrow = rows[0]; require(equal_pin(mmrow['trial'], p['trial']), 'MM trial lineage differs')
    require(set(mmrow['components_kcal_mol']) == set(GROUPS) and mmrow['components_kcal_mol'] == p['native_energies']['MM_groups_kcal_mol'], 'Native MM components differ')
    mmids = mmrow['source_ids']; require(len(mmids) == len(set(mmids)) == 2049 and set(mmids) == set(ids), 'Native MM source inventory differs')
    order = [mmids.index(s) for s in ids]
    mg = {k: array(mmrow['component_gradients_kcal_mol_A'][k])[order] for k in GROUPS}
    require(np.max(np.abs(sum(mg.values())-cg['MM_ledger'])) <= 1e-10, 'Mapped native MM gradient closure differs')
    close(math.fsum(mmrow['components_kcal_mol'].values()), mmrow['total_kcal_mol'], 1e-7, 'MM scalar closure differs')
    qm = {ids.index(a['source_id']) for a in mapping['core'] if a['kind'] == 'source'}
    targets = [a for a in mapping['core'] if a['kind'] == 'metal']
    require(len(targets) == 1, 'Unique local target required')
    t = targets[0]['source_site']; target_id = '/'.join((t['chain'], str(t['resid']), t['insertion_code'], t['name']))
    qm.add(ids.index(target_id))
    boundaries = []
    for cap in mapping['caps']:
        require(ids[cap['retained_index']] == cap['retained_id'] and ids[cap['omitted_index']] == cap['omitted_id']
            and cap['retained_index'] in qm and cap['omitted_index'] not in qm, 'Boundary physical mapping differs')
        boundaries.append((cap['retained_id'], cap['omitted_id']))
    require(len(qm) == len(mapping['core'])-len(mapping['caps']), 'Physical electronic core inventory differs')
    return dict(raw=p, request=request, physical=physical, ids=ids, x=x, gradient=g,
        components=cg, mm_gradients=mg, mm_energies=mmrow['components_kcal_mol'], qm=qm,
        boundaries=boundaries, target_id=target_id, provenance=provenance)


def load_search(record, metal):
    if record is None:
        return {'status': 'unavailable', 'metal': metal, 'reason': 'Actual completed search result not supplied'}
    r = read(record); cfg = read(r['config']); verify(r['implementation'])
    require(equal_pin(r['implementation'], cfg['implementation']) and r['stage'] == 'run' and r['metal'] == metal, 'Search implementation/state differs')
    row = {'status': 'unavailable', 'metal': metal, 'result': record, 'execution_status': r['status'], 'failure': r.get('failure')}
    if r['status'] != 'complete' or 'search' not in r:
        row['reason'] = 'Search execution unavailable'; return row
    s = read(r['search']); row.update(search=r['search'], termination=s['termination'], converged=s['converged'],
        full_force=s['final_full_force'], accepted_steps=s['accepted_steps'], distinct_evaluations=s['distinct_evaluations'])
    if s['status'] != 'complete' or not s['physical_state_admitted']:
        row['reason'] = 'Search not admitted in original bound state'; return row
    start, end = point(s['start_point'], metal), point(s['terminal_point'], metal)
    require(start['ids'] == end['ids'] and start['qm'] == end['qm'] and start['boundaries'] == end['boundaries']
        and s['start_point']['identity_sha256'] == s['terminal_point']['identity_sha256']
        and s['start_point']['state'] == s['terminal_point']['state'], 'Endpoint topology/recipe/state changed')
    guard = read(r['guard']); require(guard['physical_source_ids'] == start['ids'], 'Guard physical ordering differs')
    source = read(guard['source_inputs']); exported = read(source['export']); atoms = read(exported['atoms'])
    metadata = {a['id']: a for a in atoms}; ids = start['ids']
    protein = [i for i, aid in enumerate(ids) if metadata.get(aid, {}).get('kind') == 'protein_source']
    water = [i for i, aid in enumerate(ids) if metadata.get(aid, {}).get('kind') == 'retained_crystal_water']
    require(protein == guard['protein_indices'] and len(protein) == 1544 and len(water) == 501, 'Peptide/water inventory differs')
    ions = [v['index'] for v in guard['state_specs']]
    require(len(ions) == 4 and set(protein+water+ions) == set(range(2049)), 'Physical group partition differs')
    heavy = [i for i in protein if start['physical'][i]['element'] != 'H']
    groups = {'peptide_all': protein, 'peptide_heavy': heavy,
        'peptide_QM_heavy': [i for i in heavy if i in start['qm']],
        'peptide_outside_QM_heavy': [i for i in heavy if i not in start['qm']],
        'water_all': water, 'water_oxygen': [i for i in water if start['physical'][i]['element'] == 'O'],
        'target': [ids.index(start['target_id'])], 'Nd_spectators': [i for i in ions if ids[i] != start['target_id']]}
    require(groups['water_oxygen'] == guard['water_oxygen_indices'], 'Water oxygen inventory differs')
    works = {k: end['mm_energies'][k]-start['mm_energies'][k] for k in GROUPS}
    close(math.fsum(works.values()), s['work_from_start']['MM_ledger'], 1e-7, 'Global MM work replay differs')
    row.update(status='available', initial_point=s['start_point'], terminal_point=s['terminal_point'],
        MM_group_initial_kcal_mol=start['mm_energies'], MM_group_terminal_kcal_mol=end['mm_energies'],
        MM_group_work_kcal_mol=works, MM_total_work_kcal_mol=math.fsum(works.values()),
        source_metadata=source['export'], native_bonds=exported['bonds'], state=s['start_point']['state'],
        identity_sha256=s['start_point']['identity_sha256'])
    return dict(row=row, start=start, end=end, metadata=metadata, groups=groups,
                bonds=read(exported['bonds']), search=s)


def align(data):
    x, y = data['start']['x'], data['end']['x']; indices = data['groups']['peptide_heavy']
    xc, yc = x[indices].mean(axis=0), y[indices].mean(axis=0)
    u, singular, vt = np.linalg.svd((y[indices]-yc).T@(x[indices]-xc))
    sign = 1.0 if np.linalg.det(u@vt) >= 0 else -1.0
    rotation = u@np.diag([1., 1., sign])@vt
    require(abs(np.linalg.det(rotation)-1.) <= 1e-10 and np.max(np.abs(rotation.T@rotation-np.eye(3))) <= 1e-10, 'Improper alignment')
    translation = xc-yc@rotation; aligned = y@rotation+translation
    delta = aligned-x; raw = y-x
    fit_before = float(np.sum(raw[indices]**2)); fit_after = float(np.sum(delta[indices]**2))
    require(fit_after <= fit_before+1e-10, 'Alignment increases protein-heavy displacement')
    data.update(rotation=rotation, delta=delta, raw_delta=raw, aligned=aligned)
    return {'initial_point': data['row']['initial_point'], 'terminal_point': data['row']['terminal_point'],
        'fit_atom_ids': [data['start']['ids'][i] for i in indices], 'method': 'equal-weight whole-protein-heavy proper Kabsch; row-vector convention',
        'terminal_to_initial_rotation': rotation.tolist(), 'terminal_to_initial_translation_A': translation.tolist(),
        'initial_centroid_A': xc.tolist(), 'terminal_centroid_A': yc.tolist(), 'singular_values_A2': singular.tolist(),
        'reflection_removed': sign < 0, 'raw_fit_squared_displacement_A2': fit_before,
        'aligned_fit_squared_displacement_A2': fit_after, 'reference_scoring_uses_raw_coordinates': True}


def motion(delta, indices):
    if not indices:
        return {'atom_count': 0, 'squared_displacement_A2': 0., 'RMS_displacement_A': None, 'maximum_displacement_A': None}
    values = np.sum(delta[indices]**2, axis=1)
    return {'atom_count': len(indices), 'squared_displacement_A2': float(values.sum()),
        'RMS_displacement_A': float(np.sqrt(values.mean())), 'maximum_displacement_A': float(np.sqrt(values.max()))}


def force_descriptors(data):
    result = {}
    for endpoint in ('start', 'end'):
        p = data[endpoint]; ids = p['ids']; components = dict(total=p['gradient'], **p['components'])
        components.update({'MM/'+k: v for k, v in p['mm_gradients'].items()})
        result[endpoint] = {}
        for name, g in components.items():
            norms = np.linalg.norm(g, axis=1); order = sorted(range(2049), key=lambda i: (-float(norms[i]), ids[i]))
            residues = {}
            for i in data['groups']['peptide_all']:
                rid = residue_label(residue(ids[i])); row = residues.setdefault(rid, {'squared_norm': 0., 'maximum': 0., 'atom_count': 0})
                row['squared_norm'] += float(norms[i]**2); row['maximum'] = max(row['maximum'], float(norms[i])); row['atom_count'] += 1
            result[endpoint][name] = {'maximum_atom_kcal_mol_A': float(norms.max()),
                'RMS_atom_vector_kcal_mol_A': float(np.sqrt(np.mean(norms**2))),
                'top_atoms': [{'id': ids[i], 'norm_kcal_mol_A': float(norms[i]), 'in_real_QM': i in p['qm']} for i in order[:20]],
                'all_peptide_residues': residues}
        if endpoint == 'end':
            for k in ('maximum_atom_kcal_mol_A', 'RMS_atom_vector_kcal_mol_A'):
                close(result[endpoint]['total'][k], data['search']['final_full_force'][k], 1e-10, 'Search residual force replay differs')
    return result


def topology(data):
    ids = data['start']['ids']; protein = data['groups']['peptide_all']; byres = {}
    for i in protein:
        byres.setdefault(residue(ids[i]), []).append(i)
    edges = set()
    for a, b in data['bonds']:
        if a not in data['metadata'] or b not in data['metadata']:
            continue
        if data['metadata'][a]['kind'] != 'protein_source' or data['metadata'][b]['kind'] != 'protein_source':
            continue
        if a.rsplit('/', 1)[1] == 'N' and b.rsplit('/', 1)[1] == 'C':
            a, b = b, a
        if a.rsplit('/', 1)[1] == 'C' and b.rsplit('/', 1)[1] == 'N' and residue(a) != residue(b):
            require(residue(a)[0] == residue(b)[0], 'Interchain peptide bond needs explicit handling')
            edges.add((residue(a), residue(b)))
    successors = {}; predecessors = {}
    for a, b in edges:
        require(a not in successors and b not in predecessors, 'Branched peptide topology')
        successors[a] = b; predecessors[b] = a
    windows = [(a, successors[a], successors[successors[a]]) for a in sorted(byres)
               if a in successors and successors[a] in successors]
    require(windows, 'No actual contiguous three-residue windows')
    return byres, edges, windows


def selected_context(selected, data, byres, edges):
    residues = [tuple(k) for k in selected['residue_keys']]; ids = data['start']['ids']
    window = [i for r in residues for i in byres[r]]
    outside = [i for i in window if i in data['groups']['peptide_outside_QM_heavy']]
    heavy = [i for i in window if data['start']['physical'][i]['element'] != 'H']
    qheavy = [i for i in data['start']['qm'] if data['start']['physical'][i]['element'] != 'H']
    qres = {residue(ids[i]) for i in data['groups']['peptide_all'] if i in data['start']['qm']}
    adjacency = {r: set() for r in byres}
    for a, b in edges:
        adjacency[a].add(b); adjacency[b].add(a)
    frontier = set(residues); visited = set(frontier); depth = 0
    while not frontier.intersection(qres) and frontier:
        frontier = set().union(*(adjacency[r] for r in frontier))-visited
        visited.update(frontier); depth += 1
    sequence_distance = depth if frontier else None
    remaining = qres|set(residues); count = 0
    while remaining:
        todo = {min(remaining)}; remaining -= todo; count += 1
        while todo:
            todo = set().union(*(adjacency[r] for r in todo))&remaining; remaining -= todo
    endpoints = {}
    for endpoint in ('start', 'end'):
        p = data[endpoint]
        distances = np.linalg.norm(p['x'][outside, None, :]-p['x'][qheavy][None, :, :], axis=2)
        gradients = dict(p['components'], **{'MM/'+k: v for k, v in p['mm_gradients'].items()})
        gradients['total'] = p['gradient']
        endpoints[endpoint] = {'minimum_external_window_to_QM_heavy_A': float(distances.min()),
            'median_nearest_external_window_to_QM_heavy_A': float(np.median(distances.min(axis=1))),
            'gradient_dot_aligned_window_displacement_kcal_mol': {k: float(np.sum(
                (g@data['rotation'] if endpoint == 'end' else g)[outside]*data['delta'][outside])) for k, g in gradients.items()}}
    return {'all_window_atom_ids': [ids[i] for i in window], 'ranked_outside_QM_heavy_atom_ids': [ids[i] for i in outside],
        'whole_window_heavy_atom_count': len(heavy), 'existing_QM_peptide_residues': [residue_label(r) for r in sorted(qres)],
        'minimum_peptide_edges_to_QM_residue': sequence_distance,
        'QM_residue_and_window_union_peptide_components': count, 'endpoint_descriptors': endpoints,
        'projection_is_integrated_work_or_local_energy': False,
        'feasibility': 'Candidate only; charged-boundary and local hydrogen-bond-network review still required.'}


def selection(cells):
    la, dy = cells['La'], cells['Dy']
    require(la['start']['ids'] == dy['start']['ids'] and la['start']['qm'] == dy['start']['qm']
        and la['start']['boundaries'] == dy['start']['boundaries'] and la['bonds'] == dy['bonds']
        and la['groups'] == dy['groups'] and la['row']['identity_sha256'] == dy['row']['identity_sha256'], 'Paired model/physical topology differs')
    byres, edges, window_keys = topology(la); ids = la['start']['ids']; outside = set(la['groups']['peptide_outside_QM_heavy'])
    rrows = []
    for rid, indices in sorted(byres.items()):
        eligible = [i for i in indices if i in outside]
        values = {m: float(np.sum(cells[m]['delta'][eligible]**2)) for m in METALS}
        rrows.append({'residue_key': list(rid), 'residue_id': residue_label(rid), 'outside_QM_heavy_atoms': len(eligible),
                      'per_metal_squared_displacement_A2': values, 'pooled_squared_displacement_A2': math.fsum(values.values())})
    total = math.fsum(r['pooled_squared_displacement_A2'] for r in rrows); bykey = {tuple(r['residue_key']): r for r in rrows}
    boundaries = {i for pair in la['start']['boundaries'] for i in pair}; windows = []
    for keys in window_keys:
        indices = [i for k in keys for i in byres[k]]; values = {m: math.fsum(bykey[k]['per_metal_squared_displacement_A2'][m] for k in keys) for m in METALS}
        weight = math.fsum(values.values())
        windows.append({'residue_keys': [list(k) for k in keys], 'residue_ids': [residue_label(k) for k in keys],
            'per_metal_squared_displacement_A2': values, 'pooled_squared_displacement_A2': weight,
            'fraction_of_all_outside_QM_peptide_motion': weight/total if total else None,
            'outside_QM_heavy_atom_count': sum(i in outside for i in indices),
            'QM_atom_ids': [ids[i] for i in indices if i in la['start']['qm']],
            'boundary_anchor_ids': sorted(ids[i] for i in indices if ids[i] in boundaries)})
    windows.sort(key=lambda r: (-r['pooled_squared_displacement_A2'], tuple(tuple(k) for k in r['residue_keys'])))
    for rank, row in enumerate(windows, 1): row['rank'] = rank
    ordered = sorted(rrows, key=lambda r: (-r['pooled_squared_displacement_A2'], tuple(r['residue_key'])))
    cumulative = 0.; coverage = {}
    for n, r in enumerate(ordered, 1):
        cumulative += r['pooled_squared_displacement_A2']; r['cumulative_fraction'] = cumulative/total if total else None
        for level in (.5, .8, .9):
            if total and cumulative >= level*total and str(level) not in coverage: coverage[str(level)] = n
    selected = windows[0] if windows[0]['pooled_squared_displacement_A2'] > 0 else None
    chosen = {'status': 'candidate_for_feasibility_review' if selected else 'unavailable_zero_motion',
        'selected_window': selected, 'selection_uses_both_metals': True,
        'reference_outcomes_used': False, 'energy_or_preference_sign_used': False,
        'total_outside_QM_peptide_squared_motion_A2': total,
        'effective_participating_residues': total*total/math.fsum(r['pooled_squared_displacement_A2']**2 for r in rrows) if total else None,
        'residue_counts_for_descriptive_coverage': coverage, 'automatic_localization_pass_threshold': None,
        'one_region_explains_global_MM_release': 'not established; MM group work is global and spatial motion is not a residue energy decomposition',
        'single_expanded_region_feasibility': 'unresolved pending complete chemistry and hydrogen-bond-network review'}
    if selected:
        chosen['per_metal_context'] = {m: selected_context(selected, cells[m], byres, edges) for m in METALS}
    return chosen, ordered, windows


def run(cfg, out, result):
    require(equal_pin(cfg['implementation'], pin(__file__)) and set(cfg['search_results']) == set(METALS), 'Source or two-metal denominator differs')
    verify(cfg['specification']); verify(cfg['plan']); cells = {}; summaries = result['searches']
    for metal in METALS:
        loaded = load_search(cfg['search_results'][metal], metal)
        if 'row' in loaded:
            summaries[metal] = loaded['row']; cells[metal] = loaded
        else: summaries[metal] = loaded
    alignments = {}; mechanics = {}; atom_rows = []
    for metal, data in cells.items():
        alignments[metal] = align(data)
        mechanics[metal] = {'MM_initial_kcal_mol': data['row']['MM_group_initial_kcal_mol'],
            'MM_terminal_kcal_mol': data['row']['MM_group_terminal_kcal_mol'],
            'MM_work_kcal_mol': data['row']['MM_group_work_kcal_mol'], 'MM_total_work_kcal_mol': data['row']['MM_total_work_kcal_mol'],
            'force_residuals': force_descriptors(data),
            'group_motion': {k: {'raw': motion(data['raw_delta'], indices), 'aligned': motion(data['delta'], indices)} for k, indices in data['groups'].items()},
            'global_group_energies_are_not_residue_energy_decomposition': True}
        for i, aid in enumerate(data['start']['ids']):
            atom_rows.append({'metal': metal, 'atom_id': aid, 'element': data['start']['physical'][i]['element'],
                'in_real_QM': i in data['start']['qm'], 'raw_displacement_A': float(np.linalg.norm(data['raw_delta'][i])),
                'aligned_displacement_A': float(np.linalg.norm(data['delta'][i])),
                'initial_total_gradient_norm_kcal_mol_A': float(np.linalg.norm(data['start']['gradient'][i])),
                'terminal_total_gradient_norm_kcal_mol_A': float(np.linalg.norm(data['end']['gradient'][i]))})
    result['alignments'] = write(out/'ALIGNMENTS.json', alignments)
    result['mechanics'] = write(out/'MECHANICS.json', mechanics)
    if len(cells) == 2:
        chosen, residues, windows = selection(cells)
    else:
        chosen = {'status': 'unavailable', 'reason': 'Both admitted metals required; no La-only selection', 'selected_window': None}
        residues = []; windows = []
    result.update(selection=chosen, residues=write(out/'RESIDUES.json', residues), windows=write(out/'WINDOWS.json', windows))
    fields = ['metal', 'atom_id', 'element', 'in_real_QM', 'raw_displacement_A', 'aligned_displacement_A',
        'initial_total_gradient_norm_kcal_mol_A', 'terminal_total_gradient_norm_kcal_mol_A']
    with (out/'ATOM_MOTION.tsv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t'); writer.writeheader(); writer.writerows(atom_rows)
    lines = ['# Mex saved motion and mechanics', '', 'Segment selection: '+chosen['status']+'.', '',
        '| Metal | Bonded | MM LJ | MM Coulomb | QM/MM LJ | Total MM work (kcal/mol) |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for m in METALS:
        if m not in cells:
            lines.append('| '+m+' | unavailable | unavailable | unavailable | unavailable | unavailable |'); continue
        row = mechanics[m]; values = [row['MM_work_kcal_mol'][k] for k in GROUPS]+[row['MM_total_work_kcal_mol']]
        lines.append('| '+m+' | '+' | '.join(f'{v:.6f}' for v in values)+' |')
    lines += ['', '| Motion group | La aligned RMS / maximum (Å) | Dy aligned RMS / maximum (Å) |',
        '| --- | ---: | ---: |']
    for group in ('peptide_outside_QM_heavy', 'peptide_QM_heavy', 'water_oxygen', 'target', 'Nd_spectators'):
        values = []
        for metal in METALS:
            if metal not in mechanics:
                values.append('unavailable'); continue
            value = mechanics[metal]['group_motion'][group]['aligned']
            values.append(f"{value['RMS_displacement_A']:.4f} / {value['maximum_displacement_A']:.4f}")
        lines.append('| '+group+' | '+' | '.join(values)+' |')
    if chosen.get('selected_window'):
        w = chosen['selected_window']
        lines += ['', 'Candidate window: **'+', '.join(w['residue_ids'])+'**.',
            f"It accounts for {100*w['fraction_of_all_outside_QM_peptide_motion']:.2f}% of pooled outside-QM peptide heavy-atom squared motion after whole-protein alignment.",
            'La/Dy contributions (Å²): '+str(w['per_metal_squared_displacement_A2'])+'.',
            'Effective participating residues: '+str(chosen['effective_participating_residues'])+'.',
            'Residues covering 50/80/90% of motion: '+str(chosen['residue_counts_for_descriptive_coverage'])+'.']
    lines += ['', 'Global MM group lowering is not a residue energy decomposition. Water/peptide motion and full component force locations are recorded separately in MECHANICS.json.', '',
        'The selected segment is one local stiffness-test candidate. Its coverage and QM connectivity do not establish that one expanded region can explain the full release; complete chemistry, charge boundaries and the local hydrogen-bond network still need review. No localization threshold, preference sign or new reference outcome selected it.', '',
        'Searches remain finite and their raw force residuals remain explicit. No model evaluations or new reference launches were performed.', '']
    (out/'REPORT.md').write_text('\n'.join(lines)); result['report'] = pin(out/'REPORT.md')
    result['status'] = 'complete' if len(cells) == 2 else 'unavailable'


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--stage', choices=('bind', 'run'), required=True)
    p.add_argument('--config', type=Path, required=True); p.add_argument('--output', type=Path, required=True); args = p.parse_args()
    cfg = json.loads(args.config.read_text())
    if args.stage == 'bind':
        cfg.update(protocol=PROTOCOL, implementation=pin(__file__), specification=pin(Path(__file__).with_name('SPEC.md')))
        for k in METALS:
            value = cfg['search_results'][k]
            if isinstance(value, str): cfg['search_results'][k] = pin(value)
            elif value is not None: verify(value)
        if isinstance(cfg['plan'], str): cfg['plan'] = pin(cfg['plan'])
        else: verify(cfg['plan'])
        write(args.output, cfg); return 0
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Saved numerical analysis requires Slurm')
    global np
    import numpy as np
    require(cfg['protocol'] == PROTOCOL, 'Protocol differs')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic(); cpu = time.process_time()
    result = {'protocol': PROTOCOL, 'status': 'execution_failure', 'implementation': pin(__file__),
        'config': pin(args.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'new_molecular_calls': 0,
        'new_reference_calls': 0, 'new_optimization_steps': 0, 'declared_searches': 2,
        'searches': {m: {'metal': m, 'status': 'unavailable_not_reviewed'} for m in METALS},
        'selection': {'status': 'unavailable', 'selected_window': None}}
    try:
        run(cfg, out, result); verify(result['config'])
    except Exception as exc:
        result.update(failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick, process_CPU_seconds=time.process_time()-cpu, numpy_version=np.__version__)
    write(out/'RESULT.json', result); print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
