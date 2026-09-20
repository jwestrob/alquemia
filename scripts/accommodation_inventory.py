"""Read-only donor organization diagnosis of pinned PQQ reference/PLM structures."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import math
from pathlib import Path
import statistics

from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz

SCHEMA = 'pqq_accommodation_geometry_inventory_v1'
DONORS = {'anchor_glutamate': ('OE1', 'OE2'), 'anchor_asparagine': ('OD1',),
          'catalytic_aspartate': ('OD1', 'OD2'),
          'extra_acidic_ligand_homolog': ('OD1', 'OD2', 'OE1', 'OE2'),
          'catalytic_asp_cationic_partner': ('NH1', 'NH2', 'NE', 'NZ'),
          'pqq_cofactor': ('N6', 'O5', 'O7A', 'O7B')}


def geometry(path):
    p = Path(path).resolve(); m = read_json(p)
    pin = m['outputs']['La_xyz']; pin = dict(pin) if isinstance(pin, dict) else record(p.parent / pin)
    if not Path(pin['path']).is_absolute(): pin['path'] = str(p.parent / pin['path'])
    rows = xyz(verify(pin)); metal = rows[0][1:]
    if rows[0][0] != 'La': raise InvalidArtifact('expected source La endpoint')
    fragments = m['qm_fragments']; distances = {}; roles = {}
    flat = [a for f in fragments for a in f['atom_records']]
    if len(rows) != len(flat) + 1: raise InvalidArtifact('source fragment/endpoint count differs')
    for row, a in zip(rows[1:], flat):
        if row[0] != a['element'] or math.dist(row[1:], a['xyz_A']) > 2e-6:
            raise InvalidArtifact('fragment mapping does not match actual endpoint')
    naming = m['pqq']['residues'][0]['atom_naming_schema']
    if naming == 'pdb_ccd_pqq_v1': pqq_names = {n:n for n in DONORS['pqq_cofactor']}
    elif naming == 'protenix_generic_pqq_v1': pqq_names = {'N2':'N6', 'O1':'O5', 'O5':'O7A', 'O6':'O7B'}
    else: raise InvalidArtifact('unsupported PQQ atom naming')
    for f in fragments:
        role = f['role']; allowed = DONORS.get(role, ())
        if role == 'pqq_cofactor': allowed = tuple(pqq_names)
        contacts = []
        for a in f['atom_records']:
            if a['name'] in allowed:
                d = math.dist(metal, a['xyz_A'])
                name = pqq_names[a['name']] if role == 'pqq_cofactor' else a['name']
                distances[role + '.' + name] = d
                contacts.append(d)
        roles[role] = {'source_id': f['id'], 'formal_charge': f['formal_charge'],
                       'donor_distances_A': sorted(contacts)}
        if contacts:
            distances[role + '.nearest'] = min(contacts)
            if len(contacts) > 1: distances[role + '.second'] = sorted(contacts)[1]
    pqq = [distances['pqq_cofactor.' + a] for a in ('N6', 'O5')]
    pqq.append(min(distances['pqq_cofactor.' + a] for a in ('O7A', 'O7B')))
    distances['PQQ_three_donor_mean'] = statistics.mean(pqq)
    distances['PQQ_three_donor_max'] = max(pqq)
    distances['CN_3p1'] = m['coordination']['coordination_number']
    h = m.get('hydrogen_geometry_audit', {})
    distances['source_H_force_rms_kJ_mol_nm'] = h.get('original_H_force_rms_kJ_mol_nm')
    c = m.get('confidence', {})
    distances['protein_ion_iptm'] = c.get('value')
    return {'carve': record(p), 'xyz': pin, 'features': distances, 'roles': roles,
            'formula': dict(Counter(r[0] for r in rows[1:])),
            'charge_Ca': m['charge_ledger']['Ca_total'],
            'charge_La': m['charge_ledger']['La_total'],
            'water_policy': m['fixed_core']['water_policy'], 'pqq_source_atom_naming': naming,
            'source_structure': m['source_structure'], 'source_cif': m.get('source_cif')}


def run(plm_results, reference_results, plan, output):
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    plm = read_json(plm_results); ref = read_json(reference_results)
    prep = read_json(verify(ref['preparation']))
    sources = {t['panel_id']: t for t in prep['targets']}
    rows = []
    for r in ref['scores']:
        rows.append({'id': r['panel_id'], 'cohort': 'reference', 'label': r['class'],
                     'score': r['S_aquo_gauge_kcal_mol'], 'status': 'available',
                     **geometry(verify(sources[r['panel_id']]['manifest']))})
    for r in plm['per_target']:
        base = {'id': r['target_id'], 'cohort': 'PLM_prediction', 'label': None,
                'predicted_band': r.get('protocol_band'), 'score': r.get('S_kcal_mol')}
        if r['status'] != 'scored':
            rows.append({**base, 'status': 'unavailable', 'reason': r.get('reason', r['status'])})
            continue
        try:
            rows.append({**base, 'status': 'available', **geometry(r['carve_manifest'])})
        except (InvalidArtifact, KeyError, ValueError, OSError) as exc:
            rows.append({**base, 'status': 'unavailable', 'reason': str(exc)})
    available = [r for r in rows if r['status'] == 'available']
    keys = sorted({k for r in available for k in r['features']})
    grouped = defaultdict(list)
    for r in available:
        grouped[r['cohort'] + ':' + (r['label'] or r['predicted_band'])].append(r)
    distributions = {}
    for name, group in grouped.items():
        distributions[name] = {'n': len(group), 'features': {}}
        for key in keys:
            values = [r['features'][key] for r in group if r['features'].get(key) is not None]
            if values:
                distributions[name]['features'][key] = {'n': len(values), 'min': min(values),
                    'median': statistics.median(values), 'max': max(values)}
    from scipy.stats import spearmanr
    correlations = {}
    for cohort in ('reference', 'PLM_prediction'):
        correlations[cohort] = {}
        for key in keys:
            vals = [(r['features'].get(key), r['score']) for r in available if r['cohort'] == cohort
                    and r['features'].get(key) is not None]
            if len(vals) >= 3 and len({x for x, y in vals}) > 1:
                rho = float(spearmanr([x for x, y in vals], [y for x, y in vals]).statistic)
                correlations[cohort][key] = {'n': len(vals), 'rho': rho}
    result = {'schema': SCHEMA, 'plan': record(plan), 'implementation': record(__file__),
              'inputs': {'PLM': record(plm_results), 'reference': record(reference_results)},
              'counts': dict(Counter(r['cohort'] + ':' + r['status'] for r in rows)),
              'new_molecular_calls': 0, 'labels_fitted': False, 'baseline_changed': False,
              'interpretation': 'descriptive geometry/score association; PLM has no accuracy labels',
              'distributions': distributions, 'correlations': correlations, 'rows': rows}
    write_new(out / 'result.json', result)
    with (out / 'features.tsv').open('w') as f:
        fields = ['id', 'cohort', 'label', 'predicted_band', 'score', 'status'] + keys
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore'); w.writeheader()
        for r in rows: w.writerow({**r, **r.get('features', {})})
    return {'counts': result['counts'], 'result': record(out / 'result.json')}


if __name__ == '__main__':
    import json
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('plm-results', 'reference-results', 'plan', 'output'): p.add_argument('--' + name, required=True)
    args = p.parse_args()
    print(json.dumps(run(args.plm_results, args.reference_results, args.plan, args.output), indent=2))
