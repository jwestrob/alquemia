"""Summarize actual complete-context PQQ endpoints without fitting a new score."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import math
from pathlib import Path
import statistics

from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL


def summarize(result_path):
    result = read_json(result_path)
    manifest = read_json(verify(result['manifest']))
    collection = read_json(verify(result['collection']))
    endpoints = [read_json(verify(p)) for p in collection['results']]
    energies = {r['task_id']: r['energy_eV'] for r in endpoints}
    energies.update({r['task']['task_id']: r['accepted']['energy_eV'] for r in manifest['reused']})
    rows = []
    for old in result['rows']:
        case = old['case']
        item = next(c for c in manifest['cases'] if c['source']['case'] == case)
        prep = read_json(verify(item['preparation']))
        endpoint_eV = {metal: energies[case+'__expanded__'+metal] for metal in ('Ca', 'La')}
        contrast = (endpoint_eV['Ca'] - endpoint_eV['La']) * EV_TO_KCAL
        if contrast != old['context_R_model_kcal_mol']:
            raise InvalidArtifact('reported contrast differs from actual endpoint algebra')
        core = xyz(verify(item['source']['endpoints']['La']['xyz']))
        task = next(t for t in manifest['tasks'] + [r['task'] for r in manifest['reused']]
                    if t['case_id'] == case and t['metal'] == 'La')
        context = xyz(verify(task['xyz']))
        for i, j in prep['core_to_context'].items():
            if core[int(i)] != context[j]:
                raise InvalidArtifact('retained original coordinate changed')
        new_atoms = [a for a in prep['mapping']['source_to_qm'] if a['kind'] == 'source'
                     and a['qm_index'] not in set(prep['core_to_context'].values())]
        charged = [a for a in prep['added_fragments'] if a['formal_charge']]
        charged_polar_distances = []
        for atom in new_atoms:
            s = atom['source']
            if s['element'] not in ('N', 'O', 'S'):
                continue
            if any(all(s[k] == f['source'][k] for k in ('chain', 'resnum', 'insertion_code', 'resname')) for f in charged):
                charged_polar_distances.append(math.dist(context[atom['qm_index']][1:], context[0][1:]))
        rows.append(dict(old, endpoint_energies_eV=endpoint_eV,
                         delta_R_model_kcal_mol=contrast-old['core_R_model_kcal_mol'],
                         core_charge_Ca=item['source']['endpoints']['Ca']['charge'],
                         context_charge_Ca=item['source']['endpoints']['Ca']['charge']+prep['added_formal_charge'],
                         context_elements=dict(Counter(a[0] for a in context)),
                         original_atoms=len(core),
                         retained_original_coordinates_exact=True,
                         removed_caps=len(prep['removed_caps']),
                         charged_added_fragments=charged,
                         charged_added_polar_metal_distances_A=charged_polar_distances))
    calibration = [r for r in rows if r['role'] == 'calibration']
    if len(calibration) != 25:
        raise InvalidArtifact('canonical calibration denominator differs')
    groups = defaultdict(list)
    for r in calibration:
        groups[(r['expected_class'], r['added_formal_charge'])].append(r)
    descriptive = []
    for (label, charge), members in sorted(groups.items()):
        changes = [r['delta_R_model_kcal_mol'] for r in members]
        descriptive.append({'expected_class': label, 'added_formal_charge': charge, 'n_structures': len(members),
                            'cases': [r['case'] for r in members], 'shift_min': min(changes),
                            'shift_max': max(changes), 'shift_mean': statistics.mean(changes)})
    class_ranges = {}
    for label in ('Ca', 'La'):
        members = [r for r in calibration if r['expected_class'] == label]
        class_ranges[label] = {kind: {'min': min(r[kind] for r in members),
                                     'max': max(r[kind] for r in members),
                                     'width': max(r[kind] for r in members)-min(r[kind] for r in members)}
                               for kind in ('core_R_model_kcal_mol', 'context_R_model_kcal_mol')}
    gap = class_ranges['La']['context_R_model_kcal_mol']['min']-class_ranges['Ca']['context_R_model_kcal_mol']['max']
    if gap != result['gap_model_kcal_mol']:
        raise InvalidArtifact('frozen calibration gap differs')
    return {'result': record(result_path), 'protocol_id': manifest['protocol_id'], 'rows': rows,
            'bands': {'Ca_max_model_kcal_mol': result['Ca_band_max'], 'La_min_model_kcal_mol': result['La_band_min'],
                      'inconclusive_between': True, 'calibration_only': True},
            'core_gap_model_kcal_mol': class_ranges['La']['core_R_model_kcal_mol']['min']-class_ranges['Ca']['core_R_model_kcal_mol']['max'],
            'context_gap_model_kcal_mol': gap, 'class_ranges': class_ranges, 'descriptive_charge_strata': descriptive,
            'component_attribution': 'unresolved: no charge/density or atomic-energy components archived; geometry fixed, composition and charge coupled',
            'cost': {'new_endpoints': len(endpoints), 'reused_context_endpoints': len(manifest['reused']),
                     'new_DFT_endpoints': 0, 'collection_wall_seconds': collection['wall_seconds'],
                     'sum_endpoint_evaluation_seconds': sum(r['evaluation_seconds'] for r in endpoints),
                     'median_endpoint_evaluation_seconds': statistics.median(r['evaluation_seconds'] for r in endpoints),
                     'sum_endpoint_wall_seconds': sum(r['wall_seconds'] for r in endpoints),
                     'peak_cuda_allocated_bytes': max(r['peak_cuda_allocated_bytes'] for r in endpoints),
                     'peak_worker_host_RSS_KiB': max(r['peak_host_RSS_KiB'] for r in endpoints)}}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--result', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    write_new(a.output, summarize(a.result))


if __name__ == '__main__':
    main()
