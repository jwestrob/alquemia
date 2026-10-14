#!/usr/bin/env python3
"""Saved benchmark case review; no new endpoints, selection rules or calibration."""
import argparse
from collections import Counter, defaultdict
import csv
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import time

EV = 23.06054783061903
HA = 627.509474
ORDER = ('origin', 'adaptive_Ca', 'adaptive_La')


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


@lru_cache(None)
def read_checked(path, sha):
    assert pin(path)['sha256'] == sha, path
    return json.loads(Path(path).read_text())


def checked(p):
    return read_checked(p['path'], p['sha256'])


def call(r, bands):
    return 'Ca' if r <= bands['Ca_max'] else 'La' if r >= bands['La_min'] else 'inconclusive'


def margin(r, expected, bands):
    return r-bands['La_min'] if expected == 'La' else bands['Ca_max']-r


def require_close(a, b):
    assert math.isfinite(a) and math.isfinite(b) and abs(a-b) <= 1e-8, (a, b)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--sha', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Slurm allocation required')
    start = time.monotonic()
    inp = {'path': str(Path(args.input).resolve()), 'sha256': args.sha}
    config = checked(inp)
    checked(config['completed_scores'])
    inventory = checked(config['inventory'])
    tasks = {k: v for field in ('canonical_native_collection', 'transfer_native_collection')
             for k, v in checked(config[field])['tasks'].items()}
    diagnostic = checked(config['selection_diagnosis'])
    diagnosis = {(r['record_id'], r['policy']): r for r in diagnostic['rows']}
    bands = config['calibrations']['operational']['bands']
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    rows, flat = [], []
    for case in config['cases']:
        rec, score = case['record'], case['score']
        selected, origin = score['pool']['operational'], score['pool']['origin']
        row = {k: case[k] for k in ('case_key', 'record_id', 'aliases', 'target')}
        row.update({k: rec[k] for k in ('source_case_id', 'biological_group', 'expected_class', 'source_conditioning_metal', 'family')})
        row.update(status='unavailable', reason=None)
        if selected is None or origin is None:
            row['reason'] = 'required_saved_pool_unavailable'
            rows.append(row)
            continue
        matrix = {z: {} for z in ('Ca', 'La')}
        for z in matrix:
            for q in ORDER:
                cell = rec['matrix'][z][q]
                native = tasks[cell['native_task_id']]
                assert native['status'] == 'available'
                detail = checked(native['result'])
                assert detail['status'] == 'available'
                g = detail['G_qvac_hartree']
                require_close(g, native['G_qvac_hartree'])
                raw = {}
                for medium in ('water', 'vacuum'):
                    old = checked(detail['input']['endpoints'][medium]['point'])
                    raw[medium] = old['backends']['xtb']['energy_hartree']
                matrix[z][q] = {'MACE_eV': cell['MACE_eV'], 'G_qvac_hartree': g,
                                'U_self_consistent_hartree': raw['water']-raw['vacuum'],
                                'native_components_qvac_hartree': detail['vacuum']['components_Eh'],
                                'native_result': native['result']}
        def parts(ca_q, la_q):
            a, b = matrix['Ca'][ca_q], matrix['La'][la_q]
            mace = (a['MACE_eV']-b['MACE_eV'])*EV
            solvent = (a['G_qvac_hartree']-b['G_qvac_hartree'])*HA
            sc = (a['U_self_consistent_hartree']-b['U_self_consistent_hartree'])*HA
            return {'MACE_contrast': mace, 'fixed_charge_solvent_contrast': solvent,
                    'self_consistent_solvent_contrast_at_same_geometries': sc,
                    'additional_redistribution_plus_intrinsic_contrast': sc-solvent,
                    'R': mace+solvent}
        original = parts('origin', 'origin')
        chosen = parts(selected['Ca']['candidate'], selected['La']['candidate'])
        require_close(original['R'], origin['R_model_kcal_mol'])
        require_close(chosen['R'], selected['R_model_kcal_mol'])
        shift = chosen['R']-original['R']
        require_close(shift, selected['Ca']['work_kcal_mol']-selected['La']['work_kcal_mol'])
        by_candidate = {q: parts(q, q) for q in ORDER}
        for q, v in by_candidate.items():
            v['same_operational_band_diagnostic_call'] = call(v['R'], bands)
        state = checked(rec['state_signature'])
        state_subset = {k: state[k] for k in ('core_microstates', 'endpoint_charges', 'multiplicities', 'water_inventory')}
        state_subset['state_signature'] = rec['state_signature']
        selected_call = score['decisions']['operational']['decision']
        assert call(chosen['R'], bands) == selected_call
        row.update(status='available', matrix=matrix, origin=original, selected=chosen,
                   geometry_choices={z: selected[z]['candidate'] for z in ('Ca', 'La')},
                   works={z: selected[z]['all_works_kcal_mol'] for z in ('Ca', 'La')},
                   paired_candidates=by_candidate, state=state_subset,
                   same_band_origin_call=call(original['R'], bands),
                   separately_calibrated_origin_call=score['decisions']['origin']['decision'],
                   selected_call=selected_call,
                   correct_class_margin_kcal_mol=margin(chosen['R'], rec['expected_class'], bands),
                   origin_same_band_correct_class_margin_kcal_mol=margin(original['R'], rec['expected_class'], bands),
                   distance_above_Ca_band_kcal_mol=chosen['R']-bands['Ca_max'],
                   distance_below_La_band_kcal_mol=bands['La_min']-chosen['R'],
                   accommodation_R_shift=shift,
                   accommodation_MACE_shift=chosen['MACE_contrast']-original['MACE_contrast'],
                   accommodation_solvent_shift=chosen['fixed_charge_solvent_contrast']-original['fixed_charge_solvent_contrast'],
                   baseline_call=score['baseline']['decisions']['operational']['decision'],
                   prior_model_selection_band_diagnosis=diagnosis[(rec['record_id'], 'operational')])
        rows.append(row)
        flat.append({**{k: row[k] for k in ('source_case_id', 'biological_group', 'expected_class', 'source_conditioning_metal', 'target',
                    'selected_call', 'same_band_origin_call', 'separately_calibrated_origin_call', 'correct_class_margin_kcal_mol',
                    'distance_above_Ca_band_kcal_mol', 'distance_below_La_band_kcal_mol', 'baseline_call',
                    'accommodation_R_shift', 'accommodation_MACE_shift', 'accommodation_solvent_shift')},
                     'R_origin': original['R'], 'R_selected': chosen['R'],
                     'Ca_candidate': selected['Ca']['candidate'], 'La_candidate': selected['La']['candidate'],
                     'Ca_work': selected['Ca']['work_kcal_mol'], 'La_work': selected['La']['work_kcal_mol']})
    peer_groups = defaultdict(list)
    for row in rows:
        peer_groups[(row['biological_group'], row['source_conditioning_metal'])].append(row)
    peers = []
    for (group, condition), cases in peer_groups.items():
        good = [r for r in cases if r['status'] == 'available']
        values = [r['selected']['R'] for r in good]
        peers.append({'biological_group': group, 'conditioning': condition, 'declared': len(cases),
                      'available': len(good), 'source_ids': [r['source_case_id'] for r in cases],
                      'calls': dict(Counter(r['selected_call'] for r in good)),
                      'minimum_R': min(values) if values else None, 'maximum_R': max(values) if values else None,
                      'median_R_description_only': statistics.median(values) if values else None})
    result = {'status': 'complete', 'input': inp, 'implementation': pin(__file__), 'bands': bands, 'rows': rows,
              'peer_groups': peers, 'inconclusive_triples': config['inconclusive_triples'],
              'original_unavailable': config['original_unavailable'], 'new_molecular_calls': 0,
              'production_changed': False, 'new_calibration': False,
              'limitations': ['Origin/same-candidate calls on pooled bands are mechanistic counterfactuals, not new released classifiers.',
                              'Peer medians describe scores only; no new grouping/voting rule.',
                              'Distances to bands are not uncertainties or affinities.'],
              'job_id': os.environ['SLURM_JOB_ID'], 'elapsed_seconds': time.monotonic()-start}
    with (output/'RESULT.json').open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    with (output/'SOURCES.tsv').open('x') as f:
        writer = csv.DictWriter(f, fieldnames=list(flat[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(flat)
    print(json.dumps({'targets': [r for r in flat if r['target']], 'peer_groups': peers}, indent=2))


if __name__ == '__main__':
    main()
