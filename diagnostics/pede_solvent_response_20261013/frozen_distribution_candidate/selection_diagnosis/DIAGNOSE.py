#!/usr/bin/env python3
"""Decompose completed candidate changes using saved matrices only, on Slurm."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import time

EXPECTED_SHA = '82981de62dfa2f2ff5eaedd273e38d8b20bf36034dc39110c104908599470d25'
POLICIES = ('operational', 'mathematical', 'origin')
TOL = 1e-8


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def checked(ref):
    if pin(ref['path']) != ref:
        raise ValueError('changed artifact: '+ref['path'])
    return json.loads(Path(ref['path']).read_text())


def dump(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def finite(value):
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError('nonfinite or absent available scalar')
    return value


def closure(residual, label):
    if abs(finite(residual)) > TOL:
        raise ValueError(f'{label}: residual {residual} exceeds {TOL} kcal/mol')
    return residual


def decision(score, bands, expected):
    label = 'Ca' if score <= bands['Ca_max'] else 'La' if score >= bands['La_min'] else 'inconclusive'
    return {'decision': label, 'outcome': 'inconclusive' if label == 'inconclusive'
            else 'correct' if label == expected else 'wrong'}


def diagnose(source, policy, old_calibration, new_calibration):
    out = {k: source[k] for k in ('record_id', 'family', 'source_case_id', 'selection_id',
                                  'biological_group', 'role', 'expected_class', 'source_conditioning_metal')}
    old = source['baseline']['pool'][policy]
    new = source['pool'][policy]
    old_bands = old_calibration[policy]['bands']
    new_bands = new_calibration[policy]['bands']
    out.update(policy=policy, status='unavailable', reason=None,
               old_decision=source['baseline']['decisions'][policy]['decision'],
               old_outcome=source['baseline']['decisions'][policy]['outcome'],
               new_decision=source['decisions'][policy]['decision'],
               new_outcome=source['decisions'][policy]['outcome'],
               old_R_kcal_mol=old['R_model_kcal_mol'] if old else None,
               new_R_kcal_mol=new['R_model_kcal_mol'] if new else None,
               fixed_old_selection_R_kcal_mol=None,
               fixed_old_selection_new_band_decision=None,
               fixed_old_selection_new_band_outcome=None)
    if old is None or new is None:
        out['reason'] = 'baseline_policy_unavailable' if old is None else 'candidate_policy_unavailable'
        return out
    if old_bands is None or new_bands is None:
        out['reason'] = 'frozen_calibration_unavailable'
        return out
    expected = source['expected_class']
    if expected not in ('Ca', 'La'):
        raise ValueError('unsupported expected class')
    old_R, new_R = finite(out['old_R_kcal_mol']), finite(out['new_R_kcal_mol'])
    if decision(old_R, old_bands, expected) != source['baseline']['decisions'][policy]:
        raise ValueError('baseline decision does not reproduce')
    if decision(new_R, new_bands, expected) != source['decisions'][policy]:
        raise ValueError('candidate decision does not reproduce')
    works = {}
    for metal in ('Ca', 'La'):
        o, n = old[metal]['candidate'], new[metal]['candidate']
        available = new[metal]['all_works_kcal_mol']
        out[metal+'_old_candidate'], out[metal+'_new_candidate'] = o, n
        if o not in available or n not in available:
            out['reason'] = 'required_candidate_work_unavailable'
            return out
        wo, wn = finite(available[o]), finite(available[n])
        closure(wn-finite(new[metal]['work_kcal_mol']), metal+' selected-work')
        works[metal] = (wo, wn)
        out[metal+'_new_model_work_at_old_selection_kcal_mol'] = wo
        out[metal+'_new_model_work_at_new_selection_kcal_mol'] = wn
        out[metal+'_reselection_energy_change_kcal_mol'] = wn-wo
    fixed_R = new_R+(works['Ca'][0]-works['Ca'][1])-(works['La'][0]-works['La'][1])
    new_origin = source['pool']['origin']
    old_origin = source['baseline']['pool']['origin']
    if new_origin is None or old_origin is None:
        raise ValueError('available policy lacks required paired origin')
    independent_fixed = finite(new_origin['R_model_kcal_mol'])+works['Ca'][0]-works['La'][0]
    selected_new = finite(new_origin['R_model_kcal_mol'])+works['Ca'][1]-works['La'][1]
    selected_old = (finite(old_origin['R_model_kcal_mol'])
                    +finite(old['Ca']['work_kcal_mol'])-finite(old['La']['work_kcal_mol']))
    model_delta, selection_delta = fixed_R-old_R, new_R-fixed_R
    raw_delta = new_R-old_R
    boundary = 'Ca_max' if expected == 'Ca' else 'La_min'
    sign = -1 if expected == 'Ca' else 1
    old_boundary, new_boundary = finite(old_bands[boundary]), finite(new_bands[boundary])
    band_shift = new_boundary-old_boundary
    old_margin, new_margin = sign*(old_R-old_boundary), sign*(new_R-new_boundary)
    margin_delta = new_margin-old_margin
    fixed_call = decision(fixed_R, new_bands, expected)
    out.update(status='available', reason=None, fixed_old_selection_R_kcal_mol=fixed_R,
        fixed_old_selection_new_band_decision=fixed_call['decision'],
        fixed_old_selection_new_band_outcome=fixed_call['outcome'],
        fixed_old_selection_old_band_decision=decision(fixed_R, old_bands, expected)['decision'],
        raw_R_change_kcal_mol=raw_delta, model_effect_at_old_selection_kcal_mol=model_delta,
        geometry_reselection_effect_kcal_mol=selection_delta,
        expected_class_boundary=boundary, old_boundary_kcal_mol=old_boundary,
        new_boundary_kcal_mol=new_boundary, band_shift_kcal_mol=band_shift,
        old_expected_class_margin_kcal_mol=old_margin,
        fixed_selection_new_band_expected_class_margin_kcal_mol=sign*(fixed_R-new_boundary),
        new_expected_class_margin_kcal_mol=new_margin,
        expected_class_margin_change_kcal_mol=margin_delta,
        model_margin_contribution_kcal_mol=sign*model_delta,
        reselection_margin_contribution_kcal_mol=sign*selection_delta,
        band_margin_contribution_kcal_mol=-sign*band_shift,
        independent_fixed_R_residual_kcal_mol=closure(fixed_R-independent_fixed, 'fixed R from origin'),
        selected_new_R_residual_kcal_mol=closure(new_R-selected_new, 'new selected R from origin'),
        selected_old_R_residual_kcal_mol=closure(old_R-selected_old, 'old selected R from origin'),
        raw_change_closure_kcal_mol=closure(raw_delta-model_delta-selection_delta, 'raw change'),
        margin_change_closure_kcal_mol=closure(margin_delta-sign*model_delta-sign*selection_delta
                                             +sign*band_shift, 'margin change'))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--collection', required=True)
    parser.add_argument('--sha', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Real Slurm allocation required for saved-matrix scientific analysis')
    if args.sha != EXPECTED_SHA:
        raise SystemExit('This declared analysis pins the completed transfer collection')
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    start, cpu_start = time.monotonic(), time.process_time()
    receipt = {'status': 'failed', 'job_id': os.environ['SLURM_JOB_ID'], 'driver': pin(__file__),
               'plan': pin(Path(__file__).with_name('PLAN.md')), 'new_molecular_calls': 0}
    try:
        collection_pin = {'path': str(Path(args.collection).resolve()), 'sha256': args.sha}
        data = checked(collection_pin)
        if data['status'] != 'complete' or data['stage'] != 'transfer':
            raise ValueError('requires completed transfer scorer output')
        inventory = checked(data['inventory'])
        canonical = checked(data['canonical_calibration_receipt'])
        baseline = checked(inventory['baseline_canonical'])
        if (canonical['calibrations'] != data['calibrations']
                or baseline['calibration'] != data['baseline_calibrations']
                or canonical['inventory'] != data['inventory']
                or canonical['native_manifest'] != data['native_manifest']
                or canonical['scorer'] != data['scorer']):
            raise ValueError('frozen calibration identity differs')
        if len({r['record_id'] for r in data['rows']}) != len(data['rows']):
            raise ValueError('duplicate declared record')
        rows = [diagnose(r, p, data['baseline_calibrations'], data['calibrations'])
                for r in data['rows'] for p in POLICIES]
        summary = {}
        for family in sorted({r['family'] for r in rows}):
            summary[family] = {p: dict(Counter(r['status'] for r in rows
                               if r['family'] == family and r['policy'] == p)) for p in POLICIES}
        result = {'status': 'complete', 'protocol': 'saved_candidate_selection_diagnosis_v1',
            'collection': collection_pin, 'candidate_canonical': data['canonical_calibration_receipt'],
            'baseline_canonical': inventory['baseline_canonical'], 'inventory': data['inventory'],
            'units': 'kcal/mol; already converted saved scalars', 'closure_tolerance_kcal_mol': TOL,
            'declared_records': len(data['rows']), 'declared_policy_rows': len(rows),
            'rows': rows, 'availability_by_family_policy': summary,
            'cohorts': inventory['cohorts'],
            'current_ca125_annotated_overlay': data['current_ca125_annotated_overlay'],
            'new_molecular_calls': 0, 'new_calibration': False, 'production_changed': False,
            'interpretation': 'Arithmetic counterfactual decomposition, not a new classifier or causal electronic diagnosis.'}
        dump(output/'RESULT.json', result)
        fields = list(dict.fromkeys(k for r in rows for k in r))
        with (output/'ROWS.tsv').open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        lines = ['# Saved-matrix selection diagnosis', '',
            'Completed without new molecular calls or recalibration. All declared records and policies remain visible.', '',
            'ROWS.tsv separates energy-model change at old selections, geometry reselection and frozen-band movement.',
            'Counterfactual calls use the explicitly stated bands and are not a proposed production decision rule.', '',
            '| Family | Policy | Available | Unavailable |', '|---|---|---:|---:|']
        for family, policies in summary.items():
            for policy, counts in policies.items():
                lines.append(f"| {family} | {policy} | {counts.get('available', 0)} | {counts.get('unavailable', 0)} |")
        lines += ['', 'Missing TSV measurements are blank, not zero. Positive expected-class margin reaches that class band.',
                  'Cohorts overlap; consult the preserved membership and 125-source recovery mapping before summarizing.']
        (output/'REPORT.md').write_text('\n'.join(lines)+'\n')
        receipt.update(status='complete', result=pin(output/'RESULT.json'), table=pin(output/'ROWS.tsv'))
    except BaseException as exc:
        receipt['failure_reason'] = type(exc).__name__+': '+str(exc)
        raise
    finally:
        receipt.update(elapsed_seconds=time.monotonic()-start, CPU_seconds=time.process_time()-cpu_start,
                       maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        dump(output/'EXECUTION.json', receipt)


if __name__ == '__main__':
    main()
