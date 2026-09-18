"""Read-only decomposition of actual saved endpoint scores; no fitting or inference."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from mace_hybrid import EV_TO_KCAL

COMPONENTS = ('interaction_energy', 'electrostatic_energy', 'electron_energy')


def cohort(path, canonical):
    d = read_json(path)
    if d['status'] != 'complete':
        raise InvalidArtifact('readout diagnostic requires the completed declared cohort')
    m = read_json(verify(d['collection']['manifest'])); prep = read_json(verify(m['preparation']))
    tasks = {t['task_id']: t for t in m['tasks']}; rows = d['collection']['rows']; output = {}
    for name, s in d['scores'].items():
        if canonical and s['evaluation_role'] != 'calibration':
            continue
        ca_id = name + '_Ca_primary' if canonical else s['Ca_task']
        la_id = name + '_La_primary' if canonical else s['La_task']
        ca, la = rows[ca_id], rows[la_id]
        if ca['status'] != 'computed' or la['status'] != 'computed':
            raise InvalidArtifact('missing declared real endpoint')
        group_row = read_json(verify(prep['cases'][name])); group = group_row['primary']
        source = read_json(verify(group_row['case']['preparation']))
        total = (ca['energy_eV'] - la['energy_eV']) * EV_TO_KCAL
        if total != s['R_model_kcal']:
            raise InvalidArtifact('score/endpoint algebra differs')
        parts = {key: (ca['energy_components_eV'][key] - la['energy_components_eV'][key]) * EV_TO_KCAL for key in COMPONENTS}
        if not np.isfinite(list(parts.values())).all():
            raise InvalidArtifact('nonfinite saved energy readout')
        parts['unassigned_remainder'] = total - sum(parts.values()); parts['total'] = total
        output[name] = {'contrasts_model_kcal': parts,
                        'full_Ca_charge_e': tasks[ca_id]['charge'],
                        'atom_count': len(source['physical_atoms']),
                        'selected_group_Ca_charge_e': group['endpoint_group_charges_e']['Ca'][group['selected_group_index']],
                        'baseline_S_kcal_mol': s['baseline']['published_S_kcal_mol'] if canonical else None,
                        'expected_class': s['expected_class'] if canonical else None,
                        'source_preparation': group_row['case']['preparation'],
                        'Ca_task': ca_id, 'La_task': la_id}
    if len(output) != (25 if canonical else 34):
        raise InvalidArtifact('declared readout cohort changed')
    residual = [r['contrasts_model_kcal']['unassigned_remainder'] for r in output.values()]
    return {'source_report': record(path), 'rows': output,
            'unassigned_remainder_range_model_kcal': max(residual) - min(residual),
            'remainder_common_within_0_01_model_kcal': max(residual) - min(residual) <= .01}


def summarize(canonical):
    rows = list(canonical['rows'].values()); result = {}
    for name in (*COMPONENTS, 'unassigned_remainder', 'total'):
        values = np.array([r['contrasts_model_kcal'][name] for r in rows])
        ca = [r['contrasts_model_kcal'][name] for r in rows if r['expected_class'] == 'Ca']
        la = [r['contrasts_model_kcal'][name] for r in rows if r['expected_class'] == 'La']
        correlations = {}
        for field in ('full_Ca_charge_e', 'atom_count', 'selected_group_Ca_charge_e', 'baseline_S_kcal_mol'):
            other = np.array([r[field] for r in rows])
            resolved = np.ptp(values) > (.01 if name == 'unassigned_remainder' else 0)
            correlations[field] = float(np.corrcoef(values, other)[0, 1]) if resolved and np.ptp(other) > 0 else None
        result[name] = {'Pearson_descriptive_only': correlations, 'min_La_minus_max_Ca_model_kcal': min(la) - max(ca),
                        'Ca_range_model_kcal': [min(ca), max(ca)], 'La_range_model_kcal': [min(la), max(la)],
                        'calibrated_classifier': None}
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('canonical-report', 'transfer-report', 'agreement', 'output'):
        p.add_argument('--' + name, required=True)
    a = p.parse_args()
    canonical = cohort(a.canonical_report, True); transfer = cohort(a.transfer_report, False)
    out = Path(a.output).resolve(); out.mkdir(parents=True, exist_ok=False)
    src = out / 'audit_readouts.py'; src.write_bytes(Path(__file__).read_bytes())
    result = {'status': 'complete', 'agreement': record(a.agreement), 'implementation': record(src),
              'canonical': canonical, 'transfer': transfer, 'canonical_descriptive_summary': summarize(canonical),
              'new_energy_evaluations': 0, 'fitted_parameters': 0, 'new_classifier': None,
              'baseline_changed': False, 'interpretation': 'algebraic readout diagnostic; no unique causal or independent-validation claim'}
    write_new(out / 'result.json', result)
    print(json.dumps(result['canonical_descriptive_summary'], indent=2))


if __name__ == '__main__':
    main()
