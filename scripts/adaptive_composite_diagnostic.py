"""Normalize existing seven-context analytic derivatives; no solver calls."""
from __future__ import annotations
import argparse
import json
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from adaptive_force_diagnostic import preview
from mace_site_kinematics import Kinematics


def analyze(source, agreement, output):
    result = read_json(source); design = read_json(verify(result['design'])); rows = []
    for case in result['cases']:
        cid = case['case_id']; geometry = {}
        for metal in ('Ca', 'La'):
            data = read_json(verify(design['maps'][cid+'__'+metal]))['context']
            kin = Kinematics(data); physical, _, j, _ = kin.evaluate(np.zeros(len(kin.modes)))
            vectors = j[:, kin.heavy].reshape(len(kin.modes), -1)
            lengths = np.linalg.norm(vectors, axis=1)
            if np.any(lengths <= 1e-12): raise InvalidArtifact('zero physical mode')
            if [m['id'] for m in kin.modes] != case['mode_ids']: raise InvalidArtifact('mode identities differ')
            geometry[metal] = (vectors/lengths[:, None], lengths, kin.modes, physical[kin.heavy])
        ca, la = geometry['Ca'], geometry['La']
        if not all(np.allclose(ca[i], la[i], atol=1e-12, rtol=0) for i in (0, 1, 3)):
            raise InvalidArtifact('endpoint physical heavy-atom measures differ')
        angular = [i for i,m in enumerate(ca[2]) if m['kind'] != 'metal_translation']
        methods = {}
        for model in ('vacuum_MACE', 'composite'):
            loads = {z: np.asarray(case['endpoints'][z][model])/geometry[z][1] for z in ('Ca', 'La')}
            delta = loads['Ca']-loads['La']
            if not np.allclose(delta, np.asarray(case['Ca_minus_La_projected_derivatives'][model])/ca[1], atol=1e-10, rtol=0):
                raise InvalidArtifact('differential derivative does not replay')
            highest = sorted(angular, key=lambda i: (-abs(delta[i]), ca[2][i]['id']))[0]
            methods[model] = {'preview': preview(ca[2], ca[0], loads['Ca'], loads['La']),
                              'highest_angular_differential': ca[2][highest]['id'],
                              'highest_angular_delta_kcal_mol_A': float(delta[highest]),
                              'normalized_derivatives': {z: a.tolist() for z,a in loads.items()},
                              'normalized_differential': delta.tolist()}
        a,b = [set(s['id'] for s in methods[key]['preview']['selected']) for key in ('vacuum_MACE','composite')]
        rows.append({'case_id': cid, 'methods': methods, 'four_mode_selection_overlap': len(a & b),
                     'mode_ids': case['mode_ids'], 'mode_kinds': [m['kind'] for m in ca[2]],
                     'state_scope': 'common physical heavy donor measure; original endpoint-specific fixed H retained'})
    out = {'source': record(source), 'agreement': record(agreement), 'implementation': record(__file__),
           'selector_implementation': record(__file__.replace('adaptive_composite_diagnostic.py','adaptive_force_diagnostic.py')),
           'cases': rows, 'denominator': len(rows), 'new_molecular_calls': 0,
           'affinity_score': None, 'labels_used_for_selection': False}
    write_new(output, out)
    return [{'case_id': r['case_id'], 'four_mode_selection_overlap': r['four_mode_selection_overlap'],
             **{k: v['highest_angular_differential'] for k,v in r['methods'].items()}} for r in rows]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('source','agreement','output'): p.add_argument('--'+k,required=True)
    print(json.dumps(analyze(**vars(p.parse_args())),indent=2))


if __name__ == '__main__': main()
