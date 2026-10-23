#!/usr/bin/env python3
"""No native calls: validate the new Urey correction and fail-closed comparisons."""
import argparse
import copy
import importlib.util
import json
import math
from pathlib import Path
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--si-xyz', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('mapping_run', Path(__file__).with_name('RUN.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    start = time.monotonic()
    rows = mod.xyz(args.si_xyz)
    for row in rows:
        row['type'] = {36: 402, 37: 403, 500: 900, 501: 901}[row['type']]
    correction = mod.urey_delta(rows)
    checks = []

    def check(name, passed):
        if not passed:
            raise ValueError(name)
        checks.append(name)

    def direct_energy(atoms, length):
        return sum(-7.6*(math.dist(atoms[i]['position'], atoms[j]['position'])-length)**2
                   for i, j in correction['hydrogen_pairs_zero_based'])

    i, _ = correction['hydrogen_pairs_zero_based'][0]
    step = 1e-5
    for axis in range(3):
        plus, minus = copy.deepcopy(rows), copy.deepcopy(rows)
        plus[i]['position'][axis] += step
        minus[i]['position'][axis] -= step
        fd = ((direct_energy(plus, 1.5537)-direct_energy(plus, 1.5337))-
              (direct_energy(minus, 1.5537)-direct_energy(minus, 1.5337)))/(2*step)
        check('independent_energy_difference_derivative_axis_'+str(axis),
              abs(fd-correction['gradient_kcal_mol_A'][i][axis]) < 1e-7)
    check('zero_net_force', all(abs(sum(v[a] for v in correction['gradient_kcal_mol_A'])) < 1e-12 for a in range(3)))
    shifted = copy.deepcopy(rows)
    for row in shifted:
        row['position'] = [x+y for x, y in zip(row['position'], [3.1, -7.2, 0.5])]
    check('translation_invariance', abs(mod.urey_delta(shifted)['energy_kcal_mol']-correction['energy_kcal_mol']) < 1e-12)
    reference = {'total_energy_kcal_mol': 0, 'gradient_kcal_mol_A': [[0.0]*3 for _ in rows],
                 'components': {'Urey-Bradley': {'energy_kcal_mol': 0, 'interactions': 8},
                                'Polarization': {'energy_kcal_mol': 0, 'interactions': 276}}}
    actual = copy.deepcopy(reference)
    actual['total_energy_kcal_mol'] = correction['energy_kcal_mol']
    actual['gradient_kcal_mol_A'] = correction['gradient_kcal_mol_A']
    actual['components']['Urey-Bradley']['energy_kcal_mol'] = correction['energy_kcal_mol']
    check('known_correction_passes', mod.comparison(actual, reference, correction)['passed'])
    for name, field in [('changed_energy_fails', 'total_energy_kcal_mol'), ('changed_force_fails', 'gradient_kcal_mol_A')]:
        bad = copy.deepcopy(actual)
        if field == 'total_energy_kcal_mol':
            bad[field] += .001
        else:
            bad[field][0][0] += .01
        check(name, not mod.comparison(bad, reference, correction)['passed'])
    bad = copy.deepcopy(actual)
    bad['components']['Polarization']['energy_kcal_mol'] += .001
    check('changed_component_fails', not mod.comparison(bad, reference, correction)['passed'])
    for name, mutate in [('missing_component_rejected', lambda d: d['components'].pop('Urey-Bradley')),
                         ('missing_gradient_rejected', lambda d: d['gradient_kcal_mol_A'].pop()),
                         ('wrong_interaction_count_rejected', lambda d: d['components']['Urey-Bradley'].update(interactions=7))]:
        bad = copy.deepcopy(actual)
        mutate(bad)
        try:
            mod.comparison(bad, reference, correction)
        except ValueError:
            check(name, True)
        else:
            check(name, False)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    mod.save(out/'RESULT.json', {'status': 'passed', 'checks': checks, 'count': len(checks),
             'native_calls': 0, 'source': mod.pin(mod.__file__), 'test_source': mod.pin(__file__),
             'input': mod.pin(args.si_xyz), 'elapsed_seconds': time.monotonic()-start})
    print(json.dumps({'status': 'passed', 'checks': len(checks), 'result': mod.pin(out/'RESULT.json')}))


if __name__ == '__main__':
    main()
