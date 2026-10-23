#!/usr/bin/env python3
"""Saved conditional z-susceptibility comparison; no molecular execution."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time

STATES = ('La_dimer', 'La_trimer', 'Dy_dimer', 'Dy_trimer')
STEPS = {'h': .001, 'half': .0005}
LIMITATIONS = [
    'The shared numerical consistency limits and printed-digit bounds are not physical equivalence thresholds or confidence intervals.',
    'This is one lab-z susceptibility column at fixed consumed W9 far-pair geometry. The added-water increment is trimer minus dimer, not alpha3 or full-shell B.',
    'The provisional frozen-f PBE0/common-parent ghost basis, SCF stopping evidence, finite-field truncation and original baseline model limitations remain.',
    'Curvature supplies response information beyond the already measured zero-field dipole. Agreement in this mode cannot identify a unique cause of the zero-field energy or dipole discrepancy.',
    'No fitting, target changes, native calls, field-force qualification or physical pass/fail threshold are introduced.'
]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(value):
    require(isinstance(value, dict) and pin(value['path']) == value, 'Missing or changed pinned input')
    return Path(value['path'])


def read(value):
    return json.loads(verify(value).read_text())


def module(value, name):
    spec = importlib.util.spec_from_file_location(name, verify(value))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def compare(cfg):
    require(cfg['implementation'] == pin(__file__), 'Comparator source differs')
    verify(cfg['declaration']); verify(cfg['specification']); verify(cfg['source_inventory'])
    m = read(cfg['model_manifest']); model = read(cfg['model_result']); review = read(cfg['model_review'])
    require(review['status'] == 'passed_independent_actual_model_susceptibility_review' and
            review['input_result'] == cfg['model_result'] and review['manifest'] == cfg['model_manifest'] and
            review['implementation'] == cfg['model_implementation'] and
            review['response_checks_passed'] == review['response_checks_declared'] == 36 and
            review['full_native_denominator'] == 19 and all(x['passed'] is True for x in review['checks']),
            'Exact independent model admission required')
    verify(review['review_source'])
    require(model['status'] == 'conditional_model_susceptibility_complete' and model['error'] is None and
            model['global_stop'] is False and model['controls_passed'] is True and
            model['native_program_invocations'] == 19 and model['manifest'] == cfg['model_manifest'] and
            model['implementation'] == m['implementation'] == cfg['model_implementation'], 'Model result incomplete or changed')
    require(m['sources']['shared_specification'] == cfg['specification'] and
            m['sources']['shared_inventory'] == cfg['source_inventory'], 'Model shared inputs differ')
    model_api = module(cfg['model_implementation'], 'conditional_comparison_model')
    ms = model_api.summarize(m, model['rows'])
    require(ms == model['summary'] and ms['numerical_checks_passed'] is True, 'Model estimator replay differs')

    ref_api = module(cfg['reference_implementation'], 'conditional_comparison_reference')
    p = ref_api.load(verify(cfg['reference_panel']))
    ref = read(cfg['reference_collection'])
    ids = [c['cell_id'] for c in p['cells']]
    require(p['adapter'] == cfg['reference_implementation'] and p['specification'] == cfg['specification'] and
            p['source_inventory'] == cfg['source_inventory'], 'Reference shared inputs differ')
    require(ref['status'] == 'conditional_susceptibility_reference_complete' and ref['stage'] == 'full' and
            ref['implementation'] == cfg['reference_implementation'] and ref['panel'] == cfg['reference_panel'] and
            ref['specification'] == cfg['specification'] and ref['complete'] == ref['declared'] == 16 and
            ref['cell_ids'] == ids and len(ids) == len(set(ids)) == 16 and set(ref['rows']) == set(ids) and
            ref['numerical_checks_passed'] is True and ref['new_native_calls'] == 0 and
            ref['native_programs_in_admitted_rows'] == 16, 'Full qualified sixteen-field collector required')
    inventory = read(ref['inventory'])
    require(set(inventory) == set(ids), 'Reference inventory denominator differs')
    for cell in p['cells']:
        cid = cell['cell_id']; row = ref['rows'][cid]
        require(row['status'] == 'complete' and row['error'] is None and row['returncode'] == 0 and
                row['cell'] == cell and row['panel'] == cfg['reference_panel'] and
                row['adapter'] == cfg['reference_implementation'] and row['native_program_invocations'] == 1,
                'Reference collector row lineage differs: '+cid)
        inv = inventory[cid]
        source = verify(inv) if isinstance(inv, dict) else Path(inv)
        require(pin(source) == row['result'], 'Reference inventory/result join differs')
        require(dict(read(row['result']), result=row['result']) == row, 'Collector/native result lineage differs')
        for key in ('input', 'output', 'start', 'release'):
            verify(row[key])
        verify(row['properties']['property_file'])
        require(row['properties']['origin_Bohr'] == [0., 0., 0.], 'Reference moment origin differs')
    rs = ref_api.score(ref['rows'], p, 'full')
    require(rs == ref['response'] and rs['numerical_checks_passed'] is True, 'Reference estimator replay differs')

    joins = {}
    for state in STATES:
        row = m['rows'][state]; base = row['quantum_zero_field_id']
        require(row['quantum_zero_field_source'] == p['zero_field_records'][base]['result'], 'Zero-field authority differs')
        positions = [[float(x) for x in line.split()[2:5]] for line in verify(row['input']['xyz']).read_text().splitlines()[1:] if line.strip()]
        for cell in (c for c in p['cells'] if c['zero_id'] == base):
            require(row['parent_atom_indices'] == cell['real_parent_atom_indices'] and
                    positions == [cell['positions_A'][i] for i in row['parent_atom_indices']] and
                    row['charge'] == cell['charge'] == 3, 'Exact real geometry/mask/charge differs')
            job = m['jobs'][state+'__'+cell['field_tag']]
            require(job['field_au'] == [0., 0., cell['field_au']], 'Matched applied field differs')
        joins[state] = base

    table = []
    groups = [(s, 'state', ms['states'][s], rs['per_state'][joins[s]]) for s in STATES]
    groups += [(s, 'added_water' if s != 'Dy_minus_La' else 'paired', ms['increments'][s],
                rs['paired'] if s == 'Dy_minus_La' else rs['increments'][s]) for s in ('La', 'Dy', 'Dy_minus_La')]
    for name, kind, a, b in groups:
        for step, field in STEPS.items():
            av, bv = a['steps'][step], b['estimates'][step]
            values = {'energy': (av['chi_E_a0_cubed'], bv['chi_energy'], av['energy_curvature_print_bound_a0_cubed'], bv['roundoff_bounds']['chi_energy']),
                      'moment': (av['chi_mu_column_a0_cubed'][2], bv['chi_moment'], av['moment_slope_print_bound_a0_cubed'], bv['roundoff_bounds']['chi_moment'])}
            for estimator, (mv, rv, mb, rb) in values.items():
                require(all(math.isfinite(x) for x in (mv, rv, mb, rb)) and min(mb, rb) >= 0, 'Nonfinite susceptibility or invalid rounding bound')
                table.append({'quantity': name, 'kind': kind, 'field_size_au': field, 'step': step, 'estimator': estimator,
                              'model_a0_cubed': mv, 'reference_a0_cubed': rv, 'model_minus_reference_a0_cubed': mv-rv,
                              'model_print_bound_a0_cubed': mb, 'reference_print_bound_a0_cubed': rb,
                              'difference_print_bound_a0_cubed': mb+rb})
    return {'status': 'conditional_susceptibility_comparison_complete', 'configuration': cfg,
            'model_denominator': 19, 'reference_denominator': 16, 'matched_state_ids': joins,
            'table': table, 'model_response': ms, 'reference_response': rs,
            'model_controls': {k: model['rows'][k]['control'] for k in m['controls']},
            'model_induction': {k: v['moments']['induction'] for k, v in model['rows'].items()},
            'reference_native_evidence': {k: {'result': v['result'], 'output': v['output'],
                'property_file': v['properties']['property_file'], 'parsed': v['parsed']} for k, v in ref['rows'].items()},
            'reference_zero_native_evidence': {k: v['native'] for k, v in p['zero_replays'].items()},
            'numerical_checks_passed': True, 'physical_equivalence_threshold': None,
            'limitations': LIMITATIONS, 'new_native_calls': 0,
            'source_costs': {'model': model['cost'], 'reference_collection': ref['cost'],
                           'reference_admitted_allocated_CPU_seconds_proxy': ref['admitted_native_allocated_CPU_seconds_proxy']}}


def report(result):
    lines = ['# Conditional susceptibility: model and reference', '',
             'All 19 model calls and 16 reference fields are admitted. Both existing numerical panels pass. '
             'Values and model minus reference differences are in a0³; fields are in atomic units.', '',
             '| Quantity | Field | Estimator | Model | Reference | Model − reference | Sum of print bounds |',
             '|---|---:|---|---:|---:|---:|---:|']
    for row in result['table']:
        lines.append('| {quantity} | {field_size_au:g} | {estimator} | {model_a0_cubed:.8f} | {reference_a0_cubed:.8f} | {model_minus_reference_a0_cubed:+.8f} | {difference_print_bound_a0_cubed:.8g} |'.format(**row))
    lines += ['', '## Interpretation limits', ''] + ['- '+x for x in LIMITATIONS]
    lines += ['', 'RESULT.json retains both complete estimator/gate records, odd-slope values and their rounding bounds, '
              'all native SCF evidence, model induction records and exact input/result pins. New native calls: 0.', '']
    return '\n'.join(lines)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--config', type=Path, required=True); ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args(); start, cpu = time.monotonic(), time.process_time()
    cfg = json.loads(args.config.read_text()); args.out.mkdir(parents=True, exist_ok=False)
    try:
        result = compare(cfg); result['configuration_pin'] = pin(args.config); text = report(result); exitcode = 0
    except Exception as exc:
        result = {'status': 'conditional_susceptibility_comparison_unavailable', 'configuration_pin': pin(args.config),
                  'error': repr(exc), 'table': None, 'new_native_calls': 0}
        text = '# Conditional susceptibility comparison unavailable\n\n'+repr(exc)+'\n'; exitcode = 1
    result['cost'] = {'wall_seconds': time.monotonic()-start, 'CPU_seconds': time.process_time()-cpu, 'local_threads': 1}
    (args.out/'RESULT.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    (args.out/'REPORT.md').write_text(text)
    print(json.dumps({'result': pin(args.out/'RESULT.json'), 'status': result['status']}))
    return exitcode


if __name__ == '__main__':
    raise SystemExit(main())
