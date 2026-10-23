#!/usr/bin/env python3
"""Saved-token decomposition and finite source joins; no native/FFT/factor calls."""
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / 'workspaces').is_dir())
J = Path('explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1')
W = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021' / J
S2 = W / 'periodic_native_qualification_v1/stage2_certificate_v1'
NATIVE = W / 'periodic_stability_definition_v1/pme_export_v1/source_v3/source'

def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def checked(p):
    assert pin(p['path']) == p, p
    return Path(p['path'])

def read(p):
    return json.loads(Path(p).read_text())

start, cpu = time.monotonic(), time.process_time()
authority = pin(S2 / 'full_collection_v1/execution/RESULT.json')
result = read(authority['path'])
assert result['complete'] and not result['numerical_checks_passed'] and result['failures'] == []
assert result['counts']['admitted_gradients'] == 36
assert result['counts']['admitted_new_certificates'] == 32
rows = result['rows']
output = {'status': 'completed_saved_periodic_response_derivative_diagnosis',
          'implementation': pin(__file__), 'failed_full_collection': authority,
          'candidate': result['candidate'], 'configuration': result['configuration'],
          'build': result['build'], 'source_pins': {}, 'rows': {}, 'directions': {},
          'new_native_calls': 0, 'new_FFT_replays': 0, 'new_factorizations': 0,
          'new_solves': 0, 'author_control_reruns': 0}
for state, row in rows.items():
    checked(row['result'])
    raw = checked(row['gradient_output'])
    with raw.open() as stream:
        first = [next(stream).strip() for _ in range(5)]
    assert first[2].split() == ['TOTAL_ENERGY', row['energy_token']]
    assert first[3].split() == ['POLARIZATION_ENERGY', row['polarization_token']]
    output['rows'][state] = {k: row[k] for k in ('result', 'gradient_output', 'energy_token', 'polarization_token')}

def energy(state, key):
    return Decimal(rows[state][key])

for direction, actual in result['numerical']['directions'].items():
    values = []
    for step in ('1e-04', '5e-05'):
        minus, plus = direction + '__' + step + '__minus', direction + '__' + step + '__plus'
        denominator = 2 * Decimal(step)
        total = (energy(plus, 'energy_token') - energy(minus, 'energy_token')) / denominator
        polar = (energy(plus, 'polarization_token') - energy(minus, 'polarization_token')) / denominator
        remainder = total - polar
        values.append({'half_step_A': step, 'FD_total': str(total), 'FD_polarization': str(polar),
                       'FD_nonpolar_remainder': str(remainder)})
    differences = {k: str(Decimal(values[1][k]) - Decimal(values[0][k]))
                   for k in ('FD_total', 'FD_polarization', 'FD_nonpolar_remainder')}
    assert Decimal(differences['FD_total']) == Decimal(differences['FD_polarization']) + Decimal(differences['FD_nonpolar_remainder'])
    output['directions'][direction] = {'units': 'kcal/mol/angstrom', 'steps': values,
        'half_minus_full_step_slope_difference': differences,
        'original_numerical_checks_unchanged': actual}

active_minus_baseline = {}
for key in ('energy_token', 'polarization_token'):
    active_minus_baseline[key] = str(energy('central', key) - energy('BYPASS_OLD_ABSENT', key))
assert active_minus_baseline['energy_token'] == active_minus_baseline['polarization_token']
output['central_active_minus_baseline'] = active_minus_baseline

build = read(checked(result['build']))
manifest = read(checked(build['source_manifest']))
clone = checked(manifest['files']['induce_qualification.f'])
output['source_manifest'] = build['source_manifest']
output['source_pins']['actual_induce_clone'] = pin(clone)
output['source_pins']['actual_gradient_reporter'] = manifest['files']['GRADIENT.f90']
checked(manifest['files']['GRADIENT.f90'])
text = (NATIVE / 'induce.f').read_text()
clone_text = clone.read_text()
def routine(text, name):
    start = text.index('      subroutine ' + name + '\n')
    end = text.index('\n      end\n', start) + len('\n      end\n')
    return text[start:end]
for name in ('periodicresponsegeometry', 'periodicresponsegradient', 'localresponsegradient'):
    assert routine(text, name) == routine(clone_text, name), name
vacuum = routine(text, 'localresponsegradient')
periodic = routine(text, 'periodicresponsegradient')
assert vacuum[vacuum.index('      fac='):].replace('type(lrm)-499', 'type(lrm)-899') == periodic[periodic.index('      fac='):]
output['unchanged_executed_response_bodies'] = True
output['periodic_gradient_arithmetic_matches_qualified_vacuum_after_namespace_mapping'] = True
for name in ('induce.f', 'epolar1.f', 'epolar3.f', 'gradient.f'):
    output['source_pins'][name] = pin(NATIVE / name)
output['central_polarization_contraction'] = result['numerical']['central_polarization_contraction']
output['central_shared_h'] = result['numerical']['shared_h_checks']['central']
output['central_PCG_pre_peek'] = rows['central']['PCG_pre_peek']
output['interpretation'] = {
    'supported': 'Large step sensitivity exists in Etotal-Epolar, independently of the custom h/H energy. A custom-only omitted derivative cannot explain all failures.',
    'not_established': 'This decomposition does not identify a particular baseline cutoff/PME/coordinate defect, prove asymptotic convergence, or exclude an additional response-dependent discrepancy.',
    'gradient_component_limit': 'Only total gradient/virial and polarization energy were reported. No saved analytic polarization/nonpolar force or virial split is inferred.',
    'central_energy_limit': 'Passing stationary contraction checks energy convention only; it does not independently qualify coordinate or strain derivatives.'}
output['prospective_root_control'] = {
    'directions': ['whole_box', 'strain_1'], 'half_steps_A': ['1e-04', '5e-05'],
    'signs': ['minus', 'plus'], 'new_gradient_programs': 8,
    'central_reuse': rows['BYPASS_OLD_ABSENT']['result'],
    'change': 'Only remove the two response records; retain exact displaced geometries/cells and native numerical settings.',
    'no_new_initializers_or_certificates': True,
    'comparisons': ['Original total FD/analytic/two-step gates for baseline.',
                    'Same-step active-minus-baseline energies and analytic slopes, to measure the derivative residual of the response-induced increment.'],
    'limitation': 'A response-increment discrepancy includes response-modified native induction/PME errors and is not by itself proof of a local h/H algebra bug.',
    'executor_owner': 'Root and Planck; this audit creates no execution adapter or release.'}
output['cost'] = {'wall_seconds': time.monotonic() - start, 'CPU_seconds': time.process_time() - cpu, 'threads': 1}
with (HERE / 'AUDIT.json').open('x') as stream:
    stream.write(json.dumps(output, indent=2, allow_nan=False) + '\n')
print(json.dumps({'audit': pin(HERE / 'AUDIT.json'), 'cost': output['cost']}))
