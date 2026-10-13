"""Read-only audit/report for all 125 frozen Ca-conditioned sources; Slurm only."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import time
import traceback
from collections import Counter

POLICIES = ('operational', 'mathematical', 'origin')
ORDER = ['origin', 'adaptive_Ca', 'adaptive_La']
EV, HA = 23.06054783061903, 627.509474


# Same pinned-record and decision semantics as expanded_transfer_review/REVIEW.py.
def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def read(path):
    if isinstance(path, dict):
        assert pin(path['path']) == path, path
        path = path['path']
    return json.loads(Path(path).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def decision(value, bands, expected):
    if value is None:
        return {'decision': None, 'outcome': 'unavailable'}
    label = 'Ca' if value <= bands['Ca_max'] else 'La' if value >= bands['La_min'] else 'inconclusive'
    return {'decision': label, 'outcome': 'inconclusive' if label == 'inconclusive'
            else 'correct' if label == expected else 'wrong'}


def counts(rows):
    return {policy: dict(Counter(row['decisions'][policy]['outcome'] for row in rows))
            for policy in POLICIES}


def audit(plan):
    for item in plan['pins']:
        assert pin(item['path']) == item, item['path']
    spec = importlib.util.spec_from_file_location('qualified_raw_parser', plan['parser_helper'])
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)  # Pure parser helpers; no executable/backend calls.
    collection = read(plan['collection'])
    manifest = read(collection['manifest'])
    scope = read(manifest['scope'])
    native = read(manifest['MACE_collection'])
    canonical = read(collection['reference']['source'])
    previous = read(plan['previous_review'])
    calibration = collection['reference']['frozen_calibration']
    assert canonical['calibration'] == calibration == scope['reference']['frozen_calibration']
    assert manifest['reference'] == collection['reference']
    assert collection['reference']['new_fit'] is False and collection['new_calibration'] is False
    assert collection['new_Ca_triples'] is False and collection['production_changed'] is False
    assert manifest['pool_settings']['candidate_order'] == ORDER
    assert manifest['pool_settings']['origin_selection_tolerance_kcal_mol'] == .1
    assert manifest['settings']['component_parity_Eh'] == 1e-6
    assert manifest['settings']['population_parity_e'] == 1e-5
    assert native['terminal'] and len(scope['rows']) == len(collection['rows']) == 125
    declared = {r['case_id']: r for r in scope['rows']}
    cases = {c['case_id']: c for c in manifest['cases']}
    assert len(declared) == len(cases) == 125 and set(declared) == set(cases)
    assert {r['case_id'] for r in collection['rows']} == set(cases)
    raw_cache, absolute_errors, population_errors, component_errors = {}, [], [], []

    def raw(record):
        if record['sha256'] in raw_cache:
            return raw_cache[record['sha256']]
        data = read(record)
        task = data['task']
        helper.verify_raw(data, task | {'reuse': record}, manifest)
        energies, charges = {}, {}
        for backend in ('xtb', 'tblite'):
            stored = data['backends'][backend]
            receipt = read(stored['execution'])
            if stored['status'] != 'available':
                raw_cache[record['sha256']] = None
                return None
            assert receipt['returncode'] == 0 and receipt['owned_process_group_empty']
            folder = Path(stored['execution']['path']).parent
            if backend == 'xtb':
                text = Path(receipt['stdout']['path']).read_text()
                assert re.search(r'convergence criteria satisfied|SCC converged', text, re.I)
                assert 'FAILED TO CONVERGE' not in text.upper()
                energy = float(re.findall(r'TOTAL ENERGY\s+(-?\d+\.\d+)', text)[-1])
                charge = [float(v) for v in (folder / 'charges').read_text().split()]
            else:
                science = read(folder / 'scientific.json')
                assert science['library_returned_without_error']
                assert abs(sum(science['orbital_occupations_e']) - task['valence_electrons']) <= 1e-5
                energy, charge = science['energy_hartree'], science['charges_e']
            assert energy == stored['energy_hartree'] and charge == stored['charges_e']
            assert math.isfinite(energy) and len(charge) == task['atom_count']
            assert all(math.isfinite(v) for v in charge)
            assert abs(sum(charge) - task['charge']) <= 1e-5
            energies[backend], charges[backend] = energy, charge
        population = max(abs(a-b) for a, b in zip(charges['xtb'], charges['tblite']))
        assert abs(population - data['parity']['maximum_population_difference_e']) < 1e-12
        absolute_errors.append(energies['xtb'] - energies['tblite'])
        population_errors.append(population)
        result = {'energies': energies, 'population_error': population, 'task': task}
        raw_cache[record['sha256']] = result
        return result

    rows = []
    old_cases = {}
    for original in collection['rows']:
        case = cases[original['case_id']]
        row = {k: original[k] for k in ('case_id', 'source_case_id', 'selection_id', 'protein',
               'biological_group', 'expected_class', 'source_conditioning_metal', 'execution_status',
               'preparation_status', 'preparation_reason')}
        assert row['source_conditioning_metal'] == 'Ca'
        for key in ('source_case_id', 'selection_id', 'biological_group', 'expected_class'):
            assert row[key] == declared[row['case_id']][key]
        reasons = []
        if case['execution_status'] == 'reuse_actual_corrected_pool':
            prior = read(case['prior_collection'])
            prior_row = next(r for r in prior['rows'] if r['case_id'] == row['case_id'])
            assert prior_row == case['prior_row']
            key = prior['manifest']['sha256']
            if key not in old_cases:
                old_cases[key] = {c['case_id']: c for c in read(prior['manifest'])['cases']}
            scientific_case = old_cases[key][row['case_id']]
        else:
            scientific_case = case
        matrix = {z: {} for z in ('Ca', 'La')}
        if case['execution_status'] == 'excluded':
            assert not original['components']
            reasons.append({'kind': 'preparation', 'reason': row['preparation_reason']})
        else:
            signature = read(scientific_case['state_signature'])
            for metal in ('Ca', 'La'):
                for candidate in ORDER:
                    cell = scientific_case['matrix'][metal][candidate]
                    component = original['components'][metal + '__' + candidate]
                    matrix[metal][candidate] = None
                    if cell is None:
                        assert not component['accepted']
                        reasons.append({'kind': 'MACE_cell', 'metal': metal, 'candidate': candidate,
                                        'reason': component.get('reason')})
                        continue
                    mace = read(cell['MACE'])
                    assert mace['energy_eV'] == cell['MACE_eV'] and math.isfinite(mace['energy_eV'])
                    state = mace.get('state', mace.get('input_state_check'))
                    assert state['head'] == 'omol' and state['spin_multiplicity'] == cell['multiplicity'] == 1
                    assert state['charge'] == cell['charge'] == signature['endpoint_charges'][metal]
                    assert pin(cell['xyz']['path']) == cell['xyz']
                    values = {q: raw(component['inputs'][q]) if component['inputs'].get(q) else None
                              for q in ('vacuum', 'water')}
                    if any(v is None for v in values.values()):
                        assert not component['accepted']
                        reasons.append({'kind': 'solvent', 'metal': metal, 'candidate': candidate,
                                        'reason': component.get('reason')})
                        continue
                    vac, wet = values['vacuum'], values['water']
                    assert all(v['task']['xyz']['sha256'] == cell['xyz']['sha256'] and v['task']['charge'] == cell['charge']
                               and v['task']['metal'] == metal for v in values.values())
                    solvents = {b: wet['energies'][b] - vac['energies'][b] for b in ('xtb', 'tblite')}
                    error = solvents['xtb'] - solvents['tblite']
                    accepted = abs(error) <= 1e-6 and max(v['population_error'] for v in values.values()) <= 1e-5
                    assert accepted == component['accepted'] and error == component['component_error_hartree']
                    component_errors.append(error)
                    assert vac['energies']['xtb'] == component['primary_vacuum_energy_hartree']
                    assert wet['energies']['xtb'] == component['primary_water_energy_hartree']
                    if accepted:
                        matrix[metal][candidate] = (mace['energy_eV'], solvents['xtb'])
                    else:
                        reasons.append({'kind': 'solvent', 'metal': metal, 'candidate': candidate,
                                        'reason': component.get('reason')})
        if case.get('native_pool_status') == 'unavailable':
            source = case['native_source_record']
            reasons.append({'kind': 'native_pool_unavailable', 'status': source['status'],
                            'failed_candidates': {k: v for k, v in source['candidates'].items()
                                                  if v['status'] != 'available'}})
        complete = all(matrix[z].get(q) is not None for z in ('Ca', 'La') for q in ORDER)
        origin_available = all(matrix[z].get('origin') is not None for z in ('Ca', 'La'))
        row.update(scores_kcal_mol={}, decisions={}, selected_candidates={}, accommodation_works={}, reasons=reasons)
        works = {z: {q: (v[0]-matrix[z]['origin'][0])*EV + (v[1]-matrix[z]['origin'][1])*HA
                     for q, v in matrix[z].items() if v is not None}
                 for z in ('Ca', 'La')} if origin_available else {}
        for policy in POLICIES:
            available = origin_available if policy == 'origin' else complete and case.get('native_pool_status') != 'unavailable'
            score, chosen = None, None
            if available:
                best = {z: min(ORDER, key=lambda q: works[z][q]) if complete else 'origin' for z in ('Ca', 'La')}
                chosen = {z: 'origin' if policy == 'origin' else best[z] if policy == 'mathematical'
                          or works[z][best[z]] < -.1 else 'origin' for z in ('Ca', 'La')}
                ca, la = (matrix[z][chosen[z]] for z in ('Ca', 'La'))
                score = (ca[0]-la[0])*EV + (ca[1]-la[1])*HA
                assert original['pool'][policy] is not None
                assert abs(score-original['pool'][policy]['R_model_kcal_mol']) < 1e-7
                for z in ('Ca', 'La'):
                    assert chosen[z] == original['pool'][policy][z]['candidate']
                    assert abs(works[z][chosen[z]]-original['pool'][policy][z]['work_kcal_mol']) < 1e-7
            else:
                assert original['pool'][policy] is None
            row['scores_kcal_mol'][policy] = score
            row['selected_candidates'][policy] = chosen
            row['decisions'][policy] = decision(score, calibration[policy]['bands'], row['expected_class'])
            assert row['decisions'][policy] == original['decisions'][policy]
        row['accommodation_works'] = works
        row['delta_R_kcal_mol'] = row['scores_kcal_mol']['operational']-row['scores_kcal_mol']['origin'] if row['scores_kcal_mol']['operational'] is not None else None
        rows.append(row)
    assert counts(rows) == collection['counts']
    common = [r for r in rows if all(r['scores_kcal_mol'][p] is not None for p in ('operational', 'origin'))]
    assert len(common) == collection['common_origin_operational_coverage']['sources']
    transitions = dict(Counter(r['decisions']['origin']['outcome']+' -> '+r['decisions']['operational']['outcome'] for r in common))
    # Fresh backend costs are summed from actual immutable receipts, without reused calls.
    measured = {'calls': 0, 'CPU_seconds': 0., 'summed_call_wall_seconds': 0., 'maximum_process_RSS_KiB': 0}
    receipt_paths = set()
    for p in collection['cost']['execution_receipts']:
        assert p['path'] not in receipt_paths
        receipt_paths.add(p['path']); receipt = read(p)
        measured['calls'] += 1
        measured['CPU_seconds'] += receipt['cpu_user_seconds'] + receipt['cpu_system_seconds']
        measured['summed_call_wall_seconds'] += receipt['wall_seconds']
        measured['maximum_process_RSS_KiB'] = max(measured['maximum_process_RSS_KiB'], receipt['peak_RSS_KiB'])
    for key in ('CPU_seconds', 'summed_call_wall_seconds', 'maximum_process_RSS_KiB'):
        assert abs(measured[key]-collection['cost'][key]) < 1e-7
    assert measured['calls'] == collection['cost']['new_primary_calls'] + collection['cost']['new_independent_calls']
    stage = read(plan['stage_receipt']); assert stage['collection'] == pin(plan['collection'])
    exceptions = [r for r in rows if r['decisions']['operational']['outcome'] != 'correct'
                  or r['decisions']['origin'] != r['decisions']['operational']]
    old_la = [r for r in previous['rows'] if r['source_conditioning_metal'] == 'La'
              and r['role'] != 'canonical_calibration']
    for policy in POLICIES:
        assert previous['calibration'][policy]['bands'] == calibration[policy]['bands']
    result = {'status': 'audited', 'collection': pin(plan['collection']), 'manifest': collection['manifest'],
              'reference': collection['reference'], 'new_fit': False, 'new_Ca_triples': False,
              'declared_sources': 125, 'biological_groups': len({r['biological_group'] for r in rows}),
              'counts': counts(rows), 'common_coverage': {'sources': len(common), 'counts': counts(common),
              'transitions': transitions}, 'by_expected_class': {z: counts([r for r in rows if r['expected_class'] == z]) for z in ('Ca','La')},
              'execution_statuses': dict(Counter(r['execution_status'] for r in rows)),
              'exceptions': exceptions, 'rows': rows,
              'raw_audit': {'unique_inputs': len(raw_cache), 'components': len(component_errors),
              'maximum_component_error_Eh': max(map(abs, component_errors)),
              'maximum_population_error_e': max(population_errors),
              'absolute_v1_gate_failures_retained': sum(abs(v)>1e-6 for v in absolute_errors),
              'absolute_v1_maximum_error_Eh': max(map(abs, absolute_errors))},
              'cost': {'fresh_solvent': measured,
              'new_primary_calls': collection['cost']['new_primary_calls'],
              'new_independent_calls': collection['cost']['new_independent_calls'],
              'solvent_cost_coverage': collection['cost']['cost_coverage'],
              'solvent_stage_elapsed_seconds': stage['stage_elapsed_seconds'],
              'solvent_allocated_CPU_seconds_through_collection': stage['allocated_CPU_seconds_through_collection'],
              'native_MACE_receipted_cost': native['costs'],
              'native_MACE_cost_source': manifest['MACE_collection'],
              'scheduler_accounting': pin(plan['accounting']),
              'note': 'MACE costs retain unknown interrupted forwards; fresh solvent cost excludes reused pools. No GPU speedup extrapolation or complete end-to-end latency claim.'},
              'previous_evidence_kept_separate': {'review': pin(plan['previous_review']),
              'canonical25': previous['counts']['canonical25'], 'all_previous_strata': previous['counts'],
              'noncalibration_La_conditioned_rows': len(old_la), 'noncalibration_La_conditioned_counts': counts(old_la),
              'existing_triples100': previous['counts']['triples100'],
              'note': 'Consumed development references; canonical calibration is training fidelity, and repeated folds/triples are correlated. No populations combined or new triples constructed.'},
              'scientific_calls_by_review': 0, 'production_changed': False}
    return result


def report(result, output):
    write(output/'RESULT.json', result)
    fields = ['case_id','protein','biological_group','expected_class','source_conditioning_metal',
              'preparation_status','execution_status','origin_score','operational_score','mathematical_score',
              'origin_outcome','operational_outcome','mathematical_outcome','delta_R_kcal_mol','reasons']
    for name, population in [('ALL125.tsv', result['rows']), ('EXCEPTIONS.tsv', result['exceptions'])]:
        with (output/name).open('x') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
            for row in population:
                flat = {k: row.get(k) for k in fields}
                flat.update({p+'_score': row['scores_kcal_mol'][p] for p in POLICIES})
                flat.update({p+'_outcome': row['decisions'][p]['outcome'] for p in POLICIES})
                flat['reasons'] = json.dumps(row['reasons'], separators=(',', ':'))
                writer.writerow(flat)
    lines = ['# Corrected Ca-conditioned transfer — final review', '',
             'All 125 declared sources are retained. Frozen canonical bands are unchanged; no new calibration or triples were created.', '',
             '| Population / policy | Correct | Wrong | Inconclusive | Unavailable |',
             '|---|---:|---:|---:|---:|']
    for label, data in [('Full125', result['counts']), ('Common'+str(result['common_coverage']['sources']), result['common_coverage']['counts'])]:
        for policy in ('origin','operational','mathematical'):
            c = data[policy]; lines.append('| '+label+' / '+policy+' | '+' | '.join(str(c.get(k,0)) for k in ('correct','wrong','inconclusive','unavailable'))+' |')
    lines += ['', 'Matched transitions: '+json.dumps(result['common_coverage']['transitions'],sort_keys=True)+'.', '',
              '## Remaining errors and abstentions', '']
    for row in result['rows']:
        if row['decisions']['operational']['outcome'] in ('wrong','inconclusive'):
            lines.append('- `'+row['source_case_id']+'`: '+row['decisions']['origin']['outcome']+' → '+row['decisions']['operational']['outcome']+'.')
    lines += ['', 'All exclusions, interrupted native searches and changed decisions are listed in EXCEPTIONS.tsv; ALL125.tsv retains the full denominator. A missing accommodated pool is never filled with its valid origin score.', '',
              '## Prior evidence remains separate', '',
              'The existing 25-member canonical panel supplied these bands. Its fidelity is calibration evidence, not new predictive accuracy. Existing La-conditioned development rows and the established100 correlated triples are carried separately in RESULT.json; no new aggregate rule is invented.', '',
              '## Components and measured cost', '',
              'The audit reparsed actual backend outputs and checked unchanged component-v2 solvent differences, metal states, row minima, origin retention tolerance and all decisions. Historical absolute-v1 parity failures remain failed; this review does not reinterpret them.', '',
              json.dumps(result['raw_audit'],sort_keys=True), '',
              'Fresh solvent: '+str(result['cost']['new_primary_calls'])+' xTB calls and '+str(result['cost']['new_independent_calls'])+' tblite checks; '+format(result['cost']['fresh_solvent']['CPU_seconds'],'.3f')+' measured process CPU-seconds. Stage wall time through collection: '+format(result['cost']['solvent_stage_elapsed_seconds'],'.3f')+' s; '+format(result['cost']['solvent_allocated_CPU_seconds_through_collection'],'.3f')+' allocated CPU-seconds. These are different measures.', '',
              'MACE recorded '+str(result['cost']['native_MACE_receipted_cost']['recorded_model_forwards'])+' forwards and '+format(result['cost']['native_MACE_receipted_cost']['recorded_inference_process_cpu_seconds'],'.3f')+' inference CPU-seconds. Interrupted forwards without terminal receipts retain unknown cost. No CPU-to-GPU extrapolation, end-to-end speedup, production promotion, or independent biological validation is claimed.', '',
              'Recommendation: retain this corrected candidate as development evidence while keeping the remaining wrong/inconclusive cases and source dependence visible. Recalibration alone does not establish broader metal specificity.', '']
    with (output/'REPORT.md').open('x') as handle: handle.write('\n'.join(lines))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'): raise SystemExit('Scientific review requires Slurm')
    plan = read(args.plan); output = Path(plan['output_directory']); started = time.time()
    receipt = {'job_id': os.environ['SLURM_JOB_ID'], 'plan': pin(args.plan), 'started_epoch': started, 'new_molecular_calls': 0}
    write(output/'START.json', receipt)
    try:
        assert int(os.environ['SLURM_CPUS_ON_NODE']) == 1
        result = audit(plan); report(result, output)
        receipt.update(status='complete', result=pin(output/'RESULT.json'), report=pin(output/'REPORT.md'))
    except Exception as error:
        receipt.update(status='failed', error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        receipt['finished_epoch'] = time.time(); receipt['wall_seconds'] = receipt['finished_epoch']-started
        write(output/'RECEIPT.json', receipt)
