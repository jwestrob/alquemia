#!/usr/bin/env python3
"""Independent raw-token arithmetic; only run after the full saved collection."""
import argparse
from collections import Counter
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import time

HERE = Path(__file__).resolve().parent
BRANCH = HERE.parent
ROOT = BRANCH.parents[3]
WORK = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/common_basis_decomposition_v1'
GIDS = ('aquo8_central_La', 'aquo8_central_Dy', 'aquo9_central_La', 'aquo9_central_Dy')
SIGN = dict(zip(GIDS, (1, -1, -1, 1)))
QUANTITIES = ('cluster_contrast', 'pair_contrast', 'cooperative_remainder')
HA = Decimal('627.509474')
NUMBER = r'[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?'


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def readpin(p):
    assert pin(p['path']) == p, 'Changed pin: ' + p['path']
    return json.loads(Path(p['path']).read_text())


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def decimal_json(value):
    if isinstance(value, Decimal):
        return str(value)
    raise TypeError(type(value).__name__)


def main(result_path, output):
    tick, cpu = time.monotonic(), time.process_time()
    assert all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'))
    result = json.loads(result_path.read_text())
    assert result['status'] == 'common_basis_collection_complete' and result['new_complete'] == 68 and result['reused_complete'] == 12 and result['error'] is None, 'Full terminal result required'
    assert not output.exists(), 'Preserve existing independent review'
    checks = {}

    def check(name, condition):
        assert condition, name
        checks[name] = True

    def close(name, actual, expected, tolerance='1e-9'):
        check(name, math.isfinite(float(expected)) and abs(Decimal(actual) - Decimal(str(expected))) <= Decimal(tolerance))

    source_review = json.loads((BRANCH / 'independent_review_v1/COLLECTION_REVIEW.json').read_text())
    check('collector_source_pin', result['implementation'] == pin(BRANCH / 'collection/COLLECT.py') == source_review['implementation'])
    check('panel_pin', result['panel'] == source_review['panel'] == pin(WORK / 'prepared_v1/PANEL.json'))
    cfg = readpin(result['configuration'])
    check('collector_config_lineage', cfg['implementation'] == result['implementation'] and cfg['panel'] == result['panel'])
    a = load_module(BRANCH / 'REFERENCE.py', 'independent_common_reference')
    collector = load_module(BRANCH / 'collection/COLLECT.py', 'independent_common_admissions')
    p = a.load(Path(result['panel']['path']))
    check('adapter_pin', p['adapter'] == pin(BRANCH / 'REFERENCE.py') == cfg['adapter'])
    expected_ids = {c['cell_id'] for c in p['cells']}
    check('full68_slots', len(expected_ids) == 68 and set(result['new_fragments']) == set(cfg['result_paths']) == expected_ids)
    check('full12_reuse_slots', len(result['reused_terms']) == 12 and set(result['reused_terms']) == set(p['saved_reuse']))
    old = a.helper.load(a.b.verify(p['cp_panel']))
    cp_result, pair_result, data = readpin(p['cp_collection']), readpin(p['pair_collection']), readpin(p['consumed_input_data'])
    pair_review_path = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/matched_pair_decomposition_v1/independent_review_v1/RESULT_REVIEW.json'
    check('inherited_pair_actual_review_pin', pin(pair_review_path)['sha256'] == 'b008e6718f303f6756092622e979b14b664f158eed9a65960963cdfeb43f561f')
    check('prior_complete', cp_result['status'] == 'coordination_cp_collection_complete' and pair_result['status'] == 'matched_pair_collection_complete')

    admissions, tokens = {}, {}

    def take(cid, admitted, stored, kind):
        check(cid + ':complete', admitted['status'] == stored['status'] == 'complete')
        check(cid + ':authority', admitted['result'] == stored['result'] and admitted['output'] == stored['output'])
        raw = readpin(admitted['result'])
        text = a.b.verify(admitted['output']).read_text()
        values = re.findall(r'^Total Energy\s+:\s+(' + NUMBER + r')\s+Eh', text, re.M)
        check(cid + ':one_electronic_token', len(values) == 1)
        electronic = Decimal(values[0])
        close(cid + ':raw_electronic_scalar', electronic, stored['electronic_energy_hartree'], '1e-12')
        close(cid + ':admitted_electronic_scalar', electronic, admitted['electronic_energy_hartree'], '1e-12')
        d4_tokens = re.findall(r'^\s*Dispersion correction\s+(' + NUMBER + r')\s*$', text, re.M)
        if cid in GIDS:
            check(cid + ':one_original_D4_token', len(d4_tokens) == 1)
            d4 = Decimal(d4_tokens[0])
        else:
            check(cid + ':no_ghost_D4_token', not d4_tokens)
            d4 = Decimal(0)
        close(cid + ':D4_scalar', d4, stored['D4_energy_hartree'], '1e-12')
        check(cid + ':normal_native_termination', 'ORCA TERMINATED NORMALLY' in text and raw['status'] == 'complete')
        check(cid + ':all_native_checks', bool(stored['checks']) and all(stored['checks'].values()))
        # The original EnGrad parser stored [] even though the raw final table exists.
        # Extract the exact final table for every row, preserving that historical omission.
        tables = text.split('SCF CONVERGENCE\n')
        check(cid + ':one_final_raw_SCF_table', len(tables) == 2)
        table = tables[1].split('ORBITAL ENERGIES', 1)[0]
        table_rows = re.findall(r'^\s*Last (Energy change|MAX-Density change|RMS-Density change|DIIS Error|Orbital Gradient|Orbital Rotation)\s+\.\.\.\s+(' + NUMBER + r')\s+Tolerance\s*:\s*(' + NUMBER + r')\s*$', table, re.M)
        criteria = [{'name': name, 'value': float(value), 'tolerance': float(tol), 'passes': abs(float(value)) <= float(tol)} for name, value, tol in table_rows]
        check(cid + ':finite_scf_metadata', len(criteria) == 6 and len({v['name'] for v in criteria}) == 6 and all(math.isfinite(v['value']) and math.isfinite(v['tolerance']) for v in criteria))
        historical = raw['printed_scf_criteria']
        check(cid + ':stored_metadata_preserved', historical == stored['printed_scf_criteria'])
        if cid in GIDS:
            check(cid + ':known_original_metadata_omission', historical == [])
        else:
            check(cid + ':raw_SCF_matches_native_metadata', criteria == historical)
        stop = raw.get('native_scf_stop_signals')
        if not stop:
            stop = re.findall(r'\*{2,}\s*(.*?)\s+signals convergence\s*\*{2,}', text)
        check(cid + ':stopping_signal_recorded', bool(stop))
        admissions[cid] = {
            'kind': kind, 'result': admitted['result'], 'output': admitted['output'],
            'input': admitted.get('input', admitted.get('source_input')), 'job_id': admitted.get('job_id'),
            'electronic_token_hartree': values[0], 'D4_token_hartree': str(d4),
            'printed_scf_criteria': criteria, 'historical_stored_scf_criteria': historical, 'raw_SCF_table_extracted_independently': True, 'native_scf_stop_signals': stop,
            'native_warning_lines': [line.strip() for line in text.splitlines() if 'WARNING:' in line],
            'native_mode': raw.get('native_scf_mode'),
            'energy_stage_metadata': raw.get('energy_stage_metadata', stored.get('energy_stages')),
            'final_COSX_correction_hartree': raw.get('final_COSX_correction_hartree', stored.get('energy_stages', {}).get('final_COSX_correction_hartree')),
            'native_elapsed_times_ranks_CPU_seconds': admitted.get('allocated_CPU_seconds'),
        }
        tokens[cid] = electronic
        return d4

    original_d4 = {}
    for cid, stored in result['reused_terms'].items():
        admitted = collector.admit_reused(cid, p, old, cp_result, data)
        original_d4[cid] = take(cid, admitted, stored, 'reused')
    for cell in p['cells']:
        cid = cell['cell_id']
        stored = result['new_fragments'][cid]
        check(cid + ':configured_result_path', Path(cfg['result_paths'][cid]).resolve() == Path(stored['result']['path']).resolve())
        admitted = a.admit_result(cfg['result_paths'][cid], cell, p, Path(result['panel']['path']))
        take(cid, admitted, stored, 'new')

    def gid_cell(gid, role, index=None):
        n = p['geometries'][gid]['water_count']
        metal = gid.rsplit('_', 1)[1]
        return f'pbe0__aquo{n}__{metal}__' + (f'water_{index:02d}__' if index is not None else '') + role

    with localcontext() as context:
        context.prec = 40
        parents = {}
        for gid in GIDS:
            n = p['geometries'][gid]['water_count']
            E, ion, shell = tokens[gid], tokens[gid_cell(gid, 'ion_ghost_shell')], tokens[gid_cell(gid, 'shell_ghost_ion')]
            d = [tokens[gid_cell(gid, 'metal_water', i)] for i in range(1, n + 1)]
            w = [tokens[gid_cell(gid, 'water_only', i)] for i in range(1, n + 1)]
            interactions = [di - ion - wi for di, wi in zip(d, w)]
            C, P = E - ion - shell, sum(d) - n * ion - sum(w)
            B = E - sum(d) + (n - 1) * ion - shell + sum(w)
            check(gid + ':exact_decimal_closure', C - P - B == 0 and sum(interactions) == P)
            saved = result['decomposition']['per_parent_electronic'][gid]
            for key, value in [('cluster_interaction_hartree', C), ('pair_sum_hartree', P), ('cooperative_remainder_hartree', B), ('expanded_remainder_hartree', B), ('N_fold_ion_subtraction_hartree', n * ion)]:
                close(gid + ':' + key, value, saved[key])
            check(gid + ':all_parent_water_slots', len(saved['pair_interactions_hartree']) == n)
            for i, value in enumerate(interactions):
                close(gid + ':pair_' + str(i + 1), value, saved['pair_interactions_hartree'][i])
            parents[gid] = {'E': E, 'ion': ion, 'shell': shell, 'd': d, 'w': w, 'pair_interactions': interactions, 'C': C, 'P': P, 'B': B, 'N': n}
        electronic = {label: sum(Decimal(SIGN[g]) * parents[g][q] for g in GIDS) for label, q in zip(QUANTITIES, ('C', 'P', 'B'))}
        check('paired_exact_decimal_closure', electronic['cluster_contrast'] - electronic['pair_contrast'] - electronic['cooperative_remainder'] == 0)
        close('prior_CP8_C_replay', electronic['cluster_contrast'], cp_result['cycle']['hartree']['CP_electronic'])
        d4_C = sum(Decimal(SIGN[g]) * original_d4[g] for g in GIDS)
        old_d4 = pair_result['reference_decomposition']['D4_energy_hartree']['hartree']
        d4_P = Decimal(str(old_d4['pair_contrast']))
        d4 = dict(zip(QUANTITIES, (d4_C, d4_P, d4_C - d4_P)))
        close('D4_C_prior_CP8_replay', d4_C, cp_result['cycle']['hartree']['original_D4_unchanged'], '1e-12')
        for q in QUANTITIES:
            close('D4_prior_pair_' + q, d4[q], old_d4[q], '1e-12')
        combined = {q: electronic[q] + d4[q] for q in QUANTITIES}
        model = {q: Decimal(str(pair_result['model_decomposition'][q])) for q in QUANTITIES}
        old_total = {q: Decimal(str(pair_result['reference_decomposition']['final_energy_hartree']['kcal_mol'][q])) for q in QUANTITIES}
        errors = {q: model[q] - combined[q] * HA for q in QUANTITIES}
        close('frozen_model_closure', model['cluster_contrast'] - model['pair_contrast'] - model['cooperative_remainder'], 0, '1e-9')
        close('model_error_closure', errors['cluster_contrast'] - errors['pair_contrast'] - errors['cooperative_remainder'], 0, '1e-9')
        for q in QUANTITIES:
            for key, value, tol in [
                ('electronic_hartree', electronic[q], '1e-9'), ('saved_original_D4_hartree', d4[q], '1e-12'),
                ('electronic_plus_original_D4_hartree', combined[q], '1e-9'),
                ('electronic_kcal_mol', electronic[q] * HA, '1e-6'), ('saved_original_D4_kcal_mol', d4[q] * HA, '1e-6'),
                ('electronic_plus_original_D4_kcal_mol', combined[q] * HA, '1e-6'),
                ('unchanged_model_kcal_mol', model[q], '1e-9'), ('old_own_basis_reference_kcal_mol', old_total[q], '1e-9'),
                ('model_minus_common_basis_reference_kcal_mol', errors[q], '1e-6'),
                ('model_minus_old_own_basis_reference_kcal_mol', model[q] - old_total[q], '1e-6'),
                ('common_minus_old_own_basis_reference_kcal_mol', combined[q] * HA - old_total[q], '1e-6'),
            ]:
                close(key + ':' + q, value, result['decomposition'][key][q], tol)

    statistics = {}
    for kind in ('new', 'reused'):
        selected = [r for r in admissions.values() if r['kind'] == kind]
        statistics[kind] = {
            'cells': len(selected),
            'stopping_signals': dict(Counter(s for r in selected for s in r['native_scf_stop_signals'])),
            'failed_printed_criteria': dict(Counter(v['name'] for r in selected for v in r['printed_scf_criteria'] if not v['passes'])),
            'all_printed_criteria_passed_cells': sum(all(v['passes'] for v in r['printed_scf_criteria']) for r in selected),
            'warning_lines': dict(Counter(s for r in selected for s in r['native_warning_lines'])),
        }
    new_costs = [r['native_elapsed_times_ranks_CPU_seconds'] for r in admissions.values() if r['kind'] == 'new']
    check('all68_native_costs_present', len(new_costs) == 68 and all(isinstance(v, (int, float)) and math.isfinite(v) and v > 0 for v in new_costs))
    close('new_cost_sum_replay', Decimal(str(math.fsum(new_costs))), result['cost']['recorded_new_allocated_CPU_seconds'], '1e-6')
    check('no_native_calls_or_target_change', result['new_native_calls'] == 0 and result['fit_targets_changed'] is False and result['physical_accuracy_qualified'] is False)
    review = {
        'status': 'passed_independent_actual_common_basis_review', 'implementation': pin(__file__), 'declaration': pin(HERE / 'DECLARATION.md'),
        'collector_result': pin(result_path), 'collector_configuration': result['configuration'], 'panel': result['panel'],
        'preserved_failed_review': pin(WORK / 'independent_result_review_v1/ATTEMPT1_FAILURE.json'),
        'reviewer_v2_change': 'Independent extraction of actual six-row SCF table; preserve original four empty historical metadata arrays. No scientific or tolerance change.',
        'native_source': p['adapter'], 'collector_source': result['implementation'], 'prior_CP_result': p['cp_collection'],
        'prior_pair_result': p['pair_collection'], 'inherited_model_and_D4_raw_review': pin(pair_review_path),
        'passed_checks': len(checks), 'checks': checks, 'admissions': admissions, 'complete_new_terms': 68, 'complete_reused_terms': 12,
        'decimal_per_parent_hartree': parents, 'decimal_paired_electronic_hartree': electronic, 'decimal_saved_D4_hartree': d4,
        'decimal_paired_total_kcal_mol': {q: combined[q] * HA for q in QUANTITIES},
        'decimal_unchanged_model_kcal_mol': model, 'decimal_model_minus_common_basis_kcal_mol': errors,
        'decimal_common_minus_old_reference_kcal_mol': {q: combined[q] * HA - old_total[q] for q in QUANTITIES},
        'scf_statistics': statistics,
        'limitations': ['Electronic common-parent-basis conditional decomposition; B is not a unique physical polarization or charge-transfer component.',
                        'Original D4 is kept separate and added once, with no ghost D4 calls.',
                        'Retained native SCF stopping warnings and failed printed criteria preclude universal convergence/accuracy claims.',
                        'No forces, new fit, target update, model promotion or validation release.'],
        'cost': {'review_wall_seconds': time.monotonic() - tick, 'review_CPU_seconds': time.process_time() - cpu, 'local_threads': 1, 'new_native_calls': 0,
                 'new_native_elapsed_times_ranks_CPU_seconds': math.fsum(new_costs),
                 'reused_native_cost_excluded_from_new': sum(r['native_elapsed_times_ranks_CPU_seconds'] or 0 for r in admissions.values() if r['kind'] == 'reused'),
                 'collector_saved_cost': result['cost'], 'scheduler_accounting': 'Retained by root separately; native elapsed×ranks is not measured utilization.'},
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(review, default=decimal_json, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': review['status'], 'review': pin(output), 'checks': len(checks), 'errors_kcal_mol': errors}, default=decimal_json))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--result', type=Path, default=WORK / 'collection_v1/result/RESULT.json')
    parser.add_argument('--output', type=Path, default=WORK / 'independent_result_review_v1/REVIEW.json')
    args = parser.parse_args()
    main(args.result, args.output)
