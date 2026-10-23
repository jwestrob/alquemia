#!/usr/bin/env python3
"""Saved six-parameter review; exact evaluator replay and native producer accounting."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import time

GEOMETRIES = ('water_minus', 'water_central', 'water_plus', 'aquo8_central')
CELLS = {q+'_'+m for q in GEOMETRIES for m in ('La', 'Dy')}


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def token(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def brief(row):
    s = row['score']
    return {k: s[k] for k in ('objective', 'blocks', 'block_weights', 'energy_errors_kcal_mol',
        'energy_checks_passed', 'controls', 'development_checks_passed')}


def fixture_summary(score):
    return {kind: {q: {k: v for k, v in row.items() if k not in
        ('raw_body_components_kcal_mol_A', 'model_body_response', 'reference_body_response')}
        for q, row in cells.items()} for kind, cells in score['per_fixture'].items()}


def main(path, output, completion=None):
    started = time.monotonic(); cpu = time.process_time(); source = Path(path).resolve()
    if not source.is_file():
        print(json.dumps({'status': 'result_absent', 'no_report_written': True})); return 2
    result = json.loads(source.read_text())
    if result['status'] != 'bounded_six_parameter_water_complete':
        print(json.dumps({'status': 'not_completed', 'producer_status': result['status'],
            'error': result.get('error'), 'no_report_written': True})); return 2
    verified = {}
    def check(p):
        if p['path'] not in verified: verified[p['path']] = pin(p['path'])
        if verified[p['path']] != p: raise ValueError('Changed evidence: '+p['path'])
        return Path(p['path'])
    def read(p): return json.loads(check(p).read_text())
    def module(p, name):
        spec = importlib.util.spec_from_file_location(name, check(p))
        obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj); return obj
    cfg = read(result['config'])
    for p in cfg.values():
        if isinstance(p, dict) and set(p) == {'path', 'sha256'}: check(p)
    assert result['implementation'] == cfg['implementation']
    data = read(cfg['training_result']); ecfg = read(cfg['evaluator_config'])
    native = module(ecfg['native_helper'], 'saved_radial_native_parser')
    objective = module(cfg['objective'], 'saved_radial_objective')
    body = module(cfg['body_helper'], 'saved_radial_body')
    center = module(cfg['center_objective'], 'saved_center_objective')
    ledger = read(result['evaluator']['ledger'])
    assert ledger['native_program_invocations'] == result['new_native_program_invocations'] <= 4096
    for p in result['protected_parents']: check(p)
    prior = read(result['prior_result']); prior_cfg = read(prior['config'])
    assert prior['status'] == 'bounded_radial_continuation_complete' and prior['error'] is None
    assert result['initial']['parameters'] == prior['selected']['parameters']+[.24429636467571708,.37332330078584813]
    assert result['initial']['score'] == prior['selected']['score']
    ancestor = read(prior_cfg['prior_result']); ancestor_cfg = read(ancestor['config'])
    epsilon = read(cfg['epsilon_result']); assert epsilon['status'] == 'epsilon_diagnostic_complete'
    producers = {'17593_native':(prior,result['prior_result'],prior_cfg['evaluator']),
                 '17589_native':(ancestor,prior_cfg['prior_result'],ancestor_cfg['evaluator'])}
    memberships = {kind:{(gid,cell['provenance']['native_authority']['path'],cell['provenance']['native_authority']['sha256'])
        for vector in parent['vectors'].values() for gid,cell in vector['model'].items()} for kind,(parent,_,_) in producers.items()}
    epsilon_members = {(gid,str(Path(cfg['epsilon_result']['path']).parent/probe/gid/'RESULT.json')):row
        for probe,record in epsilon['records'].items() for gid,row in record['cells'].items() if not row.get('reused',False)}
    vectors = result['vectors']; new_native = {}; inherited_native = {}; producer_counts = {}
    for key, row in vectors.items():
        assert row['status'] == 'complete' and token(row['parameters']) == key and set(row['model']) == CELLS
        for gid, cell in row['model'].items():
            pv = cell['provenance']; authority = pv['native_authority']; check(authority)
            assert pv['request_cache_implementation'] == cfg['evaluator']
            for group in ('source_input','effective_input'):
                for p in pv[group].values():check(p)
            for name in ('analyze_output','gradient_output'):check(cell['native'][name])
            assert cell['energy'] == cell['native']['total_energy_kcal_mol'] and cell['forces'] == cell['native']['force_kcal_mol_A']
            raw = read(authority)
            if pv['reused_baseline']:
                assert authority == data['rows'][gid]['baseline_result'] and pv['native_implementation'] is None
                assert all(raw['rows'][gid]['model'][k] == cell[k] for k in ('energy','forces','native'))
                assert raw['rows'][gid]['input']['files'] == pv['effective_input'];kind='17578_baseline'
            else:
                assert raw['status']=='complete' and len(raw['receipts'])==2 and all(x['returncode']==0 for x in raw['receipts'])
                for receipt in raw['receipts']:check(receipt['output'])
                kind = pv['inherited_source_kind'] if pv['inherited_native'] else 'new_six_parameter_native'
                if kind == '17596_epsilon':
                    assert epsilon_members[(gid,authority['path'])] == raw and pv['native_implementation']==epsilon['implementation'] and pv['native_cache_key'] is None
                    assert all(raw['model'][k]==cell[k] for k in ('energy','forces','native'))
                    for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
                        assert pin(Path(authority['path']).parent/name)==pv['effective_input'][k]
                else:
                    assert raw['cell_id']==gid and all(raw[k]==cell[k] for k in ('energy','forces','native'))
                    assert raw['input']==pv['effective_input'] and raw['cache_key']==pv['native_cache_key']
                    assert raw['cache_identity']['implementation']==pv['native_implementation']
                    if kind in producers:
                        parent,parent_pin,implementation=producers[kind]
                        assert (gid,authority['path'],authority['sha256']) in memberships[kind]
                        assert Path(authority['path']).is_relative_to(Path(parent_pin['path']).parent/'cache/cells')
                        assert pv['native_implementation']==implementation
                    else:
                        assert kind=='new_six_parameter_native' and pv['native_implementation']==cfg['evaluator']
                        assert Path(authority['path']).is_relative_to(source.parent/'cache/cells')
                if pv['inherited_native']:inherited_native[authority['path']]=raw
                else:new_native[authority['path']]=raw
            producer_counts.setdefault(kind,set()).add(authority['path'])
    assert 2*len(new_native) == result['new_native_program_invocations']
    assert len(vectors)-1 == result['new_parameter_vectors'] <= 384
    assert result['completed_outer_steps'] == len(result['iterations']) <= 24
    assert result['original_results_unchanged'] and result['candidate_frozen']
    selected = result['selected']; selected_key = token(selected['parameters'])
    assert selected == vectors[selected_key]
    assert objective.evaluate(data, selected['model'], body) == selected['score']
    out = Path(output); out.mkdir(parents=True, exist_ok=False)
    evaluator = module(cfg['evaluator'],'saved_six_exact_evaluator')
    ev = evaluator.Evaluator(ecfg,out/'replay_cache',root_release=None)
    actual_calls=[]
    def forbidden(*args,**kwargs):actual_calls.append(True);raise AssertionError('Saved review forbids native calls')
    ev.native.invoke=forbidden
    (out/'replay_cache/cells').mkdir(exist_ok=True)
    for gid,cell in selected['model'].items():
        pv=cell['provenance']
        if not pv['reused_baseline'] and not pv['inherited_native']:
            link=out/'replay_cache/cells'/pv['cache_key']
            link.symlink_to(Path(pv['native_authority']['path']).parent,target_is_directory=True)
    replayed=ev.evaluate(selected['parameters'])
    assert replayed==selected['model'] and not actual_calls and ev.snapshot()['native_program_invocations']==0
    for p in result['protected_parents']:assert pin(p['path'])==p
    assert pin(result['evaluator']['ledger']['path'])==result['evaluator']['ledger']
    selected_native = {}
    for gid, cell in selected['model'].items():
        replay = native.parse(check(cell['native']['analyze_output']), check(cell['native']['gradient_output']),
            {'atoms': len(cell['forces']), 'derivative_test': False})
        assert replay == cell['native']
        selected_native[gid] = {'authority': cell['provenance']['native_authority'],
            'input': cell['provenance']['effective_input'], 'analyze_output': replay['analyze_output'],
            'gradient_output': replay['gradient_output'], 'energy_kcal_mol': cell['energy'], 'raw_replay_exact': True,'producer_kind':cell['provenance']['inherited_source_kind'] if cell['provenance']['inherited_native'] else '17578_baseline' if cell['provenance']['reused_baseline'] else 'new_six_parameter_native','native_implementation':cell['provenance']['native_implementation']}
    center_data = dict(data, rows={k: dict(data['rows'][k], baseline=data['rows'][k]['original_baseline']) for k in center.CELLS})
    center_score = center.evaluate(center_data, {k: selected['model'][k] for k in center.CELLS}, body)
    assert center_score == result['selected_center_diagnostics']['original_four_center_score']
    iterations = []; passing_unselected = []; current = token(result['initial']['parameters'])
    for entry in result['iterations']:
        assert entry['status'] == 'complete' and entry['origin_key'] == current
        assert len(entry['derivative_keys']) == len(entry['derivative_slots']) == 12
        assert len(entry['trial_keys']) == len(entry['trial_slots']) == 4
        assert all(k in vectors for k in entry['derivative_keys']+entry['trial_keys'])
        best = current
        for k in entry['trial_keys']:
            if vectors[k]['score']['objective'] < vectors[best]['score']['objective']: best = k
        assert best == entry['selected_key']
        trials = []
        for slot, k in zip(entry['trial_slots'], entry['trial_keys']):
            assert slot['parameters'] == vectors[k]['parameters']
            row = dict(id=slot['id'], selected=k == best, **brief(vectors[k]))
            trials.append(row)
            if row['development_checks_passed'] and k != best:
                passing_unselected.append(dict(outer_step=entry['outer_step'], **row))
        diagnostic = entry['linear_directions']; singular = diagnostic['singular_values']
        residual = vectors[current]['score']['weighted_residuals']; J = entry['Jacobian_normalized']
        delta = diagnostic['unrestricted_GN_not_evaluated']
        linear = [r+sum(j*d for j, d in zip(row, delta)) for r, row in zip(residual, J)]
        origin_loss = sum(x*x for x in residual); projected_loss = sum(x*x for x in linear)
        iterations.append({'outer_step': entry['outer_step'], 'origin_key': current, 'origin': brief(vectors[current]),
            'selected_key': best, 'selected': brief(vectors[best]), 'trials': trials,
            'selected_per_fixture': fixture_summary(vectors[best]['score']),
            'origin_Jacobian_rank': diagnostic['rank'], 'origin_singular_values': singular,
            'origin_condition': singular[0]/singular[-1] if singular[-1] else None,
            'origin_unrestricted_linear_residual_loss': projected_loss,
            'origin_local_projected_fraction': 1-projected_loss/origin_loss if origin_loss else 0.,
            'derivative_actual_widths': entry['derivative_actual_widths']})
        current = best
    assert current == selected_key
    assert result['development_checks_passed'] == selected['score']['development_checks_passed']
    if result['stopping_reason'] == 'all_existing_controls_passed': assert result['development_checks_passed']
    parameters = []
    for j, p in enumerate(result['family']['parameters']):
        x = selected['parameters'][j]; z = selected['normalized_parameters'][j]
        assert p['lower'] <= x <= p['upper']
        parameters.append({'name': p['name'], 'initial': result['initial']['parameters'][j], 'selected': x,
            'move': x-result['initial']['parameters'][j], 'lower': p['lower'], 'upper': p['upper'],
            'normalized_selected': z, 'distance_to_nearest_bound_normalized': min(z, 1-z),
            'exact_bound_hit': x in (p['lower'], p['upper'])})
    failures = []
    for kind, cells in selected['score']['per_fixture'].items():
        for q, row in cells.items():
            if not row['passed']:
                raw = row['raw_body_components_kcal_mol_A']; idx = max(range(len(raw)), key=lambda j: abs(raw[j]))
                failures.append({'kind': kind, 'fixture': q, **fixture_summary(selected['score'])[kind][q],
                    'largest_component_index': idx, 'largest_signed_component_kcal_mol_A': raw[idx],
                    'model_body_response': row['model_body_response'], 'reference_body_response': row['reference_body_response']})
    event = json.loads(Path(completion).read_text()) if completion else None
    summary = {'status': 'saved_result_review_passed', 'reporter': pin(__file__), 'result': pin(source),
        'config': result['config'], 'implementation': result['implementation'], 'prior_result': result['prior_result'], 'stopping_reason': result['stopping_reason'],
        'initial': brief(result['initial']), 'selected': brief(selected), 'selected_per_fixture': fixture_summary(selected['score']),
        'remaining_failed_fixtures': failures, 'parameters': parameters, 'selected_parameters': selected['parameters'],
        'original_pair_epsilon': result['family']['original_pair_epsilon'], 'iterations': iterations,
        'passing_unselected_trials': passing_unselected, 'original_four_center_score': center_score,
        'method_sensitivity': result['selected_center_diagnostics']['method_sensitivity'],
        'selected_native': selected_native, 'all_selected_raw_replays_exact': True,
        'verified_unique_files': len(verified), 'new_native_cells': len(new_native), 'native_authority_counts_by_producer':{k:len(v) for k,v in producer_counts.items()}, 'exact_evaluator_replay':True, 'replay_evaluator':ev.snapshot(),
        'inherited_native_cells_used': len(inherited_native), 'inherited_native_cells_indexed': result['evaluator']['inherited_native_cells_indexed'],
        'new_native_program_invocations': result['new_native_program_invocations'],
        'new_parameter_vectors': result['new_parameter_vectors'], 'complete_outer_steps': len(iterations),
        'cost': {k: result[k] for k in ('elapsed_seconds', 'parent_process_CPU_seconds',
            'aggregate_native_child_user_CPU_seconds', 'aggregate_native_child_system_CPU_seconds')},
        'scheduler_records': event.get('records') if event else None,
        'last_Jacobian_at_final_selected_point': bool(iterations and iterations[-1]['origin_key'] == selected_key),
        'global_optimum_claimed': False, 'parameter_family_impossibility_claimed': False,
        'new_model_evaluations': 0, 'physical_accuracy_qualified': False, 'validation_outcomes_opened': False}
    lines = ['# Six-parameter paired water development: actual result', '',
        f"Stop: `{result['stopping_reason']}` after {len(iterations)} complete outer steps from the exact17593 selection and original epsilons. All development controls pass: **{result['development_checks_passed']}**.", '',
        f"Objective {summary['initial']['objective']:.9f} → {summary['selected']['objective']:.9f}. The exact selected raw eight-cell energies/forces and saved objective replay.", '',
        '| Fixture | La RMS / max | Dy RMS / max | Paired RMS / max | All pass |', '|---|---:|---:|---:|---|']
    for q in GEOMETRIES:
        rows = [selected['score']['per_fixture'][k][q] for k in ('La', 'Dy', 'paired')]
        text = [f"{r['rigid_RMS_kcal_mol_A']:.6f} / {r['maximum_body_component_error_kcal_mol_A']:.6f}" for r in rows]
        lines.append('| '+q+' | '+' | '.join(text)+' | '+str(all(r['passed'] for r in rows))+' |')
    lines += ['', 'Rigid response units: kcal/mol/Å, with the unchanged 1 Å torque divisor. Each fixture must satisfy RMS ≤2 and maximum component ≤5.', '',
        '| Energy control | Error (kcal/mol) | Pass ≤1 |', '|---|---:|---|']
    for block, values in selected['score']['energy_errors_kcal_mol'].items():
        for name, value in values.items(): lines.append(f'| {block}: {name} | {value:+.6f} | {abs(value)<=1} |')
    lines += ['', '| Loss block | Initial | Selected | Weight |', '|---|---:|---:|---:|']
    for k, value in selected['score']['blocks'].items(): lines.append(f"| {k} | {result['initial']['score']['blocks'][k]:.9f} | {value:.9f} | .25 |")
    lines += ['', '## Selected parameters and execution', '', '| Parameter | Start | Selected | Bound hit |', '|---|---:|---:|---|']
    for p in parameters: lines.append(f"| {p['name']} | {p['initial']:.9f} | {p['selected']:.9f} | {p['exact_bound_hit']} |")
    lines += ['', f"{result['new_parameter_vectors']}/384 additional vectors; {result['new_native_program_invocations']}/4096 additional native program invocations ({len(new_native)} new physical cells). {len(inherited_native)} inherited physical cells were reused with their original native authority. Original radius/damping bounds and declared positive epsilon intervals are unchanged.",
        f"Recorded wall {result['elapsed_seconds']:.3f} s; parent CPU {result['parent_process_CPU_seconds']:.3f} s; aggregate child user/system CPU {result['aggregate_native_child_user_CPU_seconds']:.3f}/{result['aggregate_native_child_system_CPU_seconds']:.3f} s. Concurrent receipt CPU deltas are not summed.", '',
        '## Prior center controls and interpretation', '',
        f"Original four-center objective: {center_score['objective']:.9f}; all prior center controls pass: {center_score['development_checks_passed']}. These remain separate from the expanded-domain objective.",
        f"Passing but unselected trial slots: {len(passing_unselected)}. Complete per-step controls, trial selection and remaining failed body responses are retained in REVIEW.json."]
    if iterations:
        last = iterations[-1]
        condition = 'unavailable (zero smallest singular value)' if last['origin_condition'] is None else f"{last['origin_condition']:.3f}"
        lines += ['', f"Last measured Jacobian (iteration {last['outer_step']} origin): rank {last['origin_Jacobian_rank']}, condition {condition}; unrestricted linear residual loss {last['origin_unrestricted_linear_residual_loss']:.9f}, locally projected fraction {last['origin_local_projected_fraction']:.3%}. Jacobian is at the final selected point: {summary['last_Jacobian_at_final_selected_point']}.",
            'The projection describes local weighted residual response only; it is not a tested next step or evidence of global family impossibility.']
    lines += ['', 'This is an uncorrected PBE0 representability experiment on consumed water fixtures. CC/CP method sensitivity, internal atomic forces, donor/protein transfer, bulk hydration and electronic-state accuracy remain separate. No physical qualification, new evaluation or automatic extension follows. Reserved validation remains unopened.', '',
        f"Result SHA256: `{summary['result']['sha256']}`", '']
    if event: (out/'COMPLETION_SNAPSHOT.json').write_text(json.dumps(event, indent=2)+'\n')
    summary['saved_review_elapsed_seconds'] = time.monotonic()-started
    summary['saved_review_parent_cpu_seconds'] = time.process_time()-cpu
    (out/'REVIEW.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    (out/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'review': pin(out/'REVIEW.json'), 'report': pin(out/'REPORT.md')}))
    return 0


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--result', required=True); p.add_argument('--output', required=True); p.add_argument('--completion-event')
    a = p.parse_args(); raise SystemExit(main(a.result, a.output, a.completion_event))
