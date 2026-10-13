"""One repaired-source solvent pool; preserve the complete primary denominator."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
import traceback
from collections import Counter

PROTOCOL = 'Nikasha_MMOL2048_CHrepair_corrected_GFN2_recovery_v1'
TARGET = 'mmol_2048-pqq-la_model__conditioned_Ca__seed-1_sample-0'


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(path):
    if isinstance(path, dict):
        assert pin(path['path']) == path, path['path']
        path = path['path']
    return json.loads(Path(path).read_text())


def prepare(g, config, native, base):
    f, r = g.f, g.r
    scope = read(config['scope']); primary = read(config['primary']); old = read(primary['manifest'])
    plan = read(config['native_plan']); worker = read(plan['worker_manifest'])
    prepared = read(config['preparation']); admission = read(config['target_admission'])
    assert native['terminal'] and native['denominator'] == 1 and len(native['cases']) == 1
    assert native['manifest'] == plan['worker_manifest'] and native['recovery_plan'] == config['native_plan']
    assert worker['preparation'] == plan['preparation'] == config['preparation']
    assert prepared['status'] == 'prepared' and prepared['target_source'] == TARGET
    assert admission == next(row for row in prepared['rows'] if row['case_id'] == TARGET)
    assert admission['status'] == 'prepared' and admission['state_key'] == prepared['frozen_state_anchor']
    assert admission['mapping_eligible'] and admission['paired_maps_equal'] and len(prepared['controls']) == 2
    assert all(row['status'] == 'passed' for row in prepared['controls'])
    assert worker['model'] == scope['model'] == old['model'] and worker['pool_settings'] == old['pool_settings']
    assert old['scope'] == config['scope'] and old['reference'] == primary['reference'] == scope['reference']
    assert read(config['canonical_reference'])['calibration'] == scope['reference']['frozen_calibration']
    assert prepared['canonical_reference'] == config['canonical_reference'] == scope['reference']['source']
    assert len(primary['rows']) == len(scope['rows']) == 125
    original = next(row for row in scope['rows'] if row['source_case_id'] == TARGET)
    primary_row = next(row for row in primary['rows'] if row['source_case_id'] == TARGET)
    assert original['execution_status'] == 'excluded' and primary_row['decisions']['operational']['outcome'] == 'unavailable'
    case = native['cases'][0]; cid = case['case_id']
    assert case['source_case_id'] == TARGET and case['selection_id'] == admission['selection_id'] == original['selection_id']
    assert case['prepared_admission'] == config['target_admission'] and case['source']['origin_row']['source'] == admission
    assert case['source']['origin_row']['source_preparation'] == admission['representations']['context']['preparation']
    assert case['source']['origin_row']['mapping'] == admission['maps']
    assert case['expected_class'] == original['expected_class'] and case['biological_group'] == original['biological_group']
    assert set(case['candidates']) == set(g.ORDER) and case['native_pool_status'] in ('available', 'unavailable')
    m = {k: copy.deepcopy(old[k]) for k in ('runtime', 'model', 'settings', 'pool_settings', 'reference', 'scope')}
    m.update(protocol_id=PROTOCOL, driver=f.pin(g.__file__), helpers=f.pin(g.RP), implementation=f.pin(f.__file__),
        adapter=pin(__file__), execution_adapter=config['execution_adapter'], recovery_config=config['_pin'],
        MACE_collection=pin(config['native_collection']), primary_collection=config['primary'],
        repaired_preparation=config['preparation'], target_admission=config['target_admission'],
        cases=[], tasks=[], evaluation_kind='same_state_source_hydrogen_preparation_recovery', execution_disabled=False,
        declared_full_source_denominator=125, recovery_source_denominator=1, new_calibration=False, new_Ca_triples=False)
    row = copy.deepcopy(original)
    row.update(case_id=cid, preparation_status='prepared', preparation_reason=None, execution_status='new_component_tasks',
        original_exclusion=primary_row, prepared_admission=config['target_admission'], state_signature=admission['state_signature'],
        maps=admission['maps'], source_preparation=admission['representations']['context']['preparation'],
        native_pool_status=case['native_pool_status'], native_source_record=case,
        matrix={z: {} for z in ('Ca', 'La')}, paired_coordinate_max_error_A={})
    identities = {}; geometry = {}
    for z in ('Ca', 'La'):
        charge = read(admission['state_signature'])['endpoint_charges'][z]
        for q in g.ORDER:
            cell = g.verified_cell(case['matrix'][z].get(q), z, scope); row['matrix'][z][q] = cell
            if cell is None: continue
            assert case['candidates'][q]['status'] == 'available' and cell['charge'] == charge
            if q == 'origin':
                ep = admission['representations']['context']['endpoints'][z]
                assert all(cell[k] == ep[k] for k in ('xyz', 'charge', 'multiplicity'))
            geometry[z, q] = f.coordinates(cell['xyz']); cell['inputs'] = {}
            for medium in ('vacuum', 'water'):
                key = (cell['xyz']['sha256'], charge, 0, medium)
                if key not in identities:
                    digest = hashlib.sha256(json.dumps([key, m['runtime'], {k: m['settings'][k] for k in r.SCIENTIFIC_SETTINGS}], sort_keys=True).encode()).hexdigest()
                    task = {'id': 'input_'+digest, 'protein': row['protein'], 'context': cid, 'geometry': q,
                            'metal': z, 'medium': medium, 'xyz': cell['xyz'], 'charge': charge, 'multiplicity': 1,
                            'uhf': 0, 'valence_electrons': cell['valence_electrons'], 'atom_count': cell['atom_count'],
                            'reference': None, 'scout': False, 'reuse': None}
                    identities[key] = task['id']; m['tasks'].append(task)
                cell['inputs'][medium] = identities[key]
    for q in g.ORDER:
        if all((z, q) in geometry for z in ('Ca', 'La')):
            ca, la = geometry['Ca', q], geometry['La', q]
            assert len(ca[0]) == len(la[0]) and ca[0][1:] == la[0][1:]
            residual = max(abs(a-b) for x, y in zip(ca[1], la[1]) for a, b in zip(x, y))
            assert residual <= m['pool_settings']['coordinate_copy_tolerance_A']
            row['paired_coordinate_max_error_A'][q] = residual
    for z in ('Ca', 'La'):
        if (z, 'origin') not in geometry: continue
        for q in g.ORDER:
            if (z, q) in geometry: assert geometry[z, q][0] == geometry[z, 'origin'][0]
    m['cases'].append(row)
    assert len(m['tasks']) <= 12 and len({t['id'] for t in m['tasks']}) == len(m['tasks'])
    base.mkdir(parents=True, exist_ok=False); f.write(base/'manifest.json', m)
    f.write(base/'PREPARATION.json', g.dry_run(base))
    f.write(base/'SOURCE_GATE.json', {'status': 'pass', 'new_molecular_calls': 0, 'primary_row_preserved': primary_row,
        'prepared_admission': config['target_admission'], 'frozen_state_anchor': prepared['frozen_state_anchor'],
        'maximum_fresh_inputs_per_backend': 12, 'frozen_reference': m['reference'], 'new_fit': False})


def overlay(g, config, recovery_path, output):
    f = g.f; primary = read(config['primary']); recovery = read(recovery_path)
    assert len(primary['rows']) == 125 and len(recovery['rows']) == 1
    new = recovery['rows'][0]; assert new['source_case_id'] == TARGET and primary['reference'] == recovery['reference']
    prior_technical = {}; status = 'not_available_at_collection'; mxaf_pin = None
    optional = Path(config['optional_MxaF_overlay'])
    if optional.exists():
        prior = read(optional); mxaf_pin = pin(optional)
        assert prior['primary'] == config['primary'] and prior['reference'] == primary['reference']
        assert len(prior['rows']) == 125 and [r['primary'] for r in prior['rows']] == primary['rows']
        actual = read(prior['technical_recovery']); assert actual['reference'] == primary['reference']
        actual_rows = {r['source_case_id']: r for r in actual['rows']}
        for item in prior['rows']:
            r = item['technical_recovery']
            if r is not None:
                assert r == actual_rows[r['source_case_id']] and r['source_case_id'] != TARGET
                prior_technical[r['source_case_id']] = r
        assert set(prior_technical) == {f'p12293-pqq-la_model__conditioned_Ca__seed-1_sample-{i}' for i in (1, 4)}
        for key, recovered in prior_technical.items():
            old = next(r for r in primary['rows'] if r['source_case_id'] == key)
            assert old['decisions']['operational']['outcome'] == 'unavailable' and old['expected_class'] == recovered['expected_class']
        status = 'included_actual_pinned_MxaF_recovery'
    rows = [{'primary': old, 'technical_recovery': prior_technical.get(old['source_case_id']),
             'preparation_recovery': new if old['source_case_id'] == TARGET else None} for old in primary['rows']]
    selected = [r['preparation_recovery'] or r['technical_recovery'] or r['primary'] for r in rows]
    result = {'protocol_id': PROTOCOL+'_primary_preserving_overlay', 'primary': config['primary'],
        'preparation_recovery': pin(recovery_path), 'MxaF_overlay': mxaf_pin, 'MxaF_status': status,
        'rows': rows, 'declared_sources': 125, 'primary_counts': primary['counts'],
        'annotated_recovery_counts': {p: dict(Counter(r['decisions'][p]['outcome'] for r in selected)) for p in ('operational', 'mathematical', 'origin')},
        'reference': primary['reference'], 'new_fit': False, 'production_changed': False,
        'interpretation': 'Immutable primary results plus separately named same-policy technical and source-preparation recoveries; no source omitted.'}
    f.write(output, result); return result


def run(path):
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('All scientific operations require actual Slurm')
    config = read(path); config['_pin'] = pin(path)
    for name in ('implementation', 'execution_adapter', 'scope', 'primary', 'preparation', 'native_plan', 'target_admission', 'canonical_reference', 'plan'):
        expected = config[name]
        assert pin(expected['path']) == expected, name
    assert config['implementation'] == pin(__file__)
    spec = importlib.util.spec_from_file_location('existing_mxaf_solvent_transport', config['execution_adapter']['path'])
    adapter = importlib.util.module_from_spec(spec); spec.loader.exec_module(adapter)
    stage, g, placement = adapter.stage()  # Existing real12CPU allocation check before scientific imports.
    base = Path(config['output']); base.parent.mkdir(parents=True, exist_ok=True)
    started = time.time(); receipt = {'config': config['_pin'], 'placement': placement, 'status': 'started', 'started_epoch': started}
    try:
        with (base.parent/('.'+base.name+'.RECOVERY.lock')).open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            if not (base/'manifest.json').exists(): prepare(g, config, read(config['native_collection']), base)
            m = g.f.read(base/'manifest.json')
            assert m['recovery_config'] == config['_pin'] and m['adapter'] == pin(__file__)
            assert m['MACE_collection'] == pin(config['native_collection']) and m['scope'] == config['scope']
            assert read(base/'SOURCE_GATE.json')['status'] == 'pass'
            result = stage.run(Path(config['scope']['path']), Path(config['native_collection']), base)
            collection = base/'FINAL_COLLECTION.json'
            if read(collection) != result:
                matches = [p for p in base.glob('SLURM_ATTEMPT_*_COLLECTION.json') if read(p) == result]
                assert len(matches) == 1; collection = matches[0]
            name = collection.stem+'_PRIMARY_RECOVERY_OVERLAY.json'
            combined = overlay(g, config, collection, base/name)
            receipt.update(status='complete', collection=pin(collection), overlay=pin(base/name), actual_backend_cost=result['cost'],
                           annotated_counts=combined['annotated_recovery_counts'])
    except BaseException as exc:
        receipt.update(status='failed', reason=str(exc), traceback=traceback.format_exc()); raise
    finally:
        receipt.update(finished_epoch=time.time(), elapsed_seconds_including_preparation=time.time()-started,
                       allocated_CPU_seconds_including_preparation=12*(time.time()-started))
        g.f.write(base.parent/(base.name+'_EXECUTION_'+os.environ['SLURM_JOB_ID']+'.json'), receipt)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--config', required=True)
    args = parser.parse_args(); print(json.dumps(run(args.config), indent=2))
