#!/usr/bin/env python3
"""Corrected solvent evaluation and an explicit alternative-source overlay."""
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROTOCOL = 'Nikasha_two_Ca_source_corrected_GFN2_transfer_v1'
MXAF = ROOT/'diagnostics/ca_coverage_recovery_20261012/mxaf/SOLVENT.py'


def pin(path):
    path = Path(path).resolve(); digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1048576), b''): digest.update(chunk)
    return {'path': str(path), 'sha256': digest.hexdigest()}


def read(value):
    if isinstance(value, dict):
        assert pin(value['path']) == value, value['path']; value = value['path']
    return json.loads(Path(value).read_text())


def write(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(data, handle, indent=2, allow_nan=False); handle.write('\n')


def load(name, source):
    assert pin(source['path']) == source
    spec = importlib.util.spec_from_file_location(name, source['path'])
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def freeze(args):
    master = read(args.native_master)
    config = {'protocol_id': PROTOCOL, 'implementation': pin(__file__),
              'native_adapter': pin(HERE/'NATIVE.py'), 'execution_adapter': pin(MXAF), 'plan': pin(HERE/'PLAN.md'),
              'native_master': pin(args.native_master), 'native_collection': str(Path(args.native_collection).resolve()),
              'wrappers': master['wrappers'],
              'output': str(Path(args.output).resolve()),
              **{k: master[k] for k in ('preparation', 'scope', 'primary', 'prior_overlay', 'canonical_reference')},
              'declared_sources': master['declared_sources'], 'denominator': master['denominator'],
              'maximum_inputs_per_backend': 12*len(master['admitted_sources']),
              'new_calibration': False, 'new_triples': False, 'production_changed': False}
    write(args.config, config)
    print(json.dumps({'config': pin(args.config), 'native_collection_status': 'expected_path_checked_at_execution',
                      'new_molecular_calls': 0}))


def prepare(g, n, config, native, base):
    f, r = g.f, g.r
    master = read(config['native_master']); prepared = read(config['preparation']); scope = read(config['scope'])
    primary = read(config['primary']); old = read(primary['manifest'])
    assert native['terminal'] is True and native['manifest'] == config['native_master']
    assert native['declared_sources'] == config['declared_sources'] == master['declared_sources']
    assert native['denominator'] == config['denominator'] == len(prepared['rows'])
    assert native['preparation'] == master['preparation'] == config['preparation']
    assert native['canonical_reference'] == master['canonical_reference'] == config['canonical_reference']
    assert old['scope'] == config['scope'] and old['reference'] == primary['reference'] == scope['reference']
    assert read(config['canonical_reference'])['calibration'] == scope['reference']['frozen_calibration']
    assert config['canonical_reference'] == scope['reference']['source'] == prepared['canonical_reference']
    assert len(primary['rows']) == len(scope['rows']) == 125
    for lane in master['lanes']:
        worker = read(read(lane)['worker_manifest'])
        assert worker['model'] == old['model'] == scope['model']
        assert worker['pool_settings'] == old['pool_settings']
    native_rows = {c['source_case_id']: c for c in native['cases']}
    missing = {c['source_case_id']: c for c in native['missing_native_sources']}
    assert len(native_rows) == len(native['cases']) and set(native_rows).isdisjoint(missing)
    assert set(native_rows)|set(missing) == set(master['admitted_sources'])
    m = {k: copy.deepcopy(old[k]) for k in ('runtime', 'model', 'settings', 'pool_settings', 'reference', 'scope')}
    m.update(protocol_id=PROTOCOL, driver=f.pin(g.__file__), helpers=f.pin(g.RP), implementation=f.pin(f.__file__),
        adapter=pin(__file__), execution_adapter=config['execution_adapter'], recovery_config=config['_pin'],
        MACE_collection=pin(config['native_collection']), primary_collection=config['primary'],
        preparation=config['preparation'], native_master=config['native_master'], cases=[], tasks=[],
        evaluation_kind='new_two_Ca_conditioned_source_fixed_local_policy', execution_disabled=False,
        declared_full_source_denominator=125, recovery_source_denominator=len(prepared['rows']),
        new_calibration=False, new_Ca_triples=False)
    identities = {}
    for admission in prepared['rows']:
        sid = admission['case_id']; oid = admission['original_source_id']
        original = next(row for row in scope['rows'] if row['source_case_id'] == oid)
        primary_row = next(row for row in primary['rows'] if row['source_case_id'] == oid)
        assert primary_row['decisions']['operational']['outcome'] == 'unavailable'
        receipt = n.admitted(admission, scope)
        row = copy.deepcopy(original)
        row.update(source_case_id=sid, original_source_id=oid, source_conditioning_metal='2Ca_v1',
            case_id=sid+'__'+admission['selection_id'].rsplit('__', 1)[1],
            preparation_status=admission['status'], preparation_reason=admission.get('reason'),
            original_source_record=primary_row, new_source_admission=admission,
            execution_status='excluded', execution_reason=None)
        if receipt is None:
            row['execution_reason'] = 'new_source_preparation_unavailable'
            m['cases'].append(row); continue
        if sid not in native_rows:
            assert sid in missing
            row.update(execution_reason=missing[sid]['reason'], native_missing_record=missing[sid])
            m['cases'].append(row); continue
        case = native_rows[sid]; cid = row['case_id']
        assert case['case_id'] == cid and case['original_source_id'] == oid
        assert case['selection_id'] == admission['selection_id'] == original['selection_id']
        assert case['prepared_admission'] == receipt and case['source']['origin_row']['source'] == admission
        assert case['source']['origin_row']['mapping'] == admission['maps']
        assert case['source']['origin_row']['source_preparation'] == admission['representations']['context']['preparation']
        assert case['expected_class'] == original['expected_class'] and case['biological_group'] == original['biological_group']
        assert set(case['candidates']) == set(g.ORDER) and case['native_pool_status'] in ('available', 'unavailable')
        row.update(execution_status='new_component_tasks', prepared_admission=receipt,
            state_signature=admission['state_signature'], maps=admission['maps'],
            source_preparation=admission['representations']['context']['preparation'],
            native_pool_status=case['native_pool_status'], native_source_record=case,
            matrix={metal: {} for metal in ('Ca', 'La')}, paired_coordinate_max_error_A={})
        geometry = {}
        for metal in ('Ca', 'La'):
            charge = read(admission['state_signature'])['endpoint_charges'][metal]
            for candidate in g.ORDER:
                cell = g.verified_cell(case['matrix'][metal].get(candidate), metal, scope)
                row['matrix'][metal][candidate] = cell
                if cell is None: continue
                assert case['candidates'][candidate]['status'] == 'available' and cell['charge'] == charge
                if candidate == 'origin':
                    ep = admission['representations']['context']['endpoints'][metal]
                    assert all(cell[k] == ep[k] for k in ('xyz', 'charge', 'multiplicity'))
                geometry[metal, candidate] = f.coordinates(cell['xyz']); cell['inputs'] = {}
                for medium in ('vacuum', 'water'):
                    key = (cell['xyz']['sha256'], charge, 0, medium)
                    if key not in identities:
                        digest = hashlib.sha256(json.dumps([key, m['runtime'], {k: m['settings'][k] for k in r.SCIENTIFIC_SETTINGS}], sort_keys=True).encode()).hexdigest()
                        task = {'id': 'input_'+digest, 'protein': row['protein'], 'context': cid,
                            'geometry': candidate, 'metal': metal, 'medium': medium, 'xyz': cell['xyz'],
                            'charge': charge, 'multiplicity': 1, 'uhf': 0, 'valence_electrons': cell['valence_electrons'],
                            'atom_count': cell['atom_count'], 'reference': None, 'scout': False, 'reuse': None}
                        identities[key] = task['id']; m['tasks'].append(task)
                    cell['inputs'][medium] = identities[key]
        for candidate in g.ORDER:
            if all((metal, candidate) in geometry for metal in ('Ca', 'La')):
                ca, la = geometry['Ca', candidate], geometry['La', candidate]
                assert len(ca[0]) == len(la[0]) and ca[0][1:] == la[0][1:]
                residual = max(abs(a-b) for x, y in zip(ca[1], la[1]) for a, b in zip(x, y))
                assert residual <= m['pool_settings']['coordinate_copy_tolerance_A']
                row['paired_coordinate_max_error_A'][candidate] = residual
        for metal in ('Ca', 'La'):
            if (metal, 'origin') in geometry:
                for candidate in g.ORDER:
                    if (metal, candidate) in geometry: assert geometry[metal, candidate][0] == geometry[metal, 'origin'][0]
        m['cases'].append(row)
    assert len(m['tasks']) <= config['maximum_inputs_per_backend']
    assert [row['source_case_id'] for row in m['cases']] == config['declared_sources']
    base.mkdir(parents=True, exist_ok=False); f.write(base/'manifest.json', m)
    f.write(base/'PREPARATION.json', g.dry_run(base))
    f.write(base/'SOURCE_GATE.json', {'status': 'pass', 'new_molecular_calls': 0,
        'declared_new_sources': len(prepared['rows']), 'original125_unchanged': True,
        'frozen_reference': m['reference'], 'new_fit': False, 'new_triples': False})


def overlay(g, config, recovery_path, output):
    primary = read(config['primary']); prior = read(config['prior_overlay']); recovery = read(recovery_path)
    m = read(recovery['manifest']); declared = {r['source_case_id']: r for r in m['cases']}
    assert recovery['reference'] == primary['reference'] == prior['reference']
    assert prior['primary'] == config['primary'] and [r['primary'] for r in prior['rows']] == primary['rows']
    assert len(primary['rows']) == len(prior['rows']) == 125
    assert set(declared) == {r['source_case_id'] for r in recovery['rows']} == set(config['declared_sources'])
    by_old = {}
    for result in recovery['rows']:
        metadata = declared[result['source_case_id']]; old = metadata['original_source_id']
        assert old not in by_old
        by_old[old] = {**result, 'original_source_id': old,
            'execution_reason': metadata.get('execution_reason'),
            'preparation_record': metadata['new_source_admission'],
            'native_missing_record': metadata.get('native_missing_record')}
    rows = []
    for item in prior['rows']:
        row = copy.deepcopy(item); oid = item['primary']['source_case_id']
        if oid in by_old:
            assert item['primary']['decisions']['operational']['outcome'] == 'unavailable'
            assert item.get('technical_recovery') is None and item.get('preparation_recovery') is None
        row['two_Ca_source_recovery'] = by_old.get(oid); rows.append(row)
    selected = [r['two_Ca_source_recovery'] or r.get('preparation_recovery') or r.get('technical_recovery') or r['primary'] for r in rows]
    result = {'protocol_id': PROTOCOL+'_alternative_source_overlay', 'primary': config['primary'],
        'prior_recovery_overlay': config['prior_overlay'], 'two_Ca_collection': pin(recovery_path),
        'rows': rows, 'declared_sources': 125, 'declared_new_two_Ca_sources': len(by_old),
        'primary_counts': primary['counts'], 'prior_recovery_counts': prior['annotated_recovery_counts'],
        'annotated_alternative_source_counts': {policy: dict(Counter(r['decisions'][policy]['outcome'] for r in selected)) for policy in ('operational', 'mathematical', 'origin')},
        'reference': primary['reference'], 'new_fit': False, 'new_triples': False, 'production_changed': False,
        'interpretation': 'Original125 and prior recoveries preserved. New two-Ca folds are separately identified alternative source hypotheses; no source ranking, affinity claim, or independent-protein inflation.'}
    g.f.write(output, result); return result


def run(path):
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('Actual Slurm required before scientific imports')
    config = read(path); config['_pin'] = pin(path); assert config['implementation'] == pin(__file__)
    for key in ('scope', 'primary', 'prior_overlay', 'preparation', 'native_master', 'canonical_reference'): read(config[key])
    assert pin(config['plan']['path']) == config['plan']
    for value in config['wrappers'].values(): assert pin(value['path']) == value
    adapter = load('two_ca_existing_solvent', config['execution_adapter'])
    stage, g, placement = adapter.stage()
    n = load('two_ca_native_contract', config['native_adapter'])
    base = Path(config['output']); base.parent.mkdir(parents=True, exist_ok=True)
    started = time.time(); receipt = {'config': config['_pin'], 'placement': placement, 'status': 'started', 'started_epoch': started}
    try:
        with (base.parent/('.'+base.name+'.RECOVERY.lock')).open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            if not (base/'manifest.json').exists(): prepare(g, n, config, read(config['native_collection']), base)
            manifest = g.f.read(base/'manifest.json')
            assert manifest['recovery_config'] == config['_pin'] and manifest['adapter'] == pin(__file__)
            assert manifest['MACE_collection'] == pin(config['native_collection'])
            assert read(base/'SOURCE_GATE.json')['status'] == 'pass'
            result = stage.run(Path(config['scope']['path']), Path(config['native_collection']), base)
            collection = base/'FINAL_COLLECTION.json'
            if read(collection) != result:
                matches = [p for p in base.glob('SLURM_ATTEMPT_*_COLLECTION.json') if read(p) == result]
                assert len(matches) == 1; collection = matches[0]
            dest = base/(collection.stem+'_ALTERNATIVE_SOURCE_OVERLAY.json')
            combined = overlay(g, config, collection, dest)
            receipt.update(status='complete', collection=pin(collection), overlay=pin(dest), actual_backend_cost=result['cost'],
                           annotated_counts=combined['annotated_alternative_source_counts'])
    except BaseException as exc:
        receipt.update(status='failed', reason=str(exc), traceback=traceback.format_exc()); raise
    finally:
        finished = time.time(); receipt.update(finished_epoch=finished, elapsed_seconds_including_preparation=finished-started,
                                              allocated_CPU_seconds_including_preparation=12*(finished-started))
        g.f.write(base.parent/(base.name+'_EXECUTION_'+os.environ['SLURM_JOB_ID']+'.json'), receipt)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='op', required=True)
    p = sub.add_parser('freeze')
    for field in ('native-master', 'native-collection', 'output', 'config'): p.add_argument('--'+field, required=True)
    p = sub.add_parser('run'); p.add_argument('--config', required=True)
    args = parser.parse_args()
    if args.op == 'freeze': freeze(args)
    else: print(json.dumps(run(args.config), indent=2))
