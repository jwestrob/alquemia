#!/usr/bin/env python3
"""Finite two-Ca source transport around unchanged, qualified native kernels."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DISULFIDE = ROOT/'diagnostics/ca_coverage_recovery_20261012/disulfide/NATIVE.py'
PROTOCOL = 'Nikasha_two_Ca_sources_original_four_mode_CPU_v1'
POLICIES = ('operational', 'mathematical', 'origin')


def pin(path):
    p = Path(path).resolve(); h = hashlib.sha256()
    with p.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1048576), b''): h.update(chunk)
    return {'path': str(p), 'sha256': h.hexdigest()}


def verify(value):
    assert pin(value['path']) == value, value['path']
    return Path(value['path'])


def read(value):
    return json.loads((verify(value) if isinstance(value, dict) else Path(value)).read_text())


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def state_anchor(row):
    data = read(row['state_anchor'])
    if data.get('case_id') == row['anchor_case_id']:
        anchor = data
    else:
        # Curie's actual pinned canonical envelope PREPARATION uses `cases`.
        records = data['cases'] if 'cases' in data else data['rows']
        matches = [r for r in records if r['case_id'] == row['anchor_case_id']
                   and r['selection_id'] == row['selection_id']]
        assert len(matches) == 1; anchor = matches[0]
    assert anchor['status'] == 'prepared'
    assert row['state_key'] == row['frozen_state_key'] == anchor['state_key']
    return row['state_key']


def admitted(row, scope):
    """Check recorded admission identities; no molecular parsing here."""
    old = next(r for r in scope['rows'] if r['source_case_id'] == row['original_source_id'])
    assert row['case_id'] != row['original_source_id'] and '__conditioned_2Ca_v1__' in row['case_id']
    assert row['selection_id'] == old['selection_id']
    assert row['source']['expected_class'] == old['expected_class']
    assert row['source']['biological_group'] == old['biological_group']
    assert old['biological_group'] in ('q9z4j7', 'c5axv8', 'q88jh5')
    if row['status'] != 'prepared': return None
    assert row['mapping_eligible'] and row['paired_maps_equal']
    audit = read(row['spectator_admission'])
    assert audit['status'] == 'passed' and audit['passed'] is True
    state_anchor(row)
    for field in ('state_signature', 'state_anchor', 'spectator_admission'): verify(row[field])
    for value in row['maps'].values(): verify(value)
    verify(row['representations']['context']['preparation'])
    admission = pin(Path(row['maps']['Ca']['path']).parent/'ADMISSION.json')
    assert read(admission) == row
    return admission


def prepare(args):
    """Bind only real complete preparation records; this operation is metadata-only."""
    prepared = read(args.preparation); scope = read(args.scope); old = read(args.template)
    primary = read(args.primary); prior = read(args.prior_overlay)
    rows = prepared['rows']; ids = [r['case_id'] for r in rows]
    assert prepared['terminal'] is True and prepared['denominator'] == len(rows)
    assert len(rows) in (5, 10, 15) and len(set(ids)) == len(rows)
    assert len({r['original_source_id'] for r in rows}) == len(rows)
    assert len(primary['rows']) == len(scope['rows']) == 125
    assert prepared['canonical_reference'] == scope['reference']['source']
    assert read(prepared['canonical_reference'])['calibration'] == scope['reference']['frozen_calibration']
    assert primary['reference'] == scope['reference'] == prior['reference']
    assert prior['primary'] == pin(args.primary) and len(prior['rows']) == 125
    assert [r['primary'] for r in prior['rows']] == primary['rows']
    assert old['model'] == scope['model']
    eligible = [(row, admitted(row, scope)) for row in rows]
    eligible = [(row, receipt) for row, receipt in eligible if receipt is not None]
    if eligible: assert read(prepared['regression'])['status'] == 'passed'
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False)
    shared = {'protocol_id': PROTOCOL, 'implementation': pin(__file__), 'helper': pin(DISULFIDE),
              'plan': pin(HERE/'PLAN.md'), 'template': pin(args.template),
              'wrappers': {name: pin(HERE/name) for name in ('native.sbatch', 'collect.sbatch', 'solvent.sbatch')},
              'preparation': pin(args.preparation), 'scope': pin(args.scope),
              'primary': pin(args.primary), 'prior_overlay': pin(args.prior_overlay),
              'canonical_reference': prepared['canonical_reference']}
    fixtures = [next(t for t in old['tasks'] if t['task_id'] == f['task_id']) for f in old['fixtures']]
    count = min(args.lanes, len(eligible)); lanes = []
    for index in range(count):
        assigned = eligible[index::count]; base = out/('lane'+str(index)); base.mkdir()
        tasks = []; cases = []
        for row, receipt in assigned:
            cid = row['case_id']+'__'+row['selection_id'].rsplit('__', 1)[1]
            for metal in ('Ca', 'La'):
                endpoint = row['representations']['context']['endpoints'][metal]
                tasks.append({'task_id': cid+'__context__'+metal, 'case_id': cid,
                    'source_case_id': row['case_id'], 'selection_id': row['selection_id'], 'metal': metal,
                    **{k: endpoint[k] for k in ('xyz', 'charge', 'multiplicity')},
                    'mapping': row['maps'][metal], 'state_signature': row['state_signature'],
                    'context_preparation': row['representations']['context']['preparation'], 'prepared_admission': receipt})
            cases.append({'case_id': cid, 'source_case_id': row['case_id'], 'original_source_id': row['original_source_id'],
                'selection_id': row['selection_id'], 'lane': 0, 'prepared_admission': receipt,
                'state_signature': row['state_signature'], 'source_preparation': row['representations']['context']['preparation']})
        worker = copy.deepcopy(old)
        for field in ('scientific_scope', 'origin_inputs'): worker.pop(field, None)
        worker.update(protocol_id=PROTOCOL+'_transport', parent_manifest=pin(args.template),
            preparation=pin(args.preparation), plan=shared['plan'], tasks=tasks+fixtures,
            origin_task_count=len(tasks), cases=cases, lanes=1, maximum_active_numerical_threads=7,
            canonical_reference=prepared['canonical_reference'],
            finite_scope={'sources': len(cases), 'new_origin_calls': 2*len(cases), 'optimizer_starts': 2*len(cases),
                          'maximum_cross_calls': 2*len(cases), 'qualification_calls': 2, 'GFN2_calls': 0, 'DFT_calls': 0},
            authorization='User-authorized finite two-Ca source recovery; actual Slurm only; unchanged four-mode kernel')
        assert len(tasks) <= 30 and {t['case_id'] for t in fixtures}.isdisjoint(c['case_id'] for c in cases)
        write(base/'worker_manifest.json', worker)
        plan = {**shared, 'lane_index': index, 'worker_manifest': pin(base/'worker_manifest.json'),
                'resources': {'CPUs': 7, 'model_workers': 1, 'model_threads': 7,
                              'parent_threads': 1, 'parent_and_model_numerics_sequential': True},
                'source_ids': [c['source_case_id'] for c in cases]}
        write(base/'manifest.json', plan); lanes.append(pin(base/'manifest.json'))
    master = {**shared, 'lanes': lanes, 'declared_sources': ids, 'denominator': len(rows),
              'admitted_sources': [r['case_id'] for r, _ in eligible],
              'unavailable_preparations': [r for r in rows if r['status'] != 'prepared'],
              'new_fit': False, 'new_triples': False, 'production_changed': False}
    write(out/'manifest.json', master)
    write(out/'READY.json', {'status': 'metadata_prepared_scientific_preflight_UNRUN',
                            'manifest': pin(out/'manifest.json'), 'lanes': lanes, 'new_molecular_calls': 0})
    print(json.dumps({'manifest': pin(out/'manifest.json'), 'lanes': lanes, 'denominator': len(rows)}))


def execute(path):
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('Actual Slurm required before scientific imports')
    plan = read(path); assert plan['implementation'] == pin(__file__)
    helper = load('two_ca_existing_native', verify(plan['helper'])); placement = helper.allocation()
    mp = verify(plan['worker_manifest']); m = read(mp); base = Path(path).resolve().parent
    for key in ('plan', 'preparation', 'template', 'scope', 'canonical_reference'): verify(plan[key])
    for value in plan['wrappers'].values(): verify(value)
    with (base/'RUN.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        attempt = base/'attempts'/'001'; attempt.mkdir(parents=True, exist_ok=False); (attempt/'lane0').mkdir()
        write(attempt/'PLACEMENT.json', placement)
        b = load('two_ca_frozen_batch', verify(m['implementation'])); b.modules(m)
        old_handlers = b.handlers(); transport = None; error = None; started = time.monotonic()
        try:
            checks = []
            for case in m['cases']:
                row = read(case['prepared_admission']); admitted(row, read(plan['scope']))
                view = copy.deepcopy(m)
                view.update(cases=[case], origin_task_count=2,
                    tasks=[t for t in m['tasks'][:m['origin_task_count']] if t['case_id'] == case['case_id']]+m['tasks'][m['origin_task_count']:])
                local = {**plan, 'target_admission': case['prepared_admission']}
                checks.append({'source_case_id': case['source_case_id'], 'result': helper.preflight(local, view, b)})
            write(attempt/'PREFLIGHT.json', {'status': 'pass', 'cases': checks, 'new_model_calls': 0})
            transport = b.WarmCPU(mp, m, attempt, 0)
            import numpy as np
            checks = []
            for fixture in m['fixtures']:
                task = next(t for t in m['tasks'] if t['task_id'] == fixture['task_id'])
                request = {'manifest': pin(mp), 'task_id': task['task_id'], 'xyz': fixture['xyz'],
                           'charge': task['charge'], 'multiplicity': task['multiplicity']}
                result, rp = b.cell(attempt/'lane0'/'qualification'/fixture['id'], request, transport, m['model'])
                check = {'fixture': fixture['id'], 'receipt': rp, 'pass': False}
                if result['status'] == 'complete':
                    delta = (result['energy_eV']-fixture['energy_eV'])*b.EV_TO_KCAL
                    force = float(np.max(np.abs(np.load(verify(result['forces']), allow_pickle=False)-np.load(verify(fixture['forces']), allow_pickle=False))))
                    check.update(energy_error_kcal_mol=delta, force_max_error_eV_A=force,
                                 **{'pass': abs(delta) <= m['gates']['energy_kcal_mol'] and force <= m['gates']['force_max_eV_A']})
                checks.append(check)
            write(attempt/'QUALIFICATION.json', {'checks': checks, 'pass': all(c['pass'] for c in checks)})
            assert all(c['pass'] for c in checks), 'existing native force qualification failed'
            for case in m['cases']: b.run_source(mp, m, case, transport)
        except BaseException as exc:
            error = type(exc).__name__+': '+str(exc); traceback.print_exc()
        finally:
            if transport:
                try:
                    if error: transport.abort()
                    else: transport.close()
                except BaseException as exc:
                    error = (error+'; ' if error else '')+str(exc); transport.abort()
            for sig, handler in old_handlers.items(): signal.signal(sig, handler)
            paths = [base/'sources'/c['case_id']/'SOURCE_RESULT.json' for c in m['cases']]
            cases = [read(p) for p in paths if p.exists()]
            calls = [read(p) for p in sorted(base.glob('**/mace_result.json'))]
            elapsed = time.monotonic()-started
            write(base/'FINAL_COLLECTION.json', {'protocol_id': PROTOCOL, 'manifest': pin(mp),
                'recovery_plan': pin(path), 'terminal': len(cases) == len(m['cases']), 'denominator': len(m['cases']),
                'cases': cases, 'source_results': [pin(p) for p in paths if p.exists()], 'error': error,
                'missing_source_ids': [c['source_case_id'] for c, p in zip(m['cases'], paths) if not p.exists()],
                'costs': {'new_recorded_model_calls': len(calls), 'actual_forwards': sum(r['actual_model_forwards'] for r in calls),
                    'recorded_wall_seconds_sum': sum(r['wall_seconds'] for r in calls),
                    'recorded_process_CPU_seconds': sum(r['process_cpu_seconds'] for r in calls),
                    'started_without_terminal_receipt': [str(p) for p in base.glob('**/CPU_CALL_STARTED.json') if not (p.parent/'mace_result.json').exists()],
                    'elapsed_seconds': elapsed, 'allocated_CPU_seconds': elapsed*7, 'placement': placement},
                'canonical_reference': plan['canonical_reference'], 'new_GFN2_calls': 0, 'new_DFT_calls': 0,
                'original_source_results_preserved': True, 'production_changed': False})
        if error: raise RuntimeError(error)


def collect(path, output):
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('Actual Slurm required for collection')
    master = read(path); assert master['implementation'] == pin(__file__)
    verify(master['plan'])
    for value in master['wrappers'].values(): verify(value)
    lanes = []; cases = []; missing = []
    for lp in master['lanes']:
        lane = read(lp); worker = read(lane['worker_manifest']); base = Path(lp['path']).parent
        final = base/'FINAL_COLLECTION.json'
        entry = {'plan': lp, 'collection': pin(final) if final.exists() else None}
        if final.exists():
            prior = read(final); assert prior['manifest'] == lane['worker_manifest']
            entry.update(error=prior['error'], costs=prior['costs'])
        for case in worker['cases']:
            source = base/'sources'/case['case_id']/'SOURCE_RESULT.json'
            if source.exists():
                result = read(source); assert result['batch_manifest'] == lane['worker_manifest']
                assert result['source_case_id'] == case['source_case_id']
                cases.append(result)
            else:
                missing.append({'source_case_id': case['source_case_id'], 'original_source_id': case['original_source_id'],
                                'reason': 'declared_lane_has_no_terminal_source_record', 'lane': lp,
                                'preserved_source_directory': str(source.parent)})
        lanes.append(entry)
    assert {c['source_case_id'] for c in cases}.isdisjoint(r['source_case_id'] for r in missing)
    assert {c['source_case_id'] for c in cases}|{r['source_case_id'] for r in missing} == set(master['admitted_sources'])
    result = {'protocol_id': PROTOCOL, 'manifest': pin(path), 'terminal': True,
        'terminal_meaning': 'finite collection snapshot after dependent lane jobs; missing chemistry remains unavailable',
        'denominator': master['denominator'], 'declared_sources': master['declared_sources'], 'cases': cases,
        'missing_native_sources': missing, 'unavailable_preparations': master['unavailable_preparations'],
        'lanes': lanes, 'preparation': master['preparation'], 'canonical_reference': master['canonical_reference'],
        'new_molecular_calls_in_collection': 0, 'production_changed': False}
    write(output, result); return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='op', required=True)
    p = sub.add_parser('prepare')
    for field in ('template', 'preparation', 'scope', 'primary', 'prior-overlay', 'output'): p.add_argument('--'+field, required=True)
    p.add_argument('--lanes', type=int, choices=(2, 3), required=True)
    e = sub.add_parser('execute'); e.add_argument('--manifest', required=True)
    c = sub.add_parser('collect'); c.add_argument('--manifest', required=True); c.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.op == 'prepare': prepare(args)
    elif args.op == 'execute': execute(args.manifest)
    else: print(json.dumps(collect(args.manifest, args.output), indent=2))
