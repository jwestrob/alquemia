"""One-source CPU transport; call the pinned existing source/optimizer kernels."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time
import traceback

PROTOCOL = 'Nikasha_disulfide_repaired_source_original_four_mode_CPU_v1'
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(value):
    assert pin(value['path']) == value, value['path']
    return Path(value['path'])


def read(path):
    return json.loads((verify(path) if isinstance(path, dict) else Path(path)).read_text())


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def prepare(template, preparation, output):
    """Metadata-only task binding. Scientific preflight is deferred to Slurm."""
    old = read(template); prepared = read(preparation)
    assert prepared['status'] == 'prepared' and len(prepared['controls']) == 2
    row = next(r for r in prepared['rows'] if r['case_id'] == prepared['target_source'])
    assert row['status'] == 'prepared' and row['state_key'] == prepared['frozen_state_anchor']
    assert row['paired_maps_equal'] and row['mapping_eligible']
    cid = row['case_id'] + '__' + row['selection_id'].rsplit('__', 1)[1]
    admission = pin(Path(row['maps']['Ca']['path']).parent/'ADMISSION.json')
    assert read(admission) == row
    tasks = []
    for metal in ('Ca', 'La'):
        ep = row['representations']['context']['endpoints'][metal]
        tasks.append({'task_id': cid+'__context__'+metal, 'case_id': cid,
            'source_case_id': row['case_id'], 'selection_id': row['selection_id'], 'metal': metal,
            **{k: ep[k] for k in ('xyz', 'charge', 'multiplicity')}, 'mapping': row['maps'][metal],
            'state_signature': row['state_signature'], 'context_preparation': row['representations']['context']['preparation'],
            'prepared_admission': admission})
    fixtures = [next(t for t in old['tasks'] if t['task_id'] == f['task_id']) for f in old['fixtures']]
    manifest = copy.deepcopy(old)
    for key in ('scientific_scope', 'origin_inputs'): manifest.pop(key, None)
    manifest.update(protocol_id=PROTOCOL+'_transport', parent_manifest=pin(template), preparation=pin(preparation),
        plan=pin(Path(__file__).with_name('NATIVE_PLAN.md')), tasks=tasks+fixtures, origin_task_count=2,
        cases=[{'case_id': cid, 'source_case_id': row['case_id'], 'selection_id': row['selection_id'], 'lane': 0,
                'prepared_admission': admission, 'state_signature': row['state_signature'],
                'source_preparation': row['representations']['context']['preparation']}],
        lanes=1, maximum_active_numerical_threads=7,
        finite_scope={'sources': 1, 'new_origin_calls': 2, 'optimizer_starts': 2, 'maximum_cross_calls': 2,
                      'qualification_calls': 2, 'DFT_calls': 0, 'GFN2_calls': 0},
        authorization='User-authorized unavailable-source recovery; actual Slurm allocation only; original four-mode proposal policy',
        canonical_reference=prepared['canonical_reference'])
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    write(out/'worker_manifest.json', manifest)
    plan = {'protocol_id': PROTOCOL, 'implementation': pin(__file__), 'worker_manifest': pin(out/'worker_manifest.json'),
            'template': pin(template), 'preparation': pin(preparation), 'target_admission': admission,
            'canonical_reference': prepared['canonical_reference'], 'plan': manifest['plan'],
            'resources': {'CPUs': 7, 'model_workers': 1, 'model_threads': 7, 'parent_threads': 1,
                          'parent_and_model_numerical_work_sequential': True},
            'new_controls_prepared_not_rescored': 2, 'scope': manifest['finite_scope'],
            'production_changed': False, 'old_unavailable_result_preserved': True}
    write(out/'manifest.json', plan)
    write(out/'PREPARATION.json', {'status': 'metadata_prepared_scientific_preflight_UNRUN',
                                 'manifest': pin(out/'manifest.json'), 'new_molecular_calls': 0})
    print(json.dumps({'manifest': pin(out/'manifest.json'), 'scope': plan['scope']}))


def allocation():
    job = os.environ.get('SLURM_JOB_ID'); assert job and job.isdigit(), 'actual Slurm allocation required'
    text = subprocess.run(['scontrol', 'show', 'job', '-o', job], capture_output=True, text=True, check=True).stdout
    fields = dict(t.split('=', 1) for t in text.split() if '=' in t)
    assert fields['JobState'] == 'RUNNING' and int(fields['NumCPUs']) == 7 and int(fields['NumNodes']) == 1
    assert int(os.environ['SLURM_CPUS_ON_NODE']) == 7 and len(os.sched_getaffinity(0)) >= 7
    assert all(os.environ.get(k) == '1' for k in THREADS)
    return {'job_id': job, 'hostname': socket.gethostname(), 'scontrol': text,
            'allocated_CPUs': 7, 'affinity': sorted(os.sched_getaffinity(0))}


def preflight(plan, m, b):
    import numpy as np
    import scipy
    import scipy.optimize._slsqplib as kernel
    import scipy.optimize._slsqp_py as wrapper
    from mace_site_kinematics import Kinematics
    from accommodation_torsion_profiles import geometry_check
    from affordable_common import paired
    import slsqp_precision
    old = read(plan['template']); prepared = read(plan['target_admission'])
    assert plan['canonical_reference'] == read(plan['preparation'])['canonical_reference']
    verify(plan['canonical_reference']); verify(plan['plan'])
    for key in ('implementation', 'source_modules', 'source_snapshot', 'model', 'settings', 'pool_settings',
                'optimizer_software', 'optimizer_executable', 'model_executable', 'optimizer_python',
                'model_python', 'calculator_source', 'fixtures', 'gates', 'threads_per_model'):
        assert m[key] == old[key], 'unchanged scientific/runtime field differs: '+key
    for key in ('implementation', 'calculator_source', 'optimizer_executable', 'model_executable'): verify(m[key])
    for p in m['source_modules'].values(): verify(p)
    verify(m['model']['checkpoint'])
    assert scipy.__version__ == m['optimizer_software']['version']
    assert pin(kernel.__file__) == m['optimizer_software']['kernel'] and pin(wrapper.__file__) == m['optimizer_software']['wrapper']
    assert slsqp_precision.engine()[0].completion.SETTINGS == m['settings']
    assert len(m['cases']) == 1 and m['origin_task_count'] == 2 and len(m['tasks']) == 4
    assert m['lanes'] == 1 and m['maximum_active_numerical_threads'] == 7
    ts = m['tasks'][:2]; assert [t['metal'] for t in ts] == ['Ca', 'La']
    checks = []
    for t in ts:
        ep = prepared['representations']['context']['endpoints'][t['metal']]
        assert all(t[k] == ep[k] for k in ('xyz', 'charge', 'multiplicity'))
        assert t['mapping'] == prepared['maps'][t['metal']] and t['prepared_admission'] == plan['target_admission']
        for k in ('mapping', 'state_signature', 'context_preparation'): verify(t[k])
        a = b.xyz(verify(t['xyz'])); kin = Kinematics(read(t['mapping'])['context']); q = np.zeros(len(kin.modes))
        residual = float(np.max(np.abs(kin.evaluate(q)[1]-np.asarray([v[1:] for v in a]))))
        assert residual <= 1e-12 and geometry_check(kin, q, [v[0] for v in a])['pass']
        checks.append({'metal': t['metal'], 'q0_max_coordinate_error_A': residual, 'physical_map_pass': True})
    paired(verify(ts[1]['xyz']), verify(ts[0]['xyz']), ts[1]['charge'], ts[0]['charge'])
    assert read(ts[0]['mapping']) == read(ts[1]['mapping'])
    return {'status': 'pass', 'checks': checks, 'changed_method_fields': [], 'new_model_calls': 0}


def execute(path):
    placement = allocation(); plan = read(path); assert plan['implementation'] == pin(__file__)
    mp = verify(plan['worker_manifest']); m = read(mp); base = Path(path).parent
    with (base/'RUN.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        attempt = base/'attempts'/'001'; attempt.mkdir(parents=True, exist_ok=False); (attempt/'lane0').mkdir()
        write(attempt/'PLACEMENT.json', placement)
        spec = importlib.util.spec_from_file_location('frozen_disulfide_cpu', verify(m['implementation']))
        b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b); b.modules(m)
        old_handlers = b.handlers(); transport = None; error = None; results = []; started = time.monotonic()
        try:
            write(attempt/'PREFLIGHT.json', preflight(plan, m, b))
            transport = b.WarmCPU(mp, m, attempt, 0)
            import numpy as np
            checks = []
            for fixture in m['fixtures']:
                t = next(t for t in m['tasks'] if t['task_id'] == fixture['task_id'])
                request = {'manifest': pin(mp), 'task_id': t['task_id'], 'xyz': fixture['xyz'],
                           'charge': t['charge'], 'multiplicity': t['multiplicity']}
                r, rp = b.cell(attempt/'lane0'/'qualification'/fixture['id'], request, transport, m['model'])
                check = {'fixture': fixture['id'], 'receipt': rp, 'pass': False}
                if r['status'] == 'complete':
                    delta = (r['energy_eV']-fixture['energy_eV'])*b.EV_TO_KCAL
                    force = float(np.max(np.abs(np.load(verify(r['forces']), allow_pickle=False)-np.load(verify(fixture['forces']), allow_pickle=False))))
                    check.update(energy_error_kcal_mol=delta, force_max_error_eV_A=force,
                                 **{'pass': abs(delta) <= m['gates']['energy_kcal_mol'] and force <= m['gates']['force_max_eV_A']})
                checks.append(check)
            write(attempt/'QUALIFICATION.json', {'checks': checks, 'pass': all(c['pass'] for c in checks)})
            assert all(c['pass'] for c in checks), 'original CPU force qualification failed'
            for case in m['cases']: results.append(b.run_source(mp, m, case, transport))
        except BaseException as exc:
            error = type(exc).__name__+': '+str(exc); traceback.print_exc()
        finally:
            if transport:
                try:
                    if error: transport.abort()
                    else: transport.close()
                except BaseException as exc:
                    error = (error+'; ' if error else '')+str(exc); transport.abort()
            for s, handler in old_handlers.items(): signal.signal(s, handler)
            calls = [read(p) for p in sorted(base.glob('**/mace_result.json'))]
            cases = [read(p) for p in results]
            elapsed = time.monotonic()-started
            collection = {'protocol_id': PROTOCOL, 'manifest': pin(mp), 'recovery_plan': pin(path),
                'terminal': len(cases) == 1, 'denominator': 1, 'cases': cases, 'error': error,
                'missing_source_ids': [c['source_case_id'] for c in m['cases'] if c['case_id'] not in {r['case_id'] for r in cases}],
                'costs': {'new_recorded_model_calls': len(calls), 'actual_forwards': sum(r['actual_model_forwards'] for r in calls),
                    'recorded_wall_seconds_sum': sum(r['wall_seconds'] for r in calls),
                    'recorded_process_CPU_seconds': sum(r['process_cpu_seconds'] for r in calls),
                    'started_without_terminal_receipt': [str(p) for p in base.glob('**/CPU_CALL_STARTED.json') if not (p.parent/'mace_result.json').exists()],
                    'elapsed_seconds': elapsed, 'allocated_CPU_seconds': elapsed*7, 'placement': placement},
                'canonical_reference': plan['canonical_reference'], 'new_GFN2_calls': 0, 'new_DFT_calls': 0,
                'old_unavailable_result_preserved': True, 'production_changed': False}
            write(base/'FINAL_COLLECTION.json', collection)
        if error: raise RuntimeError(error)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare'); p.add_argument('--template', required=True); p.add_argument('--preparation', required=True); p.add_argument('--output', required=True)
    e = sub.add_parser('execute'); e.add_argument('--manifest', required=True)
    args = parser.parse_args()
    if args.command == 'prepare': prepare(args.template, args.preparation, args.output)
    else: execute(args.manifest)
