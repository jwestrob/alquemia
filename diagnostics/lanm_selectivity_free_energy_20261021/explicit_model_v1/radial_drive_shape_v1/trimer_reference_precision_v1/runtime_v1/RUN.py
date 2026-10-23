#!/usr/bin/env python3
"""Serial, root-released login execution of the fixed fourteen-state protocol."""
import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import re
import shutil
import signal
import socket
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import RESIDUALS as residuals
require = residuals.require
THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
           'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')
OLD_SCF = '%scf\n Guess PModel\n MaxIter 300\nend\n'
NEW_SCF = ('%scf\n Guess PModel\n MaxIter 300\n ConvCheckMode 0\n ConvForced 1\n'
           ' TolE 1e-11\n TolMAXP 1e-8\n TolRMSP 1e-9\n TolErr 1e-8\n'
           ' TolG 1e-8\n TolX 1e-8\nend\n')


def pin(path):
    p = Path(path).resolve()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(4*1024*1024), b''):
            h.update(chunk)
    return {'path': str(p), 'sha256': h.hexdigest()}


def verify(p):
    require(pin(p['path']) == p, 'Changed pin: '+p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def save(path, value):
    p = Path(path)
    tmp = p.with_name(p.name+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    os.replace(tmp, p)


def module(p, name):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def deck(old):
    require(old.count(OLD_SCF) == 1 and old.count('%scf') == 1, 'Original SCF block differs')
    new = old.replace(OLD_SCF, NEW_SCF, 1)
    require(new.replace(NEW_SCF, OLD_SCF, 1) == old, 'Deck inverse differs')
    return new


def context(cfg):
    require(cfg['status'] == 'bound_reference_precision_login_configuration' and
            cfg['implementation'] == pin(__file__), 'Configuration/source differs')
    m = read(cfg['manifest'])
    require(m['status'] == 'prepared_reference_precision_login_inputs' and
            m['implementation'] == pin(__file__) and m['residual_parser'] == pin(HERE/'RESIDUALS.py'),
            'Prepared runtime source differs')
    d = read(m['declaration'])
    require(d['status'] == 'frozen_termination_only_protocol_source_implementation_authorized' and
            d['maximum_real_reference_attempts_in_this_declaration'] == 0, 'Source declaration differs')
    pr = read(m['protocol_review'])
    require(pr['status'] == 'admitted_reference_precision_protocol_design' and pr['declaration'] == m['declaration'],
            'Protocol admission differs')
    read(pr['independent_review'])
    inv = read(m['inventory'])
    require(d['inventory'] == m['inventory'] and len(inv['rows']) == len(m['rows']) == 14 and
            m['maximum_attempts'] == 14 and m['orca'] == inv['orca'] == m['methods']['orca'], 'Finite input denominator differs')
    require(cfg['phase'] in ('scout', 'remaining') and cfg['campaign_root'] == m['campaign_root'], 'Phase/root differs')
    require(Path(cfg['output_root']).resolve() == Path(m['campaign_root'])/cfg['phase'], 'Phase output differs')
    for row, old in zip(m['rows'], inv['rows']):
        require(row['old'] == old and row['cell_id'] == old['cell_id'] and
                verify(row['input']).read_text() == deck(verify(old['input']).read_text()), 'Exact prepared deck/lineage differs')
    for p in m['methods'].values():
        verify(p)
    require(m['methods']['mpi'] == pin(shutil.which('mpirun')), 'Resolved MPI launcher changed')
    common = module(m['methods']['common_parser'], 'precision_common_parser')
    for name, p in common.HELPERS.items():
        require(m['methods'][name] == pin(p), 'Inherited native helper changed: '+name)
    collector = module(m['methods']['moment_collector'], 'precision_moment_collector')
    qhelper = module(m['methods']['property_helper'], 'precision_property_helper')
    return m, inv, d, common, collector, qhelper


def release(cfg, cfg_pin, release_pin, m, d):
    r = read(release_pin)
    ids = d['staging']['scout_order' if cfg['phase'] == 'scout' else 'continuation_order']
    require(r['status'] == 'approved_reference_precision_login_phase' and
            r['configuration'] == cfg_pin and r['implementation'] == pin(__file__) and
            r['manifest'] == cfg['manifest'] and r['declaration'] == m['declaration'] and
            r['protocol_review'] == m['protocol_review'] and r['phase'] == cfg['phase'] and
            r['cell_ids'] == ids and r['maximum_phase_attempts'] == len(ids) and
            r['maximum_total_attempts'] == 14 and r['campaign_root'] == cfg['campaign_root'],
            'Execution release differs')
    review = read(r['source_review'])
    require(review['status'] == 'passed_reference_precision_runtime_source_review' and
            review['implementation'] == pin(__file__) and review['manifest'] == cfg['manifest'] and
            review['residual_parser'] == m['residual_parser'], 'Runtime source review differs')
    require(r['execution_host'] == 'biotite' and r['scheduler_job_id'] is None and
            r['reserved_login_numerical_threads'] == 8 and r['total_owned_login_limit'] == 16 and
            r['login_reservation_confirmed'] is True and
            0 <= r['other_owned_login_numerical_threads_at_release'] <= 8, 'Login reservation differs')
    prior = None
    if cfg['phase'] == 'remaining':
        prior = read(cfg['scout_result'])
        admission = read(cfg['scout_root_review'])
        require(prior['status'] == 'reference_precision_scout_complete' and prior['manifest'] == cfg['manifest'] and
                prior['implementation'] == pin(__file__) and prior['native_attempts_this_phase'] == 2 and
                admission['status'] == 'admitted_reference_precision_scout' and admission['result'] == cfg['scout_result'] and
                admission['manifest'] == cfg['manifest'] and r['scout_result'] == cfg['scout_result'] and
                r['scout_root_review'] == cfg['scout_root_review'], 'Fresh admitted two-state scout required')
    else:
        require(cfg['scout_result'] is cfg['scout_root_review'] is None and
                r['scout_result'] is r['scout_root_review'] is None, 'Scout must not import results')
    return r, ids, prior


def process_table():
    rows = {}
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():
            continue
        try:
            if p.stat().st_uid != os.getuid():
                continue
            stat = (p/'stat').read_text()
            comm, tail = stat.rsplit(')', 1)
            v = tail.split()
            rows[int(p.name)] = {'PID': int(p.name), 'PPID': int(v[1]), 'name': comm.split('(', 1)[1],
                'state': v[0], 'user_ticks': int(v[11]), 'system_ticks': int(v[12]),
                'OS_threads': int(v[17]), 'RSS_pages': int(v[21])}
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            pass
    return rows


def runtime(m):
    require(socket.gethostname().split('.')[0] == 'biotite' and not os.environ.get('SLURM_JOB_ID'),
            'Truthful biotite login execution required')
    require(all(os.environ.get(k) == '1' for k in THREADS), 'Thread environment must already be1')
    require(os.environ.get('LD_LIBRARY_PATH') == m['runtime_environment']['LD_LIBRARY_PATH'] and
            os.environ.get('LD_PRELOAD') == m['runtime_environment']['LD_PRELOAD'], 'Loader environment differs')
    values = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    available_kib = int(values['MemAvailable'].split()[0])
    require(available_kib >= 16000*1024, 'MemAvailable below16000MiB')
    return {'host': socket.gethostname(), 'PID': os.getpid(), 'scheduler_job_id': None,
            'MemAvailable_KiB': available_kib, 'MaxCore_MiB_per_rank': 3000, 'MPI_ranks': 4,
            'threads': {k: os.environ[k] for k in THREADS}, 'reserved_login_numerical_threads': 8,
            'owned_process_snapshot': list(process_table().values()),
            'note': 'OS process/thread observations include launchers and idle MPI helpers; they are not a numerical-thread count.'}


def invoke(argv, folder, env):
    tick = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    samples, seen = [], {}
    code = None
    proc = None
    try:
        with (folder/'endpoint.out').open('x') as stdout, (folder/'endpoint.err').open('x') as stderr:
            proc = subprocess.Popen(argv, cwd=folder, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
            save(folder/'PROCESS_START.json', {'argv': argv, 'PID': proc.pid, 'controller_PID': os.getpid(),
                 'host': socket.gethostname(), 'scheduler_job_id': None})
            while True:
                table = process_table()
                owned = {proc.pid}
                while True:
                    more = {pid for pid, row in table.items() if row['PPID'] in owned}
                    if more <= owned:
                        break
                    owned |= more
                rows = [table[p] for p in sorted(owned) if p in table]
                for row in rows:
                    seen[row['PID']] = row
                samples.append({'elapsed_seconds': time.monotonic()-tick, 'processes': rows})
                code = proc.poll()
                if code is not None:
                    break
                time.sleep(1)
    except BaseException:
        if proc is not None and proc.poll() is None:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            code = proc.returncode
        raise
    finally:
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt = {'argv': argv, 'cwd': str(folder), 'returncode': code, 'PID': None if proc is None else proc.pid,
            'scheduler_job_id': None, 'host': socket.gethostname(),
            'threads': {k: env[k] for k in THREADS}, 'PATH': env['PATH'],
            'output': pin(folder/'endpoint.out') if (folder/'endpoint.out').exists() else None,
            'stderr': pin(folder/'endpoint.err') if (folder/'endpoint.err').exists() else None,
            'wall_seconds': time.monotonic()-tick, 'child_user_CPU_seconds': after.ru_utime-before.ru_utime,
            'child_system_CPU_seconds': after.ru_stime-before.ru_stime,
            'child_max_RSS_KiB_cumulative_controller_highwater': after.ru_maxrss,
            'process_observations': {'sampling_seconds': 1, 'sample_count': len(samples),
                'unique_observed_processes': list(seen.values()),
                'peak_observed_OS_threads': max((sum(r['OS_threads'] for r in s['processes']) for s in samples), default=0),
                'peak_observed_RSS_pages_sum': max((sum(r['RSS_pages'] for r in s['processes']) for s in samples), default=0),
                'page_size_bytes': os.sysconf('SC_PAGE_SIZE'), 'clock_ticks_per_second': os.sysconf('SC_CLK_TCK'),
                'limitation': 'Sampled lower coverage; very short helpers may be missed. OS threads are not all numerical workers.'}}
        save(folder/'PROCESS_OBSERVATIONS.json', samples)
        receipt['process_samples'] = pin(folder/'PROCESS_OBSERVATIONS.json')
        save(folder/'NATIVE_RECEIPT.json', receipt)
    return receipt


def admit(row, folder, receipt, common, collector, qhelper):
    text = (folder/'endpoint.out').read_text()
    cell = row['old']['cell']
    parsed = common.helper.parse(text, cell, receipt['returncode'])
    counts, aux_ok = common.auxiliary_admission(text, 990 if cell['metal'] == 'La' else 1008)
    parsed['native_aux_J_basis_function_counts'] = counts
    parsed['checks']['full_parent_aux_J_basis'] = aux_ok
    native_text = residuals.native(text)
    final_grid = re.findall(r'Recomputing exchange energy using gridx3\s*\.{3,}\s*done', native_text)
    final_integration = re.findall(r'Total energy after final integration\s*:\s*('+residuals.N+r')', native_text)
    total = re.findall(r'^Total Energy\s*:\s*('+residuals.N+r')\s+Eh', native_text, re.M)
    correction = parsed['final_COSX_correction_hartree']
    parsed['checks']['final_gridx3_completed'] = len(final_grid) == len(correction) == len(final_integration) == len(total) == 1 and all(math.isfinite(v) for v in correction)
    parsed['checks']['final_integration_energy_stage'] = bool(len(final_integration) == len(total) == 1 and residuals.print_join(final_integration[0], total[0]))
    parsed['status'] = 'complete' if all(parsed['checks'].values()) else 'failed_admission'
    save(folder/'NATIVE_READBACK.json', parsed)
    require(parsed['status'] == 'complete' and all(parsed['checks'].values()), 'Inherited native identity admission failed')
    require(residuals.settings(text) == row['baseline_grid_integral_settings'], 'Actual grid/integral settings changed')
    require(residuals.occupations(text) == row['baseline_occupations'], 'Restricted occupation inventory changed')
    current = residuals.current(text)
    # Save current fields before a strict residual rejection; failed evidence stays inspectable.
    save(folder/'CURRENT_RESIDUALS.json', current)
    require(current['all_required_passed'], 'Declared current SCF residual gate failed')
    moment = collector.quantum_moment(pin(folder/'endpoint.out'), cell['charge'], qhelper)
    require(moment['labels'] == cell['symbols'] and moment['ghost_mask'] == cell['ghost_mask'], 'Property state identity differs')
    return parsed, current, moment


def run(config_path, release_path):
    cfg_pin, release_pin = pin(config_path), pin(release_path)
    cfg = read(cfg_pin)
    out = Path(cfg['output_root']).resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick, cpu = time.monotonic(), time.process_time()
    result = {'status': 'started', 'implementation': pin(__file__), 'configuration': cfg_pin,
              'release': release_pin, 'manifest': cfg['manifest'], 'phase': cfg['phase'],
              'native_attempts_this_phase': 0, 'cells': {}, 'error': None, 'scheduler_job_id': None}
    save(out/'RESULT.json', result)
    lock, locked = None, False
    try:
        m, inv, d, common, collector, qhelper = context(cfg)
        rel, ids, prior = release(cfg, cfg_pin, release_pin, m, d)
        campaign = Path(m['campaign_root'])
        lock = (campaign/'.run.lock').open('a')
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        locked = True
        require(not (campaign/'STOP.json').exists(), 'Campaign stopped after prior failure')
        cells = campaign/'cells'
        cells.mkdir(exist_ok=True)
        existing = {p.name for p in cells.iterdir()}
        expected_existing = set() if prior is None else set(d['staging']['scout_order'])
        require(existing == expected_existing, 'Attempt inventory differs; no duplicate/recovery allowed')
        all_ids = [r['cell_id'] for r in m['rows']]
        result['cells'] = {k: {'cell_id': k, 'status': 'not_attempted'} for k in all_ids}
        if prior:
            for cid in expected_existing:
                old = prior['cells'][cid]
                require(old['status'] == 'complete' and old['result'] == pin(cells/cid/'RESULT.json') and
                        read(old['result'])['status'] == 'complete', 'Scout row unavailable')
                result['cells'][cid] = old
        for cid in ids:
            row = next(x for x in m['rows'] if x['cell_id'] == cid)
            rt = runtime(m)
            folder = cells/cid
            folder.mkdir(exist_ok=False)
            shutil.copyfile(verify(row['input']), folder/'endpoint.inp')
            require(pin(folder/'endpoint.inp')['sha256'] == row['input']['sha256'], 'Native deck copy differs')
            env = os.environ.copy()
            env['PATH'] = str(Path(m['orca']['path']).parent)+os.pathsep+env['PATH']
            require(pin(shutil.which('mpirun', path=env['PATH'])) == m['methods']['mpi'], 'Child MPI resolution differs')
            argv = [m['orca']['path'], 'endpoint.inp', '--bind-to none']
            cell_result = {'status': 'started', 'cell_id': cid, 'implementation': pin(__file__), 'manifest': cfg['manifest'],
                'configuration': cfg_pin, 'release': release_pin, 'cell': row['old']['cell'], 'old_reference': row['old']['result'],
                'input': pin(folder/'endpoint.inp'), 'runtime': rt, 'native_invocation_attempts': 1, 'error': None}
            # Intent is durable before Popen; a failed launch still consumes this state's sole attempt.
            save(folder/'START.json', cell_result)
            result['native_attempts_this_phase'] += 1
            result['cells'][cid] = {'cell_id': cid, 'status': 'started', 'start': pin(folder/'START.json')}
            save(out/'RESULT.json', result)
            try:
                receipt = invoke(argv, folder, env)
                cell_result['native_receipt'] = pin(folder/'NATIVE_RECEIPT.json')
                cell_result['parsed'], cell_result['current_residuals'], cell_result['moment'] = admit(row, folder, receipt, common, collector, qhelper)
                cell_result['status'] = 'complete'
            except BaseException as exc:
                cell_result.update(status='failed', error=repr(exc))
            finally:
                for name, filename in [('output', 'endpoint.out'), ('stderr', 'endpoint.err'),
                                       ('native_receipt', 'NATIVE_RECEIPT.json'), ('native_readback', 'NATIVE_READBACK.json'),
                                       ('property_file', 'endpoint.property.txt'), ('current_residuals_file', 'CURRENT_RESIDUALS.json')]:
                    if (folder/filename).exists():
                        cell_result[name] = pin(folder/filename)
                save(folder/'RESULT.json', cell_result)
                result['cells'][cid] = {'cell_id': cid, 'status': cell_result['status'], 'result': pin(folder/'RESULT.json')}
                save(out/'RESULT.json', result)
            require(cell_result['status'] == 'complete', cid+': '+str(cell_result['error']))
        result['status'] = 'reference_precision_scout_complete' if cfg['phase'] == 'scout' else 'reference_precision_all14_complete'
        result['native_attempts_total'] = len(list(cells.iterdir()))
        require(result['native_attempts_total'] == (2 if cfg['phase'] == 'scout' else 14), 'Terminal attempt denominator differs')
        result['full_cycle_comparison_pending'] = True
        result['scope'] = 'Admitted endpoint results; complete seven-term comparison is a later saved-only operation after all14.'
    except BaseException as exc:
        result.update(status='reference_precision_phase_failed', error=repr(exc))
        if locked:
            stop = Path(cfg['campaign_root'])/'STOP.json'
            if not stop.exists():
                save(stop, {'configuration': cfg_pin, 'phase': cfg['phase'], 'error': result['error']})
    finally:
        result['cost'] = {'wall_seconds': time.monotonic()-tick, 'controller_CPU_seconds': time.process_time()-cpu,
                          'native_costs': 'Each immutable NATIVE_RECEIPT; no allocated-Slurm cost exists.',
                          'host': socket.gethostname(), 'PID': os.getpid(), 'reserved_aggregate_login_threads': 8}
        save(out/'RESULT.json', result)
        if lock is not None:
            lock.close()
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--release', required=True)
    args = ap.parse_args()
    result = run(args.config, args.release)
    print(json.dumps({'status': result['status'], 'error': result['error'], 'attempts': result['native_attempts_this_phase']}))
    raise SystemExit(0 if result['error'] is None else 1)
