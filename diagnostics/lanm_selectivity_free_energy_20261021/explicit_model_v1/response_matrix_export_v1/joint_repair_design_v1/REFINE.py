#!/usr/bin/env python3
"""One declared joint development fit; native calls remain inside its evaluator."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import time

import numpy as np

THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
           'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
LIMITS = {'outer_steps': 24, 'new_vectors': 672, 'native_programs': 24000, 'workers': 12}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path']) == record, 'Changed artifact: '+record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def token(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def fingerprint(config):
    return token({k: v for k, v in config.items()
                  if k not in ('root_release', 'preparation_result')})


def module(record, name):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def save(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    os.replace(temporary, path)


def run(stage, config_path, output):
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out/'vectors').mkdir()
    tick, cpu = time.monotonic(), time.process_time()
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    cfg = json.loads(Path(config_path).read_text())
    ev = None
    records, values_by_key = {}, {}
    result = {'status': 'started', 'stage': stage, 'implementation': pin(__file__),
              'configuration': pin(config_path), 'scientific_fingerprint': fingerprint(cfg),
              'limits': LIMITS, 'vectors': records, 'iterations': [], 'error': None,
              'physical_accuracy_qualified': False, 'validation_outcomes_opened': False,
              'automatic_extension': False, 'job_id': os.environ.get('SLURM_JOB_ID')}

    def persist():
        if ev is not None:
            result['evaluator'] = ev.snapshot()
        result['new_parameter_vectors'] = max(0, len(records)-1)
        save(out/'RESULT.json', result)

    def retain(key, value):
        path = out/'vectors'/(key+'.json')
        save(path, value)
        values_by_key[key] = value
        score = value.get('score')
        records[key].update(status=value['status'], artifact=pin(path))
        if score is not None:
            records[key].update(objective=score['objective'],
                                block_objectives=score['block_objectives'],
                                original_EF_checks_passed=score['original_EF_checks_passed'])
        if value.get('error'):
            records[key]['error'] = value['error']

    try:
        require(cfg['implementation'] == pin(__file__), 'Driver implementation differs')
        require(all(os.environ.get(k) == '1' for k in THREADS), 'One numerical thread per worker')
        for name in ('protocol', 'root_declaration', 'parameter_domain', 'directions',
                     'evaluator', 'evaluator_config', 'objective_ready'):
            verify(cfg[name])
        require(cfg['limits'] == LIMITS, 'Declared finite limits differ')
        domain = read(cfg['parameter_domain'])
        parameters = domain['parameters']
        require(len(parameters) == 12, 'Exactly12 box parameters required')
        lower = np.array([p['lower'] for p in parameters], dtype=float)
        upper = np.array([p['upper'] for p in parameters], dtype=float)
        width = upper-lower
        initial = [p['initial'] for p in parameters]
        require(np.all(width > 0), 'Nonempty finite parameter boxes')
        release_pin = cfg.get('root_release')
        if stage == 'run':
            require(os.environ.get('SLURM_JOB_ID', '').isdigit()
                    and int(os.environ.get('SLURM_NTASKS', '1')) == 1
                    and int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) == 12,
                    'One12-CPU Slurm allocation required')
            release = read(release_pin)
            prep = read(cfg['preparation_result'])
            require(release['status'] == 'released_bounded_joint_response_fit'
                    and release['limits'] == LIMITS
                    and release['scientific_fingerprint'] == fingerprint(cfg),
                    'Root release differs')
            require(prep['status'] == 'prepared_joint_response_fit'
                    and prep['scientific_fingerprint'] == fingerprint(cfg)
                    and prep['evaluator']['native_program_invocations'] == 0,
                    'Saved preparation unavailable or different')
        else:
            require(release_pin is None and cfg.get('preparation_result') is None,
                    'Preparation must be unbound to execution')
        evaluator = module(cfg['evaluator'], 'joint_response_evaluator')
        ev = evaluator.Evaluator(read(cfg['evaluator_config']), out/'cache', root_release=release_pin)
        require(ev.initial_vector == initial and ev.family['parameters'] == parameters,
                'Evaluator/domain starting vector differs')
        math_helper = module(cfg['directions'], 'joint_response_directions')
        first = ev.evaluate(initial)
        require(ev.snapshot()['native_program_invocations'] == 0, 'Initial replay made native calls')
        require(len(first['score']['weighted_residuals']) == 531, 'Declared objective dimension')
        current_key = token(initial)
        records[current_key] = {'status': 'started', 'role': 'initial', 'parameters': initial}
        retain(current_key, {'status': 'complete', 'parameters': initial, **first})
        result.update(initial_key=current_key, current_key=current_key, parameter_domain=cfg['parameter_domain'])
        persist()
        if stage == 'prepare':
            result['status'] = 'prepared_joint_response_fit'
        else:
            require(prep['vectors'][prep['initial_key']]['parameters'] == initial,
                    'Prepared initial vector differs')

            def batch(points, role, iteration):
                keys = [token(p) for p in points]
                fresh = list(dict.fromkeys(k for k in keys if k not in records))
                require(len(records)+len(fresh) <= 673, 'Vector ceiling')
                lookup = dict(zip(keys, points))
                for key in fresh:
                    records[key] = {'status': 'scheduled', 'role': role,
                                    'outer_step': iteration, 'parameters': lookup[key]}
                persist()
                failed = []
                with ThreadPoolExecutor(max_workers=12) as pool:
                    futures = {pool.submit(ev.evaluate, lookup[k]): k for k in fresh}
                    for future in as_completed(futures):
                        key = futures[future]
                        try:
                            answer = future.result()
                            require(len(answer['score']['weighted_residuals']) == 531
                                    and np.isfinite(answer['score']['objective']), 'Finite531-vector result')
                            retain(key, {'status': 'complete', 'parameters': lookup[key], **answer})
                        except Exception as exc:
                            ev.stop_event.set()
                            failed.append(key)
                            retain(key, {'status': 'unavailable', 'parameters': lookup[key], 'error': repr(exc)})
                            for other in futures:
                                if other is not future:
                                    other.cancel()
                        persist()
                require(not failed, 'Native/vector failure; no subsequent batch')
                require(all(values_by_key[k]['status'] == 'complete' for k in keys),
                        'Unavailable cached vector; no automatic retry')
                return keys

            stop = 'outer_step_limit'
            for iteration in range(1, 25):
                current = values_by_key[current_key]
                center = np.asarray(current['parameters'], dtype=float)
                z = (center-lower)/width
                entry = {'outer_step': iteration, 'origin_key': current_key, 'status': 'derivatives'}
                result['iterations'].append(entry)
                points = []
                for j in range(12):
                    for sign in (-1, 1):
                        target = z.copy()
                        target[j] = np.clip(target[j]+sign*.005, 0, 1)
                        value = center.copy()
                        value[j] = np.clip(center[j]+(target[j]-z[j])*width[j], lower[j], upper[j])
                        points.append(value.tolist())
                keys = batch(points, 'derivative', iteration)
                columns, spans = [], []
                for j in range(12):
                    a, b = [values_by_key[k] for k in keys[2*j:2*j+2]]
                    span = (b['parameters'][j]-a['parameters'][j])/width[j]
                    require(span > 0, 'Collapsed derivative interval')
                    spans.append(span)
                    columns.append((np.asarray(b['score']['weighted_residuals'])
                                    -np.asarray(a['score']['weighted_residuals']))/span)
                J = np.column_stack(columns)
                residual = np.asarray(current['score']['weighted_residuals'])
                scaled, bounded, diagnostic = math_helper.directions(J, residual, z)
                entry.update(status='trials', derivative_keys=keys, derivative_spans=spans,
                             Jacobian_normalized=J.tolist(), linear_directions=diagnostic)
                trials = []
                for name, direction in (('scaled_GN', scaled), ('bounded_LS', bounded)):
                    for fraction in (.5, 1.):
                        value = np.clip(center+fraction*direction*width, lower, upper)
                        delta = (value-center)/width
                        require(np.max(abs(delta)) <= .025+1e-12, 'Trial exceeded trust box')
                        trials.append({'method': name, 'fraction': fraction, 'parameters': value.tolist(),
                                       'predicted_objective': float(np.linalg.norm(residual+J@delta)**2)})
                entry['trial_points'] = trials
                trial_keys = batch([p['parameters'] for p in trials], 'trial', iteration)
                selected_key = current_key
                for key in trial_keys:
                    if values_by_key[key]['score']['objective'] < values_by_key[selected_key]['score']['objective']:
                        selected_key = key
                improved = selected_key != current_key
                current_key = selected_key
                entry.update(status='complete', trial_keys=trial_keys,
                             selected_key=current_key, actual_improvement=improved)
                result['current_key'] = current_key
                persist()
                if not improved:
                    stop = 'no_actual_improvement_at_complete_trial_batch'
                    break
            selected = values_by_key[current_key]
            result.update(status='bounded_joint_response_fit_complete', stopping_reason=stop,
                          selected_key=current_key, selected=records[current_key]['artifact'],
                          selected_objective=selected['score']['objective'],
                          selected_blocks=selected['score']['block_objectives'],
                          original_EF_checks_passed=selected['score']['original_EF_checks_passed'],
                          completed_outer_steps=len(result['iterations']), candidate_frozen=True,
                          last_Jacobian_origin_key=result['iterations'][-1]['origin_key'],
                          bound_contacts=[{'name': p['name'], 'lower': x == p['lower'], 'upper': x == p['upper']}
                                          for p, x in zip(parameters, selected['parameters'])])
    except Exception as exc:
        result.update(status='failed_or_unavailable', error=repr(exc), candidate_frozen=False)
        if ev is not None:
            ev.stop_event.set()
    finally:
        end = resource.getrusage(resource.RUSAGE_CHILDREN)
        result['cost'] = {'wall_seconds': time.monotonic()-tick,
                          'parent_CPU_seconds': time.process_time()-cpu,
                          'child_CPU_seconds': end.ru_utime+end.ru_stime-children.ru_utime-children.ru_stime}
        persist()
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] in ('prepared_joint_response_fit', 'bounded_joint_response_fit_complete') else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('prepare', 'run'))
    parser.add_argument('--config', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.stage, args.config, args.out))
