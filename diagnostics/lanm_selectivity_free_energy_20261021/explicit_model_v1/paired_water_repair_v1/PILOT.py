#!/usr/bin/env python3
"""Saved preparation or one root-released four-parameter paired-water pilot."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import time

THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path']) == p, 'Changed pin: '+p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def fingerprint(cfg):
    payload = {k: v for k, v in cfg.items() if k not in ('root_release', 'preparation_result')}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def save(path, obj):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')
    os.replace(temporary, path)


def module(p, name):
    s = importlib.util.spec_from_file_location(name, verify(p))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def prepare(cfg, out):
    require(cfg['implementation'] == pin(__file__), 'Driver differs')
    require(all(os.environ.get(k) == '1' for k in THREADS), 'Numerical threads must be one per worker')
    verify(cfg['specification'])
    verify(cfg['search'])
    evaluator_cfg = read(cfg['evaluator_config'])
    require(evaluator_cfg['implementation'] == cfg['evaluator'] and evaluator_cfg['specification'] == cfg['specification'],
            'Evaluator source/specification binding')
    require(evaluator_cfg.get('root_release') is None, 'Evaluator config must remain release-free')
    native = module(cfg['evaluator'], 'paired_water_native_evaluator')
    ev = native.Evaluator(evaluator_cfg, out/'cache', root_release=cfg.get('root_release'))
    obj = module(cfg['objective'], 'paired_water_objective')
    helper = module(cfg['body_helper'], 'paired_water_body_helper')
    require(cfg['body_helper']['sha256'] == 'bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819',
            'Body weighting helper changed')
    require(cfg['water48']['sha256'] == '7be6763c6eddaa1501c0e895851433ab7117b6e039472df5b0bf50561a90e818',
            'Water48 source changed')
    water48 = read(cfg['water48'])
    panel = read(water48['panel'])
    require(water48['status'] == 'complete', 'Water48 unavailable')
    model = ev.evaluate(ev.family['initial_vector'])
    score = obj.evaluate(ev.data, model, helper)
    sensitivity = obj.method_sensitivity(ev.data, model, water48, panel)
    require(ev.snapshot()['native_program_invocations'] == 0, 'Preparation made native calls')
    return ev, obj, helper, water48, panel, model, score, sensitivity


def main(config_path, out, stage):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    start, cpu = time.monotonic(), time.process_time()
    children_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    cfg = json.loads(Path(config_path).read_text())
    result = {'status': 'unavailable', 'stage': stage, 'implementation': pin(__file__), 'config': pin(config_path),
        'scientific_config_sha256': fingerprint(cfg), 'new_native_program_invocations': 0,
        'cells': {k: {'status': 'not_attempted'} for k in ('water_central_La', 'water_central_Dy', 'aquo8_central_La', 'aquo8_central_Dy')},
        'candidate_frozen': False, 'validation_outcomes_opened': False, 'physical_accuracy_qualified': False,
        'production_released': False, 'job_id': os.environ.get('SLURM_JOB_ID'), 'error': None}
    ev = None
    try:
        if stage == 'run':
            require(os.environ.get('SLURM_JOB_ID', '').isdigit() and int(os.environ.get('SLURM_NTASKS', '1')) == 1 and
                    int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) == 8, 'Run requires one eight-CPU Slurm task')
            preparation = read(cfg['preparation_result'])
            require(preparation['status'] == 'prepared' and preparation['stage'] == 'prepare' and
                    preparation['new_native_program_invocations'] == 0 and
                    preparation['scientific_config_sha256'] == fingerprint(cfg), 'Actual zero-call preparation differs')
            release = read(cfg['root_release'])
            require(release['status'] == 'approved_paired_water_representability', 'Wrong root release')
            for k in ('implementation', 'specification', 'objective', 'search', 'evaluator', 'evaluator_config', 'preparation_result'):
                require(release[k] == cfg[k], 'Root release differs: '+k)
            require(release['driver_config_sha256'] == fingerprint(cfg) and release['workers'] == 8 and
                    release['max_parameter_vectors'] == 241 and release['max_native_program_invocations'] == 1920,
                    'Root execution limits differ')
            require(release['validation_scoring_authorized'] is False and release['physical_accuracy_claimed'] is False,
                    'Root release exceeds development scope')
        else:
            require(cfg.get('root_release') is None and cfg.get('preparation_result') is None,
                    'Preparation must use the unreleased config')
        ev, obj, helper, water48, panel, baseline, score, sensitivity = prepare(cfg, out)
        result.update(family=ev.family, baseline_score=score, baseline_method_sensitivity=sensitivity,
                      baseline_model=baseline, evaluator=ev.snapshot())
        result['cells'] = {k: {'status': 'baseline_admitted', 'provenance': v['provenance']} for k, v in baseline.items()}
        if stage == 'prepare':
            result['status'] = 'prepared'
        else:
            search = module(cfg['search'], 'paired_water_search')
            require(search.VECTOR_CAP == 241 and search.SEARCH_CAP == 233, 'Search caps changed')
            def evaluate(values):
                model = ev.evaluate(values)
                return {'model': model, 'score': obj.evaluate(ev.data, model, helper)}
            def update(state):
                result.update(search=state, evaluator=ev.snapshot(),
                              new_native_program_invocations=ev.snapshot()['native_program_invocations'])
                save(out/'RESULT.json', result)
            state = search.search(ev.family['parameters'], evaluate, ev.stop_event, update, workers=8)
            best = state['evaluations'][state['best_key']]
            result.update(status='finite_pilot_complete', candidate_parameters=state['candidate_parameters'],
                candidate_frozen=True, selected_model=best['model'], selected_score=best['score'],
                selected_method_sensitivity=obj.method_sensitivity(ev.data, best['model'], water48, panel),
                numerical_optimization_complete=state['numerical_optimization_complete'],
                optimizer_budget_limited=state['optimizer_budget_limited'],
                development_checks_passed=state['development_checks_passed'],
                identifiability_and_bounds_passed=state['identifiability_and_bounds_passed'],
                automatic_validation_release=False)
            result['candidate_keys'] = {}
            for metal in ('La', 'Dy'):
                key, rendered = ev.render(state['candidate_parameters'], 'water_central_'+metal)
                path = out/('CANDIDATE_'+metal+'.key')
                path.write_text(key)
                result['candidate_keys'][metal] = {'key': pin(path), 'rendered': rendered,
                    'scope': 'small-molecule SI water only; not a periodic protein key'}
    except Exception as exc:
        result.update(status='failed_or_unavailable', error=repr(exc), candidate_frozen=False)
    finally:
        if ev is not None:
            result['evaluator'] = ev.snapshot()
            result['new_native_program_invocations'] = result['evaluator']['native_program_invocations']
        children_end = resource.getrusage(resource.RUSAGE_CHILDREN)
        result.update(elapsed_seconds=time.monotonic()-start, parent_process_CPU_seconds=time.process_time()-cpu,
            aggregate_native_child_user_CPU_seconds=children_end.ru_utime-children_start.ru_utime,
            aggregate_native_child_system_CPU_seconds=children_end.ru_stime-children_start.ru_stime,
            child_CPU_accounting='One start/end RUSAGE_CHILDREN difference; do not sum overlapping worker receipts')
        save(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'new_native_program_invocations': result['new_native_program_invocations'],
                      'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] in ('prepared', 'finite_pilot_complete') else 1


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('stage', choices=('prepare', 'run'))
    p.add_argument('--config', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    raise SystemExit(main(a.config, a.output, a.stage))
