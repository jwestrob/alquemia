#!/usr/bin/env python3
"""Run the three unchanged remaining LBT3 searches using the frozen scout code.

CLI: python LBT3_ACCOMMODATE_REMAINING.py --stage remaining --config CONFIG
     --output FRESH [--repo-root REPO]
Config contains only path/SHA256 pins: plan, implementation, scout_result.
Requires an admitted scout and its passed finite-derivative gate. The only scope
changes are SOURCE/METAL and the cache recipe's physical multiplicity before
the origin is written. Search, forces, bounds, contacts and energy stay intact.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import socket
import sys
import time
import traceback

REPO = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
SCOUT_SHA = 'ee1db81bdf3db90c303f824d7a55b52b00c2bb64f2422903aed60c25931d8ca1'
TASKS = (('7CCO', 'Dy'), ('7CCN', 'La'), ('7CCN', 'Dy'))
PROTOCOL = 'lbt3_whole_own_composite_remaining_three_v1'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(record):
    path = Path(record['path']).resolve(strict=True)
    require(path.is_file() and pin(path)['sha256'] == record['sha256'], 'Missing/changed input: ' + record['path'])
    return path


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['remaining'], default='remaining')
    parser.add_argument('--repo-root', type=Path, default=REPO)
    for key in ('config', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    cpus = int(os.environ['SLURM_CPUS_PER_TASK']); require(cpus >= 2, 'Two allocated CPU workers required')
    require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'Allocated GPU required')
    cfg = json.loads(args.config.read_text()); verify(cfg['plan'])
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Wrong continuation adapter')
    scout = read(cfg['scout_result']); scout_config = read(scout['config']); admitted = read(scout['search'])
    require(scout['status'] == 'complete' and scout['stage'] == 'scout'
            and scout['continuation_execution_and_physical_admission']
            and scout['source_id'] == admitted['source_id'] == '7CCO'
            and scout['metal'] == admitted['metal'] == 'La'
            and admitted['status'] == 'finite_proposal_available' and admitted['candidate_admitted']
            and admitted['final_geometry']['pass'], 'Execution/physical scout gate did not pass')
    require(scout['implementation']['sha256'] == SCOUT_SHA
            and scout['implementation'] == scout_config['implementation']
            and scout['plan'] == scout_config['plan'] == cfg['plan'], 'Frozen scout code/plan differs')
    base = module('lbt3_frozen_scout_for_remaining', scout['implementation'])
    require(base.SETTINGS == scout['settings'] and (base.SOURCE, base.METAL) == ('7CCO', 'La'), 'Scout defaults differ')
    profile = read(scout['profile_result']); fd_result = read(scout['fd_result'])
    require(profile['status'] == 'complete' and fd_result['status'] == 'passed' and fd_result['composite_passed_all_states']
            and fd_result['design']['profile_result'] == scout['profile_result']
            and scout_config['profile_result'] == scout['profile_result']
            and scout_config['fd_result'] == scout['fd_result'], 'Qualified profile/derivative lineage differs')
    require(profile['active_candidate_basis']['mode_names_selection_order'] == base.NAMES, 'Common active basis differs')
    verify(fd_result['config']); fd_module = module('lbt3_remaining_native', fd_result['implementation'])
    profile_module = module('lbt3_remaining_profile', profile['implementation'])
    modes = read(profile['modes_result']); solvent = read(profile['solvent_result'])
    mforce = read(profile['mace_forces']); mcfg = read(mforce['config']); scfg = read(solvent['config'])
    require(modes['status'] == solvent['status'] == mforce['status'] == 'complete'
            and modes['preparation'] == solvent['preparation'] == mcfg['preparation']
            and mforce['model'] == mcfg['model'] and mforce['versions'] == mcfg['expected_versions'], 'Origin model/source lineage differs')
    solvent_module = module('lbt3_remaining_solvent', solvent['implementation'])
    require(scfg['tolerances'] == solvent_module.TOLERANCES, 'Solvent admission differs')
    for name, record in profile['helpers'].items():
        verify(record); require(pin(args.repo_root / 'scripts' / name)['sha256'] == record['sha256'], 'Helper changed: ' + name)
    sys.path.insert(0, str(args.repo_root.resolve() / 'scripts'))
    import numpy as np
    from mace_site_kinematics import Kinematics
    from scaffold_accommodation import Geometry, SETTINGS as GEOMETRY_SETTINGS, lbfgs_direction
    base.np = np  # The scout's own main performs exactly this module initialization.
    prep = read(modes['preparation']); mode_rows = {r['source_id']: r for r in modes['rows']}
    prep_rows = {r['source_id']: r for r in prep['rows']}
    corrections = {(r['source_id'], r['metal']): r for r in solvent['corrections']}
    raw = {(r['task']['origin_task_id'], r['task']['medium']): r for r in solvent['rows']}
    force_by = {r['task_id']: r for r in mforce['rows']}
    require(len(raw) == 8 and len(force_by) == len(corrections) == 4
            and all(r['status'] == 'accepted' for r in raw.values()), 'Origin matrix incomplete')
    world = dict(profile=profile, modes=modes, mcfg=mcfg, scfg=scfg, raw=raw, force_by=force_by,
        solvent_module=solvent_module, fd_module=fd_module, manifest=read(solvent['manifest']))
    active = profile['active_candidate_basis']['mode_indices_selection_order']

    class ScopedEvaluator(base.CompositeEvaluator):
        def save_origin(self):
            expected = 1 if base.METAL == 'La' else 6
            require((base.SOURCE, base.METAL) in TASKS and self.base['source_id'] == base.SOURCE
                    and self.base['metal'] == base.METAL and self.base['charge'] == -2
                    and self.base['multiplicity'] == self.base['spin_multiplicity'] == expected,
                    'Prepared physical source/spin scope differs')
            for medium in ('vacuum', 'ALPB(Water)'):
                task = self.raw[(self.base['task_id'], medium)]['task']
                require(task['source_id'] == base.SOURCE and task['metal'] == base.METAL
                        and task['physical_multiplicity'] == expected and task['multiplicity'] == 1
                        and task['uhf'] == 0, 'GFN2 physical/effective state differs')
            # Invoked by the inherited initializer before any origin cache key/file.
            self.recipe['physical_multiplicity'] = expected
            return super().save_origin()

    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    result = {'protocol_id': PROTOCOL, 'underlying_search_protocol': base.PROTOCOL,
        'stage': 'remaining', 'status': 'unavailable', 'implementation': pin(__file__),
        'config': pin(args.config), 'plan': cfg['plan'], 'scout_result': cfg['scout_result'],
        'frozen_scout_implementation': scout['implementation'], 'profile_result': scout['profile_result'],
        'fd_result': scout['fd_result'], 'native_helper': fd_result['implementation'],
        'solvent_helper': solvent['implementation'], 'profile_helper': profile['implementation'],
        'helpers': profile['helpers'], 'settings': base.SETTINGS,
        'declared_tasks': [{'source_id': s, 'metal': m} for s, m in TASKS],
        'scope_patch': {'module_globals': ['SOURCE', 'METAL'],
            'cache_recipe_before_origin': 'physical_multiplicity from exact prepared endpoint: La1/Dy6',
            'GFN2_effective_multiplicity': 1, 'execution': 'strictly sequential states; fresh evaluator, guards, gradient and history per state'},
        'rows': [], 'points': [], 'pool_comparison_performed': False, 'stationary_minimum_claimed': False,
        'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname()}
    write(out / 'START.json', result)
    warm = None
    for source_id, metal in TASKS:
        directory = out / (source_id + '_' + metal); directory.mkdir()
        base.SOURCE, base.METAL = source_id, metal
        prior_calls = warm.calls if warm is not None else 0
        evaluator = None
        row = {'source_id': source_id, 'metal': metal, 'status': 'unavailable', 'candidate_admitted': False,
               'physical_multiplicity': 1 if metal == 'La' else 6, 'GFN2_effective_multiplicity': 1}
        try:
            # Fresh source objects prevent stale source/guard/gradient state across tasks.
            source = profile_module.source_profile(source_id, mode_rows[source_id], prep_rows[source_id], corrections,
                force_by, raw, solvent['gradient_units'], np, Kinematics)
            require([source['kin'].modes[k]['name'] for k in active] == base.NAMES, 'Active index/source map differs')
            evaluator = ScopedEvaluator(world, source, directory)
            evaluator.native = warm
            search = base.search(source, active, evaluator, directory, Geometry, GEOMETRY_SETTINGS, lbfgs_direction)
            write(directory / 'SEARCH.json', search)
            row.update(status='complete' if search['candidate_admitted'] else 'unavailable',
                       candidate_admitted=search['candidate_admitted'], search=pin(directory / 'SEARCH.json'))
        except Exception as exc:
            row.update(failure=str(exc), traceback=traceback.format_exc())
        if evaluator is not None:
            warm = evaluator.native
            row.update(actual_cost=evaluator.cost, points=evaluator.records, requests=evaluator.requests)
            result['points'].extend(evaluator.records)
            if warm is not None:
                row['native_runtime'] = {'versions': warm.versions, 'device': warm.device,
                    'worker_prior_forwards': prior_calls, 'worker_cumulative_forwards': warm.calls,
                    'actual_model_forwards_this_search': warm.calls - prior_calls,
                    'torch_threads': warm.torch.get_num_threads()}
                require(row['native_runtime']['actual_model_forwards_this_search'] == evaluator.cost['MACE_calls_attempted'],
                        'Per-search warm-model accounting differs')
        write(directory / 'RESULT.json', row); row['result'] = pin(directory / 'RESULT.json')
        result['rows'].append(row)
        print(json.dumps({'source_id': source_id, 'metal': metal, 'status': row['status']}), flush=True)
    result['status'] = 'complete' if all(r['candidate_admitted'] for r in result['rows']) else 'unavailable'
    result['actual_cost'] = {key: sum(row.get('actual_cost', {}).get(key, 0) for row in result['rows']) for key in
        ('new_composite_evaluations', 'MACE_calls_attempted', 'GFN2_calls_attempted', 'MACE_call_seconds',
         'GFN2_call_seconds', 'model_load_seconds', 'cache_hits')}
    if warm is not None:
        result['native_runtime'] = {'actual_model_forwards': warm.calls,
            'model_loaded_once': True, 'peak_cuda_allocated_bytes': warm.torch.cuda.max_memory_allocated()}
    require(result['actual_cost']['new_composite_evaluations'] <= 72
            and result['actual_cost']['MACE_calls_attempted'] <= 72 and result['actual_cost']['GFN2_calls_attempted'] <= 144,
            'Remaining finite call ceiling exceeded')
    result['elapsed_seconds'] = time.monotonic() - started
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds'] * cpus
    result['peak_host_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['allocation'] = {k: os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_CPU', 'CUDA_VISIBLE_DEVICES')}
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}), flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
