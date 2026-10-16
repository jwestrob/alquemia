#!/usr/bin/env python3
"""One frozen Dy continuation after the admitted La scout; Slurm only.

CLI --stage remaining|mm --config CONFIG --output FRESH.
Remaining config: implementation, plan, pilot_result, search_pattern,
frozen_scout_implementation, scout_result. Root alone submits.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import socket
import time
import traceback

BASE = {
    'path': '/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_coupled_response_20261016/MEX_CONTINUATION.py',
    'sha256': 'c571f20ec76bc7bf7ecadb7b3f2c402caf0f1bdfcb57b6a84aaf4c3a2dd15712',
}


def load_base():
    path = Path(BASE['path'])
    if hashlib.sha256(path.read_bytes()).hexdigest() != BASE['sha256']:
        raise ValueError('Frozen continuation implementation differs')
    spec = importlib.util.spec_from_file_location('mex_remaining_frozen_base', path)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    return base


def scout_gate(base, cfg):
    require = base.require
    base.verify(cfg['frozen_scout_implementation'])
    require(cfg['frozen_scout_implementation']['sha256'] == BASE['sha256'], 'Configured scout implementation differs')
    scout = base.read(cfg['scout_result']); scfg = base.read(scout['config'])
    require(scout['stage'] == 'scout' and scout['status'] == 'complete', 'La scout unavailable')
    for value in (scout['implementation'], scfg['implementation']):
        base.verify(value)
        require(value['sha256'] == BASE['sha256'], 'Executed scout implementation differs')
    for key in ('plan', 'pilot_result', 'search_pattern'):
        base.verify(cfg[key]); base.verify(scfg[key])
        require(scfg[key] == cfg[key], 'Scout shared input differs: ' + key)
    require(scout['plan'] == cfg['plan'] and scout['pilot_result'] == cfg['pilot_result']
            and scout['settings'] == base.SETTINGS, 'Scout declaration differs')
    pre = base.read(scout['preflight'])
    require(pre['status'] == 'passed' and pre['new_molecular_calls'] == 0
            and pre['pilot_result'] == cfg['pilot_result'] and pre['modes'] == scout['modes'], 'Scout preflight differs')
    require(len(pre['starts']) == 4 and {(r['proposer'], r['metal']) for r in pre['starts']}
            == {(p, m) for p in ('La', 'Dy') for m in ('La', 'Dy')}
            and all(r['passed'] and r['geometry']['passed'] for r in pre['starts']), 'Scout start denominator/admission differs')
    require(len(scout['rows']) == 1 and scout['rows'][0]['metal'] == 'La', 'Scout search denominator differs')
    row = scout['rows'][0]; search = base.read(row['search'])
    require(row['status'] == search['status'] == 'finite_proposal_available'
            and row['candidate_admitted'] is True and search['candidate_admitted'] is True
            and search['metal'] == 'La' and search['settings'] == base.SETTINGS
            and search['pilot_result'] == cfg['pilot_result'] and search['initial_history_length'] == 0
            and search['final_geometry']['passed'] and search['accepted_steps'] <= 12
            and search['new_own_composite_evaluations'] <= 24
            and search['work_from_continuation_start']['composite_kcal_mol'] <= 1e-8,
            'La search admission differs')
    require(search['points'] == row['points'] == scout['points']
            and search['requests'] == scout['requests'], 'Scout point/request receipts differ')
    point = base.read(search['candidate']); cell = base.read(point['cell'])
    require(point['status'] == cell['status'] == 'complete' and point['metal'] == cell['metal'] == 'La'
            and point['native_energies'] == cell['native_energies']
            and point['assembled_gradient'] == cell['assembled_gradient'] and point['trial'] == cell['trial'],
            'Scout terminal component binding differs')
    for value in (point['physical'], point['assembled_gradient'], point['trial'],
                  *(cell[k] for k in ('MACE', 'GFN2_embedded', 'GFN2_isolated', 'ALPB', 'MM_result'))):
        base.verify(value)
    return scout, pre


def run(base, cfg, out, result):
    base.require(int(os.environ['SLURM_CPUS_PER_TASK']) >= 2, 'Two GFN workers require allocated CPUs')
    base.require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', 'NoDevFiles'), 'Allocated GPU required')
    scout, prior_preflight = scout_gate(base, cfg)
    loaded = base.inputs(cfg, out); response, fcfg, world, mech, pilot, points = loaded
    pre = base.read(base.pin(out / 'PREFLIGHT.json'))
    base.require(pre['pilot_result'] == prior_preflight['pilot_result'] and pre['modes'] == prior_preflight['modes']
                 and pre['status'] == 'passed' and pre['new_molecular_calls'] == 0, 'Remaining preflight lineage differs')
    result.update(preflight=base.pin(out / 'PREFLIGHT.json'), modes=pilot['modes'], pilot_result=cfg['pilot_result'],
                  scout_result=cfg['scout_result'], frozen_scout_implementation=cfg['frozen_scout_implementation'],
                  original_common3=pilot['comparison']['common3'], force_check_result=pilot['force_check_result'],
                  physical_model='MACE_local_vacuum + GFN2_embedded_minus_isolated + existing_MM + fixed_FF_global_ALPB')
    evaluator = base.Evaluator(cfg, out, result, loaded)
    searchdir = out / 'Dy_search'; searchdir.mkdir()
    sr = base.search(evaluator, 'Dy', searchdir)
    result.update(rows=[{'metal': 'Dy', 'status': sr['status'], 'candidate_admitted': sr['candidate_admitted'],
                         'search': base.pin(searchdir / 'SEARCH.json'), 'points': evaluator.points}],
                  points=evaluator.points, requests=evaluator.requests, accepted_steps=sr['accepted_steps'],
                  status='complete' if sr['candidate_admitted'] else 'unavailable',
                  continuation_denominator=[{'metal': 'La', 'status': scout['rows'][0]['status'],
                                             'search': scout['rows'][0]['search'], 'source': cfg['scout_result']},
                                            {'metal': 'Dy', 'status': sr['status']}])
    if evaluator.native is not None:
        result.update(native_model_calls=evaluator.native.calls, native_device=evaluator.native.device,
                      model=evaluator.mcfg['model'], versions=evaluator.native.versions)
    base.require(result['new_composite_evaluations_attempted'] <= 24 and result['new_MACE_calls_attempted'] <= 24
                 and result['new_GFN2_calls_attempted'] <= 48 and result['new_MM_configurations_attempted'] <= 24
                 and result['new_MM_group_queries_attempted'] <= 96 and result['new_ALPB_calls_attempted'] <= 24,
                 'Remaining call ceiling exceeded')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('remaining', 'mm'), default='remaining')
    parser.add_argument('--config', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); base = load_base()
    base.require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    cfg = json.loads(args.config.read_text())
    base.require(base.pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation differs')
    base.verify(cfg['plan']); out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    tick = time.monotonic(); cpu = time.process_time()
    result = {'status': 'running', 'stage': args.stage, 'implementation': base.pin(__file__), 'config': base.pin(args.config),
              'plan': cfg['plan'], 'base': BASE, 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
              'settings': base.SETTINGS, 'new_composite_evaluations_attempted': 0, 'new_MACE_calls_attempted': 0,
              'new_GFN2_calls_attempted': 0, 'new_MM_configurations_attempted': 0, 'new_MM_group_queries_attempted': 0,
              'new_ALPB_calls_attempted': 0, 'new_native_model_loads': 0, 'observed_GFN2_child_CPU_seconds': 0.,
              'new_DFT_calls': 0, 'new_wholeprotein_SCC_calls': 0}
    try:
        if args.stage == 'mm': base.mm_stage(cfg, out, result)
        else: run(base, cfg, out, result)
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic() - tick; result['process_CPU_seconds'] = time.process_time() - cpu
    if args.stage == 'remaining':
        result['application_allocated_CPU_seconds'] = result['elapsed_seconds'] * int(os.environ['SLURM_CPUS_PER_TASK'])
        result['application_allocated_GPU_seconds'] = result['elapsed_seconds']
    base.write(out / 'RESULT.json', result)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
