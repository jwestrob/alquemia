#!/usr/bin/env python3
"""One missing regional MACE forward, then saved-only four-cell assembly."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import socket
import sys
import time
import traceback
from types import SimpleNamespace

PROTOCOL = 'lbt_whole_GFN_regional_vacuum_delta_v1'
CPU_SHA = '46c8d49c8a0fd25adc52853ce8a4bd218c42a331070a05cf7a72a2c864df7bdf'
MACE_SHA = '686eba84fd3dfcf28724d2d4e7957528c987445d7137b18c88a07dfec055f661'
IDS = tuple(s+'_'+t for s in ('7CCO', '7CCN') for t in ('before_La', 'after_La', 'before_Dy', 'after_Dy'))
MISSING = '7CCO_after_Dy'
COMPONENTS = ('GFN_full', 'MACE_region', 'negative_GFN_region', 'fixed_FF_ALPB')
HA = 627.509474
EV = 23.06054783061903
BOHR = .529177210903


def require(ok, message):
    if not ok: raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as stream: sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: '+record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')
    return pin(path)


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def binding(cfg):
    return {k: cfg[k]['sha256'] if k == 'implementation' else cfg[k] for k in
            ('protocol', 'implementation', 'specification', 'design', 'inventory',
             'cpu_implementation', 'cpu_config', 'cpu_result', 'hybrid_implementation',
             'preflight_result', 'scout_result')}


def signatures(manifest):
    return [{k: v for k, v in r.items() if k not in ('whole_xyz', 'region_xyz', 'mapping_replay')}
            | {'whole_xyz_sha256': r['whole_xyz']['sha256'], 'region_xyz_sha256': r['region_xyz']['sha256']}
            for r in manifest['cells']]


def prepare(cfg, out, result):
    require(cfg['protocol'] == PROTOCOL and cfg['implementation']['sha256'] == pin(__file__)['sha256']
            and cfg['cpu_implementation']['sha256'] == CPU_SHA
            and cfg['hybrid_implementation']['sha256'] == MACE_SHA, 'Completion source declaration differs')
    cpu = module('whole_base_cpu_frozen', cfg['cpu_implementation'])
    cc = read(cfg['cpu_config']); parent = read(cfg['cpu_result'])
    require(cc['implementation']['sha256'] == CPU_SHA and parent['implementation']['sha256'] == CPU_SHA
            and parent['config'] == cfg['cpu_config'], 'CPU source/config binding differs')
    verify(parent['implementation'])
    require(parent['protocol'] == cpu.PROTOCOL and parent['stage'] == 'run-development-cpu'
            and parent['status'] == 'development_cpu_complete', 'Complete actual CPU development panel required')
    for key in ('design', 'inventory', 'preflight_result', 'scout_result'):
        require(cfg[key] == cc[key], 'CPU scientific binding differs: '+key)
    for key in ('inventory', 'preflight_result', 'scout_result'):
        require(parent[key] == cfg[key], 'CPU result lineage differs: '+key)
    require(parent['new_molecular_calls'] == 4 and parent['whole_GFN_calls_attempted'] == 3
            and parent['regional_GFN_calls_attempted'] == 1 and parent['regional_MACE_calls_attempted'] == 0
            and not parent['model_loaded'] and parent['cumulative_new_molecular_calls_including_scout'] == 5,
            'Exact CPU call budget differs')
    require(all(parent[k] == 0 for k in ('ALPB_calls', 'MM_calls', 'embedded_GFN_calls', 'DFT_calls', 'optimizer_steps')),
            'Unexpected CPU molecular component')
    require([r['cell_id'] for r in parent['rows']] == list(IDS)
            and all(r['status'] == 'not_requested' and r['whole_GFN'] is None and r['regional_GFN'] is None
                    for r in parent['rows'][4:]), 'Eight CPU slots or unrequested challenge differ')
    requests = parent['requests']
    require([(r['cell_id'], r['scope']) for r in requests] == list(cpu.CALLS)
            and all(r['status'] == 'accepted' for r in requests), 'Four default-solver requests differ')
    d = out/'saved_cpu_preparation'; d.mkdir(); scratch = {}
    base, cells, h, rows = cpu.prepare(cc, d, scratch)  # Saved only; CPU.evaluate is never called.
    base.verify_tree(cfg)
    require(signatures(read(scratch['manifest'])) == signatures(read(parent['manifest'])),
            'CPU/current eight-state manifests differ')
    result['manifest'] = scratch['manifest']
    require(parent['scout_reused_exactly'] == scratch['scout_reused_exactly']
            == parent['rows'][0]['whole_GFN'], 'Original scout was not reused exactly')
    native = {}
    for r in parent['rows'][:4]:
        cid = r['cell_id']; m = rows[cid]
        require(all(r[k] == m[k] for k in ('source_id', 'geometry', 'metal', 'state')),
                'CPU row state differs')
        require(r['status'] == ('whole_and_regional_GFN_available' if cid == MISSING else 'whole_GFN_available')
                and (r['regional_GFN'] is None) == (cid != MISSING), 'CPU component availability differs')
        native[cid] = {}
        for kind in ('whole', 'region') if cid == MISSING else ('whole',):
            record = r['whole_GFN' if kind == 'whole' else 'regional_GFN']; n = 256 if kind == 'whole' else 172
            expected = {'source_id': m['source_id'], 'metal': m['metal'], 'kind': 'isolated',
                'region': 'whole' if kind == 'whole' else 'baseline', 'atom_count': n, 'charge': -2,
                'uhf': 0, 'multiplicity': 1, 'physical_multiplicity': m['state'][kind]['MACE_physical_multiplicity'],
                'valence_electrons': m['state'][kind]['GFN_valence_electrons'], 'pointcharge_count': 0,
                'field_charge_e': 0., 'physical_source': m['physical'], 'physical_ids': m['physical_ids']}
            raw = base.native_replay(record, expected, h)
            task = read(raw['task'])
            require(task['xyz']['sha256'] == m[kind+'_xyz']['sha256'], 'Exact CPU native geometry differs')
            if kind == 'region': require(task['mapping'] == m['mapping'], 'Native regional map differs')
            if cid != IDS[0]:
                request = next(x for x in requests if (x['cell_id'], x['scope']) == (cid, kind))
                require(request['native_result'] == record and request['task'] == raw['task']
                        and request['execution'] == raw['execution']
                        and task['config'] == parent['config'], 'CPU request/native receipt differs')
            native[cid][kind] = {'pin': record, 'raw': raw}
    result['rows'] = [{k: r[k] for k in ('cell_id', 'source_id', 'geometry', 'metal', 'physical', 'positions',
                                        'mapping', 'physical_ids', 'state')} | {'status': 'not_requested'}
                      for r in rows.values()]
    result['parent_results'] = {'preflight': cfg['preflight_result'], 'scout': cfg['scout_result'], 'CPU': cfg['cpu_result']}
    result['parent_costs'] = {'scout_native': parent['scout_cost'], 'CPU_result': cfg['cpu_result'],
        'CPU_native_executions': [read(x['execution']) for x in requests],
        'CPU_elapsed_seconds': parent['elapsed_seconds'], 'CPU_allocation': parent['allocation'],
        'scheduler_utilization_not_inferred': True}
    return base, cells, h, rows, native


def mace_task(row, result):
    return {'task_id': MISSING+'__region_MACE', 'cell_id': MISSING, 'source_id': '7CCO',
        'region': 'baseline', 'metal': 'Dy', 'atom_count': 172, 'charge': -2, 'multiplicity': 6,
        'target_index': 0, 'physical_source': row['physical'], 'physical_ids': row['physical_ids'],
        'mapping': row['mapping'], 'xyz': row['region_xyz'], 'protocol': PROTOCOL,
        'config': result['config'], 'cpu_result': result['parent_results']['CPU'],
        'physical_state_qualified': False, 'derivative_qualified': False}


def run_gpu(cfg, cells, h, rows, out, result):
    require(cfg.get('gpu_result') is None, 'A prior GPU result cannot be rerun')
    cell = cells[MISSING]; row = rows[MISSING]
    require(cell['local'] is None and row['local_components'] is None
            and row['state']['region']['MACE_physical_multiplicity'] == 6, 'Only missing after-Dy regional MACE is allowed')
    d = out/'regional_MACE'; d.mkdir()
    (d/'core.xyz').write_bytes(verify(row['region_xyz']).read_bytes())
    task = mace_task(row, result); task['xyz'] = pin(d/'core.xyz')
    task_pin = write(d/'TASK.json', task)
    request = {'task': task_pin, 'cell_id': MISSING, 'status': 'unavailable', 'attempted': False}
    result['request'] = request
    hybrid = module('hybrid_forward_frozen', cfg['hybrid_implementation'])
    warm = module('native_MACE_loader_frozen', h['pins']['native_MACE_helper'])
    mcfg = read(h['pins']['MACE_config'])
    context = SimpleNamespace(native=None, warm=warm, mcfg=mcfg, alpb_helper=h['solvent'],
        cfg={'model': h['pins']['model']}, cost={'MACE_model_loads': 0, 'MACE_calls_attempted': 0})
    raw = {'status': 'unavailable', 'attempted': False, 'task': task_pin,
           'implementation': result['implementation'], 'native_forward': cfg['hybrid_implementation'],
           'native_loader': h['pins']['native_MACE_helper'], 'native_config': h['pins']['MACE_config']}
    try:
        # Exactly the frozen 172-atom native call; no HybridEvaluator initialization or mixed evaluator.
        value = hybrid.HybridEvaluator.mace(context, task, cell['core_symbols'], cell['core_x'])
        require(context.native.calls == context.cost['MACE_calls_attempted'] == 1, 'One-forward ceiling differs')
        raw.update(value, status='accepted', attempted=True)
        raw['forces_artifact'] = write(d/'FORCES_eV_A.json', raw['forces_eV_A'])
        raw['native_runtime'] = {'versions': context.native.versions, 'device': context.native.device,
            'parameter_versions_unchanged': True, 'native_forward_count': context.native.calls,
            'peak_allocated_bytes': context.native.torch.cuda.max_memory_allocated(),
            'peak_reserved_bytes': context.native.torch.cuda.max_memory_reserved(),
            'total_device_bytes': context.native.torch.cuda.get_device_properties(0).total_memory}
        request.update(status='accepted', attempted=True)
    except Exception as exc:
        raw.update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
        raise
    finally:
        count = context.cost['MACE_calls_attempted']
        raw['attempted'] = request['attempted'] = bool(count)
        result['new_molecular_calls'] = result['regional_MACE_calls_attempted'] = count
        result['model_loaded'] = context.native is not None
        result['MACE_cost'] = dict(context.cost)
        raw['cost'] = dict(context.cost)
        raw['artifacts'] = {p.name: pin(p) for p in d.iterdir() if p.is_file()}
        result['regional_MACE'] = request['native_result'] = write(d/'RESULT.json', raw)
        result['rows'][3].update(status=('regional_MACE_available' if raw['status'] == 'accepted' else
            'unavailable_attempted' if count else 'unavailable_before_forward'), regional_MACE=result['regional_MACE'])
    result['rows'][3].update(status='regional_MACE_available', regional_MACE=result['regional_MACE'])
    result.update(status='development_gpu_complete', cumulative_new_molecular_calls=6,
                  cumulative_whole_GFN_calls=4, cumulative_regional_GFN_calls=1, cumulative_regional_MACE_calls=1)


def saved_mace(cfg, rows, h, result):
    parent = read(cfg['gpu_result']); pc = read(parent['config'])
    require(parent['implementation']['sha256'] == cfg['implementation']['sha256']
            and pc['implementation']['sha256'] == cfg['implementation']['sha256']
            and binding(pc) == binding(cfg) == parent['scientific_binding'], 'GPU source/scientific binding differs')
    verify(parent['implementation'])
    require(parent['protocol'] == PROTOCOL and parent['stage'] == 'run-development-gpu'
            and parent['status'] == 'development_gpu_complete' and parent['model_loaded']
            and parent['new_molecular_calls'] == parent['regional_MACE_calls_attempted'] == 1
            and parent['cumulative_new_molecular_calls'] == 6, 'Complete one-forward GPU result required')
    require(all(parent[k] == 0 for k in ('whole_GFN_calls_attempted', 'regional_GFN_calls_attempted',
            'ALPB_calls', 'MM_calls', 'embedded_GFN_calls', 'DFT_calls', 'optimizer_steps')), 'Unexpected GPU component calls')
    require([r['cell_id'] for r in parent['rows']] == list(IDS)
            and all(r['status'] == 'not_requested' for r in parent['rows'][4:]), 'GPU eight-slot challenge boundary differs')
    raw = read(parent['regional_MACE']); hbase = rows[MISSING]; task = read(raw['task'])
    expected = mace_task(hbase, {'config': parent['config'], 'parent_results': {'CPU': cfg['cpu_result']}})
    require({k: v for k, v in task.items() if k != 'xyz'} == {k: v for k, v in expected.items() if k != 'xyz'}
            and task['xyz']['sha256'] == expected['xyz']['sha256'], 'Saved regional MACE input/state differs')
    for p in raw['artifacts'].values(): verify(p)
    require(raw['status'] == 'accepted' and raw['attempted'] and raw['model'] == h['pins']['model']
            and raw['force_convention'] == 'negative vacuum energy derivative'
            and raw['native_forward'] == cfg['hybrid_implementation']
            and raw['native_loader'] == h['pins']['native_MACE_helper']
            and raw['native_config'] == h['pins']['MACE_config']
            and raw['cost']['MACE_calls_attempted'] == raw['native_runtime']['native_forward_count'] == 1,
            'Saved native MACE lineage/forward differs')
    require(raw['native_runtime']['versions'] == h['pins']['expected_MACE_versions']
            and raw['native_runtime']['parameter_versions_unchanged'], 'MACE software/parameter provenance differs')
    h['solvent'].state(raw['input_state_check'], 172, -2, 6, 0)
    f = h['np'].asarray(raw['forces_eV_A'], dtype=float)
    require(f.shape == (172, 3) and h['np'].isfinite(f).all() and math.isfinite(raw['energy_eV'])
            and read(raw['forces_artifact']) == raw['forces_eV_A'], 'Saved raw force/energy replay differs')
    require(parent['request']['native_result'] == parent['regional_MACE']
            and parent['request']['task'] == raw['task'] and parent['request']['status'] == 'accepted',
            'GPU request/native result differs')
    result['parent_results']['GPU'] = cfg['gpu_result']
    result['parent_costs']['GPU_result'] = cfg['gpu_result']
    result['parent_costs']['GPU_native'] = {'cost': raw['cost'], 'runtime': raw['native_runtime'],
        'forward_wall_seconds': raw['wall_seconds'], 'elapsed_seconds': parent['elapsed_seconds'], 'allocation': parent['allocation']}
    return parent['regional_MACE'], raw


def collect(cfg, cells, h, rows, native, out, result):
    new_mace_pin, new_mace = saved_mace(cfg, rows, h, result)
    np = h['np']; completed = []
    result['partial_assembled_artifacts'] = []
    for cid in IDS[:4]:
        cell = cells[cid]; m = rows[cid]; full = native[cid]['whole']; local = cell['local']
        if cid == MISSING:
            mace_pin, mace = new_mace_pin, new_mace
            isolated_pin, isolated = native[cid]['region']['pin'], native[cid]['region']['raw']
            gm = -np.asarray(mace['forces_eV_A'], dtype=float)*EV
            gi = np.asarray(isolated['qm_gradient']['hartree_per_bohr'], dtype=float)*HA/BOHR
            zeros = np.zeros((len(cell['mapping']['field']), 3))
            gm_full = h['mapper'].scatter(cell['mapping'], gm, zeros)
            minus_gi_full = -h['mapper'].scatter(cell['mapping'], gi, zeros)
        else:
            mace_pin, isolated_pin = local['MACE'], local['GFN']
            mace, isolated = read(mace_pin), read(isolated_pin)
            gm_full = np.asarray(local['component_gradients_kcal_mol_A']['MACE_region'], dtype=float)
            minus_gi_full = np.asarray(local['component_gradients_kcal_mol_A']['negative_GFN_region'], dtype=float)
        require(m.get('local_components') is None or local['solver_exception_certificate'] is None,
                'No solver exception belongs to the four 7CCO development cells')
        n = {'GFN_full_Eh': full['raw']['energy_hartree'], 'MACE_region_eV': mace['energy_eV'],
             'GFN_region_Eh': isolated['energy_hartree'], 'ALPB_Eh': cell['solvent']['energy_hartree']}
        require(all(math.isfinite(v) for v in n.values()), 'Native energy ledger nonfinite')
        components = {'GFN_full': np.asarray(full['raw']['qm_gradient']['hartree_per_bohr'], dtype=float)*HA/BOHR,
            'MACE_region': gm_full, 'negative_GFN_region': minus_gi_full,
            'fixed_FF_ALPB': np.asarray(cell['solvent']['gradient_kcal_mol_A'], dtype=float)}
        require(all(a.shape == (256, 3) and np.isfinite(a).all() for a in components.values()), 'Full256 component gradient unavailable')
        candidate = sum(components[k] for k in COMPONENTS)
        base = components['GFN_full']+components['fixed_FF_ALPB']
        energies = (n['GFN_full_Eh']*HA, n['MACE_region_eV']*EV, -n['GFN_region_Eh']*HA, n['ALPB_Eh']*HA)
        row = {k: m[k] for k in ('cell_id', 'source_id', 'geometry', 'metal', 'physical', 'positions',
                                 'mapping', 'physical_ids', 'state')}
        row.update(status='complete', component_receipts={'GFN_full': full['pin'], 'MACE_region': mace_pin,
            'GFN_region': isolated_pin, 'ALPB': cell['solvent']['receipt']}, solver_exception_certificate=None,
            native_energies=n, component_gradients_kcal_mol_A={k: a.tolist() for k, a in components.items()},
            candidate_gradient_kcal_mol_A=candidate.tolist(), base_gradient_kcal_mol_A=base.tolist(),
            candidate_force_kcal_mol_A=(-candidate).tolist(), base_force_kcal_mol_A=(-base).tolist(),
            candidate_energy_kcal_mol=math.fsum(energies), base_energy_kcal_mol=math.fsum((energies[0], energies[3])),
            physical_state_qualified=False, derivative_qualified=False, new_molecular_calls=0)
        row['assembled_artifact'] = write(out/(cid+'.json'), row)
        result['partial_assembled_artifacts'].append({'cell_id': cid, 'artifact': row['assembled_artifact']})
        completed.append(row)
    # No four-cell comparison can escape from an incomplete or failed assembly.
    result['rows'][:4] = completed
    result.update(status='complete_7CCO', new_molecular_calls=0, call_accounting={
        'development_total': 6, 'whole_GFN': 4, 'regional_GFN': 1, 'regional_MACE': 1,
        'new_collect_calls': 0, 'remaining_conditional_7CCN_whole_GFN': 4, 'overall_ceiling': 10,
        'repeated_cached_components': 0, 'ALPB': 0, 'MM': 0, 'embedded_GFN': 0, 'DFT': 0},
        comparison_status='Not computed here; frozen saved comparator owns work/force triage',
        challenge_status='not_requested; requires separate frozen work triage and root continuation',
        missing_comparators={'before_Dy_whole_MACE_force': True, 'after_Dy_old_hybrid': True})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('run-development-gpu', 'collect_7CCO'), required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required before scientific imports')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic()
    result = {'protocol': PROTOCOL, 'stage': args.stage, 'status': 'execution_failure',
        'implementation': pin(__file__), 'config': pin(args.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
        'host': socket.gethostname(), 'declared_cells': list(IDS),
        'rows': [{'cell_id': cid, 'status': 'not_requested'} for cid in IDS],
        'new_molecular_calls': 0, 'whole_GFN_calls_attempted': 0, 'regional_GFN_calls_attempted': 0,
        'regional_MACE_calls_attempted': 0, 'model_loaded': False, 'ALPB_calls': 0, 'MM_calls': 0,
        'embedded_GFN_calls': 0, 'DFT_calls': 0, 'optimizer_steps': 0,
        'physical_state_qualified': False, 'derivative_qualified': False,
        'limitations': ['Conditional terminal chemistry and physical spin states remain unqualified.',
            'GFN2 UHF0 at 300K can have fractional occupations; no ground-state or derivative qualification.',
            'Regional vacuum delta is not mutual polarization; no absolute cross-element affinity claim.',
            'The 7CCN challenge, before-Dy whole-MACE force and after-Dy old hybrid remain unavailable here.']}
    try:
        cfg = read(result['config'])
        result.update(scientific_binding=binding(cfg), inventory=cfg['inventory'],
                      preflight_result=cfg['preflight_result'], scout_result=cfg['scout_result'])
        _, cells, h, rows, native = prepare(cfg, out, result)
        if args.stage == 'run-development-gpu': run_gpu(cfg, cells, h, rows, out, result)
        else: collect(cfg, cells, h, rows, native, out, result)
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick, peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        allocation={k: os.environ.get(k) for k in ('SLURM_CPUS_PER_TASK', 'SLURM_CPUS_ON_NODE', 'SLURM_JOB_CPUS_PER_NODE',
                    'SLURM_JOB_GPUS', 'SLURM_STEP_GPUS', 'CUDA_VISIBLE_DEVICES', 'SLURM_MEM_PER_CPU', 'SLURM_MEM_PER_NODE')})
    write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json'), 'new_molecular_calls': result['new_molecular_calls']}))
    return 0 if result['status'] in ('development_gpu_complete', 'complete_7CCO') else 1


if __name__ == '__main__':
    sys.exit(main())
