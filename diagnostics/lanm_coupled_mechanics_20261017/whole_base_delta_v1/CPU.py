#!/usr/bin/env python3
"""Four sequential default-solver GFN calls after the admitted fixed scout."""
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

PROTOCOL = 'lbt_whole_GFN_regional_vacuum_delta_CPU_v1'
BASE_SHA = 'f760cd799f326fe193d2a2809e4434000803e04438a3358de0210d829d32a5d1'
CALLS = (('7CCO_after_La', 'whole'), ('7CCO_before_Dy', 'whole'),
         ('7CCO_after_Dy', 'whole'), ('7CCO_after_Dy', 'region'))
IDS = tuple(s+'_'+t for s in ('7CCO', '7CCN') for t in ('before_La', 'after_La', 'before_Dy', 'after_Dy'))


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


def prepare(cfg, out, result):
    require(cfg['protocol'] == PROTOCOL and cfg['implementation']['sha256'] == pin(__file__)['sha256']
            and cfg['base_implementation']['sha256'] == BASE_SHA
            and cfg['declared_calls'] == [list(x) for x in CALLS], 'CPU continuation declaration differs')
    base = module('whole_base_delta_frozen', cfg['base_implementation'])
    base.verify_tree(cfg)
    bc = read(cfg['base_config'])
    require(bc['implementation']['sha256'] == BASE_SHA and cfg['design'] == bc['design']
            and cfg['inventory'] == bc['inventory'], 'Frozen parent scientific declaration differs')
    pf = read(cfg['preflight_result']); scout = read(cfg['scout_result']); sc = read(scout['config'])
    require(pf['stage'] == 'preflight' and pf['status'] == 'preflight_passed'
            and pf['scientific_binding'] == base.binding(bc) and pf['new_molecular_calls'] == 0,
            'Actual zero-call preflight unavailable')
    require(scout['protocol'] == base.PROTOCOL and scout['stage'] == 'run-scout'
            and scout['status'] == 'scout_passed' and scout['scientific_binding'] == base.binding(bc)
            and sc['preflight_result'] == cfg['preflight_result']
            and scout['scout_native_admitted'] and scout['scout_cost_gate_passed']
            and scout['scout_wall_seconds'] <= base.LIMITS['scout_wall_seconds']
            and scout['new_molecular_calls'] == scout['whole_GFN_calls_attempted'] == 1
            and not scout['model_loaded'], 'Actual one-call scout admission/cost unavailable')
    scratch = {}; d = out/'saved_preparation'; d.mkdir()
    cells, h = base.prepare(bc, d, scratch)
    result['manifest'] = scratch['manifest']
    manifest = read(result['manifest'])
    old_manifest = read(pf['manifest'])
    def signatures(value):
        return [{k: v for k, v in r.items() if k not in ('whole_xyz', 'region_xyz', 'mapping_replay')}
                | {'whole_xyz_sha256': r['whole_xyz']['sha256'], 'region_xyz_sha256': r['region_xyz']['sha256']}
                for r in value['cells']]
    require(signatures(manifest) == signatures(old_manifest), 'Eight prepared parent signatures differ')
    native = base.native_replay(scout['scout_native_result'],
        {'source_id': '7CCO', 'region': 'whole', 'metal': 'La', 'kind': 'isolated', 'atom_count': 256,
         'charge': -2, 'uhf': 0, 'multiplicity': 1, 'physical_multiplicity': 1,
         'valence_electrons': 750, 'pointcharge_count': 0, 'field_charge_e': 0.}, h)
    nt = read(native['task']); m = manifest['cells'][0]
    require(nt['physical_source'] == m['physical'] and nt['physical_ids'] == m['physical_ids']
            and nt['xyz']['sha256'] == m['whole_xyz']['sha256']
            and read(native['execution'])['wall_seconds'] == scout['scout_wall_seconds'],
            'Exact reused scout coordinates/cost differ')
    result['rows'] = [{'cell_id': r['cell_id'], 'source_id': r['source_id'], 'geometry': r['geometry'],
        'metal': r['metal'], 'state': r['state'], 'status': 'not_requested',
        'whole_GFN': None, 'regional_GFN': None} for r in manifest['cells']]
    result['rows'][0].update(status='whole_GFN_available', whole_GFN=scout['scout_native_result'])
    result['scout_reused_exactly'] = scout['scout_native_result']
    result['scout_cost'] = read(native['execution'])
    return base, cells, h, {r['cell_id']: r for r in manifest['cells']}


def evaluate(cfg, base, h, row, kind, out, result):
    cid = row['cell_id']; require((cid, kind) in CALLS, 'New native task not in four-call allowlist')
    key = cid+'__'+kind
    require(key not in [r['task_id'] for r in result['requests']] and len(result['requests']) < 4,
            'Duplicate request/four-call ceiling exceeded')
    d = out/key; d.mkdir()
    xyz = row[kind+'_xyz']; n = 256 if kind == 'whole' else 172
    (d/'core.xyz').write_bytes(verify(xyz).read_bytes())
    (d/'param_gfn2-xtb.txt').write_bytes(verify(h['runtime']['parameters']).read_bytes())
    staged = {'core.xyz': pin(d/'core.xyz'), 'parameters': pin(d/'param_gfn2-xtb.txt')}
    require(staged['core.xyz']['sha256'] == xyz['sha256'] and staged['parameters']['sha256'] == base.PARAM_SHA,
            'Staged scientific bytes differ')
    field = write(d/'FIELD.json', {'rows': [], 'external_field': False, 'physical_atom_count': 256})
    command = [str(verify(h['runtime']['executable']))]+base.ARGS
    task = {'task_id': key, 'source_id': row['source_id'], 'metal': row['metal'],
        'region': 'whole' if kind == 'whole' else 'baseline', 'kind': 'isolated',
        'atom_count': n, 'charge': -2, 'uhf': 0, 'multiplicity': 1,
        'physical_multiplicity': row['state'][kind]['MACE_physical_multiplicity'],
        'valence_electrons': row['state'][kind]['GFN_valence_electrons'], 'pointcharge_count': 0,
        'field_charge_e': 0., 'field': field, 'xyz': staged['core.xyz'],
        'physical_source': row['physical'], 'physical_ids': row['physical_ids'],
        'target_index': 255 if kind == 'whole' else 0, 'mapping': None if kind == 'whole' else row['mapping'],
        'runtime': h['runtime'], 'full_corrected_profile': h['pins']['GFN_profile'],
        'parameter_loading_source': h['pins']['parameter_loading_source'], 'command': command,
        'staged_files': staged, 'protocol': PROTOCOL, 'config': result['config'],
        'preflight_result': cfg['preflight_result'], 'scout_result': cfg['scout_result'],
        'electronic_state': 'published spin-independent GFN2 valence model; UHF0; no DFT ECP',
        'physical_state_qualified': False, 'derivative_qualified': False}
    task_pin = write(d/'TASK.json', task)
    native = {'status': 'unavailable_attempted', 'attempted': True, 'task': task_pin}
    request = {'task_id': key, 'cell_id': cid, 'scope': kind, 'status': 'unavailable_attempted', 'task': task_pin}
    result['requests'].append(request); result['new_molecular_calls'] += 1
    result['whole_GFN_calls_attempted' if kind == 'whole' else 'regional_GFN_calls_attempted'] += 1
    env = h['backend'].clean_env(staged['parameters']['path'])
    for k, v in h['runtime']['environment'].items():
        if v is None: env.pop(k, None)
        else: env[k] = v
    env['XTBPATH'] = str(d.resolve())
    try:
        execution = h['backend'].call(command, d, env)
        native['execution'] = pin(d/'EXECUTION.json'); base.verify_tree(staged)
        native.update(h['admit'].admit(d, task, execution, backend=h['backend'], gradients=h['gradient'], np=h['np']),
                      status='accepted')
        base.shape(native['qm_gradient']['hartree_per_bohr'], n, h['np'])
        request.update(status='accepted', execution=native['execution'])
    except Exception as exc:
        native.update(status='unavailable_attempted', failure=repr(exc), traceback=traceback.format_exc())
        if (d/'EXECUTION.json').is_file(): native['execution'] = pin(d/'EXECUTION.json')
    native['artifacts'] = {p.name: pin(p) for p in d.iterdir() if p.is_file()}
    request['native_result'] = write(d/'RESULT.json', native)
    target = next(r for r in result['rows'] if r['cell_id'] == cid)
    target['whole_GFN' if kind == 'whole' else 'regional_GFN'] = request['native_result']
    if native['status'] != 'accepted':
        target['status'] = 'unavailable_attempted'
        raise ValueError('First native/admission failure; remaining calls not requested: '+key)
    target['status'] = 'whole_GFN_available' if kind == 'whole' else 'whole_and_regional_GFN_available'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('run-development-cpu',), required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required before scientific imports')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic()
    result = {'protocol': PROTOCOL, 'stage': args.stage, 'status': 'execution_failure',
        'implementation': pin(__file__), 'config': pin(args.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
        'host': socket.gethostname(), 'declared_calls': [list(x) for x in CALLS], 'requests': [],
        'declared_cells': list(IDS), 'rows': [{'cell_id': cid, 'status': 'not_requested',
                                            'whole_GFN': None, 'regional_GFN': None} for cid in IDS],
        'new_molecular_calls': 0, 'whole_GFN_calls_attempted': 0, 'regional_GFN_calls_attempted': 0,
        'regional_MACE_calls_attempted': 0, 'model_loaded': False, 'ALPB_calls': 0, 'MM_calls': 0,
        'embedded_GFN_calls': 0, 'DFT_calls': 0, 'optimizer_steps': 0,
        'physical_state_qualified': False, 'derivative_qualified': False}
    try:
        cfg = read(result['config'])
        result.update(inventory=cfg['inventory'], preflight_result=cfg['preflight_result'], scout_result=cfg['scout_result'])
        base, _, h, rows = prepare(cfg, out, result)
        for cid, kind in CALLS:
            evaluate(cfg, base, h, rows[cid], kind, out, result)
        require(result['new_molecular_calls'] == 4 and result['whole_GFN_calls_attempted'] == 3
                and result['regional_GFN_calls_attempted'] == 1 and all(r['status'] == 'accepted' for r in result['requests']),
                'Exact four-call completion unavailable')
        result.update(status='development_cpu_complete', cumulative_new_molecular_calls_including_scout=5,
                      cumulative_whole_GFN_calls=4, cumulative_regional_GFN_calls=1,
                      cumulative_regional_MACE_calls=0, candidate_limit=base.LIMITS)
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick, peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  allocation={k: os.environ.get(k) for k in ('SLURM_CPUS_ON_NODE', 'SLURM_JOB_CPUS_PER_NODE', 'SLURM_JOB_GPUS')})
    write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json'), 'new_molecular_calls': result['new_molecular_calls']}))
    return 0 if result['status'] == 'development_cpu_complete' else 1


if __name__ == '__main__':
    sys.exit(main())
