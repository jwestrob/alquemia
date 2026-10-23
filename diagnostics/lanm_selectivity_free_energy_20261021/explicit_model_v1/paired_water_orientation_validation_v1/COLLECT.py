#!/usr/bin/env python3
"""Saved-only merge of four original water reference/model cells; retain all16 slots."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

def require(ok, message):
    if not ok: raise ValueError(message)

def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):
    require(pin(p['path']) == p, 'Changed input: ' + p['path'])
    return Path(p['path'])

def read(p): return json.loads(verify(p).read_text())

def module(p, name):
    s = importlib.util.spec_from_file_location(name, verify(p))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
    return m

def save(p, obj):
    with Path(p).open('x') as f: json.dump(obj, f, indent=2, allow_nan=False); f.write('\n')
    return pin(p)

def collect(config_path, out):
    start, cpu = time.monotonic(), time.process_time()
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')), 'One-thread collection required')
    cfg = json.loads(Path(config_path).read_text())
    require(cfg['implementation'] == pin(__file__), 'Collector changed')
    out = Path(out); out.mkdir(parents=True, exist_ok=False)
    meta = module(cfg['metadata_helper'], 'water_metadata')
    score = module(cfg['water_score'], 'water_fixed_score')
    require(cfg['stage_score']['sha256'] == '1c4bb0fa4b899991d3ca064480a10f199f2e8668469146389a1be403c257f104', 'Stage scorer changed')
    require(cfg['body_helper']['sha256'] == 'bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819', 'Rigid body helper changed')
    verify(cfg['stage_score']); verify(cfg['body_helper']); verify(cfg['root_declaration'])
    require(read(cfg['scoring_controls'])['status'] == 'passed', 'Subset controls unavailable')
    data, stage, body = meta.metadata()
    ref = module(cfg['reference_source'], 'water_reference_wrapper')
    model_adapter = module(cfg['model_source'], 'water_model_wrapper')
    for key in ('reference_config','reference_manifest','reference_release','reference_inventory','model_config','model_bound','model_release','training'):
        verify(cfg[key])
    refpin = ref.collect(cfg['reference_config']['path'], cfg['reference_manifest']['path'], cfg['reference_release']['path'], cfg['reference_inventory']['path'], out/'REFERENCE_ADMISSION.json')
    references = read(refpin)
    bound, modelcfg = read(cfg['model_bound']), read(cfg['model_config'])
    model_adapter.bound_inputs(modelcfg, bound)
    water_ids = set(model_adapter.CELLS)
    require(water_ids == set(references['rows']) == {g for g in stage.CELLS if g.startswith('validation_aquo8_')}, 'Exact four water IDs required')
    admission = {'status': 'saved_model_admission_partial', 'rows': {g: {'status': 'unavailable', 'model': None, 'reason': 'missing_result'} for g in water_ids}, 'complete': 0}
    result_path = Path(cfg['model_result_path'])
    model_result = None
    if result_path.exists():
        model_result = pin(result_path)
        try:
            actual = read(model_result)
            require(actual['bound'] == cfg['model_bound'] and actual['release'] == cfg['model_release'], 'Actual model release/binding differs')
            admission = model_adapter.admit_result(modelcfg, bound, actual)
        except Exception as exc:
            admission['reason'] = repr(exc)
    modelpin = save(out/'MODEL_ADMISSION.json', admission)
    model = {gid: None for gid in stage.CELLS}
    for gid in water_ids:
        row = data['rows'][gid]
        require(row['geometry'] == bound['rows'][gid]['geometry'] == references['rows'][gid]['geometry'], 'Paired model/reference geometry differs')
        require(row['bodies'] == bound['rows'][gid]['bodies'] and row['input']['bonds_zero_based'] == bound['rows'][gid]['input']['bonds_zero_based'], 'Model metadata differs')
        row['reference'] = references['rows'][gid]['reference']
        model[gid] = admission['rows'][gid]['model']
    training = read(cfg['training'])
    require(cfg['training']['sha256'] == '99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb', 'Consumed training changed')
    anchors = {'model': {}, 'reference': {}}
    anchor_provenance = {}
    # Re-admit the two consumed central reference anchors from their raw outputs.
    native = module(read(cfg['reference_config'])['native_adapter'], 'unchanged_reference_native')
    for metal in ('La', 'Dy'):
        gid = 'aquo8_central_' + metal
        old = training['rows'][gid]
        raw = read(old['reference_result'])
        require(raw['status'] == 'complete' and raw['returncode'] == 0 and all(raw['checks'].values()), 'Central reference unavailable')
        panel = native.load_panel(verify(raw['panel']))
        cell = next(c for c in panel['cells'] if c['cell_id'] == 'ecp_pbe0__' + gid)
        replay = native.collect_one(panel, cell, Path(old['reference_result']['path']).parent, 0)
        require(all(raw[k] == v for k, v in replay.items()), 'Central raw replay differs')
        require(old['reference']['energy'] == raw['energy_hartree'] * native.HA and old['reference']['forces'] == raw['forces_kcal_mol_A'], 'Training anchor differs')
        anchors['reference'][metal] = {'energy': old['reference']['energy']}
        anchors['model'][metal] = {'energy': bound['anchor_model'][gid]['energy']}
        anchor_provenance[gid] = {'reference': old['reference_result'], 'model_transfer': bound['transfer_result']}
    result = score.evaluate(data, model, body, anchors, stage)
    result.update(implementation=pin(__file__), configuration=pin(config_path), reference_admission=refpin,
                  model_admission=modelpin, model_result=model_result, anchor_provenance=anchor_provenance,
                  native_calls_this_collection=0, reference_allocated_CPU_seconds=references['allocated_native_cpu_seconds'],
                  wall_seconds=time.monotonic()-start, CPU_seconds=time.process_time()-cpu,
                  interpretation='Original reserved water-only orientation test against provisional frozen-f ECP-PBE0-D4; donor12 remain unavailable, full validation incomplete.')
    rp = save(out/'RESULT.json', result)
    lines = ['# Independent water orientation test', '',
             '**Water subset: ' + ('PASS' if result['water_subset_passed'] else 'FAIL' if result['water_failure_observed'] else 'INCOMPLETE') + '. Full16 validation remains incomplete.**', '',
             '| Check | Available | Passed | Failed |', '|---|---:|---:|---:|']
    for group, count in result['counts'].items():
        lines.append(f"| {group} | {count['available']}/{count['declared']} | {count['passed']} | {count['failed']} |")
    lines.extend(['', '| Energy comparison | Model | Reference | Error (kcal/mol) | Pass |', '|---|---:|---:|---:|---|'])
    for group in ('balanced_exchange', 'orientation_work'):
        for name, v in result['comparisons'][group].items():
            if v['status'] == 'complete': lines.append(f"| {name} | {v['model_kcal_mol']:.6f} | {v['reference_kcal_mol']:.6f} | {v['error_kcal_mol']:.6f} | {v['passed']} |")
            else: lines.append(f'| {name} | unavailable | unavailable | unavailable | — |')
    lines.extend(['', '| Rigid response | RMS error | Maximum component error | Pass |', '|---|---:|---:|---|'])
    for name, v in result['comparisons']['rigid_response'].items():
        if v['status'] == 'complete': lines.append(f"| {name} | {v['rigid_RMS_kcal_mol_A']:.6f} | {v['maximum_body_component_error_kcal_mol_A']:.6f} | {v['passed']} |")
        else: lines.append(f'| {name} | unavailable | unavailable | — |')
    lines.extend(['', 'Limits are unchanged: energy absolute error ≤1 kcal/mol; rigid RMS ≤2 and maximum component ≤5 kcal/mol/Å. Every one of11 comparisons is required.', '',
                  'The other12 Stage A cells have explicit unavailable records; all44 comparisons and original24-cell inventory remain visible. No protein sampling, refit or production release follows from this subset.', '',
                  'This tests transfer to four new water orientations against the provisional reference surface. It does not qualify bulk hydration, donor interactions, protein selectivity, electronic-state accuracy or exchange free energy.', '',
                  f"Reference allocated cost: {references['allocated_native_cpu_seconds']} CPU-seconds. This saved-only collection makes zero molecular calls; scheduler/model cost is retained in terminal receipts.", '', f"Result: `{rp['sha256']}`."])
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'result': rp, 'status': result['status'], 'water_subset_passed': result['water_subset_passed']}))

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--config', required=True); p.add_argument('--out', required=True)
    a = p.parse_args(); collect(a.config, a.out)
