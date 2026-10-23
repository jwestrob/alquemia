#!/usr/bin/env python3
"""Saved own-basis pair/cooperative comparison; no molecular calls."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time
import ARITHMETIC as arithmetic

HERE = Path(__file__).resolve().parent
HA = 627.509474
TERMS = ('electronic_energy_hartree', 'D4_energy_hartree', 'final_energy_hartree')


def require(ok, message):
    if not ok: raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path']) == p, 'Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p): return json.loads(verify(p).read_text())


def module(p, name):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def admit_reference(path, cell, panel_pin, panel, adapter):
    path = Path(path); r = json.loads(path.read_text())
    require(r['status'] == 'complete' and r['returncode'] == 0 and all(r['checks'].values()), 'Native admission failed')
    require(r['cell'] == cell and r['cell_id'] == cell['cell_id'] and r['panel'] == panel_pin and r['adapter'] == panel['adapter'], 'Reference lineage differs')
    for k in ('input', 'output', 'stderr', 'start', 'release'): verify(r[k])
    start = read(r['start']); release = read(r['release'])
    require(Path(r['start']['path']).resolve() == path.parent/'START.json', 'START location differs')
    require(start['cell'] == cell and start['panel'] == panel_pin and start['release'] == r['release'] and start['input'] == r['input'], 'START lineage differs')
    require(release['status'] == 'approved_matched_pair_36_cells' and release['panel'] == panel_pin and release['implementation'] == panel['adapter'], 'Reference release differs')
    require(release['cell_ids'] == [c['cell_id'] for c in panel['cells']] and release['maximum_native_calls'] == 36 and release['ranks'] == 4, 'Released denominator differs')
    for k in ('preparation', 'independent_review'): verify(release[k])
    require(path.resolve() == Path(release['execution_root']).resolve()/'tasks'/cell['cell_id']/'execution/RESULT.json', 'Native attempt differs')
    require(start['ranks'] == 4 and start['memory_MiB'] > 0 and start['maxcore_MiB'] == math.floor(.75*start['memory_MiB']/4), 'Allocation differs')
    require(start['command'] == [panel['orca']['path'], 'endpoint.inp', '--bind-to none'] and start['maximum_native_calls'] == r['native_program_invocations_attempted'] == 1, 'Invocation differs')
    require(verify(r['input']).read_text() == adapter.deck(panel, cell, 4, start['maxcore_MiB']), 'Executed deck differs')
    parsed = adapter.parse(verify(r['output']).read_text(), cell, r['returncode'])
    for k, v in parsed.items(): require(r[k] == v, 'Raw reference replay differs: '+k)
    require(parsed['status'] == 'complete' and all(parsed['checks'].values()), 'Raw reference not admitted')
    values = {TERMS[0]: parsed['electronic_energy_hartree'], TERMS[1]: parsed['D4_energy_hartree'], TERMS[2]: parsed['energy_hartree']}
    require(all(arithmetic.finite(v) for v in values.values()), 'Nonfinite reference term')
    require(abs(values[TERMS[0]] + values[TERMS[1]] - values[TERMS[2]]) < 1e-8, 'Reference dispersion closure differs')
    return dict(status='complete', result=pin(path), input=r['input'], output=r['output'], start=r['start'], release=r['release'],
                cell_id=cell['cell_id'], parent_geometry_id=cell['parent_geometry_id'], water_index=cell['water_index'],
                job_id=r['job_id'], allocated_CPU_seconds_proxy=r['allocated_cpu_seconds'],
                native_checks=parsed['checks'], printed_scf_criteria=parsed['printed_scf_criteria'],
                native_scf_mode=parsed['native_scf_mode'], native_scf_stop_signals=parsed['native_scf_stop_signals'],
                energy_stage_metadata=parsed['energy_stage_metadata'],
                final_COSX_correction_hartree=parsed['final_COSX_correction_hartree'], **values)


def reference_decomposition(panel, rows, originals):
    result = {}
    for term in TERMS:
        clusters = {g: originals[g][term] for g in arithmetic.GIDS}
        pairs = {g: {} for g in arithmetic.GIDS}; ions = {}
        for c in panel['cells']:
            if c['role'] == 'bare_ion': ions[c['metal']] = rows[c['cell_id']][term]
            else: pairs[c['parent_geometry_id']][c['water_index']-1] = rows[c['cell_id']][term]
        x = arithmetic.decompose(clusters, pairs, ions)
        result[term] = {'hartree': x, 'kcal_mol': {k: v*HA for k, v in x.items() if isinstance(v, (int, float))}}
    result['dispersion_closure_kcal_mol'] = {
        k: result[TERMS[2]]['kcal_mol'][k] - result[TERMS[0]]['kcal_mol'][k] - result[TERMS[1]]['kcal_mol'][k]
        for k in ('cluster_contrast', 'pair_contrast', 'cooperative_remainder')}
    return result


def unavailable(path, exc=None):
    p = Path(path)
    return {'status': 'missing' if not p.exists() else 'unavailable', 'path': str(p),
            'result': pin(p) if p.exists() else None, 'error': repr(exc) if exc else None}


def report(r):
    lines = ['# Matched pair and cooperative response', '', 'Status: **'+r['status']+'**.', '',
             f"Reference cells: {r['reference_complete']}/36; model dimers: {r['model_complete']}/34; original reference clusters: {r['original_complete']}/4.", '']
    if r['comparison'] is not None:
        q = r['reference_decomposition'][TERMS[2]]['kcal_mol']; m = r['model_decomposition']; e = r['comparison']
        lines += ['| Contrast | Reference | Model | Model minus reference |', '|---|---:|---:|---:|']
        for key, label in [('cluster_contrast', 'Full shell C'), ('pair_contrast', 'Individual pairs P'), ('cooperative_remainder', 'Cooperative remainder B')]:
            lines.append(f'| {label} | {q[key]:+.9f} | {m[key]:+.9f} | {e[key]:+.9f} |')
        lines += ['', 'All values are kcal/mol. C=P+B; the ion difference is subtracted once from the paired sum.', '']
    else: lines += ['The full comparison is unavailable; no missing terms were replaced or partial averages reported.', '']
    lines += ['This is the declared literal own-basis decomposition. The remainder includes finite-basis and dispersion nonadditivity; it is not a unique polarization or charge-transfer energy.',
              '', 'The completed metal–whole-shell CP cycle stays separate. Original targets and parameters are unchanged. No force, free-energy, model-selection or validation claim follows.',
              '', 'Endpoint evidence, electronic/D4 splits, native stopping metadata, unavailable slots and costs are in RESULT.json.']
    return '\n'.join(lines)+'\n'


def collect(config, out):
    start = time.monotonic(); cpu = time.process_time(); cfg = json.loads(Path(config).read_text())
    out = Path(out); out.mkdir(parents=True, exist_ok=False)
    r = {'status': 'unavailable', 'implementation': pin(__file__), 'configuration': pin(config),
         'reference_rows': {}, 'model_rows': {}, 'original_clusters': {}, 'reference_complete': 0,
         'model_complete': 0, 'original_complete': 0, 'reference_decomposition': None,
         'model_decomposition': None, 'comparison': None, 'error': None, 'new_native_calls': 0,
         'fit_targets_changed': False, 'model_selected': False, 'fresh_validation_released': False}
    try:
        require(cfg['implementation'] == pin(__file__) and cfg['arithmetic'] == pin(arithmetic.__file__) and cfg['specification'] == pin(HERE/'SPEC.md'), 'Collection source differs')
        adapter = module(cfg['reference_adapter'], 'matched_reference_saved')
        panel = adapter.load(verify(cfg['reference_panel'])); root = read(panel['root_inputs'])
        require(panel['adapter'] == cfg['reference_adapter'], 'Reference adapter differs')
        require(set(cfg['reference_result_paths']) == {c['cell_id'] for c in panel['cells']}, 'Reference denominator differs')
        r['reference_rows'] = {cid: unavailable(path) for cid, path in cfg['reference_result_paths'].items()}
        r['model_rows'] = {c['cell_id']: dict(unavailable(cfg['model_result_path']), cell_id=c['cell_id'])
                           for c in panel['cells'] if c['role'] == 'dimer'}
        cp_adapter = adapter.helper; cp_panel = cp_adapter.load(verify(root['cp_panel'])); data = read(root['consumed_input_data'])
        old_collector = module(cfg['original_collector'], 'matched_original_saved')
        old_cp = read(root['cp_collection']); require(old_cp['status'] == 'coordination_cp_collection_complete', 'Original admitted comparison unavailable')
        for g in arithmetic.GIDS:
            try: r['original_clusters'][g] = old_collector.admit_original(g, cp_panel, data, cp_adapter)
            except Exception as exc: r['original_clusters'][g] = unavailable(panel['saved_complexes'][g]['result']['path'], exc)
        for c in panel['cells']:
            path = cfg['reference_result_paths'][c['cell_id']]
            try: r['reference_rows'][c['cell_id']] = admit_reference(path, c, cfg['reference_panel'], panel, adapter)
            except Exception as exc: r['reference_rows'][c['cell_id']] = unavailable(path, exc)
        r['reference_complete'] = sum(v['status'] == 'complete' for v in r['reference_rows'].values())
        r['original_complete'] = sum(v['status'] == 'complete' for v in r['original_clusters'].values())
        if r['reference_complete'] == 36 and r['original_complete'] == 4:
            r['reference_decomposition'] = reference_decomposition(panel, r['reference_rows'], r['original_clusters'])
            require(abs(r['reference_decomposition'][TERMS[2]]['kcal_mol']['cluster_contrast'] - old_cp['cycle']['kcal_mol']['original_total']) < 1e-9, 'Original reference contrast differs')
        model_adapter = module(cfg['model_adapter'], 'matched_model_saved')
        model_cfg = read(cfg['model_config']); manifest = read(cfg['model_manifest'])
        require(manifest['panel'] == cfg['reference_panel'], 'Model/reference preparation differs')
        model_result = json.loads(Path(cfg['model_result_path']).read_text())
        admitted = model_adapter.admit_result(model_cfg, manifest, model_result)
        rows = admitted['rows']; ids = {c['cell_id'] for c in panel['cells'] if c['role'] == 'dimer'}
        require(set(rows) == ids, 'Model denominator differs')
        r['model_rows'] = rows; r['model_complete'] = len(rows)
        r['model_ion_zero_proof'] = admitted['analytic_bare_ion']
        require(set(admitted['analytic_bare_ion']) == {'La', 'Dy'} and all(
            v['status'] == 'analytically_zero' and v['energy_kcal_mol'] == 0. and v['missing_data_substitution'] is False
            for v in admitted['analytic_bare_ion'].values()), 'Unsupported model ion energy')
        r['model_result'] = pin(cfg['model_result_path']); r['model_native_cost'] = model_result['cost']
        require(all(admitted['saved_cluster_models'][g]['energy_kcal_mol'] == panel['baseline_model_energies'][g]
                    for g in arithmetic.GIDS), 'Frozen model cluster energies differ')
        pairs = {g: {} for g in arithmetic.GIDS}
        for c in panel['cells']:
            if c['role'] == 'dimer': pairs[c['parent_geometry_id']][c['water_index']-1] = rows[c['cell_id']]['energy_kcal_mol']
        r['model_decomposition'] = arithmetic.decompose(panel['baseline_model_energies'], pairs, {'La': 0., 'Dy': 0.})
        require(abs(r['model_decomposition']['cluster_contrast']-old_cp['cycle']['kcal_mol']['unchanged_model']) < 1e-9, 'Original model contrast differs')
        if r['reference_decomposition'] is not None:
            r['comparison'] = arithmetic.compare(r['reference_decomposition'][TERMS[2]]['kcal_mol'], r['model_decomposition'])
            require(abs(r['comparison']['closure_residual']) < 1e-8, 'Error closure differs')
            r['status'] = 'matched_pair_collection_complete'
        else: r['status'] = 'matched_pair_collection_incomplete'
    except Exception as exc: r.update(status='matched_pair_collection_unavailable', error=repr(exc))
    r['cost'] = {'reference_allocated_CPU_seconds_proxy': math.fsum(v.get('allocated_CPU_seconds_proxy', 0.) for v in r['reference_rows'].values()),
                 'proxy_definition': 'native elapsed wall seconds times allocated ranks; scheduler cost is recorded separately',
                 'unavailable_endpoint_costs_not_included': True, 'collector_wall_seconds': time.monotonic()-start,
                 'collector_CPU_seconds': time.process_time()-cpu, 'collector_new_native_calls': 0, 'local_threads': 1}
    (out/'RESULT.json').write_text(json.dumps(r, indent=2, allow_nan=False)+'\n')
    (out/'REPORT.md').write_text(report(r))
    print(json.dumps({k: r[k] for k in ('status', 'reference_complete', 'model_complete', 'error')}))
    return r


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--config', required=True); ap.add_argument('--out', required=True)
    args = ap.parse_args(); result = collect(args.config, args.out)
    raise SystemExit(0 if result['status'] == 'matched_pair_collection_complete' else 1)
