#!/usr/bin/env python3
"""Saved-fixture geometry adapter. All numerical work requires actual Slurm."""
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ORDER = ('origin', 'adaptive_Ca', 'adaptive_La')


def read(p):
    return json.loads(Path(p).read_text())


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    if pin(p['path']) != p:
        raise ValueError('changed input: ' + p['path'])
    return Path(p['path'])


def dump(p, value):
    with Path(p).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def confidence_pins(model):
    p = Path(model['path'])
    stem = p.name.removesuffix('model.cif')
    return {name: pin(p.parent / (stem + name + '.json'))
            if (p.parent / (stem + name + '.json')).exists() else None
            for name in ('confidences', 'summary_confidences')}


def prepare(review, sha, workspace):
    review = Path(review).resolve()
    assert pin(review)['sha256'] == sha
    scope = read(review)
    assert len(scope['cases']) == 48 and sum(c['target'] for c in scope['cases']) == 9
    out = Path(workspace).resolve()
    out.mkdir(parents=True, exist_ok=False)
    olddir = ROOT / 'workspaces/pede_failure_diagnosis_20261013'
    oldm = read(olddir / 'geometry_v1/manifest.json')
    recovery = read(olddir / 'geometry_local_v2/manifest.json')
    oldresult = read(olddir / 'geometry_local_v2/execution/RESULT.json')
    oldcases = {c['case_id']: c for c in oldm['cases']}
    oldrows = {c['case_id']: c for c in oldresult['rows']}
    code = oldm['code'] | recovery['code']
    cases, inputs = [], {}
    for item in scope['cases']:
        r = item['record']
        m = read(verify(r['baseline_manifest']))
        c = m['cases'][r['baseline_manifest_case_index']]
        if c.get('execution_status') == 'reuse_actual_corrected_pool':
            pc = read(verify(c['prior_collection']))
            pm = read(verify(pc['manifest']))
            c = next(a for a in pm['cases'] if a['case_id'] == c['case_id'])
        n = c.get('native_source_record')
        if n is None:
            col = read(verify(c['collection']))
            n = next(a for a in col['cases'] if a['case_id'] == c['case_id'])
        prep = n['source']['origin_row']['source']
        source = prep['source']
        raw = source.get('original_source_structure', source['source_structure'])
        metal, spectator = source['raw_source_metal'], None
        audit_pin = source.get('source_geometry_audit')
        if audit_pin:
            audit = read(verify(audit_pin))
            ar = next(a for a in audit['rows'] if a['new_source_id'] == r['source_case_id'])
            raw, metal, spectator = ar['model'], ar['active_Ca'], ar['spectator_Ca']
        proposals = n.get('proposals') or {
            z: n['aliases']['adaptive_' + z]['proposal_receipt'] for z in ('Ca', 'La')}
        selected = {z: item['score']['pool']['operational'][z]['candidate'] for z in ('Ca', 'La')}
        row = {'case_key': item['case_key'], 'record_id': item['record_id'], 'target': item['target'],
               'case_id': r['source_case_id'], 'selection_id': r['selection_id'],
               'biological_group': r['biological_group'], 'family': r['family'],
               'source_conditioning_metal': r['source_conditioning_metal'],
               'source_occupancy_state': source.get('source_occupancy_state', 'one source metal'),
               'state_signature': r['state_signature'], 'raw': raw, 'confidence': confidence_pins(raw),
               'metal': metal, 'spectator': spectator, 'pqq': source['pqq'], 'roles': source['roles'],
               'existing_admission': r['preparation_status'], 'mapping': r['maps']['Ca'],
               'geometries': {q: r['matrix']['Ca'][q]['xyz'] for q in ORDER},
               'selected': selected, 'proposals': proposals, 'raw_geometry_reuse': None}
        prior = oldcases.get(row['case_id'])
        if prior and all(row[k] == prior[k] for k in ('raw', 'metal', 'spectator', 'pqq', 'roles')):
            if 'raw' in oldrows[row['case_id']]:
                row['raw_geometry_reuse'] = {'collection': pin(olddir / 'geometry_local_v2/execution/RESULT.json'),
                                             'case_id': row['case_id']}
        for p in [raw, row['mapping'], row['state_signature'], *row['confidence'].values(),
                  *row['geometries'].values(), *proposals.values()]:
            if p:
                verify(p)
                inputs[p['path']] = p
        for p in proposals.values():
            result = read(verify(p))
            for name in ('manifest',):
                q = result[name]
                verify(q)
                inputs[q['path']] = q
        cases.append(row)
    for p in code.values():
        verify(p)
    m = {'protocol': 'inconclusive_saved_geometry_review_v1', 'driver': pin(__file__),
         'review': pin(review), 'plan': pin(HERE / 'PLAN.md'), 'parent_plan': scope['plan'],
         'helper': oldm['driver'], 'kinematic_replay_precedent': recovery['driver'],
         'python': oldm['python'], 'code': code, 'inputs': list(inputs.values()), 'cases': cases,
         'output': str(out / 'execution'), 'numerical_threads': 1, 'new_energy_calls': 0,
         'new_folds': 0, 'new_thresholds': False,
         'scope': 'all48 frozen source/context cases; nine targets; peer grouping remains explicit'}
    dump(out / 'manifest.json', m)
    print(json.dumps({'manifest': pin(out / 'manifest.json'), 'cases': len(cases),
                      'reused_raw_rows': sum(c['raw_geometry_reuse'] is not None for c in cases)}))


def raw_map(c, b):
    st = b.gemmi.read_structure(str(verify(c['raw'])))
    assert len(st) == 1
    amap, schema, metal = {}, None, None
    for chain in st[0]:
        for res in chain:
            if (chain.name, res.seqid.num, res.seqid.icode.strip()) == b.selkey(c['pqq'], '')[:3]:
                schema = b.pqq.identify_naming_schema(res)
            for a in res:
                k = (chain.name, res.seqid.num, res.seqid.icode.strip(), a.name.strip())
                assert k not in amap, ('duplicate raw atom', k)
                v = {'xyz': b.np.array([a.pos.x, a.pos.y, a.pos.z]),
                     'element': a.element.name, 'resname': res.name, 'plddt': float(a.b_iso)}
                amap[k] = v
                if k == b.selkey(c['metal'], c['metal']['atom']):
                    metal = v['xyz']
    assert schema is not None and metal is not None
    return amap, metal, schema


def closest_nonbonded(data, symbols, xyz, np):
    meta = data['source_atom_metadata']
    adjacency = {i: set() for i in range(len(xyz))}
    for i, j in data['bonds']:
        adjacency[i].add(j)
        adjacency[j].add(i)
    ids = [i for i, a in enumerate(meta) if a and symbols[i] not in ('H', 'D', 'Ca', 'La')]
    pairs = []
    for j, i in enumerate(ids):
        excluded = adjacency[i] | {k for a in adjacency[i] for k in adjacency[a]}
        for k in ids[j + 1:]:
            if k in excluded:
                continue
            pairs.append((float(np.linalg.norm(xyz[i] - xyz[k])), i, k))
    return [{'distance_A': d, 'physical_indices': [i, j], 'atoms': [meta[i], meta[j]]}
            for d, i, j in sorted(pairs)[:12]]


def execute(path):
    if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('SLURM_JOB_NODELIST'):
        raise RuntimeError('actual Slurm allocation required before scientific imports')
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
        os.environ[key] = '1'
    m = read(path)
    for p in [m['driver'], m['review'], m['plan'], m['parent_plan'], m['helper'],
              m['kinematic_replay_precedent'], m['python'], *m['code'].values(), *m['inputs']]:
        verify(p)
    out = Path(m['output'])
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    receipt = {'manifest': pin(path), 'job_id': os.environ['SLURM_JOB_ID'],
               'hostname': socket.gethostname(), 'status': 'running', 'new_energy_calls': 0,
               'new_folds': 0, 'numerical_threads': 1}
    dump(out / 'START.json', receipt)
    try:
        import numpy as np
        import gemmi
        sys.path.insert(0, str(ROOT / 'scripts'))
        from mace_site_kinematics import Kinematics
        import pqq_microstates as pqq
        import coordination_policy as policy
        import carve_with_pqq as carver
        spec = importlib.util.spec_from_file_location('archived_geometry_helper', m['helper']['path'])
        b = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(b)
        b.np, b.gemmi, b.pqq, b.policy, b.carver = np, gemmi, pqq, policy, carver
        reader, results, raw_cache = b.xyz_read, [], {}
        for c in m['cases']:
            r = {k: c[k] for k in ('case_key', 'record_id', 'target', 'case_id', 'selection_id',
                 'biological_group', 'family', 'source_conditioning_metal', 'source_occupancy_state',
                 'state_signature', 'raw', 'selected')}
            r.update(raw_source=c['raw'], status='unavailable', optimizer_evidence={})
            try:
                key = (c['raw']['sha256'], json.dumps({k: c[k] for k in
                       ('metal', 'spectator', 'pqq', 'roles')}, sort_keys=True))
                if key not in raw_cache:
                    if c['raw_geometry_reuse']:
                        rr = read(verify(c['raw_geometry_reuse']['collection']))
                        raw = next(a['raw'] for a in rr['rows'] if a['case_id'] == c['case_id'])
                        amap, metal, schema = raw_map(c, b)
                    else:
                        raw, amap, metal, schema = b.raw_case(c)
                    raw_cache[key] = raw, amap, metal, schema
                raw, amap, metal, schema = raw_cache[key]
                r['raw_metrics'], r['raw_geometry_reuse'] = raw, c['raw_geometry_reuse']
                data = read(verify(c['mapping']))['context']
                kin, cache, replay = Kinematics(data), {}, {}
                meta = data['source_atom_metadata']
                for name in ORDER:
                    rec = c['geometries'][name]
                    symbols, xyz = reader(rec)
                    if name == 'origin':
                        q = np.zeros(len(data['modes']))
                    else:
                        pr = c['proposals'][name.removeprefix('adaptive_')]
                        proposal = read(verify(pr))
                        assert proposal['status'] == 'proposal_available'
                        q = np.asarray(proposal['proposal']['full_q'])
                        r['optimizer_evidence'][name] = {
                            'receipt': pr, 'status': proposal['status'],
                            **{k: proposal.get(k) for k in ('optimizer', 'boundary_flag', 'residual',
                               'final_geometry', 'maximum_trial_heavy_extent_A', 'infeasible_trial_requests')},
                            'missing_recorded_fields': [k for k in ('boundary_flag', 'residual', 'final_geometry') if k not in proposal],
                            'full_q': proposal['proposal']['full_q'],
                            'active_q_radian': proposal['proposal'].get('active_q_radian')}
                    physical, rendered, _, _ = kin.evaluate(q)
                    replay[name] = {'saved_vs_kinematic_max_coordinate_difference_A': float(np.max(np.abs(xyz - rendered))),
                                    'physical_atom_count': len(physical), 'saved_atom_count': len(xyz)}
                    physical[kin.metal] = xyz[0]
                    psymbols = [v['element'] if v else symbols[0] for v in meta]
                    for link in data['core_links']:
                        i, kind, *rest = link
                        if kind == 'source':
                            j = rest[0]
                            assert symbols[i] == psymbols[j]
                            physical[j] = xyz[i]
                    cache[rec['path']] = psymbols, physical
                b.xyz_read = lambda rec: cache[rec['path']]
                r['local'] = b.local_case(c, schema, amap)
                r['physical_map_replay'] = replay
                r['physical_modes'] = data['modes']
                for name, metric in r['local'].items():
                    sy, pos = cache[c['geometries'][name]['path']]
                    metric['closest_mapped_nonbonded_heavy_pairs'] = closest_nonbonded(data, sy, pos, np)
                    metric['nonbonded_scope'] = 'physical source heavy atoms; excludes metal, H, and graph separation1/2; ranked distances, not clash labels'
                r['status'] = 'complete'
            except Exception as exc:
                r.update(status='partial_or_unavailable', reason=str(exc), traceback=traceback.format_exc())
            finally:
                b.xyz_read = reader
            results.append(r)
            dump(out / ('case_' + str(len(results)).zfill(3) + '.json'), r)
        result = {'protocol': m['protocol'], 'manifest': pin(path), 'terminal': True,
                  'declared_cases': len(m['cases']), 'declared_targets': 9, 'rows': results,
                  'new_energy_calls': 0, 'new_folds': 0,
                  'interpretation': 'saved descriptive geometry/confidence; no causal assignment, new quality cutoff, or recalibration'}
        dump(out / 'RESULT.json', result)
        receipt.update(status='complete' if all(r['status'] == 'complete' for r in results) else 'partial',
                       result=pin(out / 'RESULT.json'))
    except Exception as exc:
        receipt.update(status='failed', reason=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        u = resource.getrusage(resource.RUSAGE_SELF)
        receipt.update(wall_seconds=time.monotonic() - started, CPU_seconds=u.ru_utime + u.ru_stime,
                       peak_RSS_KiB=u.ru_maxrss)
        dump(out / 'EXECUTION.json', receipt)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='operation', required=True)
    f = sub.add_parser('prepare')
    f.add_argument('--review', required=True)
    f.add_argument('--sha', required=True)
    f.add_argument('--workspace', required=True)
    f = sub.add_parser('execute')
    f.add_argument('--manifest', required=True)
    a = p.parse_args()
    if a.operation == 'prepare':
        prepare(a.review, a.sha, a.workspace)
    else:
        execute(a.manifest)
