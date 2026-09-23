"""Read-only six-angle selection and exact warm-start/reuse feasibility."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from adaptive_force_diagnostic import preview, project
from mace_site_kinematics import Kinematics
from audit_c5ax_response import data, CID

PROBES = (
    'a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1',
    'a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4',
    'a8r3s4-pqq-la_model__conditioned_La__seed-1_sample-1',
    'a8r3s4-pqq-la_model__conditioned_La__seed-1_sample-3',
)


def audit(c5ax, strict225, output):
    a = read_json(c5ax); s = read_json(strict225)
    ids = [CID if i == 1 else CID+'__conditioned_La__seed-1_sample-'+str(i)
           for i in range(5)] + list(PROBES)
    sr = {r['case_id']: r for r in s['rows']}
    sp = {r['case_id']: r for r in s['actual_pools']}
    ar = {r['case_id']: r for r in a['rows']}
    rows = []
    for cid in ids:
        if cid in ar:
            source = ar[cid]['response']['pool_collection']
            strict = (ar[cid]['strict32']['branches']['fresh'] if cid == CID
                      else sp[cid])
        else:
            source = sr[cid]['precision_collection']; strict = sp[cid]
        coll = data(source); pool_manifest = data(coll['manifest'])
        mp = pool_manifest['source_manifest']; m = data(mp)
        ts = {z: next(t for t in m['tasks'] if (t['case_id'], t['metal']) == (cid, z))
              for z in ('Ca', 'La')}
        maps = [data(ts[z]['mapping'])['context'] for z in ('Ca', 'La')]
        if maps[0] != maps[1] or ts['Ca']['active_mode_ids'] != ts['La']['active_mode_ids']:
            raise InvalidArtifact('old common source/selection differs')
        mapping = maps[0]; kin = Kinematics(mapping)
        results = {}; origins = {}; receipt_pins = {}
        for z, t in ts.items():
            rp = Path(mp['path']).parent/'proposals'/t['task_id']/'result.json'
            r = read_json(rp)
            if r['manifest'] != mp or r['status'] != 'proposal_available':
                raise InvalidArtifact('exact successful prior proposal absent')
            receipt_pins[z] = record(rp); results[z] = r
            pt = r['origin']; f = np.load(verify(pt['forces']), allow_pickle=False)
            origins[z] = project(mapping, pt['full_q'], f)
        if not np.allclose(origins['Ca'][0], origins['La'][0], atol=1e-12, rtol=0):
            raise InvalidArtifact('origin nuclear coordinates differ')
        choice = preview(mapping['modes'], origins['Ca'][1],
                         origins['Ca'][3], origins['La'][3], maximum=6)
        names = [v['id'] for v in choice['selected']]
        if len(names) != 6 or names[:4] != ts['Ca']['active_mode_ids']:
            raise InvalidArtifact('six independent modes or original prefix unavailable')
        indices = [[v['id'] for v in mapping['modes']].index(n) for n in names]
        endpoints = {}; cells = []
        for z, t in ts.items():
            r = results[z]; pt = r['proposal']; q = np.asarray(pt['full_q'])
            if np.any(q[indices[4:]] != 0):
                raise InvalidArtifact('new coordinates were already moving')
            pq, context, j, jc = kin.evaluate(q)
            if not np.allclose(context, [v[1:] for v in xyz(verify(pt['coordinate']))], atol=1e-12, rtol=0):
                raise InvalidArtifact('warm-start mapping differs from archived geometry')
            f = np.load(verify(pt['forces']), allow_pickle=False)
            _, _, raw, normalized, _ = project(mapping, q, f)
            extent = float(np.linalg.norm(pq[kin.heavy]-kin.positions[kin.heavy], axis=1).max())
            if extent > .8+1e-7 or max(abs(q[indices])) > .8+1e-12:
                raise InvalidArtifact('old candidate outside unchanged admissible domain')
            endpoints[z] = {
                'task': t, 'proposal_receipt': receipt_pins[z], 'warm_start': pt,
                'full_q': q.tolist(), 'active_q_radian': q[indices].tolist(),
                'maximum_heavy_extent_A': extent, 'prior_boundary': r['boundary_flag'],
                'six_gradient_kcal_mol_rad': raw[indices].tolist(),
                'six_normalized_load_kcal_mol_A': normalized[indices].tolist(),
                'old_optimizer': r['optimizer'], 'old_wall_seconds': r['wall_seconds'],
            }
            for g in ('origin', 'adaptive_Ca', 'adaptive_La'):
                cell = strict['matrix'][z][g]
                if cell['status'] != 'complete': raise InvalidArtifact('old strict cell incomplete')
                # Strict scalar scoring must refer to the same physical MACE cell.
                old = next(c for c in coll['cases'] if c['case_id'] == cid)['matrix'][z][g]
                if cell['components']['MACE_eV'] != old['components']['MACE_eV']:
                    raise InvalidArtifact('strict baseline MACE energy differs')
                ax, bx = xyz(verify(cell['xyz'])), xyz(verify(old['xyz']))
                if [v[0] for v in ax] != [v[0] for v in bx] or not np.allclose(
                        [v[1:] for v in ax], [v[1:] for v in bx], atol=1e-12, rtol=0):
                    raise InvalidArtifact('strict baseline coordinates differ')
                low = cell.get('low', cell.get('actual_low'))
                if set(low) != {'vacuum', 'alpb'}: raise InvalidArtifact('strict medium missing')
                for medium, entry in low.items():
                    verify(entry['receipt']); verify(entry['output']); verify(entry['manifest'])
                    cells.append({'metal': z, 'candidate': g, 'medium': medium, 'source': entry})
        rows.append({'case_id': cid, 'status': 'feasible_read_only', 'source_pool': source,
                     'proposal_manifest': mp, 'selection': choice, 'selected_ids': names,
                     'active_indices': indices, 'added_ids': names[4:],
                     'endpoints': endpoints, 'strict_reuse_cells': cells,
                     'old_strict_R': strict['pool']['operational']['composite_R_model_kcal_mol']})
    out = {'status': 'proposed_not_executed', 'scope': 'five C5AX La sources and four declared probes',
           'c5ax_audit': record(c5ax), 'strict225': record(strict225), 'implementation': record(__file__),
           'rows': rows, 'declared_sources': ids, 'exact_strict_reuse_cells': sum(len(r['strict_reuse_cells']) for r in rows),
           'future_counts': {'single_start_searches': 18, 'maximum_new_cross_MACE': 18,
                             'maximum_new_strict_GFN2': 72, 'new_q0': 0, 'DFT': 0},
           'new_molecular_calls': 0, 'new_proposals': 0, 'new_thresholds': 0}
    write_new(output, out)
    return {'output': record(output), 'sources': len(rows), 'reuse_cells': out['exact_strict_reuse_cells'],
            'added_modes': {r['case_id']: r['added_ids'] for r in rows}, 'new_molecular_calls': 0}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for k in ('c5ax', 'strict225', 'output'): p.add_argument('--'+k, type=Path, required=True)
    print(json.dumps(audit(**vars(p.parse_args())), indent=2))
