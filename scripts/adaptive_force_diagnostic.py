"""Read-only physical force projections on consumed, pinned proposal archives."""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL
from mace_site_kinematics import Kinematics

PROTOCOL = 'archived_common_physical_load_diagnostic_v1'
RANK_TOL = 1e-5


def project(mapping, q, forces):
    kin = Kinematics(mapping)
    p, context, j, jc = kin.evaluate(q)
    forces = np.asarray(forces, dtype=float)
    if forces.shape != context.shape or not np.isfinite(forces).all():
        raise InvalidArtifact('force dimensions or values invalid')
    physical = j[:, kin.heavy].reshape(len(kin.modes), -1)
    lengths = np.linalg.norm(physical, axis=1)
    if np.any(lengths <= 1e-12):
        raise InvalidArtifact('zero physical heavy-atom mode')
    raw = -np.einsum('mij,ij->m', jc, forces) * EV_TO_KCAL
    return context, physical / lengths[:, None], raw, raw / lengths, lengths


def preview(modes, unit_vectors, ca, la, maximum=4):
    """Select common angular primitives; no energies, classes or thresholds enter."""
    eligible = [i for i, m in enumerate(modes) if m['kind'] in
                ('sidechain_torsion', 'peptide_crankshaft') and m['unit'] == 'radian']
    remaining = set(eligible); basis = []; covectors = []; selected = []; redundant = []
    loads = np.asarray([ca, la])
    for slot in range(min(maximum, len(eligible))):
        candidates = []
        for i in sorted(remaining, key=lambda k: modes[k]['id']):
            v = unit_vectors[i].copy(); g = loads[:, i].copy()
            # Two passes limit roundoff for nearly dependent real torsions.
            for _ in range(2):
                for b, bg in zip(basis, covectors):
                    a = float(v @ b); v -= a * b; g -= a * bg
            norm = float(np.linalg.norm(v))
            if norm <= RANK_TOL:
                redundant.append(modes[i]['id']); continue
            v /= norm; g /= norm
            metric = abs(g[0]-g[1]) if slot % 2 == 0 else float(np.max(np.abs(g)))
            candidates.append((metric, modes[i]['id'], i, v, g, norm))
        remaining -= {i for i in remaining if modes[i]['id'] in redundant}
        if not candidates: break
        metric, name, i, v, g, norm = sorted(candidates, key=lambda x: (-x[0], x[1]))[0]
        basis.append(v); covectors.append(g); remaining.remove(i)
        selected.append({'id': name, 'criterion': 'differential' if slot % 2 == 0 else 'individual',
                         'Ca_kcal_mol_A': float(g[0]), 'La_kcal_mol_A': float(g[1]),
                         'differential_kcal_mol_A': float(g[0]-g[1]),
                         'independent_geometric_fraction': norm})
    return {'selected': selected, 'redundant': sorted(set(redundant)),
            'eligible_count': len(eligible), 'rank_tolerance': RANK_TOL,
            'method_applied_to_optimization': False}


def endpoint(manifest, m, t, stage):
    rp = Path(manifest).parent / 'proposals' / t['task_id'] / 'result.json'
    row = {'task_id': t['task_id'], 'case_id': t['case_id'], 'metal': t['metal'],
           'stage': stage, 'status': 'unavailable', 'reason': None,
           'force_hamiltonian': 'native_vacuum_OMOL', 'composite_gradient': None}
    if not rp.exists(): return row | {'reason': 'proposal receipt missing'}, None
    r = read_json(rp)
    if r['manifest'] != record(manifest): raise InvalidArtifact('proposal belongs to another manifest')
    row.update(proposal_status=r['status'], boundary_flag=r.get('boundary_flag'),
               optimizer=r.get('optimizer'), proposal_receipt=record(rp))
    point = r.get(stage)
    if point is None: return row | {'reason': r.get('reason', stage+' unavailable')}, None
    native = read_json(verify(point['MACE'])); request = read_json(verify(native['request']))
    atoms = xyz(verify(point['coordinate']))
    if (native['status'] != 'complete' or native['model'] != m['model'] or
            request['xyz'] != point['coordinate'] or request['charge'] != t['charge'] or
            request['multiplicity'] != t['multiplicity']):
        raise InvalidArtifact('archived force method or state differs')
    data = read_json(verify(t['mapping']))['context']; q = np.asarray(point['full_q'])
    forces = np.load(verify(point['forces']), allow_pickle=False)
    if native['forces'] != point['forces']: raise InvalidArtifact('force receipt differs')
    context, vectors, raw, normed, lengths = project(data, q, forces)
    if not np.allclose(context, [a[1:] for a in atoms], atol=1e-9, rtol=0):
        raise InvalidArtifact('force source coordinates differ from mapping')
    if not np.allclose(raw[t['active_indices']], point['gradient_kcal_mol_rad'], atol=1e-8, rtol=0):
        raise InvalidArtifact('archived active gradient projection does not replay')
    active = set(t['active_mode_ids'])
    rows = [{'id': mode['id'], 'kind': mode['kind'], 'unit': mode['unit'],
             'active_in_prior_proposal': mode['id'] in active,
             'raw_gradient_kcal_per_unit': float(raw[i]), 'physical_norm_A_per_unit': float(lengths[i]),
             'normalized_gradient_kcal_mol_A': float(normed[i])} for i, mode in enumerate(data['modes'])]
    angular = [r for r in rows if r['kind'] != 'metal_translation']
    prior = [abs(r['normalized_gradient_kcal_mol_A']) for r in angular if r['active_in_prior_proposal']]
    omitted = [abs(r['normalized_gradient_kcal_mol_A']) for r in angular if not r['active_in_prior_proposal']]
    metal = [r['normalized_gradient_kcal_mol_A'] for r in rows if r['kind'] == 'metal_translation']
    row.update(status='available', modes=rows, forces=point['forces'], coordinate=point['coordinate'],
               q=point['full_q'], prior_active_max_load=max(prior, default=None),
               omitted_angular_max_load=max(omitted, default=None), metal_translation_load_norm=float(np.linalg.norm(metal)))
    return row, {'vectors': vectors, 'normalized': normed, 'mapping': data,
                 'coordinates': context, 'lengths': lengths}


def distribution(values):
    values = [v for v in values if v is not None]
    return {'n': len(values), 'min': min(values), 'median': float(np.median(values)),
            'max': max(values)} if values else {'n': 0, 'min': None, 'median': None, 'max': None}


def analyze(manifests, agreement, output):
    verify(record(agreement)); archives = []; all_endpoints = []; all_pairs = []
    for manifest in manifests:
        m = read_json(manifest); endpoint_rows = []; pairs = []; values = {}; inventories = Counter()
        for t in m['tasks']:
            for stage in ('origin', 'proposal'):
                row, data = endpoint(manifest, m, t, stage)
                endpoint_rows.append(row); values[t['case_id'], t['metal'], stage] = (row, data)
        for c in m['cases']:
            row = {'case_id': c['case_id'], 'status': 'unavailable', 'selection_preview': None}
            ca, la = [values[c['case_id'], z, 'origin'] for z in ('Ca', 'La')]
            if all(v[1] is not None for v in (ca, la)):
                a, b = ca[1], la[1]
                if (not np.allclose(a['coordinates'], b['coordinates'], atol=1e-12, rtol=0) or
                        a['mapping'] != b['mapping']): raise InvalidArtifact('joint origin is not common')
                modes = a['mapping']['modes']; inventories.update(v['kind'] for v in modes)
                loads = a['normalized'] - b['normalized']
                angular = [i for i, v in enumerate(modes) if v['kind'] != 'metal_translation']
                ordered = sorted(angular, key=lambda i: (-abs(loads[i]), modes[i]['id']))
                selected = preview(modes, a['vectors'], a['normalized'], b['normalized'])
                row.update(status='available', selection_preview=selected,
                           origin_mode_ids=[v['id'] for v in modes], normalized_delta=loads.tolist(),
                           highest_angular_differential=modes[ordered[0]]['id'] if ordered else None,
                           highest_differential_was_active=ca[0]['modes'][ordered[0]]['active_in_prior_proposal'] if ordered else None,
                           preview_new_modes=[v['id'] for v in selected['selected'] if v['id'] not in
                                              [x['id'] for x in ca[0]['modes'] if x['active_in_prior_proposal']]])
            pairs.append(row)
        for c in m.get('unavailable_cases', []):
            pairs.append({'case_id': c['case_id'], 'status': 'unavailable', 'selection_preview': None,
                          'reason': c.get('preparation_reason')})
        summaries = {}
        for stage in ('origin', 'proposal'):
            rows = [r for r in endpoint_rows if r['stage'] == stage and r['status'] == 'available']
            summaries[stage] = {'available': len(rows), 'endpoint_denominator': len(m['tasks']),
                'boundary_flags': sum(r['boundary_flag'] is True for r in rows),
                'active_load': distribution([r['prior_active_max_load'] for r in rows]),
                'omitted_angular_load': distribution([r['omitted_angular_max_load'] for r in rows]),
                'metal_block_load': distribution([r['metal_translation_load_norm'] for r in rows]),
                'omitted_exceeds_active': sum(r['omitted_angular_max_load'] > r['prior_active_max_load'] for r in rows)}
        good = [r for r in pairs if r['status'] == 'available']
        archives.append({'manifest': record(manifest), 'case_denominator': len(pairs),
                         'common_force_pairs': len(good), 'mode_inventory': dict(inventories), 'stages': summaries,
                         'highest_angular_differential_omitted': sum(not r['highest_differential_was_active'] for r in good),
                         'preview_adds_modes': sum(bool(r['preview_new_modes']) for r in good)})
        all_endpoints.extend(endpoint_rows); all_pairs.extend(pairs)
    result = {'protocol_id': PROTOCOL, 'agreement': record(agreement), 'implementation': record(__file__),
              'archives': archives, 'endpoints': all_endpoints, 'pairs': all_pairs, 'new_molecular_calls': 0,
              'affinity_score': None, 'optimizer_launched': False, 'labels_used_for_selection': False,
              'interpretation': 'vacuum diagnostic; normalized force is not relaxation energy or uncertainty'}
    write_new(output, result)
    return archives


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', action='append', required=True)
    p.add_argument('--agreement', required=True); p.add_argument('--output', required=True)
    a = p.parse_args(); print(json.dumps(analyze(a.manifest, a.agreement, a.output), indent=2))


if __name__ == '__main__': main()
