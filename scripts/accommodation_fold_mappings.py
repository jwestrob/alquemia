"""Geometry-only terminal donor maps for all primary reference folds."""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path
import time
import numpy as np
from affordable_common import InvalidArtifact, paired, read_json, record, verify, write_new, xyz
from accommodation_reference_geometry import mapping_case
from accommodation_torsion_profiles import role_mode
from mace_site_kinematics import Kinematics

PROTOCOL = 'primary225_terminal_carboxylate_mapping_v1'


def population(preparation, existing):
    prep = read_json(preparation); old = read_json(existing)
    if old['source_preparation'] != record(preparation):
        raise InvalidArtifact('old maps have a different prepared source')
    rows = [c for c in prep['cases'] if c['source']['primary_evaluation_pool']]
    if len(rows) != 225 or len({c['case_id'] for c in rows}) != 225:
        raise InvalidArtifact('primary reference denominator differs')
    if any(c['source']['canonical_coordinate_match'] for c in rows):
        raise InvalidArtifact('canonical replay entered primary transfer pool')
    if sum(c['status'] == 'prepared' for c in rows) != 208:
        raise InvalidArtifact('prepared source count differs from frozen plan')
    reusable = {r['case_id']: r['mapping_check'] for r in old['rows'] if r['mapping_status'] == 'supported'}
    return prep, rows, reusable


def supported_case(case, mapping):
    """Validate old/new maps on their actual common nuclear coordinate measure."""
    gp = read_json(verify(mapping['mapping'])); kin = Kinematics(gp['context'])
    rep = case['representations']['context']; origins = rep['endpoints']
    ca, la = [origins[z] for z in ('Ca', 'La')]
    paired(verify(la['xyz']), verify(ca['xyz']), la['charge'], ca['charge'])
    expected = np.asarray([a[1:] for a in xyz(verify(la['xyz']))])
    actual = kin.evaluate(np.zeros(len(kin.modes)))[1]
    delta = float(np.max(np.abs(actual-expected)))
    if delta > 1e-12:
        raise InvalidArtifact('physical map differs from context source')
    parent = read_json(verify(case['core']['parent']))
    selector = parent['fixed_core']['requested_roles']['extra_acidic_ligand_homolog']
    resname = selector['resname'] if isinstance(selector, dict) else selector.split(':')[1][:3]
    roles = ['anchor_glutamate'] + (['extra_acidic_ligand_homolog'] if resname == 'ASP' else [])
    modes = {role: role_mode(parent, role) for role in roles}
    names = [m['id'] for m in kin.modes]
    for role, mode in modes.items():
        i = names.index(mode); selected = kin.modes[i]
        required = {'OE1', 'OE2'} if role == 'anchor_glutamate' else {'OD1', 'OD2'}
        meta = gp['context']['source_atom_metadata']
        moving = {meta[j]['atom'] for j in selected['moving_indices'] if meta[j] is not None}
        if selected['unit'] != 'radian' or not required <= moving:
            raise InvalidArtifact('incomplete terminal carboxylate physical coordinate')
    src = case['source']
    return {'case_id': case['case_id'], 'source': case['core'], 'origins': origins,
            'maps': {z: mapping['mapping'] for z in ('Ca', 'La')}, 'preparation': rep['preparation'],
            'modes': modes, 'label_scope': case['label_scope'], 'root_case_id': src['root_case_id'],
            'biological_group': src['biological_group'], 'source_conditioning_metal': src['source_conditioning_metal'],
            'expected_class_for_later_report_only': src['expected_class'], 'canonical_coordinate_match': False,
            'primary_evaluation_pool': True, 'source_metadata': src, 'all_evidence_consumed': True,
            'map_origin_max_abs_difference_A': delta,
            'extra_homolog_residue': resname, 'active_coordinate_rule': 'Glu_chi3_and_actual_extra_Asp_chi2'}


def build_one(payload):
    case, old, topology, destination = payload
    if old is None:
        mapping = mapping_case(case, topology, destination); reused = False
    else:
        mapping = old; reused = True
    result = {'case_id': case['case_id'], 'status': 'mapping_unsupported', 'reused_mapping': reused,
              'mapping_check': mapping, 'design_case': None}
    try:
        if mapping['status'] != 'supported':
            raise InvalidArtifact(mapping.get('reason') or 'physical mapping unavailable')
        result.update(status='supported', design_case=supported_case(case, mapping))
    except (OSError, ValueError, KeyError) as exc:
        result['reason'] = str(exc)
    return result


def prepare(preparation, existing, engine_manifest, agreement, output, workers):
    if workers < 1 or (workers > 1 and not os.environ.get('SLURM_JOB_ID')):
        raise InvalidArtifact('parallel geometry preparation requires an allocation')
    if workers > int(os.environ.get('SLURM_CPUS_ON_NODE', workers)):
        raise InvalidArtifact('mapping workers exceed allocated CPUs')
    start = time.monotonic(); prep, selected, reusable = population(preparation, existing)
    engine = read_json(engine_manifest); verify(engine['orca']); verify(prep['config']['topology'])
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    payloads = [(c, reusable.get(c['case_id']), prep['config']['topology'], out/'maps')
                for c in selected if c['status'] == 'prepared']
    if workers == 1:
        built = list(map(build_one, payloads))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            built = list(pool.map(build_one, payloads))
    by_id = {r['case_id']: r for r in built}; rows = []; cases = []
    for case in selected:
        src = case['source']
        row = {'case_id': case['case_id'], 'root_case_id': src['root_case_id'],
               'biological_group': src['biological_group'], 'source_conditioning_metal': src['source_conditioning_metal'],
               'expected_class_for_later_report_only': src['expected_class'], 'canonical_coordinate_match': False,
               'primary_evaluation_pool': True, 'source': src, 'preparation_status': case['status'],
               'preparation_reason': case.get('reason'), 'status': 'preparation_unavailable',
               'mapping_check': None, 'reused_mapping': False}
        if case['case_id'] in by_id:
            checked = by_id[case['case_id']]
            row.update({k: v for k, v in checked.items() if k != 'design_case'})
            if checked['design_case'] is not None:
                cases.append(checked['design_case'])
        rows.append(row)
    elapsed = time.monotonic()-start
    result = {'protocol_id': PROTOCOL, 'agreement': record(agreement), 'source_preparation': record(preparation),
              'existing_mappings': record(existing), 'engine_manifest': record(engine_manifest),
              'implementation': record(__file__), 'model': prep['config']['model'], 'software': prep['config']['software'],
              'orca': engine['orca'], 'topology': prep['config']['topology'], 'cases': cases, 'rows': rows,
              'denominator': 225, 'original_supported': 208, 'original_unavailable': 17,
              'mapping_supported': len(cases), 'status_counts': dict(Counter(r['status'] for r in rows)),
              'maps_reused': sum(r['reused_mapping'] for r in rows),
              'maps_new_attempted': sum(not r['reused_mapping'] for r in built),
              'new_molecular_calls': 0, 'source_coordinates_changed': False, 'baseline_changed': False,
              'wall_seconds': elapsed, 'workers': workers, 'slurm_job_id': os.environ.get('SLURM_JOB_ID')}
    write_new(out/'design.json', result)
    return {k: result[k] for k in ('denominator', 'mapping_supported', 'status_counts', 'maps_reused', 'maps_new_attempted', 'new_molecular_calls', 'wall_seconds', 'slurm_job_id')}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('preparation', 'existing', 'engine-manifest', 'agreement', 'output'):
        p.add_argument('--'+key, required=True)
    p.add_argument('--workers', type=int, required=True)
    print(json.dumps(prepare(**vars(p.parse_args())), indent=2))


if __name__ == '__main__': main()
