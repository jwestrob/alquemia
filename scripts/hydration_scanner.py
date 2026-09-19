"""Compare contextual water preparation with the existing masked MACE scorer.

The endpoint-dependent preparation requires a four-state interaction descriptor;
the old two-call cancellation is used only for identical paired water geometry.
"""
from __future__ import annotations
import argparse
import copy
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact, HA_TO_KCAL, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, accepted_attempt, check_atoms, write_xyz
from mace_file_checks import cached_file_checks

STAGE = 'hydration_scanner_comparison'
PROTOCOL = 'masked_omol_context_prepared_water_comparison_v1'


def source_key(a):
    return tuple(a.get(k, '') for k in ('chain', 'resnum', 'insertion_code', 'atom', 'altloc'))


def transfer(preparation, core_parent, core_rows, metal):
    """Copy only source-identified water H; enforce exact water-O coordinates."""
    mapping = {}
    for a in core_parent['atom_graph']['source_to_qm']:
        if a['kind'] == 'source' and a['source']['canonical_resname'] == 'HOH':
            key = source_key(a['source'])
            if key in mapping:
                raise InvalidArtifact('duplicate source water mapping')
            mapping[key] = a['qm_index']
    rows = xyz(verify(preparation['endpoints'][metal]['xyz']))
    expected_count = 3 * len(preparation['explicit_waters'])
    water = [(i, a) for i, a in enumerate(preparation['physical_atoms'])
             if a['kind'] == 'retained_site_water']
    if len(water) != expected_count or len(mapping) != expected_count:
        raise InvalidArtifact('whole/core water inventories differ')
    result = list(rows); moved = []; seen = set()
    for i, a in water:
        key = source_key(a['source_record'])
        if key in seen or key not in mapping:
            raise InvalidArtifact('missing or duplicate source water atom')
        seen.add(key); j = mapping[key]
        if core_rows[j][0] != a['element'] or rows[i][0] != a['element']:
            raise InvalidArtifact('water source mapping changes element')
        if a['element'] == 'O':
            if rows[i] != core_rows[j]:
                raise InvalidArtifact('water oxygen moved between whole/core inputs')
        elif a['element'] == 'H':
            result[i] = core_rows[j]; moved.append({'whole_index': i, 'core_index': j, 'source': list(key)})
        else:
            raise InvalidArtifact('unsupported water element')
    if any(result[i] != rows[i] for i in range(len(rows)) if i not in {a['whole_index'] for a in moved}):
        raise InvalidArtifact('non-water-H coordinate changed')
    check_atoms(result, preparation['endpoints'][metal]['charge'])
    return result, moved


def contrasts(energies, disconnected, *, common_geometry=False):
    total = (energies['Ca', 'bound'] - energies['La', 'bound']
             - disconnected['Ca'] + disconnected['La']) * EV_TO_KCAL
    if all((metal, 'detached') in energies for metal in ('Ca', 'La')):
        interaction = ((energies['Ca', 'bound'] - energies['Ca', 'detached'])
                       - (energies['La', 'bound'] - energies['La', 'detached'])) * EV_TO_KCAL
        environment = (energies['Ca', 'detached'] - energies['La', 'detached']
                       - disconnected['Ca'] + disconnected['La']) * EV_TO_KCAL
    else:
        if not common_geometry:
            raise InvalidArtifact('missing detached endpoints without a proven common geometry')
        interaction = total; environment = 0.  # Proven common-geometry cancellation, not an unavailable correction.
    return {'interaction_R_model_kcal': interaction, 'bound_total_R_model_kcal': total,
            'detached_environment_difference_model_kcal': environment,
            'component_closure_error_model_kcal': total - interaction - environment}


@cached_file_checks
def prepare(inventory, core_collection, normalized_collection, ggr_comparison, agreement, output):
    import mace_omol as omol
    from mace_omol_ablation_run import source_task, descriptor_model
    from mace_omol_prepared import audit_preparation
    from mace_omol_intact import geometry
    data = read_json(inventory); core = read_json(core_collection); normal = read_json(normalized_collection)
    cm = read_json(verify(core['manifest'])); nm = read_json(verify(normal['manifest']))
    alpha = [r for r in data['rows'] if r['case_id'].startswith('ALPHA_')]
    if {r['case_id'] for r in alpha} != {'ALPHA_1F6S', 'ALPHA_6IP9'} or core['status'] != 'complete' or normal['status'] != 'complete':
        raise InvalidArtifact('complete declared alpha source panel required')
    _, parent = source_task(alpha[0]['endpoints']['La'])
    _, out, m = omol.common(verify(parent['inventory']), verify(parent['software']), agreement, output, STAGE)
    p = out/'implementation'/'hydration_scanner.py'; shutil.copyfile(__file__, p); m['implementation'][p.name] = record(p)
    tasks = []; transformations = []
    for row in alpha:
        audit_preparation(verify(row['preparation'])); prep = read_json(verify(row['preparation']))
        case = row['case_id'].split('_')[1]
        for variant in ('normalized', 'contextual'):
            paired_rows = []
            for metal in ('Ca', 'La'):
                ct = next(t for t in cm['tasks'] if t['case'] == case and t['metal'] == metal)
                cp = read_json(verify(ct['parent']))
                if prep['source_preparation'] != ct['parent']:
                    raise InvalidArtifact('whole protein and DFT source core differ')
                source = ct['xyz'] if variant == 'contextual' else next(t['xyz'] for t in nm['tasks'] if t['case'] == case+'__full' and t['metal'] == metal)
                coords, mapping = transfer(prep, cp, xyz(verify(source)), metal)
                paired_rows.append(coords)
                bp = out/(row['case_id']+'_'+variant+'_'+metal+'.xyz'); write_xyz(bp, coords)
                transformations.append({'case': case, 'variant': variant, 'metal': metal, 'source_core_xyz': source,
                                        'source_parent': ct['parent'], 'preparation': row['preparation'],
                                        'bound_xyz': record(bp), 'water_H_mapping': mapping})
                old, oldm = source_task(row['endpoints'][metal])
                if oldm['model'] != parent['model'] or oldm['software'] != parent['software']:
                    raise InvalidArtifact('archived scorer model differs')
                for position in (('bound', 'detached') if variant == 'contextual' else ('bound',)):
                    task = copy.deepcopy(old)
                    for key in ('cache_key', 'ablation_group', 'original_endpoint', 'reference_task_id'):
                        task.pop(key, None)
                    task.update(task_id=f'{case}_{variant}_{metal}_{position}', position=position,
                                preparation_variant=variant, variant='primary', source_xyz=record(bp),
                                xyz=record(bp), require_isolated_metal=position == 'detached',
                                source_core_xyz=source, source_parent=ct['parent'])
                    if position == 'detached':
                        xp = out/(task['task_id']+'.xyz'); write_xyz(xp, geometry(task)); task['xyz'] = record(xp)
                    tasks.append(task)
            if variant == 'normalized' and [a[1:] for a in paired_rows[0]] != [a[1:] for a in paired_rows[1]]:
                raise InvalidArtifact('two-call normalized geometry differs between metals')
    dry = []
    for row in data['rows']:
        if row['case_id'].startswith('ALPHA_'):
            continue
        p = read_json(verify(row['preparation']))
        if p['explicit_waters']:
            raise InvalidArtifact('declared PQQ/GGR identity unexpectedly contains site waters')
        dry.append({'case_id': row['case_id'], 'target': row['target'], 'label': row['label'],
                    'MACE_score': row['MACE_score'], 'baseline_class': row['baseline_class'],
                    'preparation': row['preparation'], 'operation': 'exact_dry_identity'})
    if len(dry) != 30:
        raise InvalidArtifact('expected25 canonical PQQ,2PQQ crystals,3GGR dry identities')
    m.update(protocol_id=PROTOCOL, model=descriptor_model(verify(parent['software'])), tasks=tasks,
             source_inventory=record(inventory), core_collection=record(core_collection),
             normalized_collection=record(normalized_collection), ggr_comparison=record(ggr_comparison),
             transformations=transformations, dry_identities=dry, new_quantum_calls=0,
             model_calls=12, baseline_changed=False, calibrated_class=None)
    return omol.seal(out, m)


@cached_file_checks
def validate(manifest):
    from mace_omol_ablation_run import descriptor_model
    from mace_omol_intact import geometry
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL or m['stage'] != STAGE or m['model'] != descriptor_model(verify(m['software'])):
        raise InvalidArtifact('scorer protocol/model changed')
    for pin in [m['agreement'], m['source_inventory'], m['core_collection'], m['normalized_collection'], m['ggr_comparison'], *m['implementation'].values()]:
        verify(pin)
    expected = {(case, variant, metal, position) for case in ('1F6S', '6IP9') for variant in ('normalized', 'contextual')
                for metal in ('Ca', 'La') for position in (('bound', 'detached') if variant == 'contextual' else ('bound',))}
    found = set()
    for t in m['tasks']:
        found.add((t['case_id'].split('_')[1], t['preparation_variant'], t['metal'], t['position']))
        prep = read_json(verify(t['preparation'])); parent = read_json(verify(t['source_parent']))
        wanted, _ = transfer(prep, parent, xyz(verify(t['source_core_xyz'])), t['metal'])
        if xyz(verify(t['source_xyz'])) != wanted or xyz(verify(t['xyz'])) != geometry(t):
            raise InvalidArtifact('derived physical coordinates differ')
        if t['state'] != check_atoms(wanted, t['charge']) or t['spin_multiplicity'] != 1:
            raise InvalidArtifact('paired chemical state changed')
        payload = {k: v for k, v in t.items() if k != 'cache_key'}
        if t['cache_key'] != cache_key({'task': payload, 'model': m['model'], 'software': m['software'], 'implementation': m['implementation']}):
            raise InvalidArtifact('cache identity changed')
    if found != expected or len(m['tasks']) != 12:
        raise InvalidArtifact('finite12-call comparison changed')
    for case in ('1F6S', '6IP9'):
        pair = [xyz(verify(t['source_xyz'])) for t in m['tasks']
                if t['case_id'] == 'ALPHA_'+case and t['preparation_variant'] == 'normalized']
        if len(pair) != 2 or [a[1:] for a in pair[0]] != [a[1:] for a in pair[1]]:
            raise InvalidArtifact('common-geometry cancellation is unsupported')
    return {'status': 'pass', 'tasks': 12, 'new_quantum_calls': 0, 'manifest': record(manifest)}


@cached_file_checks
def collect(manifest):
    validate(manifest); m = read_json(manifest); mp = Path(manifest); rows = {}; cases = []; comparisons = []
    inv = read_json(verify(m['source_inventory'])); core = read_json(verify(m['core_collection']))
    normal = read_json(verify(m['normalized_collection'])); ggr = read_json(verify(m['ggr_comparison']))
    for t in m['tasks']:
        valid = []
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r = accepted_attempt(a, t, mp)
            if r is not None:
                valid.append((a, r))
        if not valid:
            rows[t['task_id']] = {'status': 'unavailable'}
        else:
            a, r = valid[-1]
            rows[t['task_id']] = {'status': 'computed', 'energy_eV': r['energy_eV'], 'result': record(a/'result.json'),
                                  'evaluation_seconds': r['evaluation_seconds'], 'worker_wall_seconds': r['wall_seconds']}
    for case in ('1F6S', '6IP9'):
        source = next(r for r in inv['rows'] if r['case_id'] == 'ALPHA_'+case)
        dft = next(r for r in core['cases'] if r['case'] == case)
        normalized = next(r for r in normal['rows'] if r['case'] == case and r['state_id'] == case+'__full')
        entry = {'case': case, 'MACE_archived_R_model_kcal': source['MACE_score'],
                 'DFT_original_R_hartree': dft['original_R_hartree'], 'DFT_normalized_R_hartree': normalized['R_hartree'],
                 'DFT_contextual_R_hartree': dft['prepared_R_hartree']}
        for variant in ('normalized', 'contextual'):
            tasks = [t for t in m['tasks'] if t['case_id'] == 'ALPHA_'+case and t['preparation_variant'] == variant]
            if all(rows[t['task_id']]['status'] == 'computed' for t in tasks):
                energies = {(t['metal'], t['position']): rows[t['task_id']]['energy_eV'] for t in tasks}
                entry['MACE_'+variant] = contrasts(energies, source['disconnected'], common_geometry=variant == 'normalized')
            else:
                entry['MACE_'+variant] = None
        cases.append(entry)
        for g in ggr['controls']:
            gr = next(r for r in inv['rows'] if r['case_id'] == 'GGR_'+g['source'])
            c = {'alpha': case, 'GGR': g['source'], 'evidence': 'one_condition_qualified_biological_comparison',
                 'DFT_original': (dft['original_R_hartree'] - g['R_hartree'])*HA_TO_KCAL,
                 'DFT_normalized': (normalized['R_hartree'] - g['R_hartree'])*HA_TO_KCAL,
                 'DFT_contextual': (dft['prepared_R_hartree'] - g['R_hartree'])*HA_TO_KCAL,
                 'MACE_archived': source['MACE_score'] - gr['MACE_score']}
            for variant in ('normalized', 'contextual'):
                value = entry['MACE_'+variant]
                c['MACE_'+variant] = value['interaction_R_model_kcal'] - gr['MACE_score'] if value is not None else None
            comparisons.append(c)
    return {'status': 'complete' if all(r['status'] == 'computed' for r in rows.values()) else 'incomplete',
            'manifest': record(manifest), 'rows': rows, 'cases': cases, 'comparisons': comparisons,
            'dry_identities': m['dry_identities'], 'baseline_changed': False, 'new_quantum_calls': 0,
            'interpretation': 'positive alpha-minus-GGR is expected; model scales separate, no calibrated affinity or occupancy',
            'independent_biological_comparisons': 1}


if __name__ == '__main__':
    import json
    p = argparse.ArgumentParser(description=__doc__); s = p.add_subparsers(dest='op', required=True)
    q = s.add_parser('prepare')
    for name in ('inventory', 'core-collection', 'normalized-collection', 'ggr-comparison', 'agreement', 'output'):
        q.add_argument('--'+name, required=True)
    q = s.add_parser('validate'); q.add_argument('--manifest', required=True)
    q = s.add_parser('collect'); q.add_argument('--manifest', required=True); q.add_argument('--output', required=True)
    a = vars(p.parse_args()); op = a.pop('op'); dest = a.pop('output') if op == 'collect' else None
    r = globals()[op](**a)
    if dest: write_new(dest, r)
    print(json.dumps(r if op != 'collect' else {k: v for k, v in r.items() if k not in ('rows', 'dry_identities')}, indent=2))
