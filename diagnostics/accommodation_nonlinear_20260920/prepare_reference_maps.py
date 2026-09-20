"""Geometry-only mapping of all archived PQQ references; no energy engine."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from accommodation_torsion_profiles import role_mode
from coordination_preparation_context import geometry
from mace_hybrid import check_atoms
from mace_site_kinematics import Kinematics
from second_shell_context import parent_state


def role_resname(selector):
    if isinstance(selector, dict):
        return selector['resname']
    match = re.fullmatch(r'[^:]+:([A-Z]+)-?\d+[A-Za-z]?', selector)
    if not match:
        raise InvalidArtifact('unrecognized role selector')
    return match.group(1)


def parent_output(pin, parent_pin):
    """Canonical manifests store output names relative to their own directory."""
    path = Path(pin['path'])
    if not path.is_absolute():
        path = verify(parent_pin).parent / path
    resolved = {**pin, 'path': str(path.resolve())}
    verify(resolved)
    return resolved


def prepare(root, output, agreement, parent_plan):
    start = time.monotonic()
    root = Path(root).resolve()
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    compact_path = root / 'workspaces/compact_solvation_20260920/full_v1/manifest.json'
    compact = read_json(compact_path)
    inventory = read_json(verify(compact['inventory']))
    static_path = root / 'workspaces/second_shell_20260919/prepared_v2/manifest.json'
    static = read_json(static_path)
    config = read_json(verify(static['configuration']))
    selected = [c for c in inventory['cases']
                if c['role'] == 'calibration' or c['case_id'] in ('1H4I', '4MAE', '1KB0')]
    if len(selected) != 28 or sum(c['role'] == 'calibration' for c in selected) != 25:
        raise InvalidArtifact('reference denominator changed')
    cases, results = [], []
    for item in selected:
        tick = time.monotonic()
        case = item['case_id']
        result = {'case_id': case, 'status': 'unsupported', 'endpoint_checks': {}}
        try:
            reps = item['representations']
            parent_pin = reps['core']['preparation']
            parent = read_json(verify(parent_pin))
            prep_pin = reps['context']['preparation']
            prep = read_json(verify(prep_pin))
            source = {'case': case, 'parent': parent_pin, 'endpoints': {
                z: {'xyz': parent_output(parent['outputs'][z + '_xyz'], parent_pin),
                    'input': parent_output(parent['outputs'][z + '_input'], parent_pin),
                    'charge': parent['charge_ledger'][z + '_total'], 'multiplicity': 1}
                for z in ('Ca', 'La')}}
            state = parent_state(source, config['topology'], require_endpoint_receipts=False)
            requested = parent['fixed_core']['requested_roles']
            if role_resname(requested['anchor_glutamate']) != 'GLU':
                raise InvalidArtifact('anchor role is not GLU')
            roles = ['anchor_glutamate']
            extra = role_resname(requested['extra_acidic_ligand_homolog'])
            if extra == 'ASP':
                roles.append('extra_acidic_ligand_homolog')
            modes = {r: role_mode(parent, r) for r in roles}
            origins, maps, mapping_data, context_rows = {}, {}, {}, {}
            for z in ('Ca', 'La'):
                endpoint = reps['context']['endpoints'][z]
                task = next(t for t in compact['all_tasks'] if
                            (t['case_id'], t['representation'], t['metal'], t['medium']) ==
                            (case, 'context', z, 'vacuum'))
                rows = xyz(verify(endpoint['xyz']))
                if xyz(verify(task['xyz'])) != rows:
                    raise InvalidArtifact('archived compact context differs')
                if xyz(verify(reps['core']['endpoints'][z]['xyz'])) != state['original'][z]:
                    raise InvalidArtifact('archived original core differs')
                charge = source['endpoints'][z]['charge'] + prep['added_formal_charge']
                if charge != endpoint['charge'] or charge != task['charge'] or endpoint['multiplicity'] != 1:
                    raise InvalidArtifact('context charge or multiplicity changed')
                atom_state = check_atoms(rows, charge)
                gp = geometry(state, prep, rows, state['original'][z])
                checks = {'state': atom_state, 'representations': {}}
                for rep, target in (('core', state['original'][z]), ('context', rows)):
                    kin = Kinematics(gp[rep])
                    names = [m['id'] for m in kin.modes]
                    active = [names.index(mid) for mid in modes.values()]
                    zero = np.zeros(len(names))
                    observed = kin.evaluate(zero)[1]
                    expected = np.array([a[1:] for a in target])
                    max_error = float(np.max(np.abs(observed - expected)))
                    if max_error > 1e-10:
                        raise InvalidArtifact('q0 coordinate replay differs')
                    moving_groups = {}
                    for role, mid in modes.items():
                        matching = [f for f in parent['qm_fragments'] if f.get('role') == role]
                        expected_resname = 'GLU' if role == 'anchor_glutamate' else 'ASP'
                        if len(matching) != 1 or role_resname(matching[0]['id']) != expected_resname:
                            raise InvalidArtifact('active role differs from retained chemistry')
                        mode = kin.modes[names.index(mid)]
                        meta = kin.data['source_atom_metadata']
                        atom_names = {meta[i]['atom'] for i in mode['moving_indices'] if meta[i]}
                        required = {'OE1', 'OE2'} if role == 'anchor_glutamate' else {'OD1', 'OD2'}
                        if not required <= atom_names:
                            raise InvalidArtifact('terminal carboxylate mapping incomplete')
                        moving_groups[role] = sorted(atom_names)
                    probe = zero.copy()
                    probe[active] = .001
                    numerical = kin.check(probe)
                    if not numerical['pass']:
                        raise InvalidArtifact('active source/cap mapping check failed')
                    checks['representations'][rep] = {
                        'q0_bitwise_equal': bool(np.array_equal(observed, expected)),
                        'q0_max_absolute_error_A': max_error, 'active_mode_ids': list(modes.values()),
                        'active_mode_indices': active, 'active_terminal_atom_names': moving_groups,
                        'inactive_coordinates': 'exactly zero', 'mapping_checks': numerical,
                        'source_links': sum(a[1] == 'source' for a in kin.links),
                        'cap_links': sum(a[1] == 'cap' for a in kin.links)}
                dest = out / 'sources' / case / (z + '_mapping.json')
                write_new(dest, gp)
                maps[z] = record(dest)
                origins[z] = {'xyz': endpoint['xyz'], 'charge': charge, 'multiplicity': 1}
                mapping_data[z], context_rows[z] = gp, rows
                result['endpoint_checks'][z] = checks
            ca, la = context_rows['Ca'], context_rows['La']
            if ca[1:] != la[1:] or ca[0][1:] != la[0][1:] or ca[0][0] != 'Ca' or la[0][0] != 'La':
                raise InvalidArtifact('paired context coordinates/elements differ')
            if origins['La']['charge'] - origins['Ca']['charge'] != 1:
                raise InvalidArtifact('paired endpoint charge difference changed')
            if mapping_data['Ca'] != mapping_data['La']:
                raise InvalidArtifact('physical coordinate mapping differs between endpoints')
            cases.append({'case_id': case, 'source': source, 'preparation': prep_pin,
                          'modes': modes, 'maps': maps, 'origins': origins,
                          'label_scope': item['label_scope'], 'evidence_stratum': item['evidence_stratum'],
                          'biological_group': item['biological_group'],
                          'role_eligibility': {'anchor_glutamate': 'GLU_chi3_active',
                              'extra_acidic_ligand_homolog': 'ASP_chi2_active' if extra == 'ASP' else extra + '_fixed'},
                          'expected_class_for_later_report_only': item['expected_class']})
            result.update(status='supported', modes=modes, extra_role_resname=extra,
                          paired_coordinates_exact=True, paired_mappings_identical=True,
                          context_atoms=len(ca), original_core_atoms=len(state['original']['Ca']))
        except Exception as exc:
            result.update(reason=str(exc), exception_type=type(exc).__name__)
        result['wall_seconds'] = time.monotonic() - tick
        results.append(result)
        write_new(out / 'checks' / (case + '.json'), result)
        print(json.dumps({'case_id': case, 'status': result['status'], 'reason': result.get('reason')}), flush=True)
    design = {'protocol_id': 'PQQ_common_terminal_carboxylate_reference_mapping_v1',
              'agreement': record(agreement), 'parent_plan': record(parent_plan), 'cases': cases,
              'declared_case_ids': [c['case_id'] for c in selected], 'case_denominator': 28,
              'compact_source': record(compact_path), 'inventory': compact['inventory'],
              'static_source': record(static_path), 'model': inventory['model'],
              'software': inventory['software'], 'orca': compact['orca'], 'topology': config['topology'],
              'expected_new_calls': {'MACE': 0, 'GFN2': 0, 'DFT': 0},
              'execution_authorized_by_this_preparation': False, 'baseline_changed': False,
              'first_plan_settings': 'unchanged; see pinned parent_plan',
              'implementation': record(__file__)}
    write_new(out / 'design.json', design)
    result = {'design': record(out / 'design.json'), 'cases': results, 'denominator': 28,
              'status_counts': dict(Counter(r['status'] for r in results)),
              'new_molecular_calls': 0, 'wall_seconds': time.monotonic() - start,
              'implementation': record(__file__)}
    write_new(out / 'result.json', result)
    return {k: v for k, v in result.items() if k != 'cases'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'output', 'agreement', 'parent-plan'):
        parser.add_argument('--' + name, required=True)
    print(json.dumps(prepare(**vars(parser.parse_args())), indent=2))
