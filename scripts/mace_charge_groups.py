"""Finite grouped-charge POLAR pilot through the existing MACE execution runner."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import check_atoms, write_xyz, rotation, accepted_attempt, EV_TO_KCAL
from mace_charge_group_prepare import CASES, POLICY, groups, selected_residues
from mace_group_constraints import ADAPTER

SCHEMA = 'alquemia.mace_charge_groups.v1'
GB_SCHEMA = 'alquemia.mace_charge_groups_gb.v1'
PROTOCOL = 'intact_POLAR_medium_chemical_group_charge_frozen_OBC2_v1'
TOL = {'energy_kcal_mol': .01, 'force_max_eV_A': .001, 'total_charge_e': 1e-5,
       'grouping_kcal_mol': 2., 'direction_margin_kcal_mol': .02}


def kernel_parent_gate(manifest_data):
    from mace_global_benchmark import numerical_parent_gate
    return numerical_parent_gate(read_json(verify(manifest_data['parent_manifest'])))


def specs():
    return ([('GGR_1GLG', m, 'one_group') for m in ('Ca', 'La')]
            + [(c, m, 'primary') for c in sorted(CASES) for m in ('Ca', 'La')]
            + [('GGR_1GLG', m, 'rotate') for m in ('Ca', 'La')]
            + [(c, m, 'connected') for c in sorted(CASES) if c.startswith('GGR') for m in ('Ca', 'La')])


@cached_file_checks
def prepared(preparation):
    p = read_json(preparation)
    cfg = read_json(verify(p['config']))
    if p['policy_id'] != POLICY or p['status'] != 'prepared' or set(p['cases']) != CASES:
        raise InvalidArtifact('seven complete fixed group preparations required')
    verify(p['implementation'])
    cases = {}
    for name, pin in p['cases'].items():
        r = read_json(verify(pin)); c = read_json(verify(r['sources']['physical']))
        if r['sources'] != cfg['cases'][name] or c['case_id'] != name:
            raise InvalidArtifact('group physical source differs')
        bonds = sorted(tuple(sorted((b['atom_a_id'], b['atom_b_id']))) for b in c['preparation_details']['bonds'])
        selected = selected_residues(c, r['sources']['primary_core'], r['sources']['core_format'])
        if groups(c, r['protein_ledger'], bonds, selected) != r['primary']:
            raise InvalidArtifact('primary chemical group assignment changed')
        if name.startswith('GGR'):
            selected |= selected_residues(c, r['sources']['alternate_core'], r['sources']['core_format'])
            if groups(c, r['protein_ledger'], bonds, selected) != r['connected']:
                raise InvalidArtifact('connected chemical group assignment changed')
        cases[name] = (r, c)
    return p, cases, cfg


def task(name, metal, variant, group_row, physical, output=None):
    e = physical['endpoints'][metal]
    atoms = xyz(verify(e['xyz']))
    expected = [(metal if a['kind'] == 'selected_metal' else a['element'], *a['xyz_A']) for a in physical['physical_atoms']]
    if atoms != expected:
        raise InvalidArtifact('frozen source endpoint differs')
    g = group_row['connected' if variant == 'connected' else 'primary']
    ids, charges = g['atom_group_indices'], g['endpoint_group_charges_e'][metal]
    if variant == 'one_group':
        ids, charges = [0] * len(atoms), [e['charge']]
    matrix = rotation() if variant == 'rotate' else np.eye(3)
    if variant == 'rotate':
        center = np.array(next(a[1:] for a in atoms if a[0] == metal))
        atoms = [(a[0], *((np.array(a[1:]) - center) @ matrix.T + center)) for a in atoms]
    tid = name + '_' + metal + '_' + variant
    result = dict(task_id=tid, case_id=name, kind='full', metal=metal, variant=variant,
                  charge=e['charge'], spin_multiplicity=1, state=check_atoms(atoms, e['charge']),
                  preparation=group_row['sources']['physical'], source_xyz=e['xyz'],
                  atom_group_indices=ids, group_charges_e=charges,
                  charge_group_adapter=ADAPTER, rotation_matrix=matrix.tolist())
    if output is not None:
        path = output / (tid + '.xyz'); write_xyz(path, atoms); result['xyz'] = record(path)
    return result, atoms


def prepare(preparation, parent_manifest, output):
    from mace_global_benchmark import snapshot
    p, cases, cfg = prepared(preparation)
    parent = read_json(parent_manifest)
    if parent['checkpoint_label'] != 'medium' or parent['schema_version'] != 'alquemia.mace_global_benchmark.v1':
        raise InvalidArtifact('pinned medium whole-protein parent required')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    pins = snapshot(out, (*parent['implementation'], 'mace_charge_groups.py', 'mace_group_constraints.py',
                          'mace_charge_group_prepare.py', 'mace_file_checks.py'))
    model = copy.deepcopy(parent['model'])
    inv = read_json(verify(read_json(verify(parent['software']))['backend_source_inventory']))
    sources = [ref for ref in inv['files'] if ref['path'].endswith('/mace/modules/extensions.py')]
    if len(sources) != 1:
        raise InvalidArtifact('pinned native POLAR source unavailable')
    model.update(charge_group_adapter=ADAPTER, native_polar_source=sources[0], chemical_group_policy=POLICY)
    tasks = [task(c, m, v, *cases[c], out)[0] for c, m, v in specs()]
    result = dict(schema_version=SCHEMA, protocol_id=PROTOCOL, stage='MACE', preparation=record(preparation),
                  parent_manifest=record(parent_manifest), model=model, tasks=tasks, implementation=pins,
                  software=parent['software'], agreement=cfg['agreement'], tolerances=TOL, compute_budget=None,
                  new_DFT_calls=0, calibrated_class=None, reference=None, baseline_changed=False)
    for t in tasks:
        t['cache_key'] = cache_key({'task': t, 'model': model, 'software': result['software'], 'implementation': pins})
    write_new(out / 'manifest.json', result)
    return validate(out / 'manifest.json')


@cached_file_checks
def validate(manifest):
    m = read_json(manifest); is_gb = m['schema_version'] == GB_SCHEMA
    if m['schema_version'] not in (SCHEMA, GB_SCHEMA) or m['protocol_id'] != PROTOCOL or m['tolerances'] != TOL:
        raise InvalidArtifact('unknown group model or changed criteria')
    p, cases, cfg = prepared(verify(m['preparation']))
    if m['agreement'] != cfg['agreement'] or len(m['tasks']) != 24:
        raise InvalidArtifact('group plan/task count changed')
    parent = read_json(verify(m['parent_manifest']))
    if kernel_parent_gate(m)['status'] != 'pass':
        raise InvalidArtifact('existing analytic electrostatic kernel gate failed')
    expected_model = copy.deepcopy(parent['model'])
    if is_gb:
        from mace_gb import MODEL
        expected_model = MODEL
        source = collect(verify(m['source_mace_manifest']))
        if source['status'] != 'complete':
            raise InvalidArtifact('missing learned density; no fallback permitted')
    else:
        expected_model.update(charge_group_adapter=ADAPTER, native_polar_source=m['model']['native_polar_source'], chemical_group_policy=POLICY)
        verify(m['model']['native_polar_source'])
        if m['software'] != parent['software']:
            raise InvalidArtifact('MACE software differs')
    if expected_model != m['model']:
        raise InvalidArtifact('scientific model changed')
    software = read_json(verify(m['software']))
    for pin in [m['agreement'], *m['implementation'].values(), software['python'], software['requirements']]:
        verify(pin)
    for key in ('checkpoint', 'backend_source_inventory'):
        if key in software:
            verify(software[key])
    for t, (c, metal, v) in zip(m['tasks'], specs()):
        expected, atoms = task(c, metal, v, *cases[c])
        if is_gb:
            expected.pop('charge_group_adapter')
            row = source['rows'][t['task_id']]
            expected.update(platform='CUDA', solver='native', solvent_dielectric=78.5,
                            source_density=row['density_coefficients'], source_vacuum_energy_eV=row['energy_eV'])
        if {k: value for k, value in t.items() if k not in ('xyz', 'cache_key')} != expected:
            raise InvalidArtifact('task settings or chemical groups differ')
        actual = xyz(verify(t['xyz']))
        if ([a[0] for a in actual] != [a[0] for a in atoms] or
                not np.allclose([a[1:] for a in actual], [a[1:] for a in atoms], atol=1e-12, rtol=0)):
            raise InvalidArtifact('paired source geometry or rigid transform changed')
        payload = {k: v for k, v in t.items() if k != 'cache_key'}
        if t['cache_key'] != cache_key({'task': payload, 'model': m['model'], 'software': m['software'], 'implementation': m['implementation']}):
            raise InvalidArtifact('group task cache mismatch')
    return {'status': 'pass', 'tasks': 24, 'manifest': record(manifest)}


@cached_file_checks
def collect(manifest):
    validate(manifest)
    mp = Path(manifest).resolve(); m = read_json(mp); rows, attempts = {}, []
    for t in m['tasks']:
        accepted = []
        for a in sorted((mp.parent / 'execution' / t['task_id']).glob('attempt_*')):
            r = accepted_attempt(a, t, mp)
            attempts.append({'task_id': t['task_id'], 'directory': str(a), 'accepted': r is not None})
            if r is not None:
                accepted.append(r)
        rows[t['task_id']] = accepted[-1] if accepted else {'status': 'unavailable', 'energy_eV': None}
    return {'status': 'complete' if all(r['status'] == 'computed' for r in rows.values()) else 'incomplete',
            'manifest': record(mp), 'rows': rows, 'attempts': attempts}


def prepare_gb(manifest, solver_validation, output):
    from mace_global_benchmark import snapshot
    from mace_gb import MODEL
    source = collect(manifest); m = read_json(manifest)
    ref = read_json(solver_validation); rm = read_json(verify(ref['manifest']))
    if source['status'] != 'complete' or not ref['numerical_checks_pass']:
        raise InvalidArtifact('completed learned endpoints and qualified solvent backend required')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    pins = snapshot(out, m['implementation'])
    result = copy.deepcopy(m)
    result.update(schema_version=GB_SCHEMA, stage='OBC2', model=MODEL, software=rm['software'],
                  source_mace_manifest=record(manifest), solver_validation=record(solver_validation), implementation=pins)
    for t in result['tasks']:
        t.pop('cache_key'); t.pop('charge_group_adapter')
        r = source['rows'][t['task_id']]
        t.update(platform='CUDA', solver='native', solvent_dielectric=78.5,
                 source_density=r['density_coefficients'], source_vacuum_energy_eV=r['energy_eV'])
        t['cache_key'] = cache_key({'task': t, 'model': MODEL, 'software': result['software'], 'implementation': pins})
    write_new(out / 'manifest.json', result)
    return validate(out / 'manifest.json')


def report(manifest, native_reference, output, solvent_manifest=None):
    """Preserve raw directions separately from numerical and grouping admission."""
    import time
    start = time.monotonic()
    molecular = collect(manifest); m = read_json(manifest)
    old = read_json(native_reference); old_manifest = read_json(verify(old['manifest']))
    if old['manifest'] != m['parent_manifest'] or old['status'] != 'complete':
        raise InvalidArtifact('matching actual native comparison required')
    solvent = collect(solvent_manifest) if solvent_manifest else None
    if solvent_manifest and read_json(solvent_manifest)['source_mace_manifest'] != record(manifest):
        raise InvalidArtifact('solvent belongs to a different learned model')
    rows = molecular['rows']; checks = []; scores = {}; grouping = []
    def check(name, error, tol=.01):
        checks.append({'name': name, 'error': error, 'tolerance': tol,
                       'pass': error is not None and abs(error) <= tol})
    def computed(row):
        return row.get('status') == 'computed'
    for metal in ('Ca', 'La'):
        tid = 'GGR_1GLG_' + metal + '_primary'
        original = old['rows'][tid]
        ot = next(t for t in old_manifest['tasks'] if t['task_id'] == tid)
        archive = verify(old['manifest']).parent / 'execution' / tid
        if original not in [accepted_attempt(a, ot, verify(old['manifest'])) for a in archive.glob('attempt_*')]:
            raise InvalidArtifact('native reference lacks a matching actual receipt')
        new = rows['GGR_1GLG_' + metal + '_one_group']
        check(metal + '_one_group_native_energy',
              (new['energy_eV'] - original['energy_eV']) * EV_TO_KCAL if computed(new) else None)
        primary, rotated = [rows['GGR_1GLG_' + metal + '_' + v] for v in ('primary', 'rotate')]
        check(metal + '_grouped_rotation_energy',
              (rotated['energy_eV'] - primary['energy_eV']) * EV_TO_KCAL
              if computed(primary) and computed(rotated) else None)
        if solvent:
            a, b = [solvent['rows']['GGR_1GLG_' + metal + '_' + v] for v in ('primary', 'rotate')]
            check(metal + '_solvent_rotation_energy',
                  b['GB_reaction_kcal_mol'] - a['GB_reaction_kcal_mol'] if computed(a) and computed(b) else None)
    for tid, row in rows.items():
        error = None
        if computed(row):
            trace = read_json(verify(row['charge_groups']))
            error = max(s['maximum_channel_charge_error_e'] for s in trace['stages'])
        check(tid + '_group_charge_closure_e', error, 1e-5)
    for name in sorted(CASES):
        variants = ('primary', 'connected') if name.startswith('GGR') else ('primary',)
        scores[name] = {}
        for variant in variants:
            endpoints = [rows[name + '_' + metal + '_' + variant] for metal in ('Ca', 'La')]
            raw = (endpoints[0]['energy_eV'] - endpoints[1]['energy_eV']) * EV_TO_KCAL if all(map(computed, endpoints)) else None
            gb = None
            if solvent:
                endpoints_gb = [solvent['rows'][name + '_' + metal + '_' + variant] for metal in ('Ca', 'La')]
                if all(map(computed, endpoints_gb)):
                    gb = endpoints_gb[0]['GB_reaction_kcal_mol'] - endpoints_gb[1]['GB_reaction_kcal_mol']
            scores[name][variant] = {'R_vacuum_model_kcal': raw, 'GB_Ca_minus_La_kcal': gb,
                                     'R_candidate_model_kcal': raw + gb if raw is not None and gb is not None else None,
                                     'calibrated_class': None, 'qualified_score': None}
        if name.startswith('GGR'):
            a, b = [scores[name][v]['R_candidate_model_kcal'] for v in ('primary', 'connected')]
            delta = b - a if a is not None and b is not None else None
            grouping.append({'case': name, 'connected_minus_primary_model_kcal': delta,
                             'tolerance': 2., 'pass': delta is not None and abs(delta) <= 2.})
    pairs = [('PQQ_4MAE', 'PQQ_1H4I', 'PQQ_functional_class')]
    pairs += [(a, g, 'qualified_alpha_vs_direct_GGR_direction')
              for a in ('ALPHA_1F6S', 'ALPHA_6IP9') for g in ('GGR_1GLG', 'GGR_2FW0', 'GGR_2FVY')]
    contrasts = []
    for higher, lower, stratum in pairs:
        a, b = [scores[name]['primary']['R_candidate_model_kcal'] for name in (higher, lower)]
        delta = a - b if a is not None and b is not None else None
        contrasts.append({'higher_expected': higher, 'lower_expected': lower, 'evidence_stratum': stratum,
                          'margin_model_kcal': delta, 'pass': delta is not None and delta > .02})
    complete = molecular['status'] == 'complete' and solvent is not None and solvent['status'] == 'complete'
    numeric = complete and all(x['pass'] for x in checks)
    admitted = numeric and all(x['pass'] for x in grouping)
    result = {'status': 'complete' if complete else 'incomplete', 'protocol_id': PROTOCOL,
              'MACE': molecular, 'solvent': solvent, 'native_reference': record(native_reference),
              'scores': scores, 'checks': checks, 'grouping_checks': grouping,
              'numerical_checks_pass': numeric, 'representation_checks_pass': all(x['pass'] for x in grouping),
              'contrasts': contrasts, 'raw_pass_count': sum(x['pass'] for x in contrasts), 'direction_denominator': 7,
              'qualified_pass_count': sum(x['pass'] for x in contrasts) if admitted else 0,
              'biological_groups': 4, 'prospectively_blind': False, 'baseline_changed': False,
              'absolute_reference': None, 'calibrated_class': None, 'affinity_free_energy': None,
              'combined_solvent_gradient': None, 'response_status': 'response_model_not_validated',
              'wall_seconds': time.monotonic() - start}
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    reporter = out / 'report_source.py'; reporter.write_bytes(Path(__file__).read_bytes())
    result['implementation'] = record(reporter)
    write_new(out / 'result.json', result)
    return {key: result[key] for key in ('status', 'numerical_checks_pass', 'representation_checks_pass',
                                       'raw_pass_count', 'direction_denominator', 'qualified_pass_count', 'contrasts')}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for op, names in [('prepare', ('preparation', 'parent-manifest', 'output')),
                      ('validate', ('manifest',)), ('collect', ('manifest',)),
                      ('prepare-gb', ('manifest', 'solver-validation', 'output')),
                      ('report', ('manifest', 'native-reference', 'output'))]:
        p = sub.add_parser(op)
        for name in names:
            p.add_argument('--' + name, required=True)
        if op == 'report':
            p.add_argument('--solvent-manifest')
    args = vars(parser.parse_args()); op = args.pop('command').replace('-', '_')
    print(json.dumps(globals()[op](**args), indent=2))
