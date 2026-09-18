"""Finite intact-protein grouped-vacuum descriptor transfer using the existing runner."""
import argparse
import copy
import json
from pathlib import Path
import time
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_charge_group_prepare import protein_ledger, selected_residues
from mace_charge_groups import kernel_parent_gate, TOL
from mace_group_constraints import ADAPTER
from mace_site_groups import POLICY, build
from mace_hybrid import check_atoms, rotation, write_xyz, accepted_attempt, EV_TO_KCAL

SCHEMA = 'alquemia.mace_group_transfer.v1'
PROTOCOL = 'intact_POLAR_medium_typed31_group_vacuum_compatibility_v1'
COUNTS = {**{k: 1 for k in ('PQQ_1H4I', 'PQQ_4MAE', 'GGR_1GLG', 'GGR_2FW0', 'GGR_2FVY',
                           'ALPHA_1F6S', 'ALPHA_6IP9')},
          'A0A7': 6, 'HEW5': 8, 'RTX': 8, 'PARV_4CPV': 2, 'AEQ_1SL8': 3}


def family_cases(cfg):
    groups = {k: [] for k in COUNTS}
    for name, c in cfg['cases'].items():
        groups[c['group']].append(name)
    if {k: len(v) for k, v in groups.items()} != COUNTS or cfg['cutoff_A'] != 3.1:
        raise InvalidArtifact('fixed physical inventory or typed radius changed')
    return {k: sorted(v, key=lambda n: (cfg['cases'][n]['site_order'], n)) for k, v in groups.items()}


def paired_source(p):
    for metal in ('Ca', 'La'):
        want = [(metal if a['kind'] == 'selected_metal' else a['element'], *a['xyz_A']) for a in p['physical_atoms']]
        if xyz(verify(p['endpoints'][metal]['xyz'])) != want:
            raise InvalidArtifact('source endpoint atom/coordinate mismatch')


def prepare_inputs(config, output):
    start, cpu = time.monotonic(), time.process_time()
    cfg = read_json(config); families = family_cases(cfg)
    for k in ('agreement', 'typed_donor_source', 'parent_manifest', 'old_report', 'khoury', 'multisite'):
        verify(cfg[k])
    if record(Path(__file__).with_name('coordination_policy.py'))['sha256'] != cfg['typed_donor_source']['sha256']:
        raise InvalidArtifact('typed donor implementation changed')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    source = out / 'implementation'; source.mkdir()
    for name in ('mace_group_transfer.py', 'mace_site_groups.py'):
        (source / name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows = {}; ledgers = {}; reference = {}
    for family, names in families.items():
        for name in names:
            c = cfg['cases'][name]; p = read_json(verify(c['preparation'])); paired_source(p)
            protein = sorted((a['id'], a['element'], a.get('resname'), tuple(a['xyz_A'])) for a in p['physical_atoms']
                             if a['kind'] in ('protein_source', 'terminal_completion'))
            bonds = sorted(tuple(sorted((b['atom_a_id'], b['atom_b_id']))) for b in p['preparation_details']['bonds'])
            if family not in ledgers:
                ledger, expected_bonds = protein_ledger(p)
                ledgers[family] = (protein, bonds, ledger)
                if bonds != expected_bonds:
                    raise InvalidArtifact('template replay bond mismatch')
            else:
                prev, prev_bonds, ledger = ledgers[family]
                if protein != prev or bonds != prev_bonds:
                    raise InvalidArtifact('protein atoms/bonds differ across site substitutions')
            primary = build(p, ledger, bonds)
            common = {k: primary[k] for k in ('atom_ids', 'all_Ca_atoms', 'atom_group_indices', 'group_names', 'group_units')}
            common['Ca_group_charges'] = primary['endpoint_group_charges_e']['Ca']
            if family in reference and common != reference[family]:
                raise InvalidArtifact('all-Ca geometry or group partition depends on selected site')
            reference[family] = common
            row = {'case': c, 'protein_ledger': ledger, 'primary': primary}
            if c['alternate_core']:
                extra = selected_residues(p, c['alternate_core'], 'physical_mapping')
                row['connected'] = build(p, ledger, bonds, extra)
            path = out / (name + '.json'); write_new(path, row); rows[name] = record(path)
    result = {'policy_id': POLICY, 'status': 'prepared', 'config': record(config), 'cases': rows,
              'implementation': record(source / 'mace_group_transfer.py'), 'group_implementation': record(source / 'mace_site_groups.py'),
              'all_Ca_source_group_equivalence_pass': True, 'model_calls': 0,
              'wall_seconds': time.monotonic() - start, 'CPU_seconds': time.process_time() - cpu}
    write_new(out / 'preparation.json', result)
    return {'status': 'prepared', 'cases': len(rows), 'preparation': record(out / 'preparation.json'),
            'wall_seconds': result['wall_seconds'], 'CPU_seconds': result['CPU_seconds']}


@cached_file_checks
def prepared(path):
    p = read_json(path); cfg = read_json(verify(p['config'])); family_cases(cfg)
    verify(p['implementation']); verify(p['group_implementation']); verify(cfg['typed_donor_source'])
    if p['policy_id'] != POLICY or set(p['cases']) != set(cfg['cases']) or p['status'] != 'prepared':
        raise InvalidArtifact('group preparation inventory mismatch')
    cases = {}
    for name, pin in p['cases'].items():
        row = read_json(verify(pin)); src = read_json(verify(row['case']['preparation'])); paired_source(src)
        if row['case'] != cfg['cases'][name]:
            raise InvalidArtifact('source case changed')
        bonds = sorted(tuple(sorted((b['atom_a_id'], b['atom_b_id']))) for b in src['preparation_details']['bonds'])
        for variant in ('primary', 'connected') if row['case']['alternate_core'] else ('primary',):
            extra = selected_residues(src, row['case']['alternate_core'], 'physical_mapping') if variant == 'connected' else ()
            expected = json.loads(json.dumps(build(src, row['protein_ledger'], bonds, extra)))
            # NumPy versions differ by ~4e-16 A in norm roundoff. Membership,
            # atom/group identities and charges must still match exactly.
            stored = copy.deepcopy(row[variant])
            if len(stored['contacts']) != len(expected['contacts']):
                raise InvalidArtifact('typed contact inventory differs')
            for a, b in zip(stored['contacts'], expected['contacts']):
                if abs(a['distance_A'] - b['distance_A']) > 1e-12:
                    raise InvalidArtifact('typed contact source distance changed')
                a['distance_A'] = b['distance_A']
            if stored != expected:
                raise InvalidArtifact('source-derived group partition changed')
        cases[name] = row
    for names in family_cases(cfg).values():
        first = cases[names[0]]['primary']
        for name in names[1:]:
            row = cases[name]['primary']
            if any(row[k] != first[k] for k in ('atom_ids', 'all_Ca_atoms', 'atom_group_indices', 'group_names', 'group_units')) or row['endpoint_group_charges_e']['Ca'] != first['endpoint_group_charges_e']['Ca']:
                raise InvalidArtifact('shared all-Ca reference no longer matches all sites')
    return p, cfg, cases


def specs(cfg):
    families = family_cases(cfg)
    tasks = []
    for names in families.values():
        tasks.append((names[0], 'Ca', 'primary'))
        tasks += [(n, 'La', 'primary') for n in names]
    tasks += [(c, m, 'connected') for c in ('GGR_1GLG', 'GGR_2FW0', 'GGR_2FVY') for m in ('Ca', 'La')]
    first = families['RTX'][0]
    tasks += [(first, 'Ca', 'rotate'), (first, 'La', 'rotate'), (first, 'La', 'permute')]
    if len(tasks) != 55:
        raise InvalidArtifact('finite call count changed')
    return tasks


def task(name, metal, variant, row, out=None):
    g = row['connected' if variant == 'connected' else 'primary']
    src = read_json(verify(row['case']['preparation'])); coordinates = copy.deepcopy(g['all_Ca_atoms'])
    coordinates[g['selected_metal_index']][0] = metal
    ids = list(g['atom_group_indices']); order = list(range(len(ids))); rot = rotation() if variant == 'rotate' else np.eye(3)
    if variant == 'rotate':
        coordinates = [[a[0], *(np.asarray(a[1:]) @ rot.T)] for a in coordinates]
    if variant == 'permute':
        order.reverse(); coordinates = [coordinates[i] for i in order]; ids = [ids[i] for i in order]
    selected = order.index(g['selected_metal_index'])
    background = [order.index(i) for i in g['all_metal_indices'] if i != g['selected_metal_index']]
    charge = src['endpoints'][metal]['charge']; tid = name + '_' + metal + '_' + variant
    t = {'task_id': tid, 'case_id': name, 'kind': 'full', 'metal': metal, 'variant': variant,
         'charge': charge, 'spin_multiplicity': 1, 'state': check_atoms(coordinates, charge, background_calcium_indices=background),
         'background_calcium_indices': background, 'selected_metal_index': selected,
         'preparation': row['case']['preparation'], 'atom_group_indices': ids,
         'group_charges_e': g['endpoint_group_charges_e'][metal], 'charge_group_adapter': ADAPTER,
         'source_atom_ids': [g['atom_ids'][i] for i in order], 'source_atom_permutation': order,
         'rotation_matrix': rot.tolist()}
    if out:
        path = out / (tid + '.xyz'); write_xyz(path, coordinates); t['xyz'] = record(path)
    return t, coordinates


def prepare(preparation, output):
    from mace_global_benchmark import snapshot
    p, cfg, cases = prepared(preparation); parent = read_json(verify(cfg['parent_manifest']))
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    implementation = snapshot(out, (*parent['implementation'], 'mace_group_transfer.py', 'mace_site_groups.py', 'coordination_policy.py'))
    model = copy.deepcopy(parent['model']); model['chemical_group_policy'] = POLICY
    m = {'schema_version': SCHEMA, 'protocol_id': PROTOCOL, 'stage': 'MACE', 'preparation': record(preparation),
         'parent_manifest': parent['parent_manifest'], 'development_parent': cfg['parent_manifest'],
         'model': model, 'software': parent['software'], 'implementation': implementation,
         'agreement': cfg['agreement'], 'tolerances': TOL, 'compute_budget': None,
         'tasks': [task(c, metal, v, cases[c], out)[0] for c, metal, v in specs(cfg)],
         'reference': None, 'calibrated_class': None, 'baseline_changed': False, 'new_DFT_calls': 0}
    for t in m['tasks']:
        t['cache_key'] = cache_key({'task': t, 'model': model, 'software': m['software'], 'implementation': implementation})
    write_new(out / 'manifest.json', m)
    return validate(out / 'manifest.json')


@cached_file_checks
def validate(manifest):
    m = read_json(manifest); p, cfg, cases = prepared(verify(m['preparation']))
    parent = read_json(verify(cfg['parent_manifest'])); model = copy.deepcopy(parent['model']); model['chemical_group_policy'] = POLICY
    if (m['schema_version'] != SCHEMA or m['protocol_id'] != PROTOCOL or m['model'] != model or
        m['parent_manifest'] != parent['parent_manifest'] or m['development_parent'] != cfg['parent_manifest'] or
        m['software'] != parent['software'] or m['tolerances'] != TOL or m['agreement'] != cfg['agreement'] or len(m['tasks']) != 55):
        raise InvalidArtifact('declared scientific model or inventory changed')
    if kernel_parent_gate(m)['status'] != 'pass':
        raise InvalidArtifact('existing analytic kernel qualification failed')
    sm = read_json(verify(m['software']))
    for pin in [m['agreement'], *m['implementation'].values(), sm['python'], sm['checkpoint'], sm['requirements'], sm['backend_source_inventory']]:
        verify(pin)
    for t, (name, metal, variant) in zip(m['tasks'], specs(cfg)):
        want, atoms = task(name, metal, variant, cases[name])
        if {k: v for k, v in t.items() if k not in ('xyz', 'cache_key')} != want:
            raise InvalidArtifact('source-mapped task settings differ')
        actual = xyz(verify(t['xyz']))
        if ([a[0] for a in actual] != [a[0] for a in atoms] or
            not np.allclose([a[1:] for a in actual], [a[1:] for a in atoms], atol=1e-12, rtol=0)):
            raise InvalidArtifact('source geometry or declared rigid transform differs')
        payload = {k: v for k, v in t.items() if k != 'cache_key'}
        if t['cache_key'] != cache_key({'task': payload, 'model': m['model'], 'software': m['software'], 'implementation': m['implementation']}):
            raise InvalidArtifact('scientific cache mismatch')
    return {'status': 'pass', 'tasks': 55, 'manifest': record(manifest)}


@cached_file_checks
def collect(manifest):
    validate(manifest); mp = Path(manifest).resolve(); m = read_json(mp); rows = {}; attempts = []
    for t in m['tasks']:
        good = []
        for a in sorted((mp.parent / 'execution' / t['task_id']).glob('attempt_*')):
            r = accepted_attempt(a, t, mp)
            attempts.append({'task_id': t['task_id'], 'path': str(a), 'accepted': r is not None,
                             'receipt': record(a / 'receipt.json') if (a / 'receipt.json').exists() else None})
            if r is not None:
                good.append(r)
        rows[t['task_id']] = good[-1] if good else {'status': 'unavailable', 'energy_eV': None}
    return {'status': 'complete' if all(r['status'] == 'computed' for r in rows.values()) else 'incomplete',
            'manifest': record(mp), 'rows': rows, 'attempts': attempts}


def report(manifest, output):
    start = time.monotonic(); c = collect(manifest); m = read_json(manifest)
    p, cfg, cases = prepared(verify(m['preparation'])); families = family_cases(cfg)
    rows = c['rows']; scores = {}; checks = []; group_checks = []; comparisons = []
    def check(name, error, tolerance):
        checks.append({'name': name, 'error': error, 'tolerance': tolerance,
                       'pass': error is not None and abs(error) <= tolerance})
    def contrast(ca, la):
        a, b = rows[ca], rows[la]
        return (a['energy_eV'] - b['energy_eV']) * EV_TO_KCAL if a['status'] == b['status'] == 'computed' else None
    for family, names in families.items():
        ca = names[0] + '_Ca_primary'
        for name in names:
            scores[name] = {'R_model_kcal': contrast(ca, name + '_La_primary'),
                            'Ca_task': ca, 'La_task': name + '_La_primary', 'source': cfg['cases'][name],
                            'calibrated_class': None}
    for name in ('GGR_1GLG', 'GGR_2FW0', 'GGR_2FVY'):
        alternate = contrast(name + '_Ca_connected', name + '_La_connected'); primary = scores[name]['R_model_kcal']
        delta = alternate - primary if alternate is not None and primary is not None else None
        group_checks.append({'case': name, 'alternate_R_model_kcal': alternate, 'connected_minus_primary': delta,
                             'tolerance': 2., 'pass': delta is not None and abs(delta) <= 2.})
    first = families['RTX'][0]
    for metal, variant in (('Ca', 'rotate'), ('La', 'rotate'), ('La', 'permute')):
        a, b = rows[first + '_' + metal + '_primary'], rows[first + '_' + metal + '_' + variant]
        de = df = None
        if a['status'] == b['status'] == 'computed':
            de = (b['energy_eV'] - a['energy_eV']) * EV_TO_KCAL
            fa, fb = np.load(verify(a['forces'])), np.load(verify(b['forces']))
            want = fa @ rotation().T if variant == 'rotate' else fa[::-1]
            df = float(np.max(np.abs(want - fb)))
        check(metal + '_' + variant + '_energy_model_kcal', de, .01)
        check(metal + '_' + variant + '_force_eV_A', df, .001)
    for tid, row in rows.items():
        error = None
        if row['status'] == 'computed':
            trace = read_json(verify(row['charge_groups']))
            error = max(abs(row['total_charge_error_e']), *(s['maximum_channel_charge_error_e'] for s in trace['stages']))
        check(tid + '_charge_closure_e', error, 1e-5)
    def add(higher, lower, hi, lo, stratum):
        delta = hi - lo if hi is not None and lo is not None else None
        comparisons.append({'higher_expected': higher, 'lower_expected': lower, 'stratum': stratum,
                            'margin_model_kcal': delta, 'pass': delta is not None and delta > .02})
    get = lambda name: scores[name]['R_model_kcal']
    add('PQQ_4MAE', 'PQQ_1H4I', get('PQQ_4MAE'), get('PQQ_1H4I'), 'PQQ_functional_class_development')
    for a in ('ALPHA_1F6S', 'ALPHA_6IP9'):
        for g in ('GGR_1GLG', 'GGR_2FW0', 'GGR_2FVY'):
            add(a, g, get(a), get(g), 'qualified_alpha_vs_direct_GGR_development')
    means = {}
    for family in ('A0A7', 'HEW5', 'RTX'):
        values = [get(n) for n in families[family]]
        means[family] = sum(values) / len(values) if all(v is not None for v in values) else None
        for g in ('GGR_1GLG', 'GGR_2FW0', 'GGR_2FVY'):
            add(family + '_all_site_mean', g, means[family], get(g), 'Khoury_protein_level_supporting_proxy')
    for name in families['PARV_4CPV']:
        for g in ('GGR_1GLG', 'GGR_2FW0', 'GGR_2FVY'):
            add(name, g, get(name), get(g), 'parvalbumin_cross_study_supporting')
    numeric = c['status'] == 'complete' and all(x['pass'] for x in checks)
    represented = all(x['pass'] for x in group_checks)
    result = {'status': c['status'], 'protocol_id': PROTOCOL, 'collection': c, 'scores': scores,
              'checks': checks, 'grouping_checks': group_checks, 'numerical_checks_pass': numeric,
              'representation_checks_pass': represented, 'comparisons': comparisons, 'denominator': 22,
              'raw_pass_count': sum(x['pass'] for x in comparisons),
              'qualified_pass_count': sum(x['pass'] for x in comparisons) if numeric and represented else 0,
              'Khoury_means': means, 'aequorin_ordered': [{'site': n, 'R_model_kcal': get(n), 'site_label': None} for n in families['AEQ_1SL8']],
              'biological_groups': 9, 'directionally_labeled_groups': 8, 'prospectively_blind': False,
              'baseline_changed': False, 'absolute_reference': None, 'calibrated_class': None,
              'affinity_free_energy': None, 'response_status': 'response_model_not_validated',
              'prior_grouped_solvent_result_unchanged': cfg['old_report'], 'wall_seconds': time.monotonic() - start}
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    source = out / 'report_source.py'; source.write_bytes(Path(__file__).read_bytes()); result['implementation'] = record(source)
    write_new(out / 'result.json', result)
    return {k: result[k] for k in ('status', 'numerical_checks_pass', 'representation_checks_pass', 'raw_pass_count', 'denominator', 'qualified_pass_count', 'comparisons')}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); sub = ap.add_subparsers(dest='command', required=True)
    for op, names in [('prepare-inputs', ('config', 'output')), ('prepare', ('preparation', 'output')),
                      ('validate', ('manifest',)), ('collect', ('manifest',)), ('report', ('manifest', 'output'))]:
        p = sub.add_parser(op)
        for name in names:
            p.add_argument('--' + name, required=True)
    args = vars(ap.parse_args()); op = args.pop('command').replace('-', '_')
    print(json.dumps(globals()[op](**args), indent=2))
