"""Frozen large-checkpoint transfer of the existing typed-group physical panel."""
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import accepted_attempt, check_atoms, EV_TO_KCAL
from mace_charge_groups import kernel_parent_gate, TOL
from mace_group_constraints import ADAPTER
from mace_site_groups import POLICY
from mace_group_transfer import prepared, report_collection, validate as validate_medium

SCHEMA = 'alquemia.mace_group_large.v1'
PROTOCOL = 'intact_POLAR_large_typed31_group_vacuum_compatibility_v1'
IDENTITY = {'energy_model_kcal': .01, 'force_eV_A': .001, 'density_coefficients': 1e-8}


@cached_file_checks
def sources(config):
    cfg = read_json(config)
    for pin in cfg.values():
        verify(pin)
    small_path = verify(cfg['medium_manifest']); validate_medium(small_path); small = read_json(small_path)
    large_path = verify(cfg['large_parent']); large = read_json(large_path)
    previous = read_json(verify(cfg['medium_report'])); native = read_json(verify(cfg['large_native_collection']))
    if (previous['collection']['manifest'] != cfg['medium_manifest'] or previous['status'] != 'complete' or
            native['manifest'] != cfg['large_parent'] or native['status'] != 'complete' or
            not native['numerical_checks_pass'] or large['checkpoint_label'] != 'large'):
        raise InvalidArtifact('complete real medium and native-large sources required')
    ignore = {'checkpoint', 'pair_kernel', 'charge_group_adapter', 'chemical_group_policy', 'native_polar_source'}
    if ({k: v for k, v in small['model'].items() if k not in ignore} !=
            {k: v for k, v in large['model'].items() if k not in ignore}):
        raise InvalidArtifact('large parent differs beyond checkpoint and qualified tiling')
    old_tasks = {t['task_id']: t for t in large['tasks']}
    for metal in ('Ca', 'La'):
        tid = 'GGR_1GLG_' + metal + '_primary'; t = old_tasks[tid]
        typed = next(t for t in small['tasks'] if t['task_id'] == tid)
        if t['preparation'] != typed['preparation'] or t['charge'] != typed['charge']:
            raise InvalidArtifact('native/grouped GGR physical source differs')
        valid = [accepted_attempt(a, t, large_path) for a in sorted((large_path.parent / 'execution' / tid).glob('attempt_*'))]
        valid = [r for r in valid if r is not None]
        if not valid or valid[-1] != native['rows'][tid]:
            raise InvalidArtifact('archived native-large endpoint receipt differs')
    return cfg, small, large, previous, native


def model(small, large):
    result = copy.deepcopy(large['model'])
    for key in ('charge_group_adapter', 'chemical_group_policy', 'native_polar_source'):
        result[key] = small['model'][key]
    if result['chemical_group_policy'] != POLICY or result['charge_group_adapter'] != ADAPTER:
        raise InvalidArtifact('grouped model policy changed')
    return result


def tasks(small, large):
    result = []
    for metal in ('Ca', 'La'):
        t = copy.deepcopy(next(t for t in large['tasks'] if t['task_id'] == 'GGR_1GLG_' + metal + '_primary'))
        atoms = xyz(verify(t['xyz'])); check_atoms(atoms, t['charge'])
        t.pop('cache_key'); t['task_id'] = 'GGR_1GLG_' + metal + '_one_group'; t['variant'] = 'one_group'
        t.update(atom_group_indices=[0] * len(atoms), group_charges_e=[t['charge']], charge_group_adapter=ADAPTER)
        result.append(t)
    for t in small['tasks']:
        result.append({k: copy.deepcopy(v) for k, v in t.items() if k != 'cache_key'})
    return result


def prepare(config, output):
    from mace_global_benchmark import snapshot
    cfg, small, large, previous, native = sources(config)
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    implementation = snapshot(out, (*small['implementation'], 'mace_group_large.py', 'mace_group_canonical.py',
                                     'mace_canonical_run.py', 'mace_curvature.py', 'affordable_response.py'))
    m = {'schema_version': SCHEMA, 'protocol_id': PROTOCOL, 'stage': 'MACE', 'config': record(config),
         'preparation': small['preparation'], 'parent_manifest': cfg['large_parent'],
         'development_parent': cfg['medium_manifest'], 'native_reference': cfg['large_native_collection'],
         'model': model(small, large), 'software': large['software'], 'implementation': implementation,
         'agreement': cfg['agreement'], 'tolerances': TOL, 'native_identity_tolerances': IDENTITY,
         'tasks': tasks(small, large), 'compute_budget': None, 'reference': None, 'calibrated_class': None,
         'baseline_changed': False, 'new_DFT_calls': 0}
    for t in m['tasks']:
        t['cache_key'] = cache_key({'task': t, 'model': m['model'], 'software': m['software'], 'implementation': implementation})
    write_new(out / 'manifest.json', m)
    return validate(out / 'manifest.json')


@cached_file_checks
def validate(manifest):
    m = read_json(manifest); cfg, small, large, previous, native = sources(verify(m['config']))
    if (m['schema_version'] != SCHEMA or m['protocol_id'] != PROTOCOL or m['model'] != model(small, large) or
            m['preparation'] != small['preparation'] or m['parent_manifest'] != cfg['large_parent'] or
            m['development_parent'] != cfg['medium_manifest'] or m['native_reference'] != cfg['large_native_collection'] or
            m['agreement'] != cfg['agreement'] or m['software'] != large['software'] or
            m['tolerances'] != TOL or m['native_identity_tolerances'] != IDENTITY or len(m['tasks']) != 57):
        raise InvalidArtifact('declared large model, state or finite inventory changed')
    if kernel_parent_gate(m)['status'] != 'pass':
        raise InvalidArtifact('archived large analytic kernel qualification failed')
    sm = read_json(verify(m['software']))
    for pin in [*m['implementation'].values(), sm['checkpoint'], sm['python'], sm['requirements'], sm['backend_source_inventory']]:
        verify(pin)
    for name in ('mace_site_groups.py', 'coordination_policy.py', 'mace_group_constraints.py',
                 'mace_charge_group_prepare.py', 'mace_analytic.py', 'mace_realspace_compat.py'):
        if m['implementation'][name]['sha256'] != small['implementation'][name]['sha256']:
            raise InvalidArtifact('shared scientific implementation changed: ' + name)
    for t, expected in zip(m['tasks'], tasks(small, large)):
        payload = {k: v for k, v in t.items() if k != 'cache_key'}
        if payload != expected:
            raise InvalidArtifact('paired input coordinates/group state changed')
        verify(t['xyz'])
        if t['cache_key'] != cache_key({'task': payload, 'model': m['model'], 'software': m['software'], 'implementation': m['implementation']}):
            raise InvalidArtifact('large scientific cache mismatch')
    return {'status': 'pass', 'tasks': 57, 'native_identity_controls': 2, 'manifest': record(manifest)}


def accepted_rows(manifest):
    mp = Path(manifest).resolve(); m = read_json(mp); rows = {}; attempts = []
    for t in m['tasks']:
        good = []
        for path in sorted((mp.parent / 'execution' / t['task_id']).glob('attempt_*')):
            r = accepted_attempt(path, t, mp)
            attempts.append({'task_id': t['task_id'], 'path': str(path), 'accepted': r is not None,
                             'receipt': record(path / 'receipt.json') if (path / 'receipt.json').exists() else None})
            if r is not None:
                good.append(r)
        rows[t['task_id']] = good[-1] if good else {'status': 'unavailable', 'energy_eV': None}
    return rows, attempts


@cached_file_checks
def native_gate(manifest):
    m = read_json(manifest); native = read_json(verify(m['native_reference'])); rows, attempts = accepted_rows(manifest)
    checks = []
    for metal in ('Ca', 'La'):
        r = rows['GGR_1GLG_' + metal + '_one_group']; old = native['rows']['GGR_1GLG_' + metal + '_primary']
        values = dict.fromkeys(IDENTITY)
        if r['status'] == 'computed':
            values = {'energy_model_kcal': (r['energy_eV'] - old['energy_eV']) * EV_TO_KCAL,
                      'force_eV_A': float(np.max(np.abs(np.load(verify(r['forces'])) - np.load(verify(old['forces']))))),
                      'density_coefficients': float(np.max(np.abs(np.load(verify(r['density_coefficients'])) - np.load(verify(old['density_coefficients'])))))}
        for name, error in values.items():
            checks.append({'name': 'native_' + metal + '_' + name, 'error': error, 'tolerance': IDENTITY[name],
                           'pass': error is not None and abs(error) <= IDENTITY[name]})
    return {'status': 'pass' if all(c['pass'] for c in checks) else 'unavailable_or_failed', 'checks': checks}


@cached_file_checks
def collect(manifest):
    validate(manifest); rows, attempts = accepted_rows(manifest)
    return {'status': 'complete' if all(r['status'] == 'computed' for r in rows.values()) else 'incomplete',
            'manifest': record(manifest), 'rows': rows, 'attempts': attempts, 'native_identity': native_gate(manifest)}


def report(manifest, output):
    c = collect(manifest); m = read_json(manifest); cfg, small, large, previous, native = sources(verify(m['config']))
    prep, source_cfg, cases = prepared(verify(m['preparation']))
    result = report_collection(c, m, source_cfg, cases, output, extra_checks=c['native_identity']['checks'])
    paired = []
    for old, new in zip(previous['comparisons'], result['comparisons']):
        keys = ('higher_expected', 'lower_expected', 'stratum')
        if any(old[k] != new[k] for k in keys):
            raise InvalidArtifact('comparison identity changed')
        delta = new['margin_model_kcal'] - old['margin_model_kcal'] if new['margin_model_kcal'] is not None else None
        paired.append({**new, 'medium_margin_model_kcal': old['margin_model_kcal'], 'medium_pass': old['pass'],
                       'large_minus_medium_margin_model_kcal': delta})
    comparison = {'large_result': record(Path(output).resolve() / 'result.json'), 'medium_result': cfg['medium_report'],
                  'paired': paired, 'medium_raw_pass_count': previous['raw_pass_count'],
                  'large_raw_pass_count': result['raw_pass_count'],
                  'declared_directional_improvement': c['status'] == 'complete' and result['raw_pass_count'] > 17 and all(p['pass'] for p in paired[:7]),
                  'numerical_checks_pass': result['numerical_checks_pass'],
                  'representation_checks_pass': result['representation_checks_pass'], 'baseline_changed': False}
    src = Path(output) / 'checkpoint_report_source.py'; src.write_bytes(Path(__file__).read_bytes())
    comparison['implementation'] = record(src)
    write_new(Path(output) / 'checkpoint_comparison.json', comparison)
    return {k: v for k, v in comparison.items() if k != 'paired'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); sub = ap.add_subparsers(dest='command', required=True)
    for op, names in [('prepare', ('config', 'output')), ('validate', ('manifest',)),
                      ('collect', ('manifest',)), ('native-gate', ('manifest',)), ('report', ('manifest', 'output'))]:
        p = sub.add_parser(op)
        for name in names:
            p.add_argument('--' + name, required=True)
    args = vars(ap.parse_args()); op = args.pop('command').replace('-', '_')
    print(json.dumps(globals()[op](**args), indent=2))
