"""Calibrate the unchanged typed-group POLAR descriptor on the frozen PQQ panel."""
import argparse
import copy
import json
import math
from pathlib import Path
import time
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_canonical_run import inventory
from mace_charge_group_prepare import protein_ledger
from mace_charge_groups import kernel_parent_gate, TOL
from mace_group_transfer import PROTOCOL, paired_source, task
from mace_site_groups import POLICY, build
from mace_hybrid import accepted_attempt, EV_TO_KCAL

SCHEMA = 'alquemia.mace_group_canonical.v1'
DECISION = {'minimum_gap_model_kcal': .02, 'Ca': 'R<=max_Ca_calibration',
            'La': 'R>=min_La_calibration', 'inside_gap': 'inconclusive',
            'calibration_count': 25, 'transfer_count': 3}


@cached_file_checks
def config_data(path):
    cfg = read_json(path)
    for key in ('agreement', 'inventory', 'parent_manifest', 'parent_report', 'typed_donor_source'):
        verify(cfg[key])
    data = inventory(verify(cfg['inventory']))
    rows = {r['case_id']: r for r in data['rows'] if r['evaluation_role'] == 'calibration'}
    if (set(cfg['cases']) != set(rows) or cfg['cutoff_A'] != 3.1 or
            cfg['minimum_gap_model_kcal'] != DECISION['minimum_gap_model_kcal']):
        raise InvalidArtifact('declared PQQ inventory or calibration policy changed')
    sources = {}
    for name, c in cfg['cases'].items():
        r = rows[name]; p = read_json(verify(c['preparation'])); paired_source(p)
        if (c['expected_class'] != r['expected_class'] or c['evaluation_role'] != 'calibration' or
                c['prospectively_blind'] or p['case_id'] != 'PQQ_' + name or
                p['status'] != 'prepared' or p['evidence']['direction'] != r['expected_class'] or
                p['explicit_waters'] != [] or p['cofactor_charge_e'] != -3 or
                p['source_preparation'] != r['source_manifest'] or
                p['source_audit_row']['endpoints'] != r['endpoints']):
            raise InvalidArtifact('source identity, evidence or physical state changed: ' + name)
        for key in ('source', 'source_preparation', 'forcefield', 'water_forcefield'):
            verify(p[key])
        sources[name] = p
    if sum(c['expected_class'] == 'La' for c in cfg['cases'].values()) != 11:
        raise InvalidArtifact('eleven La and fourteen Ca reference labels required')
    return cfg, data, sources


def prepare_inputs(config, output):
    start, cpu = time.monotonic(), time.process_time()
    cfg, data, sources = config_data(config)
    if record(Path(__file__).with_name('coordination_policy.py'))['sha256'] != cfg['typed_donor_source']['sha256']:
        raise InvalidArtifact('typed donor implementation changed')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    source = out / 'implementation'; source.mkdir()
    pins = {}
    for name in ('mace_group_canonical.py', 'mace_group_transfer.py', 'mace_site_groups.py',
                 'mace_charge_group_prepare.py', 'coordination_policy.py'):
        dest = source / name; dest.write_bytes(Path(__file__).with_name(name).read_bytes()); pins[name] = record(dest)
    cases = {}
    for name, p in sources.items():
        ledger, bonds = protein_ledger(p)
        row = {'case': cfg['cases'][name], 'protein_ledger': ledger, 'primary': build(p, ledger, bonds)}
        dest = out / (name + '.json'); write_new(dest, row); cases[name] = record(dest)
    result = {'policy_id': POLICY, 'status': 'prepared', 'config': record(config), 'cases': cases,
              'implementation': pins, 'model_calls': 0, 'wall_seconds': time.monotonic() - start,
              'CPU_seconds': time.process_time() - cpu}
    write_new(out / 'preparation.json', result)
    return {'status': 'prepared', 'cases': len(cases), 'preparation': record(out / 'preparation.json'),
            'wall_seconds': result['wall_seconds'], 'CPU_seconds': result['CPU_seconds']}


@cached_file_checks
def prepared(path):
    p = read_json(path); cfg, data, sources = config_data(verify(p['config']))
    if p['status'] != 'prepared' or p['policy_id'] != POLICY or set(p['cases']) != set(cfg['cases']):
        raise InvalidArtifact('prepared inventory differs')
    for pin in p['implementation'].values():
        verify(pin)
    cases = {}
    for name, pin in p['cases'].items():
        row = read_json(verify(pin)); src = sources[name]
        if row['case'] != cfg['cases'][name]:
            raise InvalidArtifact('source case changed')
        bonds = sorted(tuple(sorted((b['atom_a_id'], b['atom_b_id']))) for b in src['preparation_details']['bonds'])
        expected = json.loads(json.dumps(build(src, row['protein_ledger'], bonds)))
        stored = copy.deepcopy(row['primary'])
        if len(stored['contacts']) != len(expected['contacts']):
            raise InvalidArtifact('contact inventory changed')
        for a, b in zip(stored['contacts'], expected['contacts']):
            if abs(a['distance_A'] - b['distance_A']) > 1e-12:
                raise InvalidArtifact('source contact distance changed')
            a['distance_A'] = b['distance_A']
        if stored != expected:
            raise InvalidArtifact('source-derived groups changed')
        cases[name] = row
    return p, cfg, data, cases


@cached_file_checks
def parent_results(cfg):
    from mace_group_transfer import validate as validate_parent
    mp = verify(cfg['parent_manifest']); validate_parent(mp); m = read_json(mp)
    r = read_json(verify(cfg['parent_report']))
    if (r['collection']['manifest'] != cfg['parent_manifest'] or r['protocol_id'] != PROTOCOL or
            r['status'] != 'complete' or not r['numerical_checks_pass']):
        raise InvalidArtifact('completed identical typed-group source required')
    rows = {}; tasks = {t['task_id']: t for t in m['tasks']}
    for case in ('1H4I', '4MAE'):
        for metal in ('Ca', 'La'):
            tid = 'PQQ_' + case + '_' + metal + '_primary'; found = []
            for a in sorted((mp.parent / 'execution' / tid).glob('attempt_*')):
                result = accepted_attempt(a, tasks[tid], mp)
                if result is not None:
                    found.append(result)
            if not found or found[-1] != r['collection']['rows'][tid]:
                raise InvalidArtifact('cached crystal receipt does not match report')
            rows[case + '_' + metal + '_primary'] = found[-1]
    return m, r, rows


def specs(cfg):
    return [(name, metal) for name in cfg['cases'] for metal in ('Ca', 'La')]


def source_support(path, cfg, parent):
    data = read_json(path); rows = {r['case_id']: r for r in data['rows']}
    if set(rows) != set(cfg['cases']) | {'1H4I', '4MAE', '1KB0'}:
        raise InvalidArtifact('full-protein source inventory changed')
    for name, c in cfg['cases'].items():
        if rows[name]['status'] != 'prepared' or rows[name]['preparation'] != c['preparation']:
            raise InvalidArtifact('source preparation differs from the existing valid inventory')
    for name in ('1H4I', '4MAE'):
        t = next(t for t in parent['tasks'] if t['task_id'] == 'PQQ_' + name + '_Ca_primary')
        if rows[name]['preparation'] != t['preparation']:
            raise InvalidArtifact('crystal source mismatch')
    missing = rows['1KB0']
    if missing['status'] != 'unsupported' or 'unsupported peptide connection' not in missing['reason']:
        raise InvalidArtifact('missing-source status changed')
    return {'1KB0': missing['reason']}


def prepare(preparation, output, source_preparations, clarification):
    from mace_global_benchmark import snapshot
    p, cfg, data, cases = prepared(preparation); parent, previous, reused = parent_results(cfg)
    unsupported = source_support(source_preparations, cfg, parent)
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    impl = snapshot(out, (*parent['implementation'], 'mace_group_canonical.py', 'mace_canonical_run.py',
                          'mace_curvature.py', 'affordable_response.py'))
    m = {'schema_version': SCHEMA, 'protocol_id': PROTOCOL, 'stage': 'MACE', 'preparation': record(preparation),
         'parent_manifest': parent['parent_manifest'], 'development_parent': cfg['parent_manifest'],
         'model': parent['model'], 'software': parent['software'], 'implementation': impl,
         'agreement': cfg['agreement'], 'tolerances': TOL, 'compute_budget': None,
         'decision_policy': DECISION, 'reference': None, 'calibrated_class': None,
         'baseline_changed': False, 'new_DFT_calls': 0,
         'tasks': [task(name, metal, 'primary', cases[name], out)[0] for name, metal in specs(cfg)],
         'reused': {k: {'source_report': cfg['parent_report'], 'source_task_id': v['task_id']} for k, v in reused.items()},
         'source_preparations': record(source_preparations), 'scope_clarification': record(clarification),
         'unsupported': unsupported}
    for t in m['tasks']:
        t['cache_key'] = cache_key({'task': t, 'model': m['model'], 'software': m['software'], 'implementation': impl})
    write_new(out / 'manifest.json', m)
    return validate(out / 'manifest.json')


@cached_file_checks
def validate(manifest):
    m = read_json(manifest); p, cfg, data, cases = prepared(verify(m['preparation']))
    parent, previous, reused = parent_results(cfg)
    unsupported = source_support(verify(m['source_preparations']), cfg, parent); verify(m['scope_clarification'])
    if (m['schema_version'] != SCHEMA or m['protocol_id'] != PROTOCOL or m['model'] != parent['model'] or
            m['parent_manifest'] != parent['parent_manifest'] or m['development_parent'] != cfg['parent_manifest'] or
            m['software'] != parent['software'] or m['tolerances'] != TOL or m['agreement'] != cfg['agreement'] or
            m['decision_policy'] != DECISION or len(m['tasks']) != 50 or
            m['unsupported'] != unsupported or
            m['reused'] != {k: {'source_report': cfg['parent_report'], 'source_task_id': v['task_id']} for k, v in reused.items()}):
        raise InvalidArtifact('model, inventory, reuse or decision policy differs')
    if kernel_parent_gate(m)['status'] != 'pass':
        raise InvalidArtifact('analytic kernel qualification failed')
    sm = read_json(verify(m['software']))
    for pin in [m['agreement'], *m['implementation'].values(), sm['python'], sm['checkpoint'], sm['requirements'], sm['backend_source_inventory']]:
        verify(pin)
    for name in ('mace_site_groups.py', 'coordination_policy.py', 'mace_group_constraints.py',
                 'mace_group_transfer.py', 'mace_realspace.py', 'mace_analytic.py'):
        if name in parent['implementation'] and m['implementation'][name]['sha256'] != parent['implementation'][name]['sha256']:
            raise InvalidArtifact('shared scientific implementation changed: ' + name)
    for t, (name, metal) in zip(m['tasks'], specs(cfg)):
        want, atoms = task(name, metal, 'primary', cases[name])
        if {k: v for k, v in t.items() if k not in ('xyz', 'cache_key')} != want:
            raise InvalidArtifact('task state or groups differ')
        actual = xyz(verify(t['xyz']))
        if ([a[0] for a in actual] != [a[0] for a in atoms] or
                not np.allclose([a[1:] for a in actual], [a[1:] for a in atoms], atol=1e-12, rtol=0)):
            raise InvalidArtifact('source coordinates changed')
        payload = {k: v for k, v in t.items() if k != 'cache_key'}
        if t['cache_key'] != cache_key({'task': payload, 'model': m['model'], 'software': m['software'], 'implementation': m['implementation']}):
            raise InvalidArtifact('scientific cache mismatch')
    return {'status': 'pass', 'tasks': 50, 'reused_endpoints': 4, 'unsupported_cases': 1, 'manifest': record(manifest)}


@cached_file_checks
def collect(manifest):
    validate(manifest); mp = Path(manifest).resolve(); m = read_json(mp)
    p, cfg, data, cases = prepared(verify(m['preparation'])); parent, previous, reused = parent_results(cfg)
    rows = {}; attempts = []
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
            'manifest': record(mp), 'rows': rows, 'reused_rows': reused, 'attempts': attempts}


def calibrate(scores, numerical_pass):
    """No thresholds unless every declared reference is available and classes separate."""
    calibration = [s for s in scores.values() if s['evaluation_role'] == 'calibration']
    complete = len(calibration) == 25 and all(s['R_model_kcal'] is not None and math.isfinite(s['R_model_kcal']) for s in calibration)
    result = {'status': 'unavailable_incomplete', 'Ca_max_R_model_kcal': None, 'La_min_R_model_kcal': None,
              'observed_class_gap_model_kcal': None, 'production_qualified': False}
    if complete:
        c = max(s['R_model_kcal'] for s in calibration if s['expected_class'] == 'Ca')
        l = min(s['R_model_kcal'] for s in calibration if s['expected_class'] == 'La')
        result['observed_class_gap_model_kcal'] = l - c
        result['status'] = 'unavailable_numerical_failure' if not numerical_pass else 'unavailable_class_overlap'
        if numerical_pass and l - c > DECISION['minimum_gap_model_kcal']:
            result.update(status='available_for_research', Ca_max_R_model_kcal=c, La_min_R_model_kcal=l)
    return result


def classify(value, calibration):
    if value is None or calibration['status'] != 'available_for_research':
        return None
    if value <= calibration['Ca_max_R_model_kcal']:
        return 'Ca'
    if value >= calibration['La_min_R_model_kcal']:
        return 'La'
    return 'inconclusive'


def report(manifest, output):
    start = time.monotonic(); collection = collect(manifest); m = read_json(manifest)
    p, cfg, data, cases = prepared(verify(m['preparation'])); parent, previous, reused = parent_results(cfg)
    rows = {**collection['rows'], **collection['reused_rows']}; scores = {}; checks = []
    for tid, r in rows.items():
        error = None
        if r['status'] == 'computed':
            trace = read_json(verify(r['charge_groups']))
            error = max(abs(r['total_charge_error_e']), *(s['maximum_channel_charge_error_e'] for s in trace['stages']))
        checks.append({'task_id': tid, 'charge_error_e': error, 'tolerance_e': 1e-5,
                       'pass': error is not None and error <= 1e-5})
    for r in data['rows']:
        name = r['case_id']; ca = rows.get(name + '_Ca_primary'); la = rows.get(name + '_La_primary')
        ready = ca is not None and la is not None and ca['status'] == la['status'] == 'computed'
        scores[name] = {'R_model_kcal': (ca['energy_eV'] - la['energy_eV']) * EV_TO_KCAL if ready else None,
                        'Ca_energy_eV': ca['energy_eV'] if ca else None, 'La_energy_eV': la['energy_eV'] if la else None,
                        'status': 'computed' if ready else m['unsupported'].get(name, 'unavailable'),
                        'expected_class': r['expected_class'], 'evaluation_role': r['evaluation_role'],
                        'evidence_stratum': r['evidence_stratum'], 'sequence_accession_group': r['sequence_accession_group'],
                        'homology_independence_claimed': r['homology_independence_claimed'],
                        'prospectively_blind': False, 'replicate_group': data['sequence_group_evidence'].get(name),
                        'baseline': r['baseline'], 'research_class': None, 'production_class': None}
    numeric = collection['status'] == 'complete' and all(c['pass'] for c in checks)
    calibration = calibrate(scores, numeric)
    for s in scores.values():
        s['research_class'] = classify(s['R_model_kcal'], calibration)
    comparisons = []
    la = [(n, s) for n, s in scores.items() if s['evaluation_role'] == 'calibration' and s['expected_class'] == 'La']
    ca = [(n, s) for n, s in scores.items() if s['evaluation_role'] == 'calibration' and s['expected_class'] == 'Ca']
    for ln, ls in la:
        for cn, cs in ca:
            margin = ls['R_model_kcal'] - cs['R_model_kcal'] if ls['R_model_kcal'] is not None and cs['R_model_kcal'] is not None else None
            comparisons.append({'La': ln, 'Ca': cn, 'margin_model_kcal': margin,
                                'pass': margin is not None and margin > .02})
    transfer = [s for s in scores.values() if s['evaluation_role'] != 'calibration']
    result = {'status': collection['status'], 'protocol_id': PROTOCOL, 'manifest': record(manifest),
              'collection': collection, 'scores': scores, 'calibration': calibration, 'decision_policy': DECISION,
              'checks': checks, 'numerical_checks_pass': numeric, 'prior_numerical_checks': cfg['parent_report'],
              'representation_checks_pass': previous['representation_checks_pass'], 'grouping_checks': previous['grouping_checks'],
              'calibration_pairwise': comparisons, 'pairwise_pass_count': sum(c['pass'] for c in comparisons),
              'pairwise_scored_count': sum(c['margin_model_kcal'] is not None for c in comparisons),
              'pairwise_denominator': 154, 'retrospective_transfer_denominator': 3,
              'retrospective_transfer_scored': sum(s['R_model_kcal'] is not None for s in transfer),
              'retrospective_transfer_classified': sum(s['research_class'] in ('Ca', 'La') for s in transfer),
              'retrospective_transfer_inconclusive': sum(s['research_class'] == 'inconclusive' for s in transfer),
              'retrospective_transfer_correct': sum(s['research_class'] == s['expected_class'] for s in transfer),
              'independent_transfer_validation': False, 'baseline_changed': False, 'absolute_aquo_reference': None,
              'score_formula': '(E_Ca_eV - E_La_eV) * eV_to_model_kcal', 'eV_to_model_kcal': EV_TO_KCAL,
              'endpoint_units': 'eV', 'score_units': 'model_kcal', 'S_kcal_mol': None,
              'source_preparations': m['source_preparations'], 'scope_clarification': m['scope_clarification'],
              'new_DFT_calls': 0, 'calibration_count': 25,
              'calibration_scored_count': sum(s['R_model_kcal'] is not None for s in scores.values() if s['evaluation_role'] == 'calibration'),
              'production_qualified': False,
              'wall_seconds': time.monotonic() - start}
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    src = out / 'report_source.py'; src.write_bytes(Path(__file__).read_bytes()); result['implementation'] = record(src)
    write_new(out / 'result.json', result)
    return {k: result[k] for k in ('status', 'calibration', 'numerical_checks_pass', 'pairwise_pass_count',
                                  'pairwise_denominator', 'retrospective_transfer_scored', 'retrospective_transfer_correct')}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); sub = ap.add_subparsers(dest='command', required=True)
    for op, names in [('prepare-inputs', ('config', 'output')), ('prepare', ('preparation', 'output', 'source_preparations', 'clarification')),
                      ('validate', ('manifest',)), ('collect', ('manifest',)), ('report', ('manifest', 'output'))]:
        p = sub.add_parser(op)
        for name in names:
            p.add_argument('--' + name.replace('_', '-'), required=True)
    args = vars(ap.parse_args()); op = args.pop('command').replace('-', '_')
    print(json.dumps(globals()[op](**args), indent=2))
