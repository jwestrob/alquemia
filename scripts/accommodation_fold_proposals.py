"""Fixed primary225 adapter for the existing physical proposal/selection engine."""
from __future__ import annotations
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import shutil
import statistics
import numpy as np

import accommodation_proposals as engine
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_site_kinematics import Kinematics
from mace_hybrid import check_atoms
from accommodation_folds_compare import source_groups, decision

SCOPE = 'primary225_fold_transfer_v1'


def validate_scope(m):
    design = read_json(verify(m['fold_design']))
    groups = source_groups(read_json(verify(m['sources'])))
    primary = {c['case_id'] for members in groups.values() for c in members if c['primary_evaluation_pool']}
    supported = {c['case_id'] for c in design['cases']}
    if (len(primary), len(supported), design['denominator']) != (225, 208, 225):
        raise InvalidArtifact('declared primary225/208 mapping population differs')
    if {r['case_id'] for r in design['rows']} != primary or not supported <= primary:
        raise InvalidArtifact('mapping source IDs differ from frozen primary pool')
    if {c['case_id'] for c in m['cases']} != supported or {c['case_id'] for c in m['unavailable_cases']} != primary - supported:
        raise InvalidArtifact('supported/unavailable population differs')
    if m['maximum_new_GFN2_singlepoints'] != 832 or m['algorithm_starts'] != 416:
        raise InvalidArtifact('finite fold scope differs')
    if m['settings'] != engine.SETTINGS or m['evaluation_scope'] != SCOPE:
        raise InvalidArtifact('fold proposal policy differs')
    ref = read_json(verify(m['proposal_reference']))
    if (ref['status'] != 'available' or ref['protocol_id'] != engine.PROTOCOL or
            ref['settings'] != engine.SETTINGS or ref['model'] != m['model'] or
            ref['orca'] != m['orca'] or ref['noncanonical_folds_used_for_calibration']):
        raise InvalidArtifact('proposal calibration incompatible or uses transfer cases')
    if sum(t['q0_status'] == 'available' for t in m['tasks']) != 415:
        raise InvalidArtifact('expected exactly one archived q0 failure; inspect discrepancy before execution')
    case_index = {c['case_id']: c for c in design['cases']}
    for t in m['tasks']:
        c = case_index[t['case_id']]; ep = c['origins'][t['metal']]
        if (t['xyz'], t['charge'], t['multiplicity'], t['mapping']) != (ep['xyz'], ep['charge'], ep['multiplicity'], c['maps'][t['metal']]):
            raise InvalidArtifact('fold endpoint changed from exact mapped source')
        if t['active_mode_ids'] != [c['modes'][role] for role in engine.ROLE_ORDER if role in c['modes']]:
            raise InvalidArtifact('fold active modes differ')


def prepare(design, comparison, reference, agreement, output, cpu_python, gpu_python):
    d, old, ref = map(read_json, (design, comparison, reference))
    inventory = read_json(verify(old['inventory']))
    if d['model'] != inventory['model'] or d['software'] != inventory['software']:
        raise InvalidArtifact('fold mapping and archived native model differ')
    source_groups(read_json(verify(old['sources'])))
    lookup = {}
    for pin in old['collections']:
        collection = read_json(verify(pin)); mp = verify(collection['manifest']); low = read_json(mp)
        for t in low['all_tasks']:
            key = (t['case_id'], t['representation'], t['metal'], t['medium'])
            if key in lookup: raise InvalidArtifact('duplicate archived solvent task')
            lookup[key] = (mp, t)
    inv_index = {c['case_id']: c for c in inventory['cases']}
    tasks, audits = [], []
    for c in d['cases']:
        parent = read_json(verify(c['source']['parent']))
        selector = parent['fixed_core']['requested_roles']['extra_acidic_ligand_homolog']
        rn = selector['resname'] if isinstance(selector, dict) else selector.split(':')[1][:3]
        roles = ['anchor_glutamate'] + (['extra_acidic_ligand_homolog'] if rn == 'ASP' else [])
        expected = [engine.role_mode(parent, role) for role in roles]
        if c['modes'] != dict(zip(roles, expected)):
            raise InvalidArtifact('fold role chemistry differs')
        for metal in ('Ca', 'La'):
            ep = c['origins'][metal]; atoms = xyz(verify(ep['xyz']))
            kin = Kinematics(read_json(verify(c['maps'][metal]))['context'])
            names = [mode['id'] for mode in kin.modes]; zero = np.zeros(len(names))
            difference = float(np.max(np.abs(kin.evaluate(zero)[1] - np.array([a[1:] for a in atoms]))))
            if difference > 1e-12: raise InvalidArtifact('fold origin mapping differs')
            checks = engine.geometry_check(kin, zero, [a[0] for a in atoms])
            if not checks['pass']: raise InvalidArtifact('fold origin geometry unsupported')
            check_atoms(atoms, ep['charge'])
            t = {'task_id': c['case_id'] + '__' + metal, 'case_id': c['case_id'], 'metal': metal,
                 'xyz': ep['xyz'], 'design_origin_xyz': ep['xyz'], 'executed_q0_difference_A': 0.0,
                 'map_origin_max_difference_A': difference, 'charge': ep['charge'], 'multiplicity': ep['multiplicity'],
                 'mapping': c['maps'][metal], 'source_preparation': c['preparation'],
                 'active_roles': roles, 'active_mode_ids': expected, 'active_indices': [names.index(x) for x in expected],
                 'mode_count': len(names), 'origin_geometry_checks': checks,
                 'label_scope': c['label_scope'], 'q0_status': 'unavailable', 'q0': None}
            try:
                nep = inv_index[c['case_id']]['representations']['context']['endpoints'][metal]
                nr = read_json(verify(nep['native_MACE_receipt']))
                nm = read_json(verify(nep['source_manifest']))
                nt = next(x for x in nm['tasks'] if x['task_id'] == nep['source_task_id'])
                if (nm['model'] != d['model'] or nr['status'] != 'computed' or nr['manifest'] != nep['source_manifest'] or
                        nr['energy_eV'] != nep['native_MACE_energy_eV'] or xyz(verify(nt['xyz'])) != atoms or
                        (nt['charge'], nt['spin_multiplicity']) != (ep['charge'], ep['multiplicity'])):
                    raise InvalidArtifact('fold archived native energy/state/coordinates differ')
                lows = {}
                for medium in ('vacuum', 'alpb'):
                    mp, lt = lookup[c['case_id'], 'context', metal, medium]
                    lows[medium] = engine.check_low_source(mp, lt, ep, d['orca'])
                if read_json(verify(lows['vacuum']['audit']['parameter_export'])) != read_json(verify(lows['alpb']['audit']['parameter_export'])):
                    raise InvalidArtifact('q0 GFN2 media parameter sets differ')
                t.update(q0_status='available', q0={'native_MACE_receipt': nep['native_MACE_receipt'], 'low': lows,
                    'components': {'MACE_eV': nr['energy_eV'], 'GFN2_vacuum_hartree': lows['vacuum']['energy_hartree'],
                                   'GFN2_ALPB_hartree': lows['alpb']['energy_hartree']}})
            except Exception as exc:
                t['q0_reason'] = str(exc)
            tasks.append(t); audits.append({k: t[k] for k in ('task_id', 'q0_status')} | {'reason': t.get('q0_reason')})
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    impl = out / 'implementation'; impl.mkdir(); pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        target = impl / p.name; shutil.copyfile(p, target); pins[p.name] = record(target)
    old_reference = verify(old['frozen_reference']); released = read_json(old_reference)
    unavailable = [{k: r[k] for k in ('case_id', 'root_case_id', 'biological_group', 'source_conditioning_metal',
                    'canonical_coordinate_match', 'expected_class_for_later_report_only', 'preparation_reason', 'status')}
                   for r in d['rows'] if r['status'] != 'supported']
    m = {'protocol_id': engine.PROTOCOL, 'evaluation_scope': SCOPE, 'settings': engine.SETTINGS,
         'tasks': tasks, 'cases': d['cases'], 'unavailable_cases': unavailable,
         'fold_design': record(design), 'fold_comparison': record(comparison), 'sources': old['sources'],
         'proposal_reference': record(reference), 'agreement': record(agreement),
         'frozen_comparison': record(old_reference), 'frozen_context_bands': released['calibration']['context']['bands'],
         'model': d['model'], 'software': d['software'], 'orca': d['orca'], 'implementation': pins,
         'cpu_python': str(Path(cpu_python).absolute()), 'gpu_python': str(Path(gpu_python).absolute()),
         'cpu_executable': record(cpu_python), 'gpu_executable': record(gpu_python),
         'algorithm_starts': 416, 'maximum_new_GFN2_singlepoints': 832, 'new_DFT_calls': 0,
         'threshold_refitted': False, 'baseline_changed': False,
         'resources': {'GPU_cpus': 32, 'GPUs': 1, 'GPU_host_mem_MiB': 200000,
                       'GFN2_cpus': 64, 'mpi_ranks': 8, 'concurrent_tasks': 8, 'GFN2_shards': 4}}
    write_new(out / 'manifest.json', m)
    write_new(out / 'Q0_AUDIT.json', {'manifest': record(out / 'manifest.json'), 'rows': audits,
        'declared_case_denominator': 225, 'prepared_endpoint_denominator': 416,
        'q0_reuse_available': sum(t['q0_status'] == 'available' for t in tasks), 'new_molecular_calls': 0})
    return engine.validate(out / 'manifest.json')


def outcome(call, expected):
    return ('unavailable' if call in (None, 'unavailable') else 'correct' if call == expected + '-supported'
            else 'inconclusive' if call == 'inconclusive' else 'wrong')


def summarize_rows(rows, methods):
    return {method: {'denominator': len(rows), **{key: sum(row['methods'][method]['outcome'] == key for row in rows)
             for key in ('correct', 'wrong', 'inconclusive', 'unavailable')}} for method in methods}


def strict_summary(members, index, method, bands, expected):
    values = [index[mid]['methods'][method]['R'] for mid in members]
    complete = all(v is not None for v in values)
    value = statistics.median(values) if complete else None
    call = decision(value, bands)
    return {'R': value, 'decision': call, 'outcome': outcome(call, expected), 'required_members': len(members),
            'missing_members': [mid for mid, v in zip(members, values) if v is None],
            'within_protein_range': max(values) - min(values) if complete else None}


def compare(manifest, selection, dft, output):
    m, selected, native, dr = read_json(manifest), read_json(selection), read_json(verify(read_json(manifest)['fold_comparison'])), read_json(dft)
    validate_scope(m)
    if selected['manifest'] != record(manifest) or selected['case_denominator'] != 225:
        raise InvalidArtifact('proposal collection scope differs')
    if dr['sources'] != m['sources']: raise InvalidArtifact('DFT source panel differs')
    ref = read_json(verify(m['proposal_reference']))
    bands = {'proposal_old': native['frozen_bands']['context']['composite'],
             'proposal_new': {'Ca_max': ref['bands']['Ca_supported_max_R_model_kcal_mol'], 'La_min': ref['bands']['La_supported_min_R_model_kcal_mol']},
             'DFT': {'Ca_max': dr['bands']['Ca_max'], 'La_min': dr['bands']['La_min']}}
    for rep in ('core', 'context'):
        for kind in ('native', 'composite'): bands[rep + '_' + kind] = native['frozen_bands'][rep][kind]
    si = {r['case_id']: r for r in selected['cases']}; ni = {(r['case_id'], r['representation']): r for r in native['rows']}; di = {r['case_id']: r for r in dr['rows']}
    sources = read_json(verify(m['sources'])); groups = source_groups(sources); rows = []
    for src in sources['cases']:
        if not src['primary_evaluation_pool']: continue
        cid = src['case_id']; s = si[cid]; score = s['R_selected']['composite_R_model_kcal_mol'] if s['R_selected'] else None
        values = {'proposal_old': score, 'proposal_new': score, 'DFT': di[cid]['R_kcal_mol']}
        for rep in ('core', 'context'):
            for kind in ('native', 'composite'): values[rep + '_' + kind] = ni[cid, rep][kind + '_R_model_kcal_mol']
        methods = {key: {'R': val, 'decision': decision(val, bands[key]), 'outcome': outcome(decision(val, bands[key]), src['expected_class'])} for key, val in values.items()}
        rows.append({k: src[k] for k in ('case_id', 'root_case_id', 'biological_group', 'expected_class', 'source_conditioning_metal')} | {'methods': methods,
                     'proposal_status': s['status'], 'proposal_delta_R': s['delta_R_model_kcal_mol']})
    index = {r['case_id']: r for r in rows}; pools, triples = [], []
    for root, members in groups.items():
        expected = members[0]['expected_class']; arms = {z: [c['case_id'] for c in members if c['primary_evaluation_pool'] and c['source_conditioning_metal'] == z] for z in ('La', 'Ca')}
        summaries = {z: {method: strict_summary(arms[z], index, method, band, expected) for method, band in bands.items()} for z in arms}
        for z, name in (('La', 'La4'), ('Ca', 'Ca5')):
            pools.append({'root_case_id': root, 'expected_class': expected, 'pool': name, 'members': arms[z], 'methods': summaries[z]})
        balanced = {}
        for method, band in bands.items():
            a, b = summaries['La'][method]['R'], summaries['Ca'][method]['R']
            v = (a + b) / 2 if a is not None and b is not None else None; call = decision(v, band)
            balanced[method] = {'R': v, 'decision': call, 'outcome': outcome(call, expected), 'Ca5_minus_La4': b - a if v is not None else None}
        pools.append({'root_case_id': root, 'expected_class': expected, 'pool': 'balanced', 'members': arms['La'] + arms['Ca'], 'methods': balanced})
        for ids in combinations(arms['La'], 3):
            triples.append({'root_case_id': root, 'expected_class': expected, 'members': list(ids),
                'methods': {method: strict_summary(ids, index, method, band, expected) for method, band in bands.items()}})
    if (len(rows), len(pools), len(triples)) != (225, 75, 100): raise InvalidArtifact('fixed comparison denominators differ')
    counts = {'single': {z: summarize_rows([r for r in rows if z == 'all225' or r['source_conditioning_metal'] == z], bands) for z in ('all225', 'La', 'Ca')},
              'pools': {name: summarize_rows([p for p in pools if p['pool'] == name], bands) for name in ('La4', 'Ca5', 'balanced')},
              'all100_triples': summarize_rows(triples, bands)}
    matched = {}
    subsets = {'all225': rows, 'La100': [r for r in rows if r['source_conditioning_metal'] == 'La'],
               'Ca125': [r for r in rows if r['source_conditioning_metal'] == 'Ca'], 'all100_triples': triples,
               **{name: [p for p in pools if p['pool'] == name] for name in ('La4', 'Ca5', 'balanced')}}
    for name, subset in subsets.items():
        matched[name] = {}
        for comparator in bands:
            if comparator == 'proposal_new': continue
            common = [r for r in subset if all(r['methods'][method]['R'] is not None for method in ('proposal_new', comparator))]
            matched[name][comparator] = {'declared_denominator': len(subset), 'matched_denominator': len(common),
                'counts': summarize_rows(common, ('proposal_new', comparator)),
                'outcome_transitions': dict(Counter(r['methods'][comparator]['outcome'] + '->' + r['methods']['proposal_new']['outcome'] for r in common))}
    result = {'manifest': record(manifest), 'selection': record(selection), 'DFT': record(dft), 'frozen_proposal_reference': m['proposal_reference'],
              'rows': rows, 'pools': pools, 'triples': triples, 'counts': counts, 'matched_comparisons': matched, 'bands': bands,
              'new_molecular_calls': 0, 'threshold_fitted_on_folds': False, 'baseline_changed': False,
              'interpretation': 'consumed protein structural repeats; strict completeness, no best-fold or best-triple selection'}
    write_new(output, result)
    lines = ['# Fixed proposal transfer comparison', '',
             'Consumed reference-protein structural repeats. All 225 sources and every declared pool/triple remain in the denominators. New bands were fixed using canonical25 only; no transfer-set fit.', '',
             '| Population | Method | Correct | Wrong | Inconclusive | Unavailable | Total |',
             '|---|---|---:|---:|---:|---:|---:|']
    tables = {**counts['single'], **counts['pools'], 'all100_triples': counts['all100_triples']}
    for name, table in tables.items():
        for method, tally in table.items():
            lines.append('| ' + name + ' | ' + method + ' | ' + ' | '.join(str(tally[k]) for k in ('correct', 'wrong', 'inconclusive', 'unavailable', 'denominator')) + ' |')
    lines += ['', 'Pairwise complete-case counts and transitions, all raw scores, missing members, within-protein ranges and exact source/result/reference pins are retained in the JSON. No fallback, best-fold selection, independent biological replication or production promotion is implied.', '']
    with Path(output).with_suffix('.md').open('x') as handle: handle.write('\n'.join(lines))
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='operation', required=True)
    p = sub.add_parser('prepare')
    for key in ('design', 'comparison', 'reference', 'agreement', 'output', 'cpu-python', 'gpu-python'): p.add_argument('--' + key, required=True)
    p = sub.add_parser('compare')
    for key in ('manifest', 'selection', 'dft', 'output'): p.add_argument('--' + key, required=True)
    args = vars(parser.parse_args()); op = args.pop('operation'); print(json.dumps(globals()[op](**args)))


if __name__ == '__main__': main()
