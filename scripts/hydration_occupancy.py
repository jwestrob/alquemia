"""Finite, paired water-arrangement experiment; no assumed occupancy free energy."""
from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path
import shutil

import numpy as np

from affordable_common import (HA_TO_KCAL, BOHR_TO_A, InvalidArtifact, cache_key,
                               read_json, record, verify, write_new, xyz)
from affordable_response import read_engrad
from hydration_adjudicate import METHOD
from hydration_mace import rotational_gradient
from hydration_network import geometry_checks
from hydration_square import endpoint
from hydration_state_analysis import exchange_ledger

PROTOCOL = 'native_r2scan3c_cpcm_mace_prepared_water_arrangements_v1'
METALS = ('Ca', 'La')


def pattern_set(n):
    return {''.join(bits) for bits in product('01', repeat=n)}


def physical_identity(preparation):
    return {k: preparation[k] for k in
            ('case', 'source_structure', 'atom_graph', 'charge_ledger',
             'water_groups', 'atoms', 'charges')}


def check_state(task, preparation):
    """Every nonmobile coordinate matches its source, including outer waters."""
    pattern = task['pattern']
    variable = [w for w in preparation['water_groups'] if w['role'] == 'variable']
    if pattern not in pattern_set(len(variable)):
        raise InvalidArtifact('unsupported occupancy pattern')
    omitted = {i for bit, water in zip(pattern, variable) if bit == '0'
               for i in water['indices']}
    keep = [i for i in range(len(preparation['atoms'])) if i not in omitted]
    remap = {old: new for new, old in enumerate(keep)}
    groups = [dict(w, indices=[remap[i] for i in w['indices']],
                   oxygen_index=remap[w['oxygen_index']],
                   hydrogen_indices=[remap[i] for i in w['hydrogen_indices']])
              for w in preparation['water_groups'] if w['oxygen_index'] not in omitted]
    mobile = sorted(i for w in groups if w['role'] == 'variable' for i in w['hydrogen_indices'])
    expected = [tuple(preparation['atoms'][i]) for i in keep]
    expected[0] = (task['metal'], *expected[0][1:])
    actual = xyz(verify(task['xyz']))
    if (task['charge'] != preparation['charges'][task['metal']] or
            task['multiplicity'] != 1 or task['groups'] != groups or
            task['mobile_indices'] != mobile):
        raise InvalidArtifact('water inventory/charge/mobile coordinates changed')
    frozen, internal = geometry_checks(expected, actual, groups, mobile)
    if frozen > 1e-9 or internal > 1e-9:
        raise InvalidArtifact('exact fixed geometry or rigid water shape changed')
    text = verify(task['input']).read_text()
    if text != METHOD + f"\n* xyzfile {task['charge']} 1 core.xyz\n":
        raise InvalidArtifact('DFT Hamiltonian changed')
    return keep


def prepare(states, proposals, full_collection, agreement, output):
    from hydration_proposal_opt import collect as collect_proposals
    source = read_json(states)
    proposal = read_json(proposals)
    pm = read_json(verify(proposal['manifest']))
    if (collect_proposals(verify(proposal['manifest'])) != proposal or
            not proposal['all_optimization_eligible'] or
            not pm.get('occupancy_enumerated') or pm['source_manifests'] != [record(states)]):
        raise InvalidArtifact('complete matched occupancy proposals required')
    full = read_json(full_collection)
    fm = read_json(verify(full['manifest']))
    if full['status'] != 'complete' or fm['orca'] != source['orca']:
        raise InvalidArtifact('compatible full-water endpoints required')
    full_mace = read_json(verify(read_json(verify(fm['mace_collection']))['manifest']))
    old_preparations = {}
    for pin in full_mace['source_manifests']:
        for pp in read_json(verify(pin))['preparations']:
            p = read_json(verify(pp)); old_preparations[p['case']] = p
    preparations = {p['case']: p for p in
                    (read_json(verify(pin)) for pin in source['preparations'])}
    for case, p in preparations.items():
        if physical_identity(p) != physical_identity(old_preparations[case]):
            raise InvalidArtifact('full-water cache physical context mismatch')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    rows = {r['task_id']: r for r in proposal['rows']}
    tasks = []
    for s in source['tasks']:
        if s['seed'] != 'source':
            continue
        tid = f"{s['case']}__{s['pattern']}__{s['metal']}"
        r = rows.get(tid)
        if s['variable_water_count'] and r is None:
            raise InvalidArtifact('missing nonempty-state proposal')
        selected = r['proposed_xyz'] if r else s['xyz']
        directory = out/'tasks'/tid; directory.mkdir(parents=True)
        xp = directory/'core.xyz'; xp.write_bytes(verify(selected).read_bytes())
        ip = directory/'endpoint.inp'
        ip.write_text(METHOD + f"\n* xyzfile {s['charge']} 1 core.xyz\n")
        t = {k: s[k] for k in ('case', 'pattern', 'metal', 'charge', 'multiplicity',
                               'groups', 'mobile_indices', 'variable_water_count')}
        t.update(task_id=tid, input=record(ip), xyz=record(xp),
                 output_path=str(directory/'endpoint.out'),
                 engrad_path=str(directory/'endpoint.engrad'),
                 selected_seed=r['selected_seed'] if r else 'empty_no_optimization',
                 mace_result=r['result'] if r else None, source_xyz=s['xyz'])
        t['retained_full_indices'] = check_state(t, preparations[s['case']])
        t['cache_key'] = cache_key(t | {'protocol_id': PROTOCOL, 'states': record(states)})
        tasks.append(t)
    reused = []
    for t in fm['tasks']:
        p = preparations[t['case']]
        pattern = '1' * sum(w['role'] == 'variable' for w in p['water_groups'])
        r = next(r for r in full['rows'] if r['task_id'] == t['task_id'])
        current = endpoint(r['result']['output'], r['result']['receipt'], t['xyz'], t['input'])
        if current != r['result']:
            raise InvalidArtifact('full-water parsed endpoint changed')
        copied = dict(t, pattern=pattern, variable_water_count=len(pattern),
                      task_id=f"{t['case']}__{pattern}__{t['metal']}",
                      reused_from=record(full_collection), reused_row_id=r['task_id'])
        copied['retained_full_indices'] = check_state(copied, p)
        reused.append(copied)
    expected = {(case, pattern, metal) for case, p in preparations.items()
                for pattern in pattern_set(sum(w['role'] == 'variable' for w in p['water_groups']))
                for metal in METALS}
    keys = [(t['case'], t['pattern'], t['metal']) for t in tasks + reused]
    if len(keys) != len(set(keys)) or set(keys) != expected:
        raise InvalidArtifact('occupancy experiment is incomplete or duplicated')
    impl = out/'implementation'; impl.mkdir(); pins = {}
    for name in ('hydration_occupancy.py', 'hydration_adjudicate.py', 'hydration_state_analysis.py',
                 'hydration_network.py', 'hydration_mace.py', 'hydration_square.py',
                 'affordable_peptide.py', 'affordable_response.py', 'affordable_common.py',
                 'mace_hybrid.py'):
        path = impl/name; shutil.copyfile(Path(__file__).with_name(name), path); pins[name] = record(path)
    manifest = {'protocol_id': PROTOCOL, 'tasks': tasks, 'reused': reused,
                'states_manifest': record(states), 'proposal_collection': record(proposals),
                'full_collection': record(full_collection), 'preparations': source['preparations'],
                'agreement': record(agreement), 'orca': source['orca'],
                'execution_policy': source['execution_policy'], 'implementation': pins,
                'execution_resources': {'mpi_ranks': 16, 'concurrent_tasks': 4},
                'new_endpoint_count': len(tasks), 'reused_endpoint_count': len(reused),
                'baseline_changed': False, 'occupancy_probabilities': None}
    write_new(out/'manifest.json', manifest)
    from affordable_workflow import dry_run
    return dry_run(out/'manifest.json')


def collect(manifest, output):
    import gemmi
    m = read_json(manifest)
    verify(m['agreement'])
    preparations = {p['case']: p for p in (read_json(verify(pin)) for pin in m['preparations'])}
    rows = []
    for t in m['tasks'] + m['reused']:
        row = {k: t[k] for k in ('task_id', 'case', 'pattern', 'metal', 'variable_water_count')}
        row['reused'] = 'reused_from' in t
        try:
            check_state(t, preparations[t['case']])
            op = Path(t['output_path']); rp = Path(str(op)+'.execution.json')
            result = endpoint(record(op), record(rp), t['xyz'], t['input'])
            gp = verify(read_json(rp)['artifacts']['engrad']); g = read_engrad(gp)
            atoms = xyz(verify(t['xyz'])); coords = np.array([a[1:] for a in atoms])
            if (not np.array_equal(g['atomic_numbers'], [gemmi.Element(a[0]).atomic_number for a in atoms]) or
                    abs(g['energy_Ha']-result['energy_hartree']) > 1e-7 or
                    not np.allclose(g['coordinates_bohr']*BOHR_TO_A, coords, atol=2e-6, rtol=0)):
                raise InvalidArtifact('native gradient/energy/geometry mismatch')
            torques = rotational_gradient(coords, g['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A, t['groups'])
            row.update(status='complete', result=result, gradient=record(gp),
                       selected_seed=t['selected_seed'],
                       rotational_gradient_kcal_mol_radian=torques.tolist(),
                       rotational_gradient_max_kcal_mol_radian=float(np.max(abs(torques))) if torques.size else None)
        except (OSError, ValueError, KeyError) as exc:
            row.update(status='unavailable', failure=str(exc))
        rows.append(row)
    result = {'manifest': record(manifest), 'rows': rows,
              'status': 'complete' if all(r['status'] == 'complete' for r in rows) else 'incomplete',
              'baseline_changed': False, 'occupancy_probabilities': None, 'implementation': record(__file__)}
    write_new(output, result)
    return {'status': result['status'], 'completed': sum(r['status'] == 'complete' for r in rows), 'total': len(rows)}


def analyze(collection, reference, output):
    c = read_json(collection); m = read_json(verify(c['manifest'])); ref = read_json(reference)
    if ref['status'] != 'approximate_bulk_reference_available':
        raise InvalidArtifact('computed bulk-water reference required')
    for e in ref['endpoints'].values():
        for key in ('output', 'receipt'):
            verify(e[key])
    expected = {t['task_id'] for t in m['tasks'] + m['reused']}
    source_tasks = {t['task_id']: t for t in m['tasks'] + m['reused']}
    if len(c['rows']) != len(expected) or {r['task_id'] for r in c['rows']} != expected:
        raise InvalidArtifact('missing/duplicate endpoint rows')
    cases = []
    for pp in m['preparations']:
        p = read_json(verify(pp)); case = p['case']
        waters = [w['source'] for w in p['water_groups'] if w['role'] == 'variable']
        states = []
        for pattern in sorted(pattern_set(len(waters))):
            pair = {r['metal']: r for r in c['rows'] if r['case'] == case and r['pattern'] == pattern}
            complete = set(pair) == set(METALS) and all(r['status'] == 'complete' for r in pair.values())
            endpoints = {}
            for metal, r in pair.items():
                if r['status'] == 'complete':
                    task = source_tasks[r['task_id']]
                    actual = endpoint(r['result']['output'], r['result']['receipt'], task['xyz'], task['input'])
                    if actual != r['result']:
                        raise InvalidArtifact('collected endpoint differs from its scientific output')
                    endpoints[metal] = r['result']
            states.append({'pattern': pattern, 'n_variable_waters': pattern.count('1'),
                           'occupied_waters': [w for bit, w in zip(pattern, waters) if bit == '1'],
                           'status': 'complete' if complete else 'unavailable', 'endpoints': endpoints,
                           'R_hartree': endpoints['Ca']['energy_hartree']-endpoints['La']['energy_hartree'] if complete else None})
        by_pattern = {s['pattern']: s for s in states}
        full = by_pattern['1'*len(waters)]; empty = by_pattern['0'*len(waters)]
        edges = []; nonadditivity = []; fixed_count = []
        for s in states:
            if s['status'] != 'complete':
                continue
            if full['status'] == 'complete':
                s['delta_R_from_full_kcal_mol'] = (s['R_hartree']-full['R_hartree'])*HA_TO_KCAL
            if empty['status'] == 'complete':
                s['exchange_ledger'] = {metal: exchange_ledger(s['endpoints'][metal]['energy_hartree'],
                    empty['endpoints'][metal]['energy_hartree'], s['n_variable_waters'], ref) for metal in METALS}
            for i, bit in enumerate(s['pattern']):
                if bit == '1':
                    continue
                to = s['pattern'][:i]+'1'+s['pattern'][i+1:]; dest = by_pattern[to]
                if dest['status'] != 'complete':
                    continue
                edges.append({'from': s['pattern'], 'to': to, 'water': waters[i],
                              'delta_R_water_addition_kcal_mol': (dest['R_hartree']-s['R_hartree'])*HA_TO_KCAL,
                              'exchange_ledger': {metal: exchange_ledger(dest['endpoints'][metal]['energy_hartree'],
                                  s['endpoints'][metal]['energy_hartree'], 1, ref) for metal in METALS}})
            if s['n_variable_waters'] >= 2 and empty['status'] == 'complete':
                singles = [by_pattern['0'*i+'1'+'0'*(len(waters)-i-1)] for i, bit in enumerate(s['pattern']) if bit == '1']
                if all(a['status'] == 'complete' for a in singles):
                    terms = {metal: (s['endpoints'][metal]['energy_hartree'] + (len(singles)-1)*empty['endpoints'][metal]['energy_hartree']
                             - sum(a['endpoints'][metal]['energy_hartree'] for a in singles))*HA_TO_KCAL for metal in METALS}
                    nonadditivity.append({'pattern': s['pattern'], 'electronic_nonadditivity_kcal_mol': terms,
                                          'nonadditivity_in_R_kcal_mol': terms['Ca']-terms['La']})
        for n in range(len(waters)+1):
            subset = [s for s in states if s['n_variable_waters'] == n]
            if not all(s['status'] == 'complete' for s in subset):
                continue
            chosen = {metal: min(subset, key=lambda s: s['endpoints'][metal]['energy_hartree']) for metal in METALS}
            fixed_count.append({'n_variable_waters': n, 'patterns': {z: s['pattern'] for z, s in chosen.items()},
                                'R_hartree': chosen['Ca']['endpoints']['Ca']['energy_hartree']-chosen['La']['endpoints']['La']['energy_hartree'],
                                'interpretation': 'conditional electronic minima over enumerated preparations, not occupancy predictions'})
        cases.append({'case': case, 'variable_water_order': waters, 'states': states,
                      'complete_occupancy_table': all(s['status'] == 'complete' for s in states),
                      'water_addition_edges': edges, 'nonadditivity': nonadditivity,
                      'fixed_count_electronic_contrasts': fixed_count})
    result = {'protocol_id': PROTOCOL, 'collection': record(collection), 'reference': record(reference),
              'cases': cases, 'status': 'complete' if all(c['complete_occupancy_table'] for c in cases) else 'incomplete',
              'occupancy_probabilities': None, 'formation_free_energies': None, 'calibrated_classification': None,
              'baseline_changed': False, 'implementation': record(__file__)}
    write_new(output, result)
    return {'status': result['status'], 'patterns': sum(len(c['states']) for c in cases),
            'edges': sum(len(c['water_addition_edges']) for c in cases)}


def motion(collection, output):
    """Project saved gradients onto physical water translations and rotations.

    This evaluates no potential and supplies no covariance, Hessian, or entropy.
    MACE is vacuum; DFT includes CPCM. Their difference is deliberately retained.
    """
    from mace_hybrid import EV_TO_KCAL
    c = read_json(collection); m = read_json(verify(c['manifest']))
    tasks = {t['task_id']: t for t in m['tasks']+m['reused']}
    rows = []
    for r in c['rows']:
        t = tasks[r['task_id']]
        if r['status'] != 'complete':
            rows.append({'task_id': r['task_id'], 'status': 'unavailable'})
            continue
        atoms = xyz(verify(t['xyz'])); coords = np.array([a[1:] for a in atoms])
        gradient = read_engrad(verify(r['gradient']))['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A
        mg = None
        if t.get('mace_result'):
            saved = read_json(verify(t['mace_result']))
            run = next(a for a in saved['optimization_runs'] if a['seed'] == t['selected_seed'])
            if xyz(verify(run['xyz'])) != atoms:
                raise InvalidArtifact('MACE/DFT water-motion coordinates differ')
            mg = -np.load(verify(run['forces']), allow_pickle=False)*EV_TO_KCAL
        groups = [w for w in t['groups'] if w['role'] == 'variable']
        torque = rotational_gradient(coords, gradient, groups)
        mt = rotational_gradient(coords, mg, groups) if mg is not None else None
        waters = []
        for i, w in enumerate(groups):
            translation = gradient[w['indices']].sum(axis=0)
            mtg = mg[w['indices']].sum(axis=0) if mg is not None else None
            cosine = (float(translation@mtg/(np.linalg.norm(translation)*np.linalg.norm(mtg)))
                      if mtg is not None and np.linalg.norm(translation)*np.linalg.norm(mtg)>1e-12 else None)
            waters.append({'source': w['source'],
                           'DFT_translation_gradient_kcal_mol_A': translation.tolist(),
                           'MACE_translation_gradient_kcal_mol_A': mtg.tolist() if mtg is not None else None,
                           'translation_gradient_cosine': cosine,
                           'translation_gradient_error_norm_kcal_mol_A': float(np.linalg.norm(translation-mtg)) if mtg is not None else None,
                           'DFT_rotation_gradient_kcal_mol_radian': torque[i].tolist(),
                           'MACE_rotation_gradient_kcal_mol_radian': mt[i].tolist() if mt is not None else None})
        rows.append({'task_id': r['task_id'], 'status': 'complete', 'waters': waters,
                     'DFT_gradient': r['gradient'], 'MACE_result': t.get('mace_result')})
    result = {'collection': record(collection), 'rows': rows, 'implementation': record(__file__),
              'new_energy_evaluations': 0, 'entropy_correction': None, 'curvature': None,
              'interpretation': 'physical projected gradients at prepared geometries, not equilibrium motion or validated curvature'}
    write_new(output, result)
    return {'completed': sum(r['status']=='complete' for r in rows), 'new_energy_evaluations': 0}


def report(result, output):
    """Export a compact table and plain-text report; no new energy evaluation."""
    import csv
    data = read_json(result); out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    text = ['# Joint water-arrangement results', '',
            'Electronic comparisons at fixed protein/water-oxygen coordinates. MACE prepares water H; DFT scores.',
            'Missing bound-state free-energy terms are not zero. No equilibrium occupancies or new classifications.', '',
            '| Structure | Water pattern | Count | R − R(full), kcal/mol |',
            '|---|---|---:|---:|']
    with (out/'states.csv').open('x', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['case', 'pattern', 'n_variable_waters', 'status', 'E_Ca_hartree',
                         'E_La_hartree', 'R_hartree', 'delta_R_from_full_kcal_mol'])
        for case in data['cases']:
            for s in case['states']:
                delta = s.get('delta_R_from_full_kcal_mol')
                text.append(f"| {case['case']} | {s['pattern']} | {s['n_variable_waters']} | {delta:+.3f} |" if delta is not None
                            else f"| {case['case']} | {s['pattern']} | {s['n_variable_waters']} | unavailable |")
                writer.writerow([case['case'], s['pattern'], s['n_variable_waters'], s['status'],
                                 s['endpoints'].get('Ca', {}).get('energy_hartree'),
                                 s['endpoints'].get('La', {}).get('energy_hartree'), s['R_hartree'], delta])
    with (out/'water_additions.csv').open('x', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['case', 'from_pattern', 'to_pattern', 'water', 'delta_R_kcal_mol',
                         'Ca_electronic_exchange_kcal_mol', 'La_electronic_exchange_kcal_mol',
                         'Ca_missing_nonE_for_neutral_kcal_mol', 'La_missing_nonE_for_neutral_kcal_mol',
                         'Ca_remaining_if_gas_internal_ZPE_kcal_mol', 'La_remaining_if_gas_internal_ZPE_kcal_mol'])
        for case in data['cases']:
            for e in case['water_addition_edges']:
                values = [e['exchange_ledger'][z] for z in METALS]
                fields = ('electronic_exchange_from_empty_kcal_mol', 'missing_bound_contribution_for_neutral_exchange_kcal_mol',
                          'remaining_contribution_if_gas_internal_ZPE_is_assumed_kcal_mol')
                writer.writerow([case['case'], e['from'], e['to'], str(e['water']), e['delta_R_water_addition_kcal_mol']]
                                + [v[k] for k in fields for v in values])
    text += ['', '## Pattern identities', '']
    for case in data['cases']:
        names = [f"{w['chain']}:{w['resnum']}{w['insertion_code']}" for w in case['variable_water_order']]
        text.append(f"- {case['case']}: bits refer in order to {', '.join(names)}; 1 means present.")
    text += ['', 'R = E_Ca − E_La. Positive changes are more La-like on this fixed-context scale.',
             'All source/receipt links, unrounded energies, conditional fixed-count minima and exchange thresholds are retained in the result JSON.',
             f'Source: {Path(result).resolve()}', '']
    (out/'TABLE.md').write_text('\n'.join(text))
    # Standalone scientific figure, with every pattern and the actual missingness.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(data['cases']), figsize=(9, 4.8), squeeze=False)
    values = [s['delta_R_from_full_kcal_mol'] for c in data['cases'] for s in c['states']
              if s.get('delta_R_from_full_kcal_mol') is not None] + [0.]
    span = max(max(values)-min(values), 1.)
    for ax, case in zip(axes[0], data['cases']):
        states = sorted(case['states'], key=lambda s: (s['n_variable_waters'], s['pattern']))
        for i, s in enumerate(states):
            value = s.get('delta_R_from_full_kcal_mol')
            if value is not None:
                ax.barh(i, value, color='#276b8e', height=.65)
                ax.plot(value, i, 'o', color='#14384b', markersize=4)
            else:
                ax.text(0, i, ' unavailable', va='center', color='#9c3333')
        ax.set_yticks(range(len(states)), [s['pattern'] for s in states])
        ax.axvline(0, color='#888888', linewidth=.8)
        ax.set_xlim(min(values)-.06*span, max(values)+.06*span)
        ax.set_title(case['case']); ax.set_ylabel('Water presence pattern')
        ax.set_xlabel('R − R(full water inventory), kcal/mol')
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Water arrangements change the electronic La/Ca contrast')
    fig.text(.5, .015, 'Fixed protein and water oxygens; MACE-prepared H; native DFT. These are not occupancy probabilities.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .055, 1, .96))
    fig.savefig(out/'water_arrangements.svg'); fig.savefig(out/'water_arrangements.png', dpi=180)
    plt.close(fig)
    write_new(out/'source.json', {'result': record(result), 'implementation': record(__file__)})
    return {'status': data['status'], 'output': str(out)}


if __name__ == '__main__':
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='op', required=True)
    p = sub.add_parser('prepare')
    for name in ('states', 'proposals', 'full-collection', 'agreement', 'output'):
        p.add_argument('--'+name, required=True)
    p = sub.add_parser('collect')
    for name in ('manifest', 'output'):
        p.add_argument('--'+name, required=True)
    p = sub.add_parser('analyze')
    for name in ('collection', 'reference', 'output'):
        p.add_argument('--'+name, required=True)
    p = sub.add_parser('motion')
    for name in ('collection', 'output'):
        p.add_argument('--'+name, required=True)
    p = sub.add_parser('report')
    for name in ('result', 'output'):
        p.add_argument('--'+name, required=True)
    args = vars(parser.parse_args()); op = args.pop('op')
    print(json.dumps(globals()[op](**args), indent=2))
