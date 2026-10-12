#!/usr/bin/env python3
"""Independent scalar/receipt review; no imports from scientific executors."""
import argparse
import hashlib
import json
import re
from pathlib import Path

EV_KCAL = 23.06054783061903
EH_KCAL = 627.509474
ORDER = ['origin', 'original_adaptive_Ca', 'original_adaptive_La',
         'expanded_adaptive_Ca', 'expanded_adaptive_La']
POOLS = {'origin': ORDER[:1], 'original_three': ORDER[:3],
         'expanded_three': [ORDER[0], *ORDER[3:]], 'common_five': ORDER}


def pin(p):
    p = Path(p)
    return {'path': str(p.resolve()), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def check(p):
    assert pin(p['path'])['sha256'] == p['sha256'], p['path']
    return Path(p['path'])


def read(p):
    return json.loads(Path(p).read_text())


def xyz(p):
    lines = Path(p).read_text().splitlines()
    assert int(lines[0]) == len(lines)-2
    entries = [l.split() for l in lines[2:]]
    assert all(len(a) == 4 for a in entries)
    return [a[0] for a in entries], [[float(x) for x in a[1:]] for a in entries]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for key in ['mace', 'solvent', 'output']:
        ap.add_argument('--'+key, type=Path, required=True)
    ap.add_argument('--comparison', type=Path)
    args = ap.parse_args()
    assert not args.output.exists()
    m, s = read(args.mace), read(args.solvent)
    assert m['status'] == 'complete'
    high = {x['task_id']: x for x in m['rows']}
    assert len(high) == len(s['rows']) == 40
    check(m['manifest']); check(s['manifest'])
    manifest = read(s['manifest']['path'])
    component_tol = manifest['settings']['component_parity_Eh']
    population_tol = manifest['settings']['population_parity_e']
    charge_tol = manifest['settings']['charge_sum_tolerance_e']
    endpoints, executions, geometry_matches, physical_pairs = {}, {}, [], {}
    max_component_error = max_population_error = 0.0
    failed_absolute_inputs = 0
    for low in s['rows']:
        h = high[low['id']]
        assert h['charge'] == low['charge'] and h['multiplicity'] == low['uhf']+1
        a, b = xyz(check(h['xyz'])), xyz(check(low['xyz']))
        assert a[0] == b[0]
        delta = max(abs(x-y) for ca, cb in zip(a[1],b[1]) for x,y in zip(ca,cb))
        assert delta <= 1e-12
        geometry_matches.append({'id': low['id'], 'mace_xyz': h['xyz'], 'solvent_xyz': low['xyz'],
                                 'maximum_component_error_A': delta})
        physical_pairs[h['protein'],h['target_context'],h['geometry_id'],h['metal']] = a
        mr = read(check(h['receipt']))
        assert mr['energy_eV'] == h['energy_eV']
        values = {}
        for medium, rp in low['inputs'].items():
            receipt = read(check(rp))
            assert receipt['task']['xyz'] == low['xyz']
            assert receipt['task']['charge'] == low['charge']
            assert receipt['task']['multiplicity'] == h['multiplicity']
            values[medium] = {}
            for backend, value in receipt['backends'].items():
                assert value['status'] == 'available'
                assert len(value['charges_e']) == len(a[0])
                assert abs(sum(value['charges_e'])-low['charge']) <= charge_tol
                exe = read(check(value['execution']))
                assert exe['returncode'] == 0 and exe['owned_process_group_empty']
                check(exe['stdout']); check(exe['stderr'])
                stdout = Path(exe['stdout']['path']).read_text()
                pattern = r'TOTAL ENERGY\s+([-+0-9.Ee]+)\s+Eh' if backend == 'xtb' else r'^total energy\s+([-+0-9.Ee]+)\s+Eh'
                raw_energy = float(re.findall(pattern,stdout,re.M)[-1])
                assert abs(raw_energy-value['energy_hartree']) < 1e-10
                if backend == 'xtb':
                    assert 'convergence criteria satisfied' in stdout
                executions[value['execution']['path']] = {'backend': backend, **exe}
                values[medium][backend] = value
            de = values[medium]['xtb']['energy_hartree']-values[medium]['tblite']['energy_hartree']
            dq = max(abs(x-y) for x,y in zip(values[medium]['xtb']['charges_e'],values[medium]['tblite']['charges_e']))
            max_population_error = max(max_population_error,dq)
            assert dq <= population_tol
            failed_absolute_inputs += abs(de) > manifest['settings']['energy_parity_Eh']
        transfer = {b: values['water'][b]['energy_hartree']-values['vacuum'][b]['energy_hartree']
                    for b in ['xtb','tblite']}
        assert transfer['xtb'] == low['primary_solvation_hartree']
        error = abs(transfer['xtb']-transfer['tblite'])
        max_component_error = max(max_component_error,error)
        assert error <= component_tol and low['accepted']
        k = (h['protein'], h['target_context'], h['geometry_id'],h['metal'])
        endpoints[k] = {'MACE': h['energy_eV']*EV_KCAL, 'solvent': transfer['xtb']*EH_KCAL}
        endpoints[k]['model'] = sum(endpoints[k].values())
    assert len(executions) == 160 and failed_absolute_inputs == 80
    proteins = sorted({k[0] for k in endpoints})
    contrasts, selections, interactions = {}, [], []
    for protein in proteins:
        for context in ['original','expanded']:
            for g in ORDER:
                pc, pl = [physical_pairs[protein,context,g,metal] for metal in ['Ca','La']]
                assert pc[0][0] == 'Ca' and pl[0][0] == 'La' and pc[0][1:] == pl[0][1:]
                assert pc[1] == pl[1], 'metal endpoint coordinates differ'
                ca, la = endpoints[protein,context,g,'Ca'], endpoints[protein,context,g,'La']
                contrasts[protein,context,g] = {k: ca[k]-la[k] for k in ['MACE','solvent','model']}
            for pool, candidates in POOLS.items():
                entry = {'protein': protein, 'context': context, 'pool': pool}
                for policy in ['mathematical','operational']:
                    chosen = {}
                    for metal in ['Ca','La']:
                        def energy(g): return endpoints[protein,context,g,metal]['model']
                        best = min(candidates,key=energy)
                        if policy == 'operational' and energy('origin')-energy(best) <= .1:
                            best = 'origin'
                        chosen[metal] = {'geometry':best, 'work_kcal_mol':energy(best)-energy('origin')}
                    gc, gl = chosen['Ca']['geometry'], chosen['La']['geometry']
                    chosen['R_model_kcal_mol'] = endpoints[protein,context,gc,'Ca']['model']-endpoints[protein,context,gl,'La']['model']
                    entry[policy] = chosen
                selections.append(entry)
        for g in ORDER:
            d = {'protein':protein,'geometry':g}
            for component in ['MACE','solvent','model']:
                original, expanded = [contrasts[protein,c,g][component] for c in ['original','expanded']]
                o0,e0 = [contrasts[protein,c,'origin'][component] for c in ['original','expanded']]
                d[component] = {'context_effect_kcal_mol':expanded-original,
                                'motion_original_kcal_mol':original-o0,'motion_expanded_kcal_mol':expanded-e0,
                                'interaction_kcal_mol':(expanded-e0)-(original-o0)}
            interactions.append(d)
    anchor_ordering = []
    lookup = {(x['protein'],x['context'],x['pool']):x for x in selections}
    for context in ['original','expanded']:
        for pool in POOLS:
            ca,la = [lookup[p,context,pool]['operational']['R_model_kcal_mol']
                     for p in ['q9z4j7-pqq-la_model','c5axv8-pqq-la_model']]
            anchor_ordering.append({'context':context,'pool':pool,'La_minus_Ca_kcal_mol':la-ca})
    compared = None
    if args.comparison:
        supplied = read(args.comparison)
        for p in supplied['pools']:
            own = lookup[p['protein'],p['context'],p['pool']]
            for policy in ['mathematical','operational']:
                assert abs(p[policy]['R_model_kcal_mol']-own[policy]['R_model_kcal_mol']) < 1e-7
                for metal in ['Ca','La']:
                    assert p[policy][metal]['geometry'] == own[policy][metal]['geometry']
                    assert abs(p[policy][metal]['work_kcal_mol']-own[policy][metal]['work_kcal_mol']) < 1e-7
        for r in supplied['contrasts']:
            own = contrasts[r['protein'],r['context'],r['geometry']]
            assert all(abs(r['R_'+k+'_kcal_mol']-own[k]) < 1e-7 for k in own)
        own_effects = {(r['protein'],r['geometry']):r for r in interactions}
        for r in supplied['effects']:
            own = own_effects[r['protein'],r['geometry']]
            for component in ['MACE','solvent','model']:
                for theirs,ours in [('context_effect','context_effect'),('motion_work_original','motion_original'),
                                    ('motion_work_expanded','motion_expanded'),('interaction','interaction')]:
                    assert abs(r[component+'_'+theirs+'_kcal_mol']-own[component][ours+'_kcal_mol']) < 1e-7
        own_anchors = {(r['context'],r['pool']):r['La_minus_Ca_kcal_mol'] for r in anchor_ordering}
        for r in supplied['anchor_ordering']:
            assert abs(r['two_anchor_La_minus_Ca_kcal_mol']-own_anchors[r['context'],r['pool']]) < 1e-7
        compared = pin(args.comparison)
    costs = {}
    for backend in ['xtb','tblite']:
        rows = [r for r in executions.values() if r['backend'] == backend]
        costs[backend] = {'calls':len(rows),'process_CPU_seconds':sum(r['cpu_user_seconds']+r['cpu_system_seconds'] for r in rows),
                          'summed_worker_wall_seconds':sum(r['wall_seconds'] for r in rows),
                          'peak_process_RSS_KiB':max(r['peak_RSS_KiB'] for r in rows)}
    result = {'schema':'independent_context_factorial_review.v1','new_scientific_calls':0,
              'sources':{'MACE':pin(args.mace),'solvent':pin(args.solvent)},'integrated_report_checked':compared,
              'declared_components':40,'accepted_components':40,'absolute_gate_failed_inputs':failed_absolute_inputs,
              'max_component_error_Eh':max_component_error,'max_population_error_e':max_population_error,
              'geometry_matches':geometry_matches,'selections':selections,'interactions':interactions,'anchor_ordering':anchor_ordering,
              'actual_factorial_execution_costs':costs,'MACE_cost_as_recorded':m['cost'],
              'limits':'Two consumed anchors; no new calibration, class claims, affinity or force qualification.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['geometry_matches','interactions','selections']},indent=2))


if __name__ == '__main__':
    main()
