"""Compare actual hydration states without inventing bound-water free energies."""
from __future__ import annotations
import argparse
from math import comb
from affordable_common import HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new


def exchange_ledger(full_energy, empty_energy, n, reference):
    """Report the missing state contribution required for neutral formation.

    DeltaG = (E_state-E_empty-n*Egas) + DeltaG_bound_nonE
             - n*(mu_water-Egas). DeltaG_bound_nonE is unknown, not zero.
    """
    egas=reference['endpoints']['gas']['energy_hartree']
    terms=reference['terms']
    electronic=(full_energy-empty_energy-n*egas)*HA_TO_KCAL
    required=n*terms['mu_liquid_minus_Egas_kcal_mol']-electronic
    return {'n_variable_waters':n,'electronic_exchange_from_empty_kcal_mol':electronic,
            'missing_bound_contribution_for_neutral_exchange_kcal_mol':required,
            'remaining_contribution_if_gas_internal_ZPE_is_assumed_kcal_mol':required-n*terms['harmonic_gas_ZPE_kcal_mol'],
            'bound_state_correction_kcal_mol':None,'formation_free_energy_kcal_mol':None,
            'interpretation':'threshold for an uncomputed state contribution; not an assigned correction or occupancy'}


def analyze(collections, output, reference=None, previous=None):
    gathered=[];mappings={};manifests=[];signature=None
    for path in collections:
        c=read_json(path);m=read_json(verify(c['manifest']));manifests.append(c['manifest'])
        cfg=read_json(verify(m['configuration']))
        current=(m['protocol_id'],m['configuration']['sha256'])
        if signature is None:signature=current
        elif current!=signature:raise InvalidArtifact('incompatible hydration protocols/configurations')
        for pin in m['preparations']:
            p=read_json(verify(pin))
            fields=('atom_graph','charge_ledger','water_groups','atoms','source_structure')
            if p['case'] in mappings and cache_key({k:p[k] for k in fields})!=cache_key({k:mappings[p['case']][k] for k in fields}):
                raise InvalidArtifact('case physical representation changed between collections')
            mappings[p['case']]=p
        by_id={t['task_id']:t for t in m['tasks']}
        if len(c['rows'])!=len(by_id) or {r['task_id'] for r in c['rows']}!=set(by_id):
            raise InvalidArtifact('collection does not contain exactly the manifested tasks')
        for row in c['rows']:
            t=by_id[row['task_id']]
            if row['status']=='complete':
                for name in ('output','receipt'):
                    verify(row['result'][name])
            gathered.append(dict(row,task=t))
    ids=[r['task_id'] for r in gathered]
    if len(ids)!=len(set(ids)):raise InvalidArtifact('duplicate states/starts from overlapping collections')
    ref=read_json(reference) if reference else None
    if ref and ref['status']!='approximate_bulk_reference_available':
        raise InvalidArtifact('unsupported water reference')
    old=read_json(previous) if previous else None
    cases=[]
    for case in sorted(mappings):
        rows=[r for r in gathered if r['case']==case];patterns=sorted({r['pattern'] for r in rows})
        states=[]
        for pattern in patterns:
            endpoints={};failures=[]
            for metal in ('Ca','La'):
                seeds=[r for r in rows if r['pattern']==pattern and r['metal']==metal]
                expected_seeds=set(cfg['orientation_seeds']) if '1' in pattern else {'source'}
                if {r['seed'] for r in seeds}!=expected_seeds:
                    failures.append(metal+':incomplete_predefined_seeds')
                    continue
                good=[r for r in seeds if r['status']=='complete']
                failed=[r for r in seeds if r['status']!='complete']
                if failed or not good:
                    failures.extend([r['task_id'] for r in failed] or [metal+':not_run'])
                    continue
                winner=min(good,key=lambda r:r['result']['energy_hartree'])
                energies=[r['result']['energy_hartree'] for r in good]
                endpoints[metal]={'energy_hartree':winner['result']['energy_hartree'],
                    'selected_task':winner['task_id'],'selected_seed':winner['seed'],
                    'seed_spread_kcal_mol':(max(energies)-min(energies))*HA_TO_KCAL,
                    'all_seeds':[{'task':r['task_id'],'seed':r['seed'],'result':r['result']} for r in good]}
            complete=not failures
            states.append({'pattern':pattern,'n_variable_waters':pattern.count('1'),
                'status':'complete' if complete else 'unavailable','failures':failures,
                'endpoints':endpoints,'R_hartree':endpoints['Ca']['energy_hartree']-endpoints['La']['energy_hartree'] if complete else None})
        full=next((s for s in states if set(s['pattern'])=={'1'}),None)
        expected=2**sum(w['role']=='variable' for w in mappings[case]['water_groups'])
        summary={'case':case,'states':states,'full_water_orientation_comparison':None,
                 'expected_occupancy_patterns':expected,'observed_occupancy_patterns':len(states),
                 'complete_occupancy_table':len(states)==expected and all(s['status']=='complete' for s in states),
                 'composition':{'atoms':mappings[case]['atom_count'],'charges':mappings[case]['charges']}}
        if full and full['status']=='complete':
            initial={}
            for metal in ('Ca','La'):
                source=next(r for r in rows if r['pattern']==full['pattern'] and r['metal']==metal and r['seed']=='source')
                initial[metal]=source['result']['initial_energy_hartree']
            initial_r=initial['Ca']-initial['La']
            summary['full_water_orientation_comparison']={'initial_source_R_hartree':initial_r,
                'minimum_found_R_hartree':full['R_hartree'],
                'water_reorientation_delta_R_kcal_mol':(full['R_hartree']-initial_r)*HA_TO_KCAL,
                'endpoint_lowerings_kcal_mol':{metal:(full['endpoints'][metal]['energy_hartree']-initial[metal])*HA_TO_KCAL for metal in ('Ca','La')},
                'selection_rule':'lowest electronic energy per metal among the two predefined converged starts',
                'missing_free_energy_terms':True}
            if old:
                prior=next(r for r in old['rows'] if r['case']==case and r['removed_water'] is None)
                summary['full_water_orientation_comparison']['representation_plus_SCF_policy_delta_R_kcal_mol']=(initial_r-prior['R_hartree'])*HA_TO_KCAL
        by_count=[]
        for n in sorted({s['n_variable_waters'] for s in states}):
            subset=[s for s in states if s['n_variable_waters']==n]
            water_count=sum(w['role']=='variable' for w in mappings[case]['water_groups'])
            if len(subset)!=comb(water_count,n) or any(s['status']!='complete' for s in subset):continue
            winners={metal:min(subset,key=lambda s:s['endpoints'][metal]['energy_hartree']) for metal in ('Ca','La')}
            by_count.append({'n_variable_waters':n,'patterns':{metal:w['pattern'] for metal,w in winners.items()},
                'R_hartree':winners['Ca']['endpoints']['Ca']['energy_hartree']-winners['La']['endpoints']['La']['energy_hartree'],
                'interpretation':'conditional fixed-water-count electronic minima, not predicted water count'})
        summary['fixed_count_electronic_contrasts']=by_count
        empty=next((s for s in states if s['n_variable_waters']==0),None)
        if ref and empty and empty['status']=='complete':
            for state in states:
                if state['status']=='complete':
                    state['exchange_ledger']={metal:exchange_ledger(state['endpoints'][metal]['energy_hartree'],
                        empty['endpoints'][metal]['energy_hartree'],state['n_variable_waters'],ref) for metal in ('Ca','La')}
        cases.append(summary)
    result={'collections':[record(p) for p in collections],'manifests':manifests,
            'reference':record(reference) if reference else None,'previous':record(previous) if previous else None,
            'cases':cases,'status':'complete' if gathered and all(s['status']=='complete' for c in cases for s in c['states']) else 'incomplete',
            'occupancy_probabilities':None,'calibrated_classification':None,'baseline_changed':False,
            'implementation':record(__file__)}
    write_new(output,result)
    return {'status':result['status'],'cases':[{'case':c['case'],'full_water_orientation_comparison':c['full_water_orientation_comparison']} for c in cases]}


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--collections',nargs='+',required=True)
    p.add_argument('--output',required=True);p.add_argument('--reference');p.add_argument('--previous')
    print(json.dumps(analyze(**vars(p.parse_args())),indent=2))
