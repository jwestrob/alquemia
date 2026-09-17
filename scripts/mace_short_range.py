"""Read-only assessment of the predeclared MACE short-range learned component."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from mace_hybrid import accepted_attempt, EV_TO_KCAL
from mace_global_prepare import CASES, LABELS
from mace_global_benchmark import CONTRASTS

PROTOCOL='mace_polar_interaction_energy_descriptor_v1'


def accepted_source(path):
    c=read_json(path);mp=verify(c['manifest']);m=read_json(mp)
    if c['status']!='complete':raise InvalidArtifact('complete actual MACE collection required')
    for t in m['tasks']:
        rows=[accepted_attempt(a,t,mp) for a in (mp.parent/'execution'/t['task_id']).glob('attempt_*')]
        if c['rows'][t['task_id']] not in rows:raise InvalidArtifact('saved output lacks matching successful receipt')
    return c,m


def paired_components(la,ca):
    required={'interaction_energy','electron_energy','electrostatic_energy'}
    if set(la['energy_components_eV'])!=required or set(ca['energy_components_eV'])!=required:
        raise InvalidArtifact('missing or additional energy components require an accounting revision')
    components={key:(ca['energy_components_eV'][key]-la['energy_components_eV'][key])*EV_TO_KCAL for key in required}
    total=(ca['energy_eV']-la['energy_eV'])*EV_TO_KCAL
    return {'R_short_kcal_mol':components['interaction_energy'],'components_kcal_mol':components,
            'R_total_vacuum_kcal_mol':total,'constant_atomic_reference_contrast_kcal_mol':total-sum(components.values()),
            'reference':None,'S_kcal_mol':None,'calibrated_class':None,'short_component_gradient':None}


def assess(medium,large,local,plan,output):
    start=time.monotonic();sources={label:accepted_source(path) for label,path in (('medium',medium),('large',large),('local',local))}
    mc,mm=sources['medium'];lc,lm=sources['large'];cc,cm=sources['local']
    if mm['cases']!=lm['cases'] or cm['global_mace_collection']!=record(medium):
        raise InvalidArtifact('component inputs have incompatible physical/model preparation')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    rows={};references={'medium':[],'large':[]};checks=[]
    for case in CASES:
        r={'evidence':LABELS[case],'full':{},'local_medium':{}}
        for label,c in (('medium',mc),('large',lc)):
            pair=paired_components(*(c['rows'][f'{case}_{metal}_primary'] for metal in ('La','Ca')))
            references[label].append(pair['constant_atomic_reference_contrast_kcal_mol']);r['full'][label]=pair
            rotated=[c['rows'].get(f'{case}_{metal}_rotate') for metal in ('La','Ca')]
            if all(rotated):
                rotation=paired_components(*rotated)
                error=abs(rotation['R_short_kcal_mol']-pair['R_short_kcal_mol'])
                checks.append({'name':label+'_'+case+'_short_rotation','error_kcal_mol':error,'pass':error<=.01})
        for state in ('archived','global_H'):
            pair=paired_components(*(cc['rows'][f'{case}_{state}_{metal}'] for metal in ('La','Ca')))
            references['medium'].append(pair['constant_atomic_reference_contrast_kcal_mol']);r['local_medium'][state]=pair
        r['short_full_minus_core_global_H_kcal_mol']=r['full']['medium']['R_short_kcal_mol']-r['local_medium']['global_H']['R_short_kcal_mol']
        rows[case]=r
    for label,values in references.items():
        span=max(values)-min(values)
        checks.append({'name':label+'_constant_atomic_reference','span_kcal_mol':span,'pass':span<=1e-5})
    if not all(c['pass'] for c in checks):raise InvalidArtifact('component accounting or rigid invariance failed')
    comparisons={label:[] for label in ('medium','large','local_archived','local_global_H')}
    for higher,lower in CONTRASTS:
        for label in comparisons:
            if label in ('medium','large'):a,b=[rows[case]['full'][label]['R_short_kcal_mol'] for case in (higher,lower)]
            else:a,b=[rows[case]['local_medium'][label.removeprefix('local_')]['R_short_kcal_mol'] for case in (higher,lower)]
            comparisons[label].append({'higher_expected':higher,'lower_expected':lower,'difference_kcal_mol':a-b,'expected_order':a>b})
    # Verify the installed source belongs to the pinned backend inventory.
    software=read_json(verify(mm['software']))
    inventory=read_json(verify(software['backend_source_inventory']))
    implementation=next(ref for ref in inventory['files'] if ref['path'].endswith('/mace/modules/extensions.py'))
    verify(implementation)
    result={'status':'complete','protocol_id':PROTOCOL,'plan':record(plan),'implementation':record(__file__),
            'saved_sources':{k:record(v) for k,v in (('medium',medium),('large',large),('local',local))},
            'MACE_energy_implementation':implementation,'rows':rows,'contrasts':comparisons,'checks':checks,
            'primary_development_ordering_screen_pass':all(r['expected_order'] for r in comparisons['medium']),
            'new_MACE_calls':0,'new_GB_calls':0,'new_DFT_calls':0,'wall_seconds':time.monotonic()-start,
            'reference':None,'S_kcal_mol':None,'calibrated_class':None,'evidence_use':'consumed_method_development',
            'physical_energy_claim':'partial_learned_component_only; not_a_binding_energy',
            'charge_and_spin_dependence':'no_explicit_charge_or_spin_input_in_this_component; atomic_geometry/species_only'}
    write_new(out/'result.json',result)
    lines=['# Whole-protein short-range learned component','',
           'This is a separate empirical descriptor, not a complete energy or binding free energy.',
           'It uses the original trained `interaction_energy` readout, before charge restoration and field updates.',
           'No new MACE, GB or DFT evaluation. Baseline and failed total-energy descriptors unchanged.','',
           '| Case | Full medium Rshort | Full large Rshort | Local corrected-H Rshort | Full−local corrected-H |',
           '|---|---:|---:|---:|---:|']
    for case,r in rows.items():
        values=[r['full']['medium']['R_short_kcal_mol'],r['full']['large']['R_short_kcal_mol'],r['local_medium']['global_H']['R_short_kcal_mol'],r['short_full_minus_core_global_H_kcal_mol']]
        lines.append('| '+case+' | '+' | '.join(f'{v:.6f}' for v in values)+' |')
    lines.extend(['','All values kcal/mol; larger was predeclared as more La-like. No inherited zero/band.','',
                  '## Predeclared comparisons','','```json',json.dumps(comparisons,indent=2),'```','',
                  f"Primary development screen passes: {result['primary_development_ordering_screen_pass']}.",
                  'No threshold/weight was fitted. Cases are consumed development; alpha is one observation.',
                  'The component has finite-neighborhood sensitivity and no explicit total-charge/spin response.',
                  'Do not report total-energy forces as its gradients. Independent evaluation remains necessary.'])
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('medium','large','local','plan','output'):p.add_argument('--'+key,required=True)
    a=p.parse_args();r=assess(a.medium,a.large,a.local,a.plan,a.output)
    print(json.dumps({'status':r['status'],'contrasts':r['contrasts'],'primary_screen':r['primary_development_ordering_screen_pass']},indent=2))
