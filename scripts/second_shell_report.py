"""Compare completed second-shell endpoint records; no fitting or new calculations."""
import argparse
import json
from pathlib import Path
from affordable_common import HA_TO_KCAL,read_json,record,verify,write_new
from hydration_square import endpoint


def report(dft,mace,output):
    d=read_json(dft);m=read_json(mace);manifest=read_json(verify(d['manifest']))
    if d['status']!='complete' or m['status']!='complete':raise ValueError('completed actual endpoint collections required')
    rows=[]
    for dr in d['rows']:
        case=dr['case'];mr=next(r for r in m['cases'] if r['case']==case)
        source=next(s for s in manifest['states'] if s['case']==case)
        prep=read_json(verify(source['preparation']));core={z:endpoint(e['output'],e['receipt'],e['xyz'],e['input'])
                   for z,e in source['source']['endpoints'].items()}
        components={}
        for term in ('CPCM_dielectric','dispersion','gCP','SCF'):
            values=[pair[z]['components_hartree'].get(term) for pair in (dr['endpoints'],core) for z in ('Ca','La')]
            components[term]=((values[0]-values[1])-(values[2]-values[3]))*HA_TO_KCAL if all(v is not None for v in values) else None
        rows.append({'case':case,'group':dr['group'],'evidence_stratum':source['source']['evidence_stratum'],
          'before_atoms':dr['atoms_before'],'after_atoms':dr['atoms_after'],'added_formal_charge':dr['added_formal_charge'],
          'DFT_delta_R_kcal_mol':dr['delta_R_kcal_mol'],'MACE_delta_R_kcal_mol':mr['delta_R_kcal_mol'],
          'delta_direction_agrees':(dr['delta_R_kcal_mol']>0)==(mr['delta_R_kcal_mol']>0),
          'component_delta_R_kcal_mol':components,'DFT_core_R_hartree':dr['core_R_hartree'],'DFT_expanded_R_hartree':dr['expanded_R_hartree'],
          'MACE_core_R_kcal_mol':mr['core_R_kcal_mol'],'MACE_expanded_R_kcal_mol':mr['expanded_R_kcal_mol'],
          'removed_caps_count':len(prep['removed_caps']),'added_fragments':prep['added_fragments'],
          'excluded_water_contacts':len(prep['excluded_contacts'])})
    comparisons=[]
    for dc in d['comparisons']:
        x=next(r for r in rows if r['case']==dc['La_like']);y=next(r for r in rows if r['case']==dc['Ca_like'])
        before=x['MACE_core_R_kcal_mol']-y['MACE_core_R_kcal_mol'];after=x['MACE_expanded_R_kcal_mol']-y['MACE_expanded_R_kcal_mol']
        comparisons.append(dict(dc,DFT_direction_before=dc['before_kcal_mol']>0,DFT_direction_after=dc['after_kcal_mol']>0,
          MACE_before_kcal_mol=before,MACE_after_kcal_mol=after,MACE_direction_before=before>0,MACE_direction_after=after>0))
    calls=[r['accepted'] for r in m['rows']]
    timings={'DFT_endpoint_wall_sum_seconds':sum(e['wall_seconds'] for r in d['rows'] for e in r['endpoints'].values()),
       'DFT_endpoint_wall_seconds':{r['case']:{z:e['wall_seconds'] for z,e in r['endpoints'].items()} for r in d['rows']},
       'MACE_model_evaluation_seconds':sum(r['evaluation_seconds'] for r in calls),
       'MACE_worker_seconds':sum(r['wall_seconds'] for r in calls),
       'MACE_peak_cuda_allocated_bytes':max(r['peak_cuda_allocated_bytes'] for r in calls),
       'MACE_peak_host_RSS_KiB':max(r['peak_host_RSS_KiB'] for r in calls),
       'baseline_input_preparation_and_prior_calculations':'reused; excluded from new execution cost'}
    result={'DFT_collection':record(dft),'MACE_collection':record(mace),'rows':rows,'comparisons':comparisons,'timings':timings,
       'new_DFT_endpoints':10,'new_MACE_endpoint_evaluations':20,'reference':None,'calibrated_decision':None,
       'baseline_changed':False,'accuracy':'retrospective development directions only; alpha replicates share one biological group',
       'component_warning':'SCF includes CPCM; CPCM_dielectric is an overlapping diagnostic, not a second additive energy',
       'method_warning':'DFT uses CPCM; native unmasked OMOL has no solvent and total-charge conditioning. No equality of Hamiltonians claimed.'}
    out=Path(output);out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    lines=['# Actual fixed-coordinate second-shell results','',
       '| Case | Atoms core→expanded | Added charge | DFT ΔR | Native MACE ΔR |',
       '|---|---:|---:|---:|---:|']
    for r in rows:lines.append(f"| {r['case']} | {r['before_atoms']}→{r['after_atoms']} | {r['added_formal_charge']:+d} | {r['DFT_delta_R_kcal_mol']:+.4f} | {r['MACE_delta_R_kcal_mol']:+.4f} |")
    lines+=['','All energies in kcal/mol. ΔR=R_expanded−R_core;R=E_Ca−E_La.','',
      '| Expected La-like minus Ca-like | DFT core | DFT expanded | MACE core | MACE expanded |',
      '|---|---:|---:|---:|---:|']
    for r in comparisons:lines.append(f"| {r['La_like']}−{r['Ca_like']} | {r['before_kcal_mol']:+.4f} | {r['after_kcal_mol']:+.4f} | {r['MACE_before_kcal_mol']:+.4f} | {r['MACE_after_kcal_mol']:+.4f} |")
    lines+=['','Positive is the expected direction within each method. No absolute threshold inherited.',
            'PQQ is functional association; alpha/GGR is condition-qualified affinity ordering.',
            'The two alpha structures are one biological group. All cases were consumed previously.',
            'No production default or canonical25casePQQ result changed.','']
    (out/'TABLES.md').write_text('\n'.join(lines));return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('dft','mace','output'):p.add_argument('--'+name,required=True)
    r=report(**vars(p.parse_args()));print(json.dumps({'comparisons':r['comparisons'],'timings':r['timings']},indent=2))
