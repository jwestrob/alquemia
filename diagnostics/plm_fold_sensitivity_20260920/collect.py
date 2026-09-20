"""Collect the two declared PLM sensitivity runs; no execution or recalibration."""
import argparse
import csv
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from affordable_common import read_json, record, verify, write_new
from pqq_fast_prepare import coordinate_comparison
from compact_solvation_compare import mix_pair


def collect(submissions, prior, scheduler, output):
    ss=read_json(submissions);manifest=read_json(verify(ss['manifest']));old=read_json(prior)
    old_manifest=read_json(verify(old['manifest']));design=read_json(verify(old_manifest['design']))
    old_low=read_json(verify(old_manifest['low_manifest']))
    old_mace=read_json(verify(old['MACE_result']))
    states={r['JobID']:r for r in csv.DictReader(Path(scheduler).open(),delimiter='|')}
    results=[];cost=[]
    for job in ss['jobs']:
        pp=verify(job['plan']);plan=read_json(pp);request=read_json(verify(plan['request']));run=pp.parent
        ep=run/('ensemble_'+job['job_id']+'.json');ensemble=read_json(ep);group=ensemble['groups'][0]
        prep=read_json(run/'source_preparation/preparation.json');prepared={c['case_id']:c for c in prep['cases']}
        members=[]
        for member in group['members']:
            source=next(c for c in request['cases'] if c['case_id']==member['case_id'])
            prov=source['source_provenance'];selection=read_json(verify(prov['selection_review']))
            previous=next(m for m in selection['models'] if m['sample']==prov['source_sample'] and m['seed']==prov['source_seed'])
            if previous['source_cif']!=source['source_structure']:raise ValueError('prior source admission identity differs')
            extra=source['roles']['extra_acidic_ligand_homolog']
            distances=[t['distance_A'] for t in previous['typed_contacts'] if t['chain']==extra['chain'] and t['resnum']==extra['resnum'] and t['element']=='O']
            members.append({**member,'prior_single_fold_admission':previous,'extra_acidic_min_O_distance_A':min(distances) if distances else None})
        comparison={'status':'unavailable','reason':'sample0 current score/preparation unavailable'}
        sample=next(c for c in request['cases'] if c['source_provenance']['source_sample']==0)
        fresh=next(m for m in members if m['case_id']==sample['case_id'])
        if fresh['status']=='available':
            source=prepared[sample['case_id']];old_source=next(c for c in design['cases'] if c['case_id']==job['protein_id'])
            endpoints=source['representations']['context']['endpoints'];pair={z:next(p for p in old['points'] if p['case_id']==job['protein_id'] and p['point']=='origin' and p['metal']==z) for z in ('Ca','La')}
            sm=read_json(run/'prepared_score/scoring/manifest.json');low=read_json(verify(sm['solvent_manifests'][sample['case_id']]))
            coords={z:coordinate_comparison(verify(endpoints[z]['xyz']),verify(old_source['origins'][z]['xyz'])) for z in ('Ca','La')}
            same_state=all(endpoints[z]['charge']==old_source['origins'][z]['charge'] and endpoints[z]['multiplicity']==old_source['origins'][z]['multiplicity'] for z in ('Ca','La'))
            inputs=[]
            for t in low['tasks']:
                prev=next(t0 for t0 in old_low['tasks'] if t0['point_id']==job['protein_id']+'__origin__'+t['metal'] and t0['medium']==t['medium'])
                inputs.append({'metal':t['metal'],'medium':t['medium'],'input_identical':verify(t['input']).read_bytes()==verify(prev['input']).read_bytes(),'old_input':prev['input'],'new_input':t['input']})
            compatible=all(v['exact'] for v in coords.values()) and same_state and all(x['input_identical'] for x in inputs) and request['config']['model']==old_mace['model'] and old_low['orca']==low['orca']
            comparison={'status':'compatible' if compatible else 'preparation_or_method_mismatch','coordinates':coords,'charges_and_multiplicities_match':same_state,'input_comparisons':inputs,'model_matches':request['config']['model']==old_mace['model'],'orca_matches':old_low['orca']==low['orca'],'delta_R_model_kcal_mol':None}
            if compatible and all(p['status']=='complete' for p in pair.values()):
                score=mix_pair({z:pair[z]['MACE_eV'] for z in pair},{z:pair[z]['GFN2_vacuum_hartree'] for z in pair},{z:pair[z]['GFN2_ALPB_hartree'] for z in pair})
                comparison.update(previous_components=score,delta_R_model_kcal_mol=fresh['composite_R_model_kcal_mol']-score['composite_R_model_kcal_mol'],endpoint_differences={z:{'MACE_eV':fresh['native_MACE'][z]['energy_eV']-pair[z]['MACE_eV'],'GFN2_vacuum_hartree':fresh['low_level_endpoints'][z]['vacuum_hartree']-pair[z]['GFN2_vacuum_hartree'],'GFN2_ALPB_hartree':fresh['low_level_endpoints'][z]['alpb_hartree']-pair[z]['GFN2_ALPB_hartree']} for z in pair})
        state=states[job['job_id']];seconds=int(state['ElapsedRaw']);cpus=int(state['AllocCPUS']);execution=run/'execution.json'
        cost.append({'job_id':job['job_id'],'state':state['State'],'elapsed_seconds':seconds,'allocated_CPUs':cpus,'allocated_GPUs':1,'allocated_core_seconds':seconds*cpus,'allocated_GPU_seconds':seconds,'source_to_score_seconds':read_json(execution)['source_to_score_seconds'] if execution.exists() else None,'execution':record(execution) if execution.exists() else None,'available_sources':group['available_members'],'prepared_sources':prep['supported']})
        results.append({'protein_id':job['protein_id'],'ensemble':record(ep),'source_request':plan['request'],'members':members,'required_members':3,'available_members':group['available_members'],'status':group['status'],'median_R_model_kcal_mol':group['median_R_model_kcal_mol'],'spread_model_kcal_mol':group['components']['composite_R_model_kcal_mol']['spread'],'developmental_band_transfer':group['developmental_frozen_band_transfer'],'prior_sample0_comparison':comparison})
    result={'scope':manifest,'submissions':record(submissions),'prior_torsion_result':record(prior),'scheduler':record(scheduler),'implementation':record(__file__),'rows':results,'costs':cost,'total_allocated_core_seconds':sum(c['allocated_core_seconds'] for c in cost),'total_allocated_GPU_seconds':sum(c['allocated_GPU_seconds'] for c in cost),'biological_labels':None,'accuracy':None,'production_changed':False,'threshold_refitted':False}
    write_new(output,result)
    return {r['protein_id']:{k:r[k] for k in ('status','available_members','median_R_model_kcal_mol','spread_model_kcal_mol','developmental_band_transfer')} for r in results}

if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('submissions','prior','scheduler','output'):p.add_argument('--'+n,required=True)
    print(json.dumps(collect(**vars(p.parse_args())),indent=2))
