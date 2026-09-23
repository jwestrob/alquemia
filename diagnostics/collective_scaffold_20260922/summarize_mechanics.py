"""Read only completed scaffold receipts; compute no molecular energies/forces."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,record,write_new


def summarize(collection,output):
    collected=read_json(collection);manifest=read_json(verify(collected['manifest']))
    if collected['terminal_searches']!=24:raise ValueError('full terminal 24-target collection required')
    rows=[];responses={}
    for entry in collected['rows']:
        row={k:entry[k] for k in ('case_id','target','status','candidate_admitted')}
        row['receipt']=entry.get('receipt')
        if entry['status']!='complete':rows.append(row);continue
        receipt=read_json(verify(entry['receipt']));case=read_json(verify(receipt['parent_case']))
        opt=read_json(verify(receipt['optimization']));atoms=read_json(verify(case['parent']['atoms']))
        source=np.asarray(read_json(verify(case['parent']['positions_A'])))
        initial=np.asarray(read_json(verify(receipt['initial_parent_positions_A'])))
        final=np.asarray(read_json(verify(receipt['final_parent_positions_A'])))
        mobile=read_json(verify(case['mobile_set']))['mobile_parent_indices']
        mapping=read_json(verify(case['context_parent_mapping']))
        local=set(mapping['local_parent_indices']);donors={a['parent_index'] for a in read_json(verify(case['donors']))['atoms']}
        heavy=[i for i in mobile if atoms[i]['element']!='H'];hydrogen=[i for i in mobile if atoms[i]['element']=='H']
        free_heavy=[i for i in heavy if i not in donors];response=final-initial
        displacement=np.linalg.norm(final-source,axis=1)
        boundary=[i for i in heavy if displacement[i]>=.8-1e-5]
        limits=[]
        for i in boundary:
            limits.append({'source_id':atoms[i]['source_id'],'resname':atoms[i]['resname'],
                           'displacement_A':float(displacement[i]),'donor_target':i in donors,'in_context':i in local})
        reasons=Counter(r['reason'] for r in opt['rejected_trials']);guards=Counter()
        for trial in opt['rejected_trials']:
            g=trial.get('geometry')
            if not g:continue
            if g['maximum_heavy_displacement_A']>.8+1e-7:guards['heavy_displacement']+=1
            if g['source_bond_max_error_A']>1e-6:guards['source_bond']+=1
            if g['fixed_max_error_A']>1e-12:guards['fixed_atoms']+=1
            if g['chirality_failures']:guards['true_stereocenter']+=1
            if g['peptide_basin_failures']:guards['peptide_cis_trans']+=1
            if g['new_severe_overlaps']:guards['new_severe_overlap']+=1
        row.update(optimizer_status=opt['status'],stationary=opt['stationary'],accepted_iterations=opt['accepted_iterations'],
                   parent_work_kcal_mol=opt['work_kcal_mol'],parent_energy_force_calls=opt['parent_energy_force_calls'],
                   wall_seconds=receipt['wall_seconds'],evaluation_seconds=receipt['parent_evaluation_seconds'],
                   first_evaluation_seconds=opt['first_declared_evaluation_seconds'],
                   projected_max_atom_kcal_mol_A=opt['final_projected_max_atom_kcal_mol_A'],
                   projected_rms_component_kcal_mol_A=opt['final_projected_rms_component_kcal_mol_A'],
                   heavy_source_max_A=float(displacement[heavy].max()),
                   heavy_response_RMS_A=float(np.sqrt(np.mean(np.sum(response[heavy]**2,axis=1)))),
                   H_response_RMS_A=float(np.sqrt(np.mean(np.sum(response[hydrogen]**2,axis=1)))),
                   boundary_atoms=limits,rejection_reasons=dict(reasons),rejection_guards=dict(guards),
                   actual_runtime=receipt['runtime'])
        responses[(entry['case_id'],entry['target'])]=(response,free_heavy,final)
        rows.append(row)
    comparisons=[]
    for cid in dict.fromkeys(r['case_id'] for r in rows):
        origin=responses.get((cid,'origin'))
        for target in ('Ca_adaptive','La_adaptive'):
            other=responses.get((cid,target));r={'case_id':cid,'target':target,'status':'unavailable'}
            if origin is not None and other is not None:
                a,ids,ofinal=origin;b,other_ids,tfinal=other
                if ids!=other_ids:raise ValueError('source response physical measure changed')
                av=a[ids].ravel();bv=b[ids].ravel();norm=np.linalg.norm(av)*np.linalg.norm(bv)
                r.update(status='available',heavy_response_cosine=None if norm==0 else float(np.dot(av,bv)/norm),
                         heavy_response_difference_RMS_A=float(np.sqrt(np.mean(np.sum((a[ids]-b[ids])**2,axis=1)))),
                         final_target_difference_RMS_A=float(np.sqrt(np.mean(np.sum((tfinal[ids]-ofinal[ids])**2,axis=1)))))
            comparisons.append(r)
    result={'collection':record(collection),'analysis_implementation':record(__file__),'rows':rows,
            'matched_origin_response_comparisons':comparisons,'new_molecular_calls':0,
            'stationary_count':sum(r.get('stationary',False) for r in rows),
            'admitted_count':sum(r['candidate_admitted'] for r in rows),
            'forcefield_energy_is_scoring_component':False}
    write_new(output,result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--collection',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();result=summarize(a.collection,a.output)
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','matched_origin_response_comparisons')},indent=2))
