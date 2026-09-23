"""Materialize bounded molecular requests; no molecular calls or submission."""
import argparse
from pathlib import Path
import json
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from compact_solvation_compare import native_endpoint
from mace_site_kinematics import Kinematics
from adaptive_force_diagnostic import project,preview
import adaptive_angular_proposals as angular
import union_adaptive as union
import union_triple_transfer as transfer


def prepare(audit,output):
    a=read_json(audit);rows=[];pending=[];ready=[]
    if a['protocol_id']!=transfer.PROTOCOL:raise InvalidArtifact('declared threefold transfer required')
    for r in a['rows']:
        if r['pool_reuse']:continue
        state=[];projections=[];endpoints=[]
        for z in ('Ca','La'):
            ep=r['source']['representations']['context']['endpoints'][z]
            kin=Kinematics(read_json(verify(r['mapping'][z]))['context'])
            t={'task_id':r['pair_id']+'__'+z,'case_id':r['pair_id'],'source_case_id':r['case_id'],
                'selection_id':r['selection_id'],'metal':z,'mapping':r['mapping'][z],'xyz':ep['xyz'],
                'charge':ep['charge'],'multiplicity':ep['multiplicity'],'mode_count':len(kin.modes),
                'origin_native_receipt':r['origin_reuse']['native_endpoints'].get(z,{}).get('receipt'),
                'origin_solvent':r['origin_reuse']['solvent_endpoints'].get(z,{}),
                'optimizer_starts':1,'settings':a['settings']}
            if t['origin_native_receipt']:
                ne=native_endpoint(read_json(verify(t['origin_native_receipt'])),a['config']['model'])
                native,force=union.origin(ne,a['config']['model'])
                _,v,raw,normed,_=project(kin.data,np.zeros(len(kin.modes)),force)
                projections.append((kin,v,normed));t['origin_status']='actual_energy_and_forces_reusable'
            else:t['origin_status']='await_finite_origin_tasks'
            endpoints.append(t)
        if len(projections)==2:
            choice=preview(projections[0][0].modes,projections[0][1],projections[0][2],projections[1][2])
            ids=[x['id'] for x in choice['selected']]
            if len(ids)!=4:raise InvalidArtifact('four independent physical modes unavailable')
            for t in endpoints:
                kin=Kinematics(read_json(verify(t['mapping']))['context'])
                t.update(active_indices=[[x['id'] for x in kin.modes].index(mid) for mid in ids],active_mode_ids=ids,selector=choice)
                angular.final_geometry(kin,t,np.zeros(4),[x[0] for x in xyz(verify(t['xyz']))])
            ready.append(r['pair_id'])
        else:
            if projections:raise InvalidArtifact('unexpected partially available force pair')
            pending.append(r['pair_id'])
        rows.extend(endpoints)
    if (len(rows),len(ready),len(pending))!=(110,36,19):raise InvalidArtifact('fixed110 search/36 ready/19 origin-pending scope differs')
    cross=[];low=[]
    for row in a['rows']:
        if row['pool_reuse']:continue
        for candidate in ('adaptive_Ca','adaptive_La'):
            for z in ('Ca','La'):
                t={'task_id':row['pair_id']+'__'+z+'__at_'+candidate,'case_id':row['pair_id'],'source_case_id':row['case_id'],
                   'selection_id':row['selection_id'],'candidate':candidate,'metal':z,'coordinate_rule':'exact final admitted physical map; same candidate for both metals',
                   'source_mapping':row['mapping'][z],'charge':row['source']['representations']['context']['endpoints'][z]['charge'],
                   'multiplicity':row['source']['representations']['context']['endpoints'][z]['multiplicity']}
                if candidate!='adaptive_'+z:cross.append(t)
                for medium in ('vacuum','alpb'):low.append({**t,'task_id':t['task_id']+'__'+medium,'medium':medium,'solver':'native','maxiter':500,'mpi_ranks':1})
    result={'protocol_id':transfer.PROTOCOL,'audit':record(audit),'reference':a['reference'],'agreement':a['agreement'],
        'counts':a['counts'],'search_tasks':rows,'cross_MACE_requests':cross,'candidate_GFN2_requests':low,
        'fully_preflighted_force_pairs':ready,'force_pairs_waiting_declared_origins':pending,
        'candidate_coordinates':'not yet generated; exact native manifests materialize only after admitted searches',
        'submission_authorized':False,'new_molecular_calls':0,'production_changed':False,'implementation':record(__file__)}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('search_tasks','cross_MACE_requests','candidate_GFN2_requests','fully_preflighted_force_pairs','force_pairs_waiting_declared_origins')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('audit','output'):p.add_argument('--'+k,required=True)
    print(json.dumps(prepare(**vars(p.parse_args())),indent=2))
