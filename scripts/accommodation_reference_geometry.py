"""Existing-distance warning and terminal-mode support on prepared PQQ folds.

No labels/energies are read to select structures; no molecular calls are made.
"""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import os
from pathlib import Path
import time
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz

CUTOFF_A=2.2
ROLES={'anchor_glutamate':'GLU','extra_acidic_ligand_homolog':'ASP'}


def scan(preparation,warning_reference,output):
    d=read_json(preparation);rows=[]
    for r in d['cases']:
        source=r['source'];row={'case_id':r['case_id'],'root_case_id':source['root_case_id'],
            'source_conditioning_metal':source['source_conditioning_metal'],
            'source_native_metal':source['raw_source_metal'],
            'source_structure':source.get('original_source_structure',source['source_structure']),
            'preparation_status':r['status'],'preparation_reason':r.get('reason'),
            'screen_status':'unavailable','roles':{},'flagged':None,'mapping_status':'not_checked'}
        if r['status']=='prepared':
            p=read_json(verify(r['core']['parent']));atoms=xyz(verify(r['core']['endpoints']['La']['xyz']))
            ca=xyz(verify(r['core']['endpoints']['Ca']['xyz']))
            if atoms[1:]!=ca[1:] or atoms[0][1:]!=ca[0][1:]:raise InvalidArtifact('source endpoint geometry differs')
            center=np.array(atoms[0][1:]);index=1
            if not np.allclose(center,source['raw_source_metal']['xyz_A'],atol=1e-9,rtol=0):raise InvalidArtifact('source selected metal coordinates differ')
            for f in p['qm_fragments']:
                items=f['atom_records'];role=f.get('role')
                if role in ROLES:
                    rn=f['id'].split(':')[1][:3]
                    if rn!=ROLES[role]:row['roles'][role]={'status':'not_applicable_residue_is_'+rn,'source_id':f['id'],'minimum_distance_A':None,'warning':False}
                    else:
                        values=[]
                        for j,a in enumerate(items):
                            if a['element']!='O':continue
                            coord=atoms[index+j]
                            if coord[0]!='O' or not np.allclose(coord[1:],a['xyz_A'],atol=1e-9,rtol=0):raise InvalidArtifact('source oxygen map differs')
                            values.append({'atom':a['name'],'qm_index':index+j,'distance_A':float(np.linalg.norm(np.array(coord[1:])-center))})
                        if len(values)!=2:raise InvalidArtifact('incomplete target carboxylate')
                        minimum=min(v['distance_A'] for v in values)
                        row['roles'][role]={'status':'available','source_id':f['id'],'oxygen_distances':values,
                            'minimum_distance_A':minimum,'warning':minimum<CUTOFF_A}
                index+=len(items)
            for role,expected in ROLES.items():
                if role in row['roles']:continue
                selector=p['fixed_core']['requested_roles'][role]
                rn=selector['resname'] if isinstance(selector,dict) else selector.split(':')[1][:3]
                if rn==expected:raise InvalidArtifact('required acidic role missing from prepared core')
                row['roles'][role]={'status':'not_applicable_residue_is_'+rn,'source_id':selector,'minimum_distance_A':None,'warning':False}
            row.update(screen_status='available',flagged=any(v['warning'] for v in row['roles'].values()),core=r['core']['parent'])
            row['mapping_status']='pending' if row['flagged'] else 'not_selected_by_geometry_rule'
        rows.append(row)
    by={}
    for metal in ('Ca','La'):
        subset=[r for r in rows if r['source_conditioning_metal']==metal]
        good=[r for r in subset if r['screen_status']=='available'];flagged=[r for r in good if r['flagged']]
        by[metal]={'source_denominator':len(subset),'prepared':len(good),'unsupported':len(subset)-len(good),'flagged':len(flagged),
            'by_role':{role:sum(r['roles'][role]['warning'] for r in good) for role in ROLES},
            'flagged_protein_groups':len({r['root_case_id'] for r in flagged})}
    result={'schema':'reference_donor_geometry_warning_v1','implementation':record(__file__),'source_preparation':record(preparation),
        'warning_reference':record(warning_reference),'threshold_A':CUTOFF_A,'comparison':'strictly_less_than',
        'interpretation':'existing La-oriented screening flag; Ca geometry is not thereby defective',
        'rows':rows,'denominator':len(rows),'prepared':sum(r['screen_status']=='available' for r in rows),
        'flagged':sum(r['flagged'] is True for r in rows),'by_source_condition':by,'energies_read':False,'labels_used':False,'new_molecular_calls':0}
    write_new(output,result);return result


def mapping_case(case,topology,destination):
    from second_shell_context import parent_state
    from coordination_preparation_context import geometry
    from accommodation_torsion_profiles import role_mode
    from mace_site_kinematics import Kinematics
    tick=time.monotonic();result={'case_id':case['case_id'],'status':'unsupported','roles':{},'reason':None}
    try:
        state=parent_state(case['core'],topology,require_endpoint_receipts=False)
        rep=case['representations']['context'];prep=read_json(verify(rep['preparation']))
        la,ca=[xyz(verify(rep['endpoints'][z]['xyz'])) for z in ('La','Ca')]
        if la[1:]!=ca[1:] or la[0][1:]!=ca[0][1:]:raise InvalidArtifact('context source pair coordinates differ')
        if rep['endpoints']['La']['charge']-rep['endpoints']['Ca']['charge']!=1:raise InvalidArtifact('paired charge differs')
        gp=geometry(state,prep,la,state['original']['La']);kin=Kinematics(gp['context'])
        origin=kin.evaluate(np.zeros(len(kin.modes)))[1]
        if not np.allclose(origin,[a[1:] for a in la],atol=1e-9,rtol=0):raise InvalidArtifact('physical origin differs')
        for role,rn in ROLES.items():
            f=next((f for f in state['parent']['qm_fragments'] if f.get('role')==role),None)
            if f is None or f['id'].split(':')[1][:3]!=rn:
                result['roles'][role]={'status':'not_applicable','mode':None};continue
            mode=role_mode(state['parent'],role);match=[m for m in kin.modes if m['id']==mode]
            if len(match)!=1:raise InvalidArtifact('physical terminal mode unavailable')
            moving=match[0]['moving_indices'];meta=kin.data['source_atom_metadata']
            names={meta[i]['atom'] for i in moving if meta[i] is not None};required={'OE1','OE2'} if rn=='GLU' else {'OD1','OD2'}
            if not required<=names:raise InvalidArtifact('incomplete physical carboxylate mode')
            result['roles'][role]={'status':'supported','mode':mode,'moving_atom_names':sorted(names)}
        path=Path(destination)/(case['case_id']+'.json');write_new(path,gp)
        result.update(status='supported',mapping=record(path),checks=gp['checks'],
            paired_coordinate_measure='identical Cartesian origins and physical modes; metal element/charge differ')
    except Exception as exc:result.update(reason=str(exc),exception_type=type(exc).__name__)
    result['wall_seconds']=time.monotonic()-tick
    return result


def mappings(inventory,output,workers):
    if workers>1 and not os.environ.get('SLURM_JOB_ID'):raise InvalidArtifact('parallel mapping checks require allocation')
    start=time.monotonic();inv=read_json(inventory);data=read_json(verify(inv['source_preparation']));out=Path(output);out.mkdir(parents=True,exist_ok=False)
    selected={r['case_id'] for r in inv['rows'] if r['flagged'] is True};cases=[r for r in data['cases'] if r['case_id'] in selected]
    checked=[]
    if workers==1:
        for case in cases:checked.append(mapping_case(case,data['config']['topology'],out/'maps'))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            jobs=[pool.submit(mapping_case,c,data['config']['topology'],out/'maps') for c in cases]
            for job in as_completed(jobs):checked.append(job.result())
    by={r['case_id']:r for r in checked};rows=[]
    for r in inv['rows']:
        rows.append({**r,'mapping_check':by.get(r['case_id']),'mapping_status':by[r['case_id']]['status'] if r['case_id'] in by else r['mapping_status']})
    result={k:v for k,v in inv.items() if k not in ('implementation','rows')}
    result.update(inventory=record(inventory),mapping_implementation=record(__file__),rows=rows,
        flagged_mapping_denominator=len(selected),mapping_status_counts=dict(Counter(r['status'] for r in checked)),
        mapping_wall_seconds=time.monotonic()-start,workers=workers,slurm_job_id=os.environ.get('SLURM_JOB_ID'),new_molecular_calls=0)
    write_new(out/'result.json',result)
    return {'mapping_status_counts':result['mapping_status_counts'],'mapping_wall_seconds':result['mapping_wall_seconds']}


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('scan')
    for k in ('preparation','warning_reference','output'):q.add_argument('--'+k.replace('_','-'),required=True)
    q=s.add_parser('mappings');q.add_argument('--inventory',required=True);q.add_argument('--output',required=True);q.add_argument('--workers',type=int,required=True)
    a=vars(p.parse_args());op=a.pop('op');r=globals()[op](**a);print(json.dumps({k:v for k,v in r.items() if k!='rows'},indent=2))
if __name__=='__main__':main()
