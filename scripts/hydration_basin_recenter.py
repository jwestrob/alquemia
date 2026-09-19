"""Recenter the same bounded response using actual new native DFT gradients."""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact,read_json,verify,record,write_new,xyz
from hydration_water_motion import checked_gradient
from hydration_square import endpoint
from mace_hybrid import check_atoms
from hydration_basin import STAGE,PROTOCOL,SETTINGS


def prepare(validation,agreement,output):
    import mace_omol as omol
    from hydration_basin_validation import collect_mace
    v=read_json(validation);p=read_json(verify(v['preparation']))
    if v['status']!='complete':raise InvalidArtifact('actual complete native validation required')
    mc=read_json(verify(v['MACE_collection']));dc=read_json(verify(v['DFT_collection']))
    if mc!=collect_mace(verify(mc['manifest'])) or dc['status']!='complete':raise InvalidArtifact('changed/partial scientific results')
    source_collection=read_json(verify(p['source_collection']));parent=read_json(verify(source_collection['manifest']))
    parent_p=read_json(verify(parent['preparation']));generation=parent_p.get('generation',0)+1
    if generation>2:raise InvalidArtifact('two declared recentering rounds already used')
    dm=read_json(verify(dc['manifest']));dt={t['task_id']:t for t in dm['tasks']}
    dr={r['task_id']:r for r in dc['rows']};mr={r['task_id']:r for r in mc['rows']};centers={};excluded={}
    for name,s in p['selection'].items():
        checks=[r for r in v['direction_checks'] if r['center_id']==name]
        if p.get('proposals_only'):
            # This later step inherits the original local checks, explicitly;
            # it must still pass its own actual native energy check below.
            old=read_json(verify(parent_p['previous_validation']))
            old_checks=[r for r in old['direction_checks'] if r['center_id']==name and r['kind']=='response']
            local_pass=len(old_checks)==1 and old_checks[0]['pass']
        else:
            response=[r for r in checks if r['kind']=='response']
            local_pass=len(response)==1 and response[0]['pass']
        row=next(r for r in v['rows'] if r['center_id']==name and r['kind'] in ('minimum','proposal'))
        if not local_pass or not row['energy_pass']:
            excluded[name]={'status':'unsupported_native_response','local_checks_pass':local_pass,'proposal_energy_pass':row['energy_pass']};continue
        if row['kind']=='minimum' and row['stationary_within_tolerance']:
            excluded[name]={'status':'native_stationary_minimum_already_available','validation':record(validation),'task_id':row['task_id']};continue
        t=dt[row['task_id']];d=dr[row['task_id']];a=mr[row['task_id']];c=s['center']
        actual=endpoint(d['result']['output'],d['result']['receipt'],t['xyz'],t['input'])
        if actual!=d['result']:raise InvalidArtifact('new anchor differs from native endpoint')
        atoms=xyz(verify(t['xyz']));checked_gradient(d['gradient'],actual,atoms)
        forces=np.load(verify(a['forces']),allow_pickle=False)
        if forces.shape!=(len(atoms),3):raise InvalidArtifact('new force ordering differs')
        origin=c.get('origin_xyz',c['xyz']);original=xyz(verify(origin))
        mobile={i for w in c['groups'] if w['role']=='variable' for i in w['indices']}
        if any(atoms[i]!=original[i] for i in range(len(atoms)) if i not in mobile):raise InvalidArtifact('scaffold/metal changed across response steps')
        centers[name]={**c,'xyz':t['xyz'],'DFT_result':actual,'DFT_gradient':d['gradient'],
            'MACE_result':a['result'],'MACE_forces':a['forces'],'MACE_energy_eV':a['energy_eV'],
            'MACE_manifest':mc['manifest'],'origin_xyz':origin,'previous_xyz':c['xyz'],
            'previous_validation':record(validation),'generation':generation,'entropy_curvature_qualified':len(checks)==2 and all(r['pass'] for r in checks),
            'total_max_atom_displacement_A':float(np.max(np.linalg.norm(np.array([a[1:] for a in atoms])-np.array([a[1:] for a in original]),axis=1))),
            'metal_oxygen_distances_A':{str(w['source']['resnum']):float(np.linalg.norm(np.array(atoms[w['oxygen_index']][1:])-np.array(atoms[0][1:]))) for w in c['groups'] if w['role']=='variable'}}
    if not centers:
        out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
        write_new(out/'no_new_tasks.json',{'status':'no_supported_recenter_needed_or_available','excluded':excluded,'validation':record(validation)})
        return {'status':'no_new_tasks','excluded':excluded}
    _,out,m=omol.common(verify(parent['inventory']),verify(parent['software']),agreement,output,STAGE)
    for name in ('hydration_basin.py','hydration_basin_coordinates.py','hydration_basin_recenter.py','hydration_basin_validation.py','hydration_water_motion.py','hydration_square.py'):
        path=out/'implementation'/name;shutil.copyfile(Path(__file__).with_name(name),path);m['implementation'][name]=record(path)
    pp=out/'preparation.json';write_new(pp,{'protocol_id':PROTOCOL,'generation':generation,'centers':centers,
        'settings':SETTINGS,'source_collection':parent_p['source_collection'],'previous_validation':record(validation),
        'excluded':excluded,'agreement':record(agreement)})
    tasks=[]
    for name,c in sorted(centers.items()):
        tasks.append({'task_id':name,'case_id':c['case'],'metal':c['metal'],'metal_index':0,'kind':'core',
            'variant':'primary','charge':c['charge'],'spin_multiplicity':1,'energy_component':omol.COMPONENT,
            'energy_only':False,'xyz':c['xyz'],'state':check_atoms(xyz(verify(c['xyz'])),c['charge']),'preparation':record(pp)})
    m.update(tasks=tasks,preparation=record(pp),settings=SETTINGS,protocol_id=PROTOCOL,
             generation=generation,evidence_use='consumed_bounded_water_recenter_development')
    return omol.seal(out,m)


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('validation','agreement','output'):p.add_argument('--'+name,required=True)
    print(json.dumps(prepare(**vars(p.parse_args())),indent=2))
