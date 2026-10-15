#!/usr/bin/env python3
"""Bind real saved inputs; prepare accommodated cavity unions only inside Slurm."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback

PROTOCOL='native_ALPB_capped_union_saved_accommodation_preparation_v1'
INVENTORY_SHA='9b38d4b184ddb69079fbd955bc517ec4abb642d1d71b011b734cd7efb01d7658'


def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def check(p):
    if pin(p['path'])!=p:raise ValueError('Changed input: '+p['path'])
    return Path(p['path'])


def read(p):
    return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())


def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')


def load(name,p):
    spec=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m


def freeze(root,output):
    root=Path(root).resolve();out=Path(output).resolve()
    ip=pin(root/'diagnostics/native_cavity_shielding_20261014/accommodation_feasibility/INVENTORY.json')
    if ip['sha256']!=INVENTORY_SHA:raise ValueError('Wrong saved-cell inventory')
    inv=read(ip);origin_pin=inv['sources']['origin_preparation'];origin=read(origin_pin)
    old_manifest=read(origin['manifest']);km_pin=pin(root/'workspaces/exaa_metal_midpoint_20261014/preparation_v1/manifest.json');km=read(km_pin)
    base=old_manifest['implementation']['PREPARE.py'];check(base)
    requests=[]
    for sample in (0,3):
        sid=f'q9z4j7-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-{sample}'
        for candidate in ('adaptive_Ca','adaptive_La'):
            motion,=[x for x in inv['saved_motion_metadata'] if x['source_case_id']==sid and x['candidate']==candidate]
            for metal in ('Ca','La'):
                cell,=[x for x in inv['rows'] if (x['source_case_id'],x['candidate'],x['metal'])==(sid,candidate,metal)]
                old,=[x for x in origin['rows'] if x['source_case_id']==sid and x['metal']==metal]
                if cell['status']!='available' or old['status']!='prepared':raise ValueError('Missing declared cell')
                requests.append({'id':f'sample{sample}_{candidate}_{metal}','sample':sample,'donor_candidate':candidate,
                                 'metal':metal,'cell':cell,'origin_row':old,'motion':motion})
    out.mkdir(parents=True,exist_ok=False);(out/'implementation').mkdir();(out/'inputs').mkdir();impl={};snapshots={}
    for name in ('PREPARE.py','RUN.sbatch'):
        d=out/'implementation'/name;shutil.copyfile(Path(__file__).with_name(name),d);impl[name]=pin(d)
    # Actual original coordinate/population files are copied, never reconstructed at freeze.
    for req in requests:
        for p in (req['cell']['xyz'],req['cell']['vacuum_charges'],req['cell']['native_result'],
                  req['cell']['state_signature'],req['cell']['source_preparation'],req['motion']['proposal'],
                  req['motion']['proposal_manifest'],req['cell']['maps']['Ca'],req['cell']['maps']['La']):
            original=check(p)
            if p['sha256'] not in snapshots:
                d=out/'inputs'/(p['sha256']+original.suffix);shutil.copyfile(original,d)
                snapshots[p['sha256']]={'original':p,'snapshot':pin(d)}
    m={'protocol':PROTOCOL,'inventory':ip,'origin_preparation':origin_pin,'origin_native_collection':inv['sources']['origin_native_collection'],
       'plan':pin(root/'diagnostics/native_cavity_accommodation_20261014/PLAN.md'),
       'base_preparation':base,'kinematics_provenance':km_pin,'source_snapshot':km['source_snapshot'],'source_modules':km['source_modules'],
       'python':km['python'],'implementation':impl,'input_snapshots':snapshots,'requests':requests,
       'settings':{'local_atoms':192,'cap_count':14,'kinematic_replay_A':1e-8,'source_coordinate_tolerance_A':1e-12,
                   'new_overlap_heavy_A':1.0,'new_overlap_H_A':0.55,
                   'contact_policy':'unchanged severe-overlap policy for newly introduced physical core/environment contacts versus origin; synthetic caps excluded',
                   'environment':'unchanged exact origin union positions and zero charges','charge_renormalization':False},
       'declared_nonorigin_cells':8,'output':str(out/'execution'),'new_energy_calls':0,'new_SCC_calls':0,'scientific_preparation':'unrun'}
    write(out/'manifest.json',m);print(json.dumps(pin(out/'manifest.json')))


def prepare(req,out,m,base,np,metal_geometry,Kinematics):
    c=req['cell'];old=req['origin_row'];motion=req['motion'];metal=req['metal'];candidate=req['donor_candidate']
    cp=read(c['source_preparation']);state=read(c['state_signature']);mapping=read(c['maps']['Ca'])['context']
    if read(c['maps']['Ca'])!=read(c['maps']['La']):raise ValueError('Paired maps changed')
    if (c['atom_count'],c['multiplicity'],c['valence_electrons'])!=(192,1,588) or c['charge']!=state['endpoint_charges'][metal]:raise ValueError('State differs')
    if state['water_inventory'] or state['endpoint_charges']!={'Ca':-3,'La':-2}:raise ValueError('Original dry charge state differs')
    if c['source_preparation']!=old['source_preparation'] or c['state_signature']!=old['state_signature'] or cp['source']!=old['protonated_source']:
        raise ValueError('Origin/source provenance differs')
    native=read(c['native_result']);ni=native['input'];proposal=read(motion['proposal']);pm=read(motion['proposal_manifest'])
    if native['status']!='available' or ni['xyz']!=c['xyz'] or ni['state_signature']!=c['state_signature'] or not c['component_v2']['accepted']:
        raise ValueError('Missing candidate solvent qualification')
    if not all(x['passed'] and x['informative'] for x in c['native_diagonal_checks']) or not all(x['informative'] for x in c['vacuum_precision']):
        raise ValueError('Candidate native numerical gate failed')
    for p in (c['xyz'],c['vacuum_charges'],c['MACE'],c['MACE_request'],c['vacuum_SCC']['point'],c['vacuum_SCC']['execution']):check(p)
    if c['vacuum_charges']['sha256']!=native['vacuum']['charges']['sha256'] or ni['endpoints']['vacuum']['executed_xyz']['sha256']!=c['xyz']['sha256']:
        raise ValueError('Candidate charge/coordinate mismatch')
    if proposal['status']!='proposal_available' or not proposal['proposal']['geometry_checks']['pass']:raise ValueError('Original physical proposal unavailable')
    task,=[copy.deepcopy(t) for t in pm['tasks'] if t['metal']==candidate.rsplit('_',1)[1]]
    kin=Kinematics(mapping);q=np.asarray(proposal['proposal']['full_q']);moving={i for k in task['active_indices'] for i in kin.modes[k]['moving_indices']}
    local_by_physical={x[2]:x[0] for x in kin.links if x[1]=='source'}
    if any(i not in local_by_physical for i in moving):raise ValueError('Moving source atom is omitted from local block')
    if any(x[1]=='cap' and (x[2] in moving or x[3] in moving) for x in kin.links):raise ValueError('Unsupported moving cap/source boundary')
    midx=[next(i for i,v in enumerate(kin.modes) if v['id']=='metal_'+axis) for axis in 'xyz']
    if not np.array_equal(q[midx],np.zeros(3)):raise ValueError('Saved metal is not fixed')
    task['active_indices']=midx+task['active_indices']
    if not np.array_equal(metal_geometry.full_q(task,q[task['active_indices']]),q):raise ValueError('Saved q outside original selected modes')
    a=base.xyz_atoms(check(c['xyz']));xyz=np.asarray([v['xyz'] for v in a]);physical,replayed,_,_=kin.evaluate(q)
    residual=float(np.max(np.abs(replayed-xyz)))
    if residual>m['settings']['kinematic_replay_A']:raise ValueError('Saved candidate does not replay')
    gate=metal_geometry.final_geometry(kin,task,q[task['active_indices']],[v['element'] for v in a])
    if not gate['physical_feasible'] or not gate['mapping_checks']['pass'] or not gate['all_other_physical_atoms_fixed']:
        raise ValueError('Original physical/cap/bond/overlap admission failed')
    uv=old['variants']['union'];um=read(uv['mapping']);ua=base.xyz_atoms(check(uv['xyz']));uq=check(uv['charges']).read_text().split()
    if len(ua)!=uv['atom_count'] or len(uq)!=len(ua) or len(um['centers'])!=len(ua):raise ValueError('Origin union inventory differs')
    if [x['local_index'] for x in um['centers'][:192]]!=list(range(192)):raise ValueError('Origin union core order differs')
    origin=np.asarray([v['xyz'] for v in ua]);symbols=[v['element'] for v in ua]
    if [v['element'] for v in a]!=symbols[:192]:raise ValueError('Candidate element order changed')
    if not np.array_equal(xyz[0],origin[0]):raise ValueError('Metal coordinate moved')
    qtext=check(c['vacuum_charges']).read_text().split()
    if len(qtext)!=192:raise ValueError('Candidate population count differs')
    changed=np.where(np.max(np.abs(xyz-origin[:192]),axis=1)>m['settings']['source_coordinate_tolerance_A'])[0].tolist()
    allowed={local_by_physical[i] for i in moving}
    if set(changed)-allowed:raise ValueError('Nonselected core atom changed')
    caps=[x['index'] for x in um['centers'] if x['cap']]
    if len(caps)!=14 or not np.array_equal(xyz[caps],origin[caps]):raise ValueError('Original fixed caps changed')
    centers=[]
    for i,(identity,atom,qt) in enumerate(zip(um['centers'],ua,uq)):
        identity={k:v for k,v in identity.items() if k!='index'}
        if i<192:atom=a[i];qt=qtext[i]
        elif float(qt)!=0 or identity['uncertain'] or identity['cap'] or identity['local_index'] is not None:raise ValueError('Invalid appended environment charge/mask')
        centers.append({**identity,**atom,'charge_text':qt})
    # Verify every represented kinematic physical atom against the union, including fixed flanks.
    source_index={base.key(x['source']):i for i,x in enumerate(centers) if 'source' in x}
    for i,meta in enumerate(mapping['source_atom_metadata']):
        if not meta:continue
        k=base.key(meta)
        if k not in source_index:raise ValueError('Missing physical source atom in union')
        if not np.allclose(physical[i],centers[source_index[k]]['xyz'],atol=m['settings']['kinematic_replay_A'],rtol=0):
            raise ValueError('Moving/fixed physical source reconstruction disagrees')
    if [x['line'] for x in centers[192:]]!=[x['line'] for x in ua[192:]]:raise ValueError('Appended source coordinates changed')
    heavy_local=[i for i in range(192) if symbols[i]!='H'];heavy_env=[i for i in range(192,len(centers)) if symbols[i]!='H']
    ep=np.asarray([centers[i]['xyz'] for i in heavy_env]);nearest=None
    for i in heavy_local:
        distances=np.linalg.norm(ep-xyz[i],axis=1);k=int(np.argmin(distances));d=float(distances[k]);j=heavy_env[k]
        if nearest is None or d<nearest['distance_A']:
            nearest={'local_index':i,'environment_union_index':j,'distance_A':d,'local_identity':um['centers'][i],
                     'environment_identity':um['centers'][j],'origin_distance_A':float(np.linalg.norm(origin[i]-origin[j]))}
    # Same extreme-overlap policy as original local geometry_check, now also
    # checking newly introduced physical contacts to the fixed appended source.
    # No bonded topology or general nonbonded-clash radius is invented here.
    ep_all=np.asarray([v['xyz'] for v in centers[192:]]);new_overlaps=[]
    for i in changed:
        if centers[i]['cap']:continue
        olddist=np.linalg.norm(origin[192:]-origin[i],axis=1)
        newdist=np.linalg.norm(ep_all-xyz[i],axis=1)
        limits=np.asarray([m['settings']['new_overlap_H_A'] if 'H' in (symbols[i],symbols[j]) else m['settings']['new_overlap_heavy_A'] for j in range(192,len(centers))])
        for k in np.where((newdist<limits)&(olddist>=limits))[0].tolist():
            new_overlaps.append({'local_index':i,'environment_union_index':k+192,'before_A':float(olddist[k]),
                                 'after_A':float(newdist[k]),'existing_threshold_A':float(limits[k])})
    cap_distances=[]
    for i in caps:
        cm=centers[i]['cap_mapping'];j=source_index[base.key(cm['omitted'])]
        cap_distances.append({'local_cap_index':i,'omitted_source_union_index':j,'omitted':cm['omitted'],
                              'origin_distance_A':float(np.linalg.norm(origin[i]-origin[j])),
                              'candidate_distance_A':float(np.linalg.norm(xyz[i]-np.asarray(centers[j]['xyz'])))})
    contact_receipt={'minimum_core_environment_heavy_distance':nearest,'cap_source_distances':cap_distances,
                     'new_severe_physical_overlaps':new_overlaps,'threshold_source':'original geometry_check 1.0A heavy / 0.55A H',
                     'cap_contacts_not_physical_clashes':True}
    write(out/'CONTACT_ADMISSION.json',contact_receipt)
    if new_overlaps:raise ValueError('New severe physical core/environment overlap at original thresholds; see CONTACT_ADMISSION.json')
    variants={'local':base.emit(out,'local',centers[:192],check(c['xyz'])),
              'union':base.emit(out,'union',centers),
              'cap_removed':base.emit(out,'cap_removed',[v for v in centers if not v['cap']])}
    if variants['local']['charge_sum_e']!=variants['union']['charge_sum_e']:raise ValueError('Union charge sum differs')
    if variants['union']['masks']['sha256']!=uv['masks']['sha256'] or variants['cap_removed']['masks']['sha256']!=old['variants']['cap_removed']['masks']['sha256']:
        raise ValueError('Original charge/cap mask changed')
    saved={'native_result':c['native_result'],'native_replay':{'vacuum':c['native_vacuum_replay']},
           'native_diagonal_checks':c['native_diagonal_checks'],'vacuum_precision':c['vacuum_precision'],
           'component_v2_acceptance':c['component_v2']}
    return {'status':'prepared','id':req['id'],'sample':req['sample'],'role':c['role'],'metal':metal,'donor_candidate':candidate,
            'source_case_id':c['source_case_id'],'core_charge':c['charge'],'multiplicity':1,'core_atom_count':192,'cap_count':14,
            'origin_row_id':f"sample{req['sample']}_{metal}",'origin_preparation':m['origin_preparation'],'origin_union':uv,
            'state_signature':c['state_signature'],'source_preparation':c['source_preparation'],'protonated_source':old['protonated_source'],'maps':c['maps'],
            'archived_native_result':c['native_result'],'saved_solvent':saved,'variants':variants,
            'geometry_admission':{'pass':True,'original_physical_gate':gate,'replay_max_A':residual,'changed_local_indices':changed,
                                  'all_moving_physical_atoms_in_local_block':True,'all_appended_source_atoms_fixed':True,
                                  'fixed_cap_coordinates_exact':True,'minimum_core_environment_heavy_distance':nearest,
                                  'cap_source_distances':cap_distances,'new_severe_physical_overlaps':new_overlaps,
                                  'contact_admission':pin(out/'CONTACT_ADMISSION.json'),'exact_duplicate_centers':False},
            'zero_environment_charge_is_diagnostic_only':True,'new_energy_calls':0}


def execute(manifest,sha):
    # Allocation guard before importing any numerical/scientific module.
    if not os.environ.get('SLURM_JOB_ID','').isdigit():raise RuntimeError('Scientific preparation requires Slurm')
    mp={'path':str(Path(manifest).resolve()),'sha256':sha};m=read(mp);base=load('original_cavity_preparation',m['base_preparation']);placement=base.allocation()
    if m['protocol']!=PROTOCOL or pin(__file__)!=m['implementation']['PREPARE.py']:raise ValueError('Wrong frozen preparation')
    out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);t=time.monotonic();cpu=time.process_time()
    rows=[{'id':r['id'],'sample':r['sample'],'metal':r['metal'],'donor_candidate':r['donor_candidate'],'status':'unavailable','reason':'not prepared'} for r in m['requests']]
    receipt={'status':'running','manifest':mp,'placement':placement,'new_energy_calls':0};write(out/'START.json',receipt)
    try:
        for p in (m['inventory'],m['origin_preparation'],m['origin_native_collection'],m['plan'],m['kinematics_provenance'],m['python'],*m['implementation'].values()):check(p)
        if pin(sys.executable)!=m['python']:raise ValueError('Interpreter changed')
        for p in m['source_modules'].values():check(p)
        for s in m['input_snapshots'].values():check(s['original']);check(s['snapshot'])
        sys.path.insert(0,m['source_snapshot'])
        import numpy as np
        import adaptive_metal_proposals as metal_geometry
        import mace_site_kinematics
        for module in (metal_geometry,mace_site_kinematics):
            if pin(module.__file__)!=m['source_modules'][Path(module.__file__).name]:raise ValueError('Helper version changed')
        for i,req in enumerate(m['requests']):
            d=out/req['id'];d.mkdir()
            try:rows[i]=prepare(req,d,m,base,np,metal_geometry,mace_site_kinematics.Kinematics)
            except Exception as exc:rows[i].update(reason=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc())
        for sample in (0,3):
            for candidate in ('adaptive_Ca','adaptive_La'):
                paired=[r for r in rows if r['sample']==sample and r['donor_candidate']==candidate]
                if all(r['status']=='prepared' for r in paired):
                    for name in ('local','union','cap_removed'):
                        a,b=[base.xyz_atoms(check(r['variants'][name]['xyz'])) for r in paired]
                        if a[0]['xyz']!=b[0]['xyz'] or [(x['element'],x['xyz']) for x in a[1:]]!=[(x['element'],x['xyz']) for x in b[1:]]:
                            for r in paired:r.update(status='unavailable',reason='Paired geometry differs')
        for r in rows:write(out/r['id']/'PREPARATION.json',r)
        receipt['status']='complete' if all(r['status']=='prepared' for r in rows) else 'failed_preparation'
    except BaseException as exc:
        receipt.update(status='failed',error=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc());raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-t,process_CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out/'EXECUTION.json',receipt)
        write(out/'PREPARATION.json',{'protocol':PROTOCOL,'manifest':mp,'terminal':True,'status':'prepared' if all(r['status']=='prepared' for r in rows) else 'unavailable',
                                    'rows':rows,'declared_nonorigin_cells':8,'prepared_nonorigin_cells':sum(r['status']=='prepared' for r in rows),
                                    'origin_preparation':m['origin_preparation'],'execution':pin(out/'EXECUTION.json'),'new_energy_calls':0})
    if receipt['status']!='complete':raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='command',required=True)
    f=s.add_parser('freeze');f.add_argument('--root',required=True);f.add_argument('--output',required=True)
    e=s.add_parser('execute');e.add_argument('--manifest',required=True);e.add_argument('--manifest-sha',required=True)
    a=p.parse_args()
    if a.command=='freeze':freeze(a.root,a.output)
    else:execute(a.manifest,a.manifest_sha)
