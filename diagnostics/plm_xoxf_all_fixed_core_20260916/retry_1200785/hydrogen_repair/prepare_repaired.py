#!/usr/bin/env python3
"""Recover retained3cad H coordinates under the exact captured objective; no ORCA."""
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','OPENMM_CPU_THREADS'):
    os.environ[key]='1'
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import prepare_batch as b
from geometry_checks import check
import gemmi
import numpy as np
import openmm as mm
from openmm import app,unit

A=b.A;UID='PQQSEQ_3cad7f010e2b9d9374c4';CID=UID+'_AF3_sample0'
BASE=A/'workspaces/plm_xoxf_all_fixed_core_20260916/retry_1200785'
OLD=BASE/'prepared_batch';OUT=BASE/'hydrogen_repair'
EQ=Path('/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue')
APPROVAL=EQ/'xoxf_all/authorization.json'
SOURCE=OLD/'original_preparation'/CID/f'{CID}_protonated.pdb'
SYSTEM=OLD/'hydrogen_audit'/CID/'system_0.xml'
PM=SOURCE.with_name(f'{CID}_protonation_manifest.json')
SOURCE_SHA='40d3ff5d44bfe9c4fbd8615927b67486a21fb3d8fc32925dc4f2b4aa89147c88'
SYSTEM_SHA='5c094007d008f08685cc46f78174b44e9940b425154e6380576dea12fde8e435'
CIF_SHA='2c26279c776d06703f5f38e163d2caf774f24820535eae2c877acee1f947f76b'
record=b.record;read=b.read;write=b.write
h=b.load('xoxf_retained_H_helpers',A/'diagnostics/plm_adh9_af3_20260916/hydrogen_repair/repair_hydrogens.py')

def atoms_by_identity(path):
    structure=gemmi.read_structure(str(path));result={}
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                key=(chain.name,residue.name,residue.seqid.num,residue.seqid.icode.strip(),atom.name,atom.element.name)
                if key in result:raise ValueError('Duplicate atom identity')
                result[key]=(atom.pos.x,atom.pos.y,atom.pos.z)
    return structure,result

def recover():
    b.verify({'path':str(SOURCE),'sha256':SOURCE_SHA});b.verify({'path':str(SYSTEM),'sha256':SYSTEM_SHA})
    pdb=app.PDBFile(str(SOURCE));atoms=list(pdb.topology.atoms())
    mask=np.asarray([a.element==app.element.hydrogen for a in atoms])
    system=mm.XmlSerializer.deserialize(SYSTEM.read_text())
    if system.getNumParticles()!=len(atoms):raise ValueError('Captured System atom count differs')
    if any((system.getParticleMass(i).value_in_unit(unit.dalton)>0)!=v for i,v in enumerate(mask)):
        raise ValueError('Captured objective does not hold all non-H atoms fixed')
    ctx=mm.Context(system,mm.VerletIntegrator(0),mm.Platform.getPlatformByName('CPU'),{'Threads':'1'})
    ctx.setPositions(pdb.positions);before,bm=h.state(ctx,mask);reporter=h.Reporter();start=time.monotonic()
    mm.LocalEnergyMinimizer.minimize(ctx,1.0,1000,reporter)
    after,am=h.state(ctx,mask);elapsed=time.monotonic()-start
    np.save(OUT/'recovered_full_precision_positions_nm.npy',after)
    identical_heavy=np.array_equal(before[~mask],after[~mask])
    path=OUT/'protonated_recovered.pdb'
    with path.open('x') as f:app.PDBFile.writeFile(pdb.topology,after*unit.nanometer,f,keepIds=True)
    serialized=app.PDBFile(str(path));ctx.setPositions(serialized.positions);rounded,rm=h.state(ctx,mask)
    _,old=atoms_by_identity(SOURCE);_,new=atoms_by_identity(path)
    if old.keys()!=new.keys():raise ValueError('Full protein atom or protonation identity changed')
    if any(v!=new[k] for k,v in old.items() if k[-1]!='H'):raise ValueError('Serialized heavy coordinates changed')
    result={'status':'PASS' if identical_heavy and am['hydrogen_force_rms_kJ_mol_nm']<=1.0 else 'FAIL',
      'approval':record(APPROVAL),'source_protonated':record(SOURCE),'captured_original_System':record(SYSTEM),
      'original_minimization':record(OLD/'hydrogen_audit'/CID/'minimization_observation.json'),
      'initial':bm,'final_full_precision':am,'post_PDB_serialization':rm,
      'tolerance':1.0,'maximum_iterations':1000,'iterations_reported':len(reporter.entries),
      'platform':'CPU','threads':1,'elapsed_seconds':elapsed,'original_random_seed':20260914,'new_random_draws':0,
      'objective':'Exact deserialized original System XML; no added/deleted/changed forces, particles or parameters',
      'all_full_protein_atom_identities_preserved':True,'full_precision_heavy_coordinates_identical':bool(identical_heavy),
      'serialized_heavy_coordinates_identical':True,'heavy_displacement_A':0.0,
      'protonation_state_changed':False,'output':record(path),'full_precision_positions':record(OUT/'recovered_full_precision_positions_nm.npy'),
      'script':record(__file__),'observer_helpers':record(h.__file__),'trace':reporter.entries,
      'preparation_deviation':'Retained failed PDB H-only restart; unchanged captured objective/tolerance; iteration ceiling50→1000; no fresh protonation',
      'calibration_compatibility_limit':'Same fixed core, state, heavy coordinates and electronic method; numerical H convergence changed explicitly. No calibration rerun or refit.'}
    write(OUT/'recovery_validation.json',result)
    print(json.dumps({k:result[k] for k in ('status','iterations_reported','elapsed_seconds','final_full_precision','post_PDB_serialization')}),flush=True)
    if result['status']!='PASS':raise ValueError('Same-objective numerical recovery failed convergence/identity')
    return path,result

def validate_frozen_atoms(w,target,carve):
    old_model,old=atoms_by_identity(SOURCE);_,new=atoms_by_identity(carve['source_structure']['path'])
    if old.keys()!=new.keys():raise ValueError('Atom states changed')
    roles,residues,partner=w.role_state(old_model[0],target)
    site,pkey,pqq,prepared,contacts,*_=w.site_state(old_model[0],target)
    expected_pqq=w.fixed.pqq_atom_records(pqq,prepared)
    if carve['qm_fragments'][0]['atom_records']!=expected_pqq:raise ValueError('PQQ coordinates/state changed')
    fixed_nonH=0;fixed_caps=0;changed_H=0;parent_max=0.0
    for fragment in carve['qm_fragments'][1:]:
        role=fragment['role'];residue=residues[role];key=roles[role]
        old_atoms={a.name:a for a in residue}
        for atom in fragment['atom_records']:
            xyz=tuple(atom['xyz_A']);origin=atom['origin']
            if origin=='source_heavy_atom':
                a=old_atoms[atom['name']]
                if xyz!=(a.pos.x,a.pos.y,a.pos.z):raise ValueError('Frozen core heavy coordinate changed')
                fixed_nonH+=1
            elif origin=='generated_Cbeta_link_cap':
                if xyz!=tuple(w.base._sidechain_link_h(residue,key)[1:]):raise ValueError('Frozen cap changed')
                fixed_caps+=1
            elif origin=='source_protonation_hydrogen':
                a=old_atoms[atom['name']]
                if a.element.name!='H':raise ValueError('New or renamed H atom')
                changed_H+=xyz!=(a.pos.x,a.pos.y,a.pos.z)
                parent=next(x for x in fragment['atom_records'] if x['name']==atom['parent_atom'])
                distance=float(np.linalg.norm(np.asarray(xyz)-np.asarray(parent['xyz_A'])))
                if distance>1.25:raise ValueError('Retained strict carver1.25A H-parent gate failed')
                parent_max=max(parent_max,distance)
            else:raise ValueError('Unexpected protein atom origin')
    geometry=check(carve['qm_fragments'])
    if geometry['status']!='PASS':raise ValueError('Repaired core geometry fails')
    return {'geometry':geometry,'unchanged_PQQ_atom_count':len(expected_pqq),'unchanged_protein_heavy_count':fixed_nonH,
       'unchanged_cap_count':fixed_caps,'changed_source_protein_H_count':changed_H,'protein_H_parent_max_A':parent_max,
       'strict_carver_H_parent_max_A':1.25,'all_full_protein_H_identities_preserved':True,'PQQ_chemistry_changed':False,
       'core_membership_changed':False,'sample_substituted':False}

def main():
    if record(b.__file__)['sha256']!='762ec317cd01c99e1940bcb4add11384b6cabcf33358d7229a10fe12cb59ec96':
        raise ValueError('Pinned batch adapter changed')
    w=b.wrapper();original_pins,original,cp=w.verify_pins(OLD/'candidate_implementation_pins.json')
    target=copy.deepcopy(next(t for t in read(OLD/'candidate_manifest.json')['targets'] if t['target_id']==UID))
    if target['rank']!=0 or target['source_cif']['sha256']!=CIF_SHA:raise ValueError('Selected source changed')
    OUT.mkdir(parents=True,exist_ok=False);(OUT/'prepared').mkdir()
    repaired,recovery=recover()
    candidate=OUT/'candidate_manifest.json'
    target['source_provenance'].update(original_failed_protonation=record(PM),H_recovery=record(OUT/'recovery_validation.json'))
    write(candidate,{'schema_version':'plm.adh9.candidate_manifest.v1','protocol_id':b.PROTOCOL,'approval':record(APPROVAL),'targets':[target]})
    candidate_pins=copy.deepcopy(original_pins);candidate_pins['candidate_manifest']=record(candidate)
    candidate_pins['authorities'].update(single_target_repair_script=record(__file__),H_recovery=record(OUT/'recovery_validation.json'))
    pins_path=OUT/'candidate_implementation_pins.json';write(pins_path,candidate_pins)
    pins,original,cp=w.verify_pins(pins_path)
    old_pm=read(PM)
    def retained_H(normalized,output,pins_path,ph):
        if ph!=old_pm['ph']:raise ValueError('Protonation pH state changed')
        _,old_heavy=atoms_by_identity(old_pm['source']['path']);_,new_heavy=atoms_by_identity(normalized)
        if old_heavy!=new_heavy:raise ValueError('Normalized source differs from failed original')
        shutil.copyfile(repaired,output)
        pm={'schema_version':1,'protocol_id':'retained_standard_only_H_same_objective_numerical_recovery_v1',
            'source':record(normalized),'output':record(output),'original_protonation':record(PM),
            'recovery':record(OUT/'recovery_validation.json'),'implementation':record(__file__),
            'software':old_pm['software'],'repaired_missing_atom_count':0,'repaired_missing_terminal_atom_count':0,
            'ph':ph,'new_protonation_performed':False,'protonation_states_changed':False,
            'original_frozen_protonator_reexecuted':False,'numerical_H_preparation_changed':True}
        path=Path(output).with_name(CID+'_protonation_manifest.json');write(path,pm);return path
    original_function=w.protonator.protonate_standard_only
    w.protonator.protonate_standard_only=retained_H
    try:raw=w.prepare_one(target,OUT/'prepared',pins,original,cp)
    finally:w.protonator.protonate_standard_only=original_function
    carve=read(raw['path']);details=validate_frozen_atoms(w,target,carve)
    for rec in carve['outputs'].values():rec['path']=str(Path(raw['path']).parent/rec['path'])
    la=Path(carve['outputs']['La_xyz']['path']).read_text().splitlines();ca=Path(carve['outputs']['Ca_xyz']['path']).read_text().splitlines()
    if la[3:]!=ca[3:] or la[2].split()[1:]!=ca[2].split()[1:]:raise ValueError('La/Ca nuclei differ')
    details.update(status='PASS',target_id=UID,selected_sample=0,source_cif=target['source_cif'],
       original_failed_PDB=record(SOURCE),recovery=record(OUT/'recovery_validation.json'),
       same_captured_System=record(SYSTEM),full_protein_heavy_coordinates_unchanged=True,protonation_states_unchanged=True,
       calibration_implementation_changed=True,preparation_mode='retained_H_same_objective_numerical_recovery',
       charges={m:carve['charge_ledger'][m+'_total'] for m in ('La','Ca')},paired_coordinates_identical=True,
       original_failed_H_parent_distance_A=1.2525230536800513,original_failed_full_precision_distance_A=1.252802393752406,
       limitation=recovery['calibration_compatibility_limit'])
    validation=OUT/'geometry_validation.json';write(validation,details)
    carve.update(geometry_validation=record(validation),hydrogen_recovery=details,
                 preparation_adapter=record(__file__),original_failed_protonation=record(PM))
    final_carve=OUT/'carve_manifest.json';write(final_carve,carve)
    outcomes=OUT/'target_outcomes.json';write(outcomes,{'per_target':[{'target_id':UID,'status':'ready','reason':None}],
        'target_count':1,'orca_executed':False,'untouched_ready_cases':134,'UNK_targets_remain_unscored':True})
    implementation=OUT/'implementation_pins.json';write(implementation,{'protocol_id':b.PROTOCOL,'approval':record(APPROVAL),
        'candidate_pins':record(pins_path),'candidate_manifest':record(candidate),'geometry_validation':record(validation),
        'target_outcomes':record(outcomes),'recovery':record(OUT/'recovery_validation.json'),'repair_script':record(__file__),
        'original_batch_adapter':record(b.__file__),'geometry_checker':record(HERE.parent/'geometry_checks.py'),
        'H_recovery_helpers':record(h.__file__),'original_calibration_implementation_pins':record(cp)})
    case={k:target[k] for k in ('case_id','target_id','rank','source_cif')}
    case.update(status='ready_for_orca',reason=None,carve_manifest=record(final_carve))
    result={'schema_version':'plm.adh9.prepared_pairs.v1','protocol_id':b.PROTOCOL,'cases':[case],
        'candidate_manifest':record(candidate),'implementation_pins':record(implementation),
        'geometry_validation':record(validation),'target_outcomes':record(outcomes),'prepared_pair_count':1,
        'target_count':1,'unsupported_count':0,'orca_executed':False}
    write(OUT/'prepared_pairs.json',result)
    print(json.dumps({'status':'PASS','prepared_pairs':record(OUT/'prepared_pairs.json'),'geometry':details['geometry'],
                      'protein_H_parent_max_A':details['protein_H_parent_max_A']}),flush=True)

if __name__=='__main__':main()
