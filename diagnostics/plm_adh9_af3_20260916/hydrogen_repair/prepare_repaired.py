#!/usr/bin/env python3
"""Rebuild exactly the admitted fixed core with explicitly repaired protein H."""
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import gemmi
import repair_hydrogens as helper
from validate_geometry import validate

ROOT=helper.ROOT;BASE=ROOT/'workspaces/plm_adh9_af3_20260916/hydrogen_repair'
OUT=BASE/'prepared';HERE=Path(__file__).resolve().parent
CASE=helper.STEM+'_Hrepair';UID='PQQSEQ_13d74836d4b7a3e02140'
APPROVAL=Path('/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/adh9/af3_comparison/hydrogen_repair/authorization.json')
PROTOCOL='pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3'
read=lambda p:json.loads(Path(p).read_text())

def load_wrapper():
    path=ROOT/'diagnostics/plm_adh9_fixed_core_20260916/prepare_candidates.py'
    if helper.record(path)['sha256']!='8d988aa7bd5e10924f879b464b58f4ad3403829e0b21b7db6e55bfefb304ec5a':raise ValueError('Frozen wrapper changed')
    spec=importlib.util.spec_from_file_location('hydrogen_repair_frozen_core',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module

def structure_atoms(path):
    s=gemmi.read_structure(str(path));assert len(s)==1
    out={}
    for c in s[0]:
        for r in c:
            for a in r:
                k=(c.name,r.name,r.seqid.num,r.seqid.icode,a.name,a.element.name)
                if k in out:raise ValueError('Duplicate structure atom')
                out[k]=(a.pos.x,a.pos.y,a.pos.z)
    return s,out

def main():
    if not APPROVAL.is_file():raise ValueError('Current repair approval missing')
    wrapper=load_wrapper();fixed=wrapper.fixed;cal=ROOT/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/implementation_pins.json'
    frozen=fixed.verify_pins(cal)
    old_manifest=helper.OLD/(helper.STEM+'_carve_manifest.json');old=read(old_manifest)
    if old['source_cif']['sha256']!='f8517d8ed8042e941c5d09d1ca9b927382e0e8712f9ce4ff683a9777e4f1e23d':raise ValueError('Wrong selected AF3 model')
    recovery=BASE/'recovered_exact_objective';minpath=recovery/'minimization.json';minrec=read(minpath)
    repaired_pdb=recovery/'protonated_recovered.pdb'
    if minrec['output']!=helper.record(repaired_pdb):raise ValueError('Recovered PDB differs from converged output')
    structure,newatoms=structure_atoms(repaired_pdb);_,oldatoms=structure_atoms(old['source_structure']['path'])
    if newatoms.keys()!=oldatoms.keys():raise ValueError('Full-protein atom/protonation identities changed')
    if any(newatoms[k]!=v for k,v in oldatoms.items() if k[-1]!='H'):raise ValueError('Full-protein heavy coordinates changed')
    target=read(ROOT/'workspaces/plm_adh9_af3_20260916/hydrogen_retry_01/candidate_manifest.json')['targets'][0]
    keys,residues,partner=wrapper.role_state(structure[0],target)
    site,pkey,pqq_residue,prepared,contacts,sulfur,untyped=wrapper.site_state(structure[0],target)
    fragments=[copy.deepcopy(old['qm_fragments'][0])]
    assert fixed.pqq_atom_records(pqq_residue,prepared)==fragments[0]['atom_records']
    atoms=list(prepared.atoms)
    for role in fixed.ROLE_ORDER:
        fn=fixed.cationic_sidechain_fragment if role=='catalytic_asp_cationic_partner' else fixed.canonical_sidechain_fragment
        # Preserve the frozen helper signatures and exact original atom order.
        if role=='catalytic_asp_cationic_partner':fa,fragment=fn(residues[role],keys[role],site.atom.pos)
        else:fa,fragment=fn(residues[role],keys[role],role,site.atom.pos)
        atoms.extend(fa);fragments.append(fragment)
    validation=validate(old['qm_fragments'],fragments,minrec)
    validation['convergence_evaluated_before_PDB_serialization']=True
    validation['post_serialization_force_measurement']=helper.record(recovery/'post_serialization.json')
    validation['post_serialization']=read(recovery/'post_serialization.json')['post_PDB_serialization']
    OUT.mkdir(parents=True,exist_ok=False);directory=OUT/CASE;directory.mkdir()
    pdb=directory/(CASE+'_protonated.pdb');shutil.copyfile(repaired_pdb,pdb)
    validation.update({'approval':helper.record(APPROVAL),'source_AF3_cif':old['source_cif'],
       'original_failed_carve':helper.record(old_manifest),'repaired_protein':helper.record(pdb),
       'full_structure_atom_identities_preserved':True,'full_structure_heavy_displacement_A':0.0,
       'original_AF3_source_and_selected_sample_preserved':True,
       'source_protonation_states_and_core_formal_charges_preserved':True,
       'calibration_compatibility':'Same core/state/heavy geometry/electronic method; explicitly changed convergence of same H objective. No refit or independent calibration-transfer validation.'})
    valpath=OUT/'geometry_validation.json';helper.write(valpath,validation)
    protonation=directory/(CASE+'_protonation_manifest.json')
    helper.write(protonation,{'protocol_id':'pdbfixer_standard_residue_protonation_v2_with_same_objective_H_recovery',
       'original_protonation':old['protonation_manifest'],'source':old['normalized_source_structure'],
       'output':helper.record(pdb),'ph':7.0,'nonstandard_hydrogens_added':0,'heavy_atoms_changed':0,
       'original_seed':20260914,'new_random_draws':0,'hydrogen_minimization':helper.record(minpath),
       'recovery_script':helper.record(HERE/'recover_captured_objective.py'),
       'changed_numerical_preparation':True,'original_implementation_replayed_exactly':False})
    heavy=directory/(CASE+'_heavy_coordinate_check.json')
    helper.write(heavy,{'passes':True,'source':old['source_cif'],'original_protonated_structure':old['source_structure'],
       'protonated':helper.record(pdb),'all_source_heavy_atoms_retained':True,'no_heavy_atoms_added':True,
       'maximum_movement_from_original_protonated_A':0.0,'original_AF3_serialization_check':old['heavy_coordinate_check']})
    candidate=OUT/'candidate_manifest.json';target=copy.deepcopy(target);target['case_id']=CASE
    helper.write(candidate,{'schema_version':'plm.adh9.candidate_manifest.v1','protocol_id':PROTOCOL,'approval':helper.record(APPROVAL),'targets':[target]})
    pins=OUT/'implementation_pins.json'
    helper.write(pins,{'schema_version':'plm.adh9.hydrogen_repair.pins.v1','protocol_id':PROTOCOL,'approval':helper.record(APPROVAL),
       'candidate_manifest':helper.record(candidate),'original_carve':helper.record(old_manifest),
       'original_calibration_implementation_pins':helper.record(cal),'geometry_validation':helper.record(valpath),
       'recovery_minimization':helper.record(minpath),'captured_original_objective':helper.record(BASE/'original_capture/system_0.xml'),
       'scripts':{p.name:helper.record(p) for p in HERE.glob('*.py')},
       'calibration_replayed_without_deviation':False,'deviation':'Same H objective; restart failed PDB, numerical iteration ceiling50→1000, convergence check.'})
    repaired=copy.deepcopy(old);repaired.update(stem=CASE,status='ready_for_orca',candidate_manifest=helper.record(candidate),
       wrapper=helper.record(__file__),source_structure=helper.record(pdb),protonation_manifest=helper.record(protonation),
       heavy_coordinate_check={**helper.record(heavy),'passes':True},qm_fragments=fragments,
       hydrogen_repair={'approval':helper.record(APPROVAL),'implementation_pins':helper.record(pins),'original_failed_carve':helper.record(old_manifest),
       'minimization':helper.record(minpath),'changed_numerical_preparation':True,'calibration_compatibility':validation['calibration_compatibility']},
       geometry_validation=helper.record(valpath))
    repaired['charge_ledger']['fragments']=fragments
    repaired['fixed_core']['hydrogen_only_preparation_recovery']=True
    tasks=[];outputs={}
    for metal in ('La','Ca'):
        all_atoms=[(metal,site.atom.pos.x,site.atom.pos.y,site.atom.pos.z),*atoms]
        charge=old['charge_ledger'][metal+'_total'];wrapper.base._validate_singlet(metal,all_atoms,charge)
        xyz=directory/(CASE+'_'+metal+'_qm.xyz');inp=directory/('sp_'+CASE+'_'+metal+'.inp')
        fixed.write_xyz(xyz,stem=CASE,label=metal,atoms=all_atoms,charge=charge)
        fixed.write_orca_input(inp,xyz_name=xyz.name,charge=charge,stem=CASE,label=metal)
        xr=helper.record(xyz);xr['path']=xyz.name;ir=helper.record(inp);ir['path']=inp.name
        outputs[metal+'_xyz']=xr;outputs[metal+'_input']=ir
        tasks.append({'task_id':metal,'input':ir,'xyz':xr,'output_path':inp.with_suffix('.out').name})
    la=(directory/(CASE+'_La_qm.xyz')).read_text().splitlines();ca=(directory/(CASE+'_Ca_qm.xyz')).read_text().splitlines()
    if la[3:]!=ca[3:] or la[2].split()[1:]!=ca[2].split()[1:]:raise ValueError('Paired coordinate invariant failed')
    repaired.update(tasks=tasks,outputs=outputs)
    carve=directory/(CASE+'_carve_manifest.json');helper.write(carve,repaired)
    helper.write(OUT/'prepared_pairs.json',{'schema_version':'plm.adh9.prepared_pairs.v1','protocol_id':PROTOCOL,
       'candidate_manifest':helper.record(candidate),'implementation_pins':helper.record(pins),'geometry_validation':helper.record(valpath),
       'cases':[{'case_id':CASE,'target_id':UID,'rank':1,'source_cif':old['source_cif'],'carve_manifest':helper.record(carve),'status':'ready_for_orca','reason':None}],
       'prepared_pair_count':1,'unsupported_count':0,'orca_executed':False})
    print(OUT/'prepared_pairs.json')

if __name__=='__main__':main()
