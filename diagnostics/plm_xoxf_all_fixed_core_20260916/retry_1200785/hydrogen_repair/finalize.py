#!/usr/bin/env python3
"""Finalize the already-recovered core using exact frozen coordinate serialization.
No minimization, new protonation, structural prediction or ORCA is performed.
"""
from continuation import *

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
                expected=w.fixed.atom_record(name=atom['name'],element='H',origin=origin,parent_atom='CB',xyz=w.base._sidechain_link_h(residue,key)[1:])
                if atom!=expected:raise ValueError('Frozen cap differs under canonical6decimal serialization')
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
    w=b.wrapper();pins,original,cp=w.verify_pins(OUT/'candidate_implementation_pins.json')
    target=read(OUT/'candidate_manifest.json')['targets'][0]
    if target['target_id']!=UID or target['rank']!=0 or target['source_cif']['sha256']!=CIF_SHA:raise ValueError('Scope changed')
    recovery=read(OUT/'recovery_validation.json')
    if recovery['status']!='PASS':raise ValueError('H recovery did not pass')
    repaired=b.verify(recovery['output']);candidate=OUT/'candidate_manifest.json';pins_path=OUT/'candidate_implementation_pins.json'
    raw=record(OUT/'prepared'/CID/f'{CID}_carve_manifest.json')
    carve=read(raw['path']);details=validate_frozen_atoms(w,target,carve)

    for rec in carve['outputs'].values():rec['path']=str(Path(raw['path']).parent/rec['path'])
    la=Path(carve['outputs']['La_xyz']['path']).read_text().splitlines();ca=Path(carve['outputs']['Ca_xyz']['path']).read_text().splitlines()
    if la[3:]!=ca[3:] or la[2].split()[1:]!=ca[2].split()[1:]:raise ValueError('La/Ca nuclei differ')
    details.update(status='PASS',target_id=UID,selected_sample=0,source_cif=target['source_cif'],
       original_failed_PDB=record(SOURCE),recovery=record(OUT/'recovery_validation.json'),
       same_system_objective=True,all_heavy_atoms_preserved=True,pqq_and_caps_preserved=True,only_protein_H_coordinates_changed=True,
       recovered_protonated_output=record(repaired),prepared_protonated_output=carve['source_structure'],
       source_System=record(SYSTEM),source_restart_positions=record(FIRST/'recovered_full_precision_positions_nm.npy'),
       same_captured_System=record(SYSTEM),full_protein_heavy_coordinates_unchanged=True,protonation_states_unchanged=True,
       calibration_implementation_changed=True,preparation_mode='retained_H_same_objective_numerical_recovery',
       charges={m:carve['charge_ledger'][m+'_total'] for m in ('La','Ca')},paired_coordinates_identical=True,
       original_failed_H_parent_distance_A=1.2525230536800513,original_failed_full_precision_distance_A=1.252802393752406,
       limitation=recovery['calibration_compatibility_limit'])
    metadata_atoms=[a for f in carve['qm_fragments'] for a in f['atom_records']]
    xyz_atoms=[(v.split()[0],tuple(map(float,v.split()[1:]))) for v in la[3:]]
    if xyz_atoms!=[(a['element'],tuple(a['xyz_A'])) for a in metadata_atoms]:raise ValueError('Core XYZ differs from canonical atom metadata')
    details.update(finalizer=record(__file__),executed_recovery_script=record(HERE/'continuation.py'),
       cap_comparison_policy='Exact canonical atom_record and XYZ six-decimal serialization; no tolerance relaxation',
       full_precision_H_force_rms_kJ_mol_nm=recovery['final_full_precision']['hydrogen_force_rms_kJ_mol_nm'],
       post_PDB_serialization_H_force_rms_kJ_mol_nm=recovery['post_PDB_serialization']['hydrogen_force_rms_kJ_mol_nm'],
       serialized_geometry_force_converged_not_claimed=True)
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
        'executed_recovery_script':record(HERE/'continuation.py'),'first_attempt_script':record(HERE/'prepare_repaired.py'),
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
