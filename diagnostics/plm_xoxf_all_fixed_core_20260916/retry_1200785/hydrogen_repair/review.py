#!/usr/bin/env python3
"""Read-only integrity checks on the one recovered preparation; writes receipt only."""
import hashlib
import json
from pathlib import Path
import sys
import finalize as p

def main():
    source=p.OUT/'prepared_pairs.json';pp=p.read(source)
    if len(pp['cases'])!=1 or pp['cases'][0]['target_id']!=p.UID:raise ValueError('One-target scope differs')
    cm=p.read(p.b.verify(pp['cases'][0]['carve_manifest']));gv=p.read(p.b.verify(cm['geometry_validation']))
    rec=p.read(p.b.verify(gv['recovery']))
    if rec['status']!='PASS' or gv['status']!='PASS':raise ValueError('Preparation not validated')
    if pp['cases'][0]['source_cif']['sha256']!=p.CIF_SHA:raise ValueError('Selected model changed')
    if rec['captured_original_System']['sha256']!=p.SYSTEM_SHA or rec['source_protonated']['sha256']!=p.SOURCE_SHA:
        raise ValueError('Captured objective/restart source changed')
    for key in ('captured_original_System','source_protonated','full_precision_positions','output'):
        p.b.verify(rec[key])
    before,old=p.atoms_by_identity(p.SOURCE);after,new=p.atoms_by_identity(cm['source_structure']['path'])
    if old.keys()!=new.keys():raise ValueError('Full protein atom identities changed')
    heavy_keys=[k for k in old if k[-1]!='H']
    if any(old[k]!=new[k] for k in heavy_keys):raise ValueError('Heavy coordinates changed')
    w=p.b.wrapper();target=p.read(p.OUT/'candidate_manifest.json')['targets'][0]
    keys,residues,_=w.role_state(before[0],target);site,*_=w.site_state(before[0],target)
    role='catalytic_asp_cationic_partner'
    old_failed=False
    try:w.fixed.cationic_sidechain_fragment(residues[role],keys[role],site.atom.pos)
    except w.fixed.FixedCoreError as exc:
        old_failed='hydrogen-parent counts' in str(exc)
    if not old_failed:raise ValueError('Original retained failure was not reproduced')
    keys,residues,_=w.role_state(after[0],target);site,*_=w.site_state(after[0],target)
    _,fragment=w.fixed.cationic_sidechain_fragment(residues[role],keys[role],site.atom.pos)
    details=p.validate_frozen_atoms(w,target,cm)
    if fragment['formal_charge']!=1 or details['protein_H_parent_max_A']>1.25:raise ValueError('Recovered Arg state/gate failed')
    sys.path.insert(0,str(p.EQ/'xoxf_all/retry_1200785'))
    import run_quantum as q
    parent=q.check_prepared(source,{'target_ids':[p.UID],'reused_target_ids':[]},q.pipeline())
    if parent['status']!='PASS':raise ValueError('Parent endpoint prerequisites fail')
    result={'status':'PASS','target_id':p.UID,'source_cif':pp['cases'][0]['source_cif'],
      'prepared_pairs':p.record(source),'geometry_validation':p.record(p.OUT/'geometry_validation.json'),
      'recovery_validation':p.record(p.OUT/'recovery_validation.json'),'same_system_objective':True,
      'all_heavy_atoms_preserved':True,'pqq_and_caps_preserved':True,'only_protein_H_coordinates_changed':True,
      'all_protein_H_identities_preserved':True,'original_carver_failure_reproduced':True,
      'frozen_H_attachment_cutoff_A':1.25,'actual_max_core_H_parent_A':details['protein_H_parent_max_A'],
      'full_protein_atom_count':len(old),'full_protein_heavy_count':len(heavy_keys),
      'changed_full_protein_H_coordinates':sum(old[k]!=new[k] for k in old if k[-1]=='H'),
      'parent_preflight':parent,'review_script':p.record(__file__),'endpoint_submission_performed':False,
      'new_PDBFixer_protonation_performed':False,'controls_or_other_targets_touched':False,
      'limitation':rec['calibration_compatibility_limit']}
    p.write(p.OUT/'independent_validation.json',result)
    print(json.dumps({k:result[k] for k in ('status','actual_max_core_H_parent_A','full_protein_atom_count','changed_full_protein_H_coordinates')}))

if __name__=='__main__':main()
