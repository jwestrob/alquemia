"""Prepare one supplemental source from an already executed fixed-state H repair."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import gemmi
import numpy as np
from scipy.spatial import cKDTree
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,paired
import pqq_fast_prepare as fast
import consistent_context as union
from environment_context_chemistry import complete_expansion
from second_shell_context import parent_state
from mace_hybrid import check_atoms,write_xyz

PROTOCOL='archived_same_potential_H_repair_PQQ_supplemental_source_v1'
CASE='mmol_1770-pqq-la_model__conditioned_La__seed-1_sample-4'


def atoms(path):
    return [(c.name,r.name,r.seqid.num,r.seqid.icode,a.name,a.element.name,tuple(a.pos))
            for c in gemmi.read_structure(str(path))[0] for r in c for a in r]


def prepare(audit,union_preparation,group_selection,agreement,output):
    prior=read_json(audit);rp=prior['repaired_source'];repair=read_json(verify(rp['continuation_receipt']));parent=read_json(verify(repair['parent_receipt']))
    if parent['case_id']!=CASE or not repair['mobile_H_RMS_force_converged'] or not repair['geometry']['passes_geometry']:raise InvalidArtifact('archived H repair unqualified')
    for p in rp.values():
        if isinstance(p,dict) and 'path' in p:verify(p)
    old=atoms(verify(rp['original_protonated']));new=atoms(verify(rp['continued_PDB']))
    if [x[:-1] for x in old]!=[x[:-1] for x in new] or len(old)!=9582 or sum(x[-2]!='H' for x in old)!=4895:raise InvalidArtifact('source atom inventory changed')
    if any(x!=y for x,y in zip(old,new) if x[-2]!='H') or cKDTree([x[-1] for x in new]).query_pairs(.45):raise InvalidArtifact('heavy geometry changed or H still overlaps')
    if repair['source_system']!=rp['captured_System_XML'] or parent['source']!=rp['original_normalized'] or parent['existing_protonated']!=rp['original_protonated']:raise InvalidArtifact('repair/source provenance differs')
    u=read_json(union_preparation);g=read_json(group_selection);before=next(c for c in u['cases'] if c['case_id']==CASE);source=before['source'];group=next(x for x in g['groups'] if x['root_case_id']==source['root_case_id'])
    if before['status']!='prior_unsupported' or before['reason']!='overlapping expanded atoms/caps':raise InvalidArtifact('wrong original failure source')
    canonical=next(c for c in u['cases'] if c['case_id']==group['canonical_case_id']);config=u['config'];out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);local=out/'local';local.mkdir()
    oldpm=verify(rp['original_protonated']).with_name('protonated_protonation_manifest.json');pm=read_json(oldpm)
    if pm['output']!=rp['original_protonated'] or pm['source']!=rp['original_normalized'] or pm['ph']!=7 or pm['repaired_missing_atom_count'] or pm['repaired_missing_terminal_atom_count']:raise InvalidArtifact('original protonation state differs')
    newpm={'protocol_id':PROTOCOL,'output':rp['continued_PDB'],'source':rp['original_normalized'],'ph':pm['ph'],
        'original_manifest':record(oldpm),'repair_receipt':rp['continuation_receipt'],'captured_same_potential_System':rp['captured_System_XML'],
        'original_protonation_metadata':pm,'new_hydrogens_added':0,'new_force_calls':0,'atom_inventory_unchanged':True,'heavy_coordinates_unchanged':True,
        'repair_changes_only_existing_H_coordinates':True,'repaired_missing_atom_count':0,'repaired_missing_terminal_atom_count':0}
    pp=local/'protonation_repair_manifest.json';write_new(pp,newpm)
    fixed,_,_,_,_=fast.helpers(config)
    core=fast.build_core(source,verify(rp['continued_PDB']),pp,verify(rp['original_normalized']),local,config,fixed)
    oldcore=read_json(verify(rp['original_protonated']).with_name('core_preparation.json'));cp=read_json(verify(core['parent']))
    if cp['fixed_core']!=oldcore['fixed_core'] or cp['charge_ledger']!=oldcore['charge_ledger']:raise InvalidArtifact('core state changed')
    for z in ('Ca','La'):
        a=xyz(verify(oldcore['outputs'][z+'_xyz']));b=xyz(verify(core['endpoints'][z]['xyz']))
        if len(a)!=len(b) or [v[0] for v in a]!=[v[0] for v in b] or any(x!=y for x,y in zip(a,b) if x[0]!='H'):raise InvalidArtifact('core heavy geometry/composition changed')
    for oldfrag,newfrag in zip(oldcore['qm_fragments'],cp['qm_fragments']):
        if (oldfrag['id'],oldfrag['formal_charge'],oldfrag['atom_count'])!=(newfrag['id'],newfrag['formal_charge'],newfrag['atom_count']):raise InvalidArtifact('fragment state changed')
        for a,b in zip(oldfrag['atom_records'],newfrag['atom_records']):
            if a['origin']!='source_protonation_hydrogen' and a!=b:raise InvalidArtifact('cofactor or cap recipe changed')
    state=parent_state(core,config['topology'],require_endpoint_receipts=False);expanded,detail=complete_expansion(state);ap=local/'context_preparation.json';write_new(ap,detail);eps={}
    for z in ('Ca','La'):
        q=core['endpoints'][z]['charge']+detail['added_formal_charge'];check_atoms(expanded[z],q);p=local/(z+'_context.xyz');write_xyz(p,expanded[z]);eps[z]={'xyz':record(p),'charge':q,'multiplicity':1,'metal_index':0}
    paired(verify(eps['La']['xyz']),verify(eps['Ca']['xyz']),eps['La']['charge'],eps['Ca']['charge'])
    row={'case_id':CASE,'status':'prepared','source':source,'core':core,'representations':{'context':{'protocol_id':fast.CONTEXT_PROTOCOL,'preparation':record(ap),'endpoints':eps}},
        'supplemental_preparation_protocol':PROTOCOL,'prior_failure':before,'archived_H_repair':rp,'new_protonation_or_optimization':False}
    write_new(local/'preparation_result.json',row)
    united=union.prepare_one(row,group,config,str(out/'union'))
    if united['status']=='prepared_awaiting_group_check':
        if united['state_key']!=canonical['state_key'] or united['protein_key']!=canonical['protein_key']:united.update(status='state_or_protein_mismatch',reason='supplemental state differs from existing canonical union')
        else:united['status']='prepared'
    result={'protocol_id':PROTOCOL,'agreement':record(agreement),'audit':record(audit),'existing_union_preparation':record(union_preparation),'existing_group_selection':record(group_selection),
        'frozen_group':group,'canonical_state_signature':canonical.get('state_signature'),'local':row,'union':united,'status':'prepared' if united['status']=='prepared' else 'unsupported',
        'original_failure_preserved':True,'historical225_changed':False,'production_changed':False,'new_molecular_calls':0,'implementation':record(__file__)}
    write_new(out/'PREPARATION.json',result);return {'status':result['status'],'preparation':record(out/'PREPARATION.json'),'local_atoms':len(expanded['Ca']),'union_atoms':united.get('atom_count'),'union_reason':united.get('reason'),'new_molecular_calls':0}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('audit','union_preparation','group_selection','agreement','output'):p.add_argument('--'+n.replace('_','-'),required=True,type=Path)
    print(json.dumps(prepare(**vars(p.parse_args())),indent=2))
