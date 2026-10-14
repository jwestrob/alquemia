#!/usr/bin/env python3
"""Reuse completed4088 protonation; exact single-OXT crystal completion only."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback

ROOT=Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
HERE=Path(__file__).resolve().parent
PROTOCOL='ExaA_1FLG_AB_reused_single_remote_terminal_OXT_v3'
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,('pin changed',p['path']);return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m


def freeze(output):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    prior=ROOT/'workspaces/exaa_source_packing_20261014/crystal_v2/manifest.json';m=read(prior)
    final=prior.parent/'execution/FINAL_PREPARATION.json';f=read(final)
    assert f['terminal'] and f['supported']==0 and f['denominator']==2
    rows=[]
    for spec in m['chains']:
        p=prior.parent/'execution'/('chain_'+spec['chain'])
        files={k:pin(p/v) for k,v in {'source_audit':'SOURCE_AUDIT.json','mapped_source':'selected_dry_mapped_source.cif',
          'normalized':'source_preparation/normalized.pdb','protonated':'source_preparation/protonated.pdb',
          'protonation':'source_preparation/protonated_protonation_manifest.json','failed_admission':'ADMISSION.json','failure':'FAILURE.txt'}.items()}
        pm=read(files['protonation']);assert pm['repaired_missing_atom_count']==0 and pm['repaired_missing_terminal_atom_count']==1
        assert pm['source']==files['normalized'] and pm['output']==files['protonated']
        rows.append({'spec':spec,'files':files})
    result={'protocol':PROTOCOL,'implementation':{n:pin(HERE/n) for n in ('RECOVER.py','PLAN.md','RUN.sbatch')},
      'prior_manifest':pin(prior),'prior_result':pin(final),'prior_implementation':m['implementation'],
      'python':m['python'],'cases':rows,'output':str(out/'execution'),'reuse_protonation':True,
      'allowed_heavy_addition':{'chain':'A','resnum':616,'icode':'','resname':'TRP','atom':'OXT','element':'O'},
      'terminal_boundary_A':4.3,'new_energy_calls':0,'new_protonation_calls':0}
    for p in m['implementation'].values():verify(p)
    write(out/'manifest.json',result);print(json.dumps({'manifest':pin(out/'manifest.json'),'status':'Slurm_ready_UNRUN'}))


def heavy(path,gemmi):
    st=gemmi.read_structure(str(path));assert len(st)==1
    result={}
    for c in st[0]:
        for r in c:
            for a in r:
                if a.element.name in ('H','D'):continue
                key=(c.name,r.seqid.num,r.seqid.icode.strip(),r.name,a.name.strip(),a.element.name,a.altloc.strip('\x00 '))
                assert key not in result,('duplicate heavy identity',key)
                result[key]=tuple(a.pos)
    return result,st


def validate_saved(m,case,out,gemmi):
    files=case['files'];pm=read(files['protonation'])
    assert pm['source']==files['normalized'] and pm['output']==files['protonated']
    assert pm['repaired_missing_atom_count']==0 and pm['repaired_missing_terminal_atom_count']==1,'unexpected repair counts'
    assert not pm['add_missing_residues'] and pm['detected_missing_residues']=={} and pm['ph']==7.0
    before,_=heavy(verify(files['normalized']),gemmi);after,st=heavy(verify(files['protonated']),gemmi)
    allowed=m['allowed_heavy_addition'];key=tuple(allowed[k] for k in ('chain','resnum','icode','resname','atom','element'))+('',)
    assert set(after)-set(before)=={key},'heavy additions are not exactly terminal A:TRP616/OXT'
    assert set(before)<=set(after),'deposited heavy atoms lost'
    assert all(before[k]==after[k] for k in before),'deposited heavy coordinates moved'
    counts={'LA':0,'CA':0,'PQQ':0};nh={k:0 for k in counts}
    for c in st[0]:
        for r in c:
            if r.name in counts:
                counts[r.name]+=1;nh[r.name]+=sum(a.element.name in ('H','D') for a in r)
                if r.name in ('LA','CA'):assert len(r)==1 and r[0].element.name=='Ca'
    assert counts=={'LA':1,'CA':1,'PQQ':1} and not any(nh.values()),'nonstandard occupancy/H changed'
    source,_=heavy(verify(files['mapped_source']),gemmi)
    assert source==before,'mapped deposited source differs from saved normalized heavy coordinates'
    record={'status':'validated_single_terminal_OXT_before_context','prior_failure':files['failed_admission'],
      'prior_protonation':files['protonation'],'prior_manifest':m['prior_manifest'],'added_atom':{**allowed,'xyz_A':list(after[key])},
      'deposited_heavy_count':len(before),'deposited_heavy_coordinates_exactly_unchanged':True,
      'other_heavy_additions':[],'removed_heavy_atoms':[],'nonstandard_residues':counts,'nonstandard_H':nh,
      'new_protonation_calls':0,'original_missing_terminal_count_preserved':1}
    write(out/'OXT_SOURCE_CHECK.json',record);return record,st


def terminal_check(record,st,row,out,gemmi):
    import numpy as np
    from affordable_common import xyz
    p=read(row['representations']['context']['preparation'])
    def contains(value):
        if isinstance(value,dict):
            if value.get('chain')=='A' and value.get('resnum')==616:return True
            return any(contains(x) for x in value.values())
        if isinstance(value,list):return any(contains(x) for x in value)
        return False
    assert not contains(p['mapping']),'terminal residue participates in context/source/bond/cap mapping'
    points=np.asarray([x[1:] for x in xyz(verify(row['representations']['context']['endpoints']['Ca']['xyz'])) if x[0] not in ('H','D')])
    terminal=[a for c in st[0] if c.name=='A' for r in c if r.seqid.num==616 for a in r if a.element.name not in ('H','D')]
    assert terminal and any(a.name.strip()=='OXT' for a in terminal)
    distances={a.name.strip():float(np.min(np.linalg.norm(points-np.asarray(list(a.pos)),axis=1))) for a in terminal}
    assert min(distances.values())>4.3,'terminal residue is within frozen context boundary'
    record.update(status='passed',terminal_absent_from_context_and_caps=True,terminal_heavy_nearest_context_A=distances,
      minimum_terminal_to_context_A=min(distances.values()),OXT_to_context_A=distances['OXT'],boundary_A=4.3,
      context_preparation=row['representations']['context']['preparation'])
    write(out/'TERMINAL_COMPLETION.json',record);return pin(out/'TERMINAL_COMPLETION.json')


def run(m,old,case,out,v2,helper):
    out.mkdir(parents=True,exist_ok=False);spec=case['spec'];files=case['files'];cid='1FLG_chain_'+spec['chain']+'__mapped_Q9_conditional_oxidized_v1'
    result={'case_id':cid,'chain':spec['chain'],'status':'unsupported','mapping_eligible':False,'paired_maps_equal':False,
      'chemical_state':old['chemical_state'],'source_redox_experimentally_qualified':False,'prior_failed_admission':files['failed_admission']}
    try:
        import gemmi,pqq_fast_prepare as fast,motion_envelope as envelope
        helper.imported_pins({'repository_script_pins':old['repository_script_pins']})
        record,st=validate_saved(m,case,out,gemmi)
        request=copy.deepcopy(old['reference_request']);source=request['source']
        for k in list(source):
            if k.startswith('archived_') or k in ('source_geometry_audit','original_source_id'):source.pop(k)
        active={'chain':'B','resnum':1,'icode':'','resname':'LA','atom':'CA','element':'Ca'}
        spectator={'chain':'D','resnum':1,'icode':'','resname':'CA','atom':'CA','element':'Ca'}
        source.update(case_id=cid,source_structure=files['mapped_source'],metal=active,raw_source_metal=active,
          original_source_structure=old['source'],source_kind='crystal_conditional_oxidized_PQQ',source_conditioning_metal=None,
          source_occupancy_state='deposited_native_2Ca',terminal_completion_source_check=pin(out/'OXT_SOURCE_CHECK.json'),
          assembly={'kind':'explicit_selected_crystal_chain_with_mature_sequence_alias_map','raw_chain':spec['chain'],'raw_source':old['source'],'atom_map':files['source_audit']})
        request.update(source=source,original_source_id=cid,audit_row={'spectator_Ca':spectator})
        fixed,*_=fast.helpers(old['config']);prepdir=out/'source_preparation';prepdir.mkdir()
        core=fast.build_core(source,verify(files['protonated']),verify(files['protonation']),verify(files['normalized']),prepdir,old['config'],fixed)
        state=fast.parent_state(core,old['config']['topology'],require_endpoint_receipts=False)
        expanded,audit=fast.complete_expansion(state);ap=prepdir/'context_preparation.json';write(ap,audit);eps={}
        for z in ('Ca','La'):
            xp=prepdir/(z+'_context.xyz');fast.write_xyz(xp,expanded[z]);charge=core['endpoints'][z]['charge']+audit['added_formal_charge'];fast.check_atoms(expanded[z],charge)
            eps[z]={'xyz':pin(xp),'charge':charge,'multiplicity':1,'metal_index':0}
        prepared={'case_id':cid,'source':source,'status':'prepared','role':'same_protein_crystal_conditional_state',
          'biological_group':'q9z4j7','expected_class':'Ca','label_scope':source['label_scope'],'all_evidence_consumed':True,
          'core':core,'representations':{'context':{'protocol_id':fast.CONTEXT_PROTOCOL,'preparation':pin(ap),'endpoints':eps}},
          'protonation_reused':files['protonation'],'original_source_audit':files['source_audit'],
          'terminal_completion_policy':PROTOCOL,'source_preservation':pin(out/'OXT_SOURCE_CHECK.json')}
        pp=prepdir/'preparation_result.json';write(pp,prepared)
        row=envelope.prepare_task(prepared,request['group'],old['config'],out/'expanded')
        result.update(row,chain=spec['chain'],source_preparation=pin(pp),selection_id=request['group']['selection_id'],
          raw_atom_map=files['source_audit'],chemical_state=old['chemical_state'],source_redox_experimentally_qualified=False)
        assert row['status']=='prepared_awaiting_group_check',row.get('reason',row['status'])
        assert row['mapping_eligible'] and row['paired_maps_equal'],'physical mapping gates failed'
        anchor=request['anchor'];result['local_state_identical_to_Q9']=row['state_key']==anchor['state_key'];result['reference_state_signature']=anchor['state_signature']
        assert result['local_state_identical_to_Q9'],'exact local state signature differs from Q9'
        actual=read(row['protein_identity']);ref=read(anchor['protein_identity'])
        positions={x['source_position'] for a in read(old['site_alignment_audit'])['mappings'] if a['chain']==spec['chain'] for x in a['mapping']}
        expected=[r for r in ref if r[0]=='A' and r[1] in positions]
        assert actual==expected and len(actual)==582,'exact mature protein identity differs'
        result['construct_compatibility']={'matched_residues':582,'mature_subset_exact':True,'full_protein_key_identical':row['protein_key']==anchor['protein_key'],
          'excluded_Q9_identity':[r for r in ref if r not in expected]}
        result['terminal_completion']=terminal_check(record,st,row,out,gemmi)
        result['spectator_admission']=helper.spectator_check(request,prepared,row,spectator,[record],out/'SPECTATOR_ADMISSION.json')
        result['omitted_contacts']=v2.omitted_contacts(old,spec,row,out,gemmi)
        result.update(status='prepared',preparation_profile=PROTOCOL,conditional_source_state=True)
    except Exception as e:
        result.update(status='unsupported',reason=str(e) or repr(e),exception_type=type(e).__name__);(out/'FAILURE.txt').write_text(traceback.format_exc())
    result.update(new_energy_calls=0,new_optimizations=0,new_protonation_calls=0)
    write(out/'ADMISSION.json',result);return result


def execute(path):
    assert os.environ.get('SLURM_JOB_ID'),'actual Slurm required'
    m=read(path);assert m['protocol']==PROTOCOL and pin(__file__)==m['implementation']['RECOVER.py']
    for p in m['implementation'].values():verify(p)
    for p in m['prior_implementation'].values():verify(p)
    assert pin(sys.executable)==m['python'];old=read(m['prior_manifest']);verify(m['prior_result'])
    for case in m['cases']:
        for p in case['files'].values():verify(p)
    helper=load('exaa_OXT_existing_spectator',verify(old['spectator_implementation']['PREPARE.py']))
    placement=helper.allocation(1);v2=load('exaa_OXT_prior_driver',verify(old['implementation']['PREPARE.py']))
    sys.path.insert(0,str(ROOT/'scripts'));out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);started=time.monotonic();write(out/'PLACEMENT.json',placement)
    rows=[run(m,old,c,out/('chain_'+c['spec']['chain']),v2,helper) for c in m['cases']]
    used=helper.imported_pins({'repository_script_pins':old['repository_script_pins']});u=resource.getrusage(resource.RUSAGE_SELF)
    final={'protocol':PROTOCOL,'manifest':pin(path),'terminal':True,'status':'complete','rows':rows,'denominator':2,
      'supported':sum(r['status']=='prepared' for r in rows),'config':old['config'],'prior_failed_result':m['prior_result'],
      'receipt':{'placement':placement,'wall_seconds':time.monotonic()-started,'CPU_seconds':u.ru_utime+u.ru_stime,'peak_RSS_KiB':u.ru_maxrss},
      'imported_source_pins':used,'new_energy_calls':0,'new_optimizations':0,'new_protonation_calls':0,'production_changed':False,
      'conditional_state':old['chemical_state'],'source_redox_experimentally_qualified':False}
    write(out/'FINAL_PREPARATION.json',final);print(json.dumps({'result':pin(out/'FINAL_PREPARATION.json'),'supported':final['supported'],'denominator':2}))


if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='command',required=True)
    a=s.add_parser('freeze');a.add_argument('--output',required=True)
    a=s.add_parser('execute');a.add_argument('--manifest',required=True)
    a=p.parse_args();freeze(a.output) if a.command=='freeze' else execute(a.manifest)
