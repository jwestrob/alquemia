#!/usr/bin/env python3
"""Two pinned crystal chains; existing preparation kernels, Slurm-only science."""
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
import types

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
HERE = Path(__file__).resolve().parent
PROTOCOL = 'ExaA_1FLG_AB_exact_mature_map_conditional_oxidized_PQQ_frozen_union_v1'


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    assert pin(p['path']) == p, ('pin changed', p['path'])
    return Path(p['path'])


def read(p):
    return json.loads((verify(p) if isinstance(p, dict) else Path(p)).read_text())


def write(p, value):
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); sys.modules[name] = m
    spec.loader.exec_module(m); return m


def freeze(output):
    """Metadata binding only: no scientific import or coordinate analysis."""
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    oldpath = ROOT/'workspaces/ca_coverage_recovery_20261012/spectator_preparation_q9_v2/manifest.json'
    old = read(oldpath)
    auditpath = ROOT/'workspaces/ca_coverage_recovery_20261012/refold/run_v3_h200/prepared/1FLG_AUDIT.json'
    audit = read(auditpath); verify(audit['reference'])
    completed = oldpath.parent/'execution/FINAL_PREPARATION.json'
    assert read(completed)['supported'] == 5
    m = {'protocol': PROTOCOL, 'implementation': {n: pin(HERE/n) for n in ('PREPARE.py','PLAN.md','SOURCES.md','RUN.sbatch')},
         'python': old['python'], 'spectator_implementation': old['implementation'],
         'repository_script_pins': old['repository_script_pins'], 'config': old['config'],
         'reference_preparation_manifest': pin(oldpath), 'reference_preparation': pin(completed),
         'reference_request': old['requests'][0], 'source': audit['reference'], 'site_alignment_audit': pin(auditpath),
         'cofactor_metrics': pin(ROOT/'workspaces/pede_solvent_response_20261013/cofactor_review_v2/execution/RESULT.json'),
         'chains': [{'chain':'A','PQQ':701,'active_Ca':751,'spectator_Ca':761},
                    {'chain':'B','PQQ':702,'active_Ca':752,'spectator_Ca':762}],
         'output': str(out/'execution'), 'expected_denominator': 2,
         'chemical_state': 'conditional_existing_oxidized_PQQ3minus_not_experimentally_verified',
         'new_energies': 0, 'new_heavy_repairs': False, 'new_optimizations': False}
    for p in m['implementation'].values(): verify(p)
    for p in m['spectator_implementation'].values(): verify(p)
    verify(m['reference_request']['group']['union'])
    write(out/'manifest.json',m)
    print(json.dumps({'manifest':pin(out/'manifest.json'),'status':'prepared_only_Slurm_UNRUN'}))


def meta(chain, residue, atom):
    return {'chain':chain.name,'resnum':residue.seqid.num,'icode':residue.seqid.icode.strip(),
            'resname':residue.name,'atom':atom.name.strip(),'element':atom.element.name,
            'altloc':atom.altloc.strip('\x00 '),'occupancy':float(atom.occ),'B_iso':float(atom.b_iso),
            'xyz_A':list(atom.pos)}


def identity(x):
    return tuple(x[k] for k in ('chain','resnum','icode','resname','atom','element','altloc'))


def selected_source(m, spec, out, crystal, gemmi, detect):
    raw = verify(m['source']); st = gemmi.read_structure(str(raw))
    assert len(st) == 1 and [c.name for c in st[0]] == ['A','B']
    model = crystal.exact_model(st,'1'); chain = crystal.exact_chain(model,spec['chain'])
    oldmap = next(x for x in read(m['site_alignment_audit'])['mappings'] if x['chain'] == spec['chain'])
    assert oldmap['identity_percent'] == 100.0 and all(x['identity'] for x in oldmap['mapping'])
    lookup = {(x['crystal_resnum'],x['crystal_icode'],x['crystal_resname']):x['source_position'] for x in oldmap['mapping']}
    assert len(lookup) == 582
    block = gemmi.cif.read_file(str(raw)).sole_block()
    categories = {name:block.get_mmcif_category(name) for name in ('_struct_conn.','_chem_comp.','_entity.','_pdbx_unobs_or_zero_occ_residues.')}
    conn = categories['_struct_conn.']; external=[]
    for i, kind in enumerate(conn.get('conn_type_id',[])):
        row = {k:v[i] for k,v in conn.items()}
        if 'PQQ' in (row.get('ptnr1_auth_comp_id'),row.get('ptnr2_auth_comp_id')) and kind != 'metalc': external.append(row)
    audit = {'source':m['source'],'chain':spec['chain'],'declared_sites':spec,'deposited_categories':categories,
             'PQQ_external_nonmetal_connections':external,'source_redox':'unassigned_conditional_oxidized_model_only',
             'conformer_policy':crystal.ALTLOC_POLICY_ID,'excluded_residues':[],'alternate_conformers':[],
             'retained_atom_map':[],'source_nonstandard_inventory':[], 'local_preparation_status':'not_started'}
    write(out/'SOURCE_METADATA.json',audit)
    assert not external, 'deposited PQQ external connection unsupported'
    selected = gemmi.Structure(); selected.name='1FLG_'+spec['chain']+'_mapped_Q9'
    sm = gemmi.Model(1); dest = {c:gemmi.Chain(c) for c in ('A','B','C','D')}
    for residue in chain:
        key=(residue.seqid.num,residue.seqid.icode.strip(),residue.name)
        kind='protein' if key in lookup else 'PQQ' if residue.seqid.num==spec['PQQ'] and residue.name=='PQQ' else 'active' if residue.seqid.num==spec['active_Ca'] and residue.name=='CA' else 'spectator' if residue.seqid.num==spec['spectator_Ca'] and residue.name=='CA' else None
        if residue.name not in crystal.STANDARD_AA:
            audit['source_nonstandard_inventory'].append({'chain':chain.name,'resnum':residue.seqid.num,'resname':residue.name,'atoms':[meta(chain,residue,a) for a in residue]})
        if not kind:
            assert residue.name in crystal.WATER_NAMES, ('unmapped protein or unsupported heterogen',key)
            audit['excluded_residues'].append({'chain':chain.name,'resnum':residue.seqid.num,'resname':residue.name,'reason':'frozen_dry_water_exclusion','atom_count':len(residue)})
            continue
        if kind!='protein': assert not crystal.residue_has_altloc(residue), 'PQQ/metal alternate conformer unsupported'
        if crystal.residue_has_altloc(residue): audit['alternate_conformers'].append({'residue':list(key),**crystal.altloc_decision(residue)})
        nr=crystal.clone_residue_selected(residue,normalize_to_la=(kind=='active'))
        target={'protein':'A','active':'B','PQQ':'C','spectator':'D'}[kind]
        nr.seqid=gemmi.SeqId(lookup[key] if kind=='protein' else 1,' ')
        nr.subchain=target
        source_atoms=crystal.selected_residue_atoms(residue)
        assert len(source_atoms)==len(nr)
        for a,b in zip(source_atoms,nr):
            assert a.element.name not in ('H','D'), 'raw source H unsupported by fixed protonation policy'
            assert list(a.pos)==list(b.pos) and a.element.name==b.element.name and a.occ==b.occ
            audit['retained_atom_map'].append({'source':meta(chain,residue,a),'mapped':meta(dest[target],nr,b)})
        dest[target].add_residue(nr)
    for c in dest.values(): sm.add_chain(c)
    selected.add_model(sm); selected.setup_entities()
    path=out/'selected_dry_mapped_source.cif'; selected.make_mmcif_document().write_file(str(path))
    reread=gemmi.read_structure(str(path)); expected={identity(x['mapped']):x['mapped'] for x in audit['retained_atom_map']}
    actual={identity(meta(c,r,a)):meta(c,r,a) for c in reread[0] for r in c for a in r}
    assert expected==actual, 'source mapping/coordinates/occupancy changed in serialization'
    pqq=detect.inspect_pqq_candidates(reread[0]); assert len(pqq)==1
    audit['PQQ_schema']=pqq[0].prepared.metadata()
    assert len(pqq[0].residue)==24 and all(a.element.name not in ('H','D') for a in pqq[0].residue)
    audit['selected_source']=pin(path); audit['source_heavy_coordinates_exactly_preserved']=True
    audit['full_source_dimer_retained_as_evidence']=m['source']
    audit['source_sequence_map']=m['site_alignment_audit']
    write(out/'SOURCE_AUDIT.json',audit)
    return path,audit,lookup


def omitted_contacts(m, spec, row, out, gemmi):
    """Descriptive saved-source omissions only; same existing4.3A radius."""
    import numpy as np
    from affordable_common import xyz
    raw=gemmi.read_structure(str(verify(m['source'])))[0]
    points=xyz(verify(row['representations']['context']['endpoints']['Ca']['xyz']))
    heavy=[(i,a) for i,a in enumerate(points) if a[0] not in ('H','D')]
    pos=np.asarray([x[1:] for _,x in heavy]); rows=[]
    for c in raw:
        for r in c:
            if c.name==spec['chain'] and r.name not in ('HOH','WAT','DOD'): continue
            for a in r:
                if a.element.name in ('H','D'):continue
                distances=np.linalg.norm(pos-np.asarray(list(a.pos)),axis=1); j=int(np.argmin(distances))
                if float(distances[j])<=4.3: rows.append({'source':meta(c,r,a),'distance_to_retained_heavy_A':float(distances[j]),'nearest_context_index':heavy[j][0],'reason':'source_water' if r.name in ('HOH','WAT','DOD') else 'other_crystal_chain'})
    result={'description':'omitted water/other-chain source contacts; no new admission rule or energy claim','radius_A':4.3,'contacts':rows}
    write(out/'OMITTED_CONTACTS.json',result);return pin(out/'OMITTED_CONTACTS.json')


def run_chain(m,spec,out,helper):
    out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    cid='1FLG_chain_'+spec['chain']+'__mapped_Q9_conditional_oxidized_v1'
    result={'case_id':cid,'chain':spec['chain'],'status':'unsupported','mapping_eligible':False,
            'chemical_state':m['chemical_state'],'source':m['source'],'source_redox_experimentally_qualified':False}
    try:
        import gemmi,detect_pqq,pqq_fast_prepare as fast,motion_envelope as envelope
        fixed,legacy,crystal,normalizer,cp=fast.helpers(m['config'])
        path,audit,lookup=selected_source(m,spec,out,crystal,gemmi,detect_pqq)
        request=copy.deepcopy(m['reference_request']);source=request['source']
        for k in list(source):
            if k.startswith('archived_') or k in ('source_geometry_audit','original_source_id'):source.pop(k)
        active={'chain':'B','resnum':1,'icode':'','resname':'LA','atom':'CA','element':'Ca'}
        spectator={'chain':'D','resnum':1,'icode':'','resname':'CA','atom':'CA','element':'Ca'}
        source.update(case_id=cid,source_structure=pin(path),metal=active,raw_source_metal=active,
            original_source_structure=m['source'],source_kind='crystal_conditional_oxidized_PQQ',
            source_conditioning_metal=None,source_occupancy_state='deposited_native_2Ca',
            assembly={'kind':'explicit_selected_crystal_chain_with_mature_sequence_alias_map','raw_chain':spec['chain'],'raw_source':m['source'],'atom_map':pin(out/'SOURCE_AUDIT.json')})
        request.update(source=source,original_source_id=cid,audit_row={'spectator_Ca':spectator})
        wrapper=load('exaa_spectator_protonator',verify(m['spectator_implementation']['PROTONATE.py']))
        narrow=types.SimpleNamespace(protonate_standard_only=lambda src,output,*,pins_path,ph:wrapper.protonate(src,output,pins_path=pins_path,ph=ph,legacy=legacy,spectator_count=1))
        scope=dict(fast.prepare_case.__globals__);scope['helpers']=lambda config:(fixed,narrow,crystal,normalizer,cp)
        prepare=types.FunctionType(fast.prepare_case.__code__,scope,'prepare_1FLG_with_existing_spectator_wrapper')
        prepared=prepare(source,m['config'],out/'source_preparation')
        assert prepared['status']=='prepared'
        core=read(prepared['core']['parent'])
        preservation=[helper.conserved(path,verify(core[k]),0.001,True) for k in ('normalized_source_structure','source_structure')]
        row=envelope.prepare_task(prepared,request['group'],m['config'],out/'expanded')
        result.update(row,chain=spec['chain'],source_preparation=pin(out/'source_preparation/preparation_result.json'),
            selection_id=request['group']['selection_id'],source_preservation=preservation,
            raw_atom_map=pin(out/'SOURCE_AUDIT.json'),chemical_state=m['chemical_state'],source_redox_experimentally_qualified=False)
        assert row['status']=='prepared_awaiting_group_check',row.get('reason',row['status'])
        assert row['mapping_eligible'] and row['paired_maps_equal']
        anchor=request['anchor']; result['local_state_identical_to_Q9']=row['state_key']==anchor['state_key']
        result['reference_state_signature']=anchor['state_signature']
        assert result['local_state_identical_to_Q9'],'exact local state signature differs from frozen Q9'
        actual=read(row['protein_identity']);ref=read(anchor['protein_identity'])
        positions=set(lookup.values());expected=[r for r in ref if r[0]=='A' and r[1] in positions]
        assert actual==expected and len(actual)==582,'protein identity differs from exact mature aligned subset'
        result['construct_compatibility']={'matched_residues':582,'mapped_Q9_first':min(positions),'mapped_Q9_last':max(positions),
            'full_protein_key_identical':row['protein_key']==anchor['protein_key'],'mature_subset_exact':True,
            'excluded_Q9_identity':[r for r in ref if r not in expected]}
        result['spectator_admission']=helper.spectator_check(request,prepared,row,spectator,preservation,out/'SPECTATOR_ADMISSION.json')
        result['omitted_contacts']=omitted_contacts(m,spec,row,out,gemmi)
        result.update(status='prepared',preparation_profile=PROTOCOL,conditional_source_state=True)
    except Exception as exc:
        result.update(status='unsupported',reason=str(exc),exception_type=type(exc).__name__)
        (out/'FAILURE.txt').write_text(traceback.format_exc())
    result.update(wall_seconds=time.monotonic()-start,new_energy_calls=0,new_optimizations=0)
    write(out/'ADMISSION.json',result);return result


def execute(path):
    assert os.environ.get('SLURM_JOB_ID'),'scientific preparation requires actual Slurm'
    m=read(path);assert m['protocol']==PROTOCOL and pin(__file__)==m['implementation']['PREPARE.py']
    assert pin(sys.executable)==m['python']
    for group in ('implementation','spectator_implementation'):
        for p in m[group].values():verify(p)
    for k in ('source','site_alignment_audit','reference_preparation','reference_preparation_manifest','cofactor_metrics'):verify(m[k])
    helper=load('exaa_existing_spectator_helper',verify(m['spectator_implementation']['PREPARE.py']))
    placement=helper.allocation(1)
    # Existing helper verifies every actually imported repository source against frozen pins.
    imported_manifest={'repository_script_pins':m['repository_script_pins']}
    sys.path.insert(0,str(ROOT/'scripts'))
    out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    write(out/'PLACEMENT.json',placement)
    rows=[run_chain(m,s,out/('chain_'+s['chain']),helper) for s in m['chains']]
    used=helper.imported_pins(imported_manifest)
    usage=resource.getrusage(resource.RUSAGE_SELF)
    final={'protocol':PROTOCOL,'manifest':pin(path),'terminal':True,'status':'complete','rows':rows,
           'denominator':2,'supported':sum(r['status']=='prepared' for r in rows),'config':m['config'],
           'receipt':{'placement':placement,'wall_seconds':time.monotonic()-start,'CPU_seconds':usage.ru_utime+usage.ru_stime,'peak_RSS_KiB':usage.ru_maxrss},
           'imported_source_pins':used,'new_energy_calls':0,'new_optimizations':0,'production_changed':False,
           'conditional_state':m['chemical_state'],'source_redox_experimentally_qualified':False}
    write(out/'FINAL_PREPARATION.json',final)
    print(json.dumps({'result':pin(out/'FINAL_PREPARATION.json'),'supported':final['supported'],'denominator':2}))


if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='command',required=True)
    a=s.add_parser('freeze');a.add_argument('--output',required=True)
    a=s.add_parser('execute');a.add_argument('--manifest',required=True)
    a=p.parse_args();freeze(a.output) if a.command=='freeze' else execute(a.manifest)
