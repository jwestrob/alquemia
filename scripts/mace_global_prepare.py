"""Declared whole-chain preparations for the five consumed accuracy structures."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_hybrid import check_atoms, write_xyz
from mace_hydrogen import repair as repair_hydrogen
from mace_gb import RADII

POLICY='whole_chain_A_frozen_cofactor_water_inventory_ff19sb_radial_H_terminal_OXT_v1'
CASES=('GGR_1GLG','ALPHA_1F6S','ALPHA_6IP9','PQQ_1H4I','PQQ_4MAE')
FF=Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/lib/python3.11/site-packages/openmm/app/data/amber19/protein.ff19SB.xml')
WATER_FF=FF.parent.parent/'tip3p.xml'
STANDARD=set('ALA ARG ASN ASP CYS GLN GLU GLY HID HIE HIP HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
LABELS={
 'GGR_1GLG':{'evidence_stratum':'direct_same_assay_affinity','direction':'Ca','group':'GGR'},
 'ALPHA_1F6S':{'evidence_stratum':'qualified_qualitative_strong_site_affinity','direction':'La','group':'bovine_alpha'},
 'ALPHA_6IP9':{'evidence_stratum':'qualified_qualitative_strong_site_affinity','direction':'La','group':'bovine_alpha'},
 'PQQ_1H4I':{'evidence_stratum':'canonical_PQQ_functional_metal_association','direction':'Ca','group':'MxaF'},
 'PQQ_4MAE':{'evidence_stratum':'canonical_PQQ_functional_metal_association','direction':'La','group':'XoxF'}}


def atom_id(atom):
    return f'{atom.residue.chain.id}/{atom.residue.id}/{atom.residue.insertionCode.strip()}/{atom.name}'


def terminal_position(carbon, alpha, oxygen, length=1.25):
    """Reflect C->O about C->CA; independent of coordinate axes and score."""
    c,ca,o=(np.array(x,dtype=float) for x in (carbon,alpha,oxygen))
    axis=ca-c; axis/=np.linalg.norm(axis)
    vector=o-c
    if not 1. < np.linalg.norm(vector) < 1.6 or abs(np.dot(axis,vector)/np.linalg.norm(vector))>.95:
        raise InvalidArtifact('unsupported terminal carboxyl geometry')
    reflected=2*np.dot(vector,axis)*axis-vector
    return c+length*reflected/np.linalg.norm(reflected)


def protein(row, allow_terminal_completion=False, check_peptide_connectivity=False,
            retain_source_acetyl=False):
    import openmm as mm
    from openmm import app,unit
    pdb=app.PDBFile(str(verify(row['source_structure'])))
    model=app.Modeller(pdb.topology,pdb.positions)
    supported = STANDARD | ({'ACE'} if retain_source_acetyl else set())
    inventory=[{'chain':r.chain.id,'resid':r.id,'icode':r.insertionCode.strip(),'resname':r.name,
                'atoms':len(list(r.atoms())),'selected_protein':r.chain.id=='A' and r.name in supported}
                for r in model.topology.residues()]
    model.delete([a for a in model.topology.atoms() if a.residue.chain.id!='A' or a.residue.name not in supported])
    residues=list(model.topology.residues())
    if not residues or len(list(model.topology.chains()))!=1:
        raise InvalidArtifact('exactly one standard protein chain A required')
    original=np.array(model.positions.value_in_unit(unit.angstrom))
    atoms=list(model.topology.atoms());original_ids=[atom_id(a) for a in atoms]
    if len(set(original_ids))!=len(atoms):raise InvalidArtifact('duplicate protein source identity')
    if retain_source_acetyl:
        acetyl = [r for r in residues if r.name == 'ACE']
        if acetyl != residues[:1] or len(residues) < 2:
            raise InvalidArtifact('exactly one actual N-terminal ACE residue required')
        if {a.name for a in acetyl[0].atoms()} != {'C','O','CH3','H1','H2','H3'}:
            raise InvalidArtifact('unsupported source acetyl atom inventory')
        links = [(a,b) for a,b in model.topology.bonds()
                 if a.residue != b.residue and acetyl[0] in (a.residue,b.residue)]
        if len(links) != 1 or {(a.residue.index,a.name) for a in links[0]} != {
                (residues[0].index,'C'), (residues[1].index,'N')}:
            raise InvalidArtifact('source acetyl must retain its real C--N peptide connection')
    if check_peptide_connectivity:
        for a,b in model.topology.bonds():
            if a.residue!=b.residue and {a.name,b.name}=={'C','N'}:
                distance=float(np.linalg.norm(original[a.index]-original[b.index]))
                if not 1.0<=distance<=1.8:
                    raise InvalidArtifact(f'unsupported peptide connection {atom_id(a)}--{atom_id(b)}: {distance:.6f} Angstrom; missing/malformed source geometry')
    pos=original.copy();additions=[]
    last=residues[-1];named={a.name:a for a in last.atoms()}
    if 'OXT' not in named:
        if (not allow_terminal_completion and row['case_id'] not in ('PQQ_1H4I','PQQ_4MAE')) or not {'C','CA','O'}<=named.keys():
            raise InvalidArtifact('unapproved missing terminal chemistry')
        carbon=named['C']
        neighbors=[b if a==carbon else a for a,b in model.topology.bonds() if carbon in (a,b)]
        if {a.name for a in neighbors}!={'CA','O'} or any(a.residue!=last for a in neighbors):
            raise InvalidArtifact('terminal completion would alter existing connectivity')
        newpos=terminal_position(pos[carbon.index],pos[named['CA'].index],pos[named['O'].index])
        oxygen=model.topology.addAtom('OXT',app.element.oxygen,last)
        model.topology.addBond(carbon,oxygen)
        pos=np.concatenate((pos,newpos[None,:]))
        model.positions=pos*.1*unit.nanometer
        additions.append({'id':atom_id(oxygen),'element':'O','xyz_A':newpos.tolist(),
                          'bond_to':atom_id(carbon),'length_A':1.25,
                          'recipe':'reflect_carbonyl_about_C_CA_v1',
                          'metal_distance_A':float(np.linalg.norm(newpos-row['selected_site']['xyz_A']))})
    ff=app.ForceField(str(FF))
    unmatched=ff.getUnmatchedResidues(model.topology)
    if unmatched:raise InvalidArtifact('unsupported protein topology: '+str([(r.name,r.id) for r in unmatched]))
    system=ff.createSystem(model.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    nb=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce))
    q=float(sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles())))
    if abs(q-round(q))>1e-6:raise InvalidArtifact('nonintegral standard-protein charge')
    atoms=list(model.topology.atoms());physical=[]
    for a in atoms:
        elem=a.element.symbol
        if elem not in RADII:raise InvalidArtifact('unsupported protein element')
        physical.append({'id':atom_id(a),'element':elem,'xyz_A':pos[a.index].tolist(),
                         'radius_A':RADII[elem],'source_index':a.index if a.index<len(original) else None,
                         'kind':'protein_source' if a.index<len(original) else 'terminal_completion',
                         'resname':a.residue.name})
    bondforce=next(f for f in system.getForces() if isinstance(f,mm.HarmonicBondForce))
    bonds=[]
    for i in range(bondforce.getNumBonds()):
        a,b,length,k=bondforce.getBondParameters(i)
        bonds.append({'atom_a_id':physical[a]['id'],'atom_b_id':physical[b]['id'],
                      'forcefield_equilibrium_A':length.value_in_unit(unit.angstrom)})
    repaired,moves=repair_hydrogen(physical,bonds)
    for old,new in zip(physical,repaired):
        if old['element']!='H' and old['xyz_A']!=new['xyz_A']:
            raise InvalidArtifact('protein heavy atom moved during H preparation')
    # No energy computation was used to prepare this structure.
    return repaired,round(q),{'source_residue_inventory':inventory,'terminal_additions':additions,
               'protein_H_moves':moves,'bonds':bonds,'original_protein_atoms':physical,
               'ff_charge_sum_e':q,'ff_atom_charges_used_for_scoring':False,
               'source_heavy_coordinates_preserved':True}


def cofactor_and_waters(row,manifest,expected_water_count=None):
    core=xyz(verify(row['endpoints']['La']['xyz'])); result=[];water_moves=[]
    cofactor_charge=0
    if row['case_id'].startswith('PQQ'):
        if manifest['pqq']['microstate_id']!='pqq_ox_3minus_v1' or manifest['pqq']['formal_charge']!=-3:
            raise InvalidArtifact('unsupported PQQ microstate')
        fragments=manifest['qm_fragments'];offset=1
        found=0
        for fragment in fragments:
            n=fragment['atom_count']
            if fragment['kind'] in ('pqq','fixed_core_pqq'):
                found+=1
                if n!=27:raise InvalidArtifact('wrong PQQ atom inventory')
                for i,(a,source) in enumerate(zip(core[offset:offset+n],fragment['atom_records'])):
                    if a[0]!=source['element'] or not np.allclose(a[1:],source['xyz_A'],atol=1e-6,rtol=0):
                        raise InvalidArtifact('PQQ coordinate/mapping mismatch')
                    result.append({'id':'pqq/'+source['name'],'element':a[0],'xyz_A':list(a[1:]),
                                   'radius_A':RADII[a[0]],'kind':'frozen_PQQ','source_record':source})
            offset+=n
        if found!=1:raise InvalidArtifact('missing/duplicate PQQ')
        cofactor_charge=-3
    water_inventory=manifest.get('explicit_water_inventory',[])
    if expected_water_count is None:
        expected_water_count={'GGR_1GLG':0,'ALPHA_1F6S':2,'ALPHA_6IP9':3,'PQQ_1H4I':0,'PQQ_4MAE':0}[row['case_id']]
    if len(water_inventory)!=expected_water_count:
        raise InvalidArtifact('frozen water count mismatch')
    for w in water_inventory:
        selected=[a for a in manifest['atom_graph']['source_to_qm'] if a['kind']=='source' and
                  all(a['source'][key]==w[key] for key in ('chain','resnum','insertion_code','resname'))]
        if sorted(a['source']['element'] for a in selected)!=['H','H','O']:
            raise InvalidArtifact('water has incomplete source atom mapping')
        oxygen=next(np.array(core[a['qm_index']][1:]) for a in selected if a['source']['element']=='O')
        for a in selected:
            source=a['source'];element=source['element'];coords=np.array(core[a['qm_index']][1:])
            name=f"{source['chain']}/{source['resnum']}/{source['insertion_code']}/{source['atom']}"
            if core[a['qm_index']][0]!=element:raise InvalidArtifact('water element mismatch')
            if element=='H':
                displacement=coords-oxygen;length=float(np.linalg.norm(displacement))
                if not .5<length<1.5:raise InvalidArtifact('unsupported water O-H geometry')
                target=oxygen+displacement*(.9572/length)
                water_moves.append({'id':name,'before_A':coords.tolist(),'after_A':target.tolist(),
                                    'old_length_A':length,'new_length_A':.9572})
                coords=target
            result.append({'id':name,'element':element,'xyz_A':coords.tolist(),'radius_A':RADII[element],
                           'kind':'retained_site_water','source_qm_index':a['qm_index'],'source_record':source})
    return result,cofactor_charge,water_inventory,water_moves


def prepare_case(row,output,*,allow_terminal_completion=False,expected_water_count=None,
                 evidence=None,policy_id=POLICY,check_peptide_connectivity=False):
    original_manifest=read_json(verify(row['preparation_manifest']))
    p,q,details=protein(row,allow_terminal_completion=allow_terminal_completion,
                      check_peptide_connectivity=check_peptide_connectivity)
    extra,cofactor_q,waters,water_moves=cofactor_and_waters(row,original_manifest,expected_water_count)
    metal=xyz(verify(row['endpoints']['La']['xyz']))[0]
    if metal[0]!='La' or not np.allclose(metal[1:],row['selected_site']['xyz_A'],atol=1e-6,rtol=0):
        raise InvalidArtifact('selected metal coordinate mismatch')
    physical=[*p,*extra,{'id':'metal','element':'M','xyz_A':list(metal[1:]),'radius_A':1.8,'kind':'selected_metal'}]
    if len({a['id'] for a in physical})!=len(physical):raise InvalidArtifact('duplicate global physical atom')
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    endpoints={}
    for element in ('La','Ca'):
        rows=[(element if a['element']=='M' else a['element'],*a['xyz_A']) for a in physical]
        charge=q+cofactor_q+(3 if element=='La' else 2)
        state=check_atoms(rows,charge);path=out/(element+'.xyz');write_xyz(path,rows)
        endpoints[element]={'xyz':record(path),'charge':charge,'spin_multiplicity':1,'state':state}
    result={'status':'prepared','case_id':row['case_id'],'policy_id':policy_id,'source':row['source_structure'],
            'source_preparation':row['preparation_manifest'],'source_audit_row':row,
            'physical_atoms':physical,'assembly':'deposited_chain_A_with_declared_cofactor_and_site_waters',
            'explicit_waters':waters,'water_H_moves':water_moves,'protein_charge_e':q,'cofactor_charge_e':cofactor_q,
            'forcefield':record(FF),'water_forcefield':record(WATER_FF),'preparation_details':details,
            'endpoints':endpoints,'microstate':'archived_explicit_protonation_and_disulfides; PQQ3minus_where_present; neutral_waters',
            'evidence':LABELS[row['case_id']] if evidence is None else evidence,'evidence_use':'consumed_method_development',
            'reference':None,'baseline_coordinates_modified':False}
    write_new(out/'preparation.json',result)
    return record(out/'preparation.json')


def prepare(audit,plan,output):
    start=time.monotonic();source=read_json(audit);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    if tuple(row['case_id'] for row in source['rows'])!=CASES:raise InvalidArtifact('source panel inventory changed')
    rows=[]
    for row in source['rows']:
        try:rows.append({'case_id':row['case_id'],'status':'prepared','preparation':prepare_case(row,out/row['case_id'])})
        except Exception as exc:rows.append({'case_id':row['case_id'],'status':'unsupported','reason':str(exc),'exception_type':type(exc).__name__})
    result={'status':'complete' if all(r['status']=='prepared' for r in rows) else 'incomplete',
            'source_audit':record(audit),'agreement':record(plan),'policy_id':POLICY,'rows':rows,
            'implementation':{name:record(Path(__file__).with_name(name)) for name in
                              ('mace_global_prepare.py','mace_hydrogen.py','mace_gb.py','mace_hybrid.py','affordable_common.py')},
            'wall_seconds':time.monotonic()-start,'process_CPU_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,
            'new_energy_evaluations':0}
    write_new(out/'preparation_manifest.json',result)
    preserved=out/'implementation';preserved.mkdir()
    references={}
    for name,ref in result['implementation'].items():
        shutil.copyfile(verify(ref),preserved/name)
        references[name]={'recorded_source':ref,'preserved_copy':record(preserved/name)}
    # Needed by the hydrogen helper when replaying preparation from the snapshot.
    shutil.copyfile(Path(__file__).with_name('mace_response_trace.py'),preserved/'mace_response_trace.py')
    write_new(out/'implementation_preservation.json',references)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('audit','plan','output'):p.add_argument('--'+key,required=True)
    a=p.parse_args();r=prepare(a.audit,a.plan,a.output)
    print(json.dumps({'status':r['status'],'rows':r['rows']},indent=2))
