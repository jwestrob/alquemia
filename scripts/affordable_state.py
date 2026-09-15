"""Traceable whole-chain environment for the three pinned PQQ development cores."""
from __future__ import annotations

import re
import argparse
import time
import resource
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import openmm
from openmm import app, unit

from affordable_common import InvalidArtifact, read_json, record, verify, xyz, write_new
from affordable_environment import physical_boundary_key

FF=Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/lib/python3.11/site-packages/openmm/app/data/amber19/protein.ff19SB.xml')
RADII={'H':1.2,'C':1.7,'N':1.55,'O':1.52,'S':1.8,'La':1.8,'Ca':1.8}


def validate_skeleton_pair(la,ca):
    for key in ('source','assembly','microstate','explicit_waters','environment_atoms','expected_environment_charge_e'):
        if la[key]!=ca[key]: raise InvalidArtifact(f'paired environment differs: {key}')
    if physical_boundary_key(la['physical_atoms'])!=physical_boundary_key(ca['physical_atoms']):
        raise InvalidArtifact('paired physical cavity differs')
    if la['core_total_charge_e']-ca['core_total_charge_e']!=1:
        raise InvalidArtifact('paired charge difference is not one')
    left,right=la['core_atoms'],ca['core_atoms']
    if len(left)!=len(right): raise InvalidArtifact('paired core atom count differs')
    for a,b in zip(left,right):
        if any(a[k]!=b[k] for k in ('id','xyz_A','radius_A','kind')) or (a['id']!='metal' and a['element']!=b['element']):
            raise InvalidArtifact('paired physical core mapping differs')
    return {'status':'pass','physical_boundary_hash':physical_boundary_key(la['physical_atoms'])}


def source_protein(source, *, selected_chain='A'):
    p=app.PDBFile(str(source)); modeller=app.Modeller(p.topology,p.positions)
    allowed=set('ALA ARG ASN ASP CYS GLN GLU GLY HID HIE HIP HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
    # The declared physical system is chain A. Record other chains explicitly
    # instead of confusing a deposited assembly with the selected calculation.
    excluded_chains=[c.id for c in modeller.topology.chains() if c.id!=selected_chain]
    modeller.delete([a for a in modeller.topology.atoms() if a.residue.chain.id!=selected_chain])
    nonstandard=[r for r in modeller.topology.residues() if r.name not in allowed]
    if any(r.name not in ('PQQ','LA','CA','CE') for r in nonstandard):
        raise InvalidArtifact('unsupported full-source cofactor/water; no silent omission')
    modeller.delete(nonstandard)
    system=app.ForceField(str(FF)).createSystem(modeller.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    nb=next(f for f in system.getForces() if isinstance(f,openmm.NonbondedForce))
    atoms=list(modeller.topology.atoms()); pos=np.array(modeller.positions.value_in_unit(unit.angstrom))
    charges=np.array([nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(len(atoms))])
    physical=[]
    for a in atoms:
        if a.element.symbol not in RADII: raise InvalidArtifact('unsupported element radius')
        physical.append({'id':f'{a.residue.chain.id}/{a.residue.id}/{a.residue.insertionCode}/{a.name}',
                         'element':a.element.symbol,'xyz_A':pos[a.index].tolist(),'radius_A':RADII[a.element.symbol],
                         'charge_e':0.0,'source_index':a.index})
    return atoms,pos,charges,physical,excluded_chains


def prepare_skeleton(case, task, output):
    m=read_json(verify(case['source_manifest']))
    if 'source_structure' in m:
        source=verify(m['source_structure'])
    else:
        pm=read_json(verify(m['protonation_manifest'])); source=verify(pm['output'])
    atoms,pos,ffcharges,physical,excluded_chains=source_protein(source)
    tree=cKDTree(pos); core_xyz=xyz(verify(task['xyz']))
    core=[{'id':'metal','element':core_xyz[0][0],'xyz_A':list(core_xyz[0][1:]),'radius_A':1.8,'kind':'metal'}]
    # Metal cavity ID and radius are common across endpoints.
    physical.append(dict(core[0],element='M',charge_e=0.0))
    removed=set(); charges=ffcharges.copy(); ledgers=[]; offset=1
    for f in m['qm_fragments']:
        n=f['atom_count']; fragment_rows=core_xyz[offset:offset+n]
        if len(fragment_rows)!=n: raise InvalidArtifact('fragment ordering/count mismatch')
        if f['kind'] in ('protein_sidechain','fixed_core_protein_sidechain','fixed_core_cationic_sidechain'):
            match=re.fullmatch(r'([^:]+):([A-Z]+)(-?\d+)([A-Za-z]?)',f['id'])
            if not match: raise InvalidArtifact('unsupported residue identifier')
            chain,resname,resid,icode=match.groups()
            residues={a.residue for a in atoms if a.residue.chain.id==chain and a.residue.id==resid and a.residue.insertionCode.strip()==icode}
            if len(residues)!=1: raise InvalidArtifact('ambiguous/missing fragment source residue')
            residue=next(iter(residues)); byname={a.name:a for a in residue.atoms()}
            source_indices=set()
            for i,row in enumerate(fragment_rows):
                idx=offset+i
                if i==n-1:
                    if row[0]!='H': raise InvalidArtifact('missing final sidechain cap')
                    ca,cb=byname['CA'].index,byname['CB'].index
                    expected=pos[cb]+1.09*(pos[ca]-pos[cb])/np.linalg.norm(pos[ca]-pos[cb])
                    if np.linalg.norm(expected-np.array(row[1:]))>.002: raise InvalidArtifact('cap geometry differs from source bond')
                    core.append({'id':f'cap/{f["id"]}','element':'H','kind':'cap','xyz_A':list(row[1:]),'radius_A':1.2})
                else:
                    candidates=[j for j in tree.query_ball_point(row[1:],.002) if atoms[j].element.symbol==row[0] and atoms[j].residue==residue]
                    if len(candidates)!=1: raise InvalidArtifact('source-to-QM mapping is absent/ambiguous')
                    j=candidates[0];source_indices.add(j)
                    core.append(dict(physical[j],xyz_A=list(row[1:]),kind='protein_source'))
                    # The source is serialized to 0.001 A, QM to 1e-6 A. For matched
                    # grid cancellation the physical source coordinates are used
                    # only if equal at source precision; never move the QM charge.
                    physical[j]['xyz_A']=list(row[1:])
            local=source_indices|{byname['CA'].index}
            if removed & local: raise InvalidArtifact('overlapping covalent boundary')
            allres=[a.index for a in residue.atoms()]
            recipients=[byname[x].index for x in ('N','C')]
            wanted=float(ffcharges[allres].sum())-f['formal_charge']
            retained=[i for i in allres if i not in local]
            change=wanted-float(ffcharges[retained].sum())
            for i in recipients: charges[i]+=change/2
            removed.update(local)
            ledgers.append({'fragment_id':f['id'],'formal_charge':f['formal_charge'],
                            'removed_source_indices':sorted(local),'source_qm_indices':sorted(source_indices),
                            'recipients':recipients,'each_increment_e':change/2,
                            'original_residue_charge_e':float(ffcharges[allres].sum())})
        elif f['kind'] in ('pqq','fixed_core_pqq'):
            for i,row in enumerate(fragment_rows):
                a={'id':f'pqq/{i}','element':row[0],'kind':'pqq','xyz_A':list(row[1:]),'radius_A':RADII[row[0]]}
                core.append(a);physical.append(dict(a,charge_e=0.0))
        else:
            raise InvalidArtifact(f'unsupported environment fragment {f["kind"]}')
        offset+=n
    if offset!=len(core_xyz) or len(core)!=len(core_xyz): raise InvalidArtifact('unaccounted core atoms/waters')
    environment=[dict(physical[i],charge_e=float(charges[i])) for i in range(len(atoms)) if i not in removed]
    expected=float(ffcharges.sum())-sum(l['formal_charge'] for l in ledgers)
    if abs(sum(a['charge_e'] for a in environment)-expected)>1e-6: raise InvalidArtifact('local charge closure failed')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    boundary={'source':record(source),'forcefield':record(FF),'ff_chain_charge_e':float(ffcharges.sum()),
              'ledgers':ledgers,'expected_environment_charge_e':expected,'source_protein_atom_count':len(atoms),
              'nonstandard_state':'full PQQ and selected metal in QM; no other heterogens accepted',
              'microstate':'existing source protonation, PQQ(3-), no water',
              'mapping_max_tolerance_A':.002,'explicitly_excluded_other_chains':excluded_chains,
              'physical_assembly':'deposited_catalytic_chain_A'}
    write_new(output/'boundary_mapping.json',boundary)
    return {'source':record(source),'assembly':'deposited_catalytic_chain_A','microstate':'frozen_PQQ_3minus_standard_source_protonation',
            'explicit_waters':[],'boundary_mapping':record(output/'boundary_mapping.json'),
            'physical_atoms':physical,'core_atoms':core,'environment_atoms':environment,
            'core_total_charge_e':task['charge'],'expected_environment_charge_e':expected,
            'charge_quality':{'status':'unavailable'},'source_endpoint_xyz':task['xyz'],'case':case['case'],'metal':task['metal']}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pilot-manifest',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();m=read_json(a.pilot_manifest)
    a.output.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();rows=[]
    for task in m['tasks']:
        case=next(c for c in m['cases'] if c['case']==task['case'])
        before=time.monotonic();row={'task_id':task['task_id']}
        try:
            destination=a.output/task['task_id'];s=prepare_skeleton(case,task,destination)
            write_new(destination/'skeleton.json',s)
            row.update(status='prepared_unvalidated_charges',skeleton=record(destination/'skeleton.json'))
        except Exception as exc:
            row.update(status='environment_preparation_unsupported',reason=str(exc))
        row['wall_seconds']=time.monotonic()-before;rows.append(row)
    usage=resource.getrusage(resource.RUSAGE_SELF)
    write_new(a.output/'preparation_receipt.json',{'pilot':record(a.pilot_manifest),'implementation':record(__file__),
        'rows':rows,'wall_seconds':time.monotonic()-start,'process_user_seconds':usage.ru_utime,
        'process_system_seconds':usage.ru_stime,'process_peak_RSS_KiB':usage.ru_maxrss,
        'scope':'serial preparation only; process CPU/RSS include Python imports; not allocated Slurm core-seconds'})
    print(f"{sum(r['status']=='prepared_unvalidated_charges' for r in rows)}/{len(rows)} environments prepared; failures retained")


if __name__=='__main__': main()
