#!/usr/bin/env python3
"""LBT source-derived local/MM partitions. Numerical preparation requires Slurm.

prepare --config CONFIG.json --output FRESH
No Context, molecular energy, optimizer, SCC, DFT or model load occurs here.
Public update_geometry/scatter functions use full physical coordinate order.
"""
from __future__ import annotations
import argparse
from collections import Counter
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
import traceback
import xml.etree.ElementTree as ET

PROTOCOL = 'lbt_source_segment_local_MM_partition_v1'
ADMISSION_REVISION = 'common_MM_inactive_target_sigma_v2'
SOURCES = ('7CCO', '7CCN')
METALS = ('La', 'Dy')
SEED = (3, 14)
DONORS = ('A/3//OD1','A/3//OD2','A/5//OD1','A/7//OD1','A/7//OD2',
          'A/9//O','A/11//OE1','A/11//OE2','A/14//OE1','A/14//OE2')
HB_SEEDS = DONORS + ('A/5//ND2',)
HB_RULE = {'D_A_A':3.5, 'H_A_A':2.6, 'D_H_A_min_degrees':120.,
           'exclude_native_bond_distance_le':3, 'source_union':True}
TOL = {'charge_e':1e-9, 'coordinate_A':1e-12, 'map_FD_steps_A':[1e-5,5e-6],
       'map_FD_max_error':1e-6, 'map_translation_max_error':1e-12}
NUMBERS = {'H':1,'C':6,'N':7,'O':8,'La':57,'Dy':66}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve(strict=True)
    return {'path':str(p), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    p = Path(record['path']).resolve(strict=True)
    require(pin(p)['sha256'] == record['sha256'], 'Changed input: '+str(p))
    return p


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')
    return pin(path)


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def initialize(cfg):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(), 'Slurm required')
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Adapter differs')
    require(pin(sys.executable)['sha256'] == cfg['python']['sha256'], 'Python executable differs')
    verify(cfg['plan']); verify(cfg['region_policy']); verify(cfg['metal_parameters'])
    recovery=read(cfg['recovery_evidence'])
    failed=read(recovery['failed_result'])
    verify(recovery['failed_implementation']);verify(recovery['failed_config'])
    require(failed['slurm_job_id']=='12950' and failed['new_molecular_calls']==0
            and len(failed['rows'])==4 and all(r['status']=='unavailable' and
                r['reason']=='La/Dy MM common component definitions differ' for r in failed['rows']),
            'Recovery must preserve the four original comparison failures')
    for fixture in recovery['fixtures']:
        for record in fixture['systems'].values(): verify(record)
    for record in cfg['policy_references'].values():
        verify(record)
    for record in cfg['helpers'].values():
        verify(record)
    helper_dir = verify(cfg['helpers']['affordable_common.py']).parent
    require(all(verify(p).parent == helper_dir for p in cfg['helpers'].values()), 'Helper directory differs')
    sys.path.insert(0, str(helper_dir))
    global np, mm, unit, cap_jacobians, build_system, graph_distances, lj_from_amber
    import numpy as np
    import openmm as mm
    from openmm import unit
    from affordable_response import cap_jacobians
    from metal_environment_component_checks import build_system
    from metal_environment_mechanics import graph_distances, lj_from_amber
    require(mm.__version__ == cfg['openmm_version'], 'OpenMM version differs')


def native_source(source_id, preparation_pin):
    pr = read(preparation_pin)
    require(pr['status']=='prepared' and pr['source_id']==source_id and pr['atoms']==256
            and pr['protein_atoms']==255 and pr['protein_charge']==-5 and pr['charge']==-2,
            'Unexpected parent preparation/state')
    atoms = read(pr['native_peptide_atoms']); physical = read(pr['physical_atoms'])
    require(len(atoms)==255 and len(physical)==256 and physical[-1]['id']=='metal'
            and [a['id'] for a in atoms]==[a['id'] for a in physical[:-1]], 'Physical/native order differs')
    require(all(a['xyz_A']==b['xyz_A'] and a['element']==b['element'] for a,b in zip(atoms,physical)),
            'Physical/native coordinates differ')
    require(len({a['id'] for a in physical})==256 and all(a['element'] in NUMBERS for a in atoms),
            'Unsupported identity/element inventory')
    parent = mm.XmlSerializer.deserialize(verify(pr['native_peptide_system']).read_text())
    require(parent.getNumParticles()==255 and parent.getNumConstraints()==0, 'Native system inventory differs')
    force_names=[type(f).__name__ for f in parent.getForces()]
    require(len(force_names)==len(set(force_names)), 'Ledger indices require one native force of each type')
    nb = next(f for f in parent.getForces() if isinstance(f,mm.NonbondedForce))
    require(nb.getNonbondedMethod()==mm.NonbondedForce.NoCutoff, 'Finite native NoCutoff system required')
    parameters=[]
    for i,a in enumerate(atoms):
        q,s,e=nb.getParticleParameters(i)
        row=dict(a, native_charge_e=q.value_in_unit(unit.elementary_charge),
                 sigma_nm=s.value_in_unit(unit.nanometer), epsilon_kJ_mol=e.value_in_unit(unit.kilojoule_per_mole))
        require(abs(row['native_charge_e']-a['charge_e'])<=TOL['charge_e'], 'Native charge replay differs')
        parameters.append(row)
    require(abs(math.fsum(a['native_charge_e'] for a in parameters)+5)<=TOL['charge_e'], 'Peptide charge differs')
    ids={a['id']:i for i,a in enumerate(physical)}; adjacency={i:set() for i in range(255)}
    bonds=read(pr['covalent_bonds']); declared=set()
    for b in bonds:
        i,j=b['native_particle_indices']
        require([atoms[i]['id'],atoms[j]['id']]==[b['atom_a_id'],b['atom_b_id']], 'Bond ID map differs')
        adjacency[i].add(j); adjacency[j].add(i); declared.add(tuple(sorted((i,j))))
    bondforce=next(f for f in parent.getForces() if isinstance(f,mm.HarmonicBondForce))
    require(declared=={tuple(sorted(map(int,bondforce.getBondParameters(k)[:2])))
                       for k in range(bondforce.getNumBonds())}, 'Native bond graph differs')
    ff=ET.parse(verify(pr['forcefield'])).find('.//NonbondedForce')
    require(ff is not None, 'Native FF scaling declaration absent')
    return dict(id=source_id, preparation=preparation_pin, pr=pr, atoms=atoms,
                physical=physical, ids=ids, x=np.asarray([a['xyz_A'] for a in physical],float),
                parent=parent, nb=nb, parameters=parameters, adjacency=adjacency,
                bonds=declared, coulomb14=float(ff.attrib['coulomb14scale']),
                lj14=float(ff.attrib['lj14scale']))


def hydrogen_bond_inventory(source):
    """Geometry-only source-union inclusion rule, independent of metal energies."""
    atoms,x,adj,ids=source['atoms'],source['x'],source['adjacency'],source['ids']
    rows=[]
    donors=[i for i,a in enumerate(atoms) if a['element'] in ('N','O')
            and any(atoms[j]['element']=='H' for j in adj[i])]
    acceptors=[i for i,a in enumerate(atoms) if a['element']=='O']
    seed={ids[s] for s in HB_SEEDS}
    for d in donors:
        near=graph_distances(adj,d,depth=3)
        for a in acceptors:
            if (d not in seed and a not in seed) or a in near:
                continue
            da=float(np.linalg.norm(x[d]-x[a]))
            if da>HB_RULE['D_A_A']:
                continue
            for h in sorted(j for j in adj[d] if atoms[j]['element']=='H'):
                ha=float(np.linalg.norm(x[h]-x[a]));v=x[d]-x[h];w=x[a]-x[h]
                cosine=float(np.dot(v,w)/(np.linalg.norm(v)*np.linalg.norm(w)))
                angle=math.degrees(math.acos(max(-1.,min(1.,cosine))))
                if ha<=HB_RULE['H_A_A'] and angle>=HB_RULE['D_H_A_min_degrees']:
                    rows.append(dict(donor=atoms[d]['id'],hydrogen=atoms[h]['id'],acceptor=atoms[a]['id'],
                                     donor_residue=int(atoms[d]['id'].split('/')[1]),
                                     acceptor_residue=int(atoms[a]['id'].split('/')[1]),
                                     D_A_A=da,H_A_A=ha,D_H_A_degrees=angle))
    return rows


def choose_regions(sources):
    inventories={s:hydrogen_bond_inventory(v) for s,v in sources.items()}
    residues=set(range(SEED[0],SEED[1]+1))
    for rows in inventories.values():
        for r in rows:
            residues.update((r['donor_residue'],r['acceptor_residue']))
    baseline=(min(residues),max(residues)); expanded=(max(1,baseline[0]-1),min(17,baseline[1]+1))
    return {'baseline':baseline,'expanded':expanded}, dict(seed_interval=list(SEED),
        rule=HB_RULE,source_contacts=inventories,source_union_interval=list(baseline),
        expanded_interval=list(expanded),expanded_is_distinct=expanded!=baseline,
        selected_before_any_energy=True)


def region_definition(source, interval):
    first,last=interval;atoms=source['atoms']; ids=source['ids']
    selected={i for i,a in enumerate(atoms) if first<=int(a['id'].split('/')[1])<=last}
    cuts=[]
    for r,keep,retained,length in ((first-1,('C','O'),'C',1.09),(last+1,('N','H'),'N',1.01)):
        if not 1<=r<=17:
            continue
        for name in keep:
            require(f'A/{r}//{name}' in ids, 'Partial peptide boundary atom missing')
            selected.add(ids[f'A/{r}//{name}'])
        cuts.append({'retained_source_id':f'A/{r}//{retained}', 'omitted_source_id':f'A/{r}//CA',
                     'length_A':length, 'partial_residue':r, 'formal_QM_fragment_charge':0})
    require({ids[a] for a in DONORS}<=selected and ids['A/5//ND2'] in selected, 'Incomplete donor chemistry')
    actual={tuple(sorted((i,j))) for i in selected for j in source['adjacency'][i] if j not in selected}
    require(actual=={tuple(sorted((ids[c['retained_source_id']],ids[c['omitted_source_id']]))) for c in cuts},
            'Native covalent cut inventory differs')
    charges={i:a['native_charge_e'] for i,a in enumerate(source['parameters'])}
    full=[i for i,a in enumerate(atoms) if first<=int(a['id'].split('/')[1])<=last]
    expected=-sum(first<=r<=last for r in (3,7,11,13,14))+(1 if first==1 else 0)-(1 if last==17 else 0)
    require(abs(math.fsum(charges[i] for i in full)-expected)<=TOL['charge_e'], 'Full-region formal/native charge differs')
    omitted={ids[c['omitted_source_id']] for c in cuts};removed=selected|omitted;shifted=dict(charges);shift=[]
    for c in cuts:
        r=c['partial_residue'];res={i for i,a in enumerate(atoms) if int(a['id'].split('/')[1])==r}
        remain=res-removed; o=ids[c['omitted_source_id']]
        recipients=sorted(j for j in source['adjacency'][o] if j in remain and atoms[j]['element']!='H')
        require(len(recipients)==2, 'Expected two residue-local remaining heavy recipients')
        wanted=math.fsum(charges[i] for i in res);before=math.fsum(charges[i] for i in remain)
        increment=(wanted-before)/len(recipients)
        for i in recipients:shifted[i]+=increment
        require(abs(math.fsum(shifted[i] for i in remain)-wanted)<=TOL['charge_e'], 'Boundary charge closure failed')
        shift.append(dict(partial_residue=r,formal_QM_fragment_charge=0,
            original_residue_charge_e=wanted,removed_source_ids=[atoms[i]['id'] for i in sorted(res&removed)],
            recipients=[atoms[i]['id'] for i in recipients],each_increment_e=increment,
            retained_MM_charge_before_e=before,retained_MM_charge_after_e=math.fsum(shifted[i] for i in remain)))
    field=[dict(atoms[i],source_index=i,charge_e=shifted[i]) for i in range(255) if i not in removed]
    q=expected+3;fieldq=math.fsum(a['charge_e'] for a in field)
    require(abs(q+fieldq+2)<=TOL['charge_e'], 'Electronic plus field does not close to whole -2')
    return dict(interval=list(interval),qm_source_indices=sorted(selected),qm_charge=q,caps=cuts,
                removed_MM1_indices=sorted(omitted),field=field,charge_shift=shift,field_charge=fieldq)


def update_geometry(partition, physical_positions):
    """Pure current-coordinate map. Caps have no independent coordinates."""
    x=np.asarray(physical_positions,dtype=float)
    require(x.shape==(256,3) and np.isfinite(x).all(), 'Full finite physical (256,3) coordinates required')
    atoms=partition['physical_atoms'];ids={a['id']:i for i,a in enumerate(atoms)}
    core=[{'qm_index':0,'kind':'metal','source_id':'metal','physical_index':255,
           'element':'M','xyz_A':x[255].tolist()}]
    for i in partition['qm_source_indices']:
        core.append(dict(qm_index=len(core),kind='source',source_id=atoms[i]['id'],physical_index=i,
                         element=atoms[i]['element'],xyz_A=x[i].tolist()))
    caps=[]
    for c in partition['caps']:
        a,b=ids[c['retained_source_id']],ids[c['omitted_source_id']];length=c['length_A']
        ja,jb=cap_jacobians(x[a],x[b],length);v=x[b]-x[a]
        row=dict(c,qm_index=len(core),kind='cap',element='H',retained_index=a,omitted_index=b,
                 xyz_A=(x[a]+length*v/np.linalg.norm(v)).tolist(),
                 jacobian_retained=ja.tolist(),jacobian_omitted=jb.tolist())
        core.append(row);caps.append(row)
    field=[dict(a,xyz_A=x[a['source_index']].tolist()) for a in partition['field']]
    return dict(core=core,caps=caps,field=field,physical_atoms=[dict(a,xyz_A=x[i].tolist()) for i,a in enumerate(atoms)])


def scatter(mapping, qm_gradient, pc_gradient):
    """Return full physical gradient in the SAME units as both supplied gradients.

    Caller converts native Hartree/bohr or MACE eV/Å before combining components.
    This routine performs only the coordinate chain rule, never a force sign flip.
    """
    q=np.asarray(qm_gradient,float);p=np.asarray(pc_gradient,float)
    require(q.shape==(len(mapping['core']),3) and p.shape==(len(mapping['field']),3)
            and np.isfinite(q).all() and np.isfinite(p).all(), 'Gradient shape/finite gate failed')
    result=np.zeros((256,3))
    for a in mapping['core']:
        g=q[a['qm_index']]
        if a['kind']=='cap':
            result[a['retained_index']]+=np.asarray(a['jacobian_retained']).T@g
            result[a['omitted_index']]+=np.asarray(a['jacobian_omitted']).T@g
        else:result[a['physical_index']]+=g
    for i,a in enumerate(mapping['field']):result[a['source_index']]+=p[i]
    return result


def geometry_checks(partition, x):
    mapping=update_geometry(partition,x);records=[]
    for cap in mapping['caps']:
        for name,idx in [('retained',cap['retained_index']),('omitted',cap['omitted_index'])]:
            expected=np.asarray(cap['jacobian_'+name])
            for h in TOL['map_FD_steps_A']:
                observed=np.zeros((3,3))
                for axis in range(3):
                    plus=x.copy();minus=x.copy();plus[idx,axis]+=h;minus[idx,axis]-=h
                    cp=update_geometry(partition,plus)['core'][cap['qm_index']]['xyz_A']
                    cm=update_geometry(partition,minus)['core'][cap['qm_index']]['xyz_A']
                    observed[:,axis]=(np.asarray(cp)-np.asarray(cm))/(2*h)
                error=float(np.max(np.abs(observed-expected)))
                records.append(dict(cap_qm_index=cap['qm_index'],anchor=name,step_A=h,max_error=error,
                                    passed=error<=TOL['map_FD_max_error']))
        require(np.max(np.abs(np.asarray(cap['jacobian_retained'])+np.asarray(cap['jacobian_omitted'])-np.eye(3)))
                <=TOL['map_translation_max_error'], 'Cap translation closure failed')
    # Deterministic synthetic linear load checks the entire real/cap/field scatter.
    q=np.arange(len(mapping['core'])*3,dtype=float).reshape(-1,3)/max(1,len(mapping['core'])*3)
    p=np.arange(len(mapping['field'])*3,dtype=float).reshape(-1,3)/max(1,len(mapping['field'])*3)
    physical=scatter(mapping,q,p)
    closure=float(np.max(np.abs(physical.sum(axis=0)-q.sum(axis=0)-p.sum(axis=0))))
    require(closure<=1e-10, 'Complete map translation conservation failed')
    require(all(r['passed'] for r in records), 'Actual source cap two-step derivative check failed')
    return dict(cap_checks=records,synthetic_linear_load_translation_error=closure,
                molecular_force_qualification=False,energy_calls=0)


def ion_parameters(path):
    path=Path(path); require(path.name=='frcmod.ions234lm_126_tip3p', 'Pure12-6 TIP3P ion family required')
    active=False;found={}
    for line in path.read_text().splitlines():
        if line.strip()=='NONBON':active=True;continue
        fields=line.split()
        if active and fields and fields[0] in ('La3+','Dy3+'):
            require(fields[0] not in found, 'Duplicate ion parameter')
            radius,epsilon=map(float,fields[1:3]);found[fields[0]]=dict(lj_from_amber(radius,epsilon),
                source_line=line,Rmin_half_A=radius,epsilon_kcal_mol=epsilon)
    require(set(found)=={'La3+','Dy3+'}, 'Missing exact La/Dy ion parameters')
    return {m:found[m+'3+'] for m in METALS}


def interaction_ledger(source, partition, ionpars):
    qm=set(partition['qm_source_indices'])|{255};field={a['source_index']:a for a in partition['field']}
    particles=[]
    for i,a in enumerate(source['parameters']):
        particles.append(dict(a,physical_index=i,region='QM' if i in qm else 'MM',
                              charge_e=field[i]['charge_e'] if i in field else 0.))
    particles.append(dict(source['physical'][-1],kind='metal',physical_index=255,region='QM',
                          charge_e=0.,formal_charge_e=3.,endpoint_lj=ionpars))
    ids=source['ids'];omitted_to_retained={ids[c['omitted_source_id']]:ids[c['retained_source_id']] for c in partition['caps']}
    decisions=[];selected=[];counts=Counter()
    for f in source['parent'].getForces():
        kind=type(f).__name__
        if kind in ('NonbondedForce','CMMotionRemover'):continue
        require(kind in ('HarmonicBondForce','HarmonicAngleForce','PeriodicTorsionForce','CMAPTorsionForce'),
                'Unsupported native force '+kind)
        getter,num,size={'HarmonicBondForce':('getBondParameters','getNumBonds',2),
                         'HarmonicAngleForce':('getAngleParameters','getNumAngles',3),
                         'PeriodicTorsionForce':('getTorsionParameters','getNumTorsions',4),
                         'CMAPTorsionForce':('getTorsionParameters','getNumTorsions',8)}[kind]
        for k in range(getattr(f,num)()):
            raw=getattr(f,getter)(k); ordered=list(map(int,raw[1:] if kind=='CMAPTorsionForce' else raw[:size]))
            support=set(ordered);q=support&qm;m=support-qm
            if not m:keep=False;reason='wholly_QM'
            elif not q:keep=True;reason='wholly_MM'
            elif kind=='HarmonicBondForce':keep=True;reason='cross_radial_stiffness'
            elif all(i in omitted_to_retained and omitted_to_retained[i] in support for i in m):
                keep=False;reason='cross_represented_by_actual_caps'
            else:keep=True;reason='cross_unrepresented_exterior'
            row=dict(force=kind,index=k,atoms=ordered,source_ids=[particles[i]['id'] for i in ordered],
                     keep=keep,reason=reason)
            decisions.append(row);counts[kind+':'+reason]+=1
            if keep:selected.append({'force':kind,'index':k})
    require(counts['HarmonicBondForce:cross_radial_stiffness']==len(partition['caps']), 'Cross bond count differs')
    exceptions=[];distances={}
    for k in range(source['nb'].getNumExceptions()):
        i,j,q,s,e=source['nb'].getExceptionParameters(k);i,j=int(i),int(j)
        if i not in distances:distances[i]=graph_distances(source['adjacency'],i)
        d=distances[i].get(j);require(d in (1,2,3), 'Unsupported native exception topology')
        scale=source['coulomb14'] if d==3 else 0.
        nativeq=q.value_in_unit(unit.elementary_charge**2)
        require(abs(nativeq-particles[i]['native_charge_e']*particles[j]['native_charge_e']*scale)<=TOL['charge_e'],
                'Native 1-4 charge scaling differs')
        sigma=s.value_in_unit(unit.nanometer);epsilon=e.value_in_unit(unit.kilojoule_per_mole)
        expected_e=source['lj14']*math.sqrt(particles[i]['epsilon_kJ_mol']*particles[j]['epsilon_kJ_mol']) if d==3 else 0.
        require(abs(epsilon-expected_e)<=1e-10, 'Native LJ exception scaling differs')
        region='QM' if i in qm and j in qm else 'MM' if i not in qm and j not in qm else 'cross'
        exceptions.append(dict(native_index=k,atoms=[i,j],region=region,bond_distance=d,
            chargeprod_e2=scale*particles[i]['charge_e']*particles[j]['charge_e'],
            original_chargeprod_e2=nativeq,coulomb_scale=scale,sigma_nm=sigma,epsilon_kJ_mol=epsilon))
    require(len({tuple(sorted(e['atoms'])) for e in exceptions})==len(exceptions), 'Duplicate native exception')
    return {'particles.json':particles,'exceptions.json':exceptions,'bonded_terms.json':selected},decisions,dict(counts)


def native_mm_system(source, ledger, artifacts, metal):
    """Reuse the real builder; the empty external region is exactly zero MM."""
    particles=artifacts['particles.json']
    if any(a['region']=='MM' for a in particles):
        return build_system(ledger,artifacts,metal)
    require(not artifacts['bonded_terms.json'] and all(e['region']=='QM' for e in artifacts['exceptions.json']),
            'Empty external region retained MM interactions')
    from openmm.app import element
    system=mm.System()
    for i in range(255):system.addParticle(source['parent'].getParticleMass(i))
    system.addParticle(element.get_by_symbol(metal).mass)
    return system,dict(particle_count=256,real_QM=256,MM=0,bonded_counts={},MM_exception_count=0,
                       cross_replacements=0,constraints=0,energy_force_calls=0,
                       empty_MM_policy='identically zero; no empty CustomNonbonded interaction group')


def common_mm_definitions(systems,particles):
    """Exact physical equality; only inactive target sigma is normalized in copies."""
    target=255
    require(len(particles)==256 and particles[target]['id']=='metal'
            and particles[target]['kind']=='metal' and particles[target]['region']=='QM',
            'Inactive-sigma comparison requires the declared QM target at physical255')
    qset={i for i,a in enumerate(particles) if a['region']=='QM'}
    mset={i for i,a in enumerate(particles) if a['region']=='MM'}
    require(qset|mset==set(range(256)) and target not in mset, 'Physical region membership differs')
    serialized={};inactive={};raw_common={}
    for metal in METALS:
        system=mm.XmlSerializer.deserialize(verify(systems[metal]['system']).read_text())
        forces=list(system.getForces())
        serialized[metal]={str(g):[mm.XmlSerializer.serialize(f) for f in forces if f.getForceGroup()==g]
                           for g in range(3)}
        if not mset:
            require(not forces, 'Whole-Q system must have no MM forces')
            continue
        lj=[f for f in forces if f.getForceGroup()==1]
        require(len(lj)==1 and isinstance(lj[0],mm.NonbondedForce)
                and lj[0].getName()=='MM_LJ', 'Unique native MM-only LJ force required')
        xml=ET.fromstring(serialized[metal]['1'][0]);pp=xml.find('Particles')
        require(pp is not None and len(pp)==256, 'MM-only LJ particle count differs')
        atom=pp[target]
        require(set(atom.attrib)=={'q','sig','eps'} and atom.get('q')=='0' and atom.get('eps')=='0',
                'Target sigma can be ignored only when target MM charge and epsilon are exactly zero')
        require(all(len(xml.find(name))==0 for name in
                    ('GlobalParameters','ParticleOffsets','ExceptionOffsets')),
                'Parameterized MM-only LJ offsets are not supported by inactive-sigma proof')
        require(not any(target in (int(e.get('p1')),int(e.get('p2'))) for e in xml.find('Exceptions')),
                'Target MM-only LJ exception could bypass the zero-particle epsilon')
        cross=[f for f in forces if f.getForceGroup()==3 and f.getName()=='QM_MM_LJ_generic']
        require(len(cross)==1 and isinstance(cross[0],mm.CustomNonbondedForce)
                and cross[0].getNumInteractionGroups()==1, 'Unique QM-MM LJ interaction group required')
        left,right=cross[0].getInteractionGroupParameters(0)
        require(set(left)==qset and set(right)==mset, 'Actual QM-MM interaction membership differs')
        require(not any(target in tuple(map(int,cross[0].getExclusionParticles(k)))
                        for k in range(cross[0].getNumExclusions())),
                'Unexpected target QM-MM exclusion')
        replacements=[f for f in forces if f.getForceGroup()==3 and f.getName()=='QM_MM_LJ_native_exceptions']
        require(len(replacements)==1 and not any(target in tuple(map(int,replacements[0].getBondParameters(k)[:2]))
                    for k in range(replacements[0].getNumBonds())), 'Unexpected target QM-MM exception replacement')
        inactive[metal]={'physical_index':target,'id':'metal','force_group':1,'q':atom.get('q'),
            'epsilon':atom.get('eps'),'sigma_nm_serialized':atom.get('sig'),
            'target_MM_exceptions':0,'target_QMMM_exclusions':0,
            'actual_QMMM_membership':'QM set; every MM particle remains in the opposite set'}
        # Change only this comparison copy. Actual System/XML/group3 remain untouched.
        atom.set('sig','inactive_target_sigma')
        serialized[metal]['1_normalized']=[ET.tostring(xml,encoding='unicode')]
    common={}
    for g in range(3):
        key=str(g);raw_common[key]=serialized['La'][key]==serialized['Dy'][key]
        comparekey='1_normalized' if g==1 and mset else key
        common[key]=serialized['La'][comparekey]==serialized['Dy'][comparekey]
    return common,{'raw_serialized_equal':raw_common,'physical_equal':common,
        'comparison_only_inactive_target_sigma':inactive,'actual_systems_modified':False,
        'all_other_common_force_XML_and_exceptions_exact':all(common.values()),
        'group3_metal_dependent_QMMM_LJ_preserved':True}


def save_partition(source, name, interval, hb, ionpars, cfg, out):
    out.mkdir();definition=region_definition(source,interval)
    partition=dict(definition,protocol_id=PROTOCOL,source_id=source['id'],region_name=name,
                   physical_atoms=source['physical'],physical_source=source['pr']['physical_atoms'],
                   source_preparation=source['preparation'],selection=hb)
    mapping=update_geometry(partition,source['x']);checks=geometry_checks(partition,source['x'])
    complete=set(range(interval[0],interval[1]+1))
    require(all(r['donor_residue'] in complete and r['acceptor_residue'] in complete
                for rows in hb['source_contacts'].values() for r in rows), 'Source-union local H-bond coverage incomplete')
    ppin=write(out/'PARTITION.json',partition);mpin=write(out/'MAPPING.json',mapping)
    apin=write(out/'PHYSICAL.json',source['physical']);cpin=write(out/'core_mapping.json',mapping['core'])
    epin=write(out/'environment_atoms.json',mapping['field'])
    pc=out/'environment.pc';pc.write_text(str(len(mapping['field']))+'\n'+''.join(
        ' '.join(format(v,'.17g') for v in (a['charge_e'],*a['xyz_A']))+'\n' for a in mapping['field']))
    solvent=[a['native_charge_e'] for a in source['parameters']]+[3.]
    require(abs(math.fsum(solvent)+2)<=TOL['charge_e'], 'Global solvent charge closure failed')
    solventpin=write(out/'GLOBAL_SOLVENT_CHARGES.json',{'physical_ids':[a['id'] for a in source['physical']],
        'charges_e':solvent,'total_e':-2,'policy':'unshifted native ff19SB plus formal3 target; no cap charges'})
    artifacts,decisions,counts=interaction_ledger(source,definition,ionpars)
    artifactpins={n:write(out/n,v) for n,v in artifacts.items()}
    decisionpin=write(out/'BONDED_TERM_DECISIONS.json',decisions)
    ledger={'protocol_id':PROTOCOL,'native_system':source['pr']['native_peptide_system'],
        'forcefield':source['pr']['forcefield'],'source_preparation':source['preparation'],
        'partition':ppin,'artifacts':artifactpins,'bonded_decisions':decisionpin,'counts':counts,
        'ion_parameters':cfg['metal_parameters'],'boundary_policy':'fixed_cap_radial_bonds_retained_substitutable_angular_terms_omitted_v1',
        'MM_MM':'shifted-charge Coulomb and native LJ, native exceptions',
        'real_QM_MM':'native LJ and native exception replacements; no classical Coulomb',
        'QM_QM':'no MM terms','caps':'no independent MM particles or interactions',
        'C4_induction_included':False,'QM_MM_Coulomb_included':False,
        'embedding_cross_Coulomb':'all retained shifted point charges; no native cross-region electrostatic 1-4 scaling',
        'native_coulomb14scale':source['coulomb14'],'native_lj14scale':source['lj14'],
        'interaction_model_qualified':False,'energy_force_calls':0}
    lpin=write(out/'LEDGER.json',ledger);systems={};tasks=[]
    for metal,mult,ecp in [('La',1,46),('Dy',6,55)]:
        system,description=native_mm_system(source,ledger,artifacts,metal)
        sp=out/(metal+'_MM.xml');sp.write_text(mm.XmlSerializer.serialize(system));systems[metal]=dict(system=pin(sp),description=description)
        elements=[metal]+[a['element'] for a in mapping['core'][1:]]
        electrons=sum(NUMBERS[e] for e in elements)-definition['qm_charge']
        require((electrons-(mult-1))%2==0 and (electrons-ecp)%2==0, 'Physical/frozen-f electron parity differs')
        xp=out/(metal+'_Q.xyz');xp.write_text(str(len(elements))+'\n'+source['id']+' '+name+' source-derived region\n'+''.join(
            e+' '+' '.join(format(v,'.17g') for v in a['xyz_A'])+'\n' for e,a in zip(elements,mapping['core'])))
        tasks.append(dict(task_id=source['id']+'_'+name+'_'+metal,source_id=source['id'],region=name,metal=metal,
            charge=definition['qm_charge'],multiplicity=mult,physical_spin_2S=mult-1,oxidation_state=3,
            all_electron_count=electrons,frozen_f_ecp_core_electrons=ecp,frozen_f_explicit_electrons=electrons-ecp,
            native_effective_multiplicity=1,gfn2_uhf=0,physical_state_qualified=False,
            molecular_executor_or_parameters_bound=False,
            atom_count=len(elements),metal_index=0,xyz=pin(xp),pointcharges=pin(pc),core_mapping=cpin,
            physical_source=apin,partition=ppin,ledger=lpin,global_solvent_charges=solventpin))
    common,comparison=common_mm_definitions(systems,artifacts['particles.json'])
    require(all(common.values()), 'La/Dy MM common physical component definitions differ')
    recovery=read(cfg['recovery_evidence'])
    old=next(r for r in recovery['fixtures'] if r['source_id']==source['id'] and r['region']==name)
    require(all(systems[m]['system']['sha256']==old['systems'][m]['sha256'] for m in METALS),
            'Recovery changed an actual native MM XML; all physical and inactive parameters must be preserved')
    comparison['actual_native_XML_exactly_replayed_from_12950']=True
    result=dict(status='prepared',source_id=source['id'],region=name,interval=list(interval),partition=ppin,
        mapping=mpin,ledger=lpin,systems=systems,tasks=tasks,real_QM_atoms=len(definition['qm_source_indices'])+1,
        cap_atoms=len(definition['caps']),field_atoms=len(mapping['field']),charge_accounting={
            'QM_e':definition['qm_charge'],'shifted_field_e':definition['field_charge'],'whole_e':-2},
        geometry_checks=checks,common_MM_component_definition_checks=common,
        common_MM_comparison= comparison,admission_revision=ADMISSION_REVISION,
        energy_force_calls=0,molecular_gradient_qualified=False,source_coordinates_unchanged=True)
    write(out/'PREPARATION.json',result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['prepare'])
    parser.add_argument('--config',required=True,type=Path);parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();cfg=json.loads(args.config.read_text());out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    tick=time.monotonic();result={'protocol_id':PROTOCOL,'status':'execution_failure','rows':[],
        'declared_sources':list(SOURCES),'declared_regions':['baseline','expanded'],
        'declared_states':8,'config':pin(args.config),'implementation':pin(__file__),
        'admission_revision':ADMISSION_REVISION,'recovery_evidence':cfg['recovery_evidence'],
        'new_molecular_calls':0,'new_MM_energy_calls':0,'new_quantum_calls':0,'optimization_steps':0}
    try:
        initialize(cfg);parent=read(cfg['preparation_result'])
        require(parent['status']=='complete' and parent['slurm_job_id']=='12469', 'Complete original preparation12469 required')
        rows={r['source_id']:r for r in parent['rows']}
        require(set(rows)==set(SOURCES), 'Source denominator differs')
        sources={s:native_source(s,rows[s]['preparation']) for s in SOURCES}
        require([a['id'] for a in sources[SOURCES[0]]['physical']]==[a['id'] for a in sources[SOURCES[1]]['physical']],
                'Paired source physical identities differ')
        require(sources[SOURCES[0]]['pr']['native_peptide_system']['sha256']==sources[SOURCES[1]]['pr']['native_peptide_system']['sha256'],
                'Same native peptide Hamiltonian required across sources')
        regions,hb=choose_regions(sources);result['region_selection']=write(out/'REGION_SELECTION.json',hb)
        ionpars=ion_parameters(verify(cfg['metal_parameters']))
        for sid in SOURCES:
            directory=out/sid;directory.mkdir()
            for name in ('baseline','expanded'):
                row={'source_id':sid,'region':name,'status':'unavailable'}
                try:
                    require(name!='expanded' or hb['expanded_is_distinct'], 'Baseline is whole peptide; no distinct expanded region')
                    prepared=save_partition(sources[sid],name,regions[name],hb,ionpars,cfg,directory/name)
                    row.update(status='prepared',preparation=pin(directory/name/'PREPARATION.json'),tasks=prepared['tasks'])
                except Exception as exc:
                    row.update(reason=str(exc),traceback=traceback.format_exc())
                result['rows'].append(row);write(out/(sid+'_'+name+'_STATUS.json'),row)
        result['status']='complete' if len(result['rows'])==4 and all(r['status']=='prepared' for r in result['rows']) else 'unavailable'
        result.update(software={'python':sys.version,'numpy':np.__version__,'openmm':mm.__version__},
            plan=cfg['plan'],region_policy=cfg['region_policy'],helpers=cfg['helpers'],
            policy_references=cfg['policy_references'],
            complete_hybrid_qualified=False)
    except Exception as exc:
        result.update(error=str(exc),traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick,slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                  peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    write(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')},indent=2))
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':
    raise SystemExit(main())
