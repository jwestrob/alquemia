"""Collective protein-only constrained proposals; no force-field score addition."""
from __future__ import annotations

import argparse
from collections import deque
import fcntl
import json
import os
from pathlib import Path
import shutil
import time
import traceback

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsmr
from scipy.spatial import cKDTree

from affordable_common import InvalidArtifact, paired, read_json, record, verify, write_new, xyz

PROTOCOL = 'source_ff19SB_collective_matched_donor_target_proposal_v2'
TARGETS = ('origin', 'Ca_adaptive', 'La_adaptive')
SETTINGS = {
    'optimizer': 'sparse_projected_LBFGS_Armijo_bond_retraction_v1',
    'maximum_accepted_iterations': 200, 'lbfgs_history': 8,
    'maximum_trial_atom_step_A': .05, 'maximum_heavy_displacement_A': .8,
    'heavy_displacement_tolerance_A': 1e-7, 'fixed_atom_tolerance_A': 1e-12,
    'bond_retraction_tolerance_A': 1e-8, 'bond_final_tolerance_A': 1e-6,
    'retraction_iterations': 30, 'projection_atol': 1e-11,
    'projection_maxiter': 3000, 'projection_constraint_velocity_tolerance': 1e-7,
    'armijo': 1e-4, 'line_search_backtracks': 24,
    'stationarity_max_atom_kcal_mol_A': .1,
    'stationarity_rms_component_kcal_mol_A': .03,
    'maximum_accepted_energy_increase_kcal_mol': 1e-6,
    'new_overlap_heavy_A': 1., 'new_overlap_H_A': .55,
    'overlap_roundoff_A': 1e-7, 'source_chirality_volume_min_A3': 1e-5,
    'chirality_rule': 'actual_standard_amino_acid_nonGly_CA_and_Ile_Thr_CB_v1',
    'peptide_source_cosine_minimum': .1,
    'mobile_residue_radius_A': 8., 'direct_actual_peptide_neighbors': True,
    'forcefield_energy_added_to_score': False,
    'nonstationary_feasible_nonincreasing_candidates_admitted': True,
    'platform': 'OpenCL', 'precision': 'double', 'deterministic_forces': 'no_separate_OpenCL_toggle',
}


def wrap_angle(value):
    return (value + np.pi) % (2*np.pi) - np.pi


def dihedral(points):
    a,b,c,d = np.asarray(points)
    axis = c-b; axis /= np.linalg.norm(axis)
    v = a-b; v -= axis*np.dot(v,axis)
    w = d-c; w -= axis*np.dot(w,axis)
    if min(np.linalg.norm(v),np.linalg.norm(w)) < 1e-12:
        raise InvalidArtifact('undefined actual peptide dihedral')
    return float(np.arctan2(np.dot(np.cross(axis,v),w),np.dot(v,w)))


def adjacency(n, bonds):
    result = [set() for _ in range(n)]
    for a,b in bonds:
        result[a].add(b); result[b].add(a)
    return result


def covalent_guards(atoms, bonds, positions):
    """Source signs, never an assumed relaxed protein or reassigned stereochemistry."""
    neighbours = adjacency(len(atoms),bonds)
    chiral = []
    standard=set('ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
    for i,a in enumerate(atoms):
        if a['resname'] not in standard:raise InvalidArtifact('unsupported nonstandard parent stereochemistry')
        is_alpha=a['name']=='CA' and a['resname']!='GLY'
        is_beta=a['name']=='CB' and a['resname'] in ('ILE','THR')
        if not (is_alpha or is_beta):
            continue
        expected={'N','C','CB'} if is_alpha else ({'CA','CG1','CG2'} if a['resname']=='ILE' else {'CA','OG1','CG2'})
        ns=[atoms[j] for j in neighbours[i]]
        if (a['element']!='C' or len(ns)!=4 or sum(n['element']=='H' for n in ns)!=1 or
                {n['name'] for n in ns if n['element']!='H'}!=expected or
                any(n['residue_id']!=a['residue_id'] for n in ns)):
            raise InvalidArtifact('actual standard stereocenter substituents differ: '+a['source_id'])
        ids = sorted(neighbours[i]); p = positions[ids]
        vol = float(np.linalg.det(p[:3]-p[3]))
        if abs(vol) <= SETTINGS['source_chirality_volume_min_A3']:
            raise InvalidArtifact('nearly planar source tetrahedral carbon: '+str(i))
        chiral.append({'center':i,'neighbours':ids,'source_signed_volume_A3':vol})
    peptide = []
    for a,b in bonds:
        ia,ib = a,b
        if atoms[a]['name']=='N' and atoms[b]['name']=='C': ia,ib = b,a
        if (atoms[ia]['name'],atoms[ib]['name']) != ('C','N'):
            continue
        if atoms[ia]['residue_id'] == atoms[ib]['residue_id']:
            continue
        left = [i for i in neighbours[ia] if atoms[i]['name']=='CA' and
                atoms[i]['residue_id']==atoms[ia]['residue_id']]
        right = [i for i in neighbours[ib] if atoms[i]['name']=='CA' and
                 atoms[i]['residue_id']==atoms[ib]['residue_id']]
        if len(left)!=1 or len(right)!=1:
            raise InvalidArtifact('actual peptide bond lacks unique flanking alpha carbons')
        ids = [left[0],ia,ib,right[0]]; omega = dihedral(positions[ids])
        if abs(np.cos(omega)) < SETTINGS['peptide_source_cosine_minimum']:
            raise InvalidArtifact('source peptide outside a declared cis/trans basin')
        peptide.append({'indices':ids,'source_omega_radian':omega})
    return {'chirality':chiral,'peptides':peptide}


def mobile_atoms(atoms,bonds,positions,core_heavy_positions,core_parent_indices):
    selected = set()
    tree = cKDTree(np.asarray(core_heavy_positions))
    for i,a in enumerate(atoms):
        if a['element']!='H' and tree.query(positions[i])[0] <= SETTINGS['mobile_residue_radius_A']:
            selected.add(a['residue_id'])
    selected.update(atoms[i]['residue_id'] for i in core_parent_indices)
    original = set(selected)
    for a,b in bonds:
        aa,bb=atoms[a],atoms[b]
        if aa['residue_id']==bb['residue_id'] or {aa['name'],bb['name']}!={'C','N'}:
            continue
        if aa['residue_id'] in original: selected.add(bb['residue_id'])
        if bb['residue_id'] in original: selected.add(aa['residue_id'])
    return [i for i,a in enumerate(atoms) if a['residue_id'] in selected], sorted(selected)


class Constraints:
    """All source bond lengths, with fixed coordinates eliminated exactly."""
    def __init__(self, source, initial, bonds, free):
        self.source=np.asarray(source,dtype=float); self.initial=np.asarray(initial,dtype=float)
        self.free=np.asarray(free,dtype=int); self.bonds=np.asarray(bonds,dtype=int)
        self.fixed=np.array(sorted(set(range(len(source)))-set(free)),dtype=int)
        self.slot=np.full(len(source),-1,dtype=int); self.slot[self.free]=np.arange(len(self.free))
        self.lengths=np.linalg.norm(self.source[self.bonds[:,0]]-self.source[self.bonds[:,1]],axis=1)
        if np.min(self.lengths)<1e-6: raise InvalidArtifact('degenerate source bond')
        self.dynamic=np.flatnonzero(np.any(self.slot[self.bonds]>=0,axis=1))
        fixed_bonds=np.flatnonzero(np.all(self.slot[self.bonds]<0,axis=1))
        errors=np.abs(np.linalg.norm(self.initial[self.bonds[:,0]]-self.initial[self.bonds[:,1]],axis=1)-self.lengths)
        if len(fixed_bonds) and np.max(errors[fixed_bonds])>SETTINGS['bond_final_tolerance_A']:
            raise InvalidArtifact('fixed donor/outer target violates immutable source bond lengths')

    def residual_jacobian(self,positions):
        bs=self.bonds[self.dynamic];delta=positions[bs[:,0]]-positions[bs[:,1]]
        lengths=np.linalg.norm(delta,axis=1)
        if np.any(lengths<1e-6): raise InvalidArtifact('collapsed trial bond')
        u=delta/lengths[:,None]; rows=[];cols=[];vals=[]
        for endpoint,sign in ((0,1.),(1,-1.)):
            slots=self.slot[bs[:,endpoint]]; valid=np.flatnonzero(slots>=0)
            rows.extend(np.repeat(valid,3)); cols.extend((3*slots[valid,None]+np.arange(3)).ravel())
            vals.extend((sign*u[valid]).ravel())
        jac=sparse.coo_matrix((vals,(rows,cols)),shape=(len(bs),3*len(self.free))).tocsr()
        return lengths-self.lengths[self.dynamic],jac

    def project(self,positions,vector):
        _,jac=self.residual_jacobian(positions); v=np.asarray(vector).reshape(-1)
        solve=lsmr(jac.T,v,atol=SETTINGS['projection_atol'],btol=SETTINGS['projection_atol'],
                   maxiter=SETTINGS['projection_maxiter'])
        answer=v-jac.T@solve[0]
        residual=float(np.max(np.abs(jac@answer),initial=0))
        if not np.isfinite(answer).all() or residual>SETTINGS['projection_constraint_velocity_tolerance']*max(1.,np.linalg.norm(v)):
            raise InvalidArtifact('sparse tangent projection unresolved')
        return answer,{'lsmr_stop':int(solve[1]),'iterations':int(solve[2]),'max_constraint_velocity':residual}

    def retract(self,trial):
        p=np.asarray(trial).copy();p[self.fixed]=self.initial[self.fixed]
        for iteration in range(SETTINGS['retraction_iterations']+1):
            errors,jac=self.residual_jacobian(p)
            maximum=float(np.max(np.abs(errors),initial=0))
            if maximum<=SETTINGS['bond_retraction_tolerance_A']:
                return p,{'iterations':iteration,'max_bond_error_A':maximum}
            if iteration==SETTINGS['retraction_iterations']: break
            solve=lsmr(jac,-errors,atol=SETTINGS['projection_atol'],btol=SETTINGS['projection_atol'],
                       maxiter=SETTINGS['projection_maxiter'])
            step=solve[0].reshape((-1,3))
            if not np.isfinite(step).all(): break
            p[self.free]+=step
        raise InvalidArtifact('nonlinear source-bond retraction did not converge')


class Geometry:
    def __init__(self, atoms,bonds,source,initial,free,fixed_nonprotein):
        self.atoms=atoms; self.bonds=np.asarray(bonds,dtype=int)
        self.source=np.asarray(source);self.initial=np.asarray(initial)
        self.free=list(free);self.free_set=set(free)
        self.fixed=sorted(set(range(len(atoms)))-self.free_set)
        self.heavy=[i for i,a in enumerate(atoms) if a['element']!='H']
        self.fixed_nonprotein=fixed_nonprotein
        self.other=np.array([a['xyz_A'] for a in fixed_nonprotein],dtype=float).reshape((-1,3))
        self.all_source=np.concatenate((self.source,self.other))
        self.symbols=[a['element'] for a in atoms]+[a['element'] for a in fixed_nonprotein]
        self.guards=covalent_guards(atoms,bonds,self.source)
        neighbours=adjacency(len(atoms),bonds); self.excluded=set()
        for i in range(len(atoms)):
            near=set(neighbours[i])
            for j in neighbours[i]: near.update(neighbours[j])
            self.excluded.update((min(i,j),max(i,j)) for j in near if j!=i)

    def check(self,positions):
        p=np.asarray(positions)
        if p.shape!=self.source.shape or not np.isfinite(p).all():
            raise InvalidArtifact('nonfinite/incompatible full-parent coordinates')
        bond=float(np.max(np.abs(np.linalg.norm(p[self.bonds[:,0]]-p[self.bonds[:,1]],axis=1)-
                       np.linalg.norm(self.source[self.bonds[:,0]]-self.source[self.bonds[:,1]],axis=1)),initial=0))
        fixed=float(np.max(np.abs(p[self.fixed]-self.initial[self.fixed]),initial=0))
        displacement=np.linalg.norm(p[self.heavy]-self.source[self.heavy],axis=1)
        maximum=float(np.max(displacement,initial=0));chiral=[];peptides=[]
        for row in self.guards['chirality']:
            v=p[row['neighbours']]; vol=float(np.linalg.det(v[:3]-v[3]))
            if vol*row['source_signed_volume_A3']<=0 or abs(vol)<SETTINGS['source_chirality_volume_min_A3']:
                chiral.append(row['center'])
        for row in self.guards['peptides']:
            angle=dihedral(p[row['indices']])
            if np.cos(angle)*np.cos(row['source_omega_radian'])<=0: peptides.append(row['indices'])
        allp=np.concatenate((p,self.other));clashes=[]
        for a,b in sorted(cKDTree(allp).query_pairs(SETTINGS['new_overlap_heavy_A'])):
            if (a not in self.free_set and b not in self.free_set) or (a,b) in self.excluded: continue
            cutoff=SETTINGS['new_overlap_H_A'] if 'H' in (self.symbols[a],self.symbols[b]) else SETTINGS['new_overlap_heavy_A']
            distance=float(np.linalg.norm(allp[a]-allp[b]));original=float(np.linalg.norm(self.all_source[a]-self.all_source[b]))
            if distance<min(cutoff,original)-SETTINGS['overlap_roundoff_A']:
                clashes.append({'indices':[a,b],'distance_A':distance,'source_A':original,'nonprotein_contact':b>=len(p)})
        result={'source_bond_max_error_A':bond,'fixed_max_error_A':fixed,
                'maximum_heavy_displacement_A':maximum,'boundary_heavy_count':int(sum(displacement>=.8-1e-5)),
                'chirality_failures':chiral,'peptide_basin_failures':peptides,'new_severe_overlaps':clashes}
        result['pass']=bool(bond<=SETTINGS['bond_final_tolerance_A'] and fixed<=SETTINGS['fixed_atom_tolerance_A'] and
                            maximum<=.8+SETTINGS['heavy_displacement_tolerance_A'] and not chiral and not peptides and not clashes)
        return result


class ParentEnergy:
    def __init__(self, system_path,threads,output):
        import openmm as mm
        from openmm import unit
        self.mm=mm;self.unit=unit;self.calls=0;self.requests=0;self.seconds=0.;self.threads=int(threads)
        self.output=Path(output);self.output.mkdir(parents=True,exist_ok=False)
        system=mm.XmlSerializer.deserialize(Path(system_path).read_text())
        if system.getNumConstraints()!=0: raise InvalidArtifact('parent contains unexpected model constraints')
        self.integrator=mm.VerletIntegrator(.001)
        platform=mm.Platform.getPlatformByName(SETTINGS['platform'])
        properties={'Precision':'double'}
        devices=[dict(d) for d in platform.getDevices(properties)]
        if len(devices)!=1 or devices[0].get('DeviceName')!='NVIDIA H200':
            raise InvalidArtifact('expected exactly one allocated H200 OpenCL device: '+str(devices))
        self.context=mm.Context(system,self.integrator,platform,properties)
        self.runtime={'openmm_version':mm.__version__,'platform':platform.getName(),
                      'properties':{key:platform.getPropertyValue(self.context,key) for key in ('Precision','DeviceIndex','OpenCLPlatformName')},
                      'device_name':platform.getPropertyValue(self.context,'DeviceName'),
                      'visible_devices':devices,'threads_declared':int(threads)}

    def evaluate(self,p):
        self.requests+=1; start=time.monotonic()
        entry={'request':self.requests,'status':'failed'}
        artifact=self.output/f'{self.requests:06d}.npz'
        try:
            self.context.setPositions(np.asarray(p)*self.unit.angstrom)
            state=self.context.getState(getEnergy=True,getForces=True)
            energy=float(state.getPotentialEnergy().value_in_unit(self.unit.kilocalories_per_mole))
            forces=np.asarray(state.getForces(asNumpy=True).value_in_unit(self.unit.kilocalories_per_mole/self.unit.angstrom))
            self.calls+=1
            if not np.isfinite(energy) or not np.isfinite(forces).all():raise InvalidArtifact('nonfinite actual ff19SB energy/force')
            np.savez_compressed(artifact,positions_A=p,gradient_kcal_mol_A=-forces)
            entry.update(status='complete',energy_kcal_mol=energy,coordinates_gradient=record(artifact))
        except Exception as exc:
            np.savez_compressed(artifact,positions_A=p)
            entry.update(error=str(exc),coordinates=record(artifact));raise
        finally:
            elapsed=time.monotonic()-start;self.seconds+=elapsed;entry['wall_seconds']=elapsed
            with (self.output/'calls.jsonl').open('a') as handle:handle.write(json.dumps(entry,sort_keys=True)+'\n')
        self.last_evaluation=entry
        return energy,-forces

    def close(self):
        del self.context;del self.integrator


def lbfgs_direction(gradient,history):
    q=gradient.copy();alphas=[]
    for s,y,rho in reversed(history):
        alpha=rho*np.dot(s,q);alphas.append(alpha);q-=alpha*y
    scale=np.dot(history[-1][0],history[-1][1])/np.dot(history[-1][1],history[-1][1]) if history else 1.
    r=scale*q
    for (s,y,rho),alpha in zip(history,reversed(alphas)):
        r+=s*(alpha-rho*np.dot(y,r))
    return -r


def optimize(initial,constraints,geometry,evaluator,trace_path):
    p=np.asarray(initial).copy(); check=geometry.check(p)
    if not check['pass']:raise InvalidArtifact('initial physical constraints fail before force evaluation: '+json.dumps(check))
    energy,gradient=evaluator.evaluate(p);initial_energy=energy;history=[];accepted=0;status='iteration_limit'
    initial_evaluation=evaluator.last_evaluation;final_evaluation=initial_evaluation
    trial_count=0;rejected=[]; first_eval_seconds=evaluator.seconds;trace=[]
    def log(row):
        trace.append(row)
        with Path(trace_path).open('a') as handle:handle.write(json.dumps(row,sort_keys=True)+'\n')
    while True:
        pg,projection=constraints.project(p,gradient[constraints.free].ravel())
        pg3=pg.reshape((-1,3));maximum=float(np.max(np.linalg.norm(pg3,axis=1),initial=0));rms=float(np.sqrt(np.mean(pg*pg)))
        log({'accepted_iteration':accepted,'energy_kcal_mol':energy,'parent_evaluations':evaluator.calls,
             'projected_max_atom_kcal_mol_A':maximum,'projected_rms_component_kcal_mol_A':rms,
             'projection':projection,'geometry':geometry.check(p)})
        if maximum<=SETTINGS['stationarity_max_atom_kcal_mol_A'] and rms<=SETTINGS['stationarity_rms_component_kcal_mol_A']:
            status='stationary';break
        if accepted>=SETTINGS['maximum_accepted_iterations']:break
        direction,_=constraints.project(p,lbfgs_direction(pg,history))
        if np.dot(pg,direction)>=-1e-14:
            history=[];direction=-pg
        speed=float(np.max(np.linalg.norm(direction.reshape((-1,3)),axis=1),initial=0))
        if speed<1e-14:status='zero_feasible_direction';break
        direction*=min(1.,SETTINGS['maximum_trial_atom_step_A']/speed)
        slope=float(np.dot(pg,direction));step=1.;next_state=None
        for backtrack in range(SETTINGS['line_search_backtracks']):
            trial_count+=1
            trial=p.copy();trial[constraints.free]+=step*direction.reshape((-1,3))
            reason=None;trial_check=None
            try:
                trial,retraction=constraints.retract(trial);trial_check=geometry.check(trial)
                if not trial_check['pass']:raise InvalidArtifact('physical_guard')
                new_energy,new_gradient=evaluator.evaluate(trial)
                if new_energy<=energy+SETTINGS['armijo']*step*slope:
                    next_state=(trial,new_energy,new_gradient,retraction,backtrack);break
                reason='Armijo_energy'
            except InvalidArtifact as exc:reason=str(exc)
            rejected.append({'after_iteration':accepted,'backtrack':backtrack,'step':step,'reason':reason,
                             'geometry':trial_check})
            step*=.5
        if next_state is None:status='line_search_no_admissible_decrease';break
        trial,new_energy,new_gradient,retraction,backtrack=next_state
        new_pg,_=constraints.project(trial,new_gradient[constraints.free].ravel())
        s=(trial[constraints.free]-p[constraints.free]).ravel();y=new_pg-pg;sy=float(np.dot(s,y))
        if sy>1e-10*np.linalg.norm(s)*np.linalg.norm(y) and sy>1e-16:
            history.append((s,y,1./sy));history=history[-SETTINGS['lbfgs_history']:]
        p,energy,gradient=trial,new_energy,new_gradient;accepted+=1
        final_evaluation=evaluator.last_evaluation
    final_check=geometry.check(p)
    admitted=bool(final_check['pass'] and energy<=initial_energy+SETTINGS['maximum_accepted_energy_increase_kcal_mol'])
    return {'status':status,'candidate_admitted':admitted,'stationary':status=='stationary',
            'positions_A':p.tolist(),'initial_energy_kcal_mol':initial_energy,'final_energy_kcal_mol':energy,
            'initial_evaluation':initial_evaluation,'final_evaluation':final_evaluation,
            'work_kcal_mol':energy-initial_energy,'accepted_iterations':accepted,'trial_count':trial_count,
            'parent_energy_force_calls':evaluator.calls,'parent_energy_force_requests':evaluator.requests,
            'parent_evaluation_seconds':evaluator.seconds,
            'first_declared_evaluation_seconds':first_eval_seconds,'final_geometry':final_check,
            'final_projected_max_atom_kcal_mol_A':maximum,'final_projected_rms_component_kcal_mol_A':rms,
            'rejected_trials':rejected,'forcefield_energy_added_to_score':False}


def reconstruct_context(case,positions,metal='Ca'):
    """Rebuild source atoms and existing offset caps; no independent cap motion."""
    mapping=read_json(verify(case['context_parent_mapping']))
    g=read_json(verify(case['context_maps'][metal]))['context']
    coordinates=np.asarray(g['core_positions_A'],dtype=float).copy()
    for row in mapping['source_context_atoms']:
        coordinates[row['context_index']]=positions[row['parent_index']]
    for row in mapping['caps']:
        a,b=positions[[row['retained_parent_index'],row['omitted_parent_index']]]
        distance=np.linalg.norm(b-a)
        if distance<1e-6:raise InvalidArtifact('collapsed cap boundary bond')
        coordinates[row['context_index']]+=a+row['length_A']*(b-a)/distance-np.asarray(row['q0_xyz_A'])
    return coordinates


def load_parent(case):
    parent=case['parent'];topology=read_json(verify(parent['topology']))
    atoms=read_json(verify(parent['atoms']));positions=np.asarray(read_json(verify(parent['positions_A'])),dtype=float)
    ids=read_json(verify(parent['source_atom_ids']))
    if len(atoms)!=len(positions) or [a['source_id'] for a in atoms]!=ids or len(set(ids))!=len(ids):
        raise InvalidArtifact('parent source order/coordinates inconsistent')
    if atoms!=topology['atoms']:raise InvalidArtifact('parent topology atom record differs')
    atoms=[dict(a,residue_id=a['source_id'].rsplit('/',1)[0]) for a in atoms]
    bonds=[b['indices'] for b in topology['bonds']]
    lengths=np.linalg.norm(positions[np.asarray(bonds)[:,0]]-positions[np.asarray(bonds)[:,1]],axis=1)
    if not np.allclose(lengths,[b['source_length_A'] for b in topology['bonds']],atol=1e-12,rtol=0):
        raise InvalidArtifact('source bond length record differs')
    verify(parent['system']);verify(parent['forcefield']);verify(case['source']);verify(case['source_preparation'])
    mobile=read_json(verify(case['mobile_set']))
    computed,residues=mobile_atoms(atoms,bonds,positions,mobile['canonical_core_heavy_positions_A'],mobile['core_parent_indices'])
    if computed!=mobile['mobile_parent_indices']:raise InvalidArtifact('8A complete-residue shell does not replay')
    donors=read_json(verify(case['donors']))['atoms'];donorids=[d['parent_index'] for d in donors]
    if len(set(donorids))!=len(donorids) or not set(donorids)<=set(computed):
        raise InvalidArtifact('target donor mapping outside mobile source shell')
    if any(ids[d['parent_index']]!=d['source_id'] or atoms[d['parent_index']]['element']!='O' for d in donors):
        raise InvalidArtifact('donor source identity differs')
    if any(d['role'] not in ('anchor_glutamate','anchor_asparagine','extra_acidic_ligand_homolog') for d in donors):
        raise InvalidArtifact('nondeclared coordination target role')
    expected={'anchor_glutamate':{'OE1','OE2'},'anchor_asparagine':{'OD1'}}
    if any(d['role']=='extra_acidic_ligand_homolog' for d in donors):expected['extra_acidic_ligand_homolog']={'OD1','OD2'}
    if {role:{atoms[d['parent_index']]['name'] for d in donors if d['role']==role} for role in expected}!=expected:
        raise InvalidArtifact('uniform donor atom rule differs')
    nonprotein=read_json(verify(case['nonprotein_inventory']))
    fixed_other=nonprotein['fixed_context_atoms']
    # Add any source nonprotein atom outside the fixed local context, without duplicates.
    for atom in nonprotein['source_atoms']:
        if not any(atom['element']==a['element'] and np.allclose(atom['xyz_A'],a['xyz_A'],atol=1e-12,rtol=0) for a in fixed_other):
            fixed_other.append(atom)
    ca,la=case['origins']['Ca'],case['origins']['La']
    paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
    if ca['multiplicity']!=1 or la['multiplicity']!=1:raise InvalidArtifact('changed endpoint multiplicity')
    if read_json(verify(case['context_maps']['Ca']))!=read_json(verify(case['context_maps']['La'])):
        raise InvalidArtifact('paired physical source maps differ')
    for metal in ('Ca','La'):
        actual=np.array([a[1:] for a in xyz(verify(case['origins'][metal]['xyz']))])
        if not np.allclose(actual,reconstruct_context(case,positions,metal),atol=1e-12,rtol=0):
            raise InvalidArtifact('full parent q0 fails unchanged context/cap replay')
    return {'atoms':atoms,'positions':positions,'bonds':bonds,'mobile':computed,'donors':donorids,
            'free':sorted(set(computed)-set(donorids)),'fixed_nonprotein':fixed_other,
            'mobile_residues':residues}


def target_problem(case,target,parent=None):
    data=load_parent(case) if parent is None else parent
    initial=np.array(read_json(verify(case['targets'][target]['parent_positions_A'])),dtype=float)
    if initial.shape!=data['positions'].shape:raise InvalidArtifact('target source atom shape differs')
    exterior=sorted(set(range(len(initial)))-set(data['mobile']))
    if not np.allclose(initial[exterior],data['positions'][exterior],atol=1e-12,rtol=0):
        raise InvalidArtifact('declared target moves frozen exterior source atoms')
    if not np.allclose(initial[data['donors']],case['targets'][target]['donor_positions_A'],atol=1e-12,rtol=0):
        raise InvalidArtifact('target donor positions differ')
    for metal in ('Ca','La'):
        old=xyz(verify(case['targets'][target]['context_coordinates'][metal]))
        if not np.allclose(reconstruct_context(case,initial,metal),[a[1:] for a in old],atol=1e-12,rtol=0):
            raise InvalidArtifact('target source/cap reconstruction differs')
    constraints=Constraints(data['positions'],initial,data['bonds'],data['free'])
    geometry=Geometry(data['atoms'],data['bonds'],data['positions'],initial,data['free'],data['fixed_nonprotein'])
    check=geometry.check(initial)
    if not check['pass']:raise InvalidArtifact('initial physical guards fail: '+json.dumps(check))
    return data,initial,constraints,geometry


def prepare(parents,inputs,agreement,target_mapping,plan,output,technical_addendum):
    parent=read_json(parents);registry=read_json(inputs);out=Path(output).resolve()
    out.mkdir(parents=True,exist_ok=False)
    cases=parent['cases'];ids=[c['case_id'] for c in cases]
    if len(ids)!=8 or ids!=[c['case_id'] for c in registry['cases']]:raise InvalidArtifact('common eight source order differs')
    tasks=[];case_records=[]
    for c in cases:
        case_path=Path(parents).parent/c['case_id']/'CASE.json'
        if case_path.is_file():
            if read_json(case_path)!=c:raise InvalidArtifact('parent case pin and manifest differ')
        else:
            case_path=out/'unavailable_cases'/f"{c['case_id']}.json";write_new(case_path,c)
        pin=record(case_path);case_records.append(pin)
        error=None;data=None
        try:
            if c['status']!='prepared':raise InvalidArtifact(c.get('error','parent unavailable'))
            data=load_parent(c)
        except Exception as exc:error=str(exc)
        for target in TARGETS:
            task={'task_id':c['case_id']+'__'+target,'case_id':c['case_id'],'parent_case':pin,'target':target,
                  'status':'unavailable','output':str(out/'searches'/(c['case_id']+'__'+target))}
            try:
                if error:raise InvalidArtifact(error)
                _,initial,constraints,geometry=target_problem(c,target,data)
                task.update(status='prepared',initial_parent_positions_A=c['targets'][target]['parent_positions_A'],
                            initial_checks=geometry.check(initial),free_atom_count=len(data['free']),
                            mobile_atom_count=len(data['mobile']),donor_count=len(data['donors']))
            except Exception as exc:task['reason']=str(exc)
            tasks.append(task)
    implementation=out/'implementation';implementation.mkdir()
    pins={}
    for name in ('scaffold_accommodation.py','affordable_common.py'):
        dest=implementation/name;shutil.copyfile(Path(__file__).parent/name,dest);pins[name]=record(dest)
    result={'protocol_id':PROTOCOL,'settings':SETTINGS,'parent_manifest':record(parents),'inputs':record(inputs),
            'agreement':record(agreement),'target_mapping_agreement':record(target_mapping),'execution_plan':record(plan),
            'technical_addendum':record(technical_addendum),
            'case_records':case_records,'tasks':tasks,'implementation':pins,'source_count':8,'declared_search_count':24,
            'prepared_search_count':sum(t['status']=='prepared' for t in tasks),'output_root':str(out),
            'new_MACE_calls':0,'new_GFN2_calls':0,'new_DFT_calls':0,'production_changed':False}
    write_new(out/'manifest.json',result);return result


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS:raise InvalidArtifact('mechanics protocol/settings mismatch')
    for k in ('parent_manifest','inputs','agreement','target_mapping_agreement','execution_plan','technical_addendum'):verify(m[k])
    for pin in m['implementation'].values():verify(pin)
    if m['implementation']['scaffold_accommodation.py']['sha256']!=record(__file__)['sha256']:
        raise InvalidArtifact('use exact prepared implementation')
    ids=[c['case_id'] for c in read_json(verify(m['inputs']))['cases']]
    if [(t['case_id'],t['target']) for t in m['tasks']]!=[(cid,target) for cid in ids for target in TARGETS]:
        raise InvalidArtifact('finite task set/order differs')
    summaries=[]
    for pin in m['case_records']:
        c=read_json(verify(pin));data=None
        for task in [t for t in m['tasks'] if t['case_id']==c['case_id']]:
            if task['parent_case']!=pin:raise InvalidArtifact('task parent case differs')
            if task['status']!='prepared':continue
            if data is None:data=load_parent(c)
            _,initial,constraints,geometry=target_problem(c,task['target'],data)
            if geometry.check(initial)!=task['initial_checks']:raise InvalidArtifact('initial checks do not replay')
            if Path(task['output']).parent!=Path(m['output_root'])/'searches':raise InvalidArtifact('unsafe task output path')
            summaries.append({'task_id':task['task_id'],'free_atoms':len(data['free']),'initial_geometry':'pass'})
    return {'status':'ready','manifest':record(manifest),'sources':8,'declared_searches':24,
            'prepared_searches':len(summaries),'energy_force_calls':0,'checks':summaries}


def write_context(path,symbols,coordinates):
    with Path(path).open('x') as handle:
        handle.write(str(len(symbols))+'\n'+PROTOCOL+'\n')
        for symbol,p in zip(symbols,coordinates):handle.write(symbol+' '+' '.join(f'{v:.15f}' for v in p)+'\n')


def run_one(manifest,m,task):
    import resource
    out=Path(task['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();engine=None
    case=read_json(verify(task['parent_case']))
    receipt={'protocol_id':PROTOCOL,'manifest':record(manifest),'task_id':task['task_id'],'case_id':task['case_id'],
             'parent_case_id':task['case_id'],'parent_case':task['parent_case'],'parent_manifest':m['parent_manifest'],
             'target':task['target'],'status':'failed','candidate_admitted':False,'stationary':False,
             'initial_parent_positions_A':task['initial_parent_positions_A'],
             'source_atom_ids':case['parent']['source_atom_ids'],'parent_system':case['parent']['system'],
             'forcefield_energy_added_to_score':False,
             'allocation':{k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_JOB_NODELIST','SLURM_NTASKS','SLURM_CPUS_ON_NODE','SLURM_GPUS_ON_NODE','CUDA_VISIBLE_DEVICES')}}
    try:
        data,initial,constraints,geometry=target_problem(case,task['target'])
        engine=ParentEnergy(verify(case['parent']['system']),int(os.environ.get('SLURM_NTASKS','32')),out/'evaluations')
        receipt['runtime']=engine.runtime
        result=optimize(initial,constraints,geometry,engine,out/'trace.jsonl')
        positions=np.asarray(result.pop('positions_A'));write_new(out/'final_parent_positions_A.json',positions.tolist())
        write_new(out/'optimization.json',result)
        contexts={}
        for metal in ('Ca','La'):
            coords=reconstruct_context(case,positions,metal)
            symbols=[a[0] for a in xyz(verify(case['origins'][metal]['xyz']))]
            path=out/f'{metal}_context.xyz';write_context(path,symbols,coords);contexts[metal]=record(path)
        paired(verify(contexts['La']),verify(contexts['Ca']),case['origins']['La']['charge'],case['origins']['Ca']['charge'])
        receipt.update(status='complete',candidate_admitted=result['candidate_admitted'],stationary=result['stationary'],
                       optimizer_status=result['status'],checks=result['final_geometry'],
                       final_parent_positions_A=record(out/'final_parent_positions_A.json'),context_coordinates=contexts,
                       parent_work_kcal_mol=result['work_kcal_mol'],optimization=record(out/'optimization.json'),
                       initial_parent_energy_kcal_mol=result['initial_energy_kcal_mol'],
                       final_parent_energy_kcal_mol=result['final_energy_kcal_mol'],
                       initial_evaluation=result['initial_evaluation'],final_evaluation=result['final_evaluation'],
                       trace=record(out/'trace.jsonl'),unchanged_geometry=bool(np.array_equal(positions,initial)))
    except Exception as exc:
        receipt.update(error=str(exc),exception=type(exc).__name__)
        (out/'exception.txt').write_text(traceback.format_exc())
    finally:
        receipt.update(wall_seconds=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                       parent_energy_force_calls=engine.calls if engine else 0,
                       parent_energy_force_requests=engine.requests if engine else 0,
                       parent_evaluation_seconds=engine.seconds if engine else 0.)
        if engine is not None:
            if (out/'evaluations'/'calls.jsonl').exists():receipt['all_evaluations']=record(out/'evaluations'/'calls.jsonl')
            engine.close()
        write_new(out/'receipt.json',receipt)
    return receipt


def execute(manifest):
    m=read_json(manifest);validate(manifest)
    if not os.environ.get('SLURM_JOB_ID'):raise InvalidArtifact('actual mechanics requires the recorded Slurm GPU allocation')
    lock=Path(m['output_root'])/'execution.lock'
    with lock.open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        started=time.monotonic();rows=[]
        for task in m['tasks']:
            if task['status']!='prepared':continue
            path=Path(task['output'])/'receipt.json'
            if path.exists():
                receipt=read_json(path)
                if receipt['manifest']!=record(manifest) or receipt['task_id']!=task['task_id']:
                    raise InvalidArtifact('existing execution receipt belongs to another task')
                rows.append(record(path));continue
            if Path(task['output']).exists():raise InvalidArtifact('partial search retained; no automatic extra start')
            receipt=run_one(manifest,m,task);rows.append(record(path))
            print(json.dumps({'task':task['task_id'],'status':receipt['status'],
                              'optimizer_status':receipt.get('optimizer_status'),'admitted':receipt['candidate_admitted'],
                              'wall_seconds':receipt['wall_seconds'],'force_calls':receipt['parent_energy_force_calls']}),flush=True)
        out={'manifest':record(manifest),'receipts':rows,'wall_seconds':time.monotonic()-started,
             'job_id':os.environ['SLURM_JOB_ID'],'new_MACE_calls':0,'new_GFN2_calls':0}
        write_new(Path(m['output_root'])/f"execution_{os.environ['SLURM_JOB_ID']}.json",out)
    return out


def collect(manifest,output):
    m=read_json(manifest);validate(manifest);rows=[]
    for task in m['tasks']:
        row={k:task[k] for k in ('task_id','case_id','target','status')};row['candidate_admitted']=False
        path=Path(task['output'])/'receipt.json'
        if task['status']!='prepared':row['reason']=task.get('reason')
        elif not path.exists():row.update(status='missing',reason='no_terminal_receipt')
        else:
            receipt=read_json(path)
            if receipt['manifest']!=record(manifest) or receipt['task_id']!=task['task_id']:raise InvalidArtifact('collection receipt mismatch')
            row.update(status=receipt['status'],candidate_admitted=receipt['candidate_admitted'],receipt=record(path),
                       stationary=receipt['stationary'],optimizer_status=receipt.get('optimizer_status'),
                       parent_work_kcal_mol=receipt.get('parent_work_kcal_mol'),error=receipt.get('error'),
                       parent_energy_force_calls=receipt['parent_energy_force_calls'],
                       parent_energy_force_requests=receipt['parent_energy_force_requests'],wall_seconds=receipt['wall_seconds'])
            if receipt['status']=='complete':
                case=read_json(verify(task['parent_case']));_,_,_,geometry=target_problem(case,task['target'])
                positions=np.array(read_json(verify(receipt['final_parent_positions_A'])))
                if geometry.check(positions)!=receipt['checks']:raise InvalidArtifact('final physical checks do not replay')
                for metal in ('Ca','La'):
                    actual=xyz(verify(receipt['context_coordinates'][metal]))
                    if [a[0] for a in actual]!=[a[0] for a in xyz(verify(case['origins'][metal]['xyz']))]:raise InvalidArtifact('context atom state changed')
                    if not np.allclose([a[1:] for a in actual],reconstruct_context(case,positions,metal),atol=1e-12,rtol=0):raise InvalidArtifact('final cap/source reconstruction differs')
        rows.append(row)
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'rows':rows,'source_denominator':8,'search_denominator':24,
            'terminal_searches':sum(r['status'] in ('complete','failed','unavailable') for r in rows),
            'admitted_candidates':sum(r['candidate_admitted'] for r in rows),
            'fully_available_sources':[cid for cid in dict.fromkeys(r['case_id'] for r in rows)
                                      if all(r['candidate_admitted'] for r in rows if r['case_id']==cid)],
            'parent_energy_force_calls':sum(r.get('parent_energy_force_calls',0) for r in rows),
            'parent_energy_force_requests':sum(r.get('parent_energy_force_requests',0) for r in rows),
            'forcefield_energy_added_to_score':False,'classifier_gain':None,'production_changed':False}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    prep=sub.add_parser('prepare')
    for flag in ('parents','inputs','agreement','target-mapping','plan','output','technical-addendum'):prep.add_argument('--'+flag,required=True)
    for command in ('validate','execute','collect'):
        s=sub.add_parser(command);s.add_argument('--manifest',required=True)
        if command=='collect':s.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':result=prepare(a.parents,a.inputs,a.agreement,a.target_mapping,a.plan,a.output,a.technical_addendum)
    elif a.command=='validate':result=validate(a.manifest)
    elif a.command=='execute':result=execute(a.manifest)
    else:result=collect(a.manifest,a.output)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
