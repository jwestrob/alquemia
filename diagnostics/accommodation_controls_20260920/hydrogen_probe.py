#!/usr/bin/env python3
"""Isolated, declared native minimizer replay and fixed-potential H repair."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import traceback

import numpy as np
from scipy.spatial import cKDTree
import openmm as mm
from openmm import unit
from openmm.app import PDBFile


def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def save(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n')


def state(context,movable):
    s=context.getState(positions=True,energy=True,forces=True)
    x=np.asarray(s.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    f=np.asarray(s.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
    e=float(s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    if not np.isfinite(x).all() or not np.isfinite(f).all() or not np.isfinite(e):
        raise ValueError('nonfinite native preparation potential/positions/forces')
    summary={'potential_energy_kJ_mol':e,'movable_force_RMS_kJ_mol_nm':float(np.sqrt(np.mean(f[movable]**2))),
             'movable_force_max_kJ_mol_nm':float(np.max(np.linalg.norm(f[movable],axis=1)))}
    return x,f,e,summary


def overlaps(x,movable):
    return sorted((i,j) for i,j in cKDTree(x).query_pairs(.045) if movable[i] or movable[j])


class Reporter(mm.MinimizationReporter):
    def __init__(self,movable):
        super().__init__();self.rows=[];self.movable=movable
    def report(self,iteration,x,grad,args):
        g=np.asarray(grad).reshape((-1,3))[self.movable]
        self.rows.append({'iteration':int(iteration),'gradient_RMS_kJ_mol_nm':float(np.sqrt(np.mean(g*g))),
                          'args':{str(k):float(v) for k,v in args.items()}})
        return False


def numerical_repair(context,initial,movable,native_minimize,output):
    context.setPositions(initial*unit.nanometer)
    trace=[]
    for iteration in range(200):
        x,f,e,summary=state(context,movable)
        pairs=overlaps(x,movable)
        if not pairs: break
        # Componentwise atom displacement bound; direction remains downhill.
        step=np.zeros_like(x)
        norms=np.linalg.norm(f[movable],axis=1)
        step[movable]=.002*f[movable]/np.maximum(norms[:,None],1e-300)
        slope=float(np.sum(f*step))
        if not slope>0:raise ValueError('no descending movable direction')
        accepted=False
        for backtrack in range(21):
            scale=2.**(-backtrack);trial=x+scale*step
            context.setPositions(trial*unit.nanometer)
            _,_,et,_=state(context,movable)
            if et<=e-1e-4*scale*slope:
                accepted=True;break
        if not accepted:
            context.setPositions(x*unit.nanometer)
            raise ValueError('Armijo line search failed')
        trace.append({'iteration':iteration,'initial_pairs_below_0p45A':len(pairs),
                      'energy_before_kJ_mol':e,'energy_after_kJ_mol':et,
                      'backtracks':backtrack,'maximum_step_A':.02*scale})
    before,_,_,before_summary=state(context,movable)
    save(output/'preconditioning_trace.json',trace)
    np.save(output/'preconditioned_positions_nm.npy',before)
    if overlaps(before,movable):raise ValueError('200-step preconditioning did not clear overlap gate')
    reporter=Reporter(movable)
    native_minimize(context,1.,50,reporter)
    final,_,_,final_summary=state(context,movable)
    np.save(output/'repaired_positions_nm.npy',final)
    save(output/'repair_minimizer_iterations.json',reporter.rows)
    return {'preconditioning_accepted_steps':len(trace),'before_native_minimizer':before_summary,
            'after_native_minimizer':final_summary,'native_accepted_iterations':len(reporter.rows),
            'fixed_atom_max_displacement_A':float(np.max(np.abs(final[~movable]-initial[~movable])))*10,
            'maximum_H_displacement_A':float(np.max(np.linalg.norm(final[movable]-initial[movable],axis=1)))*10,
            'pairs_below_0p45A':len(overlaps(final,movable))}


def geometry(topology,positions,movable,initial):
    aa=list(topology.atoms());neighbors=[[] for a in aa]
    for a,b in topology.bonds():neighbors[a.index].append(b.index);neighbors[b.index].append(a.index)
    distances=[];bad=[]
    for a in aa:
        if not movable[a.index]:continue
        heavy=[j for j in neighbors[a.index] if aa[j].element.symbol!='H']
        if len(heavy)!=1:bad.append(a.index)
        else:distances.append(float(np.linalg.norm(positions[a.index]-positions[heavy[0]]))*10)
    maximum=float(np.max(np.linalg.norm(positions[movable]-initial[movable],axis=1)))*10
    heavy_fixed=bool(np.array_equal(positions[~movable],initial[~movable]))
    pairs=overlaps(positions,movable)
    return {'H_parent_distance_range_A':[min(distances),max(distances)],
            'H_without_exactly_one_heavy_parent':bad,'heavy_coordinates_exact':heavy_fixed,
            'maximum_H_displacement_A':maximum,'pairs_below_0p45A':len(pairs),
            'passes_geometry':bool(heavy_fixed and not bad and not pairs and min(distances)>=.8
                                   and max(distances)<=1.5 and maximum<=2.5)}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--sample',type=int,choices=[3,4],required=True)
    args=ap.parse_args();root=args.root.resolve();out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    case=f'mmol_1770-pqq-la_model__conditioned_La__seed-1_sample-{args.sample}'
    source=root/'workspaces/accommodation_goal_20260920/folds_v1/source_preparation/cases'/case
    summary={'case_id':case,'plan':pin(root/'diagnostics/accommodation_controls_20260920/HYDROGEN_DIAGNOSTIC_PLAN.md'),
             'source':pin(source/'normalized.pdb'),'existing_protonated':pin(source/'protonated.pdb'),
             'script':pin(__file__),'slurm_job_id':os.environ.get('SLURM_JOB_ID'),
             'slurm_array_task_id':os.environ.get('SLURM_ARRAY_TASK_ID'),'status':'running'}
    save(out/'receipt.json',summary)
    original_minimize=mm.LocalEnergyMinimizer.minimize;captured=[]
    def intercept(context,tolerance=10.,maxIterations=0,reporter=None):
        if captured:raise ValueError('unexpected second minimization call')
        if tolerance!=1. or maxIterations!=50:raise ValueError('unexpected native minimizer settings')
        system=context.getSystem();movable=np.array([system.getParticleMass(i).value_in_unit(unit.dalton)>0
                    for i in range(system.getNumParticles())])
        initial,_,_,initial_summary=state(context,movable)
        (out/'native_system.xml').write_text(mm.XmlSerializer.serialize(system))
        np.save(out/'initial_positions_nm.npy',initial);np.save(out/'movable.npy',movable)
        rep=Reporter(movable);native_started=time.monotonic()
        original_minimize(context,tolerance,maxIterations,rep)
        final,_,_,final_summary=state(context,movable)
        np.save(out/'native_final_positions_nm.npy',final)
        save(out/'native_minimizer_iterations.json',rep.rows)
        row={'initial':initial_summary,'final':final_summary,'accepted_iterations':len(rep.rows),
             'wall_seconds':time.monotonic()-native_started,
             'initial_overlap_pairs':len(overlaps(initial,movable)),
             'final_overlap_pairs':len(overlaps(final,movable)),
             'maximum_H_displacement_A':float(np.max(np.linalg.norm(final[movable]-initial[movable],axis=1)))*10,
             'fixed_atoms_exact':bool(np.array_equal(initial[~movable],final[~movable])),
             'system_constraints':system.getNumConstraints(),'added_H':int(movable.sum())}
        captured.append(row)
    try:
        sys.path.insert(0,str(root/'scripts'))
        wrapper=root/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/protonate_standard_only.py'
        spec=importlib.util.spec_from_file_location('native_protonator_probe',wrapper)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        summary['native_protonator']=pin(wrapper)
        mm.LocalEnergyMinimizer.minimize=intercept
        try:
            module.protonate_standard_only(source/'normalized.pdb',out/'replayed.pdb',
                pins_path=root/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/implementation_pins.json',ph=7.)
        finally:mm.LocalEnergyMinimizer.minimize=original_minimize
        if len(captured)!=1:raise ValueError('expected exactly one native hydrogen minimization')
        summary['native_replay']=captured[0]
        pdb=PDBFile(str(out/'replayed.pdb'));old=PDBFile(str(source/'protonated.pdb'))
        xx=np.asarray(pdb.positions.value_in_unit(unit.nanometer));yy=np.asarray(old.positions.value_in_unit(unit.nanometer))
        summary['replay_matches_existing_PDB_positions']=bool(np.array_equal(xx,yy))
        movable=np.load(out/'movable.npy');initial=np.load(out/'initial_positions_nm.npy')
        summary['native_geometry']=geometry(pdb.topology,np.load(out/'native_final_positions_nm.npy'),movable,initial)
        if args.sample==4:
            system=mm.XmlSerializer.deserialize((out/'native_system.xml').read_text())
            integ=mm.VerletIntegrator(0.)
            context=mm.Context(system,integ,mm.Platform.getPlatformByName('CPU'),{'Threads':'1'})
            try:
                summary['repair']=numerical_repair(context,initial,movable,original_minimize,out)
                final=np.load(out/'repaired_positions_nm.npy')
                summary['repair_geometry']=geometry(pdb.topology,final,movable,initial)
                with (out/'repaired.pdb').open('w') as f:PDBFile.writeFile(pdb.topology,final*unit.nanometer,f,keepIds=True)
                summary['repaired_output']=pin(out/'repaired.pdb')
            finally:del context,integ
        summary['status']='completed'
    except Exception as exc:
        summary['status']='failed';summary['exception']=str(exc);summary['traceback']=traceback.format_exc()
        if captured:summary['native_replay']=captured[0]
    finally:
        mm.LocalEnergyMinimizer.minimize=original_minimize
        summary['wall_seconds']=time.monotonic()-started
        save(out/'receipt.json',summary)
    print(json.dumps(summary,indent=2),flush=True)
    if summary['status']!='completed':raise SystemExit(1)


if __name__=='__main__':main()
