#!/usr/bin/env python3
"""Diagnose/recover the selected ADH9 site's hydrogen-only OpenMM objective.

No ORCA, model selection, new protonation states, force field, or heavy motion.
This is an explicitly changed numerical preparation, not the frozen wrapper.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
import openmm as mm
from openmm import app, unit

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
OLD = ROOT/'workspaces/plm_adh9_af3_20260916/preparation/pairs/PQQSEQ_13d74836d4b7a3e02140_AF3_sample1'
STEM = 'PQQSEQ_13d74836d4b7a3e02140_AF3_sample1'

def record(p):
    p=Path(p).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def write(p, value):
    with Path(p).open('x') as f: json.dump(value,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n')

def make_system(topology):
    """Exact forcefield=None objective from installed OpenMM Modeller.addHydrogens.

    All H here were newly added by the original standard-only protonation.
    Thus immobile existing atoms are precisely the non-hydrogen atoms.
    """
    system=mm.System(); nb=mm.CustomNonbondedForce('100/(r/0.1)^4')
    nb.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    nb.setCutoffDistance(1*unit.nanometer)
    bonds=mm.HarmonicBondForce();angles=mm.HarmonicAngleForce()
    for force in (nb,bonds,angles): system.addForce(force)
    atoms=list(topology.atoms());attached=[[] for a in atoms]
    for a in atoms:
        nb.addParticle([]);system.addParticle(1.0 if a.element==app.element.hydrogen else 0.0)
    for a,b in topology.bonds():
        if app.element.hydrogen in (a.element,b.element):bonds.addBond(a.index,b.index,0.1,100000.0)
        attached[a.index].append(b);attached[b.index].append(a)
    for residue in topology.residues():
        if residue.name=='HOH':raise ValueError('No waters permitted in this exact target')
        for a in residue.atoms():
            near=attached[a.index]
            if a.element==app.element.oxygen and len(near)==2 and app.element.hydrogen in [x.element for x in near]:
                angles.addAngle(near[0].index,a.index,near[1].index,1.894,460.24)
    return system

def state(context, hydrogen):
    s=context.getState(getPositions=True,getEnergy=True,getForces=True)
    xyz=np.asarray(s.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    forces=np.asarray(s.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
    e=float(s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    if not np.isfinite(xyz).all() or not np.isfinite(forces).all() or not math.isfinite(e):raise ValueError('Nonfinite OpenMM state')
    return xyz, {'energy_kJ_mol':e,'hydrogen_force_rms_kJ_mol_nm':float(np.sqrt(np.mean(forces[hydrogen]**2))),
                 'hydrogen_force_max_kJ_mol_nm':float(np.max(np.linalg.norm(forces[hydrogen],axis=1)))}

class Reporter(mm.MinimizationReporter):
    def __init__(self):super().__init__();self.entries=[]
    def report(self,iteration,x,grad,args):
        self.entries.append({'iteration':int(iteration),'args':dict(args),'gradient_norm':float(np.linalg.norm(np.asarray(grad)))})
        return False

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--iterations',type=int,default=50)
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=OLD/(STEM+'_protonated.pdb');pdb=app.PDBFile(str(source));atoms=list(pdb.topology.atoms())
    hydrogen=np.asarray([a.element==app.element.hydrogen for a in atoms]);system=make_system(pdb.topology)
    (out/'system.xml').write_text(mm.XmlSerializer.serialize(system))
    context=mm.Context(system,mm.VerletIntegrator(0.0),mm.Platform.getPlatformByName('CPU'),{'Threads':'1'})
    context.setPositions(pdb.positions);initial,initial_metrics=state(context,hydrogen);reporter=Reporter();t=time.monotonic()
    mm.LocalEnergyMinimizer.minimize(context,1.0,args.iterations,reporter)
    final,final_metrics=state(context,hydrogen);elapsed=time.monotonic()-t
    if not np.array_equal(initial[~hydrogen],final[~hydrogen]):raise ValueError('Heavy atoms moved')
    with (out/'protonated_recovered.pdb').open('x') as f:app.PDBFile.writeFile(pdb.topology,final*unit.nanometer,f,keepIds=True)
    np.save(out/'positions_nm.npy',final)
    doc={'source':record(source),'script':record(__file__),'system':record(out/'system.xml'),
         'platform':'CPU','threads':1,'tolerance':1.0,'max_iterations':args.iterations,
         'initial':initial_metrics,'final':final_metrics,'iterations_reported':len(reporter.entries),
         'heavy_displacement_A':0.0,'hydrogen_max_displacement_A':float(np.max(np.linalg.norm(final[hydrogen]-initial[hydrogen],axis=1))*10),
         'elapsed_seconds':elapsed,'trace':reporter.entries}
    write(out/'minimization.json',doc);print(json.dumps({k:v for k,v in doc.items() if k!='trace'},indent=2),flush=True)

if __name__=='__main__':main()
