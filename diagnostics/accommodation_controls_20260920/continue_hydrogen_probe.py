#!/usr/bin/env python3
"""One frozen continuation of an existing H-only numerical preparation."""
import argparse,json,os,time,traceback
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import unit
from openmm.app import PDBFile
from hydrogen_probe import pin,save,state,geometry,Reporter

p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
p.add_argument('--plan',type=Path,required=True)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
d={'plan':pin(a.plan),'script':pin(__file__),'parent_receipt':pin(a.source/'receipt.json'),
   'source_system':pin(a.source/'native_system.xml'),'source_positions':pin(a.source/'repaired_positions_nm.npy'),
   'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'status':'running'}
save(a.output/'receipt.json',d)
try:
    system=mm.XmlSerializer.deserialize((a.source/'native_system.xml').read_text())
    movable=np.load(a.source/'movable.npy');initial=np.load(a.source/'initial_positions_nm.npy')
    assert np.array_equal(movable,np.array([system.getParticleMass(i).value_in_unit(unit.dalton)>0 for i in range(system.getNumParticles())]))
    integrator=mm.VerletIntegrator(0.)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName('CPU'),{'Threads':'1'})
    context.setPositions(np.load(a.source/'repaired_positions_nm.npy')*unit.nanometer)
    _,_,_,d['before']=state(context,movable)
    rep=Reporter(movable);mm.LocalEnergyMinimizer.minimize(context,1.,500,rep)
    final,_,_,d['after']=state(context,movable)
    d['accepted_iterations']=len(rep.rows);save(a.output/'iterations.json',rep.rows)
    np.save(a.output/'positions_nm.npy',final)
    pdb=PDBFile(str(a.source/'repaired.pdb'));d['geometry']=geometry(pdb.topology,final,movable,initial)
    with (a.output/'continued.pdb').open('w') as f:PDBFile.writeFile(pdb.topology,final*unit.nanometer,f,keepIds=True)
    d['mobile_H_RMS_force_converged']=bool(d['after']['movable_force_RMS_kJ_mol_nm']<=1.)
    d['status']='completed';del context,integrator
except Exception as exc:
    d['status']='failed';d['exception']=str(exc);d['traceback']=traceback.format_exc()
d['wall_seconds']=time.monotonic()-started;save(a.output/'receipt.json',d);print(json.dumps(d,indent=2))
if d['status']!='completed':raise SystemExit(1)
