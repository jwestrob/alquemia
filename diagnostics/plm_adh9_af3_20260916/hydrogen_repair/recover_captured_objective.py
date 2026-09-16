#!/usr/bin/env python3
"""Hydrogen-only restart using the exact System captured from the frozen helper."""
import json
import time
import numpy as np
import openmm as mm
from openmm import app,unit
import repair_hydrogens as helper

BASE=helper.ROOT/'workspaces/plm_adh9_af3_20260916/hydrogen_repair'
OUT=BASE/'recovered_exact_objective'

def main():
    OUT.mkdir(exist_ok=False,parents=True)
    source=helper.OLD/(helper.STEM+'_protonated.pdb');pdb=app.PDBFile(str(source))
    system_path=BASE/'original_capture/system_0.xml'
    system=mm.XmlSerializer.deserialize(system_path.read_text())
    hydrogen=np.asarray([a.element==app.element.hydrogen for a in pdb.topology.atoms()])
    assert system.getNumParticles()==len(hydrogen)
    assert all((system.getParticleMass(i).value_in_unit(unit.dalton)>0)==h for i,h in enumerate(hydrogen))
    ctx=mm.Context(system,mm.VerletIntegrator(0.0),mm.Platform.getPlatformByName('CPU'),{'Threads':'1'})
    ctx.setPositions(pdb.positions);initial,im=helper.state(ctx,hydrogen);reporter=helper.Reporter();t=time.monotonic()
    mm.LocalEnergyMinimizer.minimize(ctx,1.0,1000,reporter)
    final,fm=helper.state(ctx,hydrogen);elapsed=time.monotonic()-t
    assert np.array_equal(initial[~hydrogen],final[~hydrogen]),'Heavy atoms moved'
    assert fm['hydrogen_force_rms_kJ_mol_nm']<=1.0,'Hydrogen convergence failed'
    np.save(OUT/'positions_nm.npy',final)
    with (OUT/'protonated_recovered.pdb').open('x') as f:app.PDBFile.writeFile(pdb.topology,final*unit.nanometer,f,keepIds=True)
    result={'status':'PASS','source':helper.record(source),'script':helper.record(__file__),
      'helper':helper.record(helper.__file__),'system':helper.record(system_path),
      'objective':'Exact serialized System from unchanged pinned Modeller.addHydrogens(forcefield=None)',
      'deviation':'Restart from its failed, retained PDB coordinates; increase numerical iteration ceiling50→1000; no chemistry or force change',
      'seed_policy':'No new random draws; original seed20260914-derived H identities/coordinates retained as start',
      'platform':'CPU','threads':1,'tolerance_kJ_mol_nm':1.0,'max_iterations':1000,
      'iterations_reported':len(reporter.entries),'initial':im,'final':fm,'heavy_displacement_A':0.0,
      'elapsed_seconds':elapsed,'output':helper.record(OUT/'protonated_recovered.pdb'),'trace':reporter.entries}
    helper.write(OUT/'minimization.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='trace'},indent=2))

if __name__=='__main__':main()
