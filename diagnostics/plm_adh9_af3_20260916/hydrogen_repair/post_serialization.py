#!/usr/bin/env python3
"""Measure forces after normal PDB serialization, without further minimization."""
import numpy as np
import openmm as mm
from openmm import app,unit
import repair_hydrogens as helper

BASE=helper.ROOT/'workspaces/plm_adh9_af3_20260916/hydrogen_repair'

def main():
    path=BASE/'recovered_exact_objective/protonated_recovered.pdb';pdb=app.PDBFile(str(path))
    sys=mm.XmlSerializer.deserialize((BASE/'original_capture/system_0.xml').read_text())
    h=np.asarray([a.element==app.element.hydrogen for a in pdb.topology.atoms()])
    ctx=mm.Context(sys,mm.VerletIntegrator(0),mm.Platform.getPlatformByName('CPU'),{'Threads':'1'})
    ctx.setPositions(pdb.positions);xyz,metrics=helper.state(ctx,h)
    original=np.load(BASE/'recovered_exact_objective/positions_nm.npy')
    out={'status':'MEASURED','pdb':helper.record(path),'script':helper.record(__file__),'post_PDB_serialization':metrics,
      'maximum_H_rounding_displacement_A':float(np.max(np.linalg.norm(xyz[h]-original[h],axis=1))*10),
      'maximum_heavy_rounding_displacement_A':float(np.max(np.linalg.norm(xyz[~h]-original[~h],axis=1))*10),
      'minimization_performed':False,'interpretation':'Convergence criterion applies to full-precision minimizer output. These forces quantify normal three-decimal-Angstrom PDB serialization; geometry is separately checked on actual carved coordinates.'}
    helper.write(BASE/'recovered_exact_objective/post_serialization.json',out);print(out)

if __name__=='__main__':main()
