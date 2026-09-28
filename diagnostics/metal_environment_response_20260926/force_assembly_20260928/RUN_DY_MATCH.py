"""Exactly two authorized original-geometry classical evaluations, no quantum."""
from pathlib import Path
import os,sys,time,json,hashlib
import numpy as np
import openmm as mm
from openmm import unit
W=Path(sys.argv[1]).resolve();sys.path.insert(0,str(W/'implementation'))
from affordable_common import read_json,verify,record,write_new
from metal_environment_component_checks import build_system,GROUPS
if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Shared CPU allocation required')
m=read_json(W/'manifest.json');[verify(x) for x in m['implementation'].values()];ledger=read_json(verify(m['ledger']));artifacts={k:read_json(verify(v)) for k,v in ledger['artifacts'].items()};start=time.monotonic()
system,description=build_system(ledger,artifacts,'Dy');integ=mm.VerletIntegrator(.001);context=mm.Context(system,integ,mm.Platform.getPlatformByName('Reference'));calls=0;rows={}
for config in ('A','B'):
 inp=read_json(verify(m['configurations'][config]));d=W/config;d.mkdir(exist_ok=False);context.setPositions(np.array(inp['coordinates_A'])*unit.angstrom);components={};gradients=[];energies={}
 for gi,name in enumerate(GROUPS):
  st=context.getState(getEnergy=True,getForces=True,groups={gi});calls+=1;e=st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole);g=-np.array(st.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom));assert np.isfinite(g).all() and np.isfinite(e)
  p=d/(name+'.json');write_new(p,dict(component=name,energy_kcal_mol=e,gradient_kcal_mol_A=g.tolist(),source_ids=inp['source_ids'],coordinates_A=inp['coordinates_A'],quantity='gradient',units='kcal/mol/angstrom'));components[name]=record(p);energies[name]=e;gradients.append(g)
 out=dict(source_ids=inp['source_ids'],coordinates_A=inp['coordinates_A'],gradient_kcal_mol_A=sum(gradients).tolist(),metal='Dy',quantity='gradient',units='kcal/mol/angstrom',QM_MM_Coulomb_included=False,C4_induction_included=False,ledger=m['ledger'],components=components,energies_kcal_mol=energies,full_hybrid_qualified=False,solvent_included=False,physical_id_alias=inp['physical_id_alias'],preparation=inp['preparation']);write_new(d/'CLASSICAL.json',out);rows[config]=record(d/'CLASSICAL.json')
del context,integ
write_new(W/'RECEIPT.json',dict(job_id=os.environ['SLURM_JOB_ID'],wall_seconds=time.monotonic()-start,allocated_CPUs=int(os.environ['SLURM_CPUS_ON_NODE']),platform='OpenMM Reference',openmm_version=mm.__version__,component_energy_force_queries=calls,configurations=2,quantum_calls=0,manifest=record(W/'manifest.json'),rows=rows,description=description))
