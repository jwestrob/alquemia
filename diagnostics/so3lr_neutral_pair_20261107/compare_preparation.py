"""One exact force bridge and capped preparation, no fitting or new physics."""
import os,time,argparse
from util import *
assert os.environ.get('SLURM_JOB_ID')
import jax
jax.config.update('jax_enable_x64',True);jax.config.update('jax_default_matmul_precision','highest')
import numpy as np
from ase.io import read as ase_read,write as ase_write
from ase.optimize import LBFGS
from ase.data import chemical_symbols
from scipy.spatial import cKDTree
from so3lr.ase_utils import make_ase_calculator
assert jax.devices()[0].platform=='gpu' and jax.devices('cpu')
ap=argparse.ArgumentParser();ap.add_argument('--ase-only-diagnostic',action='store_true');args=ap.parse_args()
out=W/('preparation_comparison_ASE_v1' if args.ase_only_diagnostic else 'preparation_comparison_v1');out.mkdir(exist_ok=False)
for p in [Path(__file__),D/'PREPARATION_COMPARISON.md']:(out/p.name).write_bytes(p.read_bytes())
prep=read(W/'hydrated_preparation_v1/RESULT.json')['rows'][0]
atoms=ase_read(prep['A']['path']);graphs={r['ligand_id']:r for r in read(W/'chemistry.json')['rows']}
calc=make_ase_calculator(model='so3lr-tauri-l-2606',dtype=np.float64,lr_cutoff=12.,dispersion_energy_cutoff_lr_damping=2.,
 kspace_electrostatics='pme',kspace_smearing=4.,kspace_spacing=1.,has_aux=True,observables=['partial_charges'])
assert calc.head_index==5
calls=[];original_calculate=calc.calculate
def calculate(*args,**kwargs):
 t=time.time();result=original_calculate(*args,**kwargs);calls.append(time.time()-t);return result
calc.calculate=calculate
if args.ase_only_diagnostic:
 bridge=read(W/'preparation_comparison_v1/BRIDGE.json')
 assert not bridge['passed'] and pin(bridge['raw_arrays']['path'])==bridge['raw_arrays']
 save(out/'UNQUALIFIED_NATIVE_BRIDGE.json',bridge)
 print('Independent ASE-only diagnostic; failed native bridge is preserved.',flush=True)
else:
 with np.load(W/'md_scout_v1/restart.npz') as checkpoint:
  assert int(checkpoint['step'])==100 and bool(checkpoint['fractional_coordinates'])
  box=np.asarray(checkpoint['box']);assert box.shape==(3,) and np.allclose(box,42)
  terminal=atoms.copy();terminal.positions=np.asarray(checkpoint['position'])*box;terminal.info.update(charge=0,multiplicity=6)
  native_force=np.asarray(checkpoint['force'])
 calc.calculate(terminal);ase_force=np.asarray(calc.results['forces']);energy=float(calc.results['energy'])
 assert ase_force.dtype==native_force.dtype==np.float64
 error=float(abs(ase_force-native_force).max())
 np.savez(out/'force_bridge.npz',native_forces_eV_A=native_force,ASE_forces_eV_A=ase_force,
  ASE_energy_eV=energy,positions_A=terminal.positions,Z=terminal.numbers,partial_charges=calc.results['aux']['partial_charges'])
 bridge=dict(max_force_difference_eV_A=error,passed=error<=1e-5,raw_arrays=pin(out/'force_bridge.npz'))
 save(out/'BRIDGE.json',bridge);print('Native/ASE checkpoint force difference',error,'eV/A',flush=True)
 assert bridge['passed'],'Static force interface mismatch; no optimizer released.'

def delta(x):return x-42*np.round(x/42)
def integrity(x):
 rows=[]
 for fragment in prep['fragments']:
  lid=int(fragment['ligand_id']);ids=np.array(fragment['ligand_indices']);bonds=ids[np.array(graphs[lid]['bonds'])]
  dd=np.linalg.norm(delta(x[bonds[:,0]]-x[bonds[:,1]]),axis=1)
  rows.append(dict(ligand_id=lid,bond_min_A=float(dd.min()),bond_max_A=float(dd.max()),
   bad_bonds=[dict(atom_a=int(a),atom_b=int(b),distance_A=float(dist)) for (a,b),dist in zip(bonds,dd) if dist<.7 or dist>2]))
 return rows
atoms.info.update(charge=0,multiplicity=6);atoms.calc=calc;history=[]
optimizer=LBFGS(atoms,maxstep=.05,alpha=70.,memory=100,damping=1.,use_line_search=False,
 logfile=str(out/'optimizer.log'),trajectory=str(out/'optimizer.traj'))
status='iteration_limit';previous=atoms.positions.copy()
for _ in optimizer.irun(fmax=0.,steps=200):
 forces=atoms.get_forces();e=atoms.get_potential_energy();checks=integrity(atoms.positions)
 displacement=float(np.linalg.norm(delta(atoms.positions-previous),axis=1).max())
 assert displacement<=.05000001
 row=dict(step=optimizer.nsteps,energy_eV=float(e),max_force_component_eV_A=float(abs(forces).max()),
  max_atom_step_A=displacement,ligands=checks)
 history.append(row);save(out/'PARTIAL.json',dict(rows=history))
 print('step',row['step'],'max force',row['max_force_component_eV_A'],'bond maxima',[r['bond_max_A'] for r in checks],flush=True)
 if any(r['bad_bonds'] for r in checks):status='chemical_integrity_failure';break
 if row['max_force_component_eV_A']<=.2:status='force_criterion_reached';break
 previous=atoms.positions.copy()
ase_write(out/'terminal.extxyz',atoms,format='extxyz',write_results=False)
np.savez(out/'terminal_arrays.npz',positions_A=atoms.positions,Z=atoms.numbers,energy_eV=float(atoms.get_potential_energy()),
 forces_eV_A=atoms.get_forces(),partial_charges=calc.results['aux']['partial_charges'])
save(out/'RESULT.json',dict(status=status,rows=history,force_bridge=bridge,job=os.environ['SLURM_JOB_ID'],
 diagnostic_scope='ASE-only' if args.ase_only_diagnostic else 'Native/ASE bridge and preparation',native_interface_qualified=bridge['passed'],
 model_calls=len(calls),evaluation_seconds=calls,initial=prep['A'],terminal=pin(out/'terminal.extxyz'),
 raw_terminal=pin(out/'terminal_arrays.npz'),protocol=pin(D/'PREPARATION_COMPARISON.md'),potential_unchanged=True))
print('Bounded preparation comparison:',status,flush=True)
