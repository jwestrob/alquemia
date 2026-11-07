"""Full neutral periodic endpoint controls and measured real-box throughput."""
import os,time,argparse
from util import *
assert os.environ.get('SLURM_JOB_ID')
ap=argparse.ArgumentParser();ap.add_argument('--cpu-scout',action='store_true');ap.add_argument('--cpu-controls',action='store_true');ap.add_argument('--attempt',default='v1');args=ap.parse_args()
assert not (args.cpu_scout and args.cpu_controls)
import jax
jax.config.update('jax_enable_x64',True);jax.config.update('jax_default_matmul_precision','highest')
import numpy as np
from ase.io import read as ase_read
from so3lr.ase_utils import make_ase_calculator
from glp import atoms_to_system
cpu_mode=args.cpu_scout or args.cpu_controls
assert jax.devices()[0].platform==('cpu' if cpu_mode else 'gpu')
tag=('periodic_cpu_controls_' if args.cpu_controls else 'periodic_cpu_' if args.cpu_scout else 'periodic_gpu_')+args.attempt
out=W/tag;out.mkdir(exist_ok=False);started=time.time();cput=time.process_time()
for p in [Path(__file__),D/'PLAN.md',D/'util.py']:(out/p.name).write_bytes(p.read_bytes())
save(out/'START.json',dict(job=os.environ['SLURM_JOB_ID'],script=pin(__file__),plan=pin(D/'PLAN.md'),
 devices=[str(d) for d in jax.devices()],preparation=pin(W/'preparation_v1/RESULT.json')))
prep=read(W/'preparation_v1/RESULT.json');calcs={};rows=[];selections=[]
def evaluate(atoms,cid):
 n=len(atoms);assert atoms.info['charge']==0 and atoms.info['multiplicity']==6
 system=atoms_to_system(atoms,dtype=np.float64)
 assert float(system.total_charge[0])==0 and float(system.num_unpaired_electrons[0])==5
 if n not in calcs:
  calcs[n]=make_ase_calculator(model='so3lr-tauri-l-2606',dtype=np.float64,lr_cutoff=12.,
    dispersion_energy_cutoff_lr_damping=2.,kspace_electrostatics='pme',kspace_smearing=4.,kspace_spacing=1.,
    has_aux=True,observables=['partial_charges'])
 calc=calcs[n];assert calc.head_index==5
 t=time.time();calc.calculate(atoms);elapsed=time.time()-t
 energy=float(calc.results['energy']);forces=np.asarray(calc.results['forces'])
 q=np.asarray(calc.results['aux']['partial_charges'])
 assert forces.dtype==np.float64 and np.isfinite(energy) and np.isfinite(forces).all() and abs(q.sum())<1e-6
 p=out/f'{cid}.npz';np.savez(p,energy_eV=energy,forces_eV_A=forces,partial_charges=q,positions_A=atoms.positions,Z=atoms.numbers)
 row=dict(id=cid,energy_eV=energy,evaluation_seconds=elapsed,natoms=n,raw_arrays=pin(p),
  max_force_eV_A=float(abs(forces).max()),total_partial_charge=float(q.sum()),state_head_verified=True)
 rows.append(row);save(out/'PARTIAL.json',dict(rows=rows));print('saved',cid,elapsed,'s',flush=True)
 return row

for r in (prep['rows'][:1] if cpu_mode else prep['rows']):
 A=ase_read(r['A']['path']);B=ase_read(r['B']['path']);s=r['id']
 assert len(A)==7288 and np.array_equal(A.positions,B.positions) and sorted(A.numbers)==sorted(B.numbers)
 if args.cpu_controls:
  # Reuse the two completed, exact CPU endpoints; no repeat of those calls.
  prior=read(W/'periodic_cpu_v2/RESULT.json')
  for state in ('A','B'):
   row=next(row for row in prior['rows'] if row['id']==f's{s}_{state}')
   assert pin(row['raw_arrays']['path'])==row['raw_arrays']
   rows.append(dict(row,reused_from_job=prior['job']))
 else:
  evaluate(A,f's{s}_A');evaluate(B,f's{s}_B')
 if args.cpu_scout:continue
 # An arbitrary common translation plus wrapping represents the same system.
 for a,state in [(A,'A'),(B,'B')]:
  t=a.copy();t.positions=(t.positions+[.371,.533,.797])%42.;t.info.update(charge=0,multiplicity=6)
  evaluate(t,f's{s}_translated_{state}')
 # Source-only geometric selection, shared by both endpoints, before inference.
 O=A.positions[88::3];heavy=np.flatnonzero((A.numbers[:88]>1))
 delta=O[:,None,:]-A.positions[heavy][None,:,:];delta-=42*np.round(delta/42)
 clearance=np.linalg.norm(delta,axis=-1).min(axis=1)
 selected=[]
 for j in np.argsort(-clearance,kind='stable'):
  if all(np.linalg.norm((O[j]-O[k])-42*np.round((O[j]-O[k])/42))>=6 for k in selected):selected.append(int(j))
  if len(selected)==3:break
 assert len(selected)==3
 selections.append(dict(start=s,water_ids=selected,solute_clearances_A=[float(clearance[j]) for j in selected]))
 save(out/'REMOTE_SELECTIONS.json',selections)
 for j in selected:
  indices=np.array([i for i in range(len(A)) if i not in [88+3*j,89+3*j,90+3*j]])
  for a,state in [(A,'A'),(B,'B')]:
   reduced=a[indices];reduced.info.update(charge=0,multiplicity=6)
   evaluate(reduced,f's{s}_remove{j}_{state}')

timings=[]
if not cpu_mode:
 A=ase_read(prep['rows'][0]['A']['path'])
 for k in range(5):
  # calculate() explicitly bypasses ASE's unchanged-coordinate cache.
  t=time.time();calcs[len(A)].calculate(A);timings.append(time.time()-t)
save(out/'RESULT.json',dict(rows=rows,remote_selections=selections,warm_seconds=timings,
 elapsed_seconds=time.time()-started,process_CPU_seconds=time.process_time()-cput,job=os.environ['SLURM_JOB_ID'],
 metadata=dict(model='L-2606',precision='float64',charge=0,multiplicity=6,PME=True,mesh_spacing_A=1.,real_cutoff_A=12.)))
print('complete',tag,len(rows),'rows; warm seconds',timings,flush=True)
