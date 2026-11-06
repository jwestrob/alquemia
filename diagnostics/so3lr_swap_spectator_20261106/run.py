"""The authorized fixed-composition spectator check; unchanged Tauri potentials."""
import os,sys,time,csv,hashlib
from pathlib import Path
D=Path(__file__).resolve().parent;ROOT=D.parents[1]
OLD=ROOT/'diagnostics/so3lr_tauri_screen_20261106'
sys.path.insert(0,str(OLD))
from common import read,save,pin,xyzwrite,EV_KCAL
W=ROOT/'workspaces/so3lr_swap_spectator_20261106'
assert os.environ.get('SLURM_JOB_ID')
attempt=W/'run_v1';attempt.mkdir(parents=True,exist_ok=False)
for p in (Path(__file__),D/'PLAN.md',D/'run.sh',D/'analyze.py',OLD/'common.py'):
 (attempt/p.name).write_bytes(p.read_bytes())
began=time.time();cpu_began=time.process_time()
oldw=ROOT/'workspaces/so3lr_tauri_screen_20261106'
software=read(OLD/'SOFTWARE_PINS.json')
for r in read(oldw/'software/RECEIPTS.json'):
 p=oldw/'software/so3lr'/r['path']
 assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],p
save(attempt/'START.json',dict(job=os.environ['SLURM_JOB_ID'],script=pin(__file__),plan=pin(D/'PLAN.md'),
 panel=pin(oldw/'PANEL.json'),environment=pin(OLD/'environment-lock.txt'),software=pin(OLD/'SOFTWARE_PINS.json'),
 started_epoch=began,CPU_affinity=sorted(os.sched_getaffinity(0))))

import jax
jax.config.update('jax_enable_x64',True)
jax.config.update('jax_default_matmul_precision','highest')
import jax.numpy as jnp
import numpy as np
from so3lr.base_calculator import make_so3lr
from so3lr.mlff.data.dataloader_sparse_ase import AseDataLoaderSparse
from so3lr.mlff.utils.jraph_utils import dynamically_batch_with_lr,graph_to_batch_fn
assert len(jax.devices())==1 and jax.devices()[0].platform=='cpu'

parents={r['id']:r for r in read(oldw/'PANEL.json')['rows']}
ac=parents['hydrated_acetate_La'];aq=parents['aquo8_central_La']
for cid in ('hydrated_acetate_La','aquo8_central_La'):
 assert parents[cid]['positions_A']==parents[cid.replace('La','Dy')]['positions_A']
 assert parents[cid]['Z'][1:]==parents[cid.replace('La','Dy')]['Z'][1:]
assert len(ac['Z'])==29 and len(aq['Z'])==25 and aq['Z'][1:4]==[8,1,1]
acetate=np.array(ac['positions_A'])+[-20.,0.,0.]
water_shell=np.array(aq['positions_A'])+[20.,0.,0.]
x=np.concatenate([acetate,water_shell]);zA=ac['Z']+aq['Z'];zA[29]=66
zB=zA.copy();zB[0]=66;zB[29]=57
assert sorted(zA)==sorted(zB) and zA[0]==57 and zB[0]==66
local=np.array(aq['positions_A'])[1:4]-np.array(aq['positions_A'])[1]
frames=[]
for state,z in [('A',zA),('B',zB)]:
 frames.append(dict(id=state,Z=z,positions_A=x.tolist(),charge=5,multiplicity=6,state=state,distance_A=None,orientation=None))
for distance in (100.,200.,400.):
 for orient,rotation in [('original',np.eye(3)),('reversed',np.diag([-1.,1.,-1.]))]:
  wat=local@rotation.T+[0.,distance,0.]
  for state,z in [('A',zA),('B',zB)]:
   frames.append(dict(id=f'{state}_water_{distance:g}_{orient}',Z=z+[8,1,1],positions_A=np.concatenate([x,wat]).tolist(),
                      charge=5,multiplicity=6,state=state,distance_A=distance,orientation=orient))
assert len(frames)==14
for r in frames:
 assert np.array_equal(np.array(r['positions_A'])[:54],x)
 mate=next(v for v in frames if v['state']!=r['state'] and v['distance_A']==r['distance_A'] and v['orientation']==r['orientation'])
 assert r['positions_A']==mate['positions_A'] and sorted(r['Z'])==sorted(mate['Z'])
save(attempt/'FRAMES.json',dict(rows=frames,parent_panel=pin(oldw/'PANEL.json'),plan=pin(D/'PLAN.md')))
functions={};rows=[];checks=[];timing=[]
components=['nn_energy','zbl_repulsion','electrostatic_energy','dispersion_energy']

def evaluate(model,cut,r):
 p=attempt/f"{r['id']}.extxyz"
 xyzwrite(p,r['Z'],r['positions_A'],f'Properties=species:S:1:pos:R:3 charge={r["charge"]} multiplicity={r["multiplicity"]} pbc="F F F"')
 data,stats=AseDataLoaderSparse(str(p)).load(cutoff=5.,cutoff_lr=1000.,calculate_neighbors_lr=True)
 data=list(data);assert len(data)==1
 graph0,lr0=data[0];n=len(r['Z'])
 pairs={(min(int(i),int(j)),max(int(i),int(j))) for i,j in zip(lr0['idx_i_lr'],lr0['idx_j_lr']) if i!=j and 0<=i<n and 0<=j<n}
 assert len(pairs)==n*(n-1)//2,(r['id'],len(pairs),n)
 graph,lr=next(dynamically_batch_with_lr(data,n_node=64,n_edge=4096,n_graph=2,n_pairs=4096))
 b=graph_to_batch_fn(graph,lr)
 b=jax.tree.map(lambda v:v.astype(jnp.float64) if hasattr(v,'dtype') and np.issubdtype(v.dtype,np.floating) else v,b)
 assert np.array_equal(np.array(b['positions'])[:n],np.array(r['positions_A']))
 assert np.array_equal(np.array(b['atomic_numbers'])[:n],r['Z'])
 assert float(b['total_charge'][0])==r['charge'] and float(b['num_unpaired_electrons'][0])==r['multiplicity']-1
 assert np.array_equal(np.flatnonzero(np.array(b['theory_mask'][0])),[5])
 assert b['residue_charge'] is None
 key=(model,cut)
 if key not in functions:
  functions[key]=jax.jit(make_so3lr(model=f'so3lr-tauri-{model}-2606',precision=jnp.float64,lr_cutoff=cut,
    dispersion_energy_cutoff_lr_damping=None if cut is None else 2.,calculate_forces=True,
    output_intermediate_quantities=['partial_charges']+components))
 t=time.time();raw=jax.block_until_ready(functions[key](b));elapsed=time.time()-t
 assert raw['energy'].dtype==jnp.float64 and raw['forces'].dtype==jnp.float64
 arrays={k:np.asarray(v) for k,v in raw.items() if v is not None}
 assert all(np.isfinite(v).all() for v in arrays.values())
 total=float(arrays['energy'][0]);part={k:float(arrays[k][:n].sum()) for k in components}
 closure=sum(part.values())-total;assert abs(closure)<1e-8,(r['id'],closure)
 q=arrays['partial_charges'][:n];assert abs(q.sum()-r['charge'])<1e-8
 fname=attempt/f'{model}_{r["id"]}_{cut}.npz';np.savez(fname,**arrays)
 v=dict(id=r['id'],model=model,cutoff_A=cut,energy_eV=total,component_sums_eV=part,
        charge=r['charge'],multiplicity=r['multiplicity'],natoms=n,evaluation_seconds=elapsed,raw_arrays=pin(fname),
        component_closure_eV=closure,total_predicted_charge=float(q.sum()))
 if n>=54:
  v['charge_on_acetate_complex']=float(q[:29].sum());v['charge_on_aquo_complex']=float(q[29:54].sum())
  v['charge_on_spectator']=float(q[54:].sum()) if n==57 else None
 rows.append(v);checks.append(dict(model=model,id=r['id'],cutoff_A=cut,coordinates_state_head_verified=True,complete_pairs=len(pairs),
                                float64_finite=True,charge_conserved=True,component_closure=True))
 timing.append(dict(model=model,id=r['id'],cutoff_A=cut,seconds=elapsed))
 save(attempt/'PARTIAL.json',dict(rows=rows,checks=checks,timing=timing))
 print('saved',model,r['id'],cut,round(elapsed,3),'s',flush=True)
 return v

for model in ('l','m'):
 replay=dict(aq,id='LaW8_replay')
 v=evaluate(model,1000.,replay)
 old=next(r for r in read(oldw/'panel_cpu_v1/RESULT.json')['rows'] if r['id']=='aquo8_central_La' and r['model']==model)
 raw=np.load(v['raw_arrays']['path']);de=abs(v['energy_eV']-old['energy_eV']);df=float(abs(raw['forces'][:25]-np.array(old['forces_eV_A'])).max())
 assert de<=1e-6 and df<=1e-6
 checks.append(dict(test='saved_real_LaW8_replay',model=model,energy_difference_eV=de,maximum_force_difference_eV_A=df,passed=True))
 for cut in (1000.,None):
  for r in frames:evaluate(model,cut,r)
  evaluate(model,cut,dict(id='free_spectator_water',Z=[8,1,1],positions_A=local.tolist(),charge=0,multiplicity=1))

save(attempt/'RESULT.json',dict(rows=rows,checks=checks,timing=timing,model_calls=len(rows),job=os.environ['SLURM_JOB_ID'],
 wall_seconds=time.time()-began,process_CPU_seconds=time.process_time()-cpu_began,devices=[str(d) for d in jax.devices()],
 frame_manifest=pin(attempt/'FRAMES.json'),plan=pin(D/'PLAN.md')))
print('finished',len(rows),'model calls',flush=True)
import runpy
runpy.run_path(str(D/'analyze.py'),run_name='__main__')
