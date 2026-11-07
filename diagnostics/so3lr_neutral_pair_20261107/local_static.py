"""Neutral first-shell controls, sampled conformers and one real force check."""
import os,time
from util import *
assert os.environ.get('SLURM_JOB_ID')
import jax
jax.config.update('jax_enable_x64',True);jax.config.update('jax_default_matmul_precision','highest')
import jax.numpy as jnp
import numpy as np
from ase.io import read as ase_read
from so3lr.base_calculator import make_so3lr
from so3lr.mlff.data.dataloader_sparse_ase import AseDataLoaderSparse
from so3lr.mlff.utils.jraph_utils import dynamically_batch_with_lr,graph_to_batch_fn
assert jax.devices()[0].platform=='cpu'
out=W/'local_v1';out.mkdir(exist_ok=False)
for p in [Path(__file__),D/'util.py',D/'PLAN.md',D/'local.sh']:(out/p.name).write_bytes(p.read_bytes())
prep=read(W/'preparation_v1/RESULT.json');frames=[];info=[]
water=ase_read(W/'preparation_v1/water.xyz').positions;water-=water[0]
for r in prep['rows']:
 a=ase_read(r['A']['path']);p=a.positions;Z=a.numbers
 O=p[88::3];dx=O[:,None,:]-p[[0,43]][None,:,:];dx-=42*np.round(dx/42)
 keep=np.flatnonzero(np.linalg.norm(dx,axis=-1).min(axis=1)<3.5)
 ii=list(range(88))+[88+3*j+k for j in keep for k in range(3)]
 z=Z[ii];x=p[ii].copy();assert len(z)+3<192
 carbs=r['fragments'][0]['carboxylate_O_indices'];near=min(carbs,key=lambda i:np.linalg.norm(x[i]-x[0]))
 u=x[near]-x[0];u/=np.linalg.norm(u);lig=r['fragments'][0]['ligand_indices']
 info.append(dict(start=r['id'],waters=len(keep),source=r['A'],radial_direction=u.tolist(),ligand_indices=lig))
 for shift in (-.1,0.,.1):
  xp=x.copy();xp[lig]+=shift*u
  for dist,ori in [(None,None),(200.,'original'),(200.,'reversed'),(400.,'original'),(400.,'reversed')]:
   xx=xp;zz=z.copy()
   if dist is not None:
    rot=np.eye(3) if ori=='original' else np.diag([-1.,1.,-1.])
    xx=np.concatenate([xp,water@rot.T+(xp[0]+xp[43])/2+[0.,dist,0.]])
    zz=np.r_[zz,[8,1,1]]
   for state in ('A','B'):
    zs=zz.copy()
    if state=='B':zs[0]=66;zs[43]=57
    cid=f's{r["id"]}_q{shift:g}_{dist}_{ori}_{state}'
    frames.append(dict(id=cid,start=r['id'],shift=shift,distance_A=dist,orientation=ori,state=state,
                       Z=zs.tolist(),positions_A=xx.tolist(),base_atoms=len(z)))
 # FD at the same first start: baseline and400-A original orientation.
 if r['id']==0:
  for shift in (-1e-4,1e-4):
   xp=x.copy();xp[lig]+=shift*u
   for dist in (None,400.):
    xx=xp if dist is None else np.concatenate([xp,water+(xp[0]+xp[43])/2+[0.,dist,0.]])
    zz=z if dist is None else np.r_[z,[8,1,1]]
    for state in ('A','B'):
     zs=zz.copy()
     if state=='B':zs[0]=66;zs[43]=57
     frames.append(dict(id=f'FD_{shift}_{dist}_{state}',start=0,shift=shift,distance_A=dist,
        orientation='original' if dist else None,state=state,Z=zs.tolist(),positions_A=xx.tolist(),base_atoms=len(z)))
assert len(frames)==98
save(out/'FRAMES.json',dict(rows=frames,preparation=pin(W/'preparation_v1/RESULT.json'),info=info,plan=pin(D/'PLAN.md')))
functions={};rows=[]
def evaluate(r,cut):
 f=out/f'{r["id"]}.xyz';xyz(f,r['Z'],r['positions_A'])
 data,stats=AseDataLoaderSparse(str(f)).load(cutoff=5,cutoff_lr=1000,calculate_neighbors_lr=True)
 g,lr=next(dynamically_batch_with_lr(data,n_node=192,n_edge=8192,n_graph=2,n_pairs=36864))
 b=graph_to_batch_fn(g,lr);b=jax.tree.map(lambda v:v.astype(jnp.float64) if hasattr(v,'dtype') and np.issubdtype(v.dtype,np.floating) else v,b)
 n=len(r['Z']);assert np.array_equal(np.array(b['positions'])[:n],r['positions_A'])
 assert float(b['total_charge'][0])==0 and float(b['num_unpaired_electrons'][0])==5 and b['residue_charge'] is None
 assert np.array_equal(np.flatnonzero(np.array(b['theory_mask'][0])),[5])
 if cut not in functions:functions[cut]=jax.jit(make_so3lr(model='so3lr-tauri-l-2606',precision=jnp.float64,lr_cutoff=cut,
   dispersion_energy_cutoff_lr_damping=None if cut is None else 2.,calculate_forces=True,output_intermediate_quantities=['partial_charges']))
 t=time.time();raw=jax.block_until_ready(functions[cut](b));arrays={k:np.asarray(v) for k,v in raw.items() if v is not None}
 assert raw['energy'].dtype==jnp.float64 and all(np.isfinite(v).all() for v in arrays.values())
 assert abs(arrays['partial_charges'][:n].sum())<1e-7
 dest=out/f'{r["id"]}_{cut}.npz';np.savez(dest,**arrays)
 row=dict(id=r['id'],cutoff_A=cut,energy_eV=float(arrays['energy'][0]),raw_arrays=pin(dest),seconds=time.time()-t)
 rows.append(row);save(out/'PARTIAL.json',dict(rows=rows));print('saved',row['id'],cut,round(row['seconds'],3),flush=True)
for r in frames:evaluate(r,1000.)
for r in frames:
 if r['shift']==0:evaluate(r,None)
assert len(rows)==128
lookup={(r['id'],r['cutoff_A']):r for r in rows}
def E(cid,cut):return lookup[cid,cut]['energy_eV']
def array(cid,cut):return np.load(lookup[cid,cut]['raw_arrays']['path'])
changes=[]
for inf in info:
 s=inf['start']
 for q in (-.1,0.,.1):
  for cut in ([1000.,None] if q==0 else [1000.]):
   ba=f's{s}_q{q:g}_None_None_A';bb=ba[:-1]+'B';bare=(E(bb,cut)-E(ba,cut))*EV_KCAL
   bf=array(bb,cut)['forces']-array(ba,cut)['forces']
   for d in (200.,400.):
    for o in ('original','reversed'):
     wa=f's{s}_q{q:g}_{d}_{o}_A';wb=wa[:-1]+'B'
     wd=(E(wb,cut)-E(wa,cut))*EV_KCAL;nb=next(r['base_atoms'] for r in frames if r['id']==wa)
     fc=(array(wb,cut)['forces'][:nb]-array(wa,cut)['forces'][:nb]-bf[:nb])*EV_KCAL
     changes.append(dict(start=s,ligand_shift_A=q,cutoff_A=cut,distance_A=d,orientation=o,
       bare_swap_kcal_mol=bare,water_effect_kcal_mol=wd-bare,max_force_difference_change_kcal_mol_A=float(abs(fc).max()),
       trans_metal_force_difference_change_kcal_mol_A=float(np.linalg.norm(fc[0]))))
def fd_delta(shift):
 b=(E(f'FD_{shift}_None_B',1000.)-E(f'FD_{shift}_None_A',1000.))*EV_KCAL
 w=(E(f'FD_{shift}_400.0_B',1000.)-E(f'FD_{shift}_400.0_A',1000.))*EV_KCAL
 return w-b
finite=(fd_delta(1e-4)-fd_delta(-1e-4))/(2e-4)
ba='s0_q0_None_None_A';bb=ba[:-1]+'B';wa='s0_q0_400.0_original_A';wb=wa[:-1]+'B'
nb=next(r['base_atoms'] for r in frames if r['id']==wa)
fc=(array(wb,1000.)['forces'][:nb]-array(wa,1000.)['forces'][:nb]-array(bb,1000.)['forces'][:nb]+array(ba,1000.)['forces'][:nb])*EV_KCAL
analytic=-float((fc[info[0]['ligand_indices']]*np.array(info[0]['radial_direction'])).sum())
assert abs(finite-analytic)<.003,(finite,analytic)
primary=[r for r in changes if r['cutoff_A']==1000.];peak=max(abs(r['water_effect_kcal_mol']) for r in primary)
save(D/'LOCAL_RESULTS.json',dict(status='pass' if peak<=.5 else 'above_neutral_extract_tolerance',max_water_effect_kcal_mol=peak,
 rows=changes,FD=dict(finite=finite,analytic=analytic,residual=finite-analytic),predictions=pin(out/'PARTIAL.json'),
 frames=pin(out/'FRAMES.json'),model_calls=len(rows),job=os.environ['SLURM_JOB_ID'],scope='Neutral first-shell extracts; full periodic context still required.'))
save(out/'RESULT.json',dict(rows=rows,checks_all_pass=True,model_calls=len(rows)))
print('NEUTRAL EXTRACT',peak,'kcal/mol; FD',finite,analytic,flush=True)
