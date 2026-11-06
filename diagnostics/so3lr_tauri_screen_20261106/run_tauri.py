"""Frozen native Tauri panel, implementation controls and isolated-cluster timing."""
import os,time,argparse
from common import *
assert os.environ.get('SLURM_JOB_ID')
import jax
jax.config.update('jax_enable_x64',True)
jax.config.update('jax_default_matmul_precision','highest')
import jax.numpy as jnp
import numpy as np
from ase import Atoms
from so3lr.base_calculator import make_so3lr
from so3lr.ase_utils import make_ase_calculator
from so3lr.mlff.data.dataloader_sparse_ase import AseDataLoaderSparse
from so3lr.mlff.utils.jraph_utils import dynamically_batch_with_lr,graph_to_batch_fn

ap=argparse.ArgumentParser();ap.add_argument('--scout',action='store_true');ap.add_argument('--cpu',action='store_true');ap.add_argument('--attempt',default='v1');ap.add_argument('--reuse');args=ap.parse_args()
start=time.time();cpu=time.process_time();devs=jax.devices();assert len(devs)==1 and devs[0].platform==('cpu' if args.cpu else 'gpu'),devs
panel=read(W/'PANEL.json')['rows'];by={r['id']:r for r in panel}
outdir=W/(('scout' if args.scout else 'panel')+('_cpu' if args.cpu else '_gpu')+'_'+args.attempt);outdir.mkdir(exist_ok=False)
(outdir/'EXECUTED_RUNNER.py').write_bytes(Path(__file__).read_bytes())
save(outdir/'START.json',dict(script=pin(__file__),panel=pin(W/'PANEL.json'),environment=pin(D/'environment-lock.txt'),job=os.environ['SLURM_JOB_ID'],devices=[str(d) for d in devs],CPU_affinity=sorted(os.sched_getaffinity(0))))
checks=[];timing=[];functions={};call=0

def batch(r,tag,fragment=None):
 p=outdir/f'{tag}.extxyz'
 comment=f'Properties=species:S:1:pos:R:3 charge={r["charge"]} multiplicity={r["multiplicity"]} pbc="F F F"'
 if fragment:comment+=' '+fragment
 xyzwrite(p,r['Z'],r['positions_A'],comment)
 data,stats=AseDataLoaderSparse(str(p)).load(cutoff=5.,cutoff_lr=1000.,calculate_neighbors_lr=True)
 graph,lr=next(dynamically_batch_with_lr(data,n_node=32,n_edge=1024,n_graph=2,n_pairs=1024))
 b=graph_to_batch_fn(graph,lr)
 b=jax.tree.map(lambda v:v.astype(jnp.float64) if hasattr(v,'dtype') and np.issubdtype(v.dtype,np.floating) else v,b)
 n=len(r['Z']);assert np.array_equal(np.asarray(b['positions'])[:n],np.asarray(r['positions_A']))
 assert int(b['total_charge'][0])==r['charge'] and int(b['num_unpaired_electrons'][0])==r['multiplicity']-1
 assert float(b['theory_mask'][0,5])==1
 checks.append(dict(id=tag,charge=r['charge'],unpaired=int(b['num_unpaired_electrons'][0]),coordinates_exact=True))
 return b

def function(model,cut):
 key=(model,cut)
 if key not in functions:
  t=time.time();fn=make_so3lr(model='so3lr-tauri-'+model+'-2606',precision=jnp.float64,lr_cutoff=cut,
                    dispersion_energy_cutoff_lr_damping=None if cut is None else 2.,calculate_forces=True)
  functions[key]=jax.jit(fn);timing.append(dict(stage='load',model=model,cutoff=cut,wall_seconds=time.time()-t))
 return functions[key]

def evaluate(model,r,tag,cut=1000.,fragment=None):
 global call
 b=batch(r,tag)
 if fragment:b=batch(r,tag+'_fragment',fragment)
 t=time.time();raw=jax.block_until_ready(function(model,cut)(b));elapsed=time.time()-t
 assert raw['energy'].dtype==jnp.float64 and raw['forces'].dtype==jnp.float64
 assert np.isfinite(np.asarray(raw['energy'])).all() and np.isfinite(np.asarray(raw['forces'])).all()
 path=outdir/f'{model}_{tag}_{cut}.npz'
 np.savez(path,**{k:np.asarray(v) for k,v in raw.items() if v is not None})
 call+=1
 r0=dict(id=tag,model=model,cutoff_A=cut,energy_eV=float(raw['energy'][0]),forces_eV_A=np.asarray(raw['forces'])[:len(r['Z'])].tolist(),
         raw_arrays=pin(path),evaluation_wall_seconds=elapsed,charge=r['charge'],multiplicity=r['multiplicity'])
 timing.append(dict(stage='evaluate',model=model,id=tag,cutoff=cut,wall_seconds=elapsed))
 print('saved',model,tag,'seconds',round(elapsed,3),flush=True)
 return r0,b

models=['l'] if args.scout else ['l','m'];results=[]
reused=read(args.reuse)['rows'] if args.reuse else []
for model in models:
 targets=[by['aquo8_central_La'],by['aquo8_central_Dy']] if args.scout else panel
 for r in targets:
  previous=next((x for x in reused if x['model']==model and x['id']==r['id'] and x['cutoff_A']==1000.),None)
  if previous:
   assert pin(previous['raw_arrays']['path'])==previous['raw_arrays']
   result=dict(previous,reused_from=pin(args.reuse));batch(r,r['id'])
  else:result,b=evaluate(model,r,r['id'])
  results.append(result)
  save(outdir/'PARTIAL.json',dict(rows=results,checks=checks,timings=timing))
 # One intact ASE/native comparison, recreated state, using exact real coordinates.
 r=by['aquo8_central_La'];a=Atoms(numbers=r['Z'],positions=r['positions_A'],pbc=False)
 a.info.update(charge=r['charge'],multiplicity=r['multiplicity'])
 calc=make_ase_calculator(model='so3lr-tauri-'+model+'-2606',dtype=np.float64,lr_cutoff=1000.)
 a.calc=calc;t=time.time();ae=a.get_potential_energy();af=a.get_forces();elapsed=time.time()-t
 ref=next(x for x in results if x['model']==model and x['id']==r['id'])
 de=abs(ae-ref['energy_eV']);df=float(np.max(abs(af-np.array(ref['forces_eV_A']))))
 checks.append(dict(test='native_ASE',model=model,energy_diff_eV=de,max_force_diff_eV_A=df,passed=de<=1e-6 and df<=1e-6))
 save(outdir/f'ASE_COMPARISON_{model}.json',checks[-1]);assert checks[-1]['passed'],checks[-1]
 if args.scout:continue
 # Native fragment graph: translated archived acetate, explicitly split charges.
 r=json.loads(json.dumps(by['dry_acetate_La']));r['positions_A']=(np.array(r['positions_A'])+np.r_[np.zeros((1,3)),np.tile([50.,0.,0.],(7,1))]).tolist()
 frag='structure_type="dimer_translated" charge_a=3 charge_b=-1 selection_a="1" selection_b="2-8"'
 v,b=evaluate(model,r,'translated_dry_fragment',fragment=frag);results.append(v)
 assert np.array_equal(np.array(b['residue_segments'])[:8],[0,1,1,1,1,1,1,1])
 assert np.array_equal(np.array(b['residue_charge'])[0],[3,-1])
 checks.append(dict(test='fragment_fields',model=model,passed=True,residue_charge=[3,-1],residue_segments=[0,1,1,1,1,1,1,1]))
 v,_=evaluate(model,r,'translated_dry_unconstrained');results.append(v)
 # Predetermined full-precision cutoff contrasts, six states at each alternative.
 for cutoff in (100.,None):
  for cid in ['aquo8_central_La','aquo8_central_Dy','aquo9_central_La','aquo9_central_Dy','hydrated_acetate_La','hydrated_acetate_Dy']:
   v,_=evaluate(model,by[cid],cid+'_cutoff_'+str(cutoff),cut=cutoff);results.append(v)
 # Distant REAL spectator and its exact free reference geometry.
 base=by['aquo8_central_La'];wat=np.asarray(base['positions_A'])[1:4];wz=base['Z'][1:4]
 wr=dict(Z=wz,positions_A=wat.tolist(),charge=0,multiplicity=1)
 v,_=evaluate(model,wr,'spectator_water_free');results.append(v)
 for distance in (50.,100.):
  r=dict(Z=base['Z']+wz,positions_A=base['positions_A']+(wat+np.array([distance,0.,0.])).tolist(),charge=3,multiplicity=1)
  v,_=evaluate(model,r,'spectator_'+str(distance));results.append(v)
 # Real-coordinate finite difference for the collective water-shell direction.
 r0=by['aquo8_central_La'];xyz=np.array(r0['positions_A']);direction=np.zeros_like(xyz)
 for i in range(1,len(xyz),3):
  u=xyz[i]-xyz[0];u/=np.linalg.norm(u);direction[i:i+3]=u
 er=[]
 for sign in (-1,1):
  r=dict(Z=r0['Z'],positions_A=(xyz+sign*1e-4*direction).tolist(),charge=3,multiplicity=1)
  v,_=evaluate(model,r,'FD_shell_'+str(sign));results.append(v);er.append(v['energy_eV'])
 fd=(er[1]-er[0])/(2e-4);ana=-float((np.array(ref['forces_eV_A'])*direction).sum())
 checks.append(dict(test='force_energy_FD',model=model,finite_difference=fd,analytic=ana,residual_eV_A=fd-ana,passed=abs(fd-ana)<=1e-4))
 # Ten repeated warm native energy/force evaluations of the REAL AcW7 system.
 b=batch(by['hydrated_acetate_La'],'throughput_input');fn=function(model,1000.)
 jax.block_until_ready(fn(b));ts=[]
 for _ in range(10):
  t=time.perf_counter();jax.block_until_ready(fn(b));ts.append(time.perf_counter()-t)
 timing.append(dict(stage='warm_throughput',model=model,natoms=29,repeats=10,seconds=ts,mean_seconds=float(np.mean(ts))))
save(outdir/'RESULT.json',dict(rows=results,checks=checks,timings=timing,model_calls=call,
 devices=[str(d) for d in devs],jax_version=jax.__version__,job_id=os.environ['SLURM_JOB_ID'],wall_seconds=time.time()-start,process_CPU_seconds=time.process_time()-cpu))
print('completed',len(results),'saved evaluations',flush=True)
