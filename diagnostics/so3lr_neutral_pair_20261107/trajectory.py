"""Analyze real initial endpoints or native scout frames without a model load."""
import os,argparse,csv,re,time
from util import *
import numpy as np
from scipy.spatial import cKDTree
from ase.io import read as ase_read
assert os.environ.get('SLURM_JOB_ID')
ap=argparse.ArgumentParser();ap.add_argument('--md-dir',type=Path);ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--preparation',type=Path,default=W/'preparation_v1/RESULT.json')
args=ap.parse_args();out=args.output_dir;out.mkdir(exist_ok=False,parents=True)
for p in [Path(__file__),D/'OBSERVABLES.md',D/'util.py']:(out/p.name).write_bytes(p.read_bytes())
prep=read(args.preparation);chem={r['ligand_id']:r for r in read(W/'chemistry.json')['rows']}
rows=[];events=[];log_rows=[];sources=[pin(args.preparation),pin(W/'chemistry.json')]
def displacement(x):return x-42*np.round(x/42)
def volume(x,ids):
 p=x[ids];return float(np.linalg.det(p[:3]-p[3]))
def analyze(x,z,preparation,state,phase,frame,t_ps):
 assert x.shape==(7288,3) and np.isfinite(x).all()
 assert np.array_equal(z[88:].reshape(-1,3),np.tile([8,1,1],(2400,1)))
 fragments=preparation['fragments'];s=preparation['id']
 original=ase_read(preparation[state]['path']).positions
 oxygen=np.arange(88,len(z),3);hydrogen=np.ravel(np.column_stack([oxygen+1,oxygen+2]))
 parent=np.repeat(oxygen,2);oh=np.linalg.norm(displacement(x[hydrogen]-x[parent]),axis=1)
 hetero=np.flatnonzero(np.isin(z,[7,8]));tree=cKDTree(x[hetero]%42,boxsize=42)
 dd,jj=tree.query(x[hydrogen]%42)
 moved=np.flatnonzero((oh>1.35)&(dd<1.25)&(hetero[jj]!=parent))
 if len(moved):events.append(dict(start=s,state=state,phase=phase,frame=frame,
  water_proton_transfer_candidates=[dict(H=int(hydrogen[j]),original_O=int(parent[j]),new_ON=int(hetero[jj[j]])) for j in moved]))
 for f in fragments:
  m=f['metal_index'];lig=np.array(f['ligand_indices']);graph=chem[int(f['ligand_id'])]
  assert np.array_equal(z[lig],graph['Z'])
  bonds=lig[np.array(graph['bonds'])];bl=np.linalg.norm(displacement(x[bonds[:,0]]-x[bonds[:,1]]),axis=1)
  on=f['potential_ON_indices'];od=np.linalg.norm(displacement(x[on]-x[m]),axis=1)
  wd=np.linalg.norm(displacement(x[oxygen]-x[m]),axis=1)
  unwrapped=x[m]+displacement(x[on]-x[m]);centroid=np.linalg.norm(unwrapped.mean(axis=0)-x[m])
  flips=[];volumes=[]
  for center in graph['stereocenters']:
   ids=lig[center['neighbors']]
   xv=displacement(x[ids]-x[lig[center['center']]])
   ref=displacement(original[ids]-original[lig[center['center']]])
   v=volume(xv,np.arange(4));v0=volume(ref,np.arange(4));assert abs(v0)>.1
   volumes.append(v)
   if v*v0<=0:flips.append(center['center'])
  contacts={str(cut):dict(ligand_ON=int(sum(od<cut)),water_O=int(sum(wd<cut)),
    acetate_arms=int(sum(min(np.linalg.norm(displacement(x[lig[pair]]-x[m]),axis=1))<cut for pair in graph['acetate_oxygen_pairs']))) for cut in [3.1,3.3,3.5]}
  kdist=np.linalg.norm(displacement(x[[86,87]]-x[m]),axis=1)
  ko=np.linalg.norm(displacement(x[[86,87],None,:]-x[f['carboxylate_O_indices']][None,:,:]),axis=-1).min(axis=1)
  rows.append(dict(start=s,state=state,phase=phase,frame=frame,time_ps=t_ps,
    ligand_id=f['ligand_id'],metal_atom_index=m,metal_atomic_number=int(z[m]),
    ligand_ON_atom_indices=on,ligand_ON_distances_A=od.tolist(),contacts=contacts,
    donor_centroid_distance_A=float(centroid),water_O_indices_within3p3=oxygen[wd<3.3].tolist(),
    bond_min_A=float(bl.min()),bond_max_A=float(bl.max()),ligand_bonds_outside_range=int(sum((bl<.7)|(bl>2))),
    stereo_signed_volumes_A3=volumes,stereo_flipped_centers=flips,
    K_metal_distances_A=kdist.tolist(),K_ligand_O_min_distances_A=ko.tolist(),
    water_OH_min_A=float(oh.min()),water_OH_max_A=float(oh.max()),water_OH_outside_range=int(sum((oh<.7)|(oh>1.35))),
    water_proton_transfer_candidate_count=len(moved)))

for r in prep['rows']:
 for state in ['A','B']:
  a=ase_read(r[state]['path']);assert a.info['charge']==0 and a.info['multiplicity']==6
  assert np.allclose(a.cell.array,np.eye(3)*42) and a.pbc.all()
  analyze(a.positions,a.numbers,r,state,'initial',0,0.)

if args.md_dir:
 import h5py
 p=args.md_dir/'trajectory.hdf5';assert p.exists()
 sources.append(pin(p))
 r=prep['rows'][0];a=ase_read(r['A']['path'])
 opt=args.md_dir/'trajectory_opt.xyz'
 if opt.exists():
  sources.append(pin(opt))
  for i,frame in enumerate(ase_read(opt,index=':')):analyze(frame.positions,frame.numbers,r,'A','relaxation',i,None)
 with h5py.File(p,'r') as h:
  assert h['positions'].dtype==np.float64 and not bool(h.attrs['fractional_coordinates'])
  assert np.array_equal(h.attrs['atomic_numbers'],a.numbers)
  assert Path(str(h.attrs['input_file'])).resolve()==Path(r['A']['path']).resolve()
  assert h.attrs['trajectory_layout']=='frame-major-v1'
  assert abs(float(h.attrs['md_dt_fs'])-.25)<1e-12 and int(h.attrs['md_steps'])==10
  assert np.allclose(h.attrs['cell_a0'],np.eye(3)*42)
  for i,x in enumerate(h['positions']):
   box=np.asarray(h['box'][i]);box=np.diag(box) if box.shape==(3,) else box
   assert box.shape==(3,3) and np.allclose(box,np.eye(3)*42)
   q=np.asarray(h['partial_charges'][i]);assert q.dtype==np.float64 and np.isfinite(q).all() and abs(q.sum())<1e-6
   analyze(np.asarray(x),a.numbers,r,'A','MD',i+1,(i+1)*.0025)
 log=args.md_dir/'native.log'
 if log.exists():
  sources.append(pin(log))
  pattern=r'(?:^|\s)(\d+)\t([-+\deE.]+)\t([-+\deE.]+)\t([-+\deE.]+)\t([-+\deE.]+)\t([-+\deE.]+)\t([-+\deE.]+)'
  for line in log.read_text(errors='replace').splitlines():
   match=re.search(pattern,line)
   if match:
    values=list(map(float,match.groups()))
    assert np.isfinite(values).all()
    log_rows.append(dict(zip(['step','rounded_total_E_eV','rounded_KE_eV','rounded_PE_eV','rounded_H_eV','temperature_K','seconds_per_step'],values)))
save(out/'RESULT.json',dict(rows=rows,proton_transfer_events=events,sources=sources,
 native_log_rows=log_rows,log_energies_for_thermodynamic_comparison=False,
 analysis_job=os.environ['SLURM_JOB_ID'],molecular_calls=0,synthetic_fixtures=0,
 unrestrained_diagnostics_only=True))
with (out/'CONTACTS.tsv').open('w') as h:
 fields=['start','state','phase','frame','time_ps','ligand_id','metal_atomic_number','donor_centroid_distance_A','contacts','bond_min_A','bond_max_A','stereo_flipped_centers']
 writer=csv.DictWriter(h,fieldnames=fields,delimiter='\t');writer.writeheader()
 for r in rows:writer.writerow({k:r[k] for k in fields})
print('Analyzed',len(rows),'real ligand/frame records; no model calls.',flush=True)
