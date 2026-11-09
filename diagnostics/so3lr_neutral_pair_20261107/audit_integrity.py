"""Audit real covalent geometry and checkpoint mapping, without a model call."""
import os,csv,argparse
from util import *
import numpy as np
import h5py
from ase.io import read as ase_read
from ase.data import chemical_symbols
assert os.environ.get('SLURM_JOB_ID')
ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
out=args.output_dir;out.mkdir(exist_ok=False,parents=True)
for p in [Path(__file__),D/'CHEMICAL_INTEGRITY_AUDIT.md']:(out/p.name).write_bytes(p.read_bytes())
run=W/'md_scout_v1';prep=read(W/'hydrated_preparation_v1/RESULT.json')['rows'][0]
initial=ase_read(prep['A']['path']);opt=ase_read(run/'trajectory_opt.xyz',index=':')
graphs={r['ligand_id']:r for r in read(W/'chemistry.json')['rows']}
def delta(x):return x-42*np.round(x/42)
assert len(opt)==21 and all(np.array_equal(a.numbers,initial.numbers) for a in opt)
assert np.max(abs(delta(opt[0].positions-initial.positions)))<1e-7
with np.load(run/'restart.npz') as c:
 assert int(c['step'])==100 and bool(c['fractional_coordinates'])
 box=np.asarray(c['box']);assert box.shape==(3,) and np.allclose(box,42)
 terminal=np.asarray(c['position'])*box
 force=np.asarray(c['force']);assert force.dtype==np.float64 and np.isfinite(force).all()
with h5py.File(run/'trajectory.hdf5','r') as h:
 assert h['positions'].shape==(100,7288,3) and not h.attrs['fractional_coordinates']
 assert np.array_equal(h.attrs['atomic_numbers'],initial.numbers)
 assert np.max(abs(delta(terminal-h['positions'][-1])))<1e-10
 first_MD=np.asarray(h['positions'][0]);last_MD=np.asarray(h['positions'][-1])
 max_entry_displacement=float(np.max(np.linalg.norm(delta(first_MD-opt[-1].positions),axis=1)))
 frames=[('initial',0,initial.positions)]+[('relaxation',i*10,a.positions) for i,a in enumerate(opt)]+[('MD',(i+1)*10,np.asarray(a)) for i,a in enumerate(h['positions'])]
 rows=[];first_bad={}
 for phase,step,x in frames:
  for fragment in prep['fragments']:
   lid=int(fragment['ligand_id']);inds=np.array(fragment['ligand_indices']);bonds=inds[np.array(graphs[lid]['bonds'])]
   dd=np.linalg.norm(delta(x[bonds[:,0]]-x[bonds[:,1]]),axis=1)
   for k,(a,b) in enumerate(bonds):
    dist=float(dd[k]);key=(lid,int(a),int(b))
    if dist>2 and key not in first_bad:
     first_bad[key]=dict(ligand_id=lid,atom_a=int(a),atom_b=int(b),elements=chemical_symbols[int(initial.numbers[a])]+'-'+chemical_symbols[int(initial.numbers[b])],
       first_bad_phase=phase,first_bad_step=step,first_bad_distance_A=dist,
       original_distance_A=float(np.linalg.norm(delta(initial.positions[a]-initial.positions[b]))))
   rows.append(dict(ligand_id=lid,phase=phase,step=step,bond_max_A=float(dd.max()),bonds_over_2_A=int(sum(dd>2))))
 details=[]
 for key,r in first_bad.items():
  _,a,b=key;r=dict(r)
  r.update(optimizer_terminal_distance_A=float(np.linalg.norm(delta(opt[-1].positions[a]-opt[-1].positions[b]))),
    MD_terminal_distance_A=float(np.linalg.norm(delta(last_MD[a]-last_MD[b]))))
  details.append(r)
 payload=dict(atom_order_verified=True,initial_coordinates_verified=True,checkpoint_HDF_coordinates_verified=True,
 optimizer_to_first_MD_max_atom_displacement_A=max_entry_displacement,
 altered_bonds=details,rows=rows,job=os.environ['SLURM_JOB_ID'],model_calls=0,
 sources=[pin(run/p) for p in ['trajectory_opt.xyz','trajectory.hdf5','restart.npz']])
 save(out/'RESULT.json',payload)
 with (out/'BONDS.tsv').open('w') as h:
  fields=['ligand_id','phase','step','bond_max_A','bonds_over_2_A'];wr=csv.DictWriter(h,fieldnames=fields,delimiter='\t');wr.writeheader();wr.writerows(rows)
 print('Verified real atom/coordinate mapping. Altered bonds:',details,flush=True)
