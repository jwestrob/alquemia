"""One declared, geometry-only hydration repair; preserve original preparations."""
import os,time,itertools,subprocess,concurrent.futures
from util import *
import numpy as np
from scipy.spatial import cKDTree
from ase import Atoms
from ase.io import read as ase_read,write as ase_write
assert os.environ.get('SLURM_JOB_ID')
out=W/'hydrated_preparation_v1';out.mkdir(exist_ok=False)
for p in [Path(__file__),D/'HYDRATION_PREPARATION.md']:(out/p.name).write_bytes(p.read_bytes())
prep=read(W/'preparation_v1/RESULT.json');pack=W/'tools/packmol_src/packmol'
rotations=[]
for perm in itertools.permutations(range(3)):
 for signs in itertools.product([-1,1],repeat=3):
  R=np.eye(3)[list(perm)]*np.array(signs)[:,None]
  if np.linalg.det(R)<.5:continue
  for degree in range(0,90,15):
   t=np.deg2rad(degree);rz=np.array([[np.cos(t),-np.sin(t),0],[np.sin(t),np.cos(t),0],[0,0,1]])
   rotations.append(rz@R)
assert len(rotations)==144
save(out/'START.json',dict(job=os.environ['SLURM_JOB_ID'],parent=pin(W/'preparation_v1/RESULT.json'),
 script=pin(__file__),protocol=pin(D/'HYDRATION_PREPARATION.md'),packmol=pin(pack)))
def delta(x):return x-42*np.round(x/42)
def prepare(r):
 s=out/f'start{r["id"]}';s.mkdir();parent=ase_read(r['A']['path']);solute=parent[:88]
 count=3 if r['id']==1 else 2;kind='Dy' if r['id']==1 else 'La'
 source=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/external_review_20261023/data/water/fixtures'/f'aquo8_central_{kind}.json'
 g=read(source)['geometry'];source_x=np.asarray(g['xyz_A']);waters=(source_x[1:]-source_x[0]).reshape(8,3,3)
 assert g['atomic_numbers'][1:]==[8,1,1]*8
 fixed=solute.positions.copy();chosen=[]
 for fragment in r['fragments']:
  center=solute.positions[fragment['metal_index']];best=None
  for rot_id,R in enumerate(rotations):
   candidates=waters@R.T+center
   clearance=np.linalg.norm(delta(candidates[:,:,None,:]-fixed[None,None,:,:]),axis=-1).min(axis=(1,2))
   pair_clear=np.full((8,8),np.inf)
   for i,j in itertools.combinations(range(8),2):
    pair_clear[i,j]=pair_clear[j,i]=np.linalg.norm(delta(candidates[i,:,None,:]-candidates[j,None,:,:]),axis=-1).min()
   for subset in itertools.combinations(range(8),count):
    score=min(float(clearance[list(subset)].min()),*(float(pair_clear[i,j]) for i,j in itertools.combinations(subset,2)))
    if best is None or score>best[0]:best=(score,rot_id,subset,candidates[list(subset)].reshape(-1,3))
  assert best[0]>=1.6,('no plausible insertion',r['id'],fragment['ligand_id'],best[0])
  score,rot_id,subset,x=best;fixed=np.concatenate([fixed,x])
  chosen.append(dict(ligand_id=fragment['ligand_id'],source=pin(source),source_water_indices=list(subset),
    rotation=rotations[rot_id].tolist(),rotation_index=rot_id,clearance_A=score,count=count))
 nseed=2*count;zs=np.r_[solute.numbers,[8,1,1]*nseed]
 assert fixed.shape==(88+3*nseed,3)
 np.savez(s/'fixed.npz',positions_A=fixed,Z=zs)
 ase_write(s/'fixed.xyz',Atoms(numbers=zs,positions=fixed),format='xyz')
 ase_write(s/'water.xyz',Atoms(numbers=[8,1,1],positions=waters[0]-waters[0,0]),format='xyz')
 text=f'''tolerance 2.0
filetype xyz
output packed.xyz
seed {55100+r['id']}
pbc 42.0 42.0 42.0
nloop 200
structure fixed.xyz
 number 1
 fixed 0.0 0.0 0.0 0.0 0.0 0.0
end structure
structure water.xyz
 number {2400-nseed}
end structure
'''
 (s/'packmol.inp').write_text(text)
 with (s/'packmol.inp').open() as inp,(s/'packmol.out').open('w') as log:
  p=subprocess.run([str(pack)],stdin=inp,stdout=log,stderr=subprocess.STDOUT,cwd=s,timeout=120,env={**os.environ,'OMP_NUM_THREADS':'1'})
 assert p.returncode==0 and 'Success!' in (s/'packmol.out').read_text()
 a=ase_read(s/'packed.xyz');assert len(a)==7288 and np.array_equal(a.numbers[:len(zs)],zs)
 assert np.array_equal(a.numbers[88:].reshape(-1,3),np.tile([8,1,1],(2400,1)))
 assert np.max(abs(a.positions[:len(zs)]-fixed))<1e-5
 a.positions[:len(zs)]=fixed;a.set_cell([42]*3);a.set_pbc(True);a.info.update(charge=0,multiplicity=6)
 assert np.array_equal(a.positions[:88],parent.positions[:88])
 pairs=cKDTree(a.positions%42,boxsize=42).query_pairs(1.6,output_type='ndarray')
 def group(i):return 0 if i<88 else 1+(i-88)//3
 assert not any(group(i)!=group(j) for i,j in pairs)
 for j in range(2400):
  start=88+3*j;v=a.positions[start:start+3];rel=delta(v-v[0]);a.positions[start:start+3]=(v[0]%42)+rel
  assert np.all((np.linalg.norm(rel[1:],axis=1)>.85)&(np.linalg.norm(rel[1:],axis=1)<1.1))
 ase_write(s/'A.extxyz',a,format='extxyz',write_results=False)
 b=a.copy();z=b.numbers.copy();z[0]=66;z[43]=57;b.set_atomic_numbers(z);b.info.update(charge=0,multiplicity=6)
 ase_write(s/'B.extxyz',b,format='extxyz',write_results=False)
 assert np.array_equal(a.positions,b.positions) and sorted(a.numbers)==sorted(b.numbers)
 result=dict(r);result.update(A=pin(s/'A.extxyz'),B=pin(s/'B.extxyz'),seeded_waters_per_site=count,
  hydration_assignments=chosen,parent=r['A'],packing_input=pin(s/'packmol.inp'),packing_output=pin(s/'packmol.out'))
 return result
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(prepare,prep['rows']))
save(out/'RESULT.json',dict(rows=rows,job=os.environ['SLURM_JOB_ID'],molecular_calls=0,
 parent=pin(W/'preparation_v1/RESULT.json'),protocol=pin(D/'HYDRATION_PREPARATION.md')))
print('Prepared three matched hydrated pairs; seeded waters per metal:',[r['seeded_waters_per_site'] for r in rows],flush=True)
