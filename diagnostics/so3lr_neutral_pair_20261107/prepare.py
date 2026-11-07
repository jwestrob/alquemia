"""Reuse stereochemical ligand poses; build three neutral periodic preparations."""
import os,sys,json,time,hashlib,subprocess,concurrent.futures
from pathlib import Path
import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree
from ase import Atoms
from ase.io import read as ase_read,write as ase_write
D=Path(__file__).resolve().parent;ROOT=D.parents[1];W=ROOT/'workspaces/so3lr_neutral_pair_20261107'
OLD=ROOT/'diagnostics/ln_structure_increment_20261030';sys.path.insert(0,str(OLD))
from structure import organic_species
assert os.environ.get('SLURM_JOB_ID')
def pin(p):return dict(path=str(Path(p).resolve()),sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest())
def save(p,j):Path(p).write_text(json.dumps(j,indent=2,allow_nan=False)+'\n')
import csv
oldw=ROOT/'workspaces/ln_structure_increment_20261030'
rows={r['lig_id']:r for r in csv.DictReader((oldw/'primary.tsv').open(),delimiter='\t') if r['lig_id'] in ('599','600')}
out=W/'preparation_v1';out.mkdir(exist_ok=False)
for p in [Path(__file__),D/'PLAN.md',D/'prepare.sh']:(out/p.name).write_bytes(p.read_bytes())
pack=W/'tools/packmol_src/packmol'
assert subprocess.check_output(['git','-C',str(pack.parent),'rev-parse','HEAD'],text=True).strip()=='7ddbd6562301aebf2a4980599ebb3ee4191020cd'
decl=dict(plan=pin(D/'PLAN.md'),script=pin(__file__),job=os.environ['SLURM_JOB_ID'],packmol=pin(pack),
 sources=[pin(oldw/'geometry'/f'{lid}.json') for lid in ['599','600']],identity_source=pin(oldw/'primary.tsv'))
save(out/'START.json',decl)
models={};traces={}
for lid in ['600','599']:
 trace=json.loads((oldw/'geometry'/f'{lid}.json').read_text());m=organic_species(rows[lid],Chem.MolFromSmiles(trace['smiles']))
 assert Chem.GetFormalCharge(m)==-4 and m.GetNumAtoms()==42
 assert sum(a.GetSymbol()=='K' for a in m.GetAtoms())==0
 models[lid]=m;traces[lid]=trace
water_source=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/external_review_20261023/data/water/fixtures/aquo8_central_La.json'
g=json.loads(water_source.read_text())['geometry'];water=np.array(g['xyz_A'])[1:4];water-=water[0]
assert g['atomic_numbers'][1:4]==[8,1,1]
ase_write(out/'water.xyz',Atoms(numbers=[8,1,1],positions=water),format='xyz')

def rotate(axis,angle):
 a=np.asarray(axis,dtype=float);a/=np.linalg.norm(a);t=np.deg2rad(angle)
 K=np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])
 return np.eye(3)*np.cos(t)+(1-np.cos(t))*np.outer(a,a)+np.sin(t)*K

prepared=[]
for seed,(cid,metal) in enumerate([(0,'La'),(1,'Dy'),(2,'La')]):
 s=out/f'start{seed}';s.mkdir();zs=[];coords=[];chem=[];fragments=[]
 for site,lid in enumerate(['600','599']):
  m=models[lid];tr=traces[lid];c=next(c for c in tr['candidates'] if c['k']==20 and c['conformer']==cid and c['assignment']==0 and c['metal']==metal)
  assert c['converged'];x=np.array(c['coordinates']);assert len(x)==m.GetNumAtoms()
  actual=Chem.Mol(m);conf=Chem.Conformer(len(x))
  for i,p in enumerate(x):conf.SetAtomPosition(i,p)
  actual.RemoveAllConformers();actual.AddConformer(conf)
  expected=[(i,code) for i,code in Chem.FindMolChiralCenters(m,includeUnassigned=True) if m.GetAtomWithIdx(i).GetSymbol()=='C']
  Chem.RemoveStereochemistry(actual);Chem.AssignAtomChiralTagsFromStructure(actual,confId=0,replaceExistingTags=True);Chem.AssignStereochemistry(actual,cleanIt=True,force=True)
  observed=[(i,code) for i,code in Chem.FindMolChiralCenters(actual,includeUnassigned=True) if actual.GetAtomWithIdx(i).GetSymbol()=='C']
  mirror=[(i,{'R':'S','S':'R'}.get(code,code)) for i,code in expected]
  assert observed==expected or observed==mirror,(lid,cid,expected,observed)
  bond_lengths=[]
  for b in m.GetBonds():
   i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx();dist=float(np.linalg.norm(x[i]-x[j]));bond_lengths.append(dist)
   assert .7<dist<2.,(lid,i,j,dist)
  R=rotate([0,0,1],seed*120) if site==0 else rotate([0,1,0],seed*90)
  center=np.array([11.+site*20,21.,21.]);x=(x-np.array(c['center']))@R.T+center
  offset=len(zs);zs+=[57 if site==0 else 66]+[a.GetAtomicNum() for a in m.GetAtoms()];coords.append(np.vstack([center,x]))
  fragments.append(dict(ligand_id=lid,metal_index=offset,ligand_indices=list(range(offset+1,offset+43)),
    donor_indices=[offset+1+i for i in tr['potential_donors']],potential_ON_indices=[offset+1+a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() in ('N','O')],
    carboxylate_O_indices=[offset+1+a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='O'],
    original_candidate={k:c[k] for k in ['conformer','assignment','k','metal','converged']},formal_ligand_charge=-4))
  chem.append(dict(ligand=lid,expected=expected,observed=observed,relative_stereo_verified=True,bond_lengths_A=bond_lengths))
 zs+=[19,19];coords.append(np.array([[21.,10.,21.],[21.,32.,21.]]));x=np.concatenate(coords)
 assert len(zs)==88 and ((sum(zs)-5)%2)==0 and np.all((x>0)&(x<42))
 atoms=Atoms(numbers=zs,positions=x,cell=[42]*3,pbc=True);atoms.info.update(charge=0,multiplicity=6)
 # Preserve exact source coordinates separately from the packing file export.
 np.savez(s/'exact_solute.npz',Z=zs,positions_A=x)
 ase_write(s/'solute.xyz',atoms,format='xyz')
 inp=f'''tolerance 2.0
filetype xyz
output packed.xyz
seed {44100+seed}
pbc 42.0 42.0 42.0
nloop 200
structure solute.xyz
 number 1
 fixed 0.0 0.0 0.0 0.0 0.0 0.0
end structure
structure ../water.xyz
 number 2400
end structure
'''
 (s/'packmol.inp').write_text(inp)
 prepared.append(dict(id=seed,directory=str(s),fragments=fragments,chemistry_checks=chem,Z=zs,solute_positions_A=x.tolist()))

def pack_one(r):
 s=Path(r['directory']);t=time.time()
 with (s/'packmol.inp').open() as inp,(s/'packmol.out').open('w') as log:
  proc=subprocess.run([str(pack)],stdin=inp,stdout=log,stderr=subprocess.STDOUT,cwd=s,timeout=480,env={**os.environ,'OMP_NUM_THREADS':'1'})
 assert proc.returncode==0 and 'Success!' in (s/'packmol.out').read_text(),(r['id'],proc.returncode)
 a=ase_read(s/'packed.xyz');assert len(a)==88+2400*3
 assert a.numbers[:88].tolist()==r['Z'] and a.numbers[88:].reshape(-1,3).tolist()==[[8,1,1]]*2400
 original=np.array(r['solute_positions_A']);err=float(abs(a.positions[:88]-original).max());assert err<1e-5,err
 a.positions[:88]=original
 a.set_cell([42.]*3);a.set_pbc(True);a.info.update(charge=0,multiplicity=6)
 xyz=a.positions.copy();wrapped=xyz%42
 pairs=cKDTree(wrapped,boxsize=42.).query_pairs(1.8,output_type='ndarray')
 def group(i):return 0 if i<88 else 1+(i-88)//3
 invalid=[(int(i),int(j)) for i,j in pairs if group(i)!=group(j)]
 assert not invalid,(r['id'],'periodic clash',invalid[:5])
 for j in range(2400):
  w=a.positions[88+j*3:91+j*3];rel=w-w[0];rel-=42*np.round(rel/42)
  assert np.max(abs(np.sort(np.linalg.norm(rel[1:],axis=1))-np.sort(np.linalg.norm(water[1:],axis=1))))<1e-5
 # Rewrap each water rigidly so its two O-H vectors remain whole.
  a.positions[88+j*3:91+j*3]=(w[0]%42)+rel
 path=s/'A.extxyz';ase_write(path,a,format='extxyz',write_results=False)
 b=a.copy();zb=b.numbers.copy();zb[0]=66;zb[43]=57;b.set_atomic_numbers(zb);b.info.update(charge=0,multiplicity=6)
 ase_write(s/'B.extxyz',b,format='extxyz',write_results=False)
 assert sorted(a.numbers.tolist())==sorted(b.numbers.tolist()) and np.array_equal(a.positions,b.positions)
 masses=float(a.get_masses().sum());density=masses*1.66053906660/(42**3)
 assert .9<density<1.1
 return dict(id=r['id'],A=pin(path),B=pin(s/'B.extxyz'),atoms=len(a),water_molecules=2400,cell_A=42,
   charge=0,multiplicity=6,density_g_mL=density,packing_seconds=time.time()-t,solute_rounding_restored_max_A=err,
   periodic_clashes_below_1p8_A=0,fragments=r['fragments'],chemistry_checks=r['chemistry_checks'],
   counterion_indices=[86,87],water_start=88,packmol_input=pin(s/'packmol.inp'),packmol_output=pin(s/'packmol.out'))

with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:result=list(pool.map(pack_one,prepared))
save(out/'RESULT.json',dict(rows=result,start=decl,water_source=pin(water_source),preparation_only=True,ML_calls=0,DFT_calls=0))
print('Prepared three chemically checked neutral periodic systems:',[r['atoms'] for r in result],flush=True)
