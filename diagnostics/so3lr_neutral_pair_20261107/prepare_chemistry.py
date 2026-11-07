"""Export the exact ligand graph for trajectory integrity checks; no new poses."""
import os,sys
from util import *
from rdkit import Chem
assert os.environ.get('SLURM_JOB_ID')
sys.path.insert(0,str(ROOT/'diagnostics/ln_structure_increment_20261030'))
from structure import organic_species
rows=[]
for lid in (599,600):
 p=ROOT/'workspaces/ln_structure_increment_20261030/geometry'/f'{lid}.json'
 tr=read(p);m=organic_species(dict(ligand_charge=-4),Chem.MolFromSmiles(tr['smiles']))
 assert m.GetNumAtoms()==42 and Chem.GetFormalCharge(m)==-4
 centers=[]
 for i,label in Chem.FindMolChiralCenters(m,includeUnassigned=True):
  a=m.GetAtomWithIdx(i)
  if a.GetSymbol()!='C':continue
  neighbors=[a.GetIdx() for a in a.GetNeighbors()]
  assert len(neighbors)==4 and label in ('R','S')
  centers.append(dict(center=i,neighbors=neighbors,CIP_label=label))
 assert len(centers)==2
 arms=[[a,c] for _,a,c in m.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=[O])-[O-]'))]
 assert len(arms)==4
 rows.append(dict(ligand_id=lid,Z=[a.GetAtomicNum() for a in m.GetAtoms()],
  bonds=[[b.GetBeginAtomIdx(),b.GetEndAtomIdx()] for b in m.GetBonds()],
  stereocenters=centers,acetate_oxygen_pairs=arms,source=pin(p)))
p=W/'chemistry.json';assert not p.exists()
save(p,dict(rows=rows,job=os.environ['SLURM_JOB_ID'],script=pin(__file__)))
print('Verified and exported both exact 42-atom ligand graphs.')
