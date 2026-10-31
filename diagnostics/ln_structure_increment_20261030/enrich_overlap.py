"""Attach released experimental conditions and other-metal/protonation context.

Metadata join by verified source row; no additional predictions or model fits.
"""
import json,os
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from curate import W,read,write,identity
assert os.environ.get('SLURM_JOB_ID')
t=pd.read_csv(W/'released_training_inventory.tsv',sep='\t',keep_default_na=False)
condition_tables={}
for source in t.source.unique():
 if source.startswith('data/'):
  path=W/'upstream'/source
  fpath=path.with_name(path.name.replace('_input.csv','_features.csv'))
  if fpath.exists():
   f=pd.read_csv(fpath)
   n=len(pd.read_csv(path))
   assert len(f)==n,(path,len(f),n)
   condition_tables[source]=f
metals=[];temp=[];ionic=[];charge=[];electrolyte=[]
cache={}
for r in t.itertuples(index=False):
 smi=r.complex_smiles
 if smi not in cache:
  try:
   m=Chem.MolFromSmiles(smi)
   mm=[a for a in m.GetAtoms() if a.GetAtomicNum() not in (1,6,7,8,9,15,16,17,35,53)]
   if not mm: mm=[a for a in m.GetAtoms() if a.GetAtomicNum()==1 and a.GetFormalCharge()>0]
   cache[smi]=';'.join(f'{a.GetSymbol()}({a.GetFormalCharge():+d})' for a in mm)
  except Exception: cache[smi]='unparsed'
 metals.append(cache[smi])
 f=condition_tables.get(r.source)
 T=I=q=e='unreported in released row'
 if f is not None:
  v=f.iloc[r.csv_row-2]
  if r.model.endswith('model_M1'): T,I=v.iloc[:2].tolist()
  elif r.model.endswith('model_M4'): T,I=v.iloc[12:14].tolist();q=v.iloc[8];e=f'ExactMW={v.iloc[14]},density={v.iloc[15]}'
  elif r.model.endswith('model_M3'): T,I=v.iloc[4:6].tolist();q=v.iloc[0];e=f'ExactMW={v.iloc[6]},density={v.iloc[7]}'
 elif r.model.endswith('model_M2'): T=25
 temp.append(T);ionic.append(I);charge.append(q);electrolyte.append(e)
t['metal_context']=metals;t['T_C']=temp;t['I_M']=ionic;t['released_ligand_charge']=charge;t['electrolyte_context']=electrolyte
t.to_csv(W/'released_training_conditions.tsv',sep='\t',index=False)
rows=read(W/'primary.tsv');out=[]
normalization_checks=[]
for r in rows:
 audit_norm=identity(r['smiles'],False)
 normalization_checks.append(dict(lig_id=r['lig_id'],analysis_key=r['normalized_identity'],
                                  chemprop_environment_key=audit_norm,equal=audit_norm==r['normalized_identity']))
 for model,g in t.groupby('model'):
  z=g[g.split.isin(['train','validation']) & (g.normalized_identity==r['normalized_identity'])]
  other=z[~z.metal_context.isin(['La(+3)','Dy(+3)','H(+1)'])]
  out.append(dict(lig_id=r['lig_id'],model=model,n_exact_normalized_records=len(z),
              other_metal_records=len(other),protonation_records=int(z.protonation_record.sum()),
              metals=';'.join(sorted(set(z.metal_context))),
              conditions=';'.join(sorted(set(f'T_C={a},I_M={b}' for a,b in zip(z.T_C,z.I_M)))),
              source_rows=';'.join(z.source+':'+z.csv_row.astype(str))))
write(W/'training_overlap_conditions.tsv',out)
write(W/'audit_normalization_checks.tsv',normalization_checks)
assert all(r['equal'] for r in normalization_checks), 'Normalization-version mismatch requires explicit audit repair'
# Check exact amidated peptide constructs as a separate illustration. No claim of
# independence from related peptides based on absence of these exact strings.
pep=json.load(open(W/'peptide_illustration.json'))
for r in pep['constructs']:
 exact=Chem.MolToSmiles(Chem.MolFromSmiles(r['smiles']))
 loose=Chem.MolFromSmiles(exact);Chem.RemoveStereochemistry(loose);loose=Chem.MolToSmiles(loose)
 heavy_count=Chem.MolFromSmiles(exact).GetNumHeavyAtoms()
 count=0;near=0
 for s in set(t.ligand_smiles):
  if not s: continue
  try:
   m=Chem.MolFromSmiles(s)
   if m.GetNumHeavyAtoms()!=heavy_count: continue
   count+=int(Chem.MolToSmiles(m)==exact)
   Chem.RemoveStereochemistry(m);near+=int(Chem.MolToSmiles(m)==loose)
  except Exception: continue
 r['released_exact_ligand_representations']=count
 r['released_stereo_agnostic_exact_graphs']=near
 r['overlap_limit']='exact strings audited; related peptide-family absence not established'
(W/'peptide_illustration_with_overlap.json').write_text(json.dumps(pep,indent=2)+'\n')
print('condition context attached',len(t),'rows; overlap records',len(out))
