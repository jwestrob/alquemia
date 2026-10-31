"""Separate fixed external illustration, exact amidated 2025 LBT3 constructs.

No peptide simulation, parameterization or fitting. Temperature of the competition
assay is unreported in SI S6–S7, so compare log preference rather than inventing
an experimental free energy at 298 K.
"""
import json,os
import numpy as np
from rdkit import Chem
from published import predict
from curate import W,identity
assert os.environ.get('SLURM_JOB_ID')
seqs={'WT':'FIDTNNDGWYEGDELLA','N5D':'FIDTDNDGWYEGDELLA'}
rows=[];smiles=[]
for label,seq in seqs.items():
 m=Chem.MolFromSequence(seq)
 oxy=[a for a in m.GetAtoms() if a.GetPDBResidueInfo().GetName().strip()=='OXT']
 assert len(oxy)==1 and oxy[0].GetPDBResidueInfo().GetResidueNumber()==17
 rw=Chem.RWMol(m);a=rw.GetAtomWithIdx(oxy[0].GetIdx());a.SetAtomicNum(7);a.SetNumExplicitHs(2);a.SetNoImplicit(True)
 m=rw.GetMol();Chem.SanitizeMol(m)
 acid=Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
 assert len(m.GetSubstructMatches(acid))==(5 if label=='WT' else 6)
 assert Chem.GetFormalCharge(m)==0
 smi=Chem.MolToSmiles(m,isomericSmiles=True)
 rows.append({'construct':label,'sequence':seq,'smiles':smi,'N_terminus':'free amine','C_terminus':'amide',
              'assay':'Tb displacement / luminescence competition','pH':6,'MES_M':.05,'NaCl_M':.1,
              'assay_T_K':None,'temperature_status':'not stated for competition assay; NMR 25 C is a separate method',
              'source':'10.1002/chem.202501318 SI S3, S6–S7, Table S3'})
 for metal in ('La','Lu'): smiles.append(smi+f'.[{metal}+3]')
p=predict('published_M4S',smiles,None,W/'peptide_M4S')
contrast=p[:,1::2]-p[:,0::2] # logK(Lu)-logK(La), within each ensemble member
for i,r in enumerate(rows):
 r['predicted_logK_Lu_minus_La']=float(contrast[:,i].mean())
 r['member_sd_logK']=float(contrast[:,i].std(ddof=1))
 r['experimental_Kd_La_nM']=14202.6 if i==0 else 864.9
 r['experimental_Kd_Lu_nM']=139. if i==0 else 353.3
 r['experimental_logK_Lu_minus_La']=float(np.log10(r['experimental_Kd_La_nM']/r['experimental_Kd_Lu_nM']))
change=contrast[:,1]-contrast[:,0]
report={'constructs':rows,'predicted_mutation_change_log_preference':float(change.mean()),
        'predicted_member_sd_log_preference':float(change.std(ddof=1)),
        'experimental_mutation_change_log_preference':rows[1]['experimental_logK_Lu_minus_La']-rows[0]['experimental_logK_Lu_minus_La'],
        'status':'separate external illustration, not small-ligand validation or protein transfer proof',
        'source_note':'SI Table S3 La WT 14202.6 nM used, rather than rounded main-text 14211; assay T not imputed',
        'interpretation':'total mutation effect includes added carboxylate/charge and coordination; no isolated H-bond attribution'}
(W/'peptide_illustration.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='constructs'},indent=2))
