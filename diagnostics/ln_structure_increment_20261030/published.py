"""Pinned Chemprop inference, real application replay and training-exposure audit.

No checkpoint refits or model selection. No experimental affinity is supplied as
an input. Exact-electrolyte M4 coverage is deliberately limited to known NaClO4.
"""
import json
import argparse
import os
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize as std

from curate import W, ROOT, KC, family, identity, write, read

UP=W/'upstream'
MODELS={'published_M1': 'stability_constant_model_M1/best_model',
        'published_M4S': 'stability_constan_protonation_constant_model_M4/SMILES_only_model',
        'published_M4': 'stability_constan_protonation_constant_model_M4/best_model'}


def ligand_part(smi):
    m=Chem.MolFromSmiles(smi)
    if m is None:
        return None
    frags=Chem.GetMolFrags(m,asMols=True)
    organic=[f for f in frags if any(a.GetAtomicNum()==6 for a in f.GetAtoms())]
    if organic:
        return Chem.MolToSmiles(max(organic,key=lambda f:f.GetNumHeavyAtoms()))
    nonmetal=[f for f in frags if all(a.GetAtomicNum() in (1,7,8,9,15,16,17,35,53) for a in f.GetAtoms())
              and not all(a.GetAtomicNum()==1 for a in f.GetAtoms())]
    return Chem.MolToSmiles(max(nonmetal,key=lambda f:f.GetNumHeavyAtoms())) if nonmetal else None


def audit():
    cache={}
    ligand_cache={}
    evalrows=pd.read_csv(W/'primary.tsv',sep='\t',dtype={'lig_id':str})
    def signature(m):
        return tuple(sorted(Counter(a.GetAtomicNum() for a in m.GetAtoms() if a.GetAtomicNum()!=1).items()))
    eligible={signature(Chem.MolFromSmiles(s)) for s in evalrows.smiles}
    def normalize_ligand(ligand):
        if ligand not in ligand_cache:
            try:
                m=Chem.MolFromSmiles(ligand)
                fam=family(m)
                # Protonation/tautomer normalization cannot change heavy-element
                # composition. Only formula-eligible ligands can be exact matches.
                if signature(m) in eligible:
                    norm,stereo=identity(ligand,False),identity(ligand,True)
                    status='full_tautomer_stereo_normalization'
                else:
                    n=std.Uncharger().uncharge(std.Cleanup(m))
                    stereo=Chem.MolToSmiles(n)
                    Chem.RemoveStereochemistry(n)
                    norm=Chem.MolToSmiles(n)
                    status='charge_normalized_formula_ineligible_for_exact_match'
                ligand_cache[ligand]=(ligand,norm,stereo,fam,status)
            except Exception:
                ligand_cache[ligand]=(ligand,'','','unparsed','parse_failure')
        return ligand_cache[ligand]
    training=[]
    roots=list((UP/'data').glob('*model_M*'))
    for folder in roots:
        for path in sorted(folder.glob('*input.csv')):
            split='train' if 'train' in path.name else 'validation' if 'validation' in path.name else 'test'
            data=pd.read_csv(path)
            for i,smi in enumerate(data.iloc[:,0]):
                if smi not in cache:
                    ligand=ligand_part(str(smi))
                    cache[smi]=normalize_ligand(ligand)
                ligand,norm,stereo,fam,norm_status=cache[smi]
                training.append(dict(model=folder.name,split=split,source=str(path.relative_to(UP)),
                                     csv_row=i+2,complex_smiles=smi,ligand_smiles=ligand,
                                     normalized_identity=norm,stereo_identity=stereo,family=fam,
                                     protonation_record=int('[H+]' in str(smi)),normalization_status=norm_status))
    # The older LOGKPREDICT release uses dative metal bonds. Remove every metal
    # before normalizing, preserving the remaining organic graph and charge.
    lw=W/'upstream_logkpredict/data'
    for path in lw.glob('*/ML_inp_train.csv'):
        for i,smi in enumerate(pd.read_csv(path).iloc[:,0]):
            try:
                mol=Chem.MolFromSmiles(str(smi))
                rw=Chem.RWMol(mol)
                for idx in reversed([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() not in (1,6,7,8,9,15,16,17,35,53)]):
                    rw.RemoveAtom(idx)
                mol=rw.GetMol(); Chem.SanitizeMol(mol)
                ligand=ligand_part(Chem.MolToSmiles(Chem.RemoveHs(mol)))
                ligand,norm,stereo,fam,norm_status=normalize_ligand(ligand)
            except Exception:
                ligand,norm,stereo,fam,norm_status='','','','unparsed','parse_failure'
            training.append(dict(model='LOGKPREDICT_'+path.parent.name,split='train',source=str(path),csv_row=i+2,
                                 complex_smiles=smi,ligand_smiles=ligand,normalized_identity=norm,
                                 stereo_identity=stereo,family=fam,protonation_record=int('[H+]' in str(smi)),normalization_status=norm_status))
    t=pd.DataFrame(training)
    t.to_csv(W/'released_training_inventory.tsv',sep='\t',index=False)
    rows=pd.read_csv(W/'primary.tsv',sep='\t',dtype={'lig_id':str})
    audit_rows=[]
    # Close analogues: radius-2 2048-bit Morgan Tanimoto >=0.6 (audit only).
    unique=t[(t.normalized_identity!='') & t.split.isin(['train','validation'])].normalized_identity.unique()
    fps=[AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s),2,nBits=2048) for s in unique]
    for r in rows.to_dict('records'):
        fp=AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(r['normalized_identity']),2,nBits=2048)
        sim=np.array(DataStructs.BulkTanimotoSimilarity(fp,fps))
        closest=str(unique[int(sim.argmax())])
        for model,sub in t.groupby('model'):
            train=sub[sub.split.isin(['train','validation'])]
            exact=train[train.normalized_identity==r['normalized_identity']]
            stereo=exact[exact.stereo_identity==r['stereo_identity']]
            samefam=train[train.family==r['family']]
            rec=dict(lig_id=r['lig_id'],family=r['family'],model=model,
                     exact_stereo_train_validation=len(stereo),normalized_train_validation=len(exact),
                     protonation_train_validation=int(exact.protonation_record.sum()),
                     stability_train_validation=int((1-exact.protonation_record).sum()),
                     same_family_train_validation=len(samefam),
                     exact_test=int(((sub.split=='test') & (sub.normalized_identity==r['normalized_identity'])).sum()),
                     normalized_representations=';'.join(exact.complex_smiles.drop_duplicates()),
                     source_rows=';'.join(exact.source+':'+exact.csv_row.astype(str)),
                     closest_any_release=closest,closest_any_release_tanimoto=float(sim.max()),
                     close_any_release_count=int((sim>=.6).sum()),
                     status='overlapping_data_reference' if len(samefam) else 'family_absence_requires_review')
            audit_rows.append(rec)
    write(W/'training_overlap.tsv',audit_rows)
    summary={'released_rows':len(t),'unique_complexes':len(cache),'unique_ligands_normalized':len(ligand_cache),
             'unparsed_rows':int((t.family=='unparsed').sum()),
             'normalization_scope':t.normalization_status.value_counts().to_dict(),
             'audit_rows':len(audit_rows),'all_primary_families_exposed':all(r['same_family_train_validation']>0 for r in audit_rows
                  if r['model'].endswith(('model_M1','model_M4'))),
             'clean_comparison':'unresolved: no verified external held-out family; no clean from-scratch refits in this bounded run'}
    (W/'overlap_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('OVERLAP',json.dumps(summary),flush=True)


def predict(name,smiles,features,dest):
    from chemprop.args import PredictArgs
    from chemprop.train import make_predictions
    dest.mkdir(parents=True,exist_ok=True)
    inp=dest/'input.csv'
    pd.DataFrame({'smiles':smiles}).to_csv(inp,index=False)
    feature_args=[]
    if features is not None:
        pd.DataFrame(features).to_csv(dest/'features.csv',index=False)
        feature_args=['--features_path',str(dest/'features.csv')]
    paths=sorted((UP/'models'/MODELS[name]).glob('fold_*/model_0/model.pt'))
    assert len(paths)==5
    out=[]
    for member,path in enumerate(paths):
        args=PredictArgs().parse_args(['--test_path',str(inp),'--smiles_columns','smiles',
               '--checkpoint_path',str(path),'--preds_path',str(dest/f'member_{member}.csv'),
               '--num_workers','0','--no_cuda',*feature_args])
        p=np.asarray(make_predictions(args=args),dtype=float).reshape(-1)
        assert len(p)==len(smiles) and np.isfinite(p).all()
        out.append(p)
    return np.array(out)


def application():
    # Exact checked-in LOGKPREDICT example: molecule, features and expected output.
    lw=W/'upstream_logkpredict'
    expected=float((lw/'logk_output').read_text().strip())
    # The released modern Python API masks to 21 features, whereas this exact
    # checkpoint expects 48. Use the released original executable unchanged.
    dest=W/'application_LOGKPREDICT'
    dest.mkdir(exist_ok=True)
    if not (dest/'logk_output').exists():
        shutil.copyfile(lw/'logk_input',dest/'logk_input')
        env=os.environ.copy();env['LOGKPREDICT_DIR']=str(lw)
        subprocess.run([sys.executable,str(lw/'LOGKPREDICT')],cwd=dest,env=env,check=True)
    actual=float((dest/'logk_output').read_text().strip())
    replay={'name':'released LOGKPREDICT Ni(II)-ethylenediamine',
            'expected_logK':expected,'observed_logK':actual,'abs_difference':abs(actual-expected),
            'pass_tolerance':.02,'passed':abs(actual-expected)<.02,
            'note':'Real released application; raw logK placeholders are masked out by upstream preprocessing.'}
    (W/'application_replay.json').write_text(json.dumps(replay,indent=2)+'\n')
    # Preserve the failed legacy replay; it is not a verification of this
    # checkpoint's published numeric output. The separate 2025 M4-S example below
    # uses its own released inputs/checkpoints and does not require that API.
    # 2025 paper Figure 8(b) first two real application inputs, using M4-S as in paper.
    app=['CC(=O)[O-].O=[U+2]=O','C(=O)(C(=O)[O-])[O-].O=[U+2]=O']
    p=predict('published_M4S',app,None,W/'application_M4S')
    pd.DataFrame({'smiles':app,'experimental_logK':[3.1,7.3],
                  'predicted_logK':p.mean(axis=0),'member_sd':p.std(axis=0,ddof=1)}).to_csv(W/'application_M4S/result.csv',index=False)
    (W/'application_M4S/preprocessing_verification.json').write_text(json.dumps({
        'source':'2025 paper Figure 8(b), released applications_DFT_vs_ML_predictions.csv',
        'inputs':app,'checkpoint_members':5,'inorganic_ligand_charge':[-1,-2],
        'metal_fragment':'O=[U+2]=O','total_graph_charge':[1,0],
        'task':'Exp LgK1','extra_features':None,'finite_predictions':True,
        'limit':'Exact numerical ML values are not tabulated in released application data; no numeric-paper replay claim.'},indent=2)+'\n')
    print('APPLICATION',json.dumps(replay),flush=True)


def evaluate_reference():
    rows=read(W/'evaluation_records.tsv')
    baseline=pd.read_csv(UP/'data/stability_constant_protonation_constant_model_M4/train_features.csv')
    # Metal-only properties copied from the released data, with constancy checked.
    atom={}
    for metal,z in [('La',57),('Dy',66)]:
        sel=baseline[(baseline.iloc[:,0]==z)&(baseline.iloc[:,9]==3)]
        cols=list(range(8))+[11]
        assert len(sel) and all(sel.iloc[:,i].nunique()==1 for i in cols)
        atom[metal]=sel.iloc[0].values.astype(float)
    salt=baseline.loc[np.isclose(baseline.iloc[:,-2],121.9382804),baseline.columns[-2:]].drop_duplicates()
    assert len(salt)==1
    saltvalues=salt.iloc[0].values.astype(float)
    output=[]
    missing=[]
    for name in MODELS:
        selected=rows
        if name=='published_M4':
            selected=[r for r in rows if r['native_annotation'].startswith('aa: NaClO4')]
            missing.extend(dict(record_id=r['record_id'],lig_id=r['lig_id'],model=name,
                                reason='exact background electrolyte not specified; 16-feature model not imputed')
                           for r in rows if r not in selected)
        smiles=[];features=[]
        for r in selected:
            for metal in ('La','Dy'):
                smiles.append(r['smiles']+f'.[{metal}+3]')
                if name=='published_M1':
                    features.append([25.,float(r['I_M'])])
                elif name=='published_M4':
                    v=atom[metal].copy()
                    q=float(r['ligand_charge'])
                    v[8:11]=[q,3,q+3]
                    v[12:]=[25.,float(r['I_M']),*saltvalues]
                    features.append(v.tolist())
        pred=predict(name,smiles,None if name=='published_M4S' else features,W/name)
        for i,r in enumerate(selected):
            contrast=-KC*(pred[:,2*i+1]-pred[:,2*i])
            out={k:r[k] for k in ('record_id','lig_id','family','T_K','I_M','primary','y_kcal','native_annotation')}
            out.update(model=name,fold='released_pretrained_ensemble',overlap_status='overlapping_data_reference',
                       prediction_kcal=float(contrast.mean()),member_sd_kcal=float(contrast.std(ddof=1)),
                       predicted_logK_La=float(pred[:,2*i].mean()),predicted_logK_Dy=float(pred[:,2*i+1].mean()),
                       error_kcal=float(contrast.mean())-float(r['y_kcal']))
            for j,x in enumerate(contrast):
                out[f'member_{j}_ddG_kcal']=float(x)
            output.append(out)
        write(W/'published_predictions.tsv',output)
    write(W/'published_predictions.tsv',output)
    write(W/'published_exclusions.tsv',missing)


def main():
    assert os.environ.get('SLURM_JOB_ID')
    torch.set_num_threads(int(os.environ.get('SLURM_CPUS_PER_TASK','1')))
    start=time.time();cpu=time.process_time()
    parser=argparse.ArgumentParser()
    parser.add_argument('--reuse-audit',action='store_true')
    args=parser.parse_args()
    if args.reuse_audit:
        assert (W/'released_training_inventory.tsv').exists() and (W/'overlap_summary.json').exists()
        a=pd.read_csv(W/'training_overlap.tsv',sep='\t',dtype={'lig_id':str})
        assert set(a.lig_id)==set(pd.read_csv(W/'primary.tsv',sep='\t',dtype={'lig_id':str}).lig_id)
        print('Reusing completed 32006 overlap inventory; no normalization repeated.',flush=True)
    else:
        audit()
    application()
    evaluate_reference()
    (W/'published_cost.json').write_text(json.dumps({'wall_seconds':time.time()-start,
        'process_cpu_seconds':time.process_time()-cpu,'job_id':os.environ['SLURM_JOB_ID']},indent=2)+'\n')


if __name__=='__main__':
    main()
