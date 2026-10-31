"""Frozen nested family-held-out ablations and family-level uncertainty."""
import hashlib
import itertools
import json
import os
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from curate import HERE, RECORDS, ROOT, W
from structure import M1, TOPOLOGY

ALPHAS=[.01,.1,1.,10.,100.]


def ridge_fit(x,y,alpha):
    return make_pipeline(StandardScaler(),Ridge(alpha=alpha)).fit(x,y)


def nested(d,cols):
    y=d.y_kcal.values
    groups=d.family.values
    pred=np.full(len(d),np.nan)
    metadata=[]
    for family in sorted(set(groups)):
        test=groups==family
        train=~test
        assert not set(d.normalized_identity[test])&set(d.normalized_identity[train]), family
        if cols is None:
            pred[test]=y[train].mean()
            metadata.append(dict(fold=family,alpha=None,n_train=int(train.sum()),n_test=int(test.sum())))
            continue
        x=d[cols].values.astype(float)
        scores=[]
        for alpha in ALPHAS:
            inner=np.full(len(d),np.nan)
            for held in sorted(set(groups[train])):
                va=train&(groups==held)
                tr=train&(groups!=held)
                assert not set(d.normalized_identity[va])&set(d.normalized_identity[tr])
                inner[va]=ridge_fit(x[tr],y[tr],alpha).predict(x[va])
            scores.append(float(np.mean((inner[train]-y[train])**2)))
        alpha=ALPHAS[int(np.argmin(scores))]
        mdl=ridge_fit(x[train],y[train],alpha)
        pred[test]=mdl.predict(x[test])
        metadata.append(dict(fold=family,alpha=alpha,n_train=int(train.sum()),n_test=int(test.sum()),
                             inner_mse=dict(zip(map(str,ALPHAS),scores)),
                             scaler_mean=mdl[0].mean_.tolist(),scaler_scale=mdl[0].scale_.tolist(),
                             coefficients=mdl[1].coef_.tolist(),intercept=float(mdl[1].intercept_),features=cols))
    assert np.isfinite(pred).all()
    return pred,metadata


def bootstrap(errors,groups):
    rng=np.random.default_rng(20261030)
    unique=sorted(set(groups))
    group_idx={g:np.flatnonzero(groups==g) for g in unique}
    rmse=[];mae=[]
    for _ in range(2000):
        idx=np.concatenate([group_idx[g] for g in rng.choice(unique,len(unique),replace=True)])
        rmse.append(np.sqrt(np.mean(errors[idx]**2,axis=0)))
        mae.append(np.mean(abs(errors[idx]),axis=0))
    return np.array(rmse),np.array(mae)


def main():
    assert os.environ.get('SLURM_JOB_ID')
    start=time.time();cpu=time.process_time()
    freeze=json.loads((RECORDS/'FREEZE.json').read_text())
    for p,sha in freeze['pins'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha, p
    d=pd.read_csv(W/'primary.tsv',sep='\t',dtype={'lig_id':str})
    f=pd.read_csv(W/'features.tsv',sep='\t',dtype={'lig_id':str})
    d=d.merge(f.drop(columns=['family','ligand_charge']),on='lig_id',validate='one_to_one')
    assert len(d)==freeze['n_primary'] and d.record_id.is_unique and d.lig_id.is_unique
    assert np.allclose(d.y_kcal,-freeze['kcal_per_log_unit']*(d.logK_Dy-d.logK_La),atol=1e-12)
    subsets={'full':d, 'size_fit':d[d.descriptor_status=='ok']}
    sensitivity_cols=['s_primary','s_soft','s_stiff','s_two_conformers']
    subsets['sensitivity_common']=d.dropna(subset=sensitivity_cols)
    predictions=[];modelmeta={};metrics=[];contrasts=[];familyrows=[];challengerows=[]
    for subset,dd in subsets.items():
        dd=dd.reset_index(drop=True)
        models={'M0':None,'M1':M1,'M2':M1+TOPOLOGY}
        if subset!='full':
            models['M3']=M1+TOPOLOGY+['s_primary']
        if subset=='sensitivity_common':
            models.update({f'M3_{s}':M1+TOPOLOGY+[f's_{s}'] for s in ('soft','stiff','two_conformers')})
        errors=[]
        for name,cols in models.items():
            pred,meta=nested(dd,cols)
            modelmeta[subset+'/'+name]=meta
            errors.append(pred-dd.y_kcal.values)
            for i,r in dd.iterrows():
                predictions.append(dict(subset=subset,model=name,fold=r.family,record_id=r.record_id,
                    lig_id=r.lig_id,family=r.family,T_K=r.T_K,I_M=r.I_M,primary=1,
                    tetraacetate_challenge=r.tetraacetate_challenge,
                    y_kcal=r.y_kcal,prediction_kcal=pred[i],error_kcal=pred[i]-r.y_kcal,
                    overlap_status='family_excluded_from_all_fit_and_model_selection'))
            for family,g in dd.groupby('family'):
                err=pred[g.index]-g.y_kcal.values
                familyrows.append(dict(subset=subset,model=name,family=family,n_ligands=len(g),
                                        rmse=float(np.sqrt(np.mean(err**2))),mae=float(abs(err).mean())))
            challenge=dd.tetraacetate_challenge.astype(bool).values
            if challenge.any():
                err=pred[challenge]-dd.y_kcal.values[challenge]
                challengerows.append(dict(subset=subset,model=name,n_ligands=int(challenge.sum()),
                    rmse=float(np.sqrt(np.mean(err**2))),mae=float(abs(err).mean()),
                    status='previously_inspected_development_challenge',
                    exclusions='whole aminocarboxylate family excluded from training and model selection'))
        errors=np.array(errors).T
        br,bm=bootstrap(errors,dd.family.values)
        names=list(models)
        for j,name in enumerate(names):
            metrics.append(dict(subset=subset,model=name,n_ligands=len(dd),n_families=dd.family.nunique(),
                rmse=float(np.sqrt(np.mean(errors[:,j]**2))),mae=float(abs(errors[:,j]).mean()),
                rmse_ci_low=float(np.quantile(br[:,j],.025)),rmse_ci_high=float(np.quantile(br[:,j],.975)),
                mae_ci_low=float(np.quantile(bm[:,j],.025)),mae_ci_high=float(np.quantile(bm[:,j],.975))))
        for before,after in [('M0','M1'),('M1','M2'),('M2','M3'),*([('M2','M3_'+s) for s in ('soft','stiff','two_conformers')] if subset=='sensitivity_common' else [])]:
            if after not in names:
                continue
            a,b=names.index(before),names.index(after)
            delta=br[:,a]-br[:,b]
            point=float(np.sqrt(np.mean(errors[:,a]**2))-np.sqrt(np.mean(errors[:,b]**2)))
            contrasts.append(dict(subset=subset,before=before,after=after,rmse_improvement=point,
                ci_low=float(np.quantile(delta,.025)),ci_high=float(np.quantile(delta,.975)),
                mean_paired_absolute_error_improvement=float((abs(errors[:,a])-abs(errors[:,b])).mean())))
    pd.DataFrame(predictions).to_csv(W/'predictions.tsv',sep='\t',index=False)
    # Row-level paired error differences, with no unpaired denominator changes.
    p=pd.DataFrame(predictions)
    paired=[]
    for (subset,lid),rr in p.groupby(['subset','lig_id']):
        rr=rr.set_index('model')
        for before,after in [('M0','M1'),('M1','M2'),('M2','M3')]:
            if after in rr.index:
                a,b=rr.loc[before],rr.loc[after]
                paired.append(dict(subset=subset,lig_id=lid,family=a.family,before=before,after=after,
                    absolute_error_improvement=abs(a.error_kcal)-abs(b.error_kcal),
                    squared_error_improvement=a.error_kcal**2-b.error_kcal**2))
    pd.DataFrame(paired).to_csv(W/'paired_errors.tsv',sep='\t',index=False)
    pd.DataFrame(metrics).to_csv(W/'metrics.tsv',sep='\t',index=False)
    pd.DataFrame(familyrows).to_csv(W/'per_family.tsv',sep='\t',index=False)
    pd.DataFrame(contrasts).to_csv(W/'paired_improvements.tsv',sep='\t',index=False)
    pd.DataFrame(challengerows).to_csv(W/'challenge.tsv',sep='\t',index=False)
    (W/'fit_metadata.json').write_text(json.dumps(modelmeta,indent=2)+'\n')
    # Alternate conditions use their ligand's primary outer fold prediction;
    # these are linked diagnostics, never added to the primary sample count.
    allrows=pd.read_csv(W/'evaluation_records.tsv',sep='\t',dtype={'lig_id':str})
    alt=allrows[allrows.primary==0][['record_id','lig_id','T_K','I_M','y_kcal']]
    alt=alt.merge(p.drop(columns=['record_id','T_K','I_M','y_kcal','error_kcal','primary']),on='lig_id')
    alt['error_kcal']=alt.prediction_kcal-alt.y_kcal
    alt['primary']=0
    alt.to_csv(W/'alternate_condition_predictions.tsv',sep='\t',index=False)
    costs={'wall_seconds':time.time()-start,'process_cpu_seconds':time.process_time()-cpu,
           'job_id':os.environ['SLURM_JOB_ID']}
    (W/'evaluation_cost.json').write_text(json.dumps(costs,indent=2)+'\n')
    print(pd.DataFrame(metrics).to_string(index=False),flush=True)
    print(pd.DataFrame(contrasts).to_string(index=False),flush=True)


if __name__=='__main__':
    main()
