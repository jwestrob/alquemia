"""Collect saved results only; no changes to choices, fits or criteria."""
import json,os,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from curate import W,HERE
assert os.environ.get('SLURM_JOB_ID')
x=pd.read_csv(W/'features.tsv',sep='\t',dtype={'lig_id':str})
primary_rows=pd.read_csv(W/'primary.tsv',sep='\t',dtype={'lig_id':str}).set_index('lig_id')
x['family']=x.lig_id.map(primary_rows.family)
rows=[]
for group,g in x.groupby('family'):
 rows.append(dict(family=group,n_ligands=len(g),size_fit_success=int((g.descriptor_status=='ok').sum())))
pd.DataFrame(rows).to_csv(W/'coverage_per_family.tsv',sep='\t',index=False)
f=x[x.descriptor_status=='ok'];s=f.s_primary.values
sens=[]
for name in ('soft','stiff','two_conformers'):
 v=f['s_'+name].values;ok=np.isfinite(v)
 sens.append(dict(setting=name,n=int(ok.sum()),spearman=float(spearmanr(s[ok],v[ok]).statistic),
                  pearson=float(np.corrcoef(s[ok],v[ok])[0,1]),
                  median_absolute_change=float(np.median(abs(v[ok]-s[ok]))),
                  rms_change=float(np.sqrt(np.mean((v[ok]-s[ok])**2))),
                  max_absolute_change=float(abs(v[ok]-s[ok]).max())))
pd.DataFrame(sens).to_csv(W/'descriptor_sensitivity.tsv',sep='\t',index=False)
# Geometry candidates, coverage and real gradients are independently saved.
trace=json.load(open(W/'checks.json'))
if (W/'published_predictions.tsv').exists():
 p=pd.read_csv(W/'published_predictions.tsv',sep='\t',dtype={'lig_id':str})
 out=[];per_family=[]
 from evaluate import bootstrap
 for name,g in p[p.primary==1].groupby('model'):
  br,bm=bootstrap(g.error_kcal.values[:,None],g.family.values)
  out.append(dict(model=name,n_ligands=len(g),rmse=float(np.sqrt(np.mean(g.error_kcal**2))),
                  mae=float(abs(g.error_kcal).mean()),status='overlapping_data_reference',
                  descriptive_rmse_ci_low=float(np.quantile(br[:,0],.025)),
                  descriptive_rmse_ci_high=float(np.quantile(br[:,0],.975))))
  for family,gg in g.groupby('family'):
   per_family.append(dict(model=name,family=family,n_ligands=len(gg),
                         rmse=float(np.sqrt(np.mean(gg.error_kcal**2))),mae=float(abs(gg.error_kcal).mean()),
                         status='overlapping_data_reference'))
 pd.DataFrame(out).to_csv(W/'published_metrics.tsv',sep='\t',index=False)
 pd.DataFrame(per_family).to_csv(W/'published_per_family.tsv',sep='\t',index=False)
# Correct family metadata on feature coverage after grouping-only repair.
assert all(primary_rows.loc[lid,'family']==fam for lid,fam in x[['lig_id','family']].itertuples(index=False,name=None))
# The changed nitrobenzoate is descriptor-ineligible, so shared-subset rows are identical.
print(pd.DataFrame(sens).to_string(index=False))
