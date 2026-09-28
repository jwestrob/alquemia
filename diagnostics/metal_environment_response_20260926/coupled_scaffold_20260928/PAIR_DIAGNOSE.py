"""Pairwise diagnosis at the three saved real configurations; no new geometry."""
from pathlib import Path
import sys,json,heapq
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,write_new
b=ROOT/'workspaces/metal_environment_response_20260926/coupled_scaffold_classical_v1';p=read_json(b/'particles.json');ex=read_json(b/'exceptions.json');xyz={k:np.array(read_json(b/f'{k}_physical.json')['xyz_A']) for k in ['origin','minus','plus']};ids={v['id']:i for i,v in enumerate(p)};x=xyz['origin'];ni,ci=ids['A/83//N'],ids['A/83//CA'];axis=x[ci]-x[ni];axis/=np.linalg.norm(axis);tangent=np.zeros_like(x);moving=np.linalg.norm(xyz['plus']-x,axis=1)>1e-12;tangent[moving]=np.cross(axis,x[moving]-x[ci]);mm=[i for i,v in enumerate(p) if v['region']=='MM'];exc={tuple(sorted(v['atoms'])):v for v in ex};h=np.deg2rad(.05);worst=[];derivsum=0.;fdsum=0.
for ii,i in enumerate(mm[:-1]):
 jj=np.array(mm[ii+1:]);sig=np.array([(p[i]['sigma_nm']+p[j]['sigma_nm'])*5 for j in jj]);eps=np.array([np.sqrt(p[i]['epsilon_kJ_mol']*p[j]['epsilon_kJ_mol'])/4.184 for j in jj])
 for k,j in enumerate(jj):
  row=exc.get(tuple(sorted((i,int(j)))))
  if row:sig[k]=row['sigma_nm']*10;eps[k]=row['epsilon_kJ_mol']/4.184
 ds={label:np.linalg.norm(xx[i]-xx[jj],axis=1) for label,xx in xyz.items()};energies={label:4*eps*((sig/r)**12-(sig/r)**6) for label,r in ds.items()};r=ds['origin'];rvec=x[i]-x[jj];dr=np.sum(rvec*(tangent[i]-tangent[jj]),axis=1)/r;s6=(sig/r)**6;analytic=24*eps/r*(s6-2*s6*s6)*dr;fd=(energies['plus']-energies['minus'])/(2*h);res=fd-analytic;derivsum+=analytic.sum();fdsum+=fd.sum()
 for k in np.argsort(abs(res))[-5:]:
  row={'ids':[p[i]['id'],p[int(jj[k])]['id']],'elements':[p[i]['element'],p[int(jj[k])]['element']],'distances_A':{label:float(v[k]) for label,v in ds.items()},'energies_kcal_mol':{label:float(v[k]) for label,v in energies.items()},'analytic_kcal_mol_rad':float(analytic[k]),'central_difference':float(fd[k]),'residual':float(res[k]),'native_exception':tuple(sorted((i,int(jj[k])))) in exc};worst.append(row)
worst=sorted(worst,key=lambda r:abs(r['residual']),reverse=True)[:10];out={'implementation':record(__file__),'ledger':record(b/'LEDGER.json'),'analytic_sum':float(derivsum),'central_difference_sum':float(fdsum),'residual_sum':float(fdsum-derivsum),'largest_residual_pairs':worst,'new_geometries':0,'new_OpenMM_calls':0};write_new(b/'PAIR_DIAGNOSIS.json',out);print(json.dumps(out,indent=2))
