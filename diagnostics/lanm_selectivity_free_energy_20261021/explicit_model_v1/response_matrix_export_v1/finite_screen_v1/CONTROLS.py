#!/usr/bin/env python3
"""Small synthetic algebra and prepared-metadata checks; no actual screening."""
import importlib.util
import json
import math
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/finite_screen_v1'
s=importlib.util.spec_from_file_location('screen_controls',HERE/'SCREEN_v2.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
tick,cpu=time.monotonic(),time.process_time();checks=[]
def check(name,value):
    assert bool(value),name
    checks.append(name)
prep=json.loads((W/'preparation_v2/MANIFEST.json').read_text());units=prep['units']
alpha=np.repeat([1.5,.8],3);A=np.diag(1/alpha);ed=np.array([.1,.2,-.1,.3,.1,.2]);ep=ed*np.array([.7,.8,.9,1.1,1.2,1.3])
g=m.geometry([500,36],[[0,0,0],[0,0,2.5]])
base=m.kernel(A,ed,ep,alpha,g,0.,0.,units);drive=m.kernel(A,ed,ep,alpha,g,.1,0.,units)
check('S2U normalized at2.5',abs(g['s']-1)<1e-15 and abs(g['g']-1)<1e-15)
check('dual closed-form energy',abs(base['polarization_energy']+.5*units['electric']*float((alpha*ed)@ep))<1e-13)
check('dual conjugate half-sum',np.max(abs(np.asarray(base['induced_conjugate_Debye'])-units['debye_per_eA']*((alpha*(ed+ep)).reshape(2,3).sum(axis=0)/2)))<1e-14)
check('drive has zero curvature derivative',base['chi_column_a0_cubed']==drive['chi_column_a0_cubed'])
check('drive changes equilibrium moment',np.linalg.norm(np.asarray(base['induced_conjugate_Debye'])-drive['induced_conjugate_Debye'])>0)
for l in (-.2,.2):
    v=m.kernel(A,ed,ep,alpha,g,0.,l,units)
    expected=units['hartree_kcal_mol']/(units['electric']*units['bohr_A']**2)*(1.5*math.exp(-l)+.8)
    check('signed stiffness closed form '+str(l),abs(v['chi_column_a0_cubed'][2]-expected)<1e-13)
c=units['hartree_kcal_mol']/(units['electric']*units['bohr_A']);z=np.tile([0.,0.,1.],2);h=1e-5
p=m.kernel(A,ed+c*h*z,ep+c*h*z,alpha,g,0.,0.,units);n=m.kernel(A,ed-c*h*z,ep-c*h*z,alpha,g,0.,0.,units)
odd=-(p['polarization_energy']-n['polarization_energy'])/(2*h*units['hartree_kcal_mol'])
check('energy-conjugate dual field moment',abs(odd-base['induced_conjugate_Debye'][2]/(units['debye_per_eA']*units['bohr_A']))<1e-10)
curv=-(p['polarization_energy']+n['polarization_energy']-2*base['polarization_energy'])/(h*h*units['hartree_kcal_mol'])
check('energy curvature consistent',abs(curv-base['chi_column_a0_cubed'][2])<1e-6)
ion=m.geometry([500],[[0.,0.,0.]])
check('bare ion corrections vanish',ion['g']==0 and ion['drive_per_a']==[0.,0.,0.])
bad=A.copy();bad[0,0]=-1
check('unstable point retained without solve',m.kernel(bad,ed,ep,alpha,g,0.,0.,units)['status']=='unstable')
check('one common geometry-only L',prep['L']==math.log(2)/max(x['geometry']['g'] for x in prep['rows'].values()))
check('full54 plus static46',len(prep['rows'])==54 and len(prep['static'])==46)
check('27 fixtures each metal',all(sum(r['metal']==a for r in prep['rows'].values())==27 for a in ('La','Dy')))
check('nine independent points and81 compositions',len(prep['grid_per_metal'])==9 and len({(x['a'],x['l']) for x in prep['grid_per_metal']})**2==81)
check('global local exponent bound',all(abs(p['l']*r['geometry']['g'])<=math.log(2)+2e-15 for p in prep['grid_per_metal'] for r in prep['rows'].values()))
out=W/'source_controls_v1';out.mkdir(exist_ok=False)
result={'status':'passed','implementation':m.pin(__file__),'screen':m.pin(HERE/'SCREEN_v2.py'),'preparation':m.pin(W/'preparation_v2/MANIFEST.json'),'checks':checks,'count':len(checks),'actual_candidate_screen_executed':False,'new_native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
m.save(out/'RESULT.json',result);print(json.dumps({'result':m.pin(out/'RESULT.json'),'count':len(checks),'cost':result['cost']}))
