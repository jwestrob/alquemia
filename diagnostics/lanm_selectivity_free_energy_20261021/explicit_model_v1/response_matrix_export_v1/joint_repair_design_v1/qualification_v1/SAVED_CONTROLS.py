#!/usr/bin/env python3
"""Narrow synthetic readback/field arithmetic controls, zero native calls."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q',HERE/'RUN.py');q=importlib.util.module_from_spec(s);s.loader.exec_module(q)
out=Path(sys.argv[1]).resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time();checks=[]
def check(name,condition):
    if not condition:raise AssertionError(name)
    checks.append(name)
def rejects(name,fn):
    try:fn()
    except (ValueError,KeyError):checks.append(name);return
    raise AssertionError('Accepted '+name)
sources=q.read(q.pin(HERE/'SOURCE_INPUTS.json'));h=q.helpers(sources);jobs=q.inventory()
check('57 unique prospective programs',len(jobs)==len({x['id']for x in jobs})==57)
atoms={g:(1 if g in q.IONS else 6 if g==q.PURE else 4 if g in q.DIMERS else 7 if g in q.NEAR else 28)for g in q.CELLS}
check('636 exact internal displacements',sum(6*atoms[j['cell']]for j in jobs if j['fd'])==636)
check('no export during any testgrad',all(not j.get('matrix')for j in jobs if j['program']=='testgrad'))
check('six matrices and four fields',sum(bool(j.get('matrix'))for j in jobs)==6 and sum('field_au'in j for j in jobs)==4)
row={'atoms':4,'types':[500,36,37,37]};job={'profile':'combined'}
text=''.join('LOCAL_RESPONSE_PARAMETERS '+str(t)+' 36 .05 -.04 .08 -.06\n'for t in (500,501))+'LOCAL_RESPONSE_PAIR 1 2 .2 .4 2.5\nLOCAL_RESPONSE_ENVIRONMENT 1 1 0 .08 .1\nLOCAL_RESPONSE_SPECTRUM .5 2 0 0 0\n'+''.join(f'LOCAL_RESPONSE_INDUCED {i} 0 0 0 0 0 0\n'for i in range(1,5))+'LOCAL_RESPONSE_DUAL 0 0\n'
f=out/'synthetic_native.txt';f.write_text(text);check('complete native readback accepted',q.response_readback(f,job,row)['passed'])
for name,bad in [('ignored keyword',''),('missing parameter',text.replace('LOCAL_RESPONSE_PARAMETERS 501 36 .05 -.04 .08 -.06\n','')),('wrong coefficient',text.replace('.05 -.04 .08 -.06','.06 -.04 .08 -.06')),('missing site',text.replace('LOCAL_RESPONSE_INDUCED 4 0 0 0 0 0 0\n','')),('nonfinite',text.replace('SPECTRUM .5 2','SPECTRUM nan 2')),('nonpositive spectrum',text.replace('SPECTRUM .5 2','SPECTRUM -1 2')),('wrong record shape',text.replace('SPECTRUM .5 2 0 0 0','SPECTRUM .5 2 0 0'))]:
    f.write_text(bad);rejects(name,lambda:q.response_readback(f,job,row))
f.write_text(text);rejects('zero branch unexpectedly activated',lambda:q.response_readback(f,{'profile':'zero'},row));f.write_text('');check('literal zero bypass',q.response_readback(f,{'profile':'zero'},row)['passed'])
# Analytically solvable artificial 4-site harmonic operator; not a scientific point.
n=4;dim=12;A=np.diag(np.linspace(2.,3.,dim));ed=np.linspace(-.01,.03,dim);u=np.linalg.solve(A,ed);units={'bohr_A':.529177210544,'hartree_kcal_mol':627.509474063,'electric':332.0637133,'debye_per_eA':4.803206802};B=units['bohr_A'];Ha=units['hartree_kcal_mol'];K=units['electric'];D=units['debye_per_eA'];z=np.tile([0.,0.,1.],n);coupling=Ha/(K*B);permanent=np.array([1.,2.,3.]);e0=-5.;g='synthetic'
export={'A':A.tolist(),'vectors':np.column_stack([ed,ed,u,u]).tolist(),'sites':[{'atom':i+1}for i in range(n)],'spectrum':{'A':np.linalg.eigvalsh(A).tolist()},'numerical_bounds':{'post_peek_amplification':1.,'native_vs_dense_u':1e-14,'native_vs_dense_v':1e-14}}
records={'matrix__'+g:{'matrix_qualification':{'exports':[export,export]},'moments':{'energy_kcal_mol':e0}}}
for sign in (-1,1):
    field=sign*.001;p=ed+coupling*field*z;v=np.linalg.solve(A,p);energy=e0-Ha*field*permanent[2]/(D*B)-.5*K*(v@p-u@ed);outfile=out/f'synthetic_field_{sign}.txt';outfile.write_text('Induced Dipoles : Iterations 2 RMS Residual 0.0000000000\n'*2)
    records[f'field__{g}__{sign:+d}']={'response_readback':{'last_induced_eA':np.column_stack([v.reshape(n,3),v.reshape(n,3)]).tolist()},'moments':{'energy_kcal_mol':float(energy),'induced_Debye':(D*v.reshape(n,3)).tolist(),'total_lab_Debye':(permanent+D*v.reshape(n,3).sum(axis=0)).tolist(),'analyze_output':q.pin(outfile)}}
row=dict(atoms=n,permanent_lab_Debye=permanent.tolist());r=q.field_checks(g,records,row,units,h);check('exact harmonic dual energy/moment/curvature',r['passed'])
bad=copy.deepcopy(records);bad['field__synthetic__+1']['moments']['energy_kcal_mol']+=1e-4;rejects('wrong field energy factor/sign',lambda:q.field_checks(g,bad,row,units,h))
bad=copy.deepcopy(records);bad['field__synthetic__+1']['response_readback']['last_induced_eA'][0][0]+=.01;rejects('wrong native field RHS solution',lambda:q.field_checks(g,bad,row,units,h))
bad=copy.deepcopy(records);del bad['field__synthetic__-1'];rejects('missing signed field',lambda:q.field_checks(g,bad,row,units,h))
q.save(out/'RESULT.json',dict(status='passed',implementation=q.pin(__file__),qualification=q.pin(HERE/'RUN.py'),checks=checks,check_count=len(checks),synthetic=True,new_native_calls=0,new_scientific_points=0,cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}));print(json.dumps({'checks':len(checks),'result':q.pin(out/'RESULT.json')}))
