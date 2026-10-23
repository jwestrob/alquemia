#!/usr/bin/env python3
"""Manufactured-only reviewer controls. No producer, CORE, solve, or data load."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import numpy as np
from scipy import sparse

spec=importlib.util.spec_from_file_location('reviewer',Path(__file__).with_name('REVIEW.py'))
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
assert all(os.environ.get(k)=='1' for k in R.THREADS)
tick,cpu=time.monotonic(),time.process_time();checks=[]

def check(name,condition):
    assert condition,name
    checks.append(name)

def rejects(name,fn):
    try:fn()
    except (ValueError,AssertionError):checks.append(name);return
    raise AssertionError('Not rejected: '+name)

def block(n):
    return {'r':[0.]*n,'J':[[0.]*12 for _ in range(n)],'denominator':1.}

p={'blocks':{name:block(n) for name,n in zip(R.OLD+R.COMPONENTS,(413,102,12,4,2,2))},
   'energy':dict(block(10),names=['energy'+str(i) for i in range(10)]),
   'forces':{str(i):dict(block(n),weights=[1./n]*n) for i,n in enumerate([9,9,9,51,57]*3)},
   'center':[0.]*12,'width':[2.]*12,'lower':[-1.]*12,'upper':[1.]*12,
   'bounds_lower':[-.025]*12+[0.],'bounds_upper':[.025]*12+[1.],'q_current':1.}

# Independent manufactured A*x+s=b construction from named affine quantities.
rows=[];rhs=[];layout=[]
def cone(name,r,J,w,t=False):
    start=len(rhs);top=np.zeros(13);top[-1]=-1. if t else 0.
    rows.append(top);rhs.append(0. if t else 1.)
    for value,j,weight in zip(r,J,w):
        rows.append(np.r_[-weight*np.asarray(j),0.]);rhs.append(weight*value)
    layout.append({'name':name,'kind':'SOC','start':start,'stop':len(rhs)})
def side(name,r,j):
    start=len(rhs);rows.append(j);rhs.append(r)
    layout.append({'name':name,'kind':'nonnegative','start':start,'stop':len(rhs)})
for name,b in p['blocks'].items():
    # Nonzero entries make sign/normalization tests meaningful.
    b['r'][0]=.125;b['J'][0][0]=.25
    cone('block/'+name,b['r'],b['J'],np.ones(len(b['r'])),name in R.COMPONENTS)
for name,b in p['forces'].items():
    b['r'][0]=.125;b['J'][0][1]=.5
    cone('force/RMS/'+name,b['r'],b['J'],np.sqrt(b['weights'])/2)
for name,r,j in zip(p['energy']['names'],p['energy']['r'],p['energy']['J']):
    side('energy/upper/'+name,1-r,np.r_[j,0.]);side('energy/lower/'+name,1+r,np.r_[-np.asarray(j),0.])
for name,b in p['forces'].items():
    for i,(r,j) in enumerate(zip(b['r'],b['J'])):
        side('force/upper/'+name+'/'+str(i),1-r/5,np.r_[np.asarray(j)/5,0.])
        side('force/lower/'+name+'/'+str(i),1+r/5,np.r_[-np.asarray(j)/5,0.])
for i in range(13):
    e=np.eye(13)[i]
    side('bound/lower/'+str(i),-p['bounds_lower'][i],-e)
    side('bound/upper/'+str(i),p['bounds_upper'][i],e)
A=sparse.csc_matrix(np.asarray(rows));b=np.asarray(rhs);c=np.r_[np.zeros(12),1.]
R.check_matrix(p,A,b,c,layout);check('Complete manufactured1817x13 named matrix accepted',True)
bad=A.copy();bad[1,0]*=-1
rejects('Reversed SOC affine sign rejected',lambda:R.check_matrix(p,bad,b,c,layout))
badb=b.copy();badb[535]=1.
rejects('Component epigraph top changed rejected',lambda:R.check_matrix(p,A,badb,c,layout))
badp=copy.deepcopy(p);badp['blocks']['old_EF']['J'][0].pop()
rejects('Truncated Jacobian rejected before zip',lambda:R.validate_problem(badp))
badp=copy.deepcopy(p);badp['forces']['0']['weights'][0]=-1.
rejects('Negative force weight rejected',lambda:R.validate_problem(badp))
badp=copy.deepcopy(p);badp['blocks']['chi']['denominator']=0.
rejects('Nonpositive normalization rejected',lambda:R.validate_problem(badp))
badp=copy.deepcopy(p);badp['center'][0]=float('nan')
rejects('Nonfinite domain rejected',lambda:R.validate_problem(badp))
badlayout=copy.deepcopy(layout);badlayout[-1]['start']=0
rejects('Overlapping row ranges rejected',lambda:R.check_matrix(p,A,b,c,badlayout))
check('SOC negative cone margin retained',R.margins([{'name':'x','kind':'SOC','start':0,'stop':2}],np.array([1.,2.]))['x']==-1.)

for name in ['Solved','AlmostSolved','MaxIterations','InsufficientProgress']:
    s={'solver_status':name,'solver_success':name=='Solved','solution':{'status':name},'reduced_accuracy_status':name.startswith('Almost')}
    check('Truthful '+name+' status',R.check_status(s) is (name=='Solved'))
    s['solver_success']=not s['solver_success']
    rejects('False success label '+name,lambda:R.check_status(s))

zero=copy.deepcopy(p)
for group in [zero['blocks'],zero['forces']]:
    for value in group.values():value['r']=[0.]*len(value['r']);value['J']=[[0.]*12 for _ in value['r']]
boundary=copy.deepcopy(zero);boundary['blocks']['old_EF']['r'][0]=1.
check('Solved exact boundary can be feasible',R.scalar_report(boundary,[0.]*12,.5,True)['eligible'])
check('Failed solver requires strict margin',not R.scalar_report(boundary,[0.]*12,.5,False)['eligible'])
boundary['blocks']['old_EF']['r'][0]=np.nextafter(1.,2.).item()
check('Small negative margin never relaxed',not R.scalar_report(boundary,[0.]*12,.5,True)['eligible'])

# Analytically zero residuals, d=0, t=.5, f=.5 => q=.625.
x=np.r_[np.zeros(12),.5];f=.5;q=.625
for success in [True,False]:
    replay=R.scalar_report(zero,[0.]*12,q,success)
    cert={'method':'independent_scalar_fsum','solver_success':success,'failed_solver_status_remains_failed':not success,
          'q_fraction':q,'q_current':1.,'strict_q_improvement':True,'in_box_and_trust':True,
          'minimum_box_margin':.025,'minimum_residual_constraint_margin':q,'minimum_epigraph_margin':q,
          'all_residual_margins_strictly_positive':True,'all25_predicted':True,'fixed_old_ceiling_pass':True,
          'residual_constraints':replay['margins'],'blocks':replay['blocks'],'original25':replay['original25']}
    slot={'fraction':f,'parameters':[0.]*12,'key':R.digest([0.]*12),'eligible':True,'certificate':cert,
          'q_selection_epsilon':1e-10,'predicted_q_improvement_exceeds_selection_epsilon':True,'native_evaluated':False}
    check('Exact fraction arithmetic status '+str(success),R.check_fraction(zero,x,slot,success)['q']==q)
    for key in ['solver_success','failed_solver_status_remains_failed','strict_q_improvement']:
        altered=copy.deepcopy(slot);altered['certificate'][key]=not altered['certificate'][key]
        rejects('Tampered certificate '+key+' '+str(success),lambda:R.check_fraction(zero,x,altered,success))
    altered=copy.deepcopy(slot);altered['key']='wrong'
    rejects('Tampered vector key '+str(success),lambda:R.check_fraction(zero,x,altered,success))
    altered=copy.deepcopy(slot);altered['q_selection_epsilon']=1e-6
    rejects('Retuned selection epsilon '+str(success),lambda:R.check_fraction(zero,x,altered,success))

out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
result={'status':'passed_manufactured_reviewer_controls','implementation':R.pin(R.__file__),
        'control_implementation':R.pin(__file__),'checks':checks,'checks_passed':len(checks),
        'new_solver_calls':0,'new_native_calls':0,'scientific_vectors_read':0,
        'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
(out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({'result':R.pin(out/'RESULT.json'),'checks':len(checks)}))
