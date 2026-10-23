"""Narrow manufactured validation controls; never builds or runs native code."""
import copy
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time
import numpy as np

source=Path(__file__).resolve().parent.parent/'EXECUTE_v2.py'
spec=importlib.util.spec_from_file_location('initialization_runner',source)
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'))
tick,cpu=time.monotonic(),time.process_time();checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def rejects(name,fn):
    try:fn()
    except (AssertionError,ValueError):checks.append(name);return
    raise AssertionError('No rejection: '+name)

norms=[1.,2.,3.,12.,.01];qstats=[0.,.25,.25,0.,.001]
keys=['c1','cinf','cF','U_R_nominal','evenness_perturbation_nominal','q_min','q_max','q_abs_max','q_zero','q_conjugate_difference_max']
replay=dict(zip(keys,norms+qstats));state=dict(sites=[None]*11958,p=5,native_norms=norms,q_stats=qstats)
report=R.norm_comparison(state,replay)
check('Exact manufactured summary accepted',report['passed'])
check('Fixed full-size positive-accumulation count',report['accumulation_operations']==4484378)
check('Small-case tolerance preserved as floor',R.norm_comparison(dict(state,sites=[None],p=1),replay)['norm_relative_allowance']==1e-12)
check('Full-size allowance equals declared positive-accumulation formula',report['norm_relative_allowance']==32*(4484378*R.EPS)/(1-4484378*R.EPS))
for i,key in enumerate(keys):
    bad=copy.deepcopy(state);arr=bad['native_norms'] if i<5 else bad['q_stats'];j=i if i<5 else i-5
    arr[j]+=100*report['checks'][i]['bound']
    check('Corrupt native summary rejected '+key,not R.norm_comparison(bad,replay)['passed'])
bad=copy.deepcopy(state);bad['native_norms'][0]=float('nan')
rejects('Nonfinite native summary rejected',lambda:R.norm_comparison(bad,replay))

for w,expected_values,expected_derivatives in [
    (0.,[1/24,11/24,11/24,1/24,0.],[-1/6,-.5,.5,1/6,0.]),
    (1.,[0.,1/24,11/24,11/24,1/24],[0.,-1/6,-.5,.5,1/6])]:
    actual=R.spline5(w)
    check('Order5 analytic values '+str(w),max(abs(a[0]-b) for a,b in zip(actual,expected_values))<1e-15)
    check('Order5 analytic first derivative '+str(w),max(abs(a[1]-b) for a,b in zip(actual,expected_derivatives))<1e-15)
for w in (.125,.5,.875):
    actual=R.spline5(w)
    check('Spline partition/derivative zero '+str(w),abs(math.fsum(a[0] for a in actual)-1)<1e-15 and abs(math.fsum(a[1] for a in actual))<1e-15)

# Degenerate manufactured sites exercise only identity/format algebra, not a box model.
expected=[dict(atom=i,type=402,xyz=[0.,0.,0.],alpha=.837) for i in range(1,11959)]
theta=np.asarray([R.spline5(1.)]*3)
sites=[dict(v,igrid=[26]*3,theta=theta) for v in expected]
s=dict(sites=sites,shape=(64,64,64),p=5,ewald=.4,fft='FFTW',chunks=[1,1,1,64,64,64,2,2,4],volume=49.324**3,trans=np.eye(3)*(64/49.324))
check('Manufactured identity/coordinate/alpha/spline replay accepted',R.check_sites(s,expected,[49.324]*3+[90.]*3)['passed'])
for field,value in [('atom',2),('type',403),('alpha',.838),('xyz',[1e-8,0.,0.]),('igrid',[27,26,26])]:
    bad=dict(s,sites=list(sites));bad['sites'][0]=dict(sites[0]);bad['sites'][0][field]=value
    rejects('Corrupt site field rejected '+field,lambda:R.check_sites(bad,expected,[49.324]*3+[90.]*3))
bad=dict(s,sites=list(sites));bad['sites'][0]=dict(sites[0],theta=theta.copy());bad['sites'][0]['theta'][0,0,0]+=.001
rejects('Corrupt spline rejected',lambda:R.check_sites(bad,expected,[49.324]*3+[90.]*3))
bad=dict(s,volume=1.)
rejects('Corrupt cell volume rejected',lambda:R.check_sites(bad,expected,[49.324]*3+[90.]*3))
with tempfile.TemporaryDirectory() as folder:
    p=Path(folder)/'artifact';p.write_text('original');authority=R.pin(p)
    check('Exact pin accepted',R.verify(authority)==p.resolve())
    p.write_text('corrupted')
    rejects('Changed file pin rejected',lambda:R.verify(authority))

out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
result=dict(status='passed_operator_runner_manufactured_controls',implementation=R.pin(source),controls=R.pin(__file__),checks=checks,checks_passed=len(checks),native_calls=0,link_calls=0,actual_box_evaluations=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1))
(out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(result=R.pin(out/'RESULT.json'),checks=len(checks))))
