#!/usr/bin/env python3
"""Bounded synthetic shape/derivative and source controls; no native calls."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
import hashlib,json,math,sys,time
from pathlib import Path
import numpy as np
start=time.monotonic(); cpu=time.process_time(); checks=[]
def check(name,ok,**values):
    checks.append(dict(name=name,passed=bool(ok),**values))
    if not ok: raise ValueError(name)
def shape(v,a,b,c,d):
    r=np.linalg.norm(v,axis=1); t=2*r; P=1+t+t*t/3
    f=(P*np.exp(-t)/((1+5+25/3)*np.exp(-5)))**2
    fp=-(4/3)*f*t*(1+t)/P; df=fp[:,None]*v/r[:,None]
    kappa=.19; z=kappa*r**3; phi=-np.expm1(-z)/r**3
    phip=3*(z*np.exp(-z)+np.expm1(-z))/r**4
    s=f.sum(); other=np.array([sum(f[j] for j in range(len(f)) if j!=i) for i in range(len(f))])
    eta_i=other/(1+other); Q=(f*f).sum(); eta=(s*s-Q)/(1+s*s); g=2*s/(1+s)
    sigma=g*(c+d*eta); beta=np.expm1(sigma)/1.6597
    h=3*f[:,None]*(a+b*eta_i)[:,None]*phi[:,None]*v
    # Chain rule through f_j, independent of native storage layout.
    de=(2*s*(1+Q)-2*f*(1+s*s))/(1+s*s)**2
    db=np.exp(sigma)/1.6597*(2*(c+d*eta)/(1+s)**2+g*d*de)
    jb=db[:,None]*df
    jh=np.zeros((len(f),3,len(f),3))
    for i in range(len(f)):
        for j in range(len(f)):
            dB=(a+b*eta_i[i]) if i==j else f[i]*b/(1+other[i])**2
            jh[i,:,j,:]=3*phi[i]*np.outer(v[i],dB*df[j])
            if i==j: jh[i,:,j,:]+=3*f[i]*(a+b*eta_i[i])*(phi[i]*np.eye(3)+phip[i]/r[i]*np.outer(v[i],v[i]))
    return h,beta,jh,jb,eta_i,eta
for label,v in [('pair',np.array([[.3,-.4,2.4]])),('two_water',np.array([[.3,-.4,2.4],[2.1,.7,-.5]])),('three_water',np.array([[.3,-.4,2.4],[2.1,.7,-.5],[-1.9,1.2,.3]]))]:
    pars=(.05,-.04,.08,-.06); h,beta,jh,jb,eta_i,eta=shape(v,*pars)
    errh=errb=0.
    for j in range(len(v)):
        for k in range(3):
            vp=v.copy();vm=v.copy();vp[j,k]+=1e-6;vm[j,k]-=1e-6
            hp,bp,*_=shape(vp,*pars);hm,bm,*_=shape(vm,*pars)
            errh=max(errh,float(np.max(np.abs((hp-hm)/2e-6-jh[:,:,j,k]))));errb=max(errb,abs((bp-bm)/2e-6-jb[j,k]))
    check(label+'_drive_coordinate_derivative',errh<2e-9,max_absolute_difference=errh)
    check(label+'_stiffness_coordinate_derivative',errb<2e-9,max_absolute_difference=errb)
    if len(v)==1: check('pair_collective_terms_exact_zero',eta==0 and eta_i[0]==0)
    perm=np.arange(len(v))[::-1];hh,bb,*_=shape(v[perm],*pars)
    check(label+'_water_permutation',np.allclose(hh,h[perm],rtol=0,atol=1e-14) and abs(bb-beta)<1e-14)
    rot=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
    hh,bb,*_=shape(v@rot.T,*pars)
    check(label+'_rotation_covariance',np.allclose(hh,h@rot.T,rtol=0,atol=1e-14) and abs(bb-beta)<1e-14)
S=math.log(4)
check('global_stiffness_domain',all(abs(S*(u+v)/4)+abs(S*(u-v)/4)<=S/2+1e-15 for u in (-1,0,1) for v in (-1,0,1)))
for z in (0.,1e-12,1e-8,1e-4):
    k=.19;r=(z/k)**(1/3) if z else 0
    phi=k*(1-z/2+z*z/6-z**3/24+z**4/120-z**5/720)
    pr=k*k*r*r*(-1.5+z-3*z*z/8+z**3/10-z**4/48)
    exact=k if z==0 else -math.expm1(-z)/r**3
    check('regularized_small_z_'+str(z),math.isfinite(phi) and math.isfinite(pr) and abs(phi-exact)<1e-15,phi=phi,phi_radial_derivative=pr)
source=Path(sys.argv[1]); baseline=Path(sys.argv[2]); out=Path(sys.argv[3]);out.mkdir(parents=True,exist_ok=False)
changed=[]
for p in sorted((source/'source').iterdir()):
    old=baseline/'source'/p.name
    if p.read_bytes()!=old.read_bytes():changed.append(p.name)
check('exact_five_changed_files',changed==['epolar1.f','epolar3.f','induce.f','moments.f','polar.f'],changed_files=changed)
for name in changed:
    bad=[i for i,line in enumerate((source/'source'/name).read_text().splitlines(),1) if len(line)>72 and line[0] not in 'cC*!']
    check(name+'_fixed_form_columns',not bad,violations=bad)
s=(source/'source/induce.f').read_text()
check('shared_direct_mutual_hooks',s.count('field(:,1:n) = field(:,1:n) + lrh')==1 and s.count('field(:,i) = field(:,i) - lrbeta(i)*uind(:,i)')==1)
check('original_exporter_text_unchanged',s[s.index('      subroutine responsematrixflag'):s.index('c     Regularized local response:')].strip()==(baseline/'source/induce.f').read_text().split('      subroutine responsematrixflag',1)[1].join(['      subroutine responsematrixflag','']).strip())
result={'status':'passed','checks':checks,'check_count':len(checks),'new_native_calls':0,'new_scientific_points':0,'builds':0,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'implementation':{'path':str(Path(__file__).resolve()),'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},'limitations':'Pure shape/source controls; native compilation, matrix/energy/force/field qualification remains pending.'}
(out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('status','check_count','wall_seconds','new_native_calls')}))
