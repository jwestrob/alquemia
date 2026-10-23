#!/usr/bin/env python3
"""Small source-preservation and synthetic dual-functional derivative controls."""
import hashlib,json,math,re,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=Path.cwd()
WS=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1'
def pin(p):return dict(path=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def read(p):return json.loads(p.read_text())
tick,cpu=time.monotonic(),time.process_time();checks={};errors={}
def ck(name,ok):
    checks[name]=bool(ok)
    if not ok:raise AssertionError(name)
mp=WS/'donor_extension_native_v1/source_v1/MANIFEST.json';m=read(mp);old=read(Path(m['baseline_source_manifest']['path']))
src=Path(m['source_root'])/'source';prior=Path(old['source_root'])/'source'
for item in m['source_outputs']:assert pin(Path(item['path']))==item
ck('exact_three_changed_files',m['changed_files']==['induce.f','polar.f','testgrad.f'])
for name in ('epolar.f','epolar1.f','epolar3.f','moments.f','rotpole.f','torque.f'):
    ck('unchanged_'+name,(src/name).read_bytes()==(prior/name).read_bytes())
a=(prior/'induce.f').read_text();b=(src/'induce.f').read_text()
def routine(text,name):
    h=re.search(r'^      subroutine '+name+r'\b[\s\S]*?^      end$',text,re.M)
    assert h is not None,name
    return h.group()
for name in ('localresponsegeometry','localresponsespectrum','dfield0a','ufield0a','exportresponsematrix'):
    ck('unchanged_'+name,routine(a,name)==routine(b,name))
ck('unchanged_water_setup_after_donor_dispatch',routine(a,'localresponsesetup')==routine(b,'localresponsesetup').replace('      call donorresponsesetup\n      if (lrdonor) return\n',''))
ck('unchanged_water_gradient_after_donor_dispatch',routine(a,'localresponsegradient')==routine(b,'localresponsegradient').replace('      if (lrdonor) then\n         call donorresponsegradient\n         return\n      end if\n',''))
helper=(HERE/'DONOR_RESPONSE.f').read_text()
ck('zero_unmatched_bypass_precedes_mode_guards',helper.index('if (all(drcoef.eq.0.0d0)) return')<helper.index('if (nmet.eq.0 .or. no.eq.0) return')<helper.index('bad=use_bounds'))
ck('both_equivalent_oxygen_sites',"if (type(i).ne.403) cycle" in helper and 'nearest' not in helper)
ck('full_gradient_reporting_after_native_gradient','call gradient (eval,derivs)\n            if (lrdonor)\n     &         call donorresponsetotalgradient (eval,derivs)' in (src/'testgrad.f').read_text())
ck('permanent_report_before_M_induced_addition',b.index('call donorresponsepermanent')<b.index('      subroutine exportresponsematrix') and 'mu=mu+rpole(1,i)*(/x(i),y(i),z(i)/)+rpole(2:4,i)' in helper)

# Synthetic three-site system: original A/Ed/Ep are independent of position.
# Only donor h/H varies, isolating the new conservative correction derivative.
xyz=np.array([[0.,0.,0.],[2.2,.2,-.3],[2.7,-.7,.8]])
alpha=np.array([1.2,1.931,1.931]);K=332.0637133
base=np.diag(np.repeat(1/alpha,3));coupling=np.arange(81,dtype=float).reshape(9,9)/100000
base+=coupling+coupling.T
ed=np.linspace(-.03,.04,9);ep=np.linspace(.02,-.05,9)
def geometry(x,aa,cc,gamma):
    h=np.zeros((3,3));beta=np.zeros(3);dh=np.zeros((3,3,3,3));db=np.zeros((3,3,3))
    fzero=1/(((1+5+25/3)*math.exp(-5))**2)
    for i in (1,2):
        v=x[i]-x[0];r=np.linalg.norm(v);t=2*r;p=1+t+t*t/3
        f=fzero*(p*math.exp(-t))**2;fr=-4/3*f*t*(1+t)/p
        kappa=gamma/math.sqrt(alpha[0]*alpha[i]);z=kappa*r**3
        if z<1e-4:
            phi=kappa*(1-z/2+z*z/6-z**3/24+z**4/120-z**5/720)
            phir=kappa*kappa*r*r*(-1.5+z-3*z*z/8+z**3/10-z**4/48)
        else:
            si=-math.expm1(-z);phi=si/r**3;phir=3*(z*math.exp(-z)-si)/r**4
        df=fr*v/r if r else np.zeros(3);w=2*f/(1+f);wr=2*fr/(1+f)**2
        beta[i]=math.expm1(cc*w)/alpha[i]
        if r:db[i,i]=math.exp(cc*w)*cc*wr*v/(alpha[i]*r)
        db[i,0]=-db[i,i]
        h[i]=aa*3*f*phi*v
        jac=aa*3*(phi*np.outer(v,df)+f*phi*np.eye(3))
        if r:jac+=aa*3*f*phir*np.outer(v,v)/r
        dh[i,i]=jac;dh[i,0]=-jac
    return h,beta,dh,db
def value(x,aa,cc,gamma):
    h,bt,dh,db=geometry(x,aa,cc,gamma);G=base+np.diag(np.repeat(bt,3));d=ed+h.ravel();p=ep+h.ravel()
    u=np.linalg.solve(G,d).reshape(3,3);v=np.linalg.solve(G,p).reshape(3,3)
    grad=np.zeros((3,3))
    for j in range(3):
        for i in range(3):grad[j]+=.5*K*(np.dot(u[i],v[i])*db[i,j]-(u[i]+v[i])@dh[i,j])
    return -.5*K*np.dot(d,np.linalg.solve(G,p)),grad,h,bt,u,v,G
for metal,gamma in [('La',.10825158943928785),('Dy',.139)]:
    for name,aa,cc in [('drive',.1,0),('stiffness',0,.1),('combined',.1,.1)]:
        e,grad,h,bt,u,v,G=value(xyz,aa,cc,gamma);fd=np.zeros_like(xyz);step=1e-5
        for i in range(3):
            for j in range(3):
                plus=xyz.copy();minus=xyz.copy();plus[i,j]+=step;minus[i,j]-=step
                fd[i,j]=(value(plus,aa,cc,gamma)[0]-value(minus,aa,cc,gamma)[0])/(2*step)
        err=float(np.max(np.abs(fd-grad)));errors[metal+'_'+name]=err
        ck(metal+'_'+name+'_dual_gradient_FD',err<2e-8)
        ck(metal+'_'+name+'_metal_reaction',np.max(np.abs(np.sum(grad,axis=0)))<2e-14)
        ck(metal+'_'+name+'_actual_unequal_dual_fields',np.max(np.abs(u-v))>.01)
        moved=geometry(xyz+np.array([2.,-4.,.6]),aa,cc,gamma)
        ck(metal+'_'+name+'_translation',np.max(np.abs(moved[0]-h))<2e-16 and np.max(np.abs(moved[1]-bt))<2e-16)
        hg,bg,_,_=geometry(xyz[[0,2,1]],aa,cc,gamma)
        ck(metal+'_'+name+'_oxygen_permutation',np.array_equal(hg,h[[0,2,1]]) and np.array_equal(bg,bt[[0,2,1]]))
        directions=np.array([0.3,-.4,.5]);directions/=np.linalg.norm(directions);B=np.tile(directions,3)
        chi=B@np.linalg.solve(G,B);basechi=B@np.linalg.solve(base,B)
        ck(metal+'_'+name+'_drive_invariance_stiffness_sign',abs(chi-basechi)<1e-14 if cc==0 else chi<basechi)
zero=geometry(xyz,0,0,.139)
ck('zero_h_beta_derivatives',all(np.count_nonzero(v)==0 for v in zero))
contact=xyz.copy();contact[1]=contact[0]
ck('finite_exact_contact_series',all(np.all(np.isfinite(v)) for v in geometry(contact,.1,.1,.139)))
out=WS/'donor_extension_native_v1/source_controls_v1';out.mkdir(parents=True,exist_ok=False)
result=dict(status='passed',implementation=pin(Path(__file__)),source_manifest=pin(mp),checks=checks,passed=sum(checks.values()),declared=len(checks),synthetic_gradient_errors=errors,
    scope='Static preservation and synthetic dual-functional math only; no compilation, molecular program, reference point or native qualification.',new_native_calls=0,builds=0,
    cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
(out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(pin(out/'RESULT.json')))
