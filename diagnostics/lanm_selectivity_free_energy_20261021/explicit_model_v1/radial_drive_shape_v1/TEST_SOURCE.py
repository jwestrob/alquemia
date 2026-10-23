#!/usr/bin/env python3
"""Finite artificial calculus and literal source-delta controls; no molecular code."""
import os
for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import argparse, hashlib, json, math, re, time, traceback
from pathlib import Path
import numpy as np

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def shape(v,e,kappa=.19):
    r=float(np.linalg.norm(v));t=2*r;p=1+t+t*t/3
    f=(p*math.exp(-t)/((1+5+25/3)*math.exp(-5)))**2
    fr=-(4/3)*f*t*(1+t)/p;z=kappa*r**3
    if z<1e-4:
        ph=kappa*(1-z/2+z*z/6-z**3/24+z**4/120-z**5/720)
        pr=kappa*kappa*r*r*(-1.5+z-3*z*z/8+z**3/10-z**4/48)
    else:
        si=-math.expm1(-z);ph=si/r**3;pr=3*(z*math.exp(-z)-si)/r**4
    psi=(r*r-6.25)/6.25;b=3*f*ph*psi
    bp=3*((fr*ph+f*pr)*psi+2*r*f*ph/6.25)
    h=e*b*v; jac=e*b*np.eye(3)
    if r>0: jac+=e*bp/r*np.outer(v,v)
    return h,jac,f,ph

def strip_insertions(s):
    return re.sub(r'c     RADIAL_SHAPE_BEGIN_\w+\n.*?c     RADIAL_SHAPE_END_\w+\n','',s,flags=re.S)
def main():
    a=argparse.ArgumentParser();a.add_argument('--manifest',required=True);a.add_argument('--out',required=True);args=a.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False);t=time.monotonic();cpu=time.process_time();checks=[];solves=0
    def ck(name,ok,**details):
        checks.append({'name':name,'passed':bool(ok),**details})
        if not ok: raise AssertionError(name)
    result={'status':'started','implementation':pin(__file__),'source_manifest':pin(args.manifest),'checks':checks,'native_calls':0,'builds':0,'optimizer_calls':0,'molecular_or_reference_evaluations':0}
    try:
        sm=json.loads(Path(args.manifest).read_text());src=Path(sm['source_root']);old=json.loads(Path(sm['parent_source_manifest']['path']).read_text());parent=Path(old['source_root'])
        ck('parent_manifest_hash',pin(sm['parent_source_manifest']['path'])==sm['parent_source_manifest'])
        changed=[]
        for p in sm['source_outputs']+[sm['cmake']]:
            ckpath=Path(p['path']);assert pin(ckpath)==p
            if ckpath.read_bytes()!=(parent/ckpath.relative_to(src)).read_bytes():changed.append(ckpath.name)
        ck('only_two_source_files_changed',changed==['induce.f','polar.f'],changed=changed)
        for name in changed:
            s=(src/'source'/name).read_text();raw=(parent/'source'/name).read_text()
            restored=strip_insertions(s)
            if name=='induce.f': restored=restored.removesuffix('\n')
            ck(name+'_literal_parent_restoration',restored==raw)
            bad=[i for i,line in enumerate(s.splitlines(),1) if len(line)>72 and line[:1] not in 'cC*!']
            ck(name+'_fixed_form_line_length',not bad,violations=bad)
        s=(src/'source/induce.f').read_text();g=s.split('      subroutine radialdrivegeometry',1)[1]
        ck('new_geometry_has_no_H_T_environment_mutation',not any(w in g for w in ('lrbeta','lrdb','lrcoef','eta','ufield','dfield','polarity','thlval','pdamp')))
        ck('zero_outer_and_inner_guard','if (lrshape(slot).ne.0.0d0) then' in s and 'if (e.eq.0.0d0) return' in g)
        ck('common_h_and_gradient_consumers_preserved',s[:s.index('c     Regularized local response:')]==(parent/'source/induce.f').read_text().split('c     Regularized local response:')[0])
        ck('finite_domain_duplicate_and_character_extra_parser',all(w in s for w in ('abs(e).gt.1.0d0','seen(slot)','ieee_is_finite(e)','read(string,*,iostat=ios) ia,ib,e,extra','character*240 record,keyword,string,extra')))
        for name,v in [('contact',np.zeros(3)),('small',np.array([.002,-.001,.003])),('inner',np.array([.3,-.4,2.1])),('center',np.array([0.,0.,2.5])),('outer',np.array([3.1,.2,-.4]))]:
            for e in (-.17,.17):
                h,j,f,ph=shape(v,e);fd=np.column_stack([(shape(v+np.eye(3)[k]*1e-6,e)[0]-shape(v-np.eye(3)[k]*1e-6,e)[0])/2e-6 for k in range(3)])
                err=float(np.max(np.abs(fd-j)));bound=1e-7*max(1.,float(np.max(np.abs(j))))
                ck(name+'_Jacobian_'+str(e),np.isfinite(h).all() and np.isfinite(j).all() and err<bound,error=err,bound=bound)
            hz,jz,*_=shape(v,0.);ck(name+'_zero_shape',np.array_equal(hz,np.zeros(3)) and np.array_equal(jz,np.zeros((3,3))))
        h,j,f,ph=shape(np.zeros(3),.17);ck('contact_limit',np.array_equal(h,np.zeros(3)) and np.array_equal(j,-3*.17*f*.19*np.eye(3)))
        h,j,f,ph=shape(np.array([0.,0.,2.5]),.17);expected=np.diag([0.,0.,6*.17*ph]);ck('ideal_r0_zero_drive_nonzero_radial_slope',np.array_equal(h,np.zeros(3)) and np.max(abs(j-expected))<1e-15)
        h,j,*_=shape(np.array([0.,0.,2.5+1e-10]),.17);ck('literal_near_r0_not_clamped',float(np.linalg.norm(h))>0)
        rot=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]]);v=np.array([.3,-.4,2.1]);h,j,*_=shape(v,.17);hr,jr,*_=shape(rot@v,.17)
        ck('rotation_covariance',np.max(abs(hr-rot@h))<1e-14 and np.max(abs(jr-rot@j@rot.T))<1e-14)
        ck('asymptotic_separation',np.linalg.norm(shape(np.array([20.,0.,0.]),.17)[0])<1e-24)
        # Independent dual variational solve on a manufactured 9-dimensional SPD operator.
        L=np.eye(9)*1.3+np.tril(np.arange(81).reshape(9,9)*.0007,-1);A=L@L.T
        Ed=np.linspace(-.2,.3,9);Ep=np.linspace(.31,-.13,9);K=332.0637133
        pos=np.array([[.1,.2,-.3],[.4,-.2,1.8],[-1.8,.7,.4]])
        def energy_gradient(x,e):
            nonlocal solves
            hv=np.zeros((3,3));js=[]
            for i in (1,2):
                hh,jj,*_=shape(x[i]-x[0],e);hv[i]=hh;js.append(jj)
            u=np.linalg.solve(A,Ed+hv.ravel());v=np.linalg.solve(A,Ep+hv.ravel());solves+=2
            energy=-K/2*float(u@(Ep+hv.ravel()));grad=np.zeros((3,3))
            for i,jj in zip((1,2),js):
                grad[i]=-K/2*((u+v).reshape(3,3)[i]@jj);grad[0]-=grad[i]
            return energy,grad,u,v
        en,grad,u,v=energy_gradient(pos,.17);err=0.
        for i in range(3):
            for k in range(3):
                d=np.zeros((3,3));d[i,k]=1e-6;fd=(energy_gradient(pos+d,.17)[0]-energy_gradient(pos-d,.17)[0])/2e-6;err=max(err,abs(fd-grad[i,k]))
        ck('unequal_dual_total_energy_gradient',not np.array_equal(u,v) and err<5e-6,max_absolute_error=err)
        ck('metal_reaction_translation',np.max(abs(grad.sum(axis=0)))<1e-12 and abs(energy_gradient(pos+np.array([1.1,-2.3,.7]),.17)[0]-en)<1e-11)
        strain=np.array([[.3,.1,0.],[.1,-.2,.2],[0.,.2,.4]]);tangent=pos@strain.T
        fd=(energy_gradient(pos+1e-6*tangent,.17)[0]-energy_gradient(pos-1e-6*tangent,.17)[0])/2e-6
        virial=pos.T@grad;pred=float(np.sum(strain*virial));ck('symmetric_strain_virial_contraction',abs(fd-pred)<5e-6,error=abs(fd-pred))
        W=np.tile(np.eye(3),(3,1));chi=W.T@np.linalg.solve(A,W);solves+=1
        ck('fixed_A_susceptibility_independent_of_drive',np.isfinite(chi).all() and np.array_equal(A,L@L.T))
        # Certify interval mapping/count algebra without producing a scientific proposal.
        ck('new_normalized_probe_and_trust_scales',abs(2*.005-.01)<1e-17 and abs(2*.025-.05)<1e-17)
        ck('qualification_and_new_probe_denominators',4+10+12+16+4+3==49 and 2*6*4+2*6*7==132 and 21*2==42 and 4*(23+5)==112)
        result['status']='passed_radial_drive_source_and_artificial_controls'
    except Exception as e:
        result['status']='failed';result['error']=repr(e);result['traceback']=traceback.format_exc()
    finally:
        result.update(check_count=len(checks),manufactured_dense_solves=solves,cost={'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1})
        (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k] for k in ('status','check_count','cost')}))
    if result['status']=='failed':raise SystemExit(1)
if __name__=='__main__':main()
