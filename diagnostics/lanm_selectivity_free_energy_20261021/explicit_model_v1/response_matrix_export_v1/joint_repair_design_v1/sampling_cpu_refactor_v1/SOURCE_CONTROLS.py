#!/usr/bin/env python3
"""Small synthetic contraction/dispatch checks, not native qualification."""
import hashlib,json,math,os,re,time
from pathlib import Path
assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[5]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1'
NEW=W/'sampling_cpu_refactor_v1/source_v1';OLD=W/'native_v1/source_v1'
def pin(p):p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def subroutine(text,name):
 return re.search(r'(?ms)^      subroutine '+name+r'\b.*?^      end\n',text)[0]
def shapes(vec):
 r=np.linalg.norm(vec,axis=1);t=2*r;p=1+t+t*t/3
 f=(p*np.exp(-t))**2/((1+5+25/3)*math.exp(-5))**2
 fr=-(4/3)*f*t*(1+t)/p
 df=np.zeros_like(vec);np.divide(fr[:,None]*vec,r[:,None],out=df,where=r[:,None]!=0)
 kappa=.27;z=kappa*r**3;phi=np.empty(len(r));pr=np.empty(len(r))
 for i in range(len(r)):
  if z[i]<1e-4:
   phi[i]=kappa*(1-z[i]/2+z[i]**2/6-z[i]**3/24+z[i]**4/120-z[i]**5/720)
   pr[i]=kappa*kappa*r[i]**2*(-1.5+z[i]-3*z[i]**2/8+z[i]**3/10-z[i]**4/48)
  else:
   si=1-math.exp(-z[i]);phi[i]=si/r[i]**3;pr[i]=3*(z[i]*math.exp(-z[i])-si)/r[i]**4
 return r,f,df,phi,pr

def geometry(coords,coef,linear):
 vec=coords[1:]-coords[0];r,f,df,phi,pr=shapes(vec);n=len(r);a,b,c,d=coef
 if linear:
  prefix=np.r_[0.,np.cumsum(f)];suffix=np.r_[np.cumsum(f[::-1])[::-1],0.]
  other=prefix[:-1]+suffix[1:];s=prefix[-1];q=sum(2*f[i]*prefix[i] for i in range(n))
 else:
  other=np.array([sum(f[j] for j in range(n) if i!=j) for i in range(n)]);s=sum(f);q=sum(2*f[i]*f[j] for i in range(n) for j in range(i+1,n))
 eta=other/(1+other);e=q/(1+s*s);g=2*s/(1+s);sig=g*(c+d*e);alpha=2.03
 beta=math.expm1(sig)/alpha
 de=(2*other*(1+s*s)-2*s*q)/(1+s*s)**2
 db=math.exp(sig)/alpha*(2*(c+d*e)/(1+s)**2+g*d*de)[:,None]*df
 pref=a+b*eta;h=3*f[:,None]*pref[:,None]*phi[:,None]*vec
 return vec,r,f,df,phi,pr,other,pref,beta,db,h

def dense(coords,coef,u,v):
 vec,r,f,df,ph,pr,other,pref,beta,db,h=geometry(coords,coef,False);n=len(r);dh=np.zeros((n+1,3,n+1,3));bg=np.zeros((n+1,3));bg[1:]=db;bg[0]=-np.sum(db,axis=0)
 for i in range(n):
  for j in range(n):
   bij=pref[i] if i==j else f[i]*coef[1]/(1+other[i])**2
   jac=3*ph[i]*bij*np.outer(vec[i],df[j])
   if i==j:
    jac+=3*f[i]*pref[i]*ph[i]*np.eye(3)
    if r[i]>0:jac+=3*f[i]*pref[i]*pr[i]/r[i]*np.outer(vec[i],vec[i])
   dh[i+1,:,j+1,:]=jac;dh[i+1,:,0,:]-=jac
 grad=.5*(np.dot(u[0],v[0])*bg-np.einsum('ik,ikjl->jl',u+v,dh))
 return h,beta,grad

def contracted(coords,coef,u,v):
 vec,r,f,df,ph,pr,other,pref,beta,db,h=geometry(coords,coef,True);w=u[1:]+v[1:];dot=np.einsum('ij,ij->i',vec,w);z=3*ph*dot;term=coef[1]*f*z/(1+other)**2
 prefix=np.r_[0.,np.cumsum(term)];suffix=np.r_[np.cumsum(term[::-1])[::-1],0.];ci=pref*z+prefix[:-1]+suffix[1:]
 prr=np.divide(pr,r,out=np.zeros_like(pr),where=r!=0)
 gg=.5*(np.dot(u[0],v[0])*db-ci[:,None]*df-3*f[:,None]*pref[:,None]*(ph[:,None]*w+prr[:,None]*dot[:,None]*vec))
 grad=np.vstack((-np.sum(gg,axis=0),gg));return h,beta,grad

def energy(coords,coef,u,v):
 *_,beta,db,h=geometry(coords,coef,True)
 return .5*(beta*np.dot(u[0],v[0])-np.sum((u[1:]+v[1:])*h))

def main(out):
 tick=time.monotonic();cpu=time.process_time();checks=[];maxerr=maxfd=0.
 def check(name,ok,detail=None):
  assert bool(ok),(name,detail);checks.append({'name':name,'passed':True,'detail':detail})
 def close(name,a,b):
  nonlocal maxerr
  err=float(np.max(abs(np.asarray(a)-b)));scale=max(1.,float(np.max(abs(a))),float(np.max(abs(b))));bound=512*np.size(a)*np.finfo(float).eps*scale;maxerr=max(maxerr,err);check(name,err<=bound,{'error':err,'arithmetic_bound':bound})
 fixtures=[[[0,0,0],[2.5,.2,-.1]],[[.2,-.1,.3],[2.5,.2,-.1],[-1.7,1.3,.7]],[[0,0,0],[2.4,.2,.3],[-1.2,2.1,.4],[.3,-1.3,2.5]],[[0,0,0],[0,0,0],[1e-7,-2e-7,3e-7],[35.,0.,0.]]]
 profiles=[(.05,-.04,0.,0.),(0.,0.,-.08,.06),(.05,-.04,.08,-.06)]
 for fi,raw in enumerate(fixtures):
  x=np.array(raw,dtype=float);u=np.arange(x.size,dtype=float).reshape(x.shape)*.013-.04;v=u[::-1]*-.37+.02
  for pi,c in enumerate(profiles):
   a=dense(x,c,u,v);b=contracted(x,c,u,v)
   for k in range(3):close(f'fixture{fi}_profile{pi}_dense_contracted_{k}',a[k],b[k])
   eps=1e-5;fd=np.zeros_like(x)
   for atom in range(len(x)):
    for axis in range(3):
     delta=np.zeros_like(x);delta[atom,axis]=eps;fd[atom,axis]=(energy(x+delta,c,u,v)-energy(x-delta,c,u,v))/(2*eps)
   err=float(np.max(abs(fd-b[2])));maxfd=max(maxfd,err);check(f'fixture{fi}_profile{pi}_fixed_dual_functional_FD',err<=2e-6,err)
   close(f'fixture{fi}_profile{pi}_net_force',np.sum(b[2],axis=0),np.zeros(3))
 x=np.asarray(fixtures[2],dtype=float);u=np.arange(x.size).reshape(x.shape)*.013-.04;v=u[::-1]*-.37+.02;c=profiles[-1];R=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
 h,beta,g=contracted(x,c,u,v);hr,br,gr=contracted(x@R.T+np.array([5.,-3.,7.]),c,u@R.T,v@R.T)
 close('rigid_rotation_translation_h',hr,h@R.T);close('rigid_rotation_translation_beta',br,beta);close('rigid_rotation_translation_gradient',gr,g@R.T)
 perm=[0,3,1,2];hp,bp,gp=contracted(x[perm],c,u[perm],v[perm]);close('water_permutation_gradient',gp,g[perm]);close('water_permutation_beta',bp,beta)
 vir=.5*((x-x[0]).T@g+g.T@(x-x[0]));strain=np.array([[.13,.07,-.02],[.07,-.11,.03],[-.02,.03,.05]]);eps=1e-5
 derivative=(energy(x@(np.eye(3)+eps*strain).T,c,u,v)-energy(x@(np.eye(3)-eps*strain).T,c,u,v))/(2*eps)
 check('symmetric_virial_fixed_dual_strain',abs(derivative-np.sum(vir*strain))<=2e-6,abs(derivative-np.sum(vir*strain)))
 old=(OLD/'source/induce.f').read_text();new=(NEW/'source/induce.f').read_text()
 oldsetup=subroutine(old,'localresponsesetup');newsetup=subroutine(new,'localresponsesetup')
 check('exact_mode_and_topology_guards',oldsetup[oldsetup.index('      if (.not.lractive) return'):oldsetup.index('      if (allocated(lrh))')]==newsetup[newsetup.index('      if (.not.lractive) return'):newsetup.index('      if (allocated(lrh))')])
 check('exact_absent_zero_bypass',newsetup.index('      if (.not.lractive) return')<newsetup.index('      if (allocated(lrh))'))
 check('dense_SPD_unchanged_except_shared_call',subroutine(old,'localresponsespectrum').replace('call ufield0a(','call ufield0(')==subroutine(new,'localresponsespectrum'))
 check('exporter_unchanged_except_shared_calls',subroutine(old,'exportresponsematrix').replace('call dfield0a ','call dfield0 ').replace('call ufield0a ','call ufield0 ')==subroutine(new,'exportresponsematrix'))
 for kind in ('d','u'):
  wrapper=subroutine(new,kind+'field0');check(kind+'_shared_dispatch_order',wrapper.index('if (use_ewald)')<wrapper.index('else if (use_mlist)')<wrapper.index('if (lractive)'))
  check(kind+'_inner_kernel_no_second_hook','lractive' not in subroutine(new,kind+'field0a'))
 check('no_dense_geometry_jacobian','lrdh' not in new and 'lrdh' not in (NEW/'source/polar.f').read_text())
 check('energy_frame_moment_routes_byte_identical',all((OLD/'source'/name).read_bytes()==(NEW/'source'/name).read_bytes() for name in ('epolar.f','epolar1.f','epolar3.f','moments.f','torque.f')))
 check('fixed_format_executable_columns',all(len(line)<=72 or not line.strip() or line[0] in 'cC*!' for name in ('induce.f','polar.f') for line in (NEW/'source'/name).read_text().splitlines()))
 result={'status':'passed_source_contraction_controls','implementation':pin(__file__),'source_manifest':pin(NEW/'MANIFEST.json'),'checks':checks,'passed':len(checks),'maximum_dense_contracted_absolute_difference':maxerr,'maximum_synthetic_FD_absolute_difference':maxfd,'native_calls':0,'builds':0,'optimizer_calls':0,'interpretation':'Synthetic algebra and source structure only; all native qualification remains pending.','cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
 out=Path(out);out.mkdir(parents=True,exist_ok=False);(out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'result':pin(out/'RESULT.json'),'passed':len(checks),'max_dense_difference':maxerr,'max_FD_difference':maxfd}))
if __name__=='__main__':main(__import__('sys').argv[1])
