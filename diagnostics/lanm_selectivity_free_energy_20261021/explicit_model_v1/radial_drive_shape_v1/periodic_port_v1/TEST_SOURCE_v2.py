#!/usr/bin/env python3
"""Small source/manufactured checks; never imports or calls a native engine."""
import importlib.util
import json
import math
from pathlib import Path
import struct
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('periodic_radial_prepare',HERE/'PREPARE_SOURCE_v2.py')
P=importlib.util.module_from_spec(spec);spec.loader.exec_module(P)

def scalar(r,kappa):
    if r>=21:return 0.,0.,0.,0.
    t=2*r;p=1+t+t*t/3;scale=1/((1+5+25/3)*math.exp(-5))**2
    f=scale*(p*math.exp(-t))**2;fr=-(4/3)*f*t*(1+t)/p
    if r>20:
        u=r-20;s=(1-u)**3*(1+3*u+6*u*u);sr=-30*u*u*(1-u)**2
        fr=fr*s+f*sr;f*=s
    z=kappa*r**3
    if z<1e-4:
        ph=kappa*(1-z/2+z*z/6-z**3/24+z**4/120-z**5/720)
        pr=kappa*kappa*r*r*(-1.5+z-3*z*z/8+z**3/10-z**4/48)
    else:
        si=-math.expm1(-z);ph=si/r**3;pr=3*(z*math.exp(-z)-si)/r**4
    return f,fr,ph,pr

def field(R,e,kappa):
    if e==0:return [0.,0.,0.]
    r=math.hypot(*R)
    if r>=21:return [0.,0.,0.]
    f,fr,ph,pr=scalar(r,kappa);psi=(r*r-6.25)/6.25
    return [3*e*f*ph*psi*x for x in R]

def contraction(R,e,kappa,u,v):
    """Independent B/B' Jacobian contraction, not native factored syntax."""
    r=math.hypot(*R);w=[a+b for a,b in zip(u,v)]
    if e==0 or r>=21:return [0.,0.,0.]
    f,fr,ph,pr=scalar(r,kappa);psi=(r*r-6.25)/6.25
    B=3*f*ph*psi;Bp=3*((fr*ph+f*pr)*psi+2*r*f*ph/6.25)
    dot=math.fsum(x*y for x,y in zip(w,R))
    return [-e*(B*w[j]+(Bp/r if r else 0)*dot*R[j]) for j in range(3)]

def native_factored(R,e,kappa,u,v):
    r=math.hypot(*R);w=[a+b for a,b in zip(u,v)]
    if e==0 or r>=21:return [0.,0.,0.]
    f,fr,ph,pr=scalar(r,kappa);psi=(r*r-6.25)/6.25
    dot=math.fsum(x*y for x,y in zip(w,R));po=pr/r if r else 0
    df=[fr*x/r if r else 0 for x in R]
    return [-3*e*(ph*psi*dot*df[j]+f*(ph*psi*w[j]
            +(psi*po+2*ph/6.25)*dot*R[j])) for j in range(3)]

def energy(R,e,k,u,v):
    return -math.fsum((a+b)*h for a,b,h in zip(u,v,field(R,e,k)))

def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    start,cpu=time.monotonic(),time.process_time();checks=[]
    result={'status':'started','checks':checks,'error':None,'native_calls':0,
            'compilations':0,'solvers':0,'factorizations':0,
            'scope':'Source identity plus manufactured explicit h derivative only; no self-consistent force admission'}
    def check(ok,name,**data):
        checks.append({'name':name,'passed':bool(ok),**data})
        if not ok:raise AssertionError(name)
    def near(a,b,rel,name,floor=1e-290):
        error=max(abs(x-y) for x,y in zip(a,b));scale=max(floor,*map(abs,a),*map(abs,b))
        check(error<=rel*scale,name,relative_error=error/scale,bound=rel)
    try:
        inp=json.loads((HERE/'INPUTS_v2.json').read_text());src=Path(inp['source_output'])
        manifest=json.loads((src/'MANIFEST.json').read_text())
        old={n:P.verify(inp['authorities']['parent_'+n[:-2]]).read_text() for n in ('polar.f','induce.f')}
        actual={n:(src/'source'/n).read_text() for n in old}
        regenerated,edits=P.transform(old,(HERE/'PERIODIC_RADIAL.f').read_text())
        check(actual==regenerated,'Exact generator output and reversible two-file overlay')
        check(manifest['changed_files']==['source/induce.f','source/polar.f'] or
              sorted(manifest['changed_files'])==['source/induce.f','source/polar.f'],'Only two native files changed')
        for name in ('localresponsesetup','localresponsegeometry','localresponsegradient',
                     'periodicresponsebackendguard','periodicresponsestability'):
            check(P.routine(old['induce.f'],name)==P.routine(actual['induce.f'],name),'Unchanged '+name)
        for name,nextname in [('dfield0 (field,fieldp)','ufield0 (field,fieldp)'),
                              ('ufield0 (field,fieldp)','localresponsesetup')]:
            def body(t):return t.split('      subroutine '+name+'\n',1)[1].split('      subroutine '+nextname+'\n',1)[0]
            check(body(old['induce.f'])==body(actual['induce.f']),'Unchanged common '+name)
        parent=json.loads(P.verify(inp['authorities']['parent_manifest']).read_text())
        mapping={str(Path(p['path']).relative_to(src)):p['sha256'] for p in manifest['source_outputs']}
        for ref in parent['source_outputs']:
            rel=str(Path(ref['path']).relative_to(parent['source_root']))
            if rel not in ('source/induce.f','source/polar.f'):
                if mapping[rel]!=ref['sha256']:raise AssertionError('Unchanged parent '+rel)
        check(True,'All613 other frozen source entries retain parent hashes')
        setup=P.routine(actual['induce.f'],'periodicresponsesetup')
        check(setup.index('lrpshape=0.0d0')<setup.index('if (.not.handled) return'),'Reset new state even on vacuum/absent dispatch')
        check(setup.index('call periodicradialsetup')<setup.index('if (.not.lractive) return'),'e-only activation precedes inactive return')
        check("trim(keyword).eq.'POLAR-RADIAL-DRIVE'" in setup,'New keyword enters BIO parser independently of local records')
        helper=(HERE/'PERIODIC_RADIAL.f').read_text()
        for guard in ("ia.ne.900.and.ia.ne.901","ib.ne.402","ieee_is_finite(e)",
                      'seen(slot) .or. abs(e).gt.1.0d0','character*240 record,keyword,string,extra',
                      'ia,ib,e,extra','if (ios.ge.0) then'):
            check(guard in helper,'Parser source guard '+guard)
        for routine,call in [('periodicresponsegeometry','periodicradialfield'),
                             ('periodicresponsegradient','periodicradialgradient')]:
            body=P.routine(actual['induce.f'],routine)
            check(body.count('call '+call)==1 and 'if (lrpshape(slot).ne.0.0d0)' in body,'Exact zero call bypass '+routine)
        grad=P.routine(actual['induce.f'],'periodicresponsegradient')
        check(grad.index('call periodicradialgradient')<grad.index('gg=fac*gg')<grad.index('dep(:,oi)=dep(:,oi)+gg'),'Single native unit/reaction/virial path')
        check('ow(nw)' in helper and 'lrdh' not in helper,'Dynamic water extent; no dense derivative allocation')
        check('lrcoef' not in helper and 'lrbeta' not in helper,'Radial helper does not change old coefficients or H')
        u=[.31,-.24,.09];v=[-.12,.07,.41];e=.37;k=.23
        vectors=[[0.,0.,0.],[.0003,-.0004,.0005],[1.2,-.7,.9],[2.5,0.,0.],
                 [2.8,0.,0.],[20.,0.,0.],[20.5,0.,0.]]
        for R in vectors:
            analytic=contraction(R,e,k,u,v);factored=native_factored(R,e,k,u,v)
            near(analytic,factored,8e-15,'Direct-J versus native-factored dual contraction '+str(R))
            h=2e-6;fd=[]
            for j in range(3):
                plus=R[:];minus=R[:];plus[j]+=h;minus[j]-=h
                fd.append((energy(plus,e,k,u,v)-energy(minus,e,k,u,v))/(2*h))
            near(analytic,fd,3e-6,'Manufactured Cartesian explicit-energy FD '+str(R))
        contact=contraction([0.,0.,0.],e,k,u,v);f=scalar(0.,k)[0]
        near(contact,[3*e*f*k*(a+b) for a,b in zip(u,v)],2e-15,'Finite exact contact Jacobian')
        check(field([2.5,0.,0.],e,k)==[0.,0.,0.],'Exact r0 zero drive')
        check(abs(contraction([2.5,0.,0.],e,k,u,v)[0])>0,'Nonzero r0 radial derivative')
        for R in ([0.,0.,0.],[2.5,0.,0.],[20.5,0.,0.],[21.,0.,0.],[22.,0.,0.]):
            near(field(R,-e,k),[-x for x in field(R,e,k)],0,'Odd e field '+str(R))
            check(field(R,0.,k)==[0.,0.,0.] and contraction(R,0.,k,u,v)==[0.,0.,0.],'Zero-e exact scalar bypass '+str(R))
        before=[-0.,1.,-2.];after=before[:] # Native branch performs no added arithmetic.
        check(struct.pack('3d',*before)==struct.pack('3d',*after),'Manufactured signed-zero preservation of no-add route')
        for R in ([21.,0.,0.],[22.,0.,0.]):
            check(field(R,e,k)==[0.,0.,0.] and contraction(R,e,k,u,v)==[0.,0.,0.],'Exact compact support '+str(R))
        left=[abs(field([21-h,0.,0.],e,k)[0]) for h in (1e-3,5e-4)]
        leftg=[abs(contraction([21-h,0.,0.],e,k,u,v)[0]) for h in (1e-3,5e-4)]
        check(.10<left[1]/left[0]<.15,'C2 switch field approaches cubic zero',ratio=left[1]/left[0])
        check(.20<leftg[1]/leftg[0]<.30,'C2 switch first derivative approaches quadratic zero',ratio=leftg[1]/leftg[0])
        R=[1.2,-.7,.9];g=contraction(R,e,k,u,v)
        check(max(abs(a-b) for a,b in zip(g,contraction(R,e,k,u,u)))>1e-4,'Unequal duals detect incorrect u-only contraction')
        near([a+(-a) for a in g],[0.,0.,0.],0,'Equal/opposite metal reaction')
        S=[[.3,.1,-.2],[.1,-.4,.05],[-.2,.05,.2]]
        SR=[math.fsum(S[i][j]*R[j] for j in range(3)) for i in range(3)]
        V=[[.5*(R[i]*g[j]+R[j]*g[i]) for j in range(3)] for i in range(3)]
        projected=math.fsum(V[i][j]*S[i][j] for i in range(3) for j in range(3))
        h=1e-6;plus=[r+h*x for r,x in zip(R,SR)];minus=[r-h*x for r,x in zip(R,SR)]
        near([projected],[(energy(plus,e,k,u,v)-energy(minus,e,k,u,v))/(2*h)],2e-7,'Symmetric image virial versus affine explicit-energy FD')
        r0=[1.25,-.75,.5];shifted=[r0[0]+44,r0[1]-44,r0[2]+88]
        wrapped=[x-44*round(x/44) for x in shifted]
        check(wrapped==r0 and field(wrapped,e,k)==field(r0,e,k),'Equivalent orthogonal image representative')
        check(22>21 and field([22.,0.,0.],e,k)==[0.,0.,0.],'Nearest-image ties outside compact support at Lmin44')
        result.update(status='passed_periodic_radial_source_manufactured_controls',
                      implementation=P.pin(__file__),source_manifest=P.pin(src/'MANIFEST.json'),
                      tested_parameter_values={'e':e,'kappa':k,'u':u,'v':v})
    except Exception as exc:
        result.update(status='failed_periodic_radial_source_manufactured_controls',error=repr(exc))
        raise
    finally:
        result['cost']={'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'threads':1}
        P.save(out/'RESULT.json',result)
    print(json.dumps({'result':P.pin(out/'RESULT.json'),'checks':len(checks),'cost':result['cost']}))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);run(parser.parse_args().out)
