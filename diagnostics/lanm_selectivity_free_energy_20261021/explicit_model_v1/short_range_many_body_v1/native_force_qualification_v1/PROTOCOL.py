"""Strict saved/live text grammar for the one-shot native force bridge."""
import math
from pathlib import Path

MODE=['MODE','VACUUM','MUTUAL','THOLE','DIELECTRIC_1','ZERO_FIELD','ALL_ACTIVE','NO_VIRTUAL_SITES']
def floats(xs):
    a=[float(x.replace('D','E')) for x in xs];assert all(math.isfinite(x) for x in a),'nonfinite native token';return a
def lines(p):return [s.split() for s in Path(p).read_text().splitlines() if s.strip()]
def baseline(path,slot,fixture):
    raw=lines(path);it=iter(raw)
    assert next(it)==['SRMB_BASELINE_V1'];t=next(it);assert t[0]=='COUNTS';n,ns,dim=map(int,t[1:]);assert n==ns==len(slot['geometry']['xyz_A']) and dim==3*n and 1<=n<=28
    t=next(it);assert t[0]=='CONSTANTS' and len(t)==5;K,dielec,D,B=floats(t[1:]);assert dielec==1 and K==fixture['native_operator']['K_kcal_A_per_mol_e2']
    t=next(it);assert t[0]=='OPERATOR_CHECKS' and len(t)==3;sym,dual=floats(t[1:]);assert sym>=0 and dual>=0
    assert next(it)==['PERMANENT_ORIGIN_A','0','0','0'];t=next(it);assert t[0]=='PERMANENT_MU_EA' and len(t)==4;mu=floats(t[1:]);t=next(it);assert t[0]=='CHARGE' and len(t)==2;charge=floats(t[1:])[0]
    assert next(it)==MODE;t=next(it);assert t[0]=='FLAGS_DPEQUAL_PCGGUESS_TCGGUESS' and len(t)==4 and set(t[1:])<={'T','F'};flags=t[1:]
    assert next(it)==['CALLS','1',str(3*n),'0']
    sites=[];A=[];ed=[];ep=[];g=slot['geometry'];types=g.get('types',g.get('tinker_types'));op=fixture['native_operator']
    for ii in range(n):
        t=next(it);assert t[0]=='SITE' and len(t)==8;pos,atom,typ=map(int,t[1:4]);alpha,x,y,z=floats(t[4:]);assert pos==ii+1 and atom==op['site_atom_indices'][ii]+1 and typ==types[atom-1] and alpha==op['alpha_A3'][ii] and [x,y,z]==g['xyz_A'][atom-1]
        sites.append({'position':pos,'atom':atom,'type':typ,'alpha':alpha,'xyz':[x,y,z]})
        for j in range(3):
            k=3*ii+j+1;t=next(it);assert t[0]=='VECTOR' and len(t)==4 and int(t[1])==k;de,pe=floats(t[2:]);ed.append(de);ep.append(pe)
            t=next(it);assert t[0]=='A' and len(t)==dim+2 and int(t[1])==k;A.append(floats(t[2:]))
    assert next(it)==['END_SRMB_BASELINE_V1'] and list(it)==[]
    return {'n':n,'sites':sites,'A':A,'Ed':ed,'Ep':ep,'K':K,'dielec':dielec,'debye':D,'bohr':B,'permanent_mu_eA':mu,'charge_e':charge,'reported_symmetry_error':sym,'reported_dual_operator_error':dual,'native_flags':dict(zip(('dpequal','pcgguess','tcgguess'),flags)),'calls':{'direct_initializations':1,'mutual_operator_applications':dim,'induced_solves':0}}
def write_uv(path,b,u,v):
    def fmt(x):return format(float(x),'.17g')
    s=['SRMB_UV_V1',f"COUNTS {b['n']} {len(b['sites'])}"]
    for i,site in enumerate(b['sites']):
        s.append('SITE '+ ' '.join(map(str,[site['position'],site['atom'],site['type']]))+' '+ ' '.join(fmt(x) for x in [site['alpha']]+site['xyz']))
        for tag,a in [('ED',b['Ed']),('EP',b['Ep']),('U',u),('V',v)]:s.append(tag+' '+ ' '.join(fmt(x) for x in a[3*i:3*i+3]))
    s.append('END_SRMB_UV_V1');p=Path(path);assert not p.exists();p.write_text('\n'.join(s)+'\n')
def result(path,b,u,v):
    it=iter(lines(path));n=b['n'];assert next(it)==['SRMB_FIXED_UV_GRADIENT_V1'];assert next(it)==['COUNTS',str(n),str(n),str(3*n)];assert next(it)==['CALLS','1',str(3*n),'1','0'];e={}
    for tag in ('OTHER_ENERGY','RAW_REDUCED_POLARIZATION','RAW_NATIVE_ENERGY'):
        t=next(it);assert t[0]==tag and len(t)==2;e[tag]=floats(t[1:])[0]
    assert next(it)==['SUPPLIED_UV_UNCHANGED','T']
    for i,site in enumerate(b['sites']):
        t=next(it);assert t[0]=='UV' and len(t)==9 and list(map(int,t[1:3]))==[i+1,site['atom']];uv=floats(t[3:]);assert uv==u[3*i:3*i+3]+v[3*i:3*i+3],'native changed exact supplied dipoles'
    gradient=[];polar=[];vir=[]
    for i in range(n):
        t=next(it);assert t[0]=='BASELINE_GRADIENT' and len(t)==5 and int(t[1])==i+1;gradient.append(floats(t[2:]));t=next(it);assert t[0]=='POLAR_CONTRACTION_GRADIENT' and len(t)==5 and int(t[1])==i+1;polar.append(floats(t[2:]))
    for j in range(3):
        t=next(it);assert t[0]=='RAW_NATIVE_VIRIAL' and len(t)==5 and int(t[1])==j+1;vir.append(floats(t[2:]))
    assert next(it)==['END_SRMB_FIXED_UV_GRADIENT_V1'] and list(it)==[]
    return {'energy':e,'baseline_gradient':gradient,'polar_contraction_gradient':polar,'raw_native_virial':vir,'u':u,'v':v,'supplied_uv_unchanged':True,'calls':{'direct_initializations':1,'mutual_operator_applications':3*n,'baseline_gradient_calls':1,'induced_solves':0}}
