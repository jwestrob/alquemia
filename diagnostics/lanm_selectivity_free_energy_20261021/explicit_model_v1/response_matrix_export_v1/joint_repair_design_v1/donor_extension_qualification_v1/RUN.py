#!/usr/bin/env python3
"""The finite donor numerical panel; root-bound execution or saved admission."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shlex
import shutil
import threading
import time
import numpy as np

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
FILES=(('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))
PROFILES={'absent':None,'zero':(0.,0.),'drive':(.1,0.),'stiffness':(0.,.1),'combined':(.1,.1)}
MAX_PROGRAMS,MAX_FD=121,576

def require(ok,msg):
    if not ok: raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,value):
    p=Path(path);tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def maxabs(v):return float(np.max(np.abs(v))) if np.asarray(v).size else 0.
def gate(checks,name,error,bound):
    e,b=float(error),float(bound);passed=math.isfinite(e) and math.isfinite(b) and 0<=e<=b
    checks[name]={'error':e,'bound':b,'passed':bool(passed)}
    return passed
def passed(checks):return all(c['passed'] for c in checks.values())

def helpers(manifest):
    s=manifest['sources'];wm=read(s['water_manifest']);old=wm['sources']
    h={k:module(s[v],'donor_'+k) for k,v in [('native','native_parser'),('gradient','gradient_parser'),('matrix','matrix_parser')]}
    h.update({k:module(old[k],'donor_'+k) for k in ('moments','induction')})
    h['water']=module(s['water_qualification_source'],'donor_old_water')
    h['dual']=module(pin(HERE/'MATRIX.py'),'donor_dual_export')
    return h

def metadata(row,job):
    r=copy.deepcopy(row);r['input']=job['input'];atoms={};alpha={};charges={};gamma={}
    for src in ('parameters','key'):
        for line in verify(r['input'][src]).read_text().splitlines():
            raw=line.split()
            if not raw or raw[0].lower() not in ('atom','polarize','multipole','polpair'):continue
            t=shlex.split(line,comments=True)
            if not t:continue
            k=t[0].lower()
            if k=='atom':atoms[int(t[1])]={'name':t[3][:3],'mass':float(t[6])}
            elif k=='polarize':alpha[int(t[1])]=float(t[2])
            elif k=='multipole':charges[int(t[1])]=float(t[-1])
            elif k=='polpair':gamma[tuple(sorted(map(int,t[1:3])))]=float(t[3])
    xyz=[s.split() for s in verify(r['input']['xyz']).read_text().splitlines()[1:]]
    require(len(xyz)==r['atoms'] and [int(t[5]) for t in xyz]==r['types'],'Runtime atom inventory')
    r['positions_A']=[[float(x) for x in t[2:5]] for t in xyz]
    r['native_names']=[atoms[t]['name'] for t in r['types']]
    masses=np.array([atoms[t]['mass'] for t in r['types']]);r['COM_A']=(np.array(r['positions_A']).T@masses/masses.sum()).tolist()
    r['alpha_A3']=[alpha[t] for t in r['types']];r['charge']=math.fsum(charges[t] for t in r['types'])
    r['energy_replay_kcal_mol']=None
    if r['kind'].startswith('donor'):require(abs(r['charge']-2.)<1e-12,'Original donor charge+2')
    if 403 in r['types'] and any(t in (500,501) for t in r['types']):
        mt=next(t for t in r['types'] if t in (500,501));r['fixed_gamma']=gamma[tuple(sorted((mt,403)))]
    return r

def regularized(row,profile):
    coef=PROFILES[profile];a,c=coef if coef else (0.,0.)
    xyz=np.array(row['positions_A']);n=row['atoms'];h=np.zeros((n,3));beta=np.zeros(n);pairs=[]
    if not a and not c:return {'h':h,'beta':beta,'pairs':pairs}
    m=next(i for i,t in enumerate(row['types']) if t in (500,501));anchor=(1+5+25/3)*math.exp(-5.)
    for i,t in enumerate(row['types']):
        if t!=403:continue
        v=xyz[i]-xyz[m];r=float(np.linalg.norm(v));z=2*r;f=(((1+z+z*z/3)*math.exp(-z))/anchor)**2
        kappa=row['fixed_gamma']/math.sqrt(row['alpha_A3'][m]*row['alpha_A3'][i])
        phi=kappa if r==0 else -math.expm1(-kappa*r**3)/r**3;w=2*f/(1+f)
        beta[i]=math.expm1(c*w)/row['alpha_A3'][i];h[i]=3*a*f*phi*v
        pairs.append([m+1,i+1,row['fixed_gamma'],kappa,r,f,w,beta[i]])
    return {'h':h,'beta':beta,'pairs':pairs}

def numeric_records(text):
    out={}
    for line in text.splitlines():
        t=line.split()
        if not t or not t[0].startswith(('DONOR_RESPONSE_','LOCAL_RESPONSE_')):continue
        try:v=[float(x.replace('D','E')) for x in t[1:]]
        except ValueError as exc:raise ValueError('Malformed native response line') from exc
        require(v and all(math.isfinite(x) for x in v),'Finite native response records')
        out.setdefault(t[0],[]).append(v)
    return out

def readback(output,job,row,h):
    text=Path(output).read_text();allrec=numeric_records(text);n=row['atoms']
    active=row['kind'].startswith('donor') and job['profile'] in ('drive','stiffness','combined')
    if not active:
        require(not any(k in allrec for k in ('DONOR_RESPONSE_PARAMETERS','DONOR_RESPONSE_PAIR','DONOR_RESPONSE_FIELD','DONOR_RESPONSE_TOTAL_ENERGY')),'Donor bypass changed response path')
        if row.get('water_profile_active'):
            h['water'].response_readback(output,{'profile':'combined'},row)
        else:require(not any(k.startswith('LOCAL_RESPONSE_') for k in allrec),'Unexpected active correction on bypass')
    else:
        params=allrec.get('DONOR_RESPONSE_PARAMETERS',[]);require(params and len(params)%2==0,'Donor coefficients per solve')
        for i,v in enumerate(params):require(v==[500+i%2,403,*PROFILES[job['profile']]],'Actual donor coefficients differ')
        spectra=allrec.get('LOCAL_RESPONSE_SPECTRUM',[]);require(len(spectra)==len(params)//2,'Spectrum per donor solve')
        for v in spectra:
            require(len(v)==5,'Spectrum grammar');lo,hi,sym,eig,orth=v
            require(hi>0 and lo>1e-6*hi and 0<=sym<=1e-10 and 0<=eig<=1e-10*max(1,hi) and 0<=orth<=1e-10,'Donor spectral/arithmetic admission')
        dual=allrec.get('LOCAL_RESPONSE_INDUCED',[]);require(len(dual)==n*len(params)//2,'Full dual site denominator')
        for start in range(0,len(dual),n):require([v[0] for v in dual[start:start+n]]==list(range(1,n+1)) and all(len(v)==7 for v in dual[start:start+n]),'Full dual site ordering')
        require(len(allrec.get('LOCAL_RESPONSE_DUAL',[]))==len(params)//2,'Dual channel diagnostic denominator')
        for d in allrec['LOCAL_RESPONSE_DUAL']:require(len(d)==2 and min(d)>=0,'Dual channel diagnostics')
    # FD energy solves follow the central analytic readback. Never use their
    # final u/v as central dipoles merely because testgrad restores coordinates.
    central=numeric_records(text.split('DONOR_RESPONSE_TOTAL_ENERGY',1)[0]) if active and job['program']=='testgrad' else allrec
    result={'active':active,'all_solve_count':len(allrec.get('DONOR_RESPONSE_PARAMETERS',[]))//2,
            'permanent_lab_eA':central.get('DONOR_RESPONSE_PERMANENT_LAB_EA',[]),
            'charge':central.get('DONOR_RESPONSE_CHARGE',[])}
    if active:
        uv=np.array(central['LOCAL_RESPONSE_INDUCED'][-n:]);result['u_eA']=uv[:,1:4].tolist();result['v_eA']=uv[:,4:7].tolist()
        if job['program']=='analyze' and job['argv_tail'][1]=='EM':
            result['energy_solve_u_eA']=[v[1:4] for v in central['LOCAL_RESPONSE_INDUCED'][:n]]
            result['energy_solve_v_eA']=[v[4:7] for v in central['LOCAL_RESPONSE_INDUCED'][:n]]
        result['spectrum']=central['LOCAL_RESPONSE_SPECTRUM'][-1]
        ex=regularized(row,job['profile']);ar=h['matrix'].arithmetic(3*n,ex['h'],ex['beta'],ex['pairs'])
        require(maxabs(np.array(central['DONOR_RESPONSE_PAIR'][-2:])-ex['pairs'])<=ar,'Actual donor pair coefficients/geometry')
        require(maxabs(np.array(central['DONOR_RESPONSE_FIELD'][-2:])[:,1:]-ex['h'][[3,4]])<=ar,'Actual donor driving field')
        require([p[0] for p in central['DONOR_RESPONSE_FIELD'][-2:]]==[4.,5.],'Both equivalent oxygens')
        if job['program']=='analyze':require(result['all_solve_count']==(2 if job['argv_tail'][1]=='EM' else 1),'Analyze solve denominator')
        if job['program']=='testgrad':
            require(result['all_solve_count']==1+job['internal_FD_energies'],'Observed Cartesian solve denominator')
            require(len(allrec.get('DONOR_RESPONSE_TOTAL_ENERGY',[]))==1,'One central full gradient report')
            gg=allrec.get('DONOR_RESPONSE_TOTAL_GRADIENT',[]);vv=allrec.get('DONOR_RESPONSE_TOTAL_VIRIAL',[])
            require(len(gg)==n and [x[0] for x in gg]==list(range(1,n+1)) and all(len(x)==4 for x in gg),'Full central gradient rows')
            require(len(vv)==3 and [x[0] for x in vv]==[1.,2.,3.] and all(len(x)==4 for x in vv),'Full virial rows')
            result.update(full_energy=allrec['DONOR_RESPONSE_TOTAL_ENERGY'][0][0],gradient=[x[1:] for x in gg],virial=[x[1:] for x in vv])
    if result['permanent_lab_eA']:
        require(all(len(v)==3 for v in result['permanent_lab_eA']) and all(v==[row['charge']] or abs(v[0]-row['charge'])<2e-12 for v in result['charge']),'Permanent moment/charge readback')
    # Preserve bounds from the literal printed tokens, not a guessed number
    # of decimal places. For FD files, use only the central pre-gradient block.
    ct=text.split('DONOR_RESPONSE_TOTAL_ENERGY',1)[0] if active and job['program']=='testgrad' else text
    def print_rows(label,tail):
        rr=[line.split()[1:] for line in ct.splitlines() if line.strip().startswith(label+' ')]
        return [[h['matrix'].half_print(t) for t in r] for r in rr[-tail:]]
    if active:
        pb=np.array(print_rows('LOCAL_RESPONSE_INDUCED',n))[:,1:]
        result['u_print_bounds_eA']=pb[:,:3].tolist();result['v_print_bounds_eA']=pb[:,3:].tolist()
    if result['permanent_lab_eA']:result['permanent_print_bounds_eA']=print_rows('DONOR_RESPONSE_PERMANENT_LAB_EA',1)[0]
    if active and job['program']=='testgrad':
        lines=[s.split()[2:] for s in text.splitlines() if s.strip().startswith('DONOR_RESPONSE_TOTAL_GRADIENT ')]
        result['gradient_print_bounds']=[[h['matrix'].half_print(x) for x in r] for r in lines]
    return result

def parse_cell(job,m,output,h):
    row=metadata(m['rows'][job['cell']],job);text=Path(output).read_text();checks={}
    rb=readback(output,job,row,h)
    # The one old-binary baseline intentionally retains its original quiet key.
    induction=None if job['binary_role']=='original_testgrad' else h['induction'].induction(output)
    if rb['active'] and job['program']=='testgrad':require(induction['solutions']==1+job['internal_FD_energies'],'Observed testgrad induction denominator')
    rec={'output':pin(output),'response':rb,'induction':induction,'checks':checks}
    if job['program']=='testgrad':
        value=h['gradient'].gradient(output,row['atoms'],job['fd_total_separation_A'] is not None,h['native'])
        value.update(fd_total_separation_A=job['fd_total_separation_A'],fd_each_side_displacement_A=job['fd_each_side_A'])
        if job['internal_FD_energies']:gate(checks,'Cartesian_gradient',value['max_component_error_kcal_mol_A'],.001)
        if rb['active']:
            gate(checks,'full_printed_energy',abs(rb['full_energy']-value['energy']),5e-9+1e-10)
            gate(checks,'full_printed_gradient',maxabs(np.array(rb['gradient'])-value['gradient']),5e-9+1e-10)
            pos=np.array(row['positions_A']);gg=np.array(rb['gradient']);vir=np.array(rb['virial'])
            expected=.5*(pos.T@gg+gg.T@pos)
            gate(checks,'virial_coordinate_gradient',maxabs(vir-expected),h['matrix'].arithmetic(3*row['atoms'],vir,expected))
            if job['field_au'] is None:
                gate(checks,'zero_field_net_gradient',maxabs(gg.sum(axis=0)),.001)
                gate(checks,'zero_field_net_torque',maxabs(np.cross(pos,gg).sum(axis=0)),.001*max(1.,float(np.linalg.norm(pos,axis=1).sum())))
        rec['value']=value
    else:
        value=h['moments'].parse_energy(output);rec['value']={'energy':value['energy_kcal_mol'],'components':value['components']}
        if job['argv_tail'][1]=='EM':
            ip=h['moments'].induced_file(output)
            moments=h['moments'].parse_moments(output,ip,row,h['induction'].induction);rec['moments']=moments
            require(len(rb['permanent_lab_eA'])==2 and len(rb['charge'])==2,'Two permanent laboratory readbacks')
            gate(checks,'two_permanent_moments',maxabs(np.array(rb['permanent_lab_eA'][0])-rb['permanent_lab_eA'][1]),h['matrix'].arithmetic(3*row['atoms'],rb['permanent_lab_eA']))
            if job['exports']:
                paths=sorted(Path(output).parent.glob('response_matrix_*.dat'));require(len(paths)==2,'Exactly two native exports')
                rec['exports']=[pin(p) for p in paths]
                rec['operator']=h['dual'].qualify_exports(h['matrix'],rec['exports'],row,m['units'],moments)
                z=rec['operator']['exports'][1];uv=np.array(z['vectors'])[:,2:].reshape(row['atoms'],3,2)
                # Export vectors are Cartesian-major; preserve atom/axis order.
                u=uv[:,:,0];v=uv[:,:,1]
                if rb['active']:
                    gate(checks,'export_full_u',maxabs(u-rb['u_eA']),z['numerical_bounds']['arithmetic'])
                    gate(checks,'export_full_v',maxabs(v-rb['v_eA']),z['numerical_bounds']['arithmetic'])
                rb.update(u_eA=u.tolist(),v_eA=v.tolist())
            require('u_eA' in rb and 'v_eA' in rb,'Conjugate duals available')
            u,v=np.array(rb['u_eA']),np.array(rb['v_eA']);perm=np.array(rb['permanent_lab_eA'][-1]);B=m['units']['bohr_A'];D=m['units']['debye_per_eA']
            mu=(perm+.5*(u+v).sum(axis=0))/B;rec['mu_lab_au']=mu.tolist()
            gate(checks,'saved_u_site_precision',maxabs(D*u-moments['induced_Debye']),.5e-6+2e-9)
            native_lab=D*(perm+(.5*(u+v) if rb['active'] else u).sum(axis=0))
            gate(checks,'native_printed_moment_origin',maxabs(native_lab-moments['total_lab_Debye']),.0005+2e-9)
    rec['numerical_checks_passed']=passed(checks)
    return rec

def find(m,records,cell,profile,stage,program,suffix=None,field=None,fd=None):
    jobs=[j for j in m['programs'] if j['cell']==cell and j['profile']==profile and j['stage']==stage and j['program']==program
          and (suffix is None or suffix in j['id']) and (field is None or j['field_au']==field)
          and (fd is None or j['fd_total_separation_A']==fd)]
    require(len(jobs)==1,'Unique planned lookup '+str((cell,profile,stage,program,suffix,field,fd)))
    return records[jobs[0]['id']]

def baseline_checks(m,records):
    checks={};exact={}
    for j in m['programs']:
        if j['stage'] not in ('A','B') or j['binary_role']=='original_testgrad':continue
        r=records[j['id']];b=m['rows'][j['cell']]['baseline']
        if b is None:b=records['old__acetate7__gradient']['value']
        gate(checks,j['id']+'_energy',abs(r['value']['energy']-b['energy']),1e-4)
        exact[j['id']]={'energy':r['value']['energy']==b['energy']}
        if j['program']=='testgrad':
            gate(checks,j['id']+'_forces',maxabs(np.array(r['value']['forces'])-b['forces']),1e-3)
            exact[j['id']]['forces']=r['value']['forces']==b['forces']
    return {'checks':checks,'exact_printed':exact,'passed':passed(checks)}

def operator_checks(m,records,h):
    checks={};details={};units=m['units'];B=units['bohr_A'];K=units['electric'];Ha=units['hartree_kcal_mol']
    for metal in ('La','Dy'):
        cid=metal+'_central';base=find(m,records,cid,'absent','A','analyze');z0=base['operator']['exports'][1]
        A0=np.array(z0['A']);v0=np.array(z0['vectors']);bt=np.tile(m['field_direction'],m['rows'][cid]['atoms'])
        chi0=Ha/(K*B*B)*float(bt@np.linalg.solve(A0,bt))
        for profile in ('drive','stiffness','combined'):
            r=find(m,records,cid,profile,'C','analyze');z=r['operator']['exports'][1];A=np.array(z['A']);v=np.array(z['vectors'])
            job=next(j for j in m['programs'] if j['id']==r['job']['id']);row=metadata(m['rows'][cid],job);ex=regularized(row,profile)
            H=np.diag(np.repeat(ex['beta'],3));hh=ex['h'].reshape(-1);ar=h['matrix'].arithmetic(len(A),A,A0,v,v0,hh)
            tag=metal+'_'+profile
            gate(checks,tag+'_G',maxabs(A-A0-H),ar)
            for k,name in enumerate(('Ed','Ep')):gate(checks,tag+'_'+name,maxabs(v[:,k]-v0[:,k]-hh),ar)
            chi=Ha/(K*B*B)*float(bt@np.linalg.solve(A,bt));ctol=h['matrix'].arithmetic(len(A),chi,chi0)
            gate(checks,tag+'_chi_structure',abs(chi-chi0) if profile=='drive' else max(0.,chi-chi0),ctol)
            g1=find(m,records,cid,profile,'C','testgrad',fd=.0001);g2=find(m,records,cid,profile,'C','testgrad',fd=.00005)
            gate(checks,tag+'_energy_route',abs(g1['value']['energy']-r['value']['energy']),1e-4)
            gate(checks,tag+'_Cartesian_two_steps',maxabs(np.array(g1['value']['numerical_gradient'])-g2['value']['numerical_gradient']),.001+1e-8)
            details[tag]={'matrix_chi_a0_cubed':chi,'baseline_chi_a0_cubed':chi0,'Ed_Ep_max_difference':maxabs(v[:,0]-v[:,1])}
    return {'checks':checks,'states':details,'passed':passed(checks)}

def field_checks(m,records,h):
    checks={};states={};units=m['units'];Ha=units['hartree_kcal_mol'];K=units['electric'];B=units['bohr_A'];D=units['debye_per_eA'];t=np.array(m['field_direction'])
    for metal in ('La','Dy'):
        cid=metal+'_central';n=m['rows'][cid]['atoms'];bt=np.tile(t,n)
        for profile in ('drive','stiffness','combined'):
            tag=metal+'_'+profile;zero=find(m,records,cid,profile,'C','analyze');z=zero['operator']['exports'][1]
            A=np.array(z['A']);ed,ep,u0,v0=np.array(z['vectors']).T;ar=z['numerical_bounds']['arithmetic'];an=math.sqrt(3*n)*ar
            ud=np.linalg.solve(A,ed);vd=np.linalg.solve(A,ep);column=np.linalg.solve(A,bt);chi=Ha/(K*B*B)*float(bt@column)
            perm=np.array(zero['response']['permanent_lab_eA'][-1]);mu0=float(t@perm/B+.5*bt@(ud+vd)/B)
            energy_zero=zero['operator']['exports'][0]
            e0=zero['value']['energy'];e0bound=5e-9+ar+.5*K*np.linalg.norm(ep)*energy_zero['numerical_bounds']['native_vs_dense_u']
            estimates={};bounds={}
            for f in (-.001,-.0005,.0005,.001):
                r=find(m,records,cid,profile,'D','analyze',field=f);rb=r['response'];u=np.array(rb['u_eA']).reshape(-1);v=np.array(rb['v_eA']).reshape(-1)
                df=ed+Ha/(K*B)*f*bt;pf=ep+Ha/(K*B)*f*bt
                q=r['induction']['maximum_printed_RMS_residual_Debye'];amp=z['numerical_bounds']['post_peek_amplification'];resbound=amp*math.sqrt(n)*(q+.5e-10)/D+an
                for name,x,rhs in [('u',u,df),('v',v,pf)]:gate(checks,f'{tag}_{f}_{name}_RHS',np.linalg.norm(A@x-rhs),resbound)
                ub=(np.linalg.norm(A@u-df)+2*an)/z['spectrum']['A'][0];vb=(np.linalg.norm(A@v-pf)+2*an)/z['spectrum']['A'][0]
                # The printed energy and moment use distinct EM solves.
                # The largest printed residual covers both; do not bound the
                # first E solve with only the last M solve's actual residual.
                energy_ub=(resbound+2*an)/z['spectrum']['A'][0]
                eb=5e-9+ar+.5*K*np.linalg.norm(pf)*energy_ub;mb=.5*math.sqrt(n)*(ub+vb)/B+ar
                bounds[f]={'energy':float(eb),'moment':float(mb)}
                energy_u=np.array(rb['energy_solve_u_eA']).reshape(-1)
                zero_u=np.array(energy_zero['vectors'])[:,2]
                contract=-Ha*f*float(t@perm/B)-.5*K*(float(energy_u@pf)-float(zero_u@ep))
                gate(checks,f'{tag}_{f}_energy_native_contraction',abs((r['value']['energy']-e0)-contract),1e-8+h['matrix'].arithmetic(3*n,contract,e0))
                gate(checks,f'{tag}_{f}_energy_quadratic',abs((r['value']['energy']-e0)+Ha*(f*mu0+.5*f*f*chi)),e0bound+eb)
                gate(checks,f'{tag}_{f}_moment_linear',abs(float(t@np.array(r['mu_lab_au']))-mu0-f*chi),mb)
            for label,step in [('h',.001),('half',.0005)]:
                plus=find(m,records,cid,profile,'D','analyze',field=step);minus=find(m,records,cid,profile,'D','analyze',field=-step)
                ep1,em1=plus['value']['energy'],minus['value']['energy']
                ce=-(ep1+em1-2*e0)/(Ha*step*step);cm=float(t@(np.array(plus['mu_lab_au'])-minus['mu_lab_au']))/(2*step);odd=-(ep1-em1)/(2*step*Ha)
                estimates[label]={'chi_energy':ce,'chi_moment':cm,'odd_energy_moment':odd,'zero_moment':float(t@np.array(zero['mu_lab_au'])),
                    'roundoff_bounds':{'chi_energy':2e-8/(Ha*step*step),'odd_energy_moment':5e-9/(Ha*step)},
                    'conditioned_energy_curvature_bound':(bounds[step]['energy']+bounds[-step]['energy']+2*e0bound)/(Ha*step*step),
                    'conditioned_moment_slope_bound':(bounds[step]['moment']+bounds[-step]['moment'])/(2*step)}
                gate(checks,tag+'_'+label+'_E_mu',abs(ce-cm),.5)
                gate(checks,tag+'_'+label+'_odd_zero',abs(odd-estimates[label]['zero_moment']),1e-4)
                gate(checks,tag+'_'+label+'_E_matrix',abs(ce-chi),estimates[label]['conditioned_energy_curvature_bound'])
                gate(checks,tag+'_'+label+'_mu_matrix',abs(cm-chi),estimates[label]['conditioned_moment_slope_bound'])
            for estimator in ('chi_energy','chi_moment'):gate(checks,tag+'_'+estimator+'_steps',abs(estimates['h'][estimator]-estimates['half'][estimator]),.5)
            states[tag]={'estimates':estimates,'matrix_chi':chi,'dense_zero_moment':mu0,'field_bounds':bounds}
    paired={}
    for profile in ('drive','stiffness','combined'):
        a,b=states['La_'+profile]['estimates'],states['Dy_'+profile]['estimates'];p={k:{v:b[k][v]-a[k][v] for v in ('chi_energy','chi_moment','odd_energy_moment','zero_moment')} for k in ('h','half')}
        for k in ('h','half'):gate(checks,profile+'_paired_'+k,abs(p[k]['chi_energy']-p[k]['chi_moment']),1.)
        for v in ('chi_energy','chi_moment'):gate(checks,profile+'_paired_'+v+'_steps',abs(p['h'][v]-p['half'][v]),1.)
        paired[profile]=p
    return {'checks':checks,'states':states,'paired':paired,'passed':passed(checks)}

def moment_bound(rec,row,units,h):
    """Conservative conditional moment error from actual normalized spectrum.

    This is not a force sensitivity or finite-coordinate-step error bound.
    """
    rb=rec['response'];n=row['atoms'];alpha,_=h['matrix'].settings(row)
    aa=[alpha[t] for t in row['types']];lo,hi=rb['spectrum'][:2]
    peek=h['matrix'].settings(row)[1]
    lam_min=lo/max(aa);amp=1+abs(peek)*hi*max(aa)/min(aa)
    ar=h['matrix'].arithmetic(3*n,rb['u_eA'],rb['v_eA'],hi/min(aa))
    q=rec['induction']['maximum_printed_RMS_residual_Debye']
    ub=(amp*math.sqrt(n)*(q+.5e-10)/units['debye_per_eA']+3*math.sqrt(3*n)*ar)/lam_min
    pb=np.array(rb['permanent_print_bounds_eA'])+.5*(np.array(rb['u_print_bounds_eA'])+rb['v_print_bounds_eA']).sum(axis=0)
    return {'solution_norm_bound_eA':float(ub),'component_solver_bound_au':math.sqrt(n)*ub/units['bohr_A'],
            'component_print_bounds_au':(pb/units['bohr_A']).tolist(),'normalized_spectrum':rb['spectrum'],
            'force_sensitivity_bound_available':False}

def mechanics_checks(m,records,h):
    checks={};rotation={};mixed={};t=np.array(m['field_direction']);Ha=m['units']['hartree_kcal_mol']
    for metal in ('La','Dy'):
        cid=metal+'_central';row=m['rows'][cid];base=find(m,records,cid,'combined','C','testgrad',fd=.0001)
        gg=np.array(base['response']['gradient']);vir=np.array(base['response']['virial']);xyz=np.array(row['positions_A']);fr=row['rotation_frame']
        tangent=np.zeros_like(xyz);inds=fr['ligand_indices'];tangent[inds]=np.cross(np.array(fr['rotation_unit']),xyz[inds]-fr['centroid_A'])/fr['rotation_length_A']
        analytic=float(np.sum(gg*tangent));rotation[metal]={'analytic':analytic,'steps':{}};mixed[metal]={}
        for q in (.01,.005):
            ends={}
            for sign in (-1,1):
                ec=next(k for k,r in m['rows'].items() if r['kind']=='donor_rotation' and r['metal']==metal and r['q_A']==sign*q)
                a=find(m,records,ec,'combined','E','analyze');g=find(m,records,ec,'combined','E','testgrad');ends[sign]=a
                job=a['job'];a['mixed_moment_bound']=moment_bound(a,metadata(m['rows'][ec],job),m['units'],h)
                gate(checks,ec+'_energy_route',abs(a['value']['energy']-g['value']['energy']),1e-4)
            slope=(ends[1]['value']['energy']-ends[-1]['value']['energy'])/(2*q)
            gate(checks,metal+f'_rotation_{q}',abs(slope-analytic),.1+5e-9/q)
            dmu=float(t@(np.array(ends[1]['mu_lab_au'])-ends[-1]['mu_lab_au']))/(2*q)
            rotation[metal]['steps'][str(q)]={'energy_slope':slope,'moment_slope':dmu}
            for step in (.001,.0005):
                gp=find(m,records,cid,'combined','D','testgrad',field=step);gm=find(m,records,cid,'combined','D','testgrad',field=-step)
                # F_q is negative energy gradient, so its field derivative
                # equals the positive geometry derivative of conjugate mu.
                df=-float(np.sum((np.array(gp['response']['gradient'])-gm['response']['gradient'])*tangent))/(2*step*Ha)
                force_print=float(np.sum((np.array(gp['response']['gradient_print_bounds'])+gm['response']['gradient_print_bounds'])*abs(tangent)))/(2*step*Ha)
                endpoint_allowance=0.
                for sign in (-1,1):
                    bnd=ends[sign]['mixed_moment_bound']
                    endpoint_allowance+=(float(abs(t)@np.array(bnd['component_print_bounds_au']))+bnd['component_solver_bound_au']*float(np.linalg.norm(t)))/(2*q)
                bound=1e-4/q+force_print+endpoint_allowance
                gate(checks,metal+f'_Maxwell_{q}_{step}',abs(df-dmu),bound)
                mixed[metal][f'{q}_{step}']={'force_field_derivative':df,'moment_geometry_derivative':dmu,'bound':bound,
                    'force_printing_bound':force_print,'endpoint_moment_bound':endpoint_allowance,
                    'force_solver_sensitivity_bound_available':False}
        rotation[metal]['q_step_difference_descriptive']=rotation[metal]['steps']['0.01']['moment_slope']-rotation[metal]['steps']['0.005']['moment_slope']
        for kind in ('translation','rotation'):
            r=find(m,records,cid,'combined','F','testgrad',suffix=kind);j=r['job'];R=np.eye(3) if kind=='translation' else np.array(j['transform']['matrix'])
            gate(checks,metal+'_'+kind+'_energy',abs(r['value']['energy']-base['value']['energy']),1e-4)
            gate(checks,metal+'_'+kind+'_gradient',maxabs(np.array(r['response']['gradient'])-gg@R.T),.001)
            gate(checks,metal+'_'+kind+'_virial',maxabs(np.array(r['response']['virial'])-R@vir@R.T),.001*max(1.,float(np.linalg.norm(xyz,axis=1).sum())))
            for channel in ('u_eA','v_eA'):gate(checks,metal+'_'+kind+'_'+channel,maxabs(np.array(r['response'][channel])-np.array(base['response'][channel])@R.T),1e-8)
        strain={}
        for label,step in [('h',.0001),('half',.00005)]:
            plus=find(m,records,cid,'combined','G','analyze',suffix='strain_'+label+'_plus');minus=find(m,records,cid,'combined','G','analyze',suffix='strain_'+label+'_minus');tr=plus['job']['transform']
            pred=float(np.sum(np.array(tr['S'])*vir))/tr['normalization_L_A'];obs=(plus['value']['energy']-minus['value']['energy'])/(2*step)
            gate(checks,metal+'_strain_'+label,abs(obs-pred),.001+5e-9/step);strain[label]=obs
        gate(checks,metal+'_strain_two_steps',abs(strain['h']-strain['half']),.001+5e-9/.0001+5e-9/.00005)
    paired={}
    for key in mixed['La']:
        a,b=mixed['La'][key],mixed['Dy'][key];df=b['force_field_derivative']-a['force_field_derivative'];dm=b['moment_geometry_derivative']-a['moment_geometry_derivative'];bound=a['bound']+b['bound']
        gate(checks,'paired_Maxwell_'+key,abs(df-dm),bound);paired[key]={'force_field_derivative':df,'moment_geometry_derivative':dm,'bound':bound}
    return {'checks':checks,'rotation':rotation,'mixed':mixed,'paired_mixed':paired,'passed':passed(checks),
            'limitation':'The inherited mixed-derivative allowance is numerical, not a bound on finite-q truncation, SCF error or physical accuracy.'}

def context(cfg,runtime=False):
    require(cfg['implementation']==pin(__file__) and cfg['matrix_adapter']==pin(HERE/'MATRIX.py'),'Runner source binding')
    m=read(cfg['manifest']);require(m['plan']==pin(HERE/'QUALIFICATION_PLAN.md') and m['implementation']==pin(HERE/'PREPARE.py'),'Exact prepared plan/source')
    require(len(m['programs'])==121 and sum(j['internal_FD_energies'] for j in m['programs'])==576,'Finite prepared denominator')
    for p in m['sources'].values():verify(p)
    for j in m['programs']:
        for p in j['input'].values():verify(p)
        source=verify(j['source_input']['key']).read_text();overlay=j['key_overlay'];expected=source+('\n'+'\n'.join(overlay)+'\n' if overlay else '')
        require(verify(j['input']['key']).read_text()==expected,'Immutable key overlay')
    for k in ('source_manifest','source_review','build_receipt','analyze','testgrad','runner_review'):verify(cfg[k])
    build=read(cfg['build_receipt']);require(build['status']=='complete' and build['manifest']==cfg['source_manifest'] and all(build['executables'][k]==cfg[k] for k in ('analyze','testgrad')),'Exact qualified build binding')
    audit=read(cfg['source_review']);require(audit['source_manifest']==cfg['source_manifest'] and audit['status'].startswith('passed_'),'Native source review')
    review=read(cfg['runner_review']);require(review['implementation']==cfg['implementation'] and review['manifest']==cfg['manifest'] and review['matrix_adapter']==cfg['matrix_adapter'] and review['source_manifest']==cfg['source_manifest'] and review['status']=='passed_donor_extension_qualification_source_review','Independent finite runner review')
    require(1<=cfg['workers']<=12 and cfg['maximum_native_programs']==121 and cfg['maximum_internal_FD_energies']==576,'Finite runtime limits')
    release=read(cfg['root_release']);require(release['status']=='approved_donor_extension_numerical_qualification' and release['implementation']==cfg['implementation'] and release['manifest']==cfg['manifest'] and release['source_manifest']==cfg['source_manifest'] and release['build_receipt']==cfg['build_receipt'],'Root exact release binding')
    require(release['maximum_native_programs']==121 and release['maximum_internal_FD_energies']==576 and release['workers']==cfg['workers'] and release['fit'] is False,'Numerical-only release')
    if runtime:require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=cfg['workers'] and all(os.environ.get(k)=='1' for k in THREADS),'Root Slurm allocation and single-thread numerical libraries')
    return m,helpers(m)

def run(cfg,out):
    m,h=context(cfg,True);out=Path(out);out.mkdir(parents=True,exist_ok=False)
    tick,cpu=time.monotonic(),time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN);lock=threading.RLock();stop=threading.Event()
    result={'status':'running','implementation':pin(__file__),'config':cfg,'manifest':cfg['manifest'],'job_id':os.environ['SLURM_JOB_ID'],
            'rows':{j['id']:{'status':'not_attempted'} for j in m['programs']},'requested_slots':0,'native_program_invocations':0,
            'internal_FD_energies_requested':0,'internal_FD_energies_completed':0,'global_stop':False,'error':None,
            'numerical_checks_passed':False,'physical_accuracy_qualified':False}
    def persist():save(out/'RESULT.json',copy.deepcopy(result))
    def worker(job):
        rec={'status':'started','job':copy.deepcopy(job),'receipts':[]};cell=out/'cells'/job['id'];claimed=False
        try:
            with lock:
                if stop.is_set():return
                require(result['requested_slots']<121,'Finite reserved-slot ceiling')
                claimed=True
                result['requested_slots']+=1;result['rows'][job['id']]=copy.deepcopy(rec);persist()
            cell.mkdir(parents=True)
            for k,name in FILES:shutil.copyfile(verify(job['input'][k]),cell/name)
            rec['runtime_input']={k:pin(cell/name) for k,name in FILES}
            executable=m['sources']['original_testgrad'] if job['binary_role']=='original_testgrad' else cfg[job['program']]
            argv=[str(verify(executable)),*job['argv_tail']]
            with lock:
                if stop.is_set():rec['status']='not_attempted';return
                require(result['native_program_invocations']<121 and result['internal_FD_energies_requested']+job['internal_FD_energies']<=576,'Finite actual invocation ceiling')
                save(cell/'START.json',{'job':job,'argv':argv,'input':rec['runtime_input'],'job_id':os.environ['SLURM_JOB_ID'],'time_ns':time.time_ns()});rec['start']=pin(cell/'START.json')
                result['native_program_invocations']+=1;result['internal_FD_energies_requested']+=job['internal_FD_energies'];persist()
            output=cell/(job['program']+'.out');h['native'].invoke(argv,cell,output,rec['receipts'])
            rec.update(parse_cell(job,m,output,h));rec['status']='complete' if rec['numerical_checks_passed'] else 'failed_numerical'
            with lock:result['internal_FD_energies_completed']+=job['internal_FD_energies']
            require(rec['numerical_checks_passed'],'Per-program numerical check failed')
        except BaseException as exc:
            rec['status']='failed';rec['error']=repr(exc);stop.set()
            with lock:
                if result['error'] is None:result['error']=repr(exc)
        finally:
            if claimed:
                with lock:result['rows'][job['id']]=copy.deepcopy(rec);result['global_stop']=stop.is_set();persist()
    try:
        for stages in (('A','B'),('C',),('D','E','F','G')):
            if stop.is_set():break
            jobs=[j for j in m['programs'] if j['stage'] in stages]
            with ThreadPoolExecutor(max_workers=cfg['workers']) as pool:list(pool.map(worker,jobs))
            if stop.is_set():break
            if stages==('A','B'):result['baseline']=baseline_checks(m,result['rows']);require(result['baseline']['passed'],'Baseline/bypass check failed')
            if stages==('C',):result['operator']=operator_checks(m,result['rows'],h);require(result['operator']['passed'],'Operator/Cartesian check failed')
        if not stop.is_set():
            result['field']=field_checks(m,result['rows'],h);result['mechanics']=mechanics_checks(m,result['rows'],h)
            require(result['field']['passed'] and result['mechanics']['passed'],'Response/mechanics numerical checks failed')
            require(result['native_program_invocations']==121 and result['internal_FD_energies_completed']==576,'Full finite completion')
            result['numerical_checks_passed']=True;result['status']='donor_extension_numerical_qualification_complete'
    except BaseException as exc:result['error']=repr(exc);stop.set()
    finally:
        if stop.is_set():result['status']='donor_extension_numerical_qualification_failed';result['global_stop']=True
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']={'wall_seconds':time.monotonic()-tick,'parent_CPU_seconds':time.process_time()-cpu,'all_native_child_CPU_seconds':ch.ru_utime+ch.ru_stime-ch0.ru_utime-ch0.ru_stime,'maximum_workers':cfg['workers'],'scheduler_allocated_CPU_seconds':None,'per_receipt_child_CPU_may_overlap_do_not_sum':True};persist()
    return result

def admit_result(cfg,result):
    m,h=context(cfg,False);r=read(result) if isinstance(result,dict) and 'path' in result else result
    require(r['config']==cfg and r['manifest']==cfg['manifest'] and r['implementation']==pin(__file__),'Saved result lineage')
    require(set(r['rows'])=={j['id'] for j in m['programs']},'Exact saved row denominator')
    require(r['status'] in ('donor_extension_numerical_qualification_complete','donor_extension_numerical_qualification_failed'),'Terminal result required')
    complete={};calls=0;fd=0
    for job in m['programs']:
        rec=r['rows'][job['id']]
        if rec['status']=='not_attempted':continue
        calls+=len(rec['receipts']);require(len(rec['receipts'])==1,'Single call per slot')
        for k,p in rec['runtime_input'].items():verify(p);require(p['sha256']==job['input'][k]['sha256'],'Actual input bytes')
        verify(rec['start']);receipt=rec['receipts'][0];verify(receipt['output'])
        require(receipt['returncode']==0,'Saved native completion')
        executable=m['sources']['original_testgrad'] if job['binary_role']=='original_testgrad' else cfg[job['program']]
        require(receipt['argv']==[executable['path'],*job['argv_tail']],'Actual binary/argv')
        parsed=parse_cell(job,m,receipt['output']['path'],h)
        # Stored raw arrays/pins must be exact; derived arithmetic is replayed
        # under its original numerical bounds, not host-dependent bit equality.
        for key in ('output','value','moments','exports'):
            if key in parsed:require(parsed[key]==rec[key],'Raw saved parse differs: '+key)
        complete[job['id']]=dict(parsed,job=job);fd+=job['internal_FD_energies']
    require(calls==r['native_program_invocations'] and fd==r['internal_FD_energies_completed'],'Actual counters')
    require(calls==121 and len(complete)==121 and fd==576,'Full qualification admission requires121/576')
    reports={'baseline':baseline_checks(m,complete),'operator':operator_checks(m,complete,h),'field':field_checks(m,complete,h),'mechanics':mechanics_checks(m,complete,h)}
    ok=all(v['numerical_checks_passed'] for v in complete.values()) and all(v['passed'] for v in reports.values())
    require(ok==r['numerical_checks_passed'],'Saved numerical status differs')
    if ok:require(r['status']=='donor_extension_numerical_qualification_complete' and r['error'] is None and r['global_stop'] is False,'Successful terminal status/error/stop')
    else:require(r['status']=='donor_extension_numerical_qualification_failed' and r['global_stop'] is True,'Failed qualification preserved')
    return {'passed':ok,'reports':reports,'new_native_calls':0,'complete':121,'internal_FD_energies':576}

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['run','admit']);p.add_argument('--config',required=True);p.add_argument('--out');p.add_argument('--result');a=p.parse_args();cfg=json.loads(Path(a.config).read_text())
    if a.action=='run':r=run(cfg,a.out)
    else:r=admit_result(cfg,pin(a.result))
    print(json.dumps({'status':r.get('status'),'passed':r.get('passed',r.get('numerical_checks_passed'))}));return 0 if r.get('passed',r.get('numerical_checks_passed')) else 1
if __name__=='__main__':raise SystemExit(main())
