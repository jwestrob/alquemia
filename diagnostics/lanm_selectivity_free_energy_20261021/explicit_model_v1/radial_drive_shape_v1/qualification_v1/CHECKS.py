#!/usr/bin/env python3
"""Narrow radial RHS/zero bridge checks around inherited native parsers."""
import copy,math,re
from pathlib import Path
import numpy as np
from PREPARE import require,pin,verify

def gate(checks,name,error,bound):
    ok=math.isfinite(float(error)) and 0<=float(error)<=float(bound)
    checks[name]={'error':float(error),'bound':float(bound),'passed':ok};require(ok,'Numerical check: '+name)

def delta(row,e):
    h=np.zeros((row['atoms'],3));metals=[i for i,t in enumerate(row['types']) if t in (500,501)]
    if not metals:return h,[]
    m=metals[0];rows=[];anchor=(1+5+25/3)*math.exp(-5)
    for i,t in enumerate(row['types']):
        if t!=36:continue
        R=np.asarray(row['positions_A'][i])-row['positions_A'][m];r=float(np.linalg.norm(R));kap=row['physical_POLPAIR']/math.sqrt(row['alpha_A3'][m]*row['alpha_A3'][i]);t=2*r;f=((1+t+t*t/3)*math.exp(-t)/anchor)**2;phi=kap if r==0 else -math.expm1(-kap*r**3)/r**3;psi=(r*r-6.25)/6.25;h[i]=3*e*f*phi*psi*R;rows.append([m+1,i+1,e,psi,*h[i]])
    return h,rows

def readback(output,job,row,m,h,ev):
    text=Path(output).read_text();n=row['atoms'];solves=1+6*n if job['fd'] else 2 if job['program']=='analyze' and job['stage']!='directions' else 1
    groups={}
    for line in text.splitlines():
        t=line.split()
        if t and t[0].startswith(('LOCAL_RESPONSE_','RADIAL_DRIVE_')):
            vals=[float(x.replace('D','E')) for x in t[1:]];require(all(math.isfinite(v) for v in vals),'Finite native readback');groups.setdefault(t[0],[]).append(vals)
    require(len(groups.get('LOCAL_RESPONSE_PARAMETERS',[]))==solves,'Observed native solve denominator')
    require(len(groups.get('LOCAL_RESPONSE_INDUCED',[]))==solves*n and len(groups.get('LOCAL_RESPONSE_DUAL',[]))==solves,'Dual solve denominator')
    require(len(groups.get('LOCAL_RESPONSE_SPECTRUM',[]))==solves,'Spectrum denominator')
    for p in groups['LOCAL_RESPONSE_PARAMETERS']:require(len(p)==6 and p[:2]==[row['metal_type'],36] and p[2:]==row['physical'][:4],'Original coefficient readback')
    for lo,hi,sym,resid,orth in groups['LOCAL_RESPONSE_SPECTRUM']:require(hi>0 and lo>1e-6*hi and 0<=sym<=1e-10 and 0<=resid<=1e-10*max(1.,hi) and 0<=orth<=1e-10,'Unchanged spectrum gates')
    for i,p in enumerate(groups['LOCAL_RESPONSE_INDUCED']):require(len(p)==7 and p[0]==1+i%n,'All native induced sites')
    for p in groups['LOCAL_RESPONSE_DUAL']:require(len(p)==2 and min(p)>=0,'Dual discrepancy readback')
    e=0. if job['profile'] in ('absent','zero') else row['e'];active=e!=0
    require(len(groups.get('RADIAL_DRIVE_PARAMETERS',[]))==(solves if active else 0),'Actual radial parameter denominator')
    for p in groups.get('RADIAL_DRIVE_PARAMETERS',[]):require(p==[row['metal_type'],36,e],'Actual radial parameter value')
    dh,expected=delta(row,e);actual=groups.get('RADIAL_DRIVE_PAIR',[]);require(len(actual)==len(expected)*solves if active else not actual,'Actual radial pair denominator')
    checks={};ar=h['matrix'].arithmetic(3*n,dh)
    if active and expected:gate(checks,'central_radial_value_readback',np.max(abs(np.asarray(actual[:len(expected)])-expected)),ar)
    if job['cell'] not in ('ion_La','ion_Dy','aquo9_central_La__near__two_water'):
        positions=[x.start() for x in re.finditer(r'(?m)^\s*LOCAL_RESPONSE_PARAMETERS ',text)]
        scope=text[:positions[1]] if job['fd'] else text
        local=ev.local_readback(scope,row,row['physical'],m['units'],solves=1 if job['fd'] else solves)
    else:local=None
    return {'native_solve_blocks':solves,'central_induced_u_v_eA':[x[1:] for x in groups['LOCAL_RESPONSE_INDUCED'][:n]],'last_induced_u_v_eA':[x[1:] for x in groups['LOCAL_RESPONSE_INDUCED'][-n:]],'local':local,'radial_checks':checks,'groups':groups}

def zero_identity(record,row,job,h):
    checks={}
    if 'baseline_final_cell' in row:
        old=row['baseline_final_cell']
        if job['program']=='testgrad':
            new=record['value'];require(new['energy']==old['model']['energy'] and new['forces']==old['model']['forces'],'Exact old/new native E/F parity');checks['exact_raw_E_F']=True
            oldout=verify(old['model']['native']['gradient_output']).read_text();newout=verify(record['actual_output']).read_text()
            patterns=(r'Total Potential Energy\s*:\s*(\S+)',r'(?m)^\s*Anlyt\s+\d+\s+([^\n]+)')
            checks['literal_native_tokens_equal']=all([x.split() for x in re.findall(p,newout)]==[x.split() for x in re.findall(p,oldout)] for p in patterns)
            require(checks['literal_native_tokens_equal'],'Exact native E/F printed tokens')
        else:
            new=record['moments'];base=old['moments'];keys=('energy_kcal_mol','components','component_closure_kcal_mol','charge_e','induced_Debye','induced_sum_Debye','total_COM_Debye')
            for k in keys:require(new[k]==base[k],'Exact raw old/new moment parity: '+k)
            require(record['response_readback']['last_induced_u_v_eA']==old['local_response_readback']['induced_u_v_eA'][-1],'Exact full-precision u/v parity')
            gate(checks,'derived_lab_moment_arithmetic',np.max(abs(np.asarray(new['total_lab_Debye'])-base['total_lab_Debye'])),h['matrix'].arithmetic(3*row['atoms'],new['total_lab_Debye'],base['total_lab_Debye']))
            checks['exact_raw_E_M_uv']=True
            oldout=verify(base['analyze_output']).read_text();newout=verify(record['actual_output']).read_text()
            patterns=(r'Total Potential Energy\s*:\s*(\S+)',r'Total Electric Charge\s*:\s*(\S+)',r'Dipole X,Y,Z-Components\s*:\s*([^\n]+)',r'(?m)^\s*LOCAL_RESPONSE_INDUCED\s+([^\n]+)')
            checks['literal_native_tokens_equal']=all([x.split() for x in re.findall(p,newout)]==[x.split() for x in re.findall(p,oldout)] for p in patterns)
            oldsites=[x.split() for x in verify(base['induced_output']).read_text().splitlines()[1:] if x.strip()];newsites=[x.split() for x in verify(new['induced_output']).read_text().splitlines()[1:] if x.strip()]
            checks['literal_site_tokens_equal']=oldsites==newsites
            require(checks['literal_native_tokens_equal'] and checks['literal_site_tokens_equal'],'Exact native E/M/dual/site printed tokens')
    else:
        require(record['moments']['energy_kcal_mol']==row['baseline_gradient']['energy'],'No-metal/no-water energy bypass');checks['exact_saved_bypass_energy']=True
    return checks

def parse(record,job,row,m,h,ev):
    output=verify(record['actual_output']);rr=copy.deepcopy(row);rr['input']=record['runtime_input'];rr['positions_A']=job.get('runtime_positions_A',row['positions_A'])
    h['gradient'].program_admission(output,job['program'],h['native']);record['induction']=h['induction'].induction(output);record['response_readback']=readback(output,job,rr,m,h,ev)
    if job['program']=='testgrad':
        record['value']=h['gradient'].gradient(output,row['atoms'],job['fd'],h['native']);require(not job['fd'] or record['value']['fd_passed'],'Original Cartesian FD gate')
    elif job['stage']=='directions':record['value']={'energy':h['moments'].parse_energy(output)['energy_kcal_mol'],'output':pin(output)}
    else:
        mom=h['moments'].parse_moments(output,h['moments'].induced_file(output),rr,h['induction'].induction);record['moments']=mom
        paths=sorted(output.parent.glob('response_matrix_*.dat'));require([p.name for p in paths]==['response_matrix_001.dat','response_matrix_002.dat'],'Exactly two native exports')
        record['exports']=[pin(p) for p in paths];record['matrix_qualification']=h['matrix'].qualify_exports(paths,rr,m['units'],mom);require(record['matrix_qualification']['passed'],'Inherited matrix/dual/contraction checks')
        uv=np.asarray(record['response_readback']['last_induced_u_v_eA']);mu=np.asarray(row['permanent_lab_Debye'])+.5*m['units']['debye_per_eA']*np.sum(uv[:,:3]+uv[:,3:],axis=0);checks={}
        gate(checks,'dual_lab_moment',np.max(abs(mu-mom['total_lab_Debye'])),.0005+row['atoms']*.5e-6+2e-9);record['moment_checks']=checks
    if job['profile'] in ('absent','zero') or job['stage']=='bypass':record['zero_identity']=zero_identity(record,row,job,h)

def operator_checks(active,zero,row,m,h):
    dh,_=delta(row,row['e']);checks={}
    for i,(a,z) in enumerate(zip(active['matrix_qualification']['exports'],zero['matrix_qualification']['exports'])):
        require(a['sites']==z['sites'],'Exact operator site map');A=np.asarray(a['A']);Z=np.asarray(z['A']);v=np.asarray(a['vectors']);v0=np.asarray(z['vectors']);order=[s['atom']-1 for s in a['sites']];dv=dh[order].ravel();ar=h['matrix'].arithmetic(len(A),A,Z,v,v0,dv)
        gate(checks,str(i)+'_unchanged_A',np.max(abs(A-Z)),ar)
        for j,label in enumerate(('Ed','Ep')):gate(checks,str(i)+'_shared_delta_'+label,np.max(abs(v[:,j]-v0[:,j]-dv)),ar)
        # This gate is an operator identity, not a fitted chi value.
        gate(checks,str(i)+'_unchanged_H_spectrum',np.max(abs(np.asarray(a['A'])-z['A'])),ar)
    return checks

def field_checks(base,minus,plus,row,m,h):
    units=m['units'];K=units['electric'];B=units['bohr_A'];Ha=units['hartree_kcal_mol'];D=units['debye_per_eA'];n=row['atoms'];dim=3*n;f=.001;c=Ha/(K*B);z=np.tile([0.,0.,1.],n);checks={}
    # E and M are distinct solves: use the E export for energy contraction/bounds.
    b=base['matrix_qualification']['exports'][0];A=np.asarray(b['A']);ed,ep,u,v=np.asarray(b['vectors']).T;ar=h['matrix'].arithmetic(dim,A,ed,ep,u,v)
    ud=np.linalg.solve(A,ed);vd=np.linalg.solve(A,ep);permanent=np.asarray(row['permanent_lab_Debye'])/(D*B);mu=permanent[2]+.5*float(z@(ud+vd))/B;chi=c/B*float(z@np.linalg.solve(A,z))
    e0bound=.5*K*np.linalg.norm(ep)*b['numerical_bounds']['native_vs_dense_u']+5e-9+ar;bounds={}
    for sign,rec in ((-1,minus),(1,plus)):
        a=rec['matrix_qualification']['exports'][0];af=np.asarray(a['A']);efd,efp,uf,vf=np.asarray(a['vectors']).T;tol=h['matrix'].arithmetic(dim,A,af,efd,efp,ed,ep)
        gate(checks,str(sign)+'_same_A',np.max(abs(af-A)),tol)
        for label,x,y in (('Ed',efd,ed),('Ep',efp,ep)):gate(checks,str(sign)+'_field_'+label,np.max(abs(x-y-c*sign*f*z)),tol)
        observed=rec['moments']['energy_kcal_mol']-base['moments']['energy_kcal_mol'];contract=-Ha*sign*f*permanent[2]-.5*K*(float(uf@efp)-float(u@ep));gate(checks,str(sign)+'_native_energy_contraction',abs(observed-contract),1e-8+h['matrix'].arithmetic(dim,observed,contract))
        eb=.5*K*np.linalg.norm(efp)*a['numerical_bounds']['native_vs_dense_u']+5e-9+ar;bounds[sign]=eb;pred=-Ha*(sign*f*mu+.5*f*f*chi);gate(checks,str(sign)+'_dense_field_energy',abs(observed-pred),e0bound+eb)
    epf=plus['moments']['energy_kcal_mol'];emf=minus['moments']['energy_kcal_mol'];e0=base['moments']['energy_kcal_mol'];odd=-(epf-emf)/(2*f*Ha);even=-(epf+emf-2*e0)/(f*f*Ha)
    gate(checks,'odd_moment',abs(odd-mu),(bounds[1]+bounds[-1])/(2*f*Ha));gate(checks,'even_curvature',abs(even-chi),(bounds[1]+bounds[-1]+2*e0bound)/(f*f*Ha))
    return {'checks':checks,'mu_lab_au':mu,'chi_a0_cubed':chi,'odd':odd,'even':even,'E_solve_bound_separate_from_M':True}
