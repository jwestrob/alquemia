#!/usr/bin/env python3
"""One finite numerical panel, reusing admitted parsers and saved baselines."""
import argparse
import concurrent.futures
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import threading
import time
import numpy as np

HERE = Path(__file__).resolve().parent
THREADS = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
PROFILES = {'absent':None,'zero':(0.,0.,0.,0.),'drive':(.05,-.04,0.,0.),'stiffness':(0.,0.,-.08,.06),'combined':(.05,-.04,.08,-.06)}
IONS = ('ion_La','ion_Dy')
PURE = 'aquo9_central_La__near__two_water'
DIMERS = tuple('aquo9_central_'+m+'__dimer_01' for m in ('La','Dy'))
NEAR = tuple('aquo9_central_'+m+'__near__metal_two_water' for m in ('La','Dy'))
CENTERS = tuple('aquo9_central_'+m for m in ('La','Dy'))
METAL = DIMERS+NEAR+CENTERS
CELLS = IONS+(PURE,)+METAL
MAX_CALLS, MAX_FD, WORKERS = 57, 636, 12
NEW_CALLS, INHERITED_CALLS, INHERITED_FD = 30, 27, 468
FILES = (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))

def require(ok,msg):
    if not ok: raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verify(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p): return json.loads(verify(p).read_text())
def save(path,value):
    p=Path(path);t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(t,p)
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def helpers(s):
    return {k:module(s[k],'local_qualification_'+k) for k in ('old_harness','native','gradient','induction','moments','field','matrix')}
def tokens(key):
    return [s.split('#',1)[0].split() for s in key.splitlines() if s.split('#',1)[0].split()]
def scientific_key(key):
    reporting={'polar-print','save-uinduce','response-matrix-export','digits'}
    return [t for t in tokens(key) if t[0].lower() not in reporting]
def clean_key(key):
    require(not any(t[0].lower() in ('polar-local-response','polar-pair-response','external-field','exfld-freq','exchpol','exchange-polar') for t in tokens(key)),'Nonbaseline key')
    return '\n'.join(s for s in key.splitlines() if not s.split() or s.split()[0].lower()!='response-matrix-export')+'\n'

def inventory():
    jobs=[]
    for profile in ('absent','zero','combined'):
        jobs += [dict(id=profile+'__'+g,cell=g,profile=profile,program='testgrad',fd=profile=='combined' and g in METAL,stage=profile) for g in CELLS]
    for profile in ('drive','stiffness'):
        jobs += [dict(id=profile+'__'+g,cell=g,profile=profile,program='testgrad',fd=True,stage='axes') for g in NEAR]
    for g in NEAR:
        for kind in ('radial','rotation'):
            for h in (1e-4,5e-5):
                for sign in (-1,1):
                    jobs.append(dict(id=f'directional__{g}__{kind}__{h:g}__{sign:+d}',cell=g,profile='combined',program='analyze',fd=False,stage='directions',kind=kind,h=h,q=sign*h))
    jobs += [dict(id='matrix__'+g,cell=g,profile='combined',program='analyze',fd=False,stage='matrices',matrix=True) for g in METAL]
    for g in DIMERS:
        for sign in (-1,1):
            jobs.append(dict(id=f'field__{g}__{sign:+d}',cell=g,profile='combined',program='analyze',fd=False,stage='fields',field_au=sign*.001))
    require(len(jobs)==MAX_CALLS and len({j['id'] for j in jobs})==MAX_CALLS,'57 program denominator')
    return jobs

def render(row,job,units):
    key=clean_key(verify(row['input']['key']).read_text())
    if PROFILES[job['profile']] is not None:
        for typ in (500,501): key+='POLAR-LOCAL-RESPONSE '+str(typ)+' 36 '+' '.join(f'{v:.17g}' for v in PROFILES[job['profile']])+'\n'
    if job.get('matrix'):key+='RESPONSE-MATRIX-EXPORT\n'
    if 'field_au' in job:key+=f"EXTERNAL-FIELD 0 0 {job['field_au']*units['MV_cm_per_field_au']:.17g}\n"
    return key

def prepare(out):
    tick,cpu=time.monotonic(),time.process_time();s=read(pin(HERE/'SOURCE_INPUTS.json'))
    for p in s.values():verify(p)
    h=helpers(s);mm,tm=read(s['moment_manifest']),read(s['trimer_manifest']);old=read(s['old_qualification']);data=read(s['development_data'])
    em,er,review=read(s['export_manifest']),read(s['export_result']),read(s['export_review'])
    require(review['status']=='passed_baseline_matrix_inventory_review' and review['result']==s['export_result'],'Full baseline matrix admission')
    require(er['status']=='baseline_matrix_inventory_complete' and er['error'] is None and er['complete_fixtures']==54 and er['admitted_snapshots']==108,'Full baseline matrix result')
    require(old['status']=='metal_pair_numerical_qualification_complete' and old['numerical_checks_passed'],'Saved original E/F authority')
    rows={}
    for g in CELLS:
        row=copy.deepcopy((tm['model_rows'] if g in NEAR+(PURE,) else mm['rows'])[g]);row['energy_replay_kcal_mol']=None
        _,typed,coords=h['old_harness'].xyz(verify(row['input']['xyz']));row['positions_A']=coords
        require(len(typed)==row['atoms'] and [int(t[5]) for t in typed]==row['types'],'Typed fixture inventory')
        if g in CENTERS:
            b=data['baseline_model'][g];p=b['native']['gradient_output'];source_input=b['provenance']['effective_input'];baseline=h['gradient'].gradient(verify(p),row['atoms'],False,h['native'])
            require(baseline['energy']==b['energy'] and baseline['forces']==b['forces'],'Original W9 native gradient replay')
            # Moment extraction changed only title/whitespace. Keep the exact
            # original native gradient XYZ bytes after proving every atom token.
            _,original_atoms,original_coords=h['old_harness'].xyz(verify(source_input['xyz']))
            require(typed==original_atoms and coords==original_coords,'W9 title-only XYZ alias differs physically')
            row['moment_XYZ_alias']=row['input']['xyz'];row['input']['xyz']=source_input['xyz']
        else:
            b=old['rows']['original__'+g];source_input=b['input'];baseline=h['gradient'].gradient(verify(b['value']['output']),row['atoms'],False,h['native'])
            require(baseline==b['value'],'Original small fixture native gradient replay')
        require(all(verify(row['input'][k]).read_bytes()==verify(source_input[k]).read_bytes() for k in ('xyz','parameters')),'Historical geometry/prm equality')
        require(scientific_key(verify(row['input']['key']).read_text())==scientific_key(verify(source_input['key']).read_text()),'Historical physical key equality')
        row['baseline_gradient']=baseline;row['baseline_gradient_input']=source_input
        if g in METAL:
            rr=er['rows'][g];q=rr['admission']['qualification'] if 'admission' in rr else rr['numerical'];pp=q['exports'][1]
            parsed=h['matrix'].parse_export(pp['source'],em['rows'][g],em['units'])
            require(parsed['A']==pp['A'] and parsed['vectors']==pp['vectors'],'Baseline raw matrix replay')
            require(all(verify(row['input'][k]).read_bytes()==verify(em['rows'][g]['input'][k]).read_bytes() for k in ('xyz','parameters')) and scientific_key(verify(row['input']['key']).read_text())==scientific_key(verify(em['rows'][g]['input']['key']).read_text()),'Baseline matrix physical input join')
            row['baseline_matrix']=parsed
            row['permanent_lab_Debye']=h['field'].permanent_moment(row)
        rows[g]=row
    jobs=inventory();require(sum(6*rows[j['cell']]['atoms'] for j in jobs if j['fd'])==MAX_FD,'636 internal FD denominator')
    for job in jobs:
        row=rows[job['cell']];directory=out/'inputs'/job['id'];directory.mkdir(parents=True)
        for k,name in FILES:
            if k!='key':shutil.copyfile(verify(row['input'][k]),directory/name)
        (directory/'system.key').write_text(render(row,job,em['units']))
        if 'q' in job:
            header,atoms,coords=h['old_harness'].xyz(directory/'system.xyz');displaced,tangent=h['old_harness'].path_geometry(coords,job['kind'],job['q'])
            (directory/'system.xyz').write_text('\n'.join([header]+[' '.join(a[:2]+[f'{v:.17g}' for v in p]+a[5:]) for a,p in zip(atoms,displaced)])+'\n')
        job['input']={k:pin(directory/name) for k,name in FILES}
    m=dict(status='prepared_no_native_calls',implementation=pin(__file__),plan=pin(HERE/'QUALIFICATION_PLAN.md'),sources=s,rows=rows,programs=jobs,units=em['units'],maximum_native_programs=MAX_CALLS,maximum_internal_FD_energies=MAX_FD,maximum_workers=WORKERS,new_native_calls=0,cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1})
    save(out/'MANIFEST.json',m)
    save(out/'CONFIG_RUN_TEMPLATE.json',dict(implementation=pin(__file__),manifest=pin(out/'MANIFEST.json'),plan=m['plan'],analyze=None,testgrad=None,source_manifest=None,build_receipt=None,source_audit=None,root_authorized=False,maximum_native_programs=MAX_CALLS,maximum_internal_FD_energies=MAX_FD,workers=WORKERS,fit=False,physical_qualification=False))
    return m

def response_readback(output,job,row):
    text=Path(output).read_text();prefix='LOCAL_RESPONSE_';profile=job['profile'];coef=PROFILES[profile]
    lines=[s.split() for s in text.splitlines() if s.strip().startswith(prefix)]
    if profile in ('absent','zero'):
        require(not lines,'Absent/zero did not use the literal bypass');return {'passed':True,'bypass':True,'records':0}
    sizes={'PARAMETERS':7,'PAIR':6,'ENVIRONMENT':6,'SPECTRUM':6,'INDUCED':8,'DUAL':3};records={k:[] for k in sizes}
    for t in lines:
        label=t[0][len(prefix):];require(label in sizes and len(t)==sizes[label],'Unknown/malformed response readback')
        vals=[float(x.replace('D','E')) for x in t[1:]];require(all(math.isfinite(v) for v in vals),'Nonfinite response readback');records[label].append(vals)
    params=records['PARAMETERS'];require(len(params)>=2 and len(params)%2==0,'Missing coefficient application readback')
    for pos,p in enumerate(params):require(p[:2]==[500+pos%2,36] and p[2:]==list(coef),'Actual coefficient record differs')
    spectra=records['SPECTRUM'];require(spectra and len(spectra)==len(params)//2,'Spectrum per solve denominator')
    for lo,hi,sym,eigres,orth in spectra:require(hi>0 and lo>1e-6*hi and 0<=sym<=1e-10 and 0<=eigres<=1e-10*max(1.,hi) and 0<=orth<=1e-10,'Native operator stability/arithmetic failed')
    require(records['DUAL'] and len(records['INDUCED'])==row['atoms']*len(records['DUAL']),'Native dual/site denominator')
    for start in range(0,len(records['INDUCED']),row['atoms']):require([x[0] for x in records['INDUCED'][start:start+row['atoms']]]==list(range(1,row['atoms']+1)),'Native induced atom order')
    for rhs,uv in records['DUAL']:require(rhs>=0 and uv>=0,'Negative native dual discrepancy')
    require(all(p[2]>0 and p[3]>0 and p[4]>=0 for p in records['PAIR']),'Invalid fixed damping readback')
    if any(t in (500,501) for t in row['types']):require(records['ENVIRONMENT'],'Missing metal environment readback')
    return {'passed':True,'bypass':False,'records':records,'last_induced_eA':[r[1:] for r in records['INDUCED'][-row['atoms']:]],'normalized_eigen_min':min(x[0] for x in spectra),'normalized_eigen_max':max(x[1] for x in spectra)}

def regularized(row,profile,matrix_helper):
    """Independent analytic value reconstruction; no geometry derivatives here."""
    a,b,c,d=PROFILES[profile];xyz=np.array(row['positions_A']);types=row['types'];metal=next(i for i,t in enumerate(types) if t in (500,501));oxy=[i for i,t in enumerate(types) if t==36]
    alpha,_=matrix_helper.settings(row);gamma=None
    for src in ('parameters','key'):
        for t in tokens(verify(row['input'][src]).read_text()):
            if t[0].lower()=='polpair' and {int(t[1]),int(t[2])}=={types[metal],36}:gamma=float(t[3])
    require(gamma is not None and gamma>0,'Resolved fixed POLPAIR missing')
    # Native pdamp=alpha^(1/6): (pdampM*pdampO)^3=sqrt(alphaM*alphaO).
    kappa=gamma/math.sqrt(alpha[types[metal]]*alpha[36]);f=[];vectors=[];radii=[]
    anchor=(1+5+25/3)*math.exp(-5)
    for i in oxy:
        v=xyz[i]-xyz[metal];r=float(np.linalg.norm(v));t=2*r;f.append(((1+t+t*t/3)*math.exp(-t)/anchor)**2);vectors.append(v);radii.append(r)
    s=math.fsum(f);other=[math.fsum(f[j] for j in range(len(f)) if j!=i) for i in range(len(f))];eta=math.fsum(f[i]*other[i] for i in range(len(f)))/(1+s*s);g=2*s/(1+s);sigma=g*(c+d*eta);beta=math.expm1(sigma)/alpha[types[metal]]
    hv=np.zeros_like(xyz)
    for i,fi,v,r,other_i in zip(oxy,f,vectors,radii,other):
        phi=kappa if r==0 else -math.expm1(-kappa*r**3)/r**3
        hv[i]=3*fi*(a+b*other_i/(1+other_i))*phi*v
    return dict(h=hv,metal=metal,beta=beta,s=s,eta=eta,sigma=sigma,gamma=gamma,kappa=kappa,radii=radii)

def gate(checks,name,error,bound):
    passed=math.isfinite(float(error)) and math.isfinite(float(bound)) and 0<=error<=bound
    checks[name]=dict(error=float(error),bound=float(bound),passed=bool(passed));require(passed,'Numerical check failed: '+name)

def matrix_checks(record,row,units,h):
    q=h['matrix'].qualify_exports(record['exports'],row,units,record['moments']);record['matrix_qualification']=q
    z=q['exports'][1];baseline=row['baseline_matrix'];v=np.array(z['vectors']);oldv=np.array(baseline['vectors']);A=np.array(z['A']);A0=np.array(baseline['A']);expected=regularized(row,'combined',h['matrix'])
    sites=[s['atom']-1 for s in z['sites']];require(sites==[s['atom']-1 for s in baseline['sites']],'Baseline/candidate export site order')
    H=np.zeros_like(A);p=sites.index(expected['metal']);H[3*p:3*p+3,3*p:3*p+3]=np.eye(3)*expected['beta'];hh=expected['h'][sites].reshape(-1);tol=h['matrix'].arithmetic(len(A),A,A0,v,oldv,hh)
    checks={};gate(checks,'actual_A_equals_baseline_plus_H',np.max(abs(A-A0-H)),tol)
    for i,name in enumerate(('Ed','Ep')):gate(checks,'actual_'+name+'_equals_baseline_plus_h',np.max(abs(v[:,i]-oldv[:,i]-hh)),tol)
    rb=record['response_readback']['records'];last=rb['ENVIRONMENT'][-1];gate(checks,'native_environment',np.max(abs(np.array(last[1:])-[expected[k] for k in ('s','eta','sigma','beta')])),tol)
    pairs=rb['PAIR'][-len(expected['radii']):]
    require(len(pairs)==len(expected['radii']),'Pair readback denominator')
    for i,(p,r) in enumerate(zip(pairs,expected['radii'])):gate(checks,'native_fixed_pair_'+str(i),np.max(abs(np.array(p[2:])-[expected['gamma'],expected['kappa'],r])),tol)
    native_uv=np.array(record['response_readback']['last_induced_eA']);mapped=np.array([native_uv[i] for i in sites]);gate(checks,'full_precision_dual_dump',np.max(abs(mapped[:,:3].reshape(-1)-v[:,2])),tol);gate(checks,'full_precision_dual_dump_p',np.max(abs(mapped[:,3:].reshape(-1)-v[:,3])),tol)
    permanent=np.array(row['permanent_lab_Debye']);mu=permanent+units['debye_per_eA']*.5*(native_uv[:,:3]+native_uv[:,3:]).sum(axis=0)
    gate(checks,'native_conjugate_M_lab',np.max(abs(mu-record['moments']['total_lab_Debye'])),.0005+2e-9)
    record['operator_application']={'checks':checks,'expected_beta':expected['beta'],'expected_h_norm':float(np.linalg.norm(hh)),'conjugate_lab_Debye':mu.tolist(),'passed':all(c['passed'] for c in checks.values())}

def field_checks(g,records,row,units,h):
    base=records['matrix__'+g];p=base['matrix_qualification']['exports'][1];A=np.array(p['A']);ed,ep,u0,v0=np.array(p['vectors']).T;n=row['atoms'];dim=3*n;K=units['electric'];B=units['bohr_A'];Ha=units['hartree_kcal_mol'];D=units['debye_per_eA'];coupling=Ha/(K*B);z=np.tile([0.,0.,1.],n)
    ar=h['matrix'].arithmetic(dim,A,ed,ep,u0,v0);an=math.sqrt(dim)*ar;lam=p['spectrum']['A'][0];amp=p['numerical_bounds']['post_peek_amplification'];permanent=np.array(row['permanent_lab_Debye'])/(D*B)
    denseu=np.linalg.solve(A,ed);densev=np.linalg.solve(A,ep);mu0=permanent[2]+.5*float(z@(denseu+densev))/B;chi=coupling/B*float(z@np.linalg.solve(A,z));checks={};bounds={}
    u0err=p['numerical_bounds']['native_vs_dense_u'];v0err=p['numerical_bounds']['native_vs_dense_v'];e0bound=.5*K*np.linalg.norm(ep)*u0err+5e-9+ar
    for sign in (-1,1):
        rec=records[f'field__{g}__{sign:+d}'];f=sign*.001;native=np.array(rec['response_readback']['last_induced_eA']);native=native[[s['atom']-1 for s in p['sites']]];uf=native[:,:3].reshape(-1);vf=native[:,3:].reshape(-1);df=ed+coupling*f*z;pf=ep+coupling*f*z
        q=max(h['matrix'].induction_values(rec['moments']['analyze_output']));resbound=amp*math.sqrt(n)*(q+.5e-10)/D+an
        gate(checks,f'{sign}_actual_Au_field_RHS',np.linalg.norm(A@uf-df),resbound);gate(checks,f'{sign}_actual_Av_field_RHS',np.linalg.norm(A@vf-pf),resbound)
        ubound=(np.linalg.norm(A@uf-df)+2*an)/lam;vbound=(np.linalg.norm(A@vf-pf)+2*an)/lam
        ordered_sites=np.array(rec['moments']['induced_Debye'])[[s['atom']-1 for s in p['sites']]]
        gate(checks,f'{sign}_saved_site_u',np.max(abs(D*uf.reshape(n,3)-ordered_sites)),.5e-6+D*ar)
        muf=permanent*D*B+D*.5*(uf.reshape(n,3)+vf.reshape(n,3)).sum(axis=0)
        gate(checks,f'{sign}_native_conjugate_M',np.max(abs(muf-rec['moments']['total_lab_Debye'])),.0005+2e-9)
        observed=rec['moments']['energy_kcal_mol']-base['moments']['energy_kcal_mol'];contract=-Ha*f*permanent[2]-.5*K*(float(uf@pf)-float(u0@ep))
        gate(checks,f'{sign}_field_energy_contract',abs(observed-contract),1e-8+h['matrix'].arithmetic(dim,observed,contract))
        predicted=-Ha*(f*mu0+.5*f*f*chi);ebound=.5*K*np.linalg.norm(pf)*ubound+5e-9+ar
        gate(checks,f'{sign}_field_energy_dense',abs(observed-predicted),e0bound+ebound);bounds[sign]=dict(energy=ebound,u=ubound,v=vbound)
    plus=records[f'field__{g}__+1']['moments']['energy_kcal_mol'];minus=records[f'field__{g}__-1']['moments']['energy_kcal_mol'];zero=base['moments']['energy_kcal_mol'];f=.001
    odd=-(plus-minus)/(2*f*Ha);even=-(plus+minus-2*zero)/(f*f*Ha)
    gate(checks,'odd_energy_conjugate_moment',abs(odd-mu0),(bounds[1]['energy']+bounds[-1]['energy'])/(2*f*Ha))
    gate(checks,'even_energy_matrix_curvature',abs(even-chi),(bounds[1]['energy']+bounds[-1]['energy']+2*e0bound)/(f*f*Ha))
    return dict(checks=checks,mu0_lab_au=mu0,odd_energy_mu_au=odd,matrix_chi_a0_cubed=chi,even_energy_chi_a0_cubed=even,passed=all(x['passed'] for x in checks.values()))

def context(cfg,runtime):
    require(cfg['implementation']==pin(__file__) and cfg['plan']==pin(HERE/'QUALIFICATION_PLAN.md'),'Qualification source/plan')
    m=read(cfg['manifest']);require(m['implementation']==cfg['implementation'] and m['plan']==cfg['plan'] and m['new_native_calls']==0,'Preparation binding')
    require(cfg['root_authorized'] is True and cfg['fit'] is False and cfg['physical_qualification'] is False,'Numerical-only root binding required')
    for k in ('analyze','testgrad','source_manifest','build_receipt','source_audit'):verify(cfg[k])
    audit=read(cfg['source_audit']);require(audit['status']=='passed_local_response_native_source_review' and audit['source_manifest']==cfg['source_manifest'],'Exact native source review join')
    build=read(cfg['build_receipt']);require(build['status']=='complete' and build['manifest']==cfg['source_manifest'] and all(build['executables'][k]==cfg[k] for k in ('analyze','testgrad')),'Completed isolated build/source/executables')
    require(cfg['maximum_native_programs']==MAX_CALLS and cfg['maximum_internal_FD_energies']==MAX_FD and cfg['workers']==WORKERS,'Fixed native limits')
    require(set(m['rows'])==set(CELLS) and [{k:v for k,v in j.items() if k!='input'} for j in m['programs']]==inventory(),'Full exact prospective inventory')
    for p in m['sources'].values():verify(p)
    for j in m['programs']:
        for p in j['input'].values():verify(p)
        require(verify(j['input']['key']).read_text()==render(m['rows'][j['cell']],j,m['units']),'Prepared key edit differs')
    if runtime:require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==WORKERS and all(os.environ.get(k)=='1' for k in THREADS),'One 12-CPU Slurm task; every numerical library one thread')
    require(cfg['maximum_new_native_programs']==NEW_CALLS and cfg['recovery_declaration']==pin(HERE/'RECOVERY_DECLARATION.md') and m['recovery_declaration']==cfg['recovery_declaration'],'Exact recovery scope')
    require(cfg['previous_result']==m['previous_result'] and len(m['reused_rows'])==INHERITED_CALLS and set(m['reused_rows'])=={j['id'] for j in m['programs'] if j['stage'] in ('absent','zero','combined')},'Exact original reused prefix')
    old=read(m['previous_result']);require(all(m['reused_rows'][k]==old['rows'][k] for k in m['reused_rows']),'Inherited records changed')
    review=read(cfg['prefix_review']);require(cfg['prefix_review']==m['prefix_review'] and review['status']=='passed_completed_native_prefix_review' and review['result']==m['previous_result'] and review['completed_programs']==INHERITED_CALLS and review['completed_internal_FD_energies']==INHERITED_FD and review['native_derivative_failure_found'] is False,'Independent prefix admission')
    require(all(review[k]==cfg[k] for k in ('source_manifest','build_receipt','analyze','testgrad')),'Prefix build/executable continuity')
    return m,helpers(m['sources'])

def run(cfg,out):
    tick,cpu=time.monotonic(),time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN);m,h=context(cfg,True);lock=threading.RLock();stop=threading.Event()
    result=dict(status='running',implementation=pin(__file__),config=cfg,manifest=cfg['manifest'],job_id=os.environ['SLURM_JOB_ID'],rows={j['id']:copy.deepcopy(m['reused_rows'].get(j['id'],{'status':'not_attempted'})) for j in m['programs']},native_program_invocations=0,inherited_native_program_invocations=INHERITED_CALLS,total_native_program_invocations=INHERITED_CALLS,previous_result=m['previous_result'],internal_FD_energies_requested=INHERITED_FD,internal_FD_energies_completed=INHERITED_FD,new_internal_FD_energies_completed=0,inherited_internal_FD_energies_completed=INHERITED_FD,directional_checks={},field_checks={},global_stop=False,error=None,numerical_checks_passed=False,physical_accuracy_qualified=False)
    def persist():
        with lock:result['global_stop']=stop.is_set();result['total_native_program_invocations']=INHERITED_CALLS+result['native_program_invocations'];save(out/'RESULT.json',result)
    def worker(job):
        row=m['rows'][job['cell']];cell=out/'cells'/job['id'];rec=dict(status='started',job=job,receipts=[])
        try:
            with lock:
                if stop.is_set():return
                fd=6*row['atoms'] if job['fd'] else 0
                require(result['native_program_invocations']<NEW_CALLS and INHERITED_CALLS+result['native_program_invocations']<MAX_CALLS and result['internal_FD_energies_requested']+fd<=MAX_FD,'Hard native ceilings')
                result['native_program_invocations']+=1;result['internal_FD_energies_requested']+=fd;result['rows'][job['id']]=copy.deepcopy(rec);persist()
            cell.mkdir(parents=True)
            for k,name in FILES:shutil.copyfile(verify(job['input'][k]),cell/name)
            rec['runtime_input']={k:pin(cell/name) for k,name in FILES}
            if job['program']=='testgrad':argv=[cfg['testgrad']['path'],'system.xyz','Y','Y','0.00001'] if job['fd'] else [cfg['testgrad']['path'],'system.xyz','Y','N']
            else:argv=[cfg['analyze']['path'],'system.xyz','EM' if job.get('matrix') or 'field_au' in job else 'E']
            output=cell/(job['program']+'.out');h['native'].invoke(argv,cell,output,rec['receipts'])
            rec['response_readback']=response_readback(output,job,row);rec['induction']=h['induction'].induction(output)
            if job['program']=='testgrad':
                value=h['gradient'].gradient(output,row['atoms'],job['fd'],h['native']);rec['value']=value
                with lock:result['internal_FD_energies_completed']+=fd;result['new_internal_FD_energies_completed']+=fd
                require(math.isfinite(value['energy']) and (not job['fd'] or value['fd_passed']),'Cartesian derivative gate')
                if job['profile'] in ('absent','zero') or job['cell'] in IONS+(PURE,):
                    rec['identity']=h['old_harness'].identity(value,row['baseline_gradient']);require(rec['identity']['passed'],'Baseline E/F identity')
            elif job.get('matrix') or 'field_au' in job:
                rr=copy.deepcopy(row);rr['input']=rec['runtime_input'];rec['moments']=h['moments'].parse_moments(output,h['moments'].induced_file(output),rr,h['induction'].induction)
                if job.get('matrix'):
                    rec['exports']=[pin(cell/f'response_matrix_{i:03d}.dat') for i in (1,2)];matrix_checks(rec,rr,m['units'],h)
                    central=result['rows']['combined__'+job['cell']]['value']['energy'];delta=abs(central-rec['moments']['energy_kcal_mol']);rec['energy_route_check']={'error':delta,'bound':1e-4,'passed':delta<=1e-4};require(delta<=1e-4,'Analyze/testgrad energy route closure')
            else:
                h['gradient'].program_admission(output,'analyze',h['native']);e=h['moments'].parse_energy(output)['energy_kcal_mol'];rec['value']={'energy':e,'output':pin(output)}
            rec['status']='complete'
        except Exception as exc:
            stop.set();rec.update(status='failed',error=repr(exc))
            with lock:
                if result['error'] is None:result['error']=repr(exc)
        finally:
            with lock:
                if cell.exists():
                    for name in ('testgrad.out','analyze.out'):
                        if (cell/name).exists():rec['actual_output']=pin(cell/name)
                    save(cell/'RESULT.json',rec)
                if rec['status']!='started' or cell.exists():result['rows'][job['id']]=rec
                persist()
    try:
        persist()
        for stage in ('absent','zero','combined','axes','directions','matrices','fields'):
            jobs=[j for j in m['programs'] if j['stage']==stage]
            with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:list(pool.map(worker,[j for j in jobs if j['id'] not in m['reused_rows']]))
            require(not stop.is_set() and all(result['rows'][j['id']]['status']=='complete' for j in jobs),'Failed stage: '+stage)
            if stage=='directions':
                for g in NEAR:
                    for kind in ('radial','rotation'):
                        _,tangent=h['old_harness'].path_geometry(m['rows'][g]['positions_A'],kind,0.);gradient=result['rows']['combined__'+g]['value']['gradient'];analytic=math.fsum(v*t for gr,tr in zip(gradient,tangent) for v,t in zip(gr,tr));steps=[]
                        for step in (1e-4,5e-5):
                            minus=result['rows'][f'directional__{g}__{kind}__{step:g}__-1']['value']['energy'];plus=result['rows'][f'directional__{g}__{kind}__{step:g}__+1']['value']['energy'];slope=(plus-minus)/(2*step);bound=1e-8/(2*step)+5e-9*sum(abs(v) for t in tangent for v in t);checks={};gate(checks,'slope',abs(slope-analytic),.001+bound);steps.append(dict(half_step_A=step,FD=slope,analytic=analytic,printing_bound=bound,checks=checks))
                        checks={};gate(checks,'two_step',abs(steps[0]['FD']-steps[1]['FD']),.001+sum(s['printing_bound'] for s in steps));result['directional_checks'][g+'__'+kind]={'steps':steps,'checks':checks,'passed':True}
            if stage=='fields':
                for g in DIMERS:result['field_checks'][g]=field_checks(g,result['rows'],m['rows'][g],m['units'],h)
            persist()
        require(result['native_program_invocations']==NEW_CALLS and result['total_native_program_invocations']==MAX_CALLS and result['new_internal_FD_energies_completed']==MAX_FD-INHERITED_FD and result['internal_FD_energies_completed']==MAX_FD and all(r['status']=='complete' for r in result['rows'].values()),'Complete finite denominator')
        result.update(status='local_response_numerical_qualification_complete',numerical_checks_passed=True)
    except Exception as exc:stop.set();result.update(status='local_response_numerical_qualification_failed',error=result['error'] or repr(exc))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']=dict(wall_seconds=time.monotonic()-tick,parent_CPU_seconds=time.process_time()-cpu,aggregate_child_user_CPU_seconds=ch.ru_utime-ch0.ru_utime,aggregate_child_system_CPU_seconds=ch.ru_stime-ch0.ru_stime,allocated_CPUs=WORKERS,per_receipt_CPU_deltas_not_summed=True);persist()
    return result


def prepare_recovery(cfg,out):
    tick,cpu=time.monotonic(),time.process_time()
    old_impl=pin(HERE/'RUN.py');old=module(old_impl,'original_local_qualification');r=read(cfg['previous_result']);m,h=old.context(r['config'],False)
    require(r['implementation']==old_impl and r['status']=='local_response_numerical_qualification_failed' and r['error']=="RuntimeError('dictionary changed size during iteration')",'Exact original technical failure required')
    require(r['native_program_invocations']==INHERITED_CALLS and r['internal_FD_energies_completed']==INHERITED_FD and r['internal_FD_energies_requested']==INHERITED_FD,'Original completed counters')
    expected=[j for j in m['programs'] if j['stage'] in ('absent','zero','combined')];require(len(expected)==INHERITED_CALLS and set(r['rows'])=={j['id'] for j in m['programs']},'Original full denominator')
    reused={};checks=[]
    for job in m['programs']:
        rec=r['rows'][job['id']]
        if job not in expected:
            require(rec=={'status':'not_attempted'},'Unexpected later original invocation');continue
        require(rec['status']=='complete' and rec['job']==job,'Original completed prefix')
        folder=Path(rec['runtime_input']['xyz']['path']).parent
        require(read(pin(folder/'RESULT.json'))==rec,'Individual versus aggregate receipt mismatch')
        require(all(verify(rec['runtime_input'][k]).read_bytes()==verify(job['input'][k]).read_bytes() for k,_ in FILES),'Original actual prepared bytes')
        receipts=rec['receipts'];argv=[r['config']['testgrad']['path'],'system.xyz','Y','Y','0.00001'] if job['fd'] else [r['config']['testgrad']['path'],'system.xyz','Y','N']
        require(len(receipts)==1 and receipts[0]['argv']==argv and receipts[0]['cwd']==str(folder) and receipts[0]['returncode']==0 and all(receipts[0]['threads'][k]=='1' for k in THREADS),'Original invocation receipt')
        output=verify(rec['value']['output']);require(receipts[0]['output']==rec['value']['output'] and rec['actual_output']==rec['value']['output'],'Original output authority')
        row=m['rows'][job['cell']];value=h['gradient'].gradient(output,row['atoms'],job['fd'],h['native']);require(value==rec['value'] and (not job['fd'] or value['fd_passed']),'Original raw E/F/FD admission')
        require(response_readback(output,job,row)==rec['response_readback'] and h['induction'].induction(output)==rec['induction'],'Original readback/induction replay')
        if 'identity' in rec:require(h['old_harness'].identity(value,row['baseline_gradient'])==rec['identity'] and rec['identity']['passed'],'Original baseline identity replay')
        reused[job['id']]=rec;checks.append({'id':job['id'],'raw_admission_passed':True,'FD_passed':value['fd_passed'] if job['fd'] else None})
    review=read(cfg['prefix_review']);require(review['status']=='passed_completed_native_prefix_review' and review['result']==cfg['previous_result'] and review['completed_programs']==INHERITED_CALLS and review['completed_internal_FD_energies']==INHERITED_FD and review['native_derivative_failure_found'] is False,'Independent prefix review binding')
    recovery=copy.deepcopy(m);recovery.update(implementation=pin(__file__),original_implementation=old_impl,original_manifest=r['manifest'],previous_result=cfg['previous_result'],prefix_review=cfg['prefix_review'],recovery_declaration=pin(HERE/'RECOVERY_DECLARATION.md'),reused_rows=reused,reuse_checks=checks,maximum_new_native_programs=NEW_CALLS,inherited_native_programs=INHERITED_CALLS,inherited_internal_FD_energies=INHERITED_FD,cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1})
    save(out/'MANIFEST.json',recovery)
    template=dict(r['config']);template.update(implementation=pin(__file__),manifest=pin(out/'MANIFEST.json'),root_authorized=False,previous_result=cfg['previous_result'],prefix_review=cfg['prefix_review'],recovery_declaration=pin(HERE/'RECOVERY_DECLARATION.md'),maximum_new_native_programs=NEW_CALLS)
    save(out/'CONFIG_RUN_TEMPLATE.json',template);return recovery

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare-recovery','run'));p.add_argument('--config');p.add_argument('--out',required=True);a=p.parse_args();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:r=prepare_recovery(json.loads(Path(a.config).read_text()),out) if a.action=='prepare-recovery' else run(json.loads(Path(a.config).read_text()),out)
    except Exception as exc:save(out/'FAILURE.json',{'error':repr(exc),'new_native_calls':0,'phase':a.action});raise
    print(json.dumps({'status':r['status'],'out':str(out)}));raise SystemExit(0 if r['status'] in ('prepared_no_native_calls','local_response_numerical_qualification_complete') else 1)
