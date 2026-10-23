#!/usr/bin/env python3
"""Actual 57-output numerical review, saved files only, one local thread."""
import copy
import argparse
import os
from decimal import Decimal
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import time
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
if not all(os.environ.get(k)=='1' for k in THREADS):raise RuntimeError('One numerical thread required')
import numpy as np

ROOT=Path.cwd()
REL=Path('explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/sampling_cpu_refactor_v1')
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/REL
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/REL
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def canon(x):return json.dumps(x,sort_keys=True,allow_nan=False)
def main(config_path,dest):
    tick,cpu=time.monotonic(),time.process_time();binding=json.loads(config_path.read_text())
    for name,p in binding.items():
        if not isinstance(p,dict) or set(p)!={'path','sha256'}:raise ValueError('Unbound/non-pin review input '+name)
        if pin(p['path'])!=p:raise ValueError('Changed review input '+name)
    sp=importlib.util.spec_from_file_location('actual_local_response',binding['runner']['path']);q=importlib.util.module_from_spec(sp);sp.loader.exec_module(q)
    path=Path(binding['result']['path']);r=json.loads(path.read_text());cfg=r['config'];m,h=q.context(cfg,False);checks=[]

    def ck(ok,name):
        checks.append(name)
        if not ok:raise ValueError(name)
    derived_differences=[]
    def derived(a,b,label,dim):
        if isinstance(a,dict):
            ck(set(a)==set(b),label+' same keys')
            for k in a:derived(a[k],b[k],label+'/'+str(k),dim)
        elif isinstance(a,(list,tuple)):
            ck(len(a)==len(b),label+' same length')
            for i,(x,y) in enumerate(zip(a,b)):derived(x,y,label+'/'+str(i),dim)
        elif isinstance(a,(float,np.floating)):
            error=abs(float(a)-float(b));bound=h['matrix'].arithmetic(dim,a,b)
            # A field curvature bound is a four-energy error combination
            # divided by Ha*h^2. Preserve that existing propagation when
            # comparing only the re-evaluated bound across BLAS hosts.
            if '/field/checks/' in label and label.endswith('/bound'):
                bound*=max(1.,4/(m['units']['hartree_kcal_mol']*.001**2))
            ck(math.isfinite(error) and error<=bound,label+' inherited arithmetic replay')
            if error:derived_differences.append({'field':label,'error':error,'arithmetic_bound':float(bound)})
        else:ck(a==b,label+' exact nonarithmetic value')
    ck(r['implementation']==binding['runner'] and r['manifest']==cfg['manifest'],'actual source/manifest binding')
    ck(q.read(binding['configuration'])==cfg and cfg['source_manifest']==binding['source_manifest'] and cfg['build_receipt']==binding['build_receipt'] and cfg['source_audit']==binding['source_review'],'exact actual configuration/build/source audit')
    login=q.read(binding['login_review']);ck(login['status']=='passed_sampling_cpu_serial_qualification_source_review' and login['implementation']==binding['runner'] and login['source_manifest']==binding['source_manifest'] and login['resource_context']==cfg['resource_context'],'independent serial-login source review')
    ck(r['job_id'] is None and r['execution_mode']=='login' and r['hostname']=='biotite' and isinstance(r['controller_pid'],int) and r['controller_pid']>0 and r['resource_context']==cfg['resource_context'],'truthful actual login identity')
    ck(r['cost']['native_workers']==1 and r['cost']['maximum_owned_numerical_threads']==2 and r['cost']['scheduler_allocation'] is None and r['cost']['scheduler_allocated_CPU_seconds']==0,'truthful serial resource and cost ledger')
    ck(all(math.isfinite(r['cost'][k]) and r['cost'][k]>=0 for k in ('wall_seconds','parent_CPU_seconds','aggregate_child_user_CPU_seconds','aggregate_child_system_CPU_seconds')),'finite nonnegative measured execution costs')
    ck(r['status']=='local_response_numerical_qualification_complete' and r['error'] is None and r['global_stop'] is False and r['numerical_checks_passed'] is True,'actual terminal numerical status')
    ck(r['native_program_invocations']==57 and not m.get('reused_rows'),'57 fresh programs without inherited prefix')
    ck(r['internal_FD_energies_requested']==r['internal_FD_energies_completed']==636,'all636 requested/completed internal FD energies')
    jobs={j['id']:j for j in m['programs']};ck(set(r['rows'])==set(jobs) and len(jobs)==57,'exact all57 inventory')
    rows={};energies={};raw_text={};fd=0;identities=0;induction_n=0;induction_max=0.;fd_max=Decimal(0);fd_rows=0;exports_n=0;route_max=0.;spectrum=[];operator_checks={}
    for cid,job in jobs.items():
        rec=r['rows'][cid];row=m['rows'][job['cell']];out=q.verify(rec['actual_output']);folder=out.parent;raw=out.read_text();raw_text[cid]=raw
        ck(rec['status']=='complete' and rec['job']==job and len(rec['receipts'])==1,cid+' complete one-call record')
        ck(json.loads((folder/'RESULT.json').read_text())==rec,cid+' per-cell/aggregate equality')
        ck(folder==path.parent/'cells'/cid,cid+' actual fresh producer directory')
        for key,name in q.FILES:
            ck(pin(folder/name)==rec['runtime_input'][key] and (folder/name).read_bytes()==q.verify(job['input'][key]).read_bytes(),cid+' literal runtime '+key)
        receipt=rec['receipts'][0]
        args=['system.xyz','Y']+(['Y','0.00001'] if job['fd'] else ['N']) if job['program']=='testgrad' else ['system.xyz','EM' if job.get('matrix') or 'field_au' in job else 'E']
        ck(receipt['argv']==[cfg[job['program']]['path']]+args and receipt['cwd']==str(folder) and receipt['returncode']==0 and receipt['output']==rec['actual_output'],cid+' exact successful native receipt')
        ck(receipt['threads']=={k:'1' for k in q.THREADS} and receipt['timeout_seconds']==300,cid+' one-thread native receipt')
        ck(not re.search(r'\b(?:NaN|Infinity)\b',raw,re.I),cid+' raw finite tokens')
        rb=q.response_readback(out,job,row);ind=h['induction'].induction(out)
        ck(canon(rb)==canon(rec['response_readback']) and ind==rec['induction'],cid+' full readback/induction replay')
        induction_n+=ind['solutions'];induction_max=max(induction_max,ind['maximum_printed_RMS_residual_Debye'])
        if job['program']=='testgrad':
            value=h['gradient'].gradient(out,row['atoms'],job['fd'],h['native']);ck(value==rec['value'],cid+' raw gradient parser replay')
            section=raw.rsplit('Cartesian Gradient Breakdown over Individual Atoms :',1)[-1];tokens=re.findall(r'^\s*(Anlyt|Numer)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s*$',section,re.M);vectors={'Anlyt':{},'Numer':{}}
            for kind,i,*v in tokens:
                ck(int(i) not in vectors[kind],cid+' unique gradient '+kind+i);vectors[kind][int(i)]=[Decimal(s.replace('D','E')) for s in v]
            ck(set(vectors['Anlyt'])==set(range(1,row['atoms']+1)),cid+' complete independent analytic sites')
            ck([[float(x) for x in vectors['Anlyt'][i]] for i in range(1,row['atoms']+1)]==value['gradient'],cid+' independent analytic token replay')
            err=None
            if job['fd']:
                ck(set(vectors['Numer'])==set(vectors['Anlyt']),cid+' complete independent numerical sites')
                err=max(abs(a-b) for i in vectors['Anlyt'] for a,b in zip(vectors['Anlyt'][i],vectors['Numer'][i]));ck(err<=Decimal('.001') and value['fd_passed'],cid+' independent Cartesian derivative gate')
                fd+=6*row['atoms'];fd_rows+=1;fd_max=max(fd_max,err)
            else:ck(not vectors['Numer'],cid+' no undeclared finite differences')
            if 'identity' in rec:
                identities+=1;ck(value['energy']==row['baseline_gradient']['energy'] and value['forces']==row['baseline_gradient']['forces'],cid+' exact original E/F identity')
            rows[cid]={'output':pin(out),'FD_error_decimal':str(err) if err is not None else None,'induction':ind}
        else:
            es=re.findall(r'Total Potential Energy\s*:\s*(\S+)\s+Kcal/mole',raw);ck(len(es)==1,cid+' exact single energy token');energies[cid]=Decimal(es[0])
            if job.get('matrix') or 'field_au' in job:
                rr=copy.deepcopy(row);rr['input']=rec['runtime_input'];mom=h['moments'].parse_moments(out,h['moments'].induced_file(out),rr,h['induction'].induction)
                for name,value in mom.items():
                    if name=='total_lab_Debye':derived(value,rec['moments'][name],cid+'/moments/'+name,3*row['atoms'])
                    else:ck(canon(value)==canon(rec['moments'][name]),cid+' exact raw moment '+name)
                ck(set(mom)==set(rec['moments']),cid+' complete moment schema')
                if job.get('matrix'):
                    ck(sorted(p.name for p in folder.glob('response_matrix_*.dat'))==['response_matrix_001.dat','response_matrix_002.dat'],'two immutable exports '+cid)
                    replay=copy.deepcopy(rec);q.matrix_checks(replay,rr,m['units'],h)
                    for fresh,saved in zip(replay['matrix_qualification']['exports'],rec['matrix_qualification']['exports']):
                        for key in ('source','sequence','sites','vectors','A','checks','restored','peek'):ck(canon(fresh[key])==canon(saved[key]),cid+' exact raw export '+key)
                    derived(replay['matrix_qualification'],rec['matrix_qualification'],cid+'/matrix',3*row['atoms'])
                    derived(replay['operator_application'],rec['operator_application'],cid+'/operator',3*row['atoms'])
                    delta=abs(r['rows']['combined__'+job['cell']]['value']['energy']-mom['energy_kcal_mol']);ck(delta<=1e-4 and delta==rec['energy_route_check']['error'],cid+' analyze/gradient energy route');route_max=max(route_max,delta)
                    independent=[]
                    for ex in rec['matrix_qualification']['exports']:
                        exports_n+=1;A=np.array(ex['A']);V=np.array(ex['vectors']);ed,ep,u,v=V.T;sites=ex['sites'];alpha=np.repeat([s['alpha_A3'] for s in sites],3);ar=ex['numerical_bounds']['arithmetic'];eig=np.linalg.eigvalsh(np.sqrt(alpha)[:,None]*A*np.sqrt(alpha)[None,:]);spectrum.append([float(eig[0]),float(eig[-1])]);ck(eig[0]>1e-6*eig[-1],cid+' independent normalized SPD')
                        L=np.linalg.cholesky(A);solve=lambda b:np.linalg.solve(L.T,np.linalg.solve(L,b));du,dv=solve(ed),solve(ep)
                        ck(np.linalg.norm(du-u)<=ex['numerical_bounds']['native_vs_dense_u'] and np.linalg.norm(dv-v)<=ex['numerical_bounds']['native_vs_dense_v'],cid+' independent Cholesky dual solve')
                        e=-.5*m['units']['electric']*math.fsum(float(x*y) for x,y in zip(u,ep));epout=mom['components']['Polarization']['energy_kcal_mol'];ck(abs(e-epout)<=5e-9+ar,cid+' stationary polarization energy')
                        # Independently reconstruct the fixed regularized drive and
                        # metal stiffness, using only raw native sites and key gamma.
                        a,b,c,d=q.PROFILES['combined'];mi=next(i for i,s in enumerate(sites) if s['type'] in (500,501));oi=[i for i,s in enumerate(sites) if s['type']==36];xyz=np.array([s['xyz_A'] for s in sites]);gamma=rec['response_readback']['records']['PAIR'][0][2];fs=[];rrs=[]
                        anchor=(1+5+25/3)*math.exp(-5)
                        for i in oi:
                            radius=float(np.linalg.norm(xyz[i]-xyz[mi]));rrs.append(radius);t=2*radius;fs.append(((1+t+t*t/3)*math.exp(-t)/anchor)**2)
                        s=math.fsum(fs);eta=math.fsum(fi*fj for i,fi in enumerate(fs) for j,fj in enumerate(fs) if i!=j)/(1+s*s);sigma=2*s/(1+s)*(c+d*eta);beta=math.expm1(sigma)/sites[mi]['alpha_A3'];H=np.zeros_like(A);H[3*mi:3*mi+3,3*mi:3*mi+3]=np.eye(3)*beta;hv=np.zeros(len(ed))
                        for k,i in enumerate(oi):
                            other=math.fsum(f for j,f in enumerate(fs) if j!=k);kap=gamma/math.sqrt(sites[mi]['alpha_A3']*sites[i]['alpha_A3']);radius=rrs[k];phi=-math.expm1(-kap*radius**3)/radius**3
                            hv[3*i:3*i+3]=3*fs[k]*(a+b*other/(1+other))*phi*(xyz[i]-xyz[mi])
                        baseline=row['baseline_matrix'];A0=np.array(baseline['A']);V0=np.array(baseline['vectors']);bound=rec['operator_application']['checks']['actual_A_equals_baseline_plus_H']['bound']
                        errA=float(np.max(abs(A-A0-H)));errd=float(np.max(abs(ed-V0[:,0]-hv)));errp=float(np.max(abs(ep-V0[:,1]-hv)))
                        ck(max(errA,errd,errp)<=bound,cid+' independent regularized A/Ed/Ep application')
                        independent.append(dict(sequence=ex['sequence'],A_error=errA,Ed_error=errd,Ep_error=errp,polarization_energy_error=abs(e-epout),normalized_spectrum=[float(eig[0]),float(eig[-1])]))
                    operator_checks[cid]=independent
            else:
                h['gradient'].program_admission(out,'analyze',h['native']);ck(float(energies[cid])==rec['value']['energy'],cid+' raw directional energy')
            rows[cid]={'output':pin(out),'energy_token':str(energies[cid]),'induction':ind}
    ck(fd==636 and fd_rows==10 and identities==21 and exports_n==12,'complete FD/identity/export denominator')
    directional={};direction_max=0.
    for g in q.NEAR:
        for kind in ('radial','rotation'):
            _,tangent=h['old_harness'].path_geometry(m['rows'][g]['positions_A'],kind,0.);grad=r['rows']['combined__'+g]['value']['gradient'];analytic=math.fsum(x*t for gr,tr in zip(grad,tangent) for x,t in zip(gr,tr));saved=r['directional_checks'][g+'__'+kind];steps=[]
            for index,step in enumerate((1e-4,5e-5)):
                p=energies[f'directional__{g}__{kind}__{step:g}__+1'];n=energies[f'directional__{g}__{kind}__{step:g}__-1'];slope=float((p-n)/(2*Decimal(str(step))));err=abs(slope-analytic);bound=saved['steps'][index]['checks']['slope']['bound'];ck(err<=bound,g+' '+kind+' independent Decimal directional slope');direction_max=max(direction_max,err);steps.append(slope)
            ck(abs(steps[0]-steps[1])<=saved['checks']['two_step']['bound'],g+' '+kind+' independent two-step derivative');directional[g+'__'+kind]={'Decimal_slopes':steps,'analytic':analytic}
    fields={};odd_max=0.;even_max=0.
    for g in q.DIMERS:
        replay=q.field_checks(g,r['rows'],m['rows'][g],m['units'],h);derived(replay,r['field_checks'][g],g+'/field',3*m['rows'][g]['atoms'])
        ex=r['rows']['matrix__'+g]['matrix_qualification']['exports'][1];A=np.array(ex['A']);ed,ep,u,v=np.array(ex['vectors']).T;z=np.tile([0.,0.,1.],m['rows'][g]['atoms']);L=np.linalg.cholesky(A);solve=lambda b:np.linalg.solve(L.T,np.linalg.solve(L,b));units=m['units'];K=units['electric'];B=units['bohr_A'];Ha=units['hartree_kcal_mol'];debye=units['debye_per_eA']
        mu0=m['rows'][g]['permanent_lab_Debye'][2]/(debye*B)+float(z@(solve(ed)+solve(ep)))/(2*B);chi=Ha/(K*B*B)*float(z@solve(z));p=energies[f'field__{g}__+1'];n=energies[f'field__{g}__-1'];zero=energies['matrix__'+g];f=Decimal('.001');odd=float(-(p-n)/(2*f*Decimal(str(Ha))));even=float(-(p+n-2*zero)/(f*f*Decimal(str(Ha))))
        ck(abs(odd-mu0)<=replay['checks']['odd_energy_conjugate_moment']['bound'] and abs(even-chi)<=replay['checks']['even_energy_matrix_curvature']['bound'],g+' independent Decimal field/Cholesky conjugacy');odd_max=max(odd_max,abs(odd-mu0));even_max=max(even_max,abs(even-chi));fields[g]=dict(odd_energy_au=odd,conjugate_moment_au=mu0,energy_curvature_a0_cubed=even,matrix_curvature_a0_cubed=chi)
    review=dict(status='passed_local_response_native_qualification',implementation=pin(__file__),result=pin(path),runner_implementation=r['implementation'],configuration=binding['configuration'],manifest=cfg['manifest'],source_manifest=cfg['source_manifest'],build_receipt=cfg['build_receipt'],executables={k:cfg[k] for k in ('analyze','testgrad')},native_source_review=binding['source_review'],login_source_review=binding['login_review'],review_bindings=pin(config_path),zero_baseline_parity_passed=True,force_checks_passed=True,moment_response_checks_passed=True,export_checks_passed=True,complete_programs=57,new_programs=57,inherited_programs=0,internal_FD_energies=fd,exact_EF_identities=identities,induction_summaries=induction_n,maximum_induction_residual_Debye=induction_max,maximum_Cartesian_error_decimal=str(fd_max),maximum_directional_error=direction_max,maximum_energy_route_error_kcal_mol=route_max,maximum_field_odd_error_au=odd_max,maximum_field_curvature_error_a0_cubed=even_max,normalized_spectrum_range=[min(x[0] for x in spectrum),max(x[1] for x in spectrum)],rows=rows,operator_checks=operator_checks,directional_checks=directional,field_checks=fields,checks=checks,new_native_calls=0,physical_accuracy_qualified=False,limitations='Numerical qualification of the declared vacuum single-metal water route and fixed profiles only; no fitting, physical-accuracy or production/GPU/periodic claim.',cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    review['derived_arithmetic_replay']={'raw_exports_exact':True,'unchanged_scientific_bounds':True,'field_bound_replay_propagation':'Four-energy arithmetic combination divided by Hartree*h^2 at the declared h=.001; used only for comparing recomputed bound metadata. Every original numerical gate is separately rerun unchanged.','differences':derived_differences,'maximum_absolute_difference':max((v['error'] for v in derived_differences),default=0.)}
    review['producer_cost']=r['cost']
    review['moment_replay_scope']='All raw moment fields exact; only COM-derived total_lab_Debye uses the existing matrix arithmetic bound, declared before actual result inspection.'
    (dest/'REVIEW.json').write_text(json.dumps(review,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:review[k] for k in ('status','complete_programs','induction_summaries','maximum_induction_residual_Debye','maximum_Cartesian_error_decimal','maximum_directional_error','maximum_field_curvature_error_a0_cubed')}));print(json.dumps({'review':pin(dest/'REVIEW.json'),'checks':len(checks)}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True,type=Path);parser.add_argument('--out',required=True,type=Path);args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    try:main(args.config.resolve(),args.out.resolve())
    except Exception as exc:
        (args.out/'FAILED_REVIEW.json').write_text(json.dumps({'status':'failed_saved_numerical_review','implementation':pin(__file__),'bindings':pin(args.config),'error':repr(exc),'new_native_calls':0},indent=2)+'\n')
        raise
