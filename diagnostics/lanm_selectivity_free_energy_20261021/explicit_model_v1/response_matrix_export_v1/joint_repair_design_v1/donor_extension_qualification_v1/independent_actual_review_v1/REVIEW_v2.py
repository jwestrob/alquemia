#!/usr/bin/env python3
"""Saved 121-slot review: disclosed producer replay plus raw independent arithmetic."""
import argparse
import collections
import csv
from decimal import Decimal, getcontext
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import time
import traceback
import numpy as np
getcontext().prec=45
N=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
DEC=lambda s:Decimal(str(s).replace('D','E'))
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def require(v,s):
    if not v:raise ValueError(s)
def verify(p):
    require(pin(p['path'])==p,'Changed pin '+p['path']);return Path(p['path'])
def raw_records(text,label):return [[DEC(t) for t in line.split()[1:]] for line in text.splitlines() if line.strip().startswith(label+' ')]
def arr(v):return np.asarray(v,dtype=float)
def maxabs(v):return float(np.max(np.abs(v)))
def check_error(error,bound):return {'error':float(error),'bound':float(bound),'passed':math.isfinite(float(error)) and abs(error)<=bound}
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')

def saved_replay_v2(author,cfg,r):
    """Producer replay, one versioned equality delta for derived lab dipoles."""
    m,h=author.context(cfg,False);complete={};calls=0;fd=0;deltas=[]
    require(r['config']==cfg and r['manifest']==cfg['manifest'] and r['implementation']==cfg['implementation'],'Saved result lineage')
    require(set(r['rows'])=={j['id'] for j in m['programs']},'Exact saved row denominator')
    for job in m['programs']:
        rec=r['rows'][job['id']]
        if rec['status']=='not_attempted':continue
        calls+=len(rec['receipts']);require(len(rec['receipts'])==1,'Single call per slot')
        for k,p in rec['runtime_input'].items():verify(p);require(p['sha256']==job['input'][k]['sha256'],'Actual input bytes')
        verify(rec['start']);receipt=rec['receipts'][0];verify(receipt['output']);require(receipt['returncode']==0,'Saved native completion')
        executable=m['sources']['original_testgrad'] if job['binary_role']=='original_testgrad' else cfg[job['program']]
        require(receipt['argv']==[executable['path'],*job['argv_tail']],'Actual binary/argv')
        parsed=author.parse_cell(job,m,receipt['output']['path'],h)
        for key in ('output','value','moments','exports'):
            if key not in parsed:continue
            if key!='moments':require(parsed[key]==rec[key],'Raw saved parse differs: '+key);continue
            fresh,stored=parsed[key],rec[key]
            require(fresh.keys()==stored.keys(),'Moment key denominator')
            require({k:v for k,v in fresh.items() if k!='total_lab_Debye'}=={k:v for k,v in stored.items() if k!='total_lab_Debye'},'Raw moment fields differ')
            difference=maxabs(arr(fresh['total_lab_Debye'])-stored['total_lab_Debye']);bound=h['matrix'].arithmetic(3*m['rows'][job['cell']]['atoms'],fresh['total_lab_Debye'],stored['total_lab_Debye'])
            require(difference<=bound,'Derived lab dipole exceeds inherited arithmetic allowance')
            if difference:deltas.append({'row':job['id'],'field':'moments.total_lab_Debye','fresh':fresh['total_lab_Debye'],'stored':stored['total_lab_Debye'],'maximum_difference_Debye':difference,'inherited_arithmetic_bound_Debye':bound})
        complete[job['id']]=dict(parsed,job=job);fd+=job['internal_FD_energies']
    require(calls==r['native_program_invocations'] and fd==r['internal_FD_energies_completed'],'Actual counters')
    require(calls==121 and len(complete)==121 and fd==576,'Full qualification requires121/576')
    reports={'baseline':author.baseline_checks(m,complete),'operator':author.operator_checks(m,complete,h),'field':author.field_checks(m,complete,h),'mechanics':author.mechanics_checks(m,complete,h)}
    ok=all(v['numerical_checks_passed'] for v in complete.values()) and all(v['passed'] for v in reports.values())
    require(ok==r['numerical_checks_passed'],'Saved numerical status differs')
    if ok:require(r['status']=='donor_extension_numerical_qualification_complete' and r['error'] is None and r['global_stop'] is False,'Successful terminal status/error/stop')
    else:require(r['status']=='donor_extension_numerical_qualification_failed' and r['global_stop'] is True,'Failed qualification preserved')
    return {'passed':ok,'reports':reports,'new_native_calls':0,'complete':121,'internal_FD_energies':576,'derived_moment_replay_differences':deltas,'only_equality_delta':'total_lab_Debye derived from exact raw COM moment, charge and literal coordinates; existing matrix arithmetic allowance only'}

def review(a):
    started,cpu=time.monotonic(),time.process_time();out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
    cfg=json.loads(Path(a.config).read_text());r=json.loads(Path(a.result).read_text())
    base={'implementation':pin(__file__),'config':pin(a.config),'result':pin(a.result),'new_native_calls':0,'new_quantum_calls':0,'author_suite_reruns':0}
    try:
        require(r['status'] in ('donor_extension_numerical_qualification_complete','donor_extension_numerical_qualification_failed'),'Terminal result required')
        require(cfg['implementation']['sha256']=='bce1a96c7f6374429508b412a3302d63abc2a5b6a76380d05a9a01e3bd7c6c79','Frozen runner')
        require(cfg['matrix_adapter']['sha256']=='1ee3364f31074dddcb5c684355cbc08dd4dfeda55c43e338759d72aac443168e','Frozen dual matrix adapter')
        require(cfg['manifest']['sha256']=='12bc8b408c01eafd34284ae91137ae09eaea60759b2fd5dbd3f5f285b7781b54','Frozen manifest')
        for k in ('implementation','matrix_adapter','manifest','source_manifest','build_receipt','source_review','runner_review','root_release','analyze','testgrad'):verify(cfg[k])
        m=json.loads(verify(cfg['manifest']).read_text());jobs={j['id']:j for j in m['programs']}
        require(r['config']==cfg and r['implementation']==cfg['implementation'] and r['manifest']==cfg['manifest'],'Actual producer lineage')
        require(len(jobs)==121 and set(r['rows'])==set(jobs),'Full 121-slot denominator')
        spec=importlib.util.spec_from_file_location('disclosed_author_runner',verify(cfg['implementation']));author=importlib.util.module_from_spec(spec);spec.loader.exec_module(author)
        author.context(cfg,False)
        reports={};data={};raw_summary={};failures={};calls=0;starts=0;fd_requested=0;fd_observed=0;all_spectra=[];induction=[];fd_errors=[];virial_errors=[];export_spectra=[]
        Ha=DEC(m['units']['hartree_kcal_mol']);B=DEC(m['units']['bohr_A']);t=[DEC(x) for x in m['field_direction']]
        for cid,job in jobs.items():
            rec=r['rows'][cid];stat=rec['status']
            if stat!='complete':failures[cid]={'status':stat,'error':rec.get('error')}
            if stat=='not_attempted':require(not rec.get('receipts') and not rec.get('start'),'Untouched row has native receipt');continue
            require(rec.get('job')==job,'Actual job metadata '+cid)
            if rec.get('start'):
                start=json.loads(verify(rec['start']).read_text());starts+=1;fd_requested+=job['internal_FD_energies']
                exe=m['sources']['original_testgrad'] if job['binary_role']=='original_testgrad' else cfg[job['program']]
                verify(exe);argv=[exe['path'],*job['argv_tail']]
                require(start['job']==job and start['argv']==argv and str(start['job_id'])==str(r['job_id']) and start['input']==rec['runtime_input'],'Raw START binding '+cid)
                for k,p in rec['runtime_input'].items():verify(p);require(p['sha256']==job['input'][k]['sha256'],'Runtime bytes '+cid)
            receipts=rec.get('receipts',[]);require(len(receipts)<=1,'No native repeats');calls+=len(receipts)
            if not receipts:continue
            receipt=receipts[0];op=verify(receipt['output']);require(receipt['argv']==argv,'Receipt argv '+cid)
            if receipt['returncode']!=0:continue
            text=op.read_text();energies=re.findall(r'Total Potential Energy\s*:\s*('+N+')',text)
            if not energies:continue
            require(len(energies)==1,'Raw energy denominator '+cid);E=DEC(energies[0])
            if 'value' in rec:require(float(E)==rec['value']['energy'],'Raw energy '+cid)
            n=m['rows'][job['cell']]['atoms'];summaries=re.findall(r'Induced Dipoles\s*:\s*Iterations\s+(\d+)\s+RMS Residual\s+(\S+)',text)
            if job['binary_role']!='original_testgrad':require(summaries,'Missing induction readback '+cid)
            for it,q in summaries:
                q=DEC(q);require(int(it)>=1 and q.is_finite() and 0<=q<=DEC('1e-8'),'Raw induction gate '+cid);induction.append(float(q))
            spectra=raw_records(text,'LOCAL_RESPONSE_SPECTRUM')
            for s in spectra:
                require(len(s)==5 and s[1]>0 and s[0]>DEC('1e-6')*s[1] and 0<=s[2]<=DEC('1e-10') and 0<=s[3]<=DEC('1e-10')*max(Decimal(1),s[1]) and 0<=s[4]<=DEC('1e-10'),'Raw spectrum gate '+cid)
                all_spectra.append([float(x) for x in s])
            dd={'energy':E};numeric={}
            if job['program']=='testgrad':
                aa=re.findall(r'^\s*Anlyt\s+(\d+)\s+('+N+r')\s+('+N+r')\s+('+N+r')\s+'+N+r'\s*$',text,re.M)
                nn=re.findall(r'^\s*Numer\s+(\d+)\s+('+N+r')\s+('+N+r')\s+('+N+r')\s+'+N+r'\s*$',text,re.M)
                require([int(v[0]) for v in aa]==list(range(1,n+1)),'Analytic atom denominator '+cid)
                ag=[[DEC(v) for v in row[1:]] for row in aa];dd['gradient']=ag
                if job['internal_FD_energies']:
                    require([int(v[0]) for v in nn]==list(range(1,n+1)),'Numerical atom denominator '+cid)
                    ng=[[DEC(v) for v in row[1:]] for row in nn];error=max(abs(x-y) for xs,ys in zip(ag,ng) for x,y in zip(xs,ys));numeric['Cartesian_gradient']=check_error(error,Decimal('.001'));fd_errors.append(float(error));fd_observed+=6*n
                full=raw_records(text,'DONOR_RESPONSE_TOTAL_GRADIENT');vir=raw_records(text,'DONOR_RESPONSE_TOTAL_VIRIAL')
                if full:
                    require([int(v[0]) for v in full]==list(range(1,n+1)) and [int(v[0]) for v in vir]==[1,2,3],'Full gradient/virial denominator '+cid)
                    dd['gradient']=[x[1:] for x in full];dd['virial']=[x[1:] for x in vir]
                    xyz=[[DEC(x) for x in line.split()[2:5]] for line in verify(rec['runtime_input']['xyz']).read_text().splitlines()[1:]]
                    predicted=[[sum((xyz[k][i]*dd['gradient'][k][j]+xyz[k][j]*dd['gradient'][k][i])/2 for k in range(n)) for j in range(3)] for i in range(3)]
                    err=max(abs(x-y) for xs,ys in zip(predicted,dd['virial']) for x,y in zip(xs,ys));bound=DEC(rec['checks']['virial_coordinate_gradient']['bound']);numeric['virial_coordinate_gradient']=check_error(err,bound);virial_errors.append(float(err))
            central=text.split('DONOR_RESPONSE_TOTAL_ENERGY',1)[0] if job['program']=='testgrad' else text
            uv=raw_records(central,'LOCAL_RESPONSE_INDUCED');perm=raw_records(central,'DONOR_RESPONSE_PERMANENT_LAB_EA')
            if uv and perm and job['program']=='analyze' and job['argv_tail'][1]=='EM':
                uv=uv[-n:];require([int(v[0]) for v in uv]==list(range(1,n+1)),'Central dual ordering '+cid)
                mu=[(perm[-1][k]+sum((v[1+k]+v[4+k])/2 for v in uv))/B for k in range(3)]
                dd['projected_mu']=sum(x*y for x,y in zip(mu,t));require(maxabs(arr(mu)-rec['mu_lab_au'])<2e-13,'Independent charged-origin dual moment '+cid)
            for ep in rec.get('exports',[]):
                lines=verify(ep).read_text().splitlines();A=arr([[float(x) for x in line.split()[2:]] for line in lines if line.strip().startswith('A ')])
                require(A.shape==(3*n,3*n),'Raw matrix dimension '+cid);np.linalg.cholesky(A);ev=np.linalg.eigvalsh(A);export_spectra.append({'path':ep,'min_eigenvalue':float(ev[0]),'max_eigenvalue':float(ev[-1]),'symmetry_max':maxabs(A-A.T)})
            data[cid]=dd;raw_summary[cid]={'status':stat,'output':receipt['output'],'raw_checks':numeric,'induction_summaries':len(summaries),'spectrum_records':len(spectra)}
        require(calls==r['native_program_invocations'] and starts==calls and fd_requested==r['internal_FD_energies_requested'],'Actual invocation/START counters')
        require(calls<=121 and fd_requested<=576 and r['requested_slots']<=121,'Fixed finite caps')
        full=(len(data)==121 and all(len(x.get('receipts',[]))==1 and x['receipts'][0]['returncode']==0 for x in r['rows'].values()))
        replay={'authorship':'Frozen producer parsing and report functions through reviewer saved_replay_v2; only derived total_lab_Debye equality uses existing arithmetic allowance. Original failed exact replay preserved.','executed':False}
        if full:
            replay.update(executed=True,result=saved_replay_v2(author,cfg,r))
            require(fd_observed==r['internal_FD_energies_completed']==576 and len(export_spectra)==20,'Full raw FD/export counters')
            reports=replay['result']['reports']
        def find(cell,profile,stage,program,field=None,q=None,fd=None,transform=None):
            found=[j for j in jobs.values() if j['cell']==cell and j['profile']==profile and j['stage']==stage and j['program']==program and (field is None or j['field_au']==field) and (fd is None or j['fd_total_separation_A']==fd) and (q is None or j['transform'] and j['transform'].get('q_A')==q) and (transform is None or j['transform'] and j['transform']['kind']==transform)]
            require(len(found)==1,'Independent unique lookup');return data[found[0]['id']]
        field={};mechanics={}
        if full:
            for metal in ('La','Dy'):
                cell=metal+'_central'
                for profile in ('drive','stiffness','combined'):
                    zero=find(cell,profile,'C','analyze');tag=metal+'_'+profile;field[tag]={}
                    for label,h in [('h',DEC('.001')),('half',DEC('.0005'))]:
                        plus=find(cell,profile,'D','analyze',field=float(h));minus=find(cell,profile,'D','analyze',field=-float(h))
                        z={'chi_energy':-(plus['energy']+minus['energy']-2*zero['energy'])/(Ha*h*h),'chi_moment':(plus['projected_mu']-minus['projected_mu'])/(2*h),'odd_energy_moment':-(plus['energy']-minus['energy'])/(2*h*Ha),'zero_moment':zero['projected_mu']}
                        reference=reports['field']['states'][tag]['estimates'][label]
                        require(all(abs(float(v)-reference[k])<1e-7 for k,v in z.items()),'Independent Decimal field arithmetic '+tag)
                        field[tag][label]={k:float(v) for k,v in z.items()}
                center=find(cell,'combined','C','testgrad',fd=.0001);fr=m['rows'][cell]['rotation_frame'];xyz=arr(m['rows'][cell]['positions_A']);tangent=np.zeros_like(xyz);tangent[1:]=np.cross(arr(fr['rotation_unit']),xyz[1:]-fr['centroid_A'])/fr['rotation_length_A'];tan=[[DEC(x) for x in v] for v in tangent]
                project=lambda grad:sum(x*y for xs,ys in zip(grad,tan) for x,y in zip(xs,ys))
                analytic=project(center['gradient']);z={'rotation_analytic':float(analytic),'rotation':{},'Maxwell':{},'strain':{}}
                for q in (DEC('.01'),DEC('.005')):
                    endpoints={}
                    for sign in (-1,1):
                        ec=next(k for k,v in m['rows'].items() if v['kind']=='donor_rotation' and v['metal']==metal and v['q_A']==float(sign*q));endpoints[sign]=find(ec,'combined','E','analyze')
                    slope=(endpoints[1]['energy']-endpoints[-1]['energy'])/(2*q);dmu=(endpoints[1]['projected_mu']-endpoints[-1]['projected_mu'])/(2*q);z['rotation'][str(q)]={'energy_slope':float(slope),'moment_slope':float(dmu),'energy_check':check_error(abs(slope-analytic),DEC('.1')+DEC('5e-9')/q)}
                    for h in (DEC('.001'),DEC('.0005')):
                        gp=find(cell,'combined','D','testgrad',field=float(h));gm=find(cell,'combined','D','testgrad',field=-float(h));mixed=-(project(gp['gradient'])-project(gm['gradient']))/(2*h*Ha);key=f'{float(q)}_{float(h)}';author_mixed=reports['mechanics']['mixed'][metal][key]
                        require(abs(float(mixed)-author_mixed['force_field_derivative'])<1e-10 and abs(float(dmu)-author_mixed['moment_geometry_derivative'])<1e-10,'Independent Maxwell arithmetic')
                        z['Maxwell'][key]={'force_field_derivative':float(mixed),'moment_geometry_derivative':float(dmu),'check':check_error(abs(mixed-dmu),DEC(author_mixed['bound']))}
                S=[[DEC(x) for x in row] for row in [[1,.2,-.3],[.2,-.7,.4],[-.3,.4,.5]]]
                for label,h in [('h',DEC('.0001')),('half',DEC('.00005'))]:
                    jp=next(j for j in jobs.values() if j['cell']==cell and j['stage']=='G' and j['transform']['q_A']==float(h));plus=data[jp['id']];minus=find(cell,'combined','G','analyze',q=-float(h));pred=sum(x*y for xs,ys in zip(S,center['virial']) for x,y in zip(xs,ys))/DEC(jp['transform']['normalization_L_A']);obs=(plus['energy']-minus['energy'])/(2*h)
                    z['strain'][label]={'analytic':float(pred),'energy_slope':float(obs),'check':check_error(abs(obs-pred),DEC('.001')+DEC('5e-9')/h)}
                mechanics[metal]=z
        independent_gates=[v for row in raw_summary.values() for v in row['raw_checks'].values()]
        for row in mechanics.values():
            independent_gates += [v['energy_check'] for v in row['rotation'].values()]
            independent_gates += [v['check'] for v in row['Maxwell'].values()]
            independent_gates += [v['check'] for v in row['strain'].values()]
        independent_passed=all(v['passed'] for v in independent_gates)
        if full and replay['result']['passed']:require(independent_passed,'Independent raw numerical gate failed')
        base['independent_raw_numerical_checks']={'passed':independent_passed,'count':len(independent_gates)}
        table=list(csv.DictReader(Path(a.sacct).read_text().splitlines(),delimiter='|'));parent=[x for x in table if x['JobIDRaw']==str(r['job_id'])];require(len(parent)==1,'Exact parent scheduler job')
        sched=parent[0];require(sched['State'].split()[0] not in ('RUNNING','PENDING','CONFIGURING','COMPLETING') and int(sched['AllocCPUS'])==12,'Terminal12CPU allocation')
        allocated=int(sched['ElapsedRaw'])*int(sched['AllocCPUS'])
        base.update(status='passed_donor_extension_numerical_actual_review' if full else 'reviewed_incomplete_donor_extension_numerical_qualification',producer_status=r['status'],qualification_passed=bool(full and replay['result']['passed']),source_manifest=cfg['source_manifest'],build_receipt=cfg['build_receipt'],executables={k:cfg[k] for k in ('analyze','testgrad')},producer_implementation=cfg['implementation'],manifest=cfg['manifest'],author_replay=replay,slot_status_counts=dict(collections.Counter(x['status'] for x in r['rows'].values())),failed_or_untouched=failures,raw_rows=raw_summary,independent_field_estimates=field,independent_mechanics=mechanics,raw_matrix_spectra=export_spectra,raw_maximum_Cartesian_FD_error=max(fd_errors,default=None),raw_maximum_virial_identity_error=max(virial_errors,default=None),raw_maximum_induction=max(induction,default=None),raw_induction_summary_count=len(induction),raw_spectrum_count=len(all_spectra),raw_minimum_normalized_eigen_ratio=min((s[0]/s[1] for s in all_spectra),default=None),counts={'native_receipts':calls,'starts':starts,'internal_FD_requested':fd_requested,'internal_FD_observed':fd_observed,'snapshots':len(export_spectra),'full_slots':121},scheduler={'receipt':pin(a.sacct),'parent':sched,'allocated_CPU_seconds':allocated,'step_rows_not_double_counted':True},producer_cost=r.get('cost'),review_cost={'wall_seconds':time.monotonic()-started,'CPU_seconds':time.process_time()-cpu,'threads':1},limitations=['Sampled numerical consistency only; no physical accuracy, fitting or target promotion.','Producer parsing/report logic is reused explicitly via saved_replay_v2; only derived lab-moment equality is bounded by inherited arithmetic. Decimal/raw arithmetic is independent.','Finite-q mixed envelope is not a truncation or physical uncertainty bound; force-solver sensitivity remains unbounded.','Failed and untouched slots are preserved; an incomplete result is not numerical qualification.'])
        save(out/'REVIEW.json',base);print(json.dumps(pin(out/'REVIEW.json')));print(json.dumps({k:base[k] for k in ('status','qualification_passed','slot_status_counts','counts','review_cost')}))
    except BaseException as exc:
        base.update(status='failed_independent_saved_review',error=repr(exc),traceback=traceback.format_exc(),review_cost={'wall_seconds':time.monotonic()-started,'CPU_seconds':time.process_time()-cpu,'threads':1});save(out/'FAILURE.json',base);raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--result',required=True);p.add_argument('--sacct',required=True);p.add_argument('--out',required=True);review(p.parse_args())
