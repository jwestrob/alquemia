#!/usr/bin/env python3
"""Saved raw-output comparison only; no native invocation or numerical libraries."""
import argparse
from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import time

HERE = Path(__file__).resolve().parent
EPS = sys.float_info.epsilon
NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'

def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def read(p):
    assert pin(p['path'])==p, 'Changed pin '+p['path']
    return json.loads(Path(p['path']).read_text())

def raw(p):
    assert pin(p['path'])==p, 'Changed raw pin '+p['path']
    return Path(p['path']).read_text()

def dec(s): return Decimal(str(s).replace('D','E'))
def half(s): return Decimal('.5')*Decimal(10)**dec(s).as_tuple().exponent
def ar(dim,*x): return 128*dim*EPS*max([1.]+[abs(float(v)) for v in x])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--result',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    start,cpu=time.monotonic(),time.process_time();inputs=json.loads((HERE/'INPUTS.json').read_text());checks=[];deltas=[];mismatches=[];rows={};report={}
    def ck(ok,label):
        checks.append(label)
        if not ok:mismatches.append({'check':label})
        return bool(ok)
    def compare(label,x,y,dim,printing=False,bound=None):
        error=abs(dec(x)-dec(y));limit=Decimal(str(ar(dim,x,y))) if bound is None else Decimal(str(bound))
        if printing:limit+=half(x)+half(y)
        entry={'field':label,'old':str(x),'new':str(y),'difference':float(error),'bound':float(limit),'passed':error<=limit}
        checks.append(label)
        if error:deltas.append(entry)
        if error>limit:mismatches.append(entry)
        return float(error)
    try:
        assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')),'One-thread context'
        old=read(inputs['old_result']);oldfirst=read(inputs['old_first_result']);oldreview=read(inputs['old_review']);build=read(inputs['new_build']);new=json.loads(a.result.read_text())
        assert a.result.resolve()==Path(inputs['expected_new_result']), 'Exact future result path'
        assert oldreview['status']=='passed_local_response_native_qualification' and oldreview['result']==inputs['old_result'],'Admitted original 57'
        om=read(old['manifest']);nm=read(new['manifest']);cfg=new['config'];assert cfg['build_receipt']==inputs['new_build'] and build['status']=='complete','Build join'
        assert new['implementation']==inputs['new_implementation'] and cfg['source_manifest']==build['manifest'],'Source join'
        assert all(cfg[k]==build['executables'][k] and pin(cfg[k]['path'])==cfg[k] for k in ('analyze','testgrad')),'Executable pins'
        ck(new['status']=='local_response_numerical_qualification_complete' and new['error'] is None and new['global_stop'] is False and new['numerical_checks_passed'] is True,'new terminal complete numerical panel')
        ck(new['native_program_invocations']==57 and new['internal_FD_energies_completed']==636,'fresh 57 programs / 636 FD')
        jobs={j['id']:j for j in om['programs']};njob={j['id']:j for j in nm['programs']}
        assert set(jobs)==set(njob)==set(old['rows'])==set(new['rows']) and len(jobs)==57,'Exact full denominator'
        assert om['units']==nm['units'],'Exact units'
        sources={'run_v1':0,'run_v2':0};counts={'gradients':0,'FD_programs':0,'energies':0,'exports':0,'active_programs':0,'printed_exact_energy_programs':0,'printed_exact_analytic_programs':0};maxima={}
        def category(name,value):maxima[name]=max(maxima.get(name,0.),value)
        for cid,job in jobs.items():
            before,after=old['rows'][cid],new['rows'][cid];dim=3*om['rows'][job['cell']]['atoms']
            assert {k:v for k,v in job.items() if k!='input'}=={k:v for k,v in njob[cid].items() if k!='input'},'Job metadata '+cid
            if before==oldfirst['rows'].get(cid):sources['run_v1']+=1
            else:sources['run_v2']+=1
            if not ck(after['status']=='complete',cid+' complete'):continue
            for k in ('xyz','key','parameters'):
                assert raw(before['runtime_input'][k])==raw(after['runtime_input'][k])==raw(job['input'][k])==raw(njob[cid]['input'][k]),'Literal input '+cid+' '+k
            bo,no=raw(before['actual_output']),raw(after['actual_output'])
            rows[cid]={'old_output':before['actual_output'],'new_output':after['actual_output'],'profile':job['profile'],'program':job['program']}
            assert before['receipts'][0]['argv'][1:]==after['receipts'][0]['argv'][1:],'Same arguments '+cid
            assert after['receipts'][0]['returncode']==0 and after['receipts'][0]['argv'][0]==cfg[job['program']]['path'],'New receipt '+cid
            be=re.findall(r'Total Potential Energy\s*:\s*('+NUMBER+r')\s+Kcal/mole',bo);ne=re.findall(r'Total Potential Energy\s*:\s*('+NUMBER+r')\s+Kcal/mole',no)
            assert len(be)==len(ne)==1,'Energy denominator '+cid
            counts['energies']+=1;counts['printed_exact_energy_programs']+=dec(be[0])==dec(ne[0]);category('energy_kcal_mol',compare(cid+'/energy',be[0],ne[0],dim,True))
            for tag in ('PARAMETERS','PAIR','ENVIRONMENT','SPECTRUM','INDUCED','DUAL'):
                left=re.findall(r'^\s*LOCAL_RESPONSE_'+tag+r'\s+(.+)$',bo,re.M);right=re.findall(r'^\s*LOCAL_RESPONSE_'+tag+r'\s+(.+)$',no,re.M)
                assert len(left)==len(right),'Response row denominator '+cid+'/'+tag
                for i,(l,r) in enumerate(zip(left,right)):
                    l,r=l.split(),r.split();assert len(l)==len(r),'Response width'
                    for j,(x,y) in enumerate(zip(l,r)):
                        if tag=='PARAMETERS':ck(dec(x)==dec(y),cid+'/exact PARAMETERS')
                        else:category('readback_'+tag,compare(f'{cid}/{tag}/{i}/{j}',x,y,dim))
            counts['active_programs']+=bool(re.search(r'^\s*LOCAL_RESPONSE_',bo,re.M))
            if job['program']=='testgrad':
                def grad(text):
                    section=text.rsplit('Cartesian Gradient Breakdown over Individual Atoms :',1)[-1]
                    return {(kind,int(site)):xyz for kind,site,*xyz in re.findall(r'^\s*(Anlyt|Numer)\s+(\d+)\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+'+NUMBER+r'\s*$',section,re.M)}
                bg,ng=grad(bo),grad(no);assert bg.keys()==ng.keys() and len(bg)==(dim//3)*(2 if job['fd'] else 1),'Gradient denominator '+cid
                counts['gradients']+=1;exact=True
                for key in bg:
                    for axis,(x,y) in enumerate(zip(bg[key],ng[key])):
                        label=f'{cid}/{key}/{axis}'
                        if key[0]=='Anlyt':
                            category('analytic_gradient_kcal_mol_A',compare(label,x,y,dim,True));exact &= dec(x)==dec(y)
                        else:
                            ax=bg[('Anlyt',key[1])][axis];ay=ng[('Anlyt',key[1])][axis]
                            ck(abs(dec(x)-dec(ax))<=Decimal('.001') and abs(dec(y)-dec(ay))<=Decimal('.001'),label+'/original FD gates')
                            category('numerical_gradient_kcal_mol_A',compare(label,x,y,dim,bound=.002+float(half(ax)+half(ay))+ar(dim,ax,ay)))
                counts['printed_exact_analytic_programs']+=exact;counts['FD_programs']+=job['fd']
            else:
                def comp(text):
                    return {name.strip():(v,int(n)) for name,v,n in re.findall(r'^\s*([A-Za-z][A-Za-z /-]+?)\s+('+NUMBER+r')\s+(\d+)\s*$',text,re.M)}
                bc,nc=comp(bo),comp(no);assert bc.keys()==nc.keys(),'Component labels '+cid
                for name,(x,count) in bc.items():
                    y,n=nc[name];ck(n==count,cid+'/component count '+name);category('component_energy_kcal_mol',compare(cid+'/component '+name,x,y,dim,True))
                for tag in ('Dipole Moment Magnitude','Dipole X,Y,Z-Components'):
                    l=re.findall(re.escape(tag)+r'\s*:\s*([^\n]+)',bo);r=re.findall(re.escape(tag)+r'\s*:\s*([^\n]+)',no);assert len(l)==len(r),'Moment denominator'
                    for xline,yline in zip(l,r):
                        x,y=re.findall(NUMBER,xline),re.findall(NUMBER,yline);assert len(x)==len(y)
                        for axis,(xx,yy) in enumerate(zip(x,y)):category('printed_moment_Debye',compare(cid+'/'+tag+str(axis),xx,yy,dim,True))
                if 'moments' in before:
                    ub,un=raw(before['moments']['induced_output']).splitlines(),raw(after['moments']['induced_output']).splitlines();assert len(ub)==len(un),'Saved uind rows'
                    for i,(x,y) in enumerate(zip(ub[1:],un[1:])):
                        xx,yy=x.split(),y.split();assert xx[:2]==yy[:2] and len(xx)==len(yy),'uind identity'
                        for j,(l,r) in enumerate(zip(xx[2:],yy[2:])):category('saved_induced_Debye',compare(f'{cid}/uind/{i}/{j}',l,r,dim,True))
            assert len(before.get('exports',[]))==len(after.get('exports',[])),'Export denominator'
            for seq,(bp,np) in enumerate(zip(before.get('exports',[]),after.get('exports',[]))):
                bl,nl=raw(bp).splitlines(),raw(np).splitlines();assert len(bl)==len(nl),'Export lines'
                for i,(l,r) in enumerate(zip(bl,nl)):
                    x,y=l.split(),r.split();assert len(x)==len(y) and x[0]==y[0],'Export grammar'
                    tag=x[0]
                    if tag in ('A','VECTOR','CHECKS'):
                        offset=2 if tag in ('A','VECTOR') else 1;assert x[:offset]==y[:offset]
                        for j,(xx,yy) in enumerate(zip(x[offset:],y[offset:])):category('export_'+tag,compare(f'{cid}/export{seq}/{i}/{j}',xx,yy,dim))
                    else:ck(x==y,f'{cid}/export{seq}/exact '+tag)
                counts['exports']+=1
        assert sources=={'run_v1':27,'run_v2':30},'Original two-attempt map'
        ck(counts['energies']==57 and counts['gradients']==31 and counts['FD_programs']==10 and counts['exports']==12,'Complete raw comparison scope')
        for g,b in old['field_checks'].items():
            n=new['field_checks'][g];dim=3*om['rows'][g]['atoms'];Ha=om['units']['hartree_kcal_mol'];f=.001
            ck(b['passed'] and n['passed'] and all(v['passed'] for v in b['checks'].values()) and all(v['passed'] for v in n['checks'].values()),g+'/unchanged field gates')
            for key,prop in [('odd_energy_mu_au',1e-8/(f*Ha)),('even_energy_chi_a0_cubed',4e-8/(f*f*Ha)),('mu0_lab_au',0),('matrix_chi_a0_cubed',0)]:
                category(key,compare(g+'/'+key,b[key],n[key],dim,bound=prop+ar(dim,b[key],n[key])))
        for g,b in old['directional_checks'].items():
            n=new['directional_checks'][g];ck(b['passed'] and n['passed'],g+'/unchanged directional gates')
            dim=3*om['rows'][g.rsplit('__',1)[0]]['atoms']
            for i,(x,y) in enumerate(zip(b['steps'],n['steps'])):
                assert x['half_step_A']==y['half_step_A']
                category('directional_slope',compare(f'{g}/{i}/FD',x['FD'],y['FD'],dim,bound=1e-8/x['half_step_A']+ar(dim,x['FD'],y['FD'])))
                category('directional_analytic',compare(f'{g}/{i}/analytic',x['analytic'],y['analytic'],dim,bound=max(0.,x['printing_bound']-5e-9/x['half_step_A'])+max(0.,y['printing_bound']-5e-9/y['half_step_A'])+ar(dim,x['analytic'],y['analytic'])))
        report.update(original_source_counts=sources,counts=counts,maxima=maxima,new_cost=new['cost'],old_cost_recovery=old['cost'],old_cost_first_attempt=oldfirst['cost'])
    except Exception as exc:mismatches.append({'exception':repr(exc)})
    report.update(status='passed_saved_cpu_functional_parity' if not mismatches else 'failed_or_incomplete_saved_cpu_functional_parity',implementation=pin(__file__),specification=pin(HERE/'SPEC.md'),inputs=pin(HERE/'INPUTS.json'),result=pin(a.result),checked_programs=len(rows),rows=rows,check_count=len(checks),nonzero_differences=deltas,mismatches=mismatches,authorship='Reviewer authored the refactor; independent actual numerical admission is separate.',new_native_calls=0,physical_accuracy_qualified=False,limitations='Fixed vacuum 9-fixture / 57-program profiles only. No periodic, MD, GPU, global domain or physical-accuracy inference.',cost={'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'threads':1})
    (a.out/'REVIEW.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report.get(k) for k in ('status','check_count','counts','maxima','cost')}));print(json.dumps({'mismatch_count':len(mismatches),'first_mismatches':mismatches[:8]}));print(json.dumps(pin(a.out/'REVIEW.json')))
    return 0 if not mismatches else 1

if __name__=='__main__':raise SystemExit(main())
