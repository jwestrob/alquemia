#!/usr/bin/env python3
"""Independent saved-output arithmetic for the fixed 56-program qualification."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/metal_pair_response_v1'
NUMBER=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?'
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
PROFILES={'zero':(0.,0.),'direct':(.05,0.),'mutual':(0.,-.05),'combined':(.05,-.05)}


def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    started,cpu=time.monotonic(),time.process_time();checks=[]
    def check(name,ok):
        checks.append({'name':name,'passed':bool(ok)})
        if not ok:raise ValueError(name)
    def verify(p):
        check('pin '+p['path'],pin(p['path'])==p);return Path(p['path'])
    def read(p):return json.loads(verify(p).read_text())
    def number(s):return float(s.replace('D','E').replace('d','e'))
    report={'status':'review_started','implementation':pin(__file__),'checks':checks,'new_native_calls':0,'physical_accuracy_qualified':False,'error':None}
    try:
        check('one thread',all(os.environ.get(k)=='1' for k in THREADS))
        rp=pin(W/'qualification_run_v1/execution/RESULT.json');r=read(rp);c=read(r['config']);b=read(r['manifest'])
        report.update(result=rp,config=r['config'],manifest=r['manifest'])
        check('terminal original job and exact source',r['status']=='metal_pair_numerical_qualification_complete' and r['error'] is None and r['job_id']=='18383' and r['implementation']==c['implementation']==b['implementation'] and r['numerical_checks_passed'] is True)
        for k in ('implementation','prototype_testgrad','prototype_analyze','prototype_source_manifest','source_audit','build_result'):verify(c[k])
        build=read(c['build_result']);audit=read(c['source_audit'])
        check('source/build/audit authority',build['status']=='complete' and build['manifest']==c['prototype_source_manifest']==audit['source_manifest'] and build['executables']=={'testgrad':c['prototype_testgrad'],'analyze':c['prototype_analyze']} and all(audit['checks'].values()) and audit['status']=='root_source_audit_complete')
        check('exact fixed counters',r['native_program_invocations']==56 and r['internal_FD_energies_requested']==r['internal_FD_energies_completed']==384 and len(r['rows'])==56 and len(r['fixtures'])==9)
        for p in b['sources'].values():verify(p)
        for fixture in b['rows'].values():
            for p in fixture['input'].values():verify(p)
        check('exact56 slots',set(r['rows'])=={j['id'] for j in b['programs']})
        values={};identities={};fd_errors={};residuals=[];wall=[];user=[];system=[];readbacks=0
        for job in b['programs']:
            label=job['id'];row=r['rows'][label];fixture=b['rows'][job['cell']];profile=job['profile'];program=job['program']
            check(label+' completion',row['status']=='complete' and row['job']==job and len(row['receipts'])==1)
            receipt=row['receipts'][0];text=verify(receipt['output']).read_text()
            for p in row['input'].values():verify(p)
            exe=b['sources']['original_testgrad'] if profile=='original' else c['prototype_'+program]
            arguments=['system.xyz','E'] if program=='analyze' else ['system.xyz','Y','Y','0.00001'] if job['fd'] else ['system.xyz','Y','N']
            check(label+' actual execution',receipt['argv']==[exe['path']]+arguments and receipt['returncode']==0 and receipt['threads']=={k:'1' for k in THREADS} and Path(receipt['cwd'])==Path(row['input']['xyz']['path']).parent)
            check(label+' native completion markers','Tinker' in text and 'Version 26.2' in text and not re.search(r'(?i)not converged|undefined|fatal|error termination|\bNaN\b|\bInfinity\b',text))
            energy_tokens=re.findall(r'Total Potential Energy\s*:\s*('+NUMBER+')',text)
            check(label+' one finite energy',len(energy_tokens)==1 and math.isfinite(number(energy_tokens[0])))
            energy=number(energy_tokens[0]);check(label+' stored raw energy',energy==row['value']['energy'])
            summaries=re.findall(r'Induced Dipoles\s*:\s*Iterations\s+(\d+)\s+RMS Residual\s+('+NUMBER+')',text)
            check(label+' every induction summary',len(summaries)>0 and len(summaries)==len(re.findall(r'Induced Dipoles\s*:\s*Iterations',text)) and all(int(i)>=1 and math.isfinite(number(v)) and 0<=number(v)<=1e-8 for i,v in summaries))
            residuals.extend(number(v) for _,v in summaries)
            check(label+' induction report',row['induction']['solutions']==len(summaries) and row['induction']['maximum_printed_RMS_residual_Debye']==max(number(v) for _,v in summaries))
            original_key=verify(fixture['input']['key']).read_text();expected_key=original_key
            if profile in PROFILES:
                dd,dm=PROFILES[profile]
                for metal in (500,501):expected_key+=f'\nPOLAR-PAIR-RESPONSE {metal} 36 {dd:.17g} {dm:.17g}\n'
            check(label+' exact key and parameter input',verify(row['input']['key']).read_text()==expected_key and row['input']['parameters']['sha256']==fixture['input']['parameters']['sha256'])
            if 'q' not in job:check(label+' exact consumed coordinates',row['input']['xyz']['sha256']==fixture['input']['xyz']['sha256'])
            actual=re.findall(r'PAIR RESPONSE TYPES/GAMMA0/D/M:\s+(\d+)\s+(\d+)\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+')',text)
            absent=re.findall(r'PAIR RESPONSE ABSENT TYPES:\s+(\d+)\s+(\d+)',text)
            if profile not in PROFILES:
                check(label+' disabled coefficient inventory',not actual and not absent)
            else:
                check(label+' enabled native coefficient inventory',len(actual)+len(absent)==2 and {(int(x[0]),int(x[1])) for x in actual+absent}=={(500,36),(501,36)})
                for a,z,*tokens in actual:
                    a,z=int(a),int(z);check(label+' represented types',a in fixture['types'] and z in fixture['types'])
                    gamma=[float(t[3]) for t in (ln.split() for ln in original_key.splitlines()) if len(t)>=4 and t[0].lower()=='polpair' and {int(t[1]),int(t[2])}=={a,z}][-1]
                    expected=(gamma,gamma*(1+dd),gamma*(1+dm));got=[number(t) for t in tokens]
                    check(label+' native coefficient application',all(math.isfinite(x) and abs(x-y)<=1e-14*max(1,abs(y)) for x,y in zip(got,expected)))
                for a,z in absent:check(label+' native absent type',int(a) not in fixture['types'] or int(z) not in fixture['types'])
                readbacks+=1
            gradient=[]
            if program=='testgrad':
                blocks={'Anlyt':[],'Numer':[]};tail=text.split('Cartesian Gradient Breakdown over Individual Atoms :',1)[1]
                for line in tail.splitlines():
                    m=re.fullmatch(r'\s*(Anlyt|Numer)\s+(\d+)\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s*',line)
                    if m:
                        check(label+' atom ordering '+m[1]+' '+m[2],int(m[2])==len(blocks[m[1]])+1)
                        blocks[m[1]].append([number(m[k]) for k in (3,4,5)])
                check(label+' full finite arrays',len(blocks['Anlyt'])==fixture['atoms'] and len(blocks['Numer'])==(fixture['atoms'] if job['fd'] else 0) and all(math.isfinite(x) for ar in blocks.values() for v in ar for x in v))
                gradient=blocks['Anlyt'];check(label+' independent gradient/force sign',gradient==row['value']['gradient'] and [[-x for x in v] for v in gradient]==row['value']['forces'])
                if job['fd']:
                    error=max(abs(x-y) for ar,nr in zip(gradient,blocks['Numer']) for x,y in zip(ar,nr));fd_errors[label]=error
                    check(label+' raw Cartesian derivative',error<=.001 and error==row['value']['max_component_error_kcal_mol_A'] and blocks['Numer']==row['value']['numerical_gradient'])
            values[label]={'energy':energy,'gradient':gradient}
            if 'identity' in row:
                origin=values['original__'+job['cell']];de=abs(energy-origin['energy']);df=max(abs(x-y) for ar,nr in zip(gradient,origin['gradient']) for x,y in zip(ar,nr))
                identities[label]={'energy_error':de,'force_max_error':df,'printed_exact':de==0 and df==0}
                check(label+' raw identity',de<=1e-4 and df<=1e-3 and row['identity']=={'energy_error_kcal_mol':de,'force_max_error_kcal_mol_A':df,'printed_exact':de==0 and df==0,'passed':True})
            wall.append(receipt['elapsed_seconds']);user.append(receipt['child_user_cpu_seconds']);system.append(receipt['child_system_cpu_seconds'])
        directions={}
        for name,record in r['directional_checks'].items():
            g,kind=name.rsplit('__',1);coords=b['rows'][g]['positions_A'];tangent=[[0.,0.,0.] for _ in coords]
            if kind=='radial':
                delta=[coords[1][i]-coords[0][i] for i in range(3)];length=math.sqrt(sum(x*x for x in delta));u=[x/length for x in delta]
                tangent[0]=[-x/2 for x in u]
                for a in (1,2,3):tangent[a]=[x/2 for x in u]
            else:
                axis=[x/math.sqrt(14) for x in (1.,2.,3.)];raw=[]
                for a in (1,2,3):
                    v=[coords[a][i]-coords[1][i] for i in range(3)];raw.append([axis[1]*v[2]-axis[2]*v[1],axis[2]*v[0]-axis[0]*v[2],axis[0]*v[1]-axis[1]*v[0]])
                length=math.sqrt(sum(x*x for v in raw for x in v))
                for a,v in zip((1,2,3),raw):tangent[a]=[x/length for x in v]
            derivative=sum(x*t for ar,tr in zip(values['combined__'+g]['gradient'],tangent) for x,t in zip(ar,tr));steps=[]
            for h,saved in zip((1e-4,5e-5),record['steps']):
                for sign in (-1,1):
                    label=f'directional__{g}__{kind}__{h:g}__{sign:+d}';row=r['rows'][label];distance=sign*h
                    original=[ln.split() for ln in verify(b['rows'][g]['input']['xyz']).read_text().splitlines()[1:] if ln.strip()]
                    displaced=[ln.split() for ln in verify(row['input']['xyz']).read_text().splitlines()[1:] if ln.strip()]
                    expected=[v[:] for v in coords]
                    if kind=='radial':expected=[[x+distance*t for x,t in zip(v,tr)] for v,tr in zip(coords,tangent)]
                    else:
                        theta=distance/length
                        for a in (1,2,3):
                            v=[coords[a][i]-coords[1][i] for i in range(3)];velocity=[axis[1]*v[2]-axis[2]*v[1],axis[2]*v[0]-axis[0]*v[2],axis[0]*v[1]-axis[1]*v[0]];dot=sum(x*y for x,y in zip(axis,v))
                            expected[a]=[coords[1][i]+v[i]*math.cos(theta)+velocity[i]*math.sin(theta)+axis[i]*dot*(1-math.cos(theta)) for i in range(3)]
                    check(label+' exact path graph and coordinates',len(displaced)==len(original) and all(a[:2]+a[5:]==z[:2]+z[5:] for a,z in zip(original,displaced)) and max(abs(float(v)-e) for line,ex in zip(displaced,expected) for v,e in zip(line[2:5],ex))<1e-14)
                    check(label+' tangent metadata',max(abs(x-y) for tr,actual in zip(tangent,row['direction']['tangent']) for x,y in zip(tr,actual))<1e-15)
                lo=values[f'directional__{g}__{kind}__{h:g}__-1']['energy'];hi=values[f'directional__{g}__{kind}__{h:g}__+1']['energy'];slope=(hi-lo)/(2*h);bound=1e-8/(2*h)+5e-9*sum(abs(x) for v in tangent for x in v);error=abs(slope-derivative)
                check(name+f' independent directional {h:g}',abs(derivative-saved['analytic_slope'])<1e-12 and slope==saved['finite_difference_slope'] and abs(bound-saved['printing_bound'])<1e-16 and error<=.001+bound and abs(error-saved['absolute_error'])<1e-12)
                steps.append({'half_step_A':h,'analytic_slope':derivative,'FD_slope':slope,'error':error,'printing_bound':bound})
            difference=abs(steps[0]['FD_slope']-steps[1]['FD_slope']);check(name+' independent two-step',difference==record['two_step_difference'] and difference<=.001+sum(x['printing_bound'] for x in steps))
            directions[name]={'steps':steps,'two_step_difference':difference}
        check('all thresholds and counters accounted',len(fd_errors)==10 and len(identities)==21 and readbacks==38 and len(directions)==4)
        check('serial child CPU sum',abs(sum(user)-r['cost']['child_user_CPU_seconds'])<1e-8 and abs(sum(system)-r['cost']['child_system_CPU_seconds'])<1e-8)
        report.update(status='actual_numerical_qualification_review_passed',raw_outputs_admitted=56,Cartesian_FD_programs=10,internal_FD_energies=384,identity_comparisons=identities,all_identity_comparisons_printed_exact=all(v['printed_exact'] for v in identities.values()),Cartesian_errors=fd_errors,maximum_Cartesian_error_kcal_mol_A=max(fd_errors.values()),directional_checks=directions,maximum_directional_error_kcal_mol_A=max(s['error'] for d in directions.values() for s in d['steps']),maximum_induction_RMS_Debye=max(residuals),induction_summaries=len(residuals),enabled_programs_with_native_coefficients=readbacks,program_wall_seconds_sum=sum(wall),native_child_user_CPU_seconds=sum(user),native_child_system_CPU_seconds=sum(system),original_application_cost=r['cost'],build_cost={'wall_seconds':build['wall_seconds'],'child_cpu_seconds':build['child_cpu_seconds']},scope='Numerical consistency only on the declared isolated CPU route, fixtures and probes; no physical efficacy or wider engine qualification.')
    except Exception as exc:report.update(status='actual_numerical_review_failed',error=repr(exc))
    finally:
        report.update(check_count=len(checks),review_cost={'wall_seconds':time.monotonic()-started,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1})
        (out/'REVIEW.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('status','check_count','error')}))
    if report['error']:raise SystemExit(1)


if __name__=='__main__':main()
