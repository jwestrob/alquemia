#!/usr/bin/env python3
"""Small synthetic field arithmetic/serial staging checks; zero native calls."""
import copy,importlib.util,json,os,re,sys,time
from decimal import Decimal
from pathlib import Path
from types import SimpleNamespace
H=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('dry_field_under_test',H/'RUN.py');a=importlib.util.module_from_spec(sp);sp.loader.exec_module(a)
mp=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time();checks={}
def check(k,x):checks[k]=bool(x);assert x,k
def rejects(k,f):
    try:f()
    except Exception:check(k,True)
    else:check(k,False)
m,comp,parser,native=a.load_manifest(mp)
check('actual10_zeros40_fieldkeys',len(m['rows'])==10 and len(m['jobs'])==40)
check('exact_charge2_all10',all(r['charge_origin']['charge_e']=='2.00000' for r in m['rows'].values()))
check('native_counts_field8_increment',all(r['expected_field_interactions']['Atomic Multipoles']==r['expected_field_interactions']['Polarization']==21 for r in m['rows'].values()))
row=next(iter(m['rows'].values()));key=a.verify(row['input']['key']).read_text();prm=a.verify(row['input']['parameters']).read_text();units=SimpleNamespace(FIELD_FACTOR=m['conversion']['FIELD_FACTOR'])
rejects('reject_preexisting_field',lambda:a.render(key+'\nexternal-field 0 0 1\n',prm,.001,m['direction'],units))
rejects('reject_unknown_field_step',lambda:a.render(key,prm,.002,m['direction'],units))
rejects('reject_nonunit_direction',lambda:a.render(key,prm,.001,[2.,0.,0.],units))
# Synthetic linear moment with cubic geometry dependence; no physical model fit.
def synthetic_token(jid,central_failure=False):
    job=m['jobs'][jid];r=m['rows'][job['geometry_id']];q=Decimal(str(r['q_A']));sign=1 if r['metal']=='La' else -1
    mu=Decimal(sign)*(Decimal(2)+3*q+10000*q**3);chi=Decimal(20 if sign==1 else 22)+q;f=Decimal(str(job['signed_field_au']));K=Decimal(str(m['conversion']['HARTREE']));e0=Decimal(r['zero_energy']['tokens']['Total'])
    e=e0-K*(mu*f+chi*f*f/2)
    if central_failure and jid=='La_central__plus_h':e+=K*Decimal('.00001')
    return format(e,'.8f'),float(mu),float(chi)
rows={jid:{'status':'complete','parsed':{'tokens':{'Total':synthetic_token(jid)[0]}}} for jid in m['jobs']}
s=a.summarize(m,rows);check('synthetic31_consistency_checks',s['numerical_check_count']==31 and s['numerical_checks_passed']);check('synthetic_central8_gate',a.central_gate(m,rows)['passed'])
for gid,v in s['states'].items():
    _,mu,chi=synthetic_token(gid+'__plus_h')
    check('signed_moment_'+gid,abs(v['steps']['half']['mu_E_au']-mu)<=v['steps']['half']['mu_print_bound_au']+1e-12)
    check('curvature_'+gid,abs(v['steps']['half']['chi_E_a0_cubed']-chi)<=v['steps']['half']['chi_print_bound_a0_cubed']+1e-9)
for q in (.01,.005):
    expected=-2*(3+10000*q*q);d=s['rotation_dipole_derivatives']['paired'][str(q)]['field_steps']['half'];check('paired_geometry_sign_'+str(q),abs(d['dmu_dq_ea0_per_A']-expected)<=d['print_bound_ea0_per_A']+1e-10)
check('q_step_difference_descriptive',abs(s['rotation_dipole_derivatives']['La']['smaller_minus_larger_q_descriptive']['half']['dmu_dq_ea0_per_A'])>.7 and s['numerical_checks_passed'])
missing=copy.deepcopy(rows);missing.pop(next(iter(missing)));rejects('missing40_no_summary',lambda:a.summarize(m,missing))
bad=copy.deepcopy(rows);bad[next(iter(bad))]['parsed']['tokens']['Total']='NaN';rejects('nonfinite_no_summary',lambda:a.summarize(m,bad))
bad=copy.deepcopy(rows);bad['La_central__plus_h']['parsed']['tokens']['Total']=synthetic_token('La_central__plus_h',True)[0];check('central_failure_detected',not a.central_gate(m,bad)['passed'])
# Exercise actual serial runner with mocked invocation only. Preserve native-format E/component tokens.
def native_text(jid,token):
    r=m['rows'][m['jobs'][jid]['geometry_id']];z=r['zero_energy'];text=Path(z['native']['analyze_output']['path']).read_text();delta=Decimal(token)-Decimal(z['tokens']['Total']);pol=format(Decimal(z['tokens']['Polarization'])+delta,'.8f')
    text=re.sub(r'(Total Potential Energy\s*:\s*)'+parser.NUMBER,lambda z:z[1]+token,text)
    for name in ('Atomic Multipoles','Polarization'):
        val=z['tokens'][name] if name=='Atomic Multipoles' else pol
        text=re.sub(r'^\s*'+re.escape(name)+r'\s+'+parser.NUMBER+r'\s+\d+\s*$',lambda _:f' {name} {val} 21',text,flags=re.M)
    return text
origload,origgate=a.load_manifest,a.release_gate;oldenv={k:os.environ.get(k) for k in (*a.THREADS,'SLURM_JOB_ID','SLURM_NTASKS','SLURM_CPUS_PER_TASK')};os.environ.update({k:'1' for k in a.THREADS});os.environ.update(SLURM_JOB_ID='999001',SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='1');rp=out/'SYNTHETIC_RELEASE.json';a.save(rp,{'status':'synthetic_only'})
try:
    for mode in ('success','central_failure','native_failure'):
        calls=[]
        def invoke(argv,cwd,output,receipts):
            calls.append(Path(cwd).name);receipts.append({'returncode':0,'elapsed_seconds':0})
            if mode=='native_failure':raise RuntimeError('Explicit synthetic native failure')
            token=synthetic_token(Path(cwd).name,mode=='central_failure')[0];output.write_text(native_text(Path(cwd).name,token))
        a.load_manifest=lambda _: (copy.deepcopy(m),comp,parser,SimpleNamespace(invoke=invoke));a.release_gate=lambda *args:{};dest=out/('mock_'+mode);rc=a.run(mp,rp,dest);r=json.loads((dest/'RESULT.json').read_text())
        if mode=='success':check('actual_driver_mock40_full_summary',len(calls)==r['new_native_calls']==r['complete_cells']==40 and rc==0 and r['summary']['numerical_checks_passed'])
        elif mode=='central_failure':check('actual_driver_mock8_then_stop32',len(calls)==r['new_native_calls']==8 and rc==1 and r['central_gate']['passed'] is False and sum(v['status']=='unattempted' for v in r['rows'].values())==32)
        else:check('actual_driver_mock_first_failure_full40',len(calls)==r['new_native_calls']==1 and rc==1 and len(r['rows'])==40 and sum(v['status']=='unattempted' for v in r['rows'].values())==39)
finally:
    a.load_manifest,a.release_gate=origload,origgate
    for k,v in oldenv.items():
        if v is None:os.environ.pop(k,None)
        else:os.environ[k]=v
review=out/'SYNTHETIC_REVIEW.json';a.save(review,{'status':'passed_dry_acetate_electrical_source_review','implementation':a.pin(H/'RUN.py'),'manifest':a.pin(mp)});release={'status':'released_dry_acetate_electrical_response','implementation':a.pin(H/'RUN.py'),'manifest':a.pin(mp),'analyze':m['analyze'],'shared_specification':m['shared_specification'],'maximum_native_calls':40,'workers':1,'cell_ids':list(m['jobs']),'output_directory':str(out/'never_run'),'independent_source_review':a.pin(review)};a.save(rp,release);check('concrete_release_binding',a.release_gate(m,mp,rp,out/'never_run')==release);release['maximum_native_calls']=41;a.save(rp,release);rejects('reject_41_scope',lambda:a.release_gate(m,mp,rp,out/'never_run'))
r={'status':'passed_saved_controls','implementation':a.pin(H/'RUN.py'),'controls_source':a.pin(__file__),'manifest':a.pin(mp),'checks':checks,'check_count':len(checks),'all_passed':all(checks.values()),'actual_new_native_calls':0,'synthetic_driver_invocations':49,'cost':{'wall_seconds':time.monotonic()-tick,'process_CPU_seconds':time.process_time()-cpu}};a.save(out/'RESULT.json',r);print(json.dumps(r,indent=2))
