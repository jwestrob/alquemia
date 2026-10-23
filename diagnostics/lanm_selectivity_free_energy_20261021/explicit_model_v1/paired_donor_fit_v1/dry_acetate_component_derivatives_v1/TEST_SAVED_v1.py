#!/usr/bin/env python3
"""Bounded saved/synthetic checks. No native executable is invoked."""
import copy
from decimal import Decimal
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace

H=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('component_under_test',H/'RUN.py');a=importlib.util.module_from_spec(sp);sp.loader.exec_module(a)
manifest_path=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time();checks={}
def check(name,ok):checks[name]=bool(ok);assert ok,name
def rejects(name,fun):
    try:fun()
    except Exception:check(name,True)
    else:check(name,False)

m,helper,_=a.load_manifest(manifest_path)
check('prepared16_same_motion_inventory',len(m['rows'])==16 and set(r['direction'] for r in m['rows'].values())=={'radial','rotation'})
for metal in a.METALS:
    e=m['central'][metal]['energy'];check('actual_'+metal+'_eight_terms_and_rounding',len(e['interactions'])==8 and abs(Decimal(e['component_sum_remainder_kcal_mol']))<=Decimal(e['component_sum_rounding_bound_kcal_mol']))
old=Path(m['central']['La']['energy']['native']['analyze_output']['path']).read_text()
for label,bad in [('nonfinite',old.replace('-367.44333349','NaN')),('component_closure',old.replace('0.96620352','1.96620352')),('duplicate_total',old+'\n Total Potential Energy : 0.00000000\n')]:
    p=out/(label+'.out');p.write_text(bad);rejects(label,lambda p=p:a.parse_energy(p,helper))
# Synthetic polynomial components, independent exact known slopes and unequal metal signs.
sm=copy.deepcopy(m);rows={}
for metal in a.METALS:
    for direction in a.DIRECTIONS:
        sm['central'][metal]['analytic'][direction]=(1 if metal=='La' else -2)*(1 if direction=='radial' else 3)*36
        sm['central'][metal]['analytic_gradient_print_bound_kcal_mol_A'][direction]=0
for cid,r in sm['rows'].items():
    q=Decimal(str(r['q_A']));factor=(1 if r['metal']=='La' else -2)*(1 if r['direction']=='radial' else 3)
    tokens={term:format(Decimal(i)+factor*i*q+Decimal(i)*q*q,'.8f') for i,term in enumerate(a.COMPONENTS,1)}
    tokens['Total']=format(sum(Decimal(v) for v in tokens.values()),'.8f');rows[cid]={'status':'complete','parsed':{'tokens':tokens}}
c=a.compare(sm,rows);check('synthetic_full18_total_checks',c['numerical_checks_passed'] and c['numerical_check_count']==18)
for direction,factor in [('radial',1),('rotation',3)]:
    v=c['directions'][direction]['paired']['steps']['0.005'];check('synthetic_'+direction+'_Dy_minus_La',v['derivatives_kcal_mol_A']['Total']==-108*factor)
    check('synthetic_'+direction+'_signed_force',v['generalized_forces_kcal_mol_A']['Atomic Multipoles']==63*factor)
    check('synthetic_'+direction+'_component_closure',v['component_sum_remainder_kcal_mol_A']==0)
check('absolute_coefficient_paired_print_bound',c['directions']['radial']['paired']['steps']['0.005']['print_bounds_kcal_mol_A']['Total']==2e-6)
bad=copy.deepcopy(sm);bad['central']['La']['analytic']['radial']+=.2
check('analytical_failure_retained',not a.compare(bad,rows)['numerical_checks_passed'])
missing=copy.deepcopy(rows);missing.pop(next(iter(missing)));rejects('missing16_denominator',lambda:a.compare(sm,missing))
nonfinite=copy.deepcopy(rows);nonfinite[next(iter(nonfinite))]['parsed']['tokens']['Total']='NaN';rejects('nonfinite_arithmetic',lambda:a.compare(sm,nonfinite))
# Exercise the actual serial driver with synthetic invocation only; native executable never called.
orig_load,orig_gate=a.load_manifest,a.release_gate
mockm=copy.deepcopy(m)
for metal in a.METALS:
    for direction in a.DIRECTIONS:mockm['central'][metal]['analytic'][direction]=0
release_path=out/'SYNTHETIC_RELEASE.json';a.save(release_path,{'status':'synthetic_control_only'})
old_env={k:os.environ.get(k) for k in (*a.THREADS,'SLURM_JOB_ID','SLURM_NTASKS','SLURM_CPUS_PER_TASK')}
os.environ.update({k:'1' for k in a.THREADS});os.environ.update(SLURM_JOB_ID='999000',SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='1')
try:
    for failure in (False,True):
        calls=[]
        def fake_invoke(argv,cwd,output,receipts):
            calls.append(str(cwd));receipts.append({'returncode':0,'elapsed_seconds':0,'threads':{k:'1' for k in a.THREADS}})
            if failure:raise RuntimeError('Explicit synthetic first-call failure')
            metal=mockm['rows'][Path(cwd).name]['metal'];output.write_text(Path(mockm['central'][metal]['energy']['native']['analyze_output']['path']).read_text())
        a.load_manifest=lambda _: (copy.deepcopy(mockm),helper,SimpleNamespace(invoke=fake_invoke));a.release_gate=lambda *args:{}
        dest=out/('mock_failure' if failure else 'mock_success');rc=a.run(manifest_path,release_path,dest);res=json.loads((dest/'RESULT.json').read_text())
        if failure:
            check('mock_first_failure_stops',len(calls)==res['new_native_calls']==1 and rc==1)
            check('mock_failure_retains16_slots',len(res['rows'])==16 and sum(v['status']=='unattempted' for v in res['rows'].values())==15 and res['comparison'] is None)
        else:
            check('mock_exact16_invocations',len(calls)==res['new_native_calls']==res['complete_cells']==16 and rc==0)
            check('mock_all18_arithmetic',res['comparison']['numerical_checks_passed'] and res['comparison']['numerical_check_count']==18)
finally:
    a.load_manifest,a.release_gate=orig_load,orig_gate
    for key,val in old_env.items():
        if val is None:os.environ.pop(key,None)
        else:os.environ[key]=val
# Concrete release rejects unbound and wrong native budget, with no reverse config hash.
rm=copy.deepcopy(m);rm['_path']=str(manifest_path);review_path=out/'SYNTHETIC_REVIEW.json';a.save(review_path,{'status':'passed_dry_acetate_component_source_review','implementation':a.pin(H/'RUN.py'),'manifest':a.pin(manifest_path)})
release={'status':'released_dry_acetate_component_derivatives','implementation':a.pin(H/'RUN.py'),'manifest':a.pin(manifest_path),'analyze':m['analyze'],'selected_review':m['selected_review'],'maximum_native_calls':16,'workers':1,'cell_ids':list(m['rows']),'output_directory':str(out/'never_run'),'independent_source_review':a.pin(review_path)}
p=out/'SYNTHETIC_BOUND_RELEASE.json';a.save(p,release);check('concrete_release_binding',a.release_gate(rm,p,out/'never_run')==release)
release['maximum_native_calls']=17;a.save(p,release);rejects('release_overcap',lambda:a.release_gate(rm,p,out/'never_run'))
release['maximum_native_calls']=16;release['status']='unreleased';a.save(p,release);rejects('unreleased_guard',lambda:a.release_gate(rm,p,out/'never_run'))
result={'status':'passed_saved_controls','implementation':a.pin(H/'RUN.py'),'controls':a.pin(__file__),'manifest':a.pin(manifest_path),'checks':checks,'check_count':len(checks),'all_passed':all(checks.values()),'synthetic_driver_invocations':17,'actual_new_native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'process_CPU_seconds':time.process_time()-cpu}}
a.save(out/'RESULT.json',result);print(json.dumps(result,indent=2))
