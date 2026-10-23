#!/usr/bin/env python3
"""One-thread saved-only actual raw admission and independent Decimal cycle."""
from pathlib import Path
from decimal import Decimal, getcontext
from collections import Counter
import hashlib
import importlib.util
import json
import math
import os
import re
import sys
import time

START=time.monotonic(); CPU=time.process_time(); getcontext().prec=40
ROOT=Path(__file__).resolve().parents[5]
D=Path(__file__).resolve().parent.parent
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/matched_pair_decomposition_v1'
OUT=W/'independent_review_v1'
HA=Decimal('627.509474')
GIDS=tuple(f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy'))
SIGNS=dict(zip(GIDS,(1,-1,-1,1)))
TERMS=('electronic_energy_hartree','D4_energy_hartree','final_energy_hartree')
CHECKS={}

def ck(name,passed):
    if not passed: raise AssertionError(name)
    CHECKS[name]=True
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):return json.loads(Path(p).read_text())
def verify(p):
    ck('hash:'+p['path'],pin(p['path'])==p);return Path(p['path'])
def module(p,name):
    spec=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def dec_token(text,pattern):
    vals=re.findall(pattern,text,re.M)
    ck('single_token:'+hashlib.sha256((text+pattern).encode()).hexdigest()[:14],len(vals)==1)
    return Decimal(vals[0].replace('D','E').replace('d','e'))
NUM=r'([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)'
def quantum_tokens(output):
    text=verify(output).read_text().split('****END OF INPUT****')[-1]
    return dict(zip(TERMS,[dec_token(text,r'^\s*Total Energy\s*:\s*'+NUM+r'\s+Eh'),
        dec_token(text,r'^\s*Dispersion correction\s+'+NUM+r'\s*$'),
        dec_token(text,r'^\s*FINAL SINGLE POINT ENERGY\s+'+NUM+r'\s*$')]))
def model_token(output):
    return dec_token(verify(output).read_text(),r'^\s*Total Potential Energy\s*:\s*'+NUM)
def decompose(clusters,pairs,ions):
    sums={g:sum(pairs[g].values(),Decimal(0)) for g in GIDS}
    c=sum((SIGNS[g]*clusters[g] for g in GIDS),Decimal(0))
    p=sum((SIGNS[g]*sums[g] for g in GIDS),Decimal(0))-(ions['Dy']-ions['La'])
    return {'cluster_contrast':c,'pair_contrast':p,'cooperative_remainder':c-p}
def numeric(x):
    if isinstance(x,Decimal):return float(x)
    if isinstance(x,dict):return {k:numeric(v) for k,v in x.items()}
    if isinstance(x,list):return [numeric(v) for v in x]
    return x

ck('one_thread',all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')))
rp=W/'collection/run_v1/result/RESULT.json';r=read(rp);cfg=read(verify(r['configuration']))
ck('terminal_complete',r['status']=='matched_pair_collection_complete' and r['error'] is None)
ck('complete_counts',(r['reference_complete'],r['model_complete'],r['original_complete'])==(36,34,4))
ck('unchanged_scope',r['new_native_calls']==0 and not r['fit_targets_changed'] and not r['model_selected'] and not r['fresh_validation_released'])
sys.path.insert(0,str(D/'collection'))
col=module(cfg['implementation'],'actual_collector');ref=module(cfg['reference_adapter'],'actual_ref');panel=ref.load(verify(cfg['reference_panel']))
ck('all36_ids',set(r['reference_rows'])=={c['cell_id'] for c in panel['cells']})
qr={}; states={}; cosine={}; replay=[]
for cell in panel['cells']:
    cid=cell['cell_id'];row=r['reference_rows'][cid]
    actual=col.admit_reference(cfg['reference_result_paths'][cid],cell,cfg['reference_panel'],panel,ref)
    ck('reference_raw_replay:'+cid,actual==row)
    q=quantum_tokens(row['output']);qr[cid]=q
    ck('reference_raw_values:'+cid,all(float(q[k])==row[k] for k in TERMS))
    ck('reference_printed_D4_closure:'+cid,abs(q[TERMS[0]]+q[TERMS[1]]-q[TERMS[2]])<Decimal('1e-8'))
    if cell['role']=='bare_ion':ck('bare_D4_zero:'+cid,q[TERMS[1]]==0)
    states[cid]={'native_stop':row['native_scf_stop_signals'],'printed_criteria':row['printed_scf_criteria'],
        'final_COSX_correction_hartree':row['final_COSX_correction_hartree'],'effective_multiplicity':cell['multiplicity'],
        'physical_multiplicity':cell['physical_multiplicity'],'expected_electrons':cell['expected_electrons'],
        'expected_basis_functions':sum(cell['expected_basis_functions_by_center'])}
    replay.append({'cell_id':cid,'result':row['result'],'input':row['input'],'output':row['output'],'start':row['start'],'release':row['release']})
root=read(verify(panel['root_inputs']));cppanel=ref.helper.load(verify(root['cp_panel']));data=read(verify(root['consumed_input_data']));old=module(cfg['original_collector'],'actual_original')
qclusters={}
for g in GIDS:
    ck('original_raw_replay:'+g,old.admit_original(g,cppanel,data,ref.helper)==r['original_clusters'][g])
    qclusters[g]=quantum_tokens(r['original_clusters'][g]['output'])
model=module(cfg['model_adapter'],'actual_model');mcfg=read(verify(cfg['model_config']));manifest=read(verify(cfg['model_manifest']));mr=read(verify(r['model_result']))
admitted=model.admit_result(mcfg,manifest,mr)
ck('model_all34_replayed',admitted['rows']==r['model_rows'] and len(admitted['rows'])==34 and mr['native_program_invocations']==34)
ck('actual_model_job',mr['job_id']=='17819')
mc={g:model_token(admitted['saved_cluster_models'][g]['native']['analyze_output']) for g in GIDS}
mp={g:{} for g in GIDS};qp={k:{g:{} for g in GIDS} for k in TERMS};qi={k:{} for k in TERMS}
for cell in panel['cells']:
    cid=cell['cell_id'];q=qr[cid]
    if cell['role']=='bare_ion':
        for k in TERMS:qi[k][cell['metal']]=q[k]
    else:
        g=cell['parent_geometry_id'];i=cell['water_index'];mp[g][i]=model_token(r['model_rows'][cid]['native']['analyze_output'])
        ck('model_printed_energy:'+cid,float(mp[g][i])==r['model_rows'][cid]['energy_kcal_mol'])
        for k in TERMS:qp[k][g][i]=q[k]
for g in GIDS:ck('literal_water_denominator:'+g,set(mp[g])==set(range(1,int(g[4])+1)))
qcycles={k:{x:v*HA for x,v in decompose({g:qclusters[g][k] for g in GIDS},qp[k],qi[k]).items()} for k in TERMS}
mcycle=decompose(mc,mp,{'La':Decimal(0),'Dy':Decimal(0)})
for k in TERMS:
    for term,value in qcycles[k].items():ck('quantum_decimal_cycle:'+k+term,abs(float(value)-r['reference_decomposition'][k]['kcal_mol'][term])<1e-8)
for term,value in mcycle.items():ck('model_decimal_cycle:'+term,abs(float(value)-r['model_decomposition'][term])<1e-8)
errors={k:mcycle[k]-qcycles[TERMS[2]][k] for k in mcycle}
ck('exact_decimal_error_closure',errors['cluster_contrast']==errors['pair_contrast']+errors['cooperative_remainder'])
for term,value in errors.items():ck('reported_errors:'+term,abs(float(value)-r['comparison'][term])<1e-8)
water=[]
for n in (8,9):
    sign=-1 if n==8 else 1;lg=f'aquo{n}_central_La';dg=f'aquo{n}_central_Dy'
    for i in range(1,n+1):
        q={k:(qp[k][dg][i]-qp[k][lg][i]-(qi[k]['Dy']-qi[k]['La']))*HA for k in TERMS}
        m=mp[dg][i]-mp[lg][i]
        water.append({'water_count':n,'water_index':i,'cycle_sign':sign,'reference_metal_exchange_kcal_mol':q,
            'model_metal_exchange_kcal_mol':m,'model_minus_reference_kcal_mol':m-q[TERMS[2]],
            'signed_pair_cycle_reference_kcal_mol':sign*q[TERMS[2]],'signed_pair_cycle_model_kcal_mol':sign*m})
ck('all17_matched_contributions',len(water)==17)
ck('matched_reference_pair_closure',sum((x['signed_pair_cycle_reference_kcal_mol'] for x in water),Decimal(0))==qcycles[TERMS[2]]['pair_contrast'])
ck('matched_model_pair_closure',sum((x['signed_pair_cycle_model_kcal_mol'] for x in water),Decimal(0))==mcycle['pair_contrast'])
stop=Counter();flags=Counter()
for x in states.values():
    stop.update(x['native_stop'])
    flags.update(z['name'] for z in x['printed_criteria'] if not z['passes'])
events={};costs={}
for label,path,expected in [('references',W/'run_v1/completion_event.json',{str(x['job_id']) for x in r['reference_rows'].values()}),('model',W/'model/run_v1/completion_event.json',{'17819'}),('collector',W/'collection/run_v1/completion_event.json',{'17822'})]:
    event=read(path);records=event['records'];ck('scheduler_inventory:'+label,set(records)==expected)
    ck('scheduler_terminal:'+label,all(x['state']=='COMPLETED' and x['exit_code']=='0:0' for x in records.values()))
    snapshot=OUT/('actual_'+label+'_scheduler_records.json');assert not snapshot.exists();snapshot.write_text(json.dumps(records,indent=2)+'\n');events[label]=pin(snapshot)
    costs[label+'_scheduler_allocated_CPU_seconds']=sum(int(x['elapsed_seconds'])*int(x['allocated_cpus']) for x in records.values())
costs.update(reference_native_wall_times_ranks_proxy=sum(x['allocated_CPU_seconds_proxy'] for x in r['reference_rows'].values()),
    model_native_wall_seconds=mr['cost']['wall_seconds'],model_aggregate_child_CPU_seconds=mr['cost']['aggregate_child_user_CPU_seconds']+mr['cost']['aggregate_child_system_CPU_seconds'],
    model_parent_CPU_seconds=mr['cost']['parent_CPU_seconds'],collector_CPU_seconds=r['cost']['collector_CPU_seconds'],collector_wall_seconds=r['cost']['collector_wall_seconds'],
    reused_reference_proxy_excluded_from_new=sum(x['allocated_CPU_seconds'] for x in r['original_clusters'].values()))
costs['new_scheduler_total_CPU_seconds']=sum(costs[k] for k in ('references_scheduler_allocated_CPU_seconds','model_scheduler_allocated_CPU_seconds','collector_scheduler_allocated_CPU_seconds'))
ck('native_cost_replay',abs(costs['reference_native_wall_times_ranks_proxy']-r['cost']['reference_allocated_CPU_seconds_proxy'])<1e-8)
release=read(W/'collection/run_v1/ROOT_RELEASE.json');ck('saved_collector_release',release['configuration']==r['configuration'] and release['new_native_calls']==0 and release['threads']==1)
result={'status':'passed_actual_complete_raw_and_independent_decimal_review','implementation':pin(__file__),'declaration':pin(D/'independent_review_v1/ACTUAL_REVIEW_DECLARATION.json'),
    'result':pin(rp),'configuration':r['configuration'],'integration_review':pin(OUT/'INTEGRATION_REVIEW.json'),'model_result':r['model_result'],
    'reference_replay_inventory':replay,'quantum_decimal_cycles_kcal_mol':qcycles,'model_decimal_cycle_kcal_mol':mcycle,'model_minus_reference_kcal_mol':errors,
    'all17_matched_water_contributions':water,'new_native_stopping_metadata':states,'new_native_stop_counts':dict(stop),'new_native_failed_printed_criterion_counts':dict(flags),
    'scheduler_snapshots':events,'costs':costs,'checks':CHECKS,'check_count':len(CHECKS),
    'conclusion':'The reference has a larger positive cooperative remainder. A +1.263316 kcal/mol pair error partly cancels the -5.601600 remainder error, producing the -4.338284 full-shell error. This own-basis decomposition does not establish a unique missing physical term.',
    'limits':['Frozen PBE0-D4 ECP reference is provisional; native admission does not imply all printed SCF thresholds pass or physical state/basis/method accuracy.',
      'Each own-basis dimer is literal and contains full D4. Cooperative B includes basis sharing and environment-dependent/nonadditive dispersion and electronic effects; it is not a unique polarization or charge-transfer observable.',
      'The completed metal-whole-shell CP correction is kept separate, not inserted into B.',
      'W9 changes orientations of the first eight waters as well as coordination number; this is not a pure ninth-water insertion or bulk/free-energy qualification.',
      'Scheduler allocation cost is distinct from actual CPU utilization; reference native wall-times-ranks is only a proxy. Historical reused work is excluded from new costs.'],
    'new_native_calls':0,'threads':1,'review_wall_seconds':time.monotonic()-START,'review_CPU_seconds':time.process_time()-CPU}
target=OUT/'RESULT_REVIEW.json';assert not target.exists();target.write_text(json.dumps(numeric(result),indent=2,allow_nan=False)+'\n')
print(json.dumps({'review':pin(target),'checks':len(CHECKS),'errors':numeric(errors),'stop':dict(stop),'failed_criteria':dict(flags),'cost':costs},indent=2))
