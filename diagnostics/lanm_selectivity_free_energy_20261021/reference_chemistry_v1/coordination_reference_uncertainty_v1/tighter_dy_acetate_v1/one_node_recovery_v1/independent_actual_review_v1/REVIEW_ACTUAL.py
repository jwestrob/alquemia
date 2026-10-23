#!/usr/bin/env python3
"""Three saved native admissions and independent raw-token CP arithmetic."""
from decimal import Decimal
import hashlib, importlib.util, json, re, time
from pathlib import Path
HERE=Path(__file__).resolve().parent; BASE=HERE.parent
ROOT=next(p for p in HERE.parents if (p/'diagnostics').is_dir() and (p/'workspaces').is_dir())
W=ROOT/'workspaces'/BASE.relative_to(ROOT/'diagnostics'); OUT=W/'independent_actual_review_v1'
tick,cpu=time.monotonic(),time.process_time();checks={}
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def check(k,v):checks[k]=bool(v);assert v,k
def verify(p):check('pin/'+p['path'],pin(p['path'])==p);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
collection_path=W/'run_v1/collection_v1/RESULT.json';result=json.loads(collection_path.read_text())
sp=importlib.util.spec_from_file_location('unchanged_recovery_admitter',verify(result['implementation']));a=importlib.util.module_from_spec(sp);sp.loader.exec_module(a)
manifest=verify(result['manifest']);m,p,cell=a.load(manifest);release_path=verify(result['release']);release=a.release_contract(m,manifest,json.loads(release_path.read_text()))
check('complete3',result['status']=='dy_acetate_tighter_attribution_complete' and result['complete']==result['declared']==3)
check('no_promotion',result['full_eighteen_qualified'] is result['cluster_release_authorized'] is result['physical_accuracy_qualified'] is False)
number=r'[-+]?\d+\.\d+(?:[EeDd][-+]?\d+)?'
patterns={'CCSD':r'^\s*E\(CCSD\)\s+\.\.\.\s*('+number+')','triples':r'^\s*(?:Triples Correction \(T\)|Triples Correction \(T1\))\s+\.\.\.\s*('+number+')','total':r'^\s*E\(CCSD\(T\)\)\s+\.\.\.\s*('+number+')'}
def components(text):
    vals={}
    for key,pattern in patterns.items():
        v=re.findall(pattern,text,re.M);assert len(v)==1,(key,v);vals[key]=Decimal(v[0].replace('D','E'))
    return vals
sets={name:{k:Decimal(0) for k in patterns} for name in ('canonical','tight','tighter')};roles={}; admitted={}
for c in p['cells']:
    cid=c['cell_id'];saved=result['rows'][cid];native_path=verify(saved['result']);actual=json.loads(native_path.read_text())
    if cid==a.CID:q=a.admit_replacement(native_path,m,p,c,manifest,release_path)
    else:q=a.old.admit_result(native_path,c,p,verify(result['panel']))
    q=json.loads(json.dumps(q));check(cid+'/complete_saved_replay',q==saved);admitted[cid]=q
    raw=verify(saved['output']).read_text();vals=components(raw);sign=1 if c['role']=='complex' else -1
    check(cid+'/raw_component_values',float(vals['CCSD'])==q['CCSD_energies_hartree'][0] and float(vals['triples'])==q['triples_hartree'][0] and float(vals['total'])==q['energy_hartree'])
    for key,value in vals.items():sets['tighter'][key]+=sign*value*Decimal(str(a.b.HA))
    previous=p['baseline_components'][c['source_tight_cell_id']]
    for kind,sourcekey in [('canonical','canonical'),('tight','local')]:
        authority=previous[sourcekey];text=verify(authority['output']).read_text();values=components(text)
        for key,value in values.items():
            check(cid+'/'+kind+'/'+key+'/raw_token',value==Decimal(authority['values'][key]));sets[kind][key]+=sign*value*Decimal(str(a.b.HA))
    st=read(actual['start']);inp=verify(actual['input']).read_text()
    check(cid+'/ranks_memory',st['ranks']==8 and st['memory_MiB']==32000 and st['maxcore_MiB']==3000)
    check(cid+'/fixed_method',inp.startswith('! RHF DLPNO-CCSD(T1) TightPNO def2-TZVP def2-TZVP/C NoFrozenCore VeryTightSCF NoAutostart PrintBasis') and all(x in inp for x in [' STol 1e-8',' TriTolE 1e-8',' TCutPairs 1e-06',' TCutPNO 1e-08',' TCutMKN 0.0001',' TCutDO 0.005']))
    check(cid+'/effective_inventory',q['basis_functions_by_center']==c['expected_basis_functions_by_center'] and set(q['total_electrons'])=={c['expected_electrons']} and q['correlated_electrons']==[c['expected_correlated_electrons']])
    check(cid+'/CC_residual',len(q['CC_terminal_residuals'])==1 and q['CC_terminal_residuals'][0]<=q['CC_printed_tolerance'][0]==1e-8)
    check(cid+'/triples_residual',len(q['triples_terminal_deltas'])==1 and abs(q['triples_terminal_deltas'][0])<=q['triples_printed_energy_tolerance'][0]==1e-8)
    check(cid+'/iterative_T1', 'Parameters used in iterative triples' in raw and q['actual_triples_settings']['complete_and_finite'])
    check(cid+'/SCF_flags_retained',q['native_scf_stop_signals']==['Energy Check'] and any(not x['passes'] for x in q['printed_scf_criteria']))
    roles[c['role']]={k:q[k] for k in ('energy_hartree','CCSD_energies_hartree','triples_hartree','RHF_energies_hartree','CC_terminal_residuals','triples_terminal_deltas','T1_diagnostic','total_electrons','correlated_electrons','ecp_core_counts','basis_functions_by_center','aux_C','printed_scf_criteria','actual_triples_settings','native_local_thresholds')}
    roles[c['role']].update(job_id=actual['job_id'],result=saved['result'],input=saved['input'],output=saved['output'],raw_CCSD_plus_triples_minus_total_hartree=str(vals['CCSD']+vals['triples']-vals['total']))
    if cid==a.CID:
        check('replacement_job_and_one_node',actual['job_id']=='18693' and st['nodes']==1 and actual['native_program_invocations_attempted']==1)
        check('canceled_input_byte_identity',verify(actual['input']).read_bytes()==verify(m['original_input']).read_bytes())
    else:check(cid+'/reused_authority',saved['result']==m['original_results'][cid] and actual['job_id'] in ('18516','18517'))
replayed=a.old.comparison(p,admitted)
for k,v in replayed.items():check('frozen_comparison/'+k,v==result[k])
maxdifference=0.
for kind,values in sets.items():
    for key,value in values.items():
        difference=abs(float(value)-result['CP_components_kcal_mol'][kind][key]);maxdifference=max(maxdifference,difference);check('decimal_CP/'+kind+'/'+key,difference<=1e-9)
errors={kind:{key:values[key]-sets['canonical'][key] for key in patterns} for kind,values in sets.items() if kind!='canonical'}
check('original_point25_gate', (abs(errors['tighter']['total'])<=Decimal('.25'))==result['diagnostic_limit_passed'])
check('original_was_outside_point25',abs(errors['tight']['total'])>Decimal('.25'))
canceled=read(result['original_canceled_receipt'])['records']['18518'];check('cancellation_and_cost',canceled['state'].startswith('CANCELLED') and int(canceled['elapsed_seconds'])*int(canceled['allocated_cpus'])==result['abandoned_attempt_scheduler_CPU_seconds']==59512)
incomplete=read(result['original_incomplete_collection']);check('original_incomplete_preserved',incomplete['status']=='unavailable' and incomplete['complete']==2 and incomplete['declared']==3 and incomplete['rows'][a.CID]['status']=='unavailable')
check('selected_native_cost',sum(x['allocated_cpu_seconds'] for x in admitted.values())==result['selected_native_allocated_CPU_seconds'])
event_text=(W/'run_v1/completion_event.json').read_text();event=json.loads(event_text);endpoint=event['records']['18693'];check('new_scheduler_cost',endpoint['state']=='COMPLETED' and endpoint['exit_code']=='0:0' and int(endpoint['elapsed_seconds'])*int(endpoint['allocated_cpus'])==2169*8)
OUT.mkdir(parents=True,exist_ok=True);(OUT/'COMPLETION_SNAPSHOT.json').write_text(event_text)
review=dict(status='passed_tighter_dy_acetate_recovery_actual_review',implementation=result['implementation'],manifest=result['manifest'],collection=pin(collection_path),release=result['release'],reviewer_source=pin(__file__),checks=checks,check_count=len(checks),roles=roles,CP_components_decimal_kcal_mol={k:{a:str(b) for a,b in v.items()} for k,v in sets.items()},CP_errors_decimal_kcal_mol={k:{a:str(b) for a,b in v.items()} for k,v in errors.items()},max_decimal_vs_float_difference_kcal_mol=maxdifference,diagnostic_limit_passed=result['diagnostic_limit_passed'],full_eighteen_qualified=False,physical_accuracy_qualified=False,cluster_release_authorized=False,cost={'selected_native_allocated_CPU_seconds':result['selected_native_allocated_CPU_seconds'],'abandoned18518_scheduler_CPU_seconds':59512,'replacement18693_scheduler_elapsed_seconds':2169,'replacement18693_scheduler_allocated_CPU_seconds':17352,'collector_application_wall_seconds':result['wall_seconds'],'collector_application_CPU_seconds':result['CPU_seconds']},completion_snapshot=pin(OUT/'COMPLETION_SNAPSHOT.json'),saved_admission_calls=3,new_native_calls=0,review_wall_seconds=time.monotonic()-tick,review_CPU_seconds=time.process_time()-cpu)
(OUT/'REVIEW.json').write_text(json.dumps(review,indent=2,allow_nan=False)+'\n');print(json.dumps({'review':pin(OUT/'REVIEW.json'),'checks':len(checks),'decimal_error':review['CP_errors_decimal_kcal_mol'],'max_float_difference':maxdifference},indent=2))
