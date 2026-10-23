#!/usr/bin/env python3
"""Literal metadata inventory only; no render, parser, cycle or solver evaluation."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import socket
import time

HERE=Path(__file__).resolve().parent
def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(set(p)=={'path','sha256'} and pin(p['path'])==p,'Changed pin '+str(p));return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def inventory(out):
    tick,cpu=time.monotonic(),time.process_time();require(not out.exists(),'Fresh immutable output');out.mkdir(parents=True)
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')),'One thread')
    require(not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID'),'Truthful metadata-only login')
    checks=[]
    def check(ok,msg):require(ok,msg);checks.append(msg)
    result={'status':'metadata_inventory_started','implementation':pin(__file__),'inputs':pin(HERE/'SOURCE_INPUTS.json'),'error':None,'checks':checks,'new_native_calls':0,'scalar_replays':0,'heldout_reads':0}
    try:
        s=read(result['inputs'])
        for k,p in s.items():
            if k!='prior_regressions':verify(p)
        dec=read(s['root_declaration']);check(dec['status']=='declared_frozen_radial14_consumed_trimer_regression','Declared finite question')
        adm=read(s['terminal_root_admission']);prod=read(s['terminal_result']);selected=read(s['terminal_candidate']);fitcfg=read(s['fit_configuration']);m=read(s['parent_manifest']);em=read(s['evaluator_manifest']);target=read(s['target_inventory']);original=read(s['original_trimer_comparison'])
        check(adm['status']=='admitted_bounded_radial14_development_fit' and adm['producer']==s['terminal_result'] and adm['comparison']==s['terminal_comparison'],'Exact admitted terminal producer/comparison')
        check(adm['selected_artifact']==prod['selected']==dec['parents']['terminal_candidate']==s['terminal_candidate'],'Frozen5d992 candidate; pending witness excluded')
        check(prod['status']=='bounded_radial14_iterative_fit_complete' and prod['error'] is None and prod['candidate_frozen'],'Terminal producer complete')
        check(prod['configuration']==s['fit_configuration'] and digest(selected['parameters'])==dec['candidate_key']==adm['selected_key']==prod['selected_key'],'Literal parameter/configuration identity')
        check(len(selected['parameters'])==14 and selected['status']=='complete' and len(selected['cells'])==46 and all(r['status']=='complete' for r in selected['cells'].values()),'Complete admitted46 cells')
        check(fitcfg['parent_evaluator']==s['parent_evaluator'] and fitcfg['parent_configuration']==s['evaluator_configuration'],'Same original evaluator authority')
        for k in ('radial_adapter','radial_runner','build','native_source_manifest','qualification_review','qualification_root_review','root_cache_bridge','environment'):check(fitcfg[k]==s[k],'Same fit source '+k)
        check(fitcfg['executables']['analyze']==s['analyze'],'Same qualified radial analyze')
        check(m['evaluator_manifest']==s['evaluator_manifest'] and m['sources']['target_inventory']==s['target_inventory'],'Original regression/evaluator/target joins')
        group=target['groups']['trimer_energy_and_dipoles']
        check(group['cycles']==m['cycles'] and len(m['cycles'])==8,'Eight unchanged literal cycles')
        ids=[f'aquo{n}_central_{metal}__{pair}__metal_two_water' for n in (8,9) for metal in ('La','Dy') for pair in (('near','far') if n==8 else ('near',))]
        need=set().union(*(set(c['model_terms']) for c in m['cycles'].values()))
        selected_ids=sorted(need&set(selected['cells']));static_ids=sorted(need&set(em['static']));missing=need-set(selected_ids)-set(static_ids)
        check(missing==set(ids)==set(m['rows']) and len(ids)==6,'Exactly six missing analyze endpoints')
        check(len(need)==42 and len(selected_ids)==14 and len(static_ids)==22 and not(set(selected_ids)&set(static_ids)),'42 disjoint terms:14selected+22static+6new')
        check(all(f'aquo9_central_{metal}__far__metal_two_water' in selected_ids for metal in ('La','Dy')),'Both W9 far trimers reused')
        for cid in ids:
            row=m['rows'][cid]
            check(row['id']==cid and row['atoms']==7 and row['charge']==3 and row['types']==[500 if row['metal']=='La' else 501,36,37,37,36,37,37],'Exact vacuum trimer roles/types '+cid)
            for p in row['source_input'].values():verify(p)
        reused={}
        for cid in selected_ids:
            record=selected['cells'][cid];row=em['rows'][cid]
            check(record['status']=='complete' and all(r['returncode']==0 for r in record['receipts']),'Admitted selected dependency '+cid)
            reused[cid]={'artifact':s['terminal_candidate'],'json_path':['cells',cid],'definition':row,'source_inputs':row['source_input'],'program_inputs':record['program_inputs'],'receipts':record['receipts'],'cache_key':record['cache_key'],'cache_identity':record['cache_identity'],'native_admission_reused':s['terminal_root_admission'],'no_new_scalar_replay':True}
        refs=set().union(*(set(c['reference_terms']) for c in m['cycles'].values()))
        check(refs<=set(original['quantum_moment_rows']),'Complete original reference property rows')
        check(group['physical_acceptance_threshold'] is None and target['policy']['physical_state_or_reference_accuracy_qualified'] is False,'No new threshold or reference promotion')
        prior=read(s['feca_result']);prior_adm=read(s['feca_root_admission'])
        check(prior_adm['result']==s['feca_result'] and prior_adm['actual_review']==s['feca_actual_review'] and prior['status']=='consumed_trimer_regression_complete' and prior['error'] is None,'Admitted feca comparator')
        check(len(prior['parameters'])==12 and all(selected['parameters'][i]!=0 for i in (12,13)),'No historical12parameter six-row import authorized')
        for cid in m['cycles']:
            a=m['cycles'][cid];b=prior['cycles'][cid]
            check(all(a[k]==b[k] for k in ('model_terms','reference_electronic_E3_kcal_mol','reference_mu3_Debye')),'Exact retained reference cycle '+cid)
        contrasts={f'aquo{n}__{pair}':{'La_cycle':f'aquo{n}_central_La__{pair}','Dy_cycle':f'aquo{n}_central_Dy__{pair}','sign':'Dy-minus-La'} for n in (8,9) for pair in ('near','far')}
        metrics=prod['selected_metrics']
        retained={k:copy.deepcopy(metrics[k]) for k in ('component_RMS','component_errors_kcal_mol','q','J_P','J_B','signed_W9_z_safeguard','electrical_residual_RMS','all25_original_physical_gates','all_declared_safeguards_pass')}
        retained['balanced_W9_exchange']=copy.deepcopy(metrics['original25']['energy/balanced_exchange/aquo9_central'])
        payload={'status':'prepared_radial14_terminal_trimer_metadata_only','implementation':pin(__file__),'source_inputs':result['inputs'],'sources':s,'root_declaration':s['root_declaration'],'parent_manifest':s['parent_manifest'],
                 'candidate':{'artifact':s['terminal_candidate'],'key':dec['candidate_key'],'parameters':selected['parameters'],'fit_result':s['terminal_result'],'fit_configuration':s['fit_configuration'],'root_admission':s['terminal_root_admission'],'comparison':s['terminal_comparison']},
                 'rows':m['rows'],'new_analyze_order':ids,'cycles':m['cycles'],'contrasts':contrasts,'reused_selected':reused,'selected46_ids':sorted(selected['cells']),
                 'static_terms':{k:em['static'][k] for k in static_ids},'historical_comparators':{'feca_result':s['feca_result'],'feca_actual_review':s['feca_actual_review'],'feca_root_admission':s['feca_root_admission'],'baseline_original_comparison':s['original_trimer_comparison'],'prior_regressions':s['prior_regressions']},
                 'reference':{'historical_reference':m['historical_reference'],'original_comparison':s['original_trimer_comparison'],'original_native_admission':s['original_trimer_native_admission'],'quantum_moment_rows':{k:original['quantum_moment_rows'][k] for k in sorted(refs)},'reference_native_metadata':original['reference_native_metadata'],'limitations':original['limitations']},
                 'terminal_metrics':retained,'units':m['units'],'counts':{'cycles':8,'paired_contrasts':4,'unique_model_terms':42,'reused_selected':14,'reused_static':22,'maximum_new_analyze':6,'new_testgrad':0,'new_reference':0,'new_solver':0,'new_probes':0},'execution':dec['execution'],'literal_radii_policy':'Use deposited XYZ bytes; do not clamp metal-O radii to nominal2.50A or regenerate geometry.','historical_six_row_import_authorized':False,'scalar_calculation_performed':False,'execution_release':None,'physical_accuracy_qualified':False}
        save(out/'INVENTORY.json',payload)
        result.update(status='prepared_radial14_terminal_trimer_inventory',inventory=pin(out/'INVENTORY.json'),counts=payload['counts'],candidate=s['terminal_candidate'])
    except BaseException as exc:result.update(status='failed_metadata_inventory',error=repr(exc));raise
    finally:
        result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'hostname':socket.gethostname(),'PID':os.getpid(),'job_id':None};save(out/'RESULT.json',result)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=inventory(a.out);print(json.dumps({'status':r['status'],'inventory':r.get('inventory'),'checks':len(r['checks']),'cost':r['cost']}))
