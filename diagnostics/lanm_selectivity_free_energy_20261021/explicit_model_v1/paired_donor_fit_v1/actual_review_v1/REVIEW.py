#!/usr/bin/env python3
"""Terminal-only selected28 replay. No native calls, fitting, or intermediate-vector audit."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import time

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
CONFIG_SHA='51f93f2b5475d16a6d6ed21a6a19e5f6e9ceee52eb1ba63941a1590fa0c9a93c'
EXPECTED_EVALUATOR='c0d3085cce6e4b7d8e74687ee0729922fe57836426e53731e2b526dc4790592a'
EXPECTED_OBJECTIVE='86e85fcb0537bdf5b459874eeda0c8bb68aead9758c712014b18b2db94757974'

def require(ok,message):
    if not ok:raise ValueError(message)

def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):
    require(pin(p['path'])==p,'Changed evidence: '+p['path']);return Path(p['path'])

def read(p):return json.loads(verify(p).read_text())

def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def save(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def terminal_scheduler(run,out):
    """Use a terminal watcher snapshot or exactly one read-only accounting query."""
    observed=[];receipt_path=run/'completion_event.json'
    terminal={'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','BOOT_FAIL','PREEMPTED','DEADLINE','REVOKED'}
    if receipt_path.exists():
        raw=receipt_path.read_bytes();snapshot=out/'COMPLETION_SNAPSHOT.json';snapshot.write_bytes(raw)
        observed.append({'kind':'watcher_snapshot','original_path':str(receipt_path),'snapshot':pin(snapshot)})
        try:
            record=json.loads(raw)['records']['17621']
            if record['state'].split()[0].rstrip('+') in terminal:return record,{'source':'watcher_snapshot','evidence':observed}
        except (ValueError,KeyError,TypeError,IndexError):pass
    argv=['sacct','--allocations','--jobs=17621','--noheader','--parsable2','--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS']
    before=time.monotonic();query={'argv':argv,'query_count':1,'read_only':True,'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    try:
        result=subprocess.run(argv,text=True,capture_output=True,timeout=30,check=False)
        query.update(returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
    except Exception as exc:
        query.update(error=repr(exc),returncode=None,stdout='',stderr='')
    query['elapsed_seconds']=time.monotonic()-before;qp=out/'SACCT_QUERY.json';save(qp,query)
    observed.append({'kind':'single_sacct_query','receipt':pin(qp)})
    require(query['returncode']==0,'Terminal scheduler evidence unavailable: one sacct query failed')
    rows=[line.strip().split('|') for line in query['stdout'].splitlines() if line.strip()]
    rows=[row for row in rows if row[0]=='17621'];require(len(rows)==1 and len(rows[0])>=5,'Terminal scheduler evidence unavailable: exact job row missing')
    job,state,exit_code,elapsed,cpus=rows[0][:5]
    require(state.split()[0].rstrip('+') in terminal,'Terminal scheduler evidence unavailable: job is not terminal')
    return {'job_id':job,'state':state,'exit_code':exit_code,'elapsed_seconds':elapsed,'allocated_cpus':cpus},{'source':'single_sacct_query','evidence':observed}

def controls(score):
    energies={group+'/'+name:{'error_kcal_mol':value,'passed':abs(value)<=1} for group,values in score['energy_errors_kcal_mol'].items() for name,value in values.items()}
    forces={kind+'/'+q:{'rigid_RMS_kcal_mol_A':value['rigid_RMS_kcal_mol_A'],'maximum_body_component_error_kcal_mol_A':value['maximum_body_component_error_kcal_mol_A'],'passed':value['passed']} for kind,values in score['per_fixture'].items() for q,value in values.items()}
    require(len(energies)==33 and len(forces)==36,'Declared control denominator differs')
    return {'energy':energies,'force':forces,'energy_failures':sum(not v['passed'] for v in energies.values()),'force_failures':sum(not v['passed'] for v in forces.values())}

def main(run,out):
    run=Path(run).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    report={'status':'unavailable','implementation':pin(__file__),'checks':{},'new_native_calls':0,'local_threads':1,'error':None}
    def check(label,yes):report['checks'][label]=bool(yes);require(yes,label)
    try:
        check('one_thread',all(os.environ.get(k)=='1' for k in THREADS))
        cp=pin(run/'CONFIG_RUN.json');check('released_config',cp['sha256']==CONFIG_SHA);cfg=read(cp)
        scheduler,completion=terminal_scheduler(run,out);report['scheduler_evidence']=completion
        check('scheduler_terminal',scheduler['state']=='COMPLETED' and scheduler['exit_code']=='0:0' and int(scheduler['allocated_cpus'])==8)
        rp=pin(run/'execution/RESULT.json');actual=read(rp)
        check('terminal_application',actual['status']=='bounded_paired_donor_complete' and actual['error'] is None and actual['candidate_frozen'] and actual['job_id']=='17621')
        check('configuration_join',actual['config']==cp and actual['implementation']==cfg['implementation'])
        for k in ('implementation','evaluator','evaluator_config','objective','body_helper','training_result','water_result','transfer_result','root_release','preparation_result','specification'):verify(cfg[k])
        check('frozen_sources',cfg['evaluator']['sha256']==EXPECTED_EVALUATOR and cfg['objective']['sha256']==EXPECTED_OBJECTIVE)
        evmod=module(cfg['evaluator'],'review_donor_evaluator');ecfg=read(cfg['evaluator_config']);data=read(cfg['training_result']);transfer=read(cfg['transfer_result']);water=read(cfg['water_result']);prep=read(cfg['preparation_result'])
        selected=actual['selected'];model=selected['model'];check('complete_selected28',selected['status']=='complete' and set(model)==set(evmod.CELLS) and all(v=='complete' for v in selected['cells'].values()))
        check('original_family',actual['family']==prep['family'] and len(actual['family']['parameters'])==8)
        check('initial_score_unchanged',actual['initial']['parameters']==prep['initial']['parameters'] and actual['initial']['score']==prep['initial']['score'])
        # Build only the saved renderer context; no evaluator constructor/cache/native entry.
        ev=evmod.Evaluator.__new__(evmod.Evaluator);ev.config=ecfg;ev.data=data;ev.family=actual['family'];ev.resolver=module(ecfg['resolver'],'review_donor_resolver');ev._base_text={};ev._base_parameters={}
        native=module(ecfg['native_helper'],'review_native_parser')
        def forbidden(*a,**k):raise AssertionError('Native calls forbidden in saved review')
        native.invoke=forbidden
        for gid in evmod.CELLS:
            files=ecfg['baseline_inputs'][gid];prm=verify(files['parameters']).read_text();key=verify(files['key']).read_text()
            ev._base_text[gid]=(prm,key);ev._base_parameters[gid]=ev.resolver.parse_parameters(prm,key)
        ev.checked_values(selected['parameters']);sources={};baseline_count=0
        for gid,value in model.items():
            row=data['rows'][gid];pv=value['provenance'];token,identity,key,render=ev.cache_identity(selected['parameters'],gid)
            check('request_'+gid,pv['requested_vector']==selected['parameters'] and pv['cache_key']==token and pv['effective_parameters']==render['effective_values'] and pv['render']==render and pv['request_cache_implementation']==cfg['evaluator'])
            check('source_'+gid,pv['source_input']==ecfg['source_inputs'][gid])
            inp=pv['effective_input']
            for p in inp.values():verify(p)
            check('physical_input_'+gid,inp['key']['sha256']==hashlib.sha256(key.encode()).hexdigest() and all(inp[k]['sha256']==row['input'][k]['sha256'] for k in ('xyz','parameters')))
            authority=read(pv['native_authority'])
            if pv['reused_baseline']:
                baseline_count+=1;check('baseline_authority_'+gid,pv['native_authority']==cfg['transfer_result'] and render['baseline_key_reused'] and pv['native_implementation'] is None and pv['native_cache_key'] is None)
                check('baseline_values_'+gid,all(value[k]==transfer['rows'][gid]['model'][k] for k in ('energy','forces','native')))
            else:
                check('new_authority_'+gid,Path(pv['native_authority']['path']).is_relative_to(run/'execution/cache/cells') and authority['status']=='complete' and authority['cell_id']==gid and authority['job_id']=='17621')
                check('new_identity_'+gid,authority['cache_key']==pv['native_cache_key']==token and authority['cache_identity']==identity and pv['native_implementation']==cfg['evaluator'])
                check('new_values_'+gid,all(value[k]==authority[k] for k in ('energy','forces','native')) and authority['input']==inp and authority['render']==render)
                check('new_receipts_'+gid,len(authority['receipts'])==2 and all(r['returncode']==0 for r in authority['receipts']))
                for r in authority['receipts']:verify(r['output'])
            n=value['native'];replay=native.parse(verify(n['analyze_output']),verify(n['gradient_output']),dict(row['input'],derivative_test=n['finite_difference'] is not None))
            check('raw_'+gid,replay==n and value['energy']==n['total_energy_kcal_mol'] and value['forces']==n['force_kcal_mol_A'])
            sources[gid]={'authority':pv['native_authority'],'effective_input':inp,'reused_baseline':pv['reused_baseline'],'native_outputs':{k:n[k] for k in ('analyze_output','gradient_output')}}
        frozen=actual['frozen_water_invariants'];check('frozen_water_evidence',frozen==prep['frozen_water_invariants'] and frozen['passed'] and set(frozen['records'])==set(evmod.WATER_CELLS))
        for gid,record in frozen['records'].items():
            original=water['selected']['model'][gid]
            check('water_selected_'+gid,all(record[k]==original[k] for k in ('energy','forces','native')) and record['native_authority']==original['provenance']['native_authority'] and record['unchanged'])
            verify(record['native_authority']);n=record['native'];verify(n['analyze_output']);verify(n['gradient_output'])
            if gid in model:check('water_invariant_'+gid,all(model[gid][k]==original[k] for k in ('energy','forces','native')))
        water_baseline={g:transfer['rows'][g]['model'] for g in evmod.WATER_CENTERS};obj=module(cfg['objective'],'review_donor_objective');body=module(cfg['body_helper'],'review_donor_body')
        score=obj.evaluate(data,model,body,water_baseline);check('objective_exact_replay',score==selected['score'])
        check('declared_pass_status',actual['development_checks_passed']==score['development_checks_passed'])
        check('frozen_water_controls',score['frozen_water_centers']['passed'] and water['selected']['score']['development_checks_passed'])
        ledger=read(actual['evaluator']['ledger']);calls=actual['new_native_program_invocations']
        check('finite_native_budget',calls==actual['evaluator']['native_program_invocations']==ledger['native_program_invocations'] and calls%2==0 and 0<=calls<=9216 and not ledger['stopped'])
        check('finite_vector_budget',len(actual['vectors'])==actual['new_parameter_vectors']+1==len(ledger['vector_keys']) and len(actual['vectors'])<=481)
        check('finite_step_budget',actual['completed_outer_steps']==len(actual['iterations'])<=24)
        check('parent_preservation',actual['original_results_unchanged'])
        for p in actual['protected_parents']:verify(p)
        bounds=[]
        for v,p in zip(selected['parameters'],actual['family']['parameters']):
            span=p['upper']-p['lower'];bounds.append(dict(p,selected=v,normalized=(v-p['lower'])/span,at_lower=v==p['lower'],at_upper=v==p['upper'],normalized_nearest_bound_distance=min(v-p['lower'],p['upper']-v)/span))
        last=None
        if actual['iterations']:
            entry=actual['iterations'][-1];origin=actual['vectors'][entry['origin_key']]
            check('last_jacobian_origin_key',evmod.digest(origin['parameters'])==entry['origin_key'] and entry['status']=='complete')
            check('last_jacobian_dimensions',len(entry['Jacobian_normalized'])==len(origin['score']['weighted_residuals']) and all(len(row)==8 for row in entry['Jacobian_normalized']) and len(entry['derivative_keys'])==16)
            last={'outer_step':entry['outer_step'],'origin_key':entry['origin_key'],'origin_parameters':origin['parameters'],'selected_parameters':selected['parameters'],
                  'origin_is_selected':origin['parameters']==selected['parameters'],'origin_objective':origin['score']['objective'],
                  'rank_at_origin':entry['linear_directions']['rank'],'singular_values_at_origin':entry['linear_directions']['singular_values'],
                  'condition_at_origin':entry['origin_condition'],'column_norms_at_origin':entry['Jacobian_column_norms'],
                  'actual_FD_widths':entry['derivative_actual_widths'],'actual_improvement_in_last_batch':entry['actual_improvement']}
        initial_controls=controls(actual['initial']['score']);selected_controls=controls(score)
        cost={k:actual[k] for k in ('elapsed_seconds','parent_process_CPU_seconds','aggregate_native_child_user_CPU_seconds','aggregate_native_child_system_CPU_seconds')}
        cost.update(scheduler_elapsed_seconds=int(scheduler['elapsed_seconds']),allocated_cpus=int(scheduler['allocated_cpus']),allocated_cpu_seconds=int(scheduler['elapsed_seconds'])*int(scheduler['allocated_cpus']),new_native_programs=calls,new_vectors=actual['new_parameter_vectors'])
        report.update(status='passed',result=rp,configuration=cp,completion_receipt=completion,selected_sources=sources,selected_baseline_reused_cells=baseline_count,
            selected_new_native_cells=28-baseline_count,initial_objective=actual['initial']['score']['objective'],selected_objective=score['objective'],
            initial_controls=initial_controls,selected_controls=selected_controls,parameter_bounds=bounds,last_Jacobian=last,
            stopping_reason=actual['stopping_reason'],development_checks_passed=score['development_checks_passed'],completed_outer_steps=actual['completed_outer_steps'],
            frozen_water_unchanged=True,native_cost=cost,reference_eligibility=cfg['reference_eligibility'],physical_accuracy_qualified=False,reserved_validation_opened=False)
    except Exception as exc:report.update(status='unavailable' if 'Terminal scheduler evidence unavailable' in str(exc) else 'failed',error=repr(exc))
    finally:
        report.update(check_count=len(report['checks']),review_wall_seconds=time.monotonic()-start,review_CPU_seconds=time.process_time()-cpu)
        save(out/'REVIEW.json',report)
    if report['status']=='passed':
        c=report['selected_controls'];old=report['initial_controls'];last=report['last_Jacobian'];cost=report['native_cost']
        lines=['# Donor development: actual17621 selected-point review','',
            f"Saved native/input/objective replay **passed** ({report['check_count']} checks; zero new molecular calls). The selected consumed-data model {'passes' if report['development_checks_passed'] else 'does not pass'} the declared development controls.",'',
            f"Objective {report['initial_objective']:.9g} → {report['selected_objective']:.9g}; stop `{report['stopping_reason']}` after {report['completed_outer_steps']} steps.",
            f"Energy failures {old['energy_failures']}/33 → {c['energy_failures']}/33; rigid-force failures {old['force_failures']}/36 → {c['force_failures']}/36. All eight water records and their qualification remain unchanged.",'',
            '## Parameter bounds','', '| Parameter | Selected | Normalized position | At bound |','|---|---:|---:|---|']
        for p in report['parameter_bounds']:lines.append(f"| {p['name']} | {p['selected']:.9g} | {p['normalized']:.6f} | {'yes' if p['at_lower'] or p['at_upper'] else 'no'} |")
        if last:lines.extend(['',f"Last Jacobian rank {last['rank_at_origin']}/8, condition {last['condition_at_origin']}; its origin {'is' if last['origin_is_selected'] else 'is not'} the selected point. These diagnostics describe that local origin."])
        lines.extend(['','## Energy controls','', '| Control | Error, kcal/mol | Pass |','|---|---:|---|'])
        for name,v in c['energy'].items():lines.append(f"| {name} | {v['error_kcal_mol']:.6f} | {v['passed']} |")
        lines.extend(['','## Rigid-force controls','', '| Fixture | RMS | Max component | Pass |','|---|---:|---:|---|'])
        for name,v in c['force'].items():lines.append(f"| {name} | {v['rigid_RMS_kcal_mol_A']:.6f} | {v['maximum_body_component_error_kcal_mol_A']:.6f} | {v['passed']} |")
        lines.extend(['','## Cost and interpretation','',f"{cost['new_native_programs']} native programs; {cost['scheduler_elapsed_seconds']} scheduler seconds × {cost['allocated_cpus']} CPUs = {cost['allocated_cpu_seconds']} allocated CPU-seconds. Measured parent CPU {cost['parent_process_CPU_seconds']:.3f}s; aggregate native child CPU {cost['aggregate_native_child_user_CPU_seconds']+cost['aggregate_native_child_system_CPU_seconds']:.3f}s. Overlapping per-invocation CPU receipts are not summed.",'',
            'These are consumed-geometry representability results against provisional uncorrected frozen-f ECP-PBE0-D4 targets. They do not establish reference accuracy, protein transfer, bulk hydration or exchange free energies. No reserved configurations were opened. A bounded or local stop does not establish global family impossibility.','',
            f"Actual result SHA `{report['result']['sha256']}`; replay receipt SHA `{pin(out/'REVIEW.json')['sha256']}`."])
        (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:report[k] for k in ('status','error','check_count','review_wall_seconds','review_CPU_seconds')}))
    return report['status']=='passed'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',required=True);p.add_argument('--output',required=True);a=p.parse_args();raise SystemExit(0 if main(a.run,a.output) else 1)
