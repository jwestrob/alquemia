#!/usr/bin/env python3
"""Finite saved receipt/flag/pin audit. No producer imports or numerical replay."""
import datetime
import hashlib
import json
import math
from pathlib import Path
import time

ROOT=Path.cwd().resolve()
REL=Path('explicit_model_v1/radial_drive_shape_v1')
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/REL
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/REL
RUN=W/'qualification_v1/recovery_v1/run_v1'
OUT=W/'qualification_v1/recovery_v1/independent_actual_review_v1/actual_v1'
EXPECTED='fe8d006571885eac09d7368e9b0a94668a5b6427e0e1b6b5be0536d8694f6f08'

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    t,c=time.monotonic(),time.process_time();checks=[];pins={};saved_gates=[]
    def check(name,ok):
        if not ok:raise ValueError(name)
        checks.append({'name':name,'passed':True})
    def verify(p):
        if p['path'] not in pins:
            check('pin:'+p['path'],pin(p['path'])==p);pins[p['path']]=p
        else:check('same_pin:'+p['path'],pins[p['path']]==p)
        return Path(p['path'])
    def read(p):return json.loads(verify(p).read_text())
    def allpins(x):
        if isinstance(x,dict):
            if set(x)=={'path','sha256'}:verify(x)
            else:
                for v in x.values():allpins(v)
        elif isinstance(x,list):
            for v in x:allpins(v)
    def gates(x,path):
        if isinstance(x,dict):
            if 'passed' in x:
                check('saved_gate:'+path,x['passed'] is True)
                if 'error' in x and 'bound' in x:
                    check('saved_error_bound:'+path,all(math.isfinite(float(x[k])) for k in ('error','bound')) and 0<=x['error']<=x['bound'])
                saved_gates.append(path)
            for k,v in x.items():gates(v,path+'/'+k)
        elif isinstance(x,list):
            for i,v in enumerate(x):gates(v,path+'/'+str(i))
    rp=pin(RUN/'execution/RESULT.json');check('exact_terminal_pin',rp['sha256']==EXPECTED)
    r=read(rp);cfgp=pin(RUN/'CONFIG_RUN.json');cfg=read(cfgp);inp=read(cfg['inputs']);old=read(inp['parent_result']);pc=read(inp['parent_configuration']);m=read(r['manifest']);release=read(cfg['root_release']);sr=read(cfg['source_review'])
    allpins(cfg);allpins(inp)
    check('terminal_success',r['status']=='radial_drive_numerical_qualification_recovered' and r['error'] is None and r['numerical_checks_passed'] and not r['global_stop'])
    check('execution_configuration',r['configuration']==cfg and r['implementation']==cfg['implementation'] and r['inputs']==cfg['inputs'] and cfg['output']==str(RUN/'execution'))
    fp=hashlib.sha256(json.dumps({k:v for k,v in cfg.items() if k!='root_release'},sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    check('root_release',release['status']=='released_radial_field_recovery' and release['configuration_fingerprint']==fp and release['source_review']==cfg['source_review'] and release['parent_result']==r['parent_result']==inp['parent_result'] and release['limits']==cfg['limits'])
    check('source_review',sr['status']=='passed_radial_field_recovery_source_review' and all(sr[k]==cfg[k] for k in ('implementation','field_adapter','inputs')))
    check('exact_native_lineage',all(r[k]==pc[k] for k in ('manifest','source_manifest','build','executables')))
    build=read(r['build']);check('actual_build',build['status']=='complete' and build['manifest']==r['source_manifest'] and build['executables']==r['executables']);allpins(r['executables'])
    check('preserved_parent_failure',old['status']=='radial_drive_numerical_qualification_failed' and old['error']==r['preserved_parent_error']=="ValueError('Zero-field export required')")
    check('parent46_recovery3_total49',old['native_program_invocation_attempts']==old['native_programs_returned_zero']==46 and r['new_native_program_invocation_attempts']==r['new_native_programs_returned_zero']==r['new_native_programs_completed']==3 and r['native_program_invocation_attempts']==r['native_programs_returned_zero']==r['native_programs_completed']==49)
    check('saved4_and_FD132',r['saved_field_qualification_attempts']==r['saved_fields_qualified']==4 and all(r[k]==old[k]==132 for k in ('internal_FD_energies_requested','internal_FD_energies_observed','internal_FD_energies_completed')) and r['new_internal_FD_energies']==0)
    check('inventory',set(r['rows'])==set(old['rows'])=={j['id'] for j in m['programs']} and len(r['rows'])==49)
    stage_counts={};raw_receipts=[];new_receipts=[];zero=[];fd=[];export_count=0
    for job in m['programs']:
        key=job['id'];row=r['rows'][key];stage=job['stage'];stage_counts[stage]=stage_counts.get(stage,0)+1
        check('complete:'+key,row['status']=='complete' and row['job']==job)
        allpins(row)
        check('literal_inputs:'+key,all(row['runtime_input'][k]['sha256']==job['input'][k]['sha256'] for k in ('xyz','key','parameters')))
        check('single_native_receipt:'+key,len(row['receipts'])==1)
        rec=row['receipts'][0];args=['system.xyz','Y','Y','0.00001'] if job['fd'] else ['system.xyz','Y','N'] if job['program']=='testgrad' else ['system.xyz','E' if stage=='directions' else 'EM']
        check('native_receipt:'+key,rec['argv']==[r['executables'][job['program']]['path']]+args and rec['returncode']==0 and rec['timeout_seconds']==120 and rec['output']==row['actual_output'] and all(v=='1' for v in rec['threads'].values()))
        raw_receipts.append((rec['cwd'],tuple(rec['argv'])))
        check('native_output_location:'+key,Path(row['actual_output']['path']).parent==Path(rec['cwd']))
        if stage not in ('fields','bypass'):check('unchanged42:'+key,row==old['rows'][key])
        elif stage=='fields':
            check('saved_field_execution_unchanged:'+key,all(row[k]==old['rows'][key][k] for k in ('job','receipts','runtime_input','actual_output','exports','induction','response_readback','moments')))
            check('new_field_qualification:'+key,'matrix_qualification' not in old['rows'][key] and row['matrix_qualification']['passed'])
            for ex in row['matrix_qualification']['exports']:
                h=ex['external_field'];check('saved_field_header:'+key,h['native_external_tokens'][:2]==['EXTERNAL','T'] and h['maximum_error']<=h['arithmetic_bound'])
        else:
            check('only_new_bypass:'+key,old['rows'][key]=={'status':'not_attempted'} and Path(rec['cwd'])==RUN/'execution/cells'/key)
            check('per_cell_receipt:'+key,json.loads((Path(rec['cwd'])/'RESULT.json').read_text())==row)
            check('bypass_exact_energy:'+key,row['zero_identity']['exact_saved_bypass_energy'])
            text=verify(row['actual_output']).read_text()
            check('raw_two_exports:'+key,text.count('Response Matrix Written To :')==2 and 'Total Potential Energy' in text and 'Induced Dipoles Written To' in text)
            new_receipts.append(rec)
        check('induction_gate:'+key,0<=row['induction']['maximum_printed_RMS_residual_Debye']<=row['induction']['tolerance_Debye']==1e-8)
        if 'matrix_qualification' in row:
            check('matrix_pass:'+key,row['matrix_qualification']['passed'] and len(row['exports'])==len(row['matrix_qualification']['exports'])==2);export_count+=2
        if job['fd']:
            check('Cartesian_FD_pass:'+key,row['value']['fd_passed'] and row['response_readback']['native_solve_blocks']==1+6*m['rows'][job['cell']]['atoms'])
            fd.append({'id':key,'maximum_component_error_kcal_mol_A':row['value']['max_component_error_kcal_mol_A']})
        if stage in ('absent','zero'):
            z=row['zero_identity'];check('literal_zero_parity:'+key,z['literal_native_tokens_equal'] and (z.get('exact_raw_E_F') or (z.get('exact_raw_E_M_uv') and z['literal_site_tokens_equal'])))
            zero.append({'id':key,'zero_identity':z,'actual_output':row['actual_output']})
        gates(row,'rows/'+key)
    check('unique_native_attempts',len(set(raw_receipts))==49)
    check('finite_program_stages',stage_counts=={'absent':4,'zero':10,'active':12,'directions':16,'fields':4,'bypass':3} and export_count==42 and len(fd)==4)
    check('serial_new_calls',all(datetime.datetime.fromisoformat(a['ended_utc'])<=datetime.datetime.fromisoformat(b['started_utc']) for a,b in zip(new_receipts,new_receipts[1:])))
    check('old_operator_and_direction_checks_exact',r['operator_checks']==old['operator_checks'] and r['directional_checks']==old['directional_checks'] and len(r['operator_checks'])==6 and len(r['directional_checks'])==4)
    check('two_field_closures',set(r['field_checks'])=={'water_central_La','water_central_Dy'} and all(v['E_solve_bound_separate_from_M'] for v in r['field_checks'].values()))
    for key in ('operator_checks','directional_checks','field_checks'):gates(r[key],key)
    ap=pin(RUN/'SACCT_SUMMARY.json');acct=read(ap);raw=pin(RUN/'SACCT.psv')
    check('scheduler20413',acct['job_id']==r['job_id']=='20413' and acct['state']=='COMPLETED' and acct['exit']=='0:0' and acct['elapsed_seconds']==9 and acct['allocated_CPUs']==12 and acct['allocated_CPU_seconds']==108 and acct['raw_sha256']==raw['sha256'])
    check('cost_separation',r['parent_scheduler_cost']['allocated_CPU_seconds']==288 and r['cost']['parent_cost_not_added'] and r['cost']['per_receipt_child_deltas_not_summed'])
    bp=pin(D/'zero_cache_bridge_assessment_v1/ASSESSMENT.json');check('frozen_conditional_assessment',bp['sha256']=='b883c4ab0daaa4d0a54ce90cb08fab02381543af2ae1ce5f38ba2391eb86e2ad');b=read(bp)
    for key in ('old_center','old24_result','old24_actual_root_admission','old24_independent_actual_review','old24_inventory','old_domain','radial_domain','old_source_manifest','new_source_manifest','source_inverse_and_ancillary_review','old_build','new_build'):verify(b[key])
    check('bridge_actual_source_binary_join',b['new_source_manifest']==r['source_manifest'] and b['new_build']==r['build'] and b['new_executables']==r['executables'] and b['preserved_failed_qualification']==r['parent_result'])
    check('same_zero_receipts',b['actual_absent_zero_parity_rows']==zero)
    check('exact_old24_order_and_zero_view',len(b['probes_in_declared_order'])==24 and [x['ordinal'] for x in b['probes_in_declared_order']]==list(range(1,25)) and all(x['extended_parameters_view']==x['old_parameters']+[0.,0.] for x in b['probes_in_declared_order']))
    check('actual_Dy_v_span',len(b['actual_normalized_spans'])==12 and b['actual_normalized_spans'][7]['parameter']=='Dy_v' and b['actual_normalized_spans'][7]['actual_normalized_span']==0.005170593745372654 and b['actual_normalized_spans'][7]['clamped_minus'])
    bridge={'status':'ready_for_root_zero_cache_equivalence_approval','proposed_root_status':b['bridge_contract']['status'],'root_approval':None,'cache_bridge_admitted':False,'qualification_result':rp,'conditional_assessment':bp,**{k:b[k] for k in ('old_center','old_center_key','old24_result','old24_actual_root_admission','old24_independent_actual_review','old24_inventory','old_domain','radial_domain','old_source_manifest','new_source_manifest','source_inverse_and_ancillary_review','old_build','new_build','old_executables','new_executables','probes_in_declared_order','actual_normalized_spans')},'rules':b['bridge_contract'],'claim_limit':b['claim_limit']}
    review={'status':'passed_radial_drive_qualification_actual_review','reviewer_implementation':pin(__file__),'result':rp,'configuration':cfgp,'implementation':r['implementation'],'source_review':cfg['source_review'],'root_release':cfg['root_release'],'parent_result':r['parent_result'],'numerical_checks_passed':True,'zero_parity_passed':True,'selected_physical_accuracy_qualified':False,'cache_bridge':bridge,'counts':{'native_programs':49,'inherited_native_programs':46,'new_native_programs':3,'saved_field_requalifications':4,'rows_reused_exactly':42,'internal_FD_energies':132,'new_internal_FD_energies':0,'matrix_exports':42,'absent_or_zero_parity_programs':14,'saved_gate_records':len(saved_gates)},'stage_counts':stage_counts,'Cartesian_FD':fd,'field_checks':r['field_checks'],'directional_checks':r['directional_checks'],'checks':checks,'unique_pins':len(pins),'scheduler_accounting':{'parent':inp['parent_sacct'],'recovery':ap,'recovery_raw':raw,'parent_allocated_CPU_seconds':288,'recovery_allocated_CPU_seconds':108,'combined_allocated_CPU_seconds':396},'application_cost':r['cost'],'limitations':['No scientific numerical replay: admitted producer functions/arrays were not re-solved; saved flags/error bounds, raw pins and receipt lineage were checked.','All four declared directional slopes are zero at these fixtures; these are zero-slope consistency checks. The four Cartesian panels supply nonzero-component force checks.','Field z odd moments are zero at these symmetric dimers; curvature and native/dense energy closures pass their original bounds.','Cache reuse requires the separate explicit root bridge. Old24 artifact bodies and raw native calculations were not re-admitted; their frozen prior assessment and admissions are reused.','Source equivalence plus representative zero tests does not claim universal binary bitwise equivalence, physical model accuracy or permission for fitting/new comparisons.'],'new_native_calls':0,'new_solver_calls':0,'new_dense_replays':0,'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-c,'numerical_threads':1}}
    (OUT/'REVIEW.json').write_text(json.dumps(review,indent=2,allow_nan=False)+'\n')
    (OUT/'REPORT.txt').write_text('PASS: finite saved radial qualification review.\n\n49 programs/132 internal FD energies are complete:46 original executions plus3 serial bypasses;42 completed row objects unchanged;4 saved field outputs newly qualified without native repeats.14 absent/zero programs retain exact native E/F or E/M/dual/site token parity. All original saved Cartesian, directional, matrix, dual, field and bypass gates pass.\n\nMaximum Cartesian component error is about2e-8 kcal/mol/Angstrom. The four declared directional slopes and field-z odd moments are zero on these fixtures; this limits their coverage. No numerical matrix/field/FD replay was performed.\n\nOriginal failed20380 cost288 allocated CPU-s; recovery20413 cost108; combined396. Recovery application7.897wall,5.784controller CPU and1.179child CPU seconds.\n\nAn explicit root zero-cache bridge is now supportable with old producer/raw provenance intact, centerfeca, exact24 derivative order and saved actual spans. Dy_v uses0.005170593745372654 normalized span. Root approval remains null; no physical accuracy, fit, or112-call comparison is admitted here.\n')
    print(json.dumps({'review':pin(OUT/'REVIEW.json'),'report':pin(OUT/'REPORT.txt'),'checks':len(checks),'saved_gates':len(saved_gates),'cost':review['cost']}))

if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=False)
    try:main()
    except Exception as e:
        (OUT/'FAILURE.json').write_text(json.dumps({'error':repr(e),'new_native_calls':0,'new_dense_replays':0})+'\n');raise
