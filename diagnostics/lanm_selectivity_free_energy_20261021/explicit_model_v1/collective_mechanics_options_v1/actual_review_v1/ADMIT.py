#!/usr/bin/env python3
"""Saved actual application admission; no native invocation or scientific cycles."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('actual_trimer_executor',HERE.parent/'RUN.py')
r=importlib.util.module_from_spec(s);s.loader.exec_module(r)


def inputs(config):
    require=r.require; require(config['implementation']==r.pin(__file__),'Admission source changed')
    cfg=r.read(config['execution_config']); release=r.read(config['root_release']);m=r.read(cfg['manifest'])
    require(release['configuration']==config['execution_config'],'Root configuration differs')
    r.release(cfg,m,release)
    return cfg,release,m


def model(config,out):
    t=time.monotonic();cpu=time.process_time();cfg,release,m=inputs(config)
    path=Path(release['model_execution_root'])/'RESULT.json';result=json.loads(path.read_text())
    r.require(result['release']==config['root_release'],'Actual model release differs')
    admitted=r.admit_model(cfg,m,result);rows={}
    for cid,x in admitted.items():
        q=x['parsed_moments'];row=m['model_rows'][cid]
        r.require(len(q['induced_Debye'])==row['atoms'] and q['induction']['solutions']==2,'Actual model moment denominator differs')
        rows[cid]={'energy_kcal_mol':q['energy_kcal_mol'],'components':q['components'],
            'induced_sum_Debye':q['induced_sum_Debye'],'total_lab_Debye':q['total_lab_Debye'],
            'induction':q['induction'],'analyze_output':q['analyze_output'],'induced_output':q['induced_output'],
            'runtime_input':x['runtime_input'],'start':x['start'],'atoms':row['atoms'],'charge_e':q['charge_e']}
    summary={'status':'passed_saved_model_application_admission','implementation':r.pin(__file__),
        'execution_config':config['execution_config'],'root_release':config['root_release'],'manifest':cfg['manifest'],
        'model_result':r.pin(path),'job_id':result['job_id'],'complete':len(rows),'declared':16,
        'exact_raw_energy_and_moment_replays':16,'rows':rows,'native_program_invocations':result['native_program_invocations'],
        'maximum_induction_residual_Debye':max(x['induction']['maximum_printed_RMS_residual_Debye'] for x in rows.values()),
        'actual_process_cost':result['cost'],'scheduler_cost':None,'new_native_calls':0,
        'admission_cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'local_threads':1},
        'independent_scientific_comparison_performed':False}
    out.mkdir(parents=True,exist_ok=False);r.save(out/'MODEL_ADMISSION.json',summary)
    return r.pin(out/'MODEL_ADMISSION.json')


def all_results(config,model_admission,out):
    t=time.monotonic();cpu=time.process_time();cfg,release,m=inputs(config)
    old=r.read(model_admission)
    r.require(old['status']=='passed_saved_model_application_admission' and old['implementation']==r.pin(__file__) and
              old['execution_config']==config['execution_config'] and old['root_release']==config['root_release'] and old['manifest']==cfg['manifest'],
              'Previous model admission differs')
    r.verify(old['model_result'])
    for row in old['rows'].values():
        for k in ('analyze_output','induced_output','start'):r.verify(row[k])
        for pin in row['runtime_input'].values():r.verify(pin)
    rows={};failures={}
    for c in m['reference_cells']:
        cid=c['cell_id'];path=Path(release['reference_execution_root'])/cid/'RESULT.json'
        try:
            x=r.admit_reference(cfg,m,path);r.require(x['release']==config['root_release'],'Actual reference release differs')
            q=x['parsed'];r.require(q['D4_energy_hartree']==0 and math.isfinite(q['electronic_energy_hartree']),'Electronic fragment nonfinite/changed')
            failed=[v for v in q.get('printed_scf_criteria',[]) if not v['passes']]
            rows[cid]={'status':'admitted','result':r.pin(path),'output':x['output'],'input':x['input'],'start':x['start'],
                'property_files':x['property_files'],'job_id':x['job_id'],'role':c['role'],'parent_geometry_id':c['parent_geometry_id'],
                'real_parent_atom_indices':c['real_parent_atom_indices'],'electronic_energy_hartree':q['electronic_energy_hartree'],
                'D4_energy_hartree':q['D4_energy_hartree'],'checks':q['checks'],
                'printed_failed_scf_criteria':failed,'native_program_invocations':x['native_program_invocations'],
                'ranks':x['ranks'],'elapsed_seconds':x['elapsed_seconds'],
                'allocated_elapsed_CPU_seconds_proxy':x['allocated_elapsed_CPU_seconds_proxy']}
        except Exception as exc:
            rows[cid]={'status':'unavailable_or_failed','result':r.pin(path) if path.exists() else None,'error':repr(exc)}
            failures[cid]=repr(exc)
    complete=[x for x in rows.values() if x['status']=='admitted']
    status='passed_saved_trimer_application_admission' if len(complete)==16 else 'failed_or_incomplete_saved_trimer_application_admission'
    result={'status':status,'implementation':r.pin(__file__),'execution_config':config['execution_config'],
        'root_release':config['root_release'],'manifest':cfg['manifest'],'model_admission':model_admission,
        'reference_declared':16,'reference_complete':len(complete),'model_declared':16,'model_complete':old['complete'],
        'reference_rows':rows,'failures':failures,'model_maximum_induction_residual_Debye':old['maximum_induction_residual_Debye'],
        'reference_native_program_invocations_admitted':sum(x['native_program_invocations'] for x in complete),
        'model_native_program_invocations':old['native_program_invocations'],
        'reference_allocated_elapsed_CPU_seconds_proxy':math.fsum(x['allocated_elapsed_CPU_seconds_proxy'] for x in complete),
        'model_process_cost':old['actual_process_cost'],'scheduler_cost':None,
        'reference_rows_retaining_failed_printed_scf_criteria':sum(bool(x['printed_failed_scf_criteria']) for x in complete),
        'new_native_calls':0,'quantum_dipole_property_arithmetic_owned_by_root_collector':True,
        'no_C_P_B_or_E3_scientific_comparison_performed':True,
        'admission_cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'local_threads':1}}
    out.mkdir(parents=True,exist_ok=False);r.save(out/'RESULT.json',result)
    lines=['# Trimer actual application admission','',
        f'Status: **{status}**. Reference {len(complete)}/16; model {old["complete"]}/16. Zero new native calls.',
        '', 'All admitted references replay exact prepared masks, states, orbital/auxiliary/ECP inventories and electronic energies. '
        'The original endpoint.property.txt pins are verified; numerical dipole extraction and eight E3/mu3 cycles belong to the separate root collector.',
        '',f'Model16 raw energies/moments and unchanged runtime inputs replay exactly. Maximum induction residual: {old["maximum_induction_residual_Debye"]:.3g} Debye.',
        '',f'{result["reference_rows_retaining_failed_printed_scf_criteria"]} admitted reference rows retain one or more failed printed SCF criteria; native admission does not establish reference accuracy.',
        '',f'Admitted reference elapsed-times-four-ranks proxy: {result["reference_allocated_elapsed_CPU_seconds_proxy"]:.3f} CPU-seconds. This is not measured CPU utilization. Scheduler cost is left to root accounting.',
        '', 'Model process cost and all exact output/source pins are in RESULT.json. No retry, target change, fit, or physics promotion follows from application admission.']
    if failures:lines += ['', 'Unavailable/failed slots:', *[f'- {k}: {v}' for k,v in failures.items()]]
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return r.pin(out/'RESULT.json')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['model','all']);ap.add_argument('--config',required=True)
    ap.add_argument('--model-admission');ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    cfg=json.loads(Path(a.config).read_text())
    result=model(cfg,a.out.resolve()) if a.mode=='model' else all_results(cfg,r.pin(a.model_admission),a.out.resolve())
    print(json.dumps({'result':result,'new_native_calls':0}))
