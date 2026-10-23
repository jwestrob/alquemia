#!/usr/bin/env python3
"""Saved-only initial admission and numerical interpretation of model job18473."""
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('actual_field_model',HERE.parent/'RUN_v2.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
W=a.ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/conditional_susceptibility_v1/model'
out=W/'initial_actual_review_v1';out.mkdir(exist_ok=False)
cfg_path=W/'run_v1/CONFIG.json';result_path=W/'run_v1/execution/RESULT.json'
cfg=json.loads(cfg_path.read_text());m=a.read(cfg['manifest']);r=json.loads(result_path.read_text());tick=time.monotonic();cpu=time.process_time()
admitted=a.admit_result(cfg,m,r)
metrics={}
for section in ('states','increments'):
    values=r['summary'][section];steps=[step for v in values.values() for step in v['steps'].values()]
    metrics[section]={
        'count':len(values),'all_passed':all(v['passed'] for v in values.values()),
        'checks_passed':sum(sum(v['checks'].values()) for v in values.values()),
        'maximum_energy_moment_discrepancy_a0_cubed':max(abs(v['chi_E_a0_cubed']-v['chi_mu_column_a0_cubed'][2]) for v in steps),
        'maximum_energy_step_difference_a0_cubed':max(abs(v['steps']['h']['chi_E_a0_cubed']-v['steps']['half']['chi_E_a0_cubed']) for v in values.values()),
        'maximum_moment_step_difference_a0_cubed':max(abs(v['steps']['h']['chi_mu_column_a0_cubed'][2]-v['steps']['half']['chi_mu_column_a0_cubed'][2]) for v in values.values()),
        'maximum_odd_slope_zero_difference_au':max(abs(v['mu_E_minus_mu0_au']) for v in steps),
        'energy_moment_agreement_within_combined_print_bounds':all(abs(v['chi_E_a0_cubed']-v['chi_mu_column_a0_cubed'][2])<=v['energy_curvature_print_bound_a0_cubed']+v['moment_slope_print_bound_a0_cubed'] for v in steps)}
value={'status':'passed_initial_saved_model_review','implementation':a.pin(__file__),'executed_adapter':cfg['implementation'],'result':a.pin(result_path),'config':a.pin(cfg_path),
       'manifest':cfg['manifest'],'release':cfg['root_release'],'saved_admission_status':admitted['status'],'native_program_invocations':19,'complete_rows':19,
       'controls':{k:r['rows'][k]['control'] for k in a.CONTROL_IDS},'metrics':metrics,'summary':r['summary'],
       'maximum_induction_residual_Debye':max(x['moments']['induction']['maximum_printed_RMS_residual_Debye'] for x in r['rows'].values()),
       'induction_solutions':sum(x['moments']['induction']['solutions'] for x in r['rows'].values()),
       'actual_application_cost':r['cost'],'new_native_calls':0,'reference_agreement_evaluated':False,
       'initial_display_note':'A preliminary saved admission succeeded; its later optional maximum-residual display referenced the nonexistent key maximum_final_residual_Debye. This receipt uses the existing maximum_printed_RMS_residual_Debye. No executed source/result/native call changed.',
       'review_cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1}}
a.save(out/'REVIEW.json',value);print(json.dumps(a.pin(out/'REVIEW.json')))
