#!/usr/bin/env python3
"""Narrow v2 reporting-delta control, reusing the saved v1 linear fixture."""
import importlib.util
import json
from pathlib import Path
import time
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('model_v2',HERE/'RUN_v2.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
W=a.ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/conditional_susceptibility_v1/model'
t=time.monotonic();cpu=time.process_time();out=W/'saved_controls_v2';out.mkdir(exist_ok=False)
cfg=json.loads((W/'preparation_v2/CONFIG_RUN_TEMPLATE.json').read_text());m,_=a.context(cfg)
old=json.loads((W/'saved_controls_v1/complete/RESULT.json').read_text());value=a.summarize(m,old['rows'])
keys={'mu_E_print_bound_au','mu0_print_bound_au','mu_E_minus_mu0_print_bound_au'}
def strip(v):
    if isinstance(v,dict):return {k:strip(x) for k,x in v.items() if k not in keys}
    if isinstance(v,list):return [strip(x) for x in v]
    return v
checks=[]
assert strip(value)==old['summary'];checks.append('all_previous_arithmetic_and_gates_exact')
for state,item in value['states'].items():
    for step,v in item['steps'].items():
        assert v['mu_E_print_bound_au']==5e-9/(a.HARTREE*v['field_au'])
        assert v['mu0_print_bound_au']==m['rows'][state]['baseline']['total_precise_print_bound_au_per_component']
        assert v['mu_E_minus_mu0_print_bound_au']==v['mu_E_print_bound_au']+v['mu0_print_bound_au']
        checks.append(state+'_'+step+'_odd_bound')
for name,item in value['increments'].items():
    for step,v in item['steps'].items():
        for key in keys:assert abs(v[key]-sum(abs(c)*value['states'][state]['steps'][step][key] for state,c in item['terms'].items()))<1e-18
        checks.append(name+'_'+step+'_absolute_coefficient_bound')
a.save(out/'RESULT.json',{'status':'passed_reporting_delta_controls','implementation':a.pin(HERE/'RUN_v2.py'),'control_source':a.pin(__file__),'manifest':a.pin(W/'preparation_v2/MANIFEST.json'),
                        'original_source':a.pin(HERE/'RUN.py'),'prior_controls':a.pin(W/'saved_controls_v1/RESULT.json'),'prior_saved_fixture':a.pin(W/'saved_controls_v1/complete/RESULT.json'),
                        'checks':checks,'passed':len(checks),'new_native_calls':0,'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'local_threads':1}})
print(json.dumps(a.pin(out/'RESULT.json')))
