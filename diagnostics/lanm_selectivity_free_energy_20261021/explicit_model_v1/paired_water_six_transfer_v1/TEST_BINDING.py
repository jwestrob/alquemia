#!/usr/bin/env python3
"""Actual saved producer joins and conditional-admission controls; no candidate bind."""
import copy
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
OUT=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_six_transfer_v1/saved_controls_v2'
OUT.mkdir(parents=True,exist_ok=False)
spec=importlib.util.spec_from_file_location('six_transfer_saved_controls',HERE/'RUN.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
start,cpu=time.monotonic(),time.process_time();checks=[]
def check(ok,name):
    checks.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)
def reject(call,reason,name):
    try:call()
    except ValueError as exc:check(reason in str(exc),name)
    else:raise AssertionError(name+' accepted unexpectedly')
data=m.read(m.PINS['pins']['prepared32']);training=m.read(m.PINS['pins']['radial_training'])
cfg=m.read(m.PINS['pins']['six_evaluator_config']);prep=m.read(m.PINS['pins']['six_preparation'])
native=m.module(data['config']['native_helper'],'six_transfer_saved_native')
native.invoke=lambda *a,**k: (_ for _ in ()).throw(AssertionError('Native call forbidden'))
sources={'17593':m.read(m.PINS['pins']['six_prior_result']),'17589':m.read(m.PINS['pins']['continuation_prior_result']),'17596':m.read(m.PINS['pins']['epsilon_result'])}
previous=m.read(m.PINS['pins']['previous_transfer_result'])
controls=m.read(m.PINS['pins']['six_saved_controls'])
requests=controls['initial_and_inherited_admission']['requests'];examples={}
for request in requests:
    for gid,pv in request['cells'].items():
        label='17578' if pv['reused_baseline'] else pv['inherited_source_kind']
        examples.setdefault(label,(request['vector'],gid,pv))
check(set(examples)=={'17578','17593_native','17589_native','17596_epsilon'},'four actual inherited producer cases available')
branches={}
candidate_path=Path(m.PINS['pins']['six_run_config']['path']).parent/'execution/RESULT.json'
for label,(values,gid,pv) in examples.items():
    authority=m.read(pv['native_authority'])
    raw=authority['rows'][gid]['model'] if label=='17578' else authority['model'] if label=='17596_epsilon' else authority
    value={k:copy.deepcopy(raw[k])for k in ('energy','forces','native')};value['provenance']=copy.deepcopy(pv)
    row=data['rows'][gid];inp=pv['effective_input']
    key,render=m.render(data['config'],row,values)
    check(key==m.verify(inp['key']).read_text() and render['water_recipient'],label+':exact six-rendered physical input')
    branches[label]=m.selected_authority(gid,value,row,inp,cfg,candidate_path,previous,sources,values)
    m.replay_native(value,row,native);check(True,label+':actual raw E/F replay')
    changed=copy.deepcopy(value);changed['provenance']['cache_key']='wrong-current-request'
    reject(lambda:m.selected_authority(gid,changed,row,inp,cfg,candidate_path,previous,sources,values),'Current six-value cache identity differs',label+':wrong current request rejected')
    changed=copy.deepcopy(value);changed['provenance']['native_implementation']=cfg['implementation']
    reason='Baseline native source relabeled' if label=='17578' else 'Epsilon native source relabeled' if label=='17596_epsilon' else 'Actual native producer relabeled'
    reject(lambda:m.selected_authority(gid,changed,row,inp,cfg,candidate_path,previous,sources,values),reason,label+':native relabeling rejected')
    if label=='17596_epsilon':
        changed=copy.deepcopy(value);changed['provenance']['native_cache_key']='invented'
        reject(lambda:m.selected_authority(gid,changed,row,inp,cfg,candidate_path,previous,sources,values),'Epsilon native source relabeled','epsilon invented native-cache identity rejected')
    if label=='17593_native':
        changed=copy.deepcopy(value);changed['provenance']['inherited_native']=False
        reject(lambda:m.selected_authority(gid,changed,row,inp,cfg,candidate_path,previous,sources,values),'New native outside exact six-fit cache','old cache cannot masquerade as new17605 cache')
m.family_gate(prep['family'],training,cfg);check(True,'actual six-family preparation admitted')
bad=copy.deepcopy(prep['family']);bad['parameters'][4]['lower']*=.5
reject(lambda:m.family_gate(bad,training,cfg),'Epsilon bounds differ','relaxed epsilon family rejected')
bad=copy.deepcopy(prep['family']);bad['parameters'][0]['upper']+=.1
reject(lambda:m.family_gate(bad,training,cfg),'Original four parameter bounds differ','changed radius family rejected')
reject(lambda:m.candidate_gate({'status':'unavailable'}),'Only a complete error-free frozen six fit may bind','running-placeholder status cannot bind')
reject(lambda:m.candidate_gate({'status':'bounded_six_parameter_water_complete','error':'failure','candidate_frozen':False}),'Only a complete error-free frozen six fit may bind','failed six fit cannot bind')
fixture={'status':'bounded_six_parameter_water_complete','error':None,'candidate_frozen':True,'implementation':m.PINS['pins']['six_source'],'config':m.PINS['pins']['six_run_config'],'job_id':'17605','original_results_unchanged':True,'development_checks_passed':False}
reject(lambda:m.candidate_gate(fixture),'All original water controls must pass','control-failing terminal status cannot bind')
check(len(m.d.CELLS)==32 and len(m.WATER)==8 and len(m.d.HYDRATED)==18 and len(m.d.DRY)==6 and len(m.d.FD)==10,'unchanged complete inventory')
check(m.PINS['maximum_new_native_program_invocations']==40 and m.PINS['maximum_internal_displaced_energy_evaluations']==1380,'unchanged finite call scope')
check(not (OUT/'BOUND.json').exists() and m.PINS['candidate_result'] is None,'no real or synthetic candidate bound')
result={'status':'passed_saved_six_transfer_controls','implementation':m.pin(HERE/'RUN.py'),'controls_source':m.pin(__file__),'pin_inventory':m.pin(HERE/'PINS.json'),'renderer':m.PINS['pins']['renderer'],'check_count':len(checks),'checks':checks,'actual_producer_branches':branches,'native_calls':0,'candidate_bound':False,'six_fit_outcomes_opened':False,'reserved_validation_opened':False,'cost':{'wall_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu,'maximum_threads':1}}
m.save(OUT/'RESULT.json',result);print(json.dumps({'result':m.pin(OUT/'RESULT.json'),'status':result['status'],'checks':len(checks),'cost':result['cost']}))
