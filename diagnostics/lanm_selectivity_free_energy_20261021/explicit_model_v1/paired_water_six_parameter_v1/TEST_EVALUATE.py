#!/usr/bin/env python3
"""Finite saved-only cache/render/limit controls; no engine invocation."""
import importlib.util
import json
import os
import time
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OUT=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_six_parameter_v1/evaluator_controls_v3'
OUT.mkdir(parents=True,exist_ok=False)
spec=importlib.util.spec_from_file_location('six_water_evaluator_controls',HERE/'EVALUATE.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
start,cpu=time.monotonic(),time.process_time();checks=[];actual_calls=0
def check(ok,name):
    checks.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)
def blocked(call,text,name):
    try:call()
    except ValueError as exc:check(text in str(exc),name)
    else:raise AssertionError(name+' was not rejected')
def forbidden(*args,**kwargs):
    global actual_calls
    actual_calls+=1
    raise AssertionError('Native execution is forbidden in saved controls')
cfg=m.default_config();m.save(OUT/'CONFIG.json',cfg)
ev=m.Evaluator(cfg,OUT/'cache');ev.native.invoke=forbidden
prior=m.read(cfg['prior_result']);epsilon=m.read(cfg['epsilon_result'])
pcfg=m.read(prior['config']);pecfg=m.read(pcfg['evaluator_config']);ancestor=m.read(pecfg['prior_result'])
protected=[cfg['prior_result'],cfg['epsilon_result'],prior['evaluator']['ledger'],ancestor['evaluator']['ledger'],cfg['training_result']]
check(ev.inherited_counts=={'17593_native':736,'17589_native':672,'17596_epsilon':16},'all completed parent and probe cells indexed')
check(ev.continuation_start==list(m.START) and len(ev.family['parameters'])==6,'six-value exact initial point')
baseline=ev.evaluate(ev.continuation_start)
for gid in m.CELLS:
    check(all(baseline[gid][k]==prior['selected']['model'][gid][k] for k in ('energy','forces','native')),gid+':exact initial E/F replay')
    check(baseline[gid]['provenance']['native_authority']==prior['selected']['model'][gid]['provenance']['native_authority'],gid+':actual initial producer')
for name,probe in epsilon['records'].items():
    values=list(m.START);metal=probe['metal'];j=4 if metal=='La' else 5;values[j]=probe['epsilon']
    model=ev.evaluate(values)
    for gid in m.CELLS:
        check(all(model[gid][k]==probe['model'][gid][k] for k in ('energy','forces','native')),name+':'+gid+':exact epsilon raw E/F replay')
        if gid.endswith('_'+metal):
            pv=model[gid]['provenance']
            check(pv['inherited_source_kind']=='17596_epsilon' and pv['native_cache_key'] is None and pv['native_implementation']==epsilon['implementation'],name+':'+gid+':epsilon authority has no invented native cache')
        else:
            check(model[gid]['provenance']['cache_key']==baseline[gid]['provenance']['cache_key'],name+':'+gid+':unchanged opposite-metal identity')
prior_members={v['provenance']['native_authority']['path'] for vector in prior['vectors'].values() for v in vector['model'].values()}
older=next(vector for vector in ancestor['vectors'].values() if any(not value['provenance']['reused_baseline'] and value['provenance']['native_authority']['path'] not in prior_members for value in vector['model'].values()))
old_values=older['parameters']+[m.EPSILON['La'],m.EPSILON['Dy']]
old_model=ev.evaluate(old_values)
for gid in m.CELLS:
    check(all(old_model[gid][k]==older['model'][gid][k] for k in ('energy','forces','native')),gid+':ancestor outside17593 requests exactly reused')
    check(old_model[gid]['provenance']['native_authority']==older['model'][gid]['provenance']['native_authority'],gid+':ancestor producer retained')
for j in range(6):
    values=list(m.START);values[j]+=1e-5
    represented='La' if j in (0,1,4) else 'Dy'
    for gid in m.CELLS:
        before=ev.cache_identity(m.START,gid);after=ev.cache_identity(values,gid)
        check((before[0]!=after[0])==gid.endswith('_'+represented),'active-metal cache identity '+str(j)+':'+gid)
        check(after[1]['effective_values']==[values[0 if represented=='La' else 2],values[1 if represented=='La' else 3],values[4 if represented=='La' else 5]] if gid.endswith('_'+represented) else before[1]==after[1],'three effective values '+str(j)+':'+gid)
for j,metal in ((4,'La'),(5,'Dy')):
    bounds=ev.family['parameters'][j]
    check(bounds['lower']==.5*m.EPSILON[metal] and bounds['upper']==1.5*m.EPSILON[metal],metal+':positive frozen epsilon interval')
    for outside in (bounds['lower']-1e-8,bounds['upper']+1e-8,float('nan')):
        values=list(m.START);values[j]=outside
        blocked(lambda:ev.checked_values(values),'outside frozen bounds',metal+':invalid epsilon rejected')
blocked(lambda:ev.checked_values(m.START[:4]),'Six parameters required','four-value request rejected')
check(ev.snapshot()['native_program_invocations']==0 and actual_calls==0,'all inherited controls made zero native calls')
admitted=ev.snapshot()
old_keys=ev._ledger['vector_keys'];ev._ledger['vector_keys']=['synthetic-limit-control']*385
values=list(m.START);values[0]+=0.000137
blocked(lambda:ev.evaluate(values),'Parameter-vector ceiling reached','385 total-vector hard ceiling')
ev._ledger['vector_keys']=old_keys
ev._ledger['native_program_invocations']=4096;ev._gate=lambda:None
os.environ['SLURM_JOB_ID']='0'  # Test-only fixture; invoke remains forbidden.
key=ev.cache_identity(values,'water_minus_La')[0]
check(key not in ev._inherited and key not in ev._cache,'synthetic call-cap fixture is uncached')
blocked(lambda:ev._cell(values,'water_minus_La'),'Native call ceiling reached','4096 native-program hard ceiling')
check(ev.stop_event.is_set() and actual_calls==0,'synthetic ceiling stops before invoke')
for p in protected:m.verify(p)
result={'status':'passed_saved_evaluator_controls','implementation':m.pin(HERE/'EVALUATE.py'),'controls_source':m.pin(__file__),'config':m.pin(OUT/'CONFIG.json'),'check_count':len(checks),'checks':checks,'initial_and_inherited_admission':admitted,'synthetic_limit_test':{'vector_count':385,'program_count':4096,'actual_new_native_calls':actual_calls,'synthetic_cache_not_for_scientific_reuse':True},'cost':{'wall_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu,'maximum_threads':1,'native_calls':0},'protected_sources_unchanged':True}
m.save(OUT/'RESULT.json',result)
print(json.dumps({'result':m.pin(OUT/'RESULT.json'),'source':result['implementation'],'checks':len(checks),'cost':result['cost']}))
