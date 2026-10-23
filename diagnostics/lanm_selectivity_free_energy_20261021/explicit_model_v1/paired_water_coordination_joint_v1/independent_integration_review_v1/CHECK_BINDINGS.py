#!/usr/bin/env python3
"""Read-only final source/preparation integration joins; no engine imports."""
import hashlib
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent;BRANCH=HERE.parent;ROOT=BRANCH.parents[3]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_coordination_joint_v1'

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,p
    return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())

def main():
    start,cpu=time.monotonic(),time.process_time();checks=[]
    def check(ok,name):
        checks.append({'name':name,'passed':bool(ok)})
        assert ok,name
    cp=pin(BRANCH/'CONFIG_PREPARE.json');cfg=read(cp)
    for value in cfg.values():
        if isinstance(value,dict) and set(value)=={'path','sha256'}:verify(value)
    prep_pin=pin(W/'preparation_v1/RESULT.json');prep=read(prep_pin)
    check(prep['status']=='prepared' and prep['error'] is None and prep['new_native_program_invocations']==0 and prep['config']==cp and prep['implementation']==cfg['implementation'],'Actual zero-call driver preparation/source/config join')
    fingerprint=hashlib.sha256(json.dumps({k:v for k,v in cfg.items() if k not in ('root_release','preparation_result')},sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    check(prep['scientific_config_sha256']==fingerprint and cfg['root_release'] is None and cfg['preparation_result'] is None,'Preparation configuration fingerprint and unreleased state')
    ecfg=read(cfg['evaluator_config']);ep_pin=pin(W/'evaluator_preparation_v2/RESULT.json');ep=read(ep_pin)
    check(ep['status']=='prepared_no_calls' and ep['new_native_program_invocations']==0 and ep['config']==cfg['evaluator_config'] and ep['implementation']==cfg['evaluator'],'Actual evaluator preparation joins final evaluator config/source')
    check(ecfg['training_result']==cfg['training_result'] and ecfg['implementation']==cfg['evaluator'],'Driver/evaluator input authority join')
    check([ecfg[k] for k in ('maximum_parameter_vectors','maximum_native_program_invocations','maximum_native_concurrency')]==[385,5120,8],'Exact evaluator limits')
    data=read(cfg['training_result']);prior=read(cfg['prior_result']);coord=read(cfg['coordination_result'])
    check(len(data['rows'])==10 and len(prep['initial']['model'])==10 and not data['orientation_in_training'],'Exact ten consumed cells and no orientation fitting')
    check(prep['initial']['parameters']==prior['selected']['parameters']==ep['initial_parameters']==ecfg['initial_vector'],'Exact17605 start across integrations')
    for gid,row in data['rows'].items():
        expected=coord['rows'][gid]['model'] if gid.startswith('aquo9_') else prior['selected']['model'][gid]
        check(all(prep['initial']['model'][gid][k]==ep['initial_model'][gid][k]==expected[k] for k in ('energy','forces','native')),gid+' exact saved native arrays through driver/evaluator')
        check(row['initial_input']==ecfg['baseline_inputs'][gid] and all(row['input'][k]==ecfg['source_inputs'][gid][k] for k in ('xyz','key','parameters')),gid+' source and initial physical-input binding')
        for record in list(row['initial_input'].values())+[row['input'][k] for k in ('xyz','key','parameters')]:verify(record)
    controls_pin=pin(W/'objective_controls_v1/RESULT.json');controls=read(controls_pin)
    check(controls['status']=='passed' and controls['checks_passed']==51 and controls['objective_source']==cfg['objective'],'Existing51 objective controls bind final objective without rerun')
    check(prep['initial']['score']==controls['actual_baseline_score'],'Actual driver baseline score matches independently controlled scorer')
    check(prep['initial_old_eight_score']==prior['selected']['score'],'Old eight-cell score exactly preserved')
    check(prep['evaluator']['native_program_invocations']==0 and prep['evaluator']['inherited_native_cells_indexed']==490 and prep['evaluator']['inherited_native_cells_admitted']==10,'Actual effective-input cache index and zero-call admission counts')
    for record in prep['protected_parents']:verify(record)
    notice=read(ecfg['coordination_provenance_notice'])
    check(data['historical_review_continuity_broken'] and notice['reconstruction_attempted'] is False and notice['science_rerun'] is False,'Historical missing root-review bytes remain explicit; no reconstruction claim')
    result={'status':'bindings_passed','implementation':pin(__file__),'driver':cfg['implementation'],'evaluator':cfg['evaluator'],'driver_config':cp,'evaluator_config':cfg['evaluator_config'],'training':cfg['training_result'],'driver_preparation':prep_pin,'evaluator_preparation':ep_pin,'objective_controls':controls_pin,'root_declaration':cfg['root_declaration'],'specification':cfg['specification'],'checks':checks,'checks_passed':len(checks),'native_calls':0,'new_objective_controls_run':False,'review_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    out=W/'independent_integration_review_v1';out.mkdir(parents=True,exist_ok=False);p=out/'BINDINGS.json';p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'result':pin(p),'checks':len(checks),'wall':result['wall_seconds']}))

if __name__=='__main__':main()
