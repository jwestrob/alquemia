#!/usr/bin/env python3
"""Saved baseline and explicitly synthetic report metadata; no fitting/native calls."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    tick,cpu=time.monotonic(),time.process_time();spec=importlib.util.spec_from_file_location('saved_fit_report_control',HERE/'REPORT.py');r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
    cfg=json.loads((HERE/'CONFIG_UNBOUND.json').read_text());domain=r.read(cfg['domain']);values=[p['initial'] for p in domain['parameters']];key=r.digest(values)
    original=json.loads(a.baseline.read_text());point={'status':'complete','parameters':values,**original};r.save(a.out/'SAVED_BASELINE_VECTOR.json',point);vp=r.pin(a.out/'SAVED_BASELINE_VECTOR.json')
    evcfg=HERE.parent/'evaluator_v1/CONFIG_EVALUATOR_UNBOUND.json';limits={'outer_steps':24,'new_vectors':672,'native_programs':24000,'workers':12}
    dc={'synthetic_report_control':True,'evaluator':cfg['evaluator'],'evaluator_config':r.pin(evcfg),'parameter_domain':cfg['domain'],'limits':limits,'root_release':None,'preparation_result':None}
    r.save(a.out/'SYNTHETIC_DRIVER_CONFIG.json',dc)
    entry={'status':'complete','parameters':values,'role':'initial','artifact':vp,'objective':point['score']['objective'],'block_objectives':point['score']['block_objectives']}
    fit={'synthetic_report_control':True,'stage':'run','status':'bounded_joint_response_fit_complete','implementation':cfg['driver'],'configuration':r.pin(a.out/'SYNTHETIC_DRIVER_CONFIG.json'),
         'scientific_fingerprint':r.digest({k:v for k,v in dc.items() if k not in ('root_release','preparation_result')}),'limits':limits,'vectors':{key:entry},'initial_key':key,'selected_key':key,'selected':vp,
         'candidate_frozen':True,'error':None,'iterations':[{'outer_step':1,'status':'complete','origin_key':key,'trial_keys':[key]*4,'selected_key':key,'actual_improvement':False}],
         'last_Jacobian_origin_key':key,'stopping_reason':'no_actual_improvement_at_complete_trial_batch','evaluator':{'native_program_invocations':0,'distinct_parameter_vectors':1,'requests':[]},
         'new_parameter_vectors':0,'cost':{'wall_seconds':0,'parent_CPU_seconds':0,'child_CPU_seconds':0},'job_id':'synthetic-not-executed'}
    r.save(a.out/'SYNTHETIC_TERMINAL.json',fit);cfg.update(fit_result=r.pin(a.out/'SYNTHETIC_TERMINAL.json'),initial_artifact=vp,selected_artifact=vp)
    checks={}
    def check(name,ok):r.require(ok,'Failed report control: '+name);checks[name]=True
    def reject(name,fn):
        try:fn()
        except (ValueError,KeyError):checks[name]=True;return
        raise AssertionError('Invalid report accepted: '+name)
    result=r.collect(cfg)
    check('full25_original_failed_gates_retained',result['selected']['EF']['energy_total']==10 and result['selected']['EF']['force_total']==15 and result['selected']['EF']['all25_passed'] is False and len(result['selected']['EF']['failed_ids'])==3)
    check('all34_pair_four_C_four_chi',len(result['selected']['pair_rows'])==34 and len(result['selected']['C'])==4 and len(result['selected']['chi'])==4 and sum(v['count'] for v in result['selected']['pair_contexts'].values())==34)
    check('four_engineering_blocks',set(result['selected']['block_objectives'])=={'old_EF','pair','C','chi'})
    check('coordination_sign_and_initial_identity',result['initial']==result['selected'] and result['selected']['paired_coordination']['model_minus_reference_kcal_mol']< -4)
    check('28_descriptive_field_rows',len(result['selected']['field_rows_descriptive'])==28)
    check('all8_regression_pending',result['consumed_trimer_regression']['status']=='pending_not_supplied' and len(result['consumed_trimer_regression']['cycles'])==8 and result['consumed_trimer_regression']['complete_cycles']==0)
    check('last_Jacobian_distinct_field',result['last_Jacobian_origin']['key']==key and result['last_Jacobian_origin']['parameters']==values)
    check('no_fake_response_pass',result['selected']['response_physical_acceptance'] is None)
    check('markdown_preserves_failed_gate_and_pending',all(x in r.markdown(result) for x in ('9/10 energy; 13/15 force','pending_not_supplied','-4.33828383')))
    # The failure route must not invent a frozen selected candidate.
    failed=copy.deepcopy(fit);failed.update(status='failed_or_unavailable',candidate_frozen=False,error='synthetic native stop',selected_key=None,selected=None)
    failed['vectors']['failed']={'status':'unavailable','error':'synthetic failed fixture','parameters':values}
    failed['evaluator']['requests']=[{'values':values,'status':'failed','error':'synthetic stop','cells':{'fixture':{'status':'failed_or_stopped','error':'synthetic native error'}}}]
    r.save(a.out/'SYNTHETIC_FAILED.json',failed);fc=dict(cfg,fit_result=r.pin(a.out/'SYNTHETIC_FAILED.json'),selected_artifact=None)
    f=r.collect(fc);check('failed_fit_never_selects',f['selected'] is None and f['error']=='synthetic native stop' and 'failed' in f['failed_or_unavailable_vectors'])
    check('failed_cell_preserved',f['failed_or_stopped_requests'][0]['cells']['fixture']['error']=='synthetic native error')
    running=copy.deepcopy(fit);running['status']='started';r.save(a.out/'SYNTHETIC_RUNNING.json',running)
    reject('running_result_not_terminal',lambda:r.collect(dict(cfg,fit_result=r.pin(a.out/'SYNTHETIC_RUNNING.json'))))
    damaged=copy.deepcopy(point);damaged['score']['objective']+=1;r.save(a.out/'SYNTHETIC_WRONG_SCORE.json',damaged)
    reject('wrong_exact_selected_artifact',lambda:r.collect(dict(cfg,selected_artifact=r.pin(a.out/'SYNTHETIC_WRONG_SCORE.json'))))
    incomplete=copy.deepcopy(original);del incomplete['score']['old_EF']['energy_controls']['balanced_exchange']['aquo9_central']
    reject('missing_original_gate',lambda:r.summarize(incomplete,r.read(cfg['development_data']),r.read(cfg['data_contract'])))
    r.save(a.out/'RESULT.json',{'status':'passed_saved_report_controls','implementation':r.pin(HERE/'REPORT.py'),'test_implementation':r.pin(__file__),'configuration':r.pin(HERE/'CONFIG_UNBOUND.json'),
         'baseline':r.pin(a.baseline),'checks':checks,'count':len(checks),'synthetic_terminal_metadata_only':True,'fit_executed':False,'new_native_calls':0,
         'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}})
    print(json.dumps({'result':r.pin(a.out/'RESULT.json'),'count':len(checks)}))


if __name__=='__main__':main()
