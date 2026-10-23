#!/usr/bin/env python3
"""One prospective eight-step component campaign; no continuation or retry."""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import resource
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('component_policy',HERE/'ADAPTER.py');A=importlib.util.module_from_spec(s);s.loader.exec_module(A)


def persist_vector(out,key,value):
    path=out/'vectors'/(key+'.json');A.require(not path.exists(),'Immutable vector artifacts');A.save(path,value);return A.pin(path)


def run(mode,config_path,out):
    tick,cpu=time.monotonic(),time.process_time();children=resource.getrusage(resource.RUSAGE_CHILDREN)
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False);(out/'vectors').mkdir()
    cfg=json.loads(Path(config_path).read_text());ev=None
    result={'status':'started','mode':mode,'implementation':A.pin(__file__),'configuration':A.pin(config_path),'configuration_fingerprint':A.fingerprint(cfg),'error':None,'vectors':{},'iterations':[],
            'limits':A.LIMITS,'job_id':os.environ.get('SLURM_JOB_ID'),'native_programs':0,'real_proposal_calls':0,'physical_accuracy_qualified':False,'reserved_outcomes_opened':False,'automatic_extension':False}
    current=None;fixed=None;core=None;values={}
    def persist():
        if ev is not None:result['evaluator']=ev.snapshot();result['native_programs']=result['evaluator']['native_program_invocations'];result['new_parameter_vectors']=result['evaluator']['new_parameter_vectors']
        A.save(out/'RESULT.json',result)
    def obtain(parameters,role,step):
        key=A.digest(parameters);ev.authorize(parameters,role,step)
        usage={'outer_step':step,'role':role}
        if key in values:
            result['vectors'][key].setdefault('requested_roles',[]).append(usage);return values[key]
        rec={'parameters':parameters[:],'status':'started','requested_roles':[usage],'inherited':key in ev.imported_vector_keys};result['vectors'][key]=rec;persist()
        try:
            v=ev.evaluate(parameters);values[key]=v;rec.update(status='complete',artifact=persist_vector(out,key,v),score=v['score']['block_objectives'])
            persist();return v
        except Exception as exc:
            rec.update(status='budget_not_attempted' if isinstance(exc,A.BudgetStop) else 'failed',error=repr(exc));persist();raise
    try:
        A.require(cfg['implementation']==A.pin(__file__),'Driver pin')
        if mode=='run':A.require(str(out)==cfg['output_directory'],'Frozen output directory')
        base,parent,inv,seeds=A.context(cfg);core=A.module(cfg['core'],'component_core')
        ev=A.evaluator_class(base)(cfg,out/'cache',mode,inv,seeds)
        current=copy.deepcopy(seeds[inv['origin_key']]);fixed=core.make_fixed(current,ev.m,A.read(cfg['decomposition']),ev.domain,A.read(cfg['declaration']))
        expected=core.probe_points(current['parameters'],ev.domain)
        A.require([seeds[k]['parameters'] for k in inv['derivative_keys']]==expected,'Exact first24 saved probe coordinates/order')
        for key,v in seeds.items():
            values[key]=v;result['vectors'][key]={'parameters':v['parameters'],'status':'complete','artifact':persist_vector(out,key,v),'inherited':True,'import_role':'optional_pilot' if key==inv['pilot_key'] else ('origin' if key==inv['origin_key'] else 'initial_probe'),'requested_roles':[]}
        result.update(origin_key=inv['origin_key'],current_key=inv['origin_key'],reuse_inventory=cfg['reuse_inventory'],initial_metrics=core.assess(current,fixed),immutable_saved_scores_preserved=True)
        A.require(result['initial_metrics']['physical_constraints_pass'],'Origin25/fixedceilings')
        # These metadata/count checks invoke no numerical proposal or new point.
        A.require(all(ev.required_programs(v['parameters'])==0 for v in seeds.values()),'All26 immutable vectors completely reusable')
        result['import_counts']={'origin':1,'initial_derivatives':24,'additional_pilot':1,'unique_cells':len(ev.imported_cell_keys)}
        if mode=='prepare':
            result.update(status='prepared_component_fit_zero_native',preparation_checks={'exact_first_stencil':True,'origin_constraints':True,'all26_reuse_zero_programs':True,'no_real_proposal':True})
        else:
            ev._gate();persist();stop=None
            for step in range(1,9):
                iteration={'outer_step':step,'origin_key':A.digest(current['parameters']),'status':'started','derivative_keys':[],'trial_slots':[]};result['iterations'].append(iteration);persist()
                points=core.probe_points(current['parameters'],ev.domain);probes=[]
                try:
                    for point in points:
                        value=obtain(point,'probe',step);probes.append(value);iteration['derivative_keys'].append(A.digest(point));persist()
                    proposal=core.propose(current,probes,fixed);result['real_proposal_calls']+=1
                    A.require(proposal['solver_calls']==1,'Exactly one solver attempt per step')
                    path=out/f'PROPOSAL_{step:02d}.json';A.save(path,proposal);iteration['proposal']=A.pin(path)
                    slots=copy.deepcopy(proposal['trial_slots']);iteration['trial_slots']=[{k:v for k,v in x.items() if k!='certificate'} for x in slots];persist()
                    if not slots or not any(x['eligible'] for x in slots):
                        iteration['status']='no_eligible_fraction';stop='no_eligible_fraction';break
                    A.require(len(slots)==4 and [x['fraction'] for x in slots]==[1.,.5,.25,.125],'Exact frozen four fractions')
                    for i,slot in enumerate(slots):
                        if not slot['eligible']:slot['execution_status']='not_eligible';continue
                        A.require(A.digest(slot['parameters'])==slot['key'],'Certified trial identity')
                        value=obtain(slot['parameters'],'trial',step)
                        slot.update(execution_status='complete',native_admitted=True,vector=value)
                        iteration['trial_slots'][i].update(execution_status='complete',native_admitted=True,artifact=result['vectors'][slot['key']]['artifact']);persist()
                    # Only this complete-batch boundary can alter the incumbent.
                    selection=core.select(current,slots,fixed);iteration['selection']=selection
                    chosen=selection['selected_trial_index']
                    if chosen is None:iteration['status']='no_actual_improvement';stop='no_actual_admitted_improvement';break
                    current=copy.deepcopy(slots[chosen]['vector']);iteration.update(status='complete',selected_key=selection['selected_key']);result['current_key']=selection['selected_key'];persist()
                    if selection['selected_metrics']['full_development_target']:stop='development_target_met';break
                except A.BudgetStop:
                    iteration['status']='hard_budget_exhausted';stop='hard_budget_exhausted';break
                except Exception:
                    iteration['status']='failed_incomplete_batch_no_selection';raise
            result.update(status='bounded_component_fit_complete',stopping_reason=stop or 'outer_step_limit',completed_outer_steps=sum(x['status']=='complete' for x in result['iterations']))
        metrics=core.assess(current,fixed);key=A.digest(current['parameters'])
        result.update(selected_key=key,selected=result['vectors'][key]['artifact'],selected_metrics=metrics,candidate_frozen=(mode=='run'),native_programs=ev.snapshot()['native_program_invocations'])
    except Exception as exc:
        if ev is not None:ev._stop(exc)
        result.update(status='failed_component_fit' if mode=='run' else 'failed_preparation',error=repr(exc),stopping_reason='failure_no_retry')
        if current is not None and core is not None and fixed is not None:
            key=A.digest(current['parameters']);result.update(last_complete_incumbent_key=key,last_complete_incumbent=result['vectors'].get(key,{}).get('artifact'),last_complete_incumbent_metrics=core.assess(current,fixed),candidate_frozen=False)
    finally:
        now=resource.getrusage(resource.RUSAGE_CHILDREN)
        result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':now.ru_utime-children.ru_utime,'child_system_CPU_seconds':now.ru_stime-children.ru_stime,'numerical_threads':1,'native_workers':12,'scheduler_allocated_CPU_seconds':0,'execution_mode':'login','scheduler_allocation':None,'inherited_native_costs_not_added_to_new_cost':True}
        persist()
    print(json.dumps({'result':A.pin(out/'RESULT.json'),'status':result['status'],'error':result['error'],'native_programs':result['native_programs']}))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('prepare','run'));p.add_argument('--config',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    r=run(a.mode,a.config,a.out);raise SystemExit(1 if r['error'] else 0)
