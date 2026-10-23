#!/usr/bin/env python3
"""Saved terminal fit reporting; no evaluator construction or native invocation."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(q):require(pin(q['path'])==q,'Changed artifact: '+q['path']);return Path(q['path'])
def read(q):return json.loads(verify(q).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def module(q,name):
    spec=importlib.util.spec_from_file_location(name,verify(q));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def summarize(value,data,contract):
    score=value['score'];old=score['old_EF'];ef=contract['development_set']['energy_force'];energies={}
    prefix={'balanced_exchange':'exchange__','individual_radial_work':'work__','paired_radial_work':'paired_work__'}
    for group,items in old['energy_controls'].items():
        for name,item in items.items():energies[prefix[group]+name]=copy.deepcopy(item)
    forces={name+'__'+group:{k:r[k] for k in ('rigid_RMS_kcal_mol_A','maximum_body_component_error_kcal_mol_A','passed')}
            for group,rows in old['per_fixture'].items() for name,r in rows.items()}
    require(set(energies)==set(ef['energy_gate_ids']) and set(forces)==set(ef['force_gate_ids']),'Exact original10+15 gate denominator')
    passed=all(r['passed'] is True for r in list(energies.values())+list(forces.values()))
    require(passed==score['original_EF_checks_passed']==old['development_checks_passed'],'Original gate summary differs')
    details=score['response']['per_observable'];parents=contract['development_set']['response']['C_parent_ids'];pairs={}
    require(set(details['pair'])==set(contract['development_set']['response']['pair_ids']) and set(details['C'])==set(parents) and set(details['chi'])==set(contract['development_set']['response']['state_ids']),'Exact34/4/4 response denominator')
    for parent in parents:
        rows={k:v for k,v in details['pair'].items() if v['parent']==parent};errors=[v['error_Debye'] for v in rows.values()]
        require(len(rows)==(8 if parent.startswith('aquo8_') else 9),'Literal pair context count')
        norms=[math.sqrt(math.fsum(x*x for x in vector)) for vector in errors]
        pairs[parent]={'count':len(rows),'component_RMS_error_Debye':math.sqrt(math.fsum(x*x for vector in errors for x in vector)/(3*len(rows))),
                       'maximum_absolute_component_error_Debye':max(abs(x) for vector in errors for x in vector),
                       'maximum_vector_error_norm_Debye':max(norms),'mean_vector_error_norm_Debye':math.fsum(norms)/len(norms),
                       'mean_signed_error_Debye':[math.fsum(v[j] for v in errors)/len(errors) for j in range(3)]}
    coefficients={'aquo9_central_Dy':1,'aquo9_central_La':-1,'aquo8_central_Dy':-1,'aquo8_central_La':1}
    C=math.fsum(c*value['model'][gid]['energy'] for gid,c in coefficients.items())
    ref=math.fsum(c*data['rows'][gid]['reference']['energy'] for gid,c in coefficients.items())
    error=energies['exchange__aquo9_central']['error_kcal_mol']
    require(abs((C-ref)-error)<=1e-10,'Paired coordination arithmetic closure')
    return {'objective':score['objective'],'block_objectives':score['block_objectives'],
            'EF':{'all25_passed':passed,'energy_passed':sum(r['passed'] for r in energies.values()),'energy_total':10,
                  'force_passed':sum(r['passed'] for r in forces.values()),'force_total':15,'energy':energies,'forces':forces,
                  'failed_ids':[k for k,r in {**energies,**forces}.items() if not r['passed']]},
            'paired_coordination':{'definition':'(Dy W9-La W9)-(Dy W8-La W8)','model_kcal_mol':C,'reference_kcal_mol':ref,'model_minus_reference_kcal_mol':error},
            'chi':copy.deepcopy(details['chi']),'C':copy.deepcopy(details['C']),'pair_contexts':pairs,'pair_rows':copy.deepcopy(details['pair']),
            'field_rows_descriptive':copy.deepcopy(score['response']['field_rows']),'response_physical_acceptance':None}


def trajectory(fit):
    records=fit['vectors'];rows=[]
    def brief(key):
        r=records.get(key,{})
        return {'key':key,'status':r.get('status','unavailable'),'objective':r.get('objective'),'block_objectives':r.get('block_objectives'),'error':r.get('error')}
    for entry in fit['iterations']:
        rows.append({'outer_step':entry['outer_step'],'status':entry['status'],'origin':brief(entry['origin_key']),
                     'derivative_count':len(entry.get('derivative_keys',[])),'trials':[brief(k) for k in entry.get('trial_keys',[])],
                     'selected':brief(entry['selected_key']) if entry.get('selected_key') else None,'actual_improvement':entry.get('actual_improvement'),
                     'linear_direction_diagnostics':entry.get('linear_directions')})
    key=fit.get('last_Jacobian_origin_key')
    if key is None:
        key=next((e['origin_key'] for e in reversed(fit['iterations']) if 'Jacobian_normalized' in e),None)
    return rows, None if key is None else {**brief(key),'parameters':records[key]['parameters'],'is_selected_point':key==fit.get('selected_key')}


def regression(cfg,fit,contract):
    ids=contract['consumed_regression_outside_proposed_fit']['cycle_ids'];q=cfg.get('regression_result')
    if q is None:return {'status':'pending_not_supplied','expected_cycles':8,'complete_cycles':0,'cycles':{cid:{'status':'pending'} for cid in ids},'fresh_validation':False}
    raw=read(q);require(raw['fit_result']==cfg['fit_result'] and raw['selected']==cfg['selected_artifact']==fit.get('selected'),'Regression selected-fit authority differs')
    rows={cid:copy.deepcopy(raw.get('cycles',{}).get(cid,{'status':'unavailable'})) for cid in ids}
    complete=sum(r.get('status')=='complete' for r in rows.values())
    require(raw.get('status')!='consumed_trimer_regression_complete' or complete==8,'Complete regression lost a cycle')
    return {'status':raw.get('status','unavailable'),'source':q,'expected_cycles':8,'complete_cycles':complete,'cycles':rows,'error':raw.get('error'),'fresh_validation':False,'new_physical_threshold':None}


def collect(cfg):
    require(cfg['implementation']==pin(__file__),'Report implementation binding')
    for name in ('driver','evaluator','objective','old_objective','body_helper','data_contract','domain','targets','development_data'):verify(cfg[name])
    fit=read(cfg['fit_result']);require(fit['stage']=='run' and fit['status'] in ('bounded_joint_response_fit_complete','failed_or_unavailable') and 'cost' in fit,'Actual terminal run required; no preparation or running placeholder')
    require(fit['implementation']==cfg['driver'],'Frozen fit driver differs')
    drivercfg=read(fit['configuration']);evcfg=read(drivercfg['evaluator_config'])
    require(fit['limits']==drivercfg['limits']=={'outer_steps':24,'new_vectors':672,'native_programs':24000,'workers':12},'Declared fit limits differ')
    require(drivercfg['evaluator']==cfg['evaluator'] and drivercfg['parameter_domain']==cfg['domain'],'Driver evaluator/domain joins')
    require(all(evcfg[k]==cfg[k] for k in ('objective','old_objective','body_helper','targets','development_data')),'Evaluator scientific source joins')
    require(fit['scientific_fingerprint']==digest({k:v for k,v in drivercfg.items() if k not in ('root_release','preparation_result')}),'Terminal scientific fingerprint differs')
    contract=read(cfg['data_contract']);data=read(cfg['development_data']);domain=read(cfg['domain']);targets=read(cfg['targets'])
    require(contract['authorities']['development_data']==cfg['development_data'] and contract['authorities']['development_objective']==cfg['old_objective'] and contract['authorities']['body_helper']==cfg['body_helper'],'Data contract source joins')
    helpers=[module(cfg[k],'fit_report_'+k) for k in ('objective','old_objective','body_helper')]
    def point(q,key):
        require(q==fit['vectors'][key]['artifact'],'Exact terminal vector artifact required');value=read(q)
        require(value['status']=='complete' and value['parameters']==fit['vectors'][key]['parameters'] and digest(value['parameters'])==key,'Complete vector/key/parameter joins')
        require(len(value['cells'])==46 and all(r['status']=='complete' for r in value['cells'].values()),'Full46 saved vector denominator')
        replay=helpers[0].evaluate(data,value['model'],value['response'],targets,helpers[1],helpers[2])
        require(replay==value['score'] and replay['objective']==fit['vectors'][key]['objective'],'Exact saved objective replay')
        return value,summarize(value,data,contract)
    initial=None;initial_summary=None
    if fit.get('initial_key'):
        require(cfg.get('initial_artifact') is not None,'Bind exact initial artifact');initial,initial_summary=point(cfg['initial_artifact'],fit['initial_key'])
        require(initial['parameters']==[p['initial'] for p in domain['parameters']],'Declared starting vector differs')
    else:require(cfg.get('initial_artifact') is None,'No initial artifact was produced')
    selected=None;selected_summary=None;bounds=[]
    if fit['status']=='bounded_joint_response_fit_complete':
        require(fit['error'] is None and fit['candidate_frozen'] is True and cfg.get('selected_artifact')==fit['selected'],'Actual complete frozen selection required')
        selected,selected_summary=point(cfg['selected_artifact'],fit['selected_key'])
        require(len(selected['parameters'])==len(domain['parameters'])==12,'Exact selected parameter denominator')
        for p,x in zip(domain['parameters'],selected['parameters']):
            require(math.isfinite(x) and p['lower']<=x<=p['upper'],'Selected parameter outside domain')
            bounds.append({'name':p['name'],'value':x,'lower':p['lower'],'upper':p['upper'],'at_lower':x==p['lower'],'at_upper':x==p['upper']})
    else:require(cfg.get('selected_artifact') is None,'Failed fit has no frozen selected candidate')
    path,lastJ=trajectory(fit);snapshot=fit.get('evaluator',{});failures={k:copy.deepcopy(r) for k,r in fit['vectors'].items() if r['status']!='complete'}
    requests=[{'values':r['values'],'status':r['status'],'error':r.get('error'),'cells':{cid:c for cid,c in r['cells'].items() if c['status']!='complete'}}
              for r in snapshot.get('requests',[]) if r['status']!='complete']
    return {'status':'saved_terminal_fit_report','fit_status':fit['status'],'fit_result':cfg['fit_result'],'initial_artifact':cfg.get('initial_artifact'),'selected_artifact':cfg.get('selected_artifact'),
            'initial':initial_summary,'selected':selected_summary,'selected_parameters':None if selected is None else selected['parameters'],
            'stopping_reason':fit.get('stopping_reason'),'error':fit.get('error'),'failed_or_unavailable_vectors':failures,'failed_or_stopped_requests':requests,
            'trajectory':path,'last_Jacobian_origin':lastJ,'parameter_bounds':bounds,'limits':fit['limits'],
            'counts':{'recorded_vectors':len(fit['vectors']),'new_parameter_vectors':fit.get('new_parameter_vectors'),'native_program_invocations':snapshot.get('native_program_invocations'),
                      'distinct_evaluator_vectors':snapshot.get('distinct_parameter_vectors'),'completed_outer_steps':sum(e['status']=='complete' for e in fit['iterations'])},
            'fit_cost':copy.deepcopy(fit['cost']),'job_id':fit.get('job_id'),'scheduler_accounting':None if cfg.get('scheduler_accounting') is None else {'source':cfg['scheduler_accounting'],'record':read(cfg['scheduler_accounting'])},
            'consumed_trimer_regression':regression(cfg,fit,contract),'new_native_calls':0,'fit_executed_by_report':False,
            'limitations':['Consumed development data; improved objective is not physical accuracy or fresh validation.',
                          'Response scales are engineering weights, not variances or physical tolerances.',
                          'Reference h=.001 energy curvature is fixed; model uses the qualified exact linear-response curvature.',
                          'Original provisional ECP-PBE0 reference and its unresolved method/state limitations are unchanged.',
                          'Eight trimer cycles are consumed regression outside the objective. Fresh radial, donor and protein outcomes remain excluded.',
                          'This report replays saved scoring, not raw native admission or a new independent force qualification.']}


def markdown(result):
    def number(x):return 'unavailable' if x is None else f'{x:.8g}'
    out=['# Bounded joint response fit report','',f"Fit status: **{result['fit_status']}**. Stop: `{result['stopping_reason']}`.",'']
    if result['error']:out+=['Recorded error: `'+result['error']+'`.','']
    s=result['selected'];initial=result['initial']
    if s is None:out+=['No frozen selected candidate is available. The terminal failure and unavailable vectors are retained in RESULT.json.','']
    else:
        ef=s['EF'];out += [f"Original gates: **{ef['energy_passed']}/10 energy; {ef['force_passed']}/15 force**. All25 pass: **{ef['all25_passed']}**.",
                           f"Paired W8→W9 model-minus-reference error: **{s['paired_coordination']['model_minus_reference_kcal_mol']:+.8f} kcal/mol**.",'',
                           '| Objective | Initial | Selected |','|---|---:|---:|',f"| Joint | {number(initial['objective'])} | {number(s['objective'])} |"]
        for name,v in s['block_objectives'].items():out.append(f"| {name} | {number(initial['block_objectives'][name])} | {number(v)} |")
        out+=['','## Four state susceptibilities','', '| State | Model | Reference | Error (a0³) |','|---|---:|---:|---:|']
        for name,r in s['chi'].items():out.append(f"| {name} | {number(r['model_a0_cubed'])} | {number(r['reference_a0_cubed'])} | {number(r['error_a0_cubed'])} |")
        out+=['','## Full-shell C dipole errors','', '| Parent | x | y | z (Debye) |','|---|---:|---:|---:|']
        for name,r in s['C'].items():out.append('| '+name+' | '+' | '.join(number(x) for x in r['error_Debye'])+' |')
        out+=['','## Literal pair dipoles','', '| Parent | Count | Component RMS | Maximum vector error (Debye) |','|---|---:|---:|---:|']
        for name,r in s['pair_contexts'].items():out.append(f"| {name} | {r['count']} | {number(r['component_RMS_error_Debye'])} | {number(r['maximum_vector_error_norm_Debye'])} |")
        out+=['','Failed original gates: '+(', '.join(ef['failed_ids']) or 'none')+'.','']
        contacts=[p['name'] for p in result['parameter_bounds'] if p['at_lower'] or p['at_upper']]
        out+=['Parameter bound contacts: '+(', '.join(contacts) or 'none')+'.','']
    out+=['## Convergence and cost','', '| Step | Status | Origin J | Selected J | Improved |','|---|---|---:|---:|---|']
    for r in result['trajectory']:out.append(f"| {r['outer_step']} | {r['status']} | {number(r['origin']['objective'])} | {number(None if r['selected'] is None else r['selected']['objective'])} | {r['actual_improvement']} |")
    cost=result['fit_cost'];counts=result['counts'];out+=['',f"Recorded vectors: {counts['recorded_vectors']}; native programs: {counts['native_program_invocations']}. "
         f"Fit wall: {number(cost.get('wall_seconds'))} s; parent CPU: {number(cost.get('parent_CPU_seconds'))} s; child CPU: {number(cost.get('child_CPU_seconds'))} s.",
         'Scheduler allocation accounting: '+('retained in RESULT.json.' if result['scheduler_accounting'] else 'not supplied; no allocation cost inferred.'),'',
         'Last Jacobian origin: `'+str(None if result['last_Jacobian_origin'] is None else result['last_Jacobian_origin']['key'])+'` (may differ from the selected point).','',
         f"Consumed trimer regression: **{result['consumed_trimer_regression']['status']}**, {result['consumed_trimer_regression']['complete_cycles']}/8 cycles.",'',
         '## Interpretation','']
    out += ['- '+s for s in result['limitations']]
    out += ['','All34 pair rows, all25 gate details, all28 descriptive field rows, parameter bounds, failures and exact source pins are retained in RESULT.json. Zero native calls were made by this report.','']
    return '\n'.join(out)


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    require(all(os.environ.get(k)=='1' for k in THREADS),'One local numerical thread')
    a.out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(a.config.read_text())
    try:result=collect(cfg)
    except Exception as exc:
        result={'status':'unavailable_saved_report','error':repr(exc),'configuration':pin(a.config),'new_native_calls':0,'report_cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu}}
        save(a.out/'RESULT.json',result);(a.out/'REPORT.md').write_text('# Saved fit report unavailable\n\n'+repr(exc)+'\n\nNo native calls; no values or passing result inferred.\n');raise
    result.update(configuration=pin(a.config),implementation=pin(__file__),report_cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1})
    save(a.out/'RESULT.json',result);(a.out/'REPORT.md').write_text(markdown(result));print(json.dumps({'result':pin(a.out/'RESULT.json'),'report':pin(a.out/'REPORT.md')}))


if __name__=='__main__':main()
