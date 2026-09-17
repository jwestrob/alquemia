"""Exact native-reference qualification of the fixed shared-neutral descriptor."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,cache_key,read_json,record,verify,write_new
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL,accepted_attempt
from mace_omol_neutral import ADAPTER,COMPONENT,PROTOCOL,RAW_SHA,INDEX
from mace_omol_gradients import ADAPTER as GRADIENT_ADAPTER
from mace_omol_response import CASES

STAGE='shared_neutral_core'
SETTINGS={'id':PROTOCOL,'charge_feature_adapter':ADAPTER,'conditioned_charge_category':0,
          'conditioned_table_index':INDEX,'raw_feature_sha256':RAW_SHA,
          'gradient_adapter':GRADIENT_ADAPTER,'chunk_size':1024,'energy_tolerance_model_kcal':.01,
          'gradient_tolerance_eV_A':.001,'partition_tolerance_kcal_scale':2.}


@cached_file_checks
def sources(response_report,native_collection,embedding_audit):
    from mace_omol import validate as validate_native, collect as collect_native, model as native_model
    from mace_omol_response import validate as validate_response,collect as collect_response
    r=read_json(response_report);mp=verify(r['manifest']);validate_response(mp);actual=collect_response(mp)
    if r['rows']!=actual['rows'] or r['status']!='complete' or r['gates']!=actual['gates'] or not r['gates']['analytic_derivative']:
        raise InvalidArtifact('actual numerically qualified response source required')
    m=read_json(mp);n=read_json(native_collection);npth=verify(n['manifest']);nm=read_json(npth)
    if (n['status']!='complete' or nm['stage']!='benchmark'
            or nm['model']!=native_model(verify(nm['software']))):
        raise InvalidArtifact('wrong native reference model/stage')
    for pin in nm['implementation'].values():verify(pin)
    # Revalidate the native qualification from exact arrays and its fixed gates.
    # Old qualified() also requires host-derived error summaries to be bit-exact.
    q=read_json(verify(nm['qualification']));qp=verify(q['manifest']);validate_native(qp)
    now=collect_native(qp);qm=read_json(qp)
    if (q['rows']!=now['rows'] or q['manifest']!=now['manifest'] or q['status']!='complete'
            or now['status']!='complete' or not q['numerical_gate_pass'] or not now['numerical_gate_pass']
            or nm['model']!=qm['model'] or nm['software']!=qm['software']):
        raise InvalidArtifact('native qualification raw artifacts or decisions differ')
    error_fields={'energy_error_kcal_mol','force_error_eV_A','error_kcal_mol'}
    decisions=lambda checks:[{k:v for k,v in c.items() if k not in error_fields} for c in checks]
    if decisions(q['checks'])!=decisions(now['checks']):raise InvalidArtifact('native qualification checks/decisions differ')
    if q['checks']!=now['checks']:
        print(json.dumps({'event':'native_reference_derived_report_roundoff',
                          'saved_checks':q['checks'],'actual_checks':now['checks']}),flush=True)
    if nm['software']!=m['software'] or nm['model']['checkpoint']!=m['model']['checkpoint']:
        raise InvalidArtifact('native reference model/software differ')
    audit=read_json(embedding_audit)
    for pin in (audit['checkpoint'],audit['arrays'],audit['native_embedding_source'],audit['implementation']):verify(pin)
    zero=next(row for row in audit['rows'] if row['charge']==0)
    if audit['checkpoint']!=m['model']['checkpoint'] or zero['raw_sha256']!=RAW_SHA or audit['specs']['total_charge']['offset']!=INDEX:
        raise InvalidArtifact('neutral category/feature differs from actual inspection')
    refs={};tasks=[]
    for t in m['tasks']:
        if t['point']!='center':continue
        task=copy.deepcopy(t)
        task.update(neutral_feature_experiment=PROTOCOL,descriptor_gradient_experiment=PROTOCOL,
                    charge_feature_adapter=ADAPTER,energy_component=COMPONENT)
        task.pop('cache_key');tasks.append(task)
        if t['metal']=='La':
            hits=[v for v in nm['tasks'] if v['xyz']['sha256']==t['xyz']['sha256']]
            if len(hits)!=1:raise InvalidArtifact('missing/ambiguous exact native charge-zero reference')
            old=hits[0];row=n['rows'][old['task_id']]
            payload={k:v for k,v in old.items() if k!='cache_key'}
            if old['cache_key']!=cache_key({'task':payload,'model':nm['model'],'software':nm['software'],'implementation':nm['implementation']}):
                raise InvalidArtifact('native endpoint cache identity differs')
            candidates=[accepted_attempt(a,old,npth) for a in (npth.parent/'execution'/old['task_id']).glob('attempt_*')]
            if row not in candidates:raise InvalidArtifact('native endpoint lacks its exact successful execution receipt')
            if (old['charge']!=0 or t['charge']!=0 or old['spin_multiplicity']!=t['spin_multiplicity']
                    or row.get('charge_feature_adapter') is not None or row.get('energy_only') or not row.get('forces')):
                raise InvalidArtifact('reference is not unmodified native charge-zero analytic output')
            verify(row['forces']);refs[t['case_id']]={'task_id':old['task_id'],'collection':record(native_collection),'result':row}
    if len(tasks)!=8 or set(refs)!=set(CASES):raise InvalidArtifact('eight-center/four-reference inventory differs')
    return m,tasks,refs


def model(software):
    from mace_omol import model as native_model
    r=native_model(software);r.update(energy_component=COMPONENT,output_semantics='energy_like_descriptor_not_quantum_endpoint',
                                    shared_neutral_settings=SETTINGS)
    return r


@cached_file_checks
def prepare(response_report,native_collection,embedding_audit,agreement,output):
    from mace_omol import common,seal
    parent,tasks,refs=sources(response_report,native_collection,embedding_audit)
    _,out,m=common(verify(parent['inventory']),verify(parent['software']),agreement,output,STAGE)
    if m['implementation']['mace_omol_gradients.py']['sha256']!=parent['implementation']['mace_omol_gradients.py']['sha256']:
        raise InvalidArtifact('qualified exact autograd adapter changed')
    m.update(protocol_id=PROTOCOL,model=model(verify(m['software'])),tasks=tasks,settings=SETTINGS,
             response_report=record(response_report),native_collection=record(native_collection),
             embedding_audit=record(embedding_audit),native_references=refs)
    return seal(out,m)


@cached_file_checks
def validate(manifest):
    m=read_json(manifest);p,tasks,refs=sources(verify(m['response_report']),verify(m['native_collection']),verify(m['embedding_audit']))
    if (m['stage']!=STAGE or m['protocol_id']!=PROTOCOL or m['model']!=model(verify(m['software']))
            or m['settings']!=SETTINGS or m['reused'] or m['native_references']!=refs or len(m['tasks'])!=8
            or m['software']!=p['software'] or m['inventory']!=p['inventory']):
        raise InvalidArtifact('neutral experiment/reference/configuration differs')
    for pin in (m['agreement'],*m['implementation'].values()):verify(pin)
    if m['implementation']['mace_omol_gradients.py']['sha256']!=p['implementation']['mace_omol_gradients.py']['sha256']:
        raise InvalidArtifact('qualified exact autograd adapter changed')
    for actual,wanted in zip(m['tasks'],tasks):
        payload={k:v for k,v in actual.items() if k!='cache_key'}
        if payload!=wanted:raise InvalidArtifact('neutral task differs from exact source')
        if actual['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('neutral scientific cache changed')
    return {'status':'pass','tasks':8,'actual_native_references':4,'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[];checks=[];cases={}
    for t in m['tasks']:
        good=[]
        for path in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(path,t,mp)
            if r is not None:good.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(path),'accepted':r is not None,
                             'receipt':record(path/'receipt.json') if (path/'receipt.json').exists() else None})
        rows[t['task_id']]=good[-1] if good else {'status':'unavailable','energy_eV':None}
    complete=all(r['status']=='computed' for r in rows.values());partition=None
    def check(name,error,tolerance):
        checks.append({'name':name,'error':float(error),'tolerance':tolerance,'pass':bool(abs(error)<=tolerance)})
    if complete:
        for key,row in rows.items():check(key+'_readout',row['native_readout']['component_sum_error_kcal_mol'],.01)
        tasks={t['task_id']:t for t in m['tasks']}
        for name in CASES:
            ref=m['native_references'][name]['result'];la=rows[name+'_center_La'];ca=rows[name+'_center_Ca']
            check(name+'_native_La_energy',(la['energy_eV']-ref['energy_eV'])*EV_TO_KCAL,.01)
            check(name+'_native_La_gradient',np.max(np.abs(np.load(verify(la['forces']))-np.load(verify(ref['forces'])))),.001)
            tc,tl=(tasks[name+'_center_'+metal] for metal in ('Ca','La'))
            d=(tc['DFT_source']['energy_hartree']-tl['DFT_source']['energy_hartree'])*HA_TO_KCAL
            learned=(ca['energy_eV']-la['energy_eV'])*EV_TO_KCAL
            jac=np.load(verify(tc['core_jacobians']));projection={}
            for metal in ('Ca','La'):
                t=tasks[name+'_center_'+metal];g=np.load(verify(rows[t['task_id']]['gradient']))*EV_TO_KCAL
                projection[metal]={'learned':[float(np.sum(j*g)) for j in jac],
                    'DFT':[float(np.sum(j*np.array(t['DFT_gradient']['gradient_kcal_mol_per_A']))) for j in jac],
                    'units':['kcal_scale_per_A','kcal_scale_per_radian']}
            cases[name]={'DFT_R_kcal_mol':d,'learned_R_model_kcal':learned,'core_anchor_R_kcal_scale':d-learned,
                         'projections':projection,'evidence':tc['evidence']}
        d=cases['GGR_connected']['DFT_R_kcal_mol']-cases['GGR_extended']['DFT_R_kcal_mol']
        t=cases['GGR_connected']['learned_R_model_kcal']-cases['GGR_extended']['learned_R_model_kcal']
        partition={'convention':'connected_minus_extended','DFT_shift_kcal_mol':d,'learned_shift_model_kcal':t,
                   'hybrid_shift_kcal_scale':d-t,'tolerance':2.,'pass':abs(d-t)<=2.}
    numerical=complete and all(c['pass'] for c in checks)
    return {'status':'complete' if complete else 'incomplete','manifest':record(mp),'protocol_id':PROTOCOL,
            'rows':rows,'attempts':attempts,'checks':checks,'numerical_gate_pass':numerical,'cases':cases,
            'partition':partition,'whole_inference_eligible':numerical and partition is not None and partition['pass'],
            'physical_charge_modified':False,'response_status':'response_model_not_validated',
            'relaxation_correction_kcal_mol':None,'calibrated_class':None,'predictive_improvement_claimed':False,'baseline_changed':False}


@cached_file_checks
def report(manifest,output):
    validate(manifest);r=collect(manifest);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    r['report_implementation']=record(__file__);write_new(out/'result.json',r);return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for k in ('response-report','native-collection','embedding-audit','agreement','output'):a.add_argument('--'+k,required=True)
    a=sub.add_parser('report');a.add_argument('--manifest',required=True);a.add_argument('--output',required=True)
    args=vars(p.parse_args());command=args.pop('command');r={'prepare':prepare,'report':report}[command](**args)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','attempts','checks','cases')},indent=2))
