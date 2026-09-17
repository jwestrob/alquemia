"""Masked-descriptor response on exact archived DFT displacement geometries."""
from __future__ import annotations
import argparse
import copy
import json
import math
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, HA_TO_KCAL, cache_key, energy, read_json, record, verify, write_new, xyz
from affordable_response import extract
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL, accepted_attempt, check_atoms
from mace_omol_ablation import ADAPTER as MASK, COMPONENT
from mace_omol_ablation_run import PROTOCOL, SEMANTICS
from mace_omol_gradient_run import CONFIG as GRADIENT, check_parent, model as gradient_model

STAGE='masked_core_response'
ID='masked_omol_archived_DFT_response_screen_v1'
PREPARED_SHA='361d2ac85b93360c8560a761e0d76d622b8f314a7c238013f8fb7faa184ae99c'
POINTS=('center','s_minus','s_plus','theta_minus','theta_plus')
CASES=('GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9')
STEPS=(.02,math.pi/180)
TOL={'own_odd_absolute':.01,'own_odd_relative':.01,'even_absolute':.005,
     'even_relative':.25,'anchored_absolute':.02,'anchored_relative_to_even':.25,
     'direct_absolute':.02,'direct_relative':.25}


@cached_file_checks
def expected(prepared,dft):
    from mace_mechanics_run import preparation
    from mace_curvature import dft_source
    if record(prepared)['sha256']!=PREPARED_SHA:raise InvalidArtifact('exact declared mechanics preparation required')
    p=preparation(prepared);dc=read_json(dft);dm=read_json(verify(dc['manifest']))
    if dc['manifest']!=p['DFT_tasks'] or dc['status']!='complete':raise InvalidArtifact('wrong/incomplete DFT collection')
    old,om=dft_source(verify(p['DFT_reference']))
    sources={'new':{t['task_id']:t for t in dm['tasks']},'old':{t['task_id']:t for t in om['tasks']}}
    tasks=[]
    for name in CASES:
        c=read_json(verify(p['cases'][name]));jac=np.load(verify(c['core_jacobians']))
        for metal in ('Ca','La'):
            for point in POINTS:
                logical=f'{name}_{point}_{metal}';e=c['grids'][point]['endpoints'][metal]
                coords=xyz(verify(e['xyz']))
                if check_atoms(coords,e['charge'])!=e['state'] or coords[0][0]!=metal:
                    raise InvalidArtifact('archived state/atom order changed')
                if jac.shape!=(2,len(coords),3) or not np.isfinite(jac).all():raise InvalidArtifact('physical Jacobian shape changed')
                if logical in p['reused_GGR']:
                    source_id=p['reused_GGR'][logical]['source_task_id'];row=old['energies'][source_id]
                    st=sources['old'][source_id];collection=p['DFT_reference'];g=old['gradients'].get(source_id)
                else:
                    source_id=logical;row=dc['rows'][source_id];st=sources['new'][source_id]
                    collection=record(dft);g=dc['gradients'].get(source_id)
                if (st['xyz']['sha256']!=e['xyz']['sha256'] or st['charge']!=e['charge']
                        or st['multiplicity']!=e['spin_multiplicity'] or row['status']!='complete'
                        or energy(verify(row['output']))!=row['energy_hartree']):
                    raise InvalidArtifact('DFT comparison geometry/state/energy mismatch')
                verify(st['input']);receipt=read_json(verify(row['receipt']))
                for pin in receipt['artifacts'].values():verify(pin)
                gradient=None
                if point=='center':
                    if g is None:raise InvalidArtifact('actual analytic DFT gradient absent')
                    a=g['artifacts'];gradient=extract(verify(a['engrad']),verify(a['output']),verify(a['input']),verify(a['xyz']))
                    if a['xyz']['sha256']!=e['xyz']['sha256'] or gradient['gradient_kcal_mol_per_A']!=g['gradient_kcal_mol_per_A']:
                        raise InvalidArtifact('parsed DFT gradient/source differs')
                tasks.append({'task_id':logical,'case_id':name,'kind':'core','variant':'primary',
                    'point':point,'metal':metal,'metal_index':0,**copy.deepcopy(e),
                    'core_mapping':p['cases'][name],'core_jacobians':c['core_jacobians'],
                    'DFT_source':{'collection':collection,'task_id':source_id,**row},'DFT_gradient':gradient,
                    'evidence':c['evidence'],'evidence_use':'consumed_method_development',
                    'energy_only':point!='center','derivative_backend':'checkpointed',
                    'descriptor_gradient_experiment':ID,'energy_component':COMPONENT,
                    'charge_feature_adapter':MASK,'output_semantics':SEMANTICS,'capture_native_readout':True})
    if len(tasks)!=40:raise InvalidArtifact('response task inventory differs')
    return tasks


def model(software):
    r=gradient_model(software);r['response_screen']={'id':ID,'tolerances':TOL,'physical_steps':STEPS}
    # JSON roundtrip preserves the manifest's list representation.
    return json.loads(json.dumps(r))


@cached_file_checks
def prepare(prepared,dft,core_report,agreement,output):
    from mace_omol import common,seal
    q=read_json(core_report);qm=read_json(verify(q['manifest']))
    tasks=expected(prepared,dft)
    _,out,m=common(verify(qm['inventory']),verify(qm['software']),agreement,output,STAGE)
    m.update(protocol_id=PROTOCOL,model=model(verify(qm['software'])),tasks=tasks,
             prepared=record(prepared),DFT_collection=record(dft),core_report=record(core_report),
             gradient_settings=GRADIENT,response_settings={'id':ID,'tolerances':TOL,'points':list(POINTS)})
    check_parent(m)
    return seal(out,m)


@cached_file_checks
def validate(manifest):
    m=read_json(manifest);tasks=expected(verify(m['prepared']),verify(m['DFT_collection']))
    if (m['stage']!=STAGE or m['protocol_id']!=PROTOCOL or m['reused'] or m['gradient_settings']!=GRADIENT
            or m['response_settings']!={'id':ID,'tolerances':TOL,'points':list(POINTS)}
            or m['model']!=model(verify(m['software'])) or len(m['tasks'])!=len(tasks)):
        raise InvalidArtifact('response method or finite inventory changed')
    for pin in [m['agreement'],*m['implementation'].values()]:verify(pin)
    check_parent(m)
    q=read_json(verify(m['core_report']));qm=read_json(verify(q['manifest']))
    if m['software']!=qm['software'] or m['inventory']!=qm['inventory']:raise InvalidArtifact('qualification backend changed')
    for actual,wanted in zip(m['tasks'],tasks):
        payload={k:v for k,v in actual.items() if k!='cache_key'}
        if payload!=wanted:raise InvalidArtifact('response task differs from exact archived input')
        if actual['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('response cache identity differs')
    return {'status':'pass','tasks':len(tasks),'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[]
    for t in m['tasks']:
        good=[]
        for path in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(path,t,mp)
            if r is not None:good.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(path),'accepted':r is not None,
                             'receipt':record(path/'receipt.json') if (path/'receipt.json').exists() else None})
        rows[t['task_id']]=good[-1] if good else {'status':'unavailable','energy_eV':None}
    complete=all(v['status']=='computed' for v in rows.values());comparisons=[]
    if complete:
        tasks={t['task_id']:t for t in m['tasks']}
        for name in CASES:
            for axis,(label,h) in enumerate(zip(('s','theta'),STEPS)):
                values={}
                for metal in ('Ca','La'):
                    key=f'{name}_center_{metal}';t=tasks[key];jac=np.load(verify(t['core_jacobians']))
                    gradient=np.load(verify(rows[key]['gradient']))*EV_TO_KCAL
                    own=float(np.sum(jac[axis]*gradient))*h
                    dg=float(np.sum(jac[axis]*np.array(t['DFT_gradient']['gradient_kcal_mol_per_A'])))*h
                    energies={};dft_energies={}
                    for sign in ('minus','plus'):
                        other=f'{name}_{label}_{sign}_{metal}'
                        energies[sign]=(rows[other]['energy_eV']-rows[key]['energy_eV'])*EV_TO_KCAL
                        dft_energies[sign]=(tasks[other]['DFT_source']['energy_hartree']-t['DFT_source']['energy_hartree'])*HA_TO_KCAL
                    values[metal]={'own_linear':own,'DFT_linear':dg,'changes':energies,'DFT_changes':dft_energies}
                values['R']={k:({s:values['Ca'][k][s]-values['La'][k][s] for s in ('minus','plus')}
                                if isinstance(values['Ca'][k],dict) else values['Ca'][k]-values['La'][k]) for k in values['Ca']}
                for metal,v in values.items():
                    even=.5*sum(v['changes'].values());deven=.5*sum(v['DFT_changes'].values())
                    odd=.5*(v['changes']['plus']-v['changes']['minus'])
                    errors={s:v['changes'][s]-v['DFT_changes'][s] for s in ('minus','plus')}
                    anchored={s:sgn*v['DFT_linear']+even-v['DFT_changes'][s] for s,sgn in (('minus',-1),('plus',1))}
                    own_tol=max(TOL['own_odd_absolute'],TOL['own_odd_relative']*abs(v['own_linear']))
                    even_tol=max(TOL['even_absolute'],TOL['even_relative']*abs(deven))
                    anchored_tol=max(TOL['anchored_absolute'],TOL['anchored_relative_to_even']*abs(deven))
                    direct_tol={s:max(TOL['direct_absolute'],TOL['direct_relative']*abs(v['DFT_changes'][s])) for s in errors}
                    comparisons.append({'case_id':name,'direction':label,'metal':metal,'physical_amplitude':h,
                        'physical_unit':'A' if label=='s' else 'radian','model_energy_scale':'eV_equivalent_converted_once_to_model_kcal',
                        **v,'even_model_kcal':even,'DFT_even_kcal_mol':deven,'even_error':even-deven,
                        'own_odd_error':odd-v['own_linear'],'direct_errors':errors,'anchored_errors':anchored,
                        'tolerances':{'own_odd':own_tol,'even':even_tol,'anchored':anchored_tol,'direct':direct_tol},
                        'gates':{'analytic_derivative':abs(odd-v['own_linear'])<=own_tol,'curvature':abs(even-deven)<=even_tol,
                                 'DFT_anchored':all(abs(x)<=anchored_tol for x in anchored.values()),
                                 'direct_response':all(abs(errors[s])<=direct_tol[s] for s in errors)}})
    return {'status':'complete' if complete else 'incomplete','manifest':record(mp),'protocol_id':PROTOCOL,
            'experiment_id':ID,'rows':rows,'attempts':attempts,'comparisons':comparisons,
            'gates':{key:complete and all(c['gates'][key] for c in comparisons) for key in ('analytic_derivative','curvature','DFT_anchored','direct_response')},
            'response_status':'response_model_not_validated','relaxation_correction_kcal_mol':None,
            'entropy_correction_kcal_mol':None,'calibrated_class':None,'baseline_changed':False,'predictive_improvement_claimed':False}


@cached_file_checks
def report(manifest,output):
    validate(manifest);r=collect(manifest);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    r['report_implementation']=record(__file__);write_new(out/'result.json',r)
    return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for key in ('prepared','dft','core-report','agreement','output'):a.add_argument('--'+key,required=True)
    a=sub.add_parser('report');a.add_argument('--manifest',required=True);a.add_argument('--output',required=True)
    args=vars(p.parse_args());command=args.pop('command');r={'prepare':prepare,'report':report}[command](**args)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','attempts','comparisons')},indent=2))
