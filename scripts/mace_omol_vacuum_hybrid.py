"""Six intact-context calls for a separately declared vacuum DFT/MACE descriptor."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,cache_key,energy,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL,accepted_attempt,check_atoms
from mace_omol_ablation import ADAPTER,COMPONENT
from mace_omol_ablation_run import descriptor_model,SEMANTICS
from mace_omol_intact import energy_task,EVALUATION
from mace_omol_panel import ADAPTER_TASK

PROTOCOL='masked_omol_subtractive_context_DFT_vacuum_v1'
STAGE='vacuum_context'
VAC_SHA='be4ac74d2fdfcd942087ad6e0f4d7f3082609297eb748e8cccac0ebd44770a0a'
CONTEXT_SHA='cf41d2741fff33421fc4786dc6db6434260596f0cc8ca5b32ddaf86dbc69efa2'
PREP_SHA='361d2ac85b93360c8560a761e0d76d622b8f314a7c238013f8fb7faa184ae99c'
CASES=('GGR_1GLG','ALPHA_1F6S','ALPHA_6IP9')
TOL={'accounting_model_kcal':.01,'algebra_kcal_scale':1e-7,'partition_kcal_scale':2.,'ordering_kcal_scale':.02}
UNAVAILABLE={'aqueous_score':None,'calibrated_class':None,'predictive_improvement_claimed':False,
 'response_status':'response_model_not_validated','relaxation_correction_kcal_mol':None,'baseline_changed':False}


@cached_file_checks
def sources(vacuum_report,context_report):
    if record(vacuum_report)['sha256']!=VAC_SHA or record(context_report)['sha256']!=CONTEXT_SHA:
        raise InvalidArtifact('exact declared vacuum and context reports required')
    v=read_json(vacuum_report);c=read_json(context_report);vm=read_json(verify(v['manifest']))
    r=read_json(verify(c['response_report']));rm=read_json(verify(r['manifest']))
    prepared=verify(rm['prepared']);p=read_json(prepared)
    if record(prepared)['sha256']!=PREP_SHA or v['status']!='complete' or not all(x['pass'] for x in v['checks']):
        raise InvalidArtifact('complete matched vacuum diagnostic required')
    residual=v['partition']['models']['raw_zero']['vacuum_residual']
    if abs(residual)>TOL['partition_kcal_scale']:raise InvalidArtifact('vacuum partition prerequisite fails')
    for pin in [*vm['implementation'].values(),*rm['implementation'].values()]:verify(pin)
    for name,case in c['rows'].items():
        for metal,end in case['endpoints'].items():
            tid=f'{name}_center_{metal}';raw=v['rows'][tid];vr=read_json(verify(raw['receipt']))
            for pin in vr['artifacts'].values():verify(pin)
            if raw['energy_scope']!='isolated_vacuum_endpoint' or energy(verify(raw['output']))!=raw['energy_hartree']:
                raise InvalidArtifact('actual vacuum endpoint differs')
            task=next(t for t in vm['tasks'] if t['task_id']==tid)
            if task['xyz']['sha256']!=end['core_xyz']['sha256'] or task['core_jacobians']!=next(t for t in rm['tasks'] if t['task_id']==tid)['core_jacobians']:
                raise InvalidArtifact('vacuum/core source mapping differs')
            learned=r['rows'][tid];st=next(t for t in rm['tasks'] if t['task_id']==tid);mp=verify(r['manifest'])
            if learned!=end['masked_result'] or not any(accepted_attempt(a,st,mp)==learned for a in (mp.parent/'execution'/tid).glob('attempt_*')):
                raise InvalidArtifact('learned core lacks identical actual receipt')
        mapping=c['mapping_checks'][name]
        for k in ('physical_atoms','source_mapping'):verify(mapping[k])
        if not mapping['original_source_H_preserved'] or mapping['maximum_coordinate_error_A']>1e-6:
            raise InvalidArtifact('source physical/cap mapping is unsupported')
    tasks=[]
    for name in CASES:
        full=p['physical_cases'][name];physical=read_json(verify(full['physical_atoms']))
        prep=read_json(verify(full['original_global_preparation']));both={}
        for metal in ('Ca','La'):
            e=full['grids']['center']['endpoints'][metal];atoms=xyz(verify(e['xyz']));both[metal]=atoms
            idx=[i for i,a in enumerate(atoms) if a[0] in ('Ca','La')]
            if len(idx)!=1 or atoms[idx[0]][0]!=metal or check_atoms(atoms,e['charge'])!=e['state']:
                raise InvalidArtifact('full paired state differs')
            if len(atoms)!=len(physical):raise InvalidArtifact('physical atom inventory differs')
            for a,b in zip(atoms,physical):
                if a[0]!=(metal if b['element']=='M' else b['element']) or not np.allclose(a[1:],b['xyz_A'],atol=1e-9,rtol=0):
                    raise InvalidArtifact('original physical coordinates differ')
            for old in [case for case in c['rows'].values() if case['global_id']==name]:
                if old['endpoints'][metal]['full_xyz']!=e['xyz']:raise InvalidArtifact('full/core mapping source differs')
            tasks.append(energy_task({'task_id':f'{name}_{metal}_full','case_id':name,'kind':'full','position':'bound',
                'variant':'primary','metal':metal,'metal_index':idx[0],**copy.deepcopy(e),'source_xyz':e['xyz'],
                'physical_atoms':full['physical_atoms'],'original_global_preparation':full['original_global_preparation'],
                'assembly':full['assembly'],'microstate':full['microstate'],'explicit_waters':prep['explicit_waters'],
                'evidence':full['evidence'],'evidence_use':'consumed_method_development','energy_component':COMPONENT,
                'charge_feature_adapter':ADAPTER,'output_semantics':SEMANTICS,'edge_adapter':copy.deepcopy(ADAPTER_TASK)}))
        ca,la=both['Ca'],both['La'];i=idx[0]
        if (ca[:i]!=la[:i] or ca[i+1:]!=la[i+1:] or ca[i][1:]!=la[i][1:]
                or full['grids']['center']['endpoints']['La']['charge']-full['grids']['center']['endpoints']['Ca']['charge']!=1):
            raise InvalidArtifact('paired full coordinates/charge differ')
    return tasks,rm,v,c


def prepare(vacuum_report,context_report,agreement,output):
    from mace_omol import common,seal
    tasks,parent,v,c=sources(vacuum_report,context_report)
    _,out,m=common(verify(parent['inventory']),verify(parent['software']),agreement,output,STAGE)
    m.update(protocol_id=PROTOCOL,model=descriptor_model(verify(parent['software'])),tasks=tasks,
             vacuum_report=record(vacuum_report),context_report=record(context_report),hybrid_tolerances=TOL,
             energy_evaluation=EVALUATION,output_semantics=SEMANTICS,**UNAVAILABLE)
    return seal(out,m)


@cached_file_checks
def validate(manifest):
    m=read_json(manifest);wanted,parent,v,c=sources(verify(m['vacuum_report']),verify(m['context_report']))
    if (m['protocol_id']!=PROTOCOL or m['stage']!=STAGE or m['hybrid_tolerances']!=TOL or m['reused']
        or m['model']!=descriptor_model(verify(m['software'])) or m['software']!=parent['software']
        or m['inventory']!=parent['inventory'] or m['energy_evaluation']!=EVALUATION or len(m['tasks'])!=6):
        raise InvalidArtifact('vacuum context method or inventory changed')
    for pin in [m['agreement'],*m['implementation'].values()]:verify(pin)
    for name in ('mace_omol_ablation.py','mace_omol_products.py','mace_omol_edges.py','mace_omol_readout.py'):
        if m['implementation'][name]['sha256']!=parent['implementation'][name]['sha256']:
            raise InvalidArtifact('scientific descriptor adapter changed')
    for actual,w in zip(m['tasks'],wanted):
        payload={k:v for k,v in actual.items() if k!='cache_key'}
        if payload!=w or actual['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('exact original-H task or cache differs')
    return {'status':'pass','tasks':6,'new_DFT_calls':0,'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[]
    _,_,v,c=sources(verify(m['vacuum_report']),verify(m['context_report']))
    for t in m['tasks']:
        good=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:good.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,
                'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        rows[t['task_id']]=good[-1] if good else {'status':'unavailable','energy_eV':None}
    complete=all(r['status']=='computed' for r in rows.values());checks=[];cases={};contrasts=[];partition=None
    if complete:
        for tid,r in rows.items():
            error=r['native_readout']['component_sum_error_kcal_mol'];checks.append({'name':tid+'_readout','error':error,'pass':abs(error)<=TOL['accounting_model_kcal']})
        for name,old in c['rows'].items():
            endpoints={}
            for metal,core in old['endpoints'].items():
                tid=f'{name}_center_{metal}';d=v['rows'][tid]['energy_hartree'];low=core['masked_energy_eV'];full=rows[old['global_id']+'_'+metal+'_full']['energy_eV']
                endpoints[metal]={'DFT_vacuum_hartree':d,'DFT_CPCM_hartree':core['DFT_energy_hartree'],
                    'learned_core_eV':low,'learned_full_eV':full,'context_model_kcal':(full-low)*EV_TO_KCAL,
                    'hybrid_kcal_scale':d*HA_TO_KCAL+(full-low)*EV_TO_KCAL}
            ca,la=endpoints['Ca'],endpoints['La'];dr=(ca['DFT_vacuum_hartree']-la['DFT_vacuum_hartree'])*HA_TO_KCAL
            cr=(ca['learned_core_eV']-la['learned_core_eV'])*EV_TO_KCAL;fr=(ca['learned_full_eV']-la['learned_full_eV'])*EV_TO_KCAL
            hr=dr+fr-cr;error=hr-(ca['hybrid_kcal_scale']-la['hybrid_kcal_scale'])
            checks.append({'name':name+'_component_algebra','error':error,'pass':abs(error)<=TOL['algebra_kcal_scale']})
            cases[name]={'endpoints':endpoints,'DFT_CPCM_R_kcal_mol':old['DFT_core_R_kcal_mol'],'DFT_vacuum_R_kcal_mol':dr,
                'learned_core_R_model_kcal':cr,'learned_full_R_model_kcal':fr,'context_R_model_kcal':fr-cr,
                'hybrid_R_kcal_scale':hr,'evidence':old['evidence'],'global_id':old['global_id']}
        partition=cases['GGR_connected']['hybrid_R_kcal_scale']-cases['GGR_extended']['hybrid_R_kcal_scale']
        error=partition-v['partition']['models']['raw_zero']['vacuum_residual']
        checks.append({'name':'partition_full_term_cancellation','error':error,'pass':abs(error)<=TOL['algebra_kcal_scale']})
        for a in ('ALPHA_1F6S','ALPHA_6IP9'):
            for g in ('GGR_extended','GGR_connected'):
                delta=cases[a]['hybrid_R_kcal_scale']-cases[g]['hybrid_R_kcal_scale']
                contrasts.append({'alpha':a,'GGR':g,'difference_kcal_scale':delta,'pass':delta>TOL['ordering_kcal_scale']})
    numerical=complete and bool(checks) and all(x['pass'] for x in checks)
    return {'manifest':record(mp),'protocol_id':PROTOCOL,'status':'complete' if complete else 'incomplete',
        'rows':rows,'attempts':attempts,'cases':cases,'checks':checks,'numerical_gate_pass':numerical,
        'partition_shift_kcal_scale':partition,'partition_gate_pass':partition is not None and abs(partition)<=TOL['partition_kcal_scale'],
        'contrasts':contrasts,'ordering_gate_pass':numerical and len(contrasts)==4 and all(x['pass'] for x in contrasts),
        'independent_biological_groups':2,'tolerances':TOL,**UNAVAILABLE}


def report(manifest,output):
    r=collect(manifest);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);r['report_implementation']=record(__file__)
    write_new(out/'result.json',r);return {k:r[k] for k in ('status','numerical_gate_pass','partition_shift_kcal_scale','partition_gate_pass','contrasts','ordering_gate_pass')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('vacuum-report','context-report','agreement','output'):a.add_argument('--'+k,required=True)
    a=s.add_parser('report');a.add_argument('--manifest',required=True);a.add_argument('--output',required=True)
    a=p.parse_args();r=prepare(a.vacuum_report,a.context_report,a.agreement,a.output) if a.op=='prepare' else report(a.manifest,a.output)
    print(json.dumps(r,indent=2))
