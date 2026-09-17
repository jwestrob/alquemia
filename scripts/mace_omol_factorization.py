"""Verify exact disconnected-component factorization from actual masked readouts."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL
from mace_omol_ablation_run import PROTOCOL,SEMANTICS,validate as validate_development,collect as collect_development
from mace_omol_readout import checked_arrays

FACTORIZATION='omol_charge_mask_disconnected_factorization_v1'
TOLERANCE=.01
ANCHOR='PQQ_1H4I'


@cached_file_checks
def calculate(collection):
    saved=read_json(collection);mp=verify(saved['manifest']);m=read_json(mp)
    validate_development(mp)
    if saved!=collect_development(mp) or not saved['canonical_extension_permitted']:
        raise InvalidArtifact('completed passing masked descriptor development required')
    reference={};reference_endpoints={};arrays={};checks=[]
    def check(name,error):
        checks.append({'name':name,'error_model_kcal':float(error),'pass':abs(error)<=TOLERANCE})
    detached=[t for t in m['tasks'] if t['kind']=='full' and t['position']=='detached']
    for t in detached:
        r=saved['rows'][t['task_id']]
        if r['input_state_check']['metal_neighbor_edge_count']!=0:
            raise InvalidArtifact('factorization requires an actually disconnected metal')
        a=checked_arrays(r,t);total=a['node_energy_eV']+a['embedding_energy_eV'];arrays[t['task_id']]=total
        if t['case_id']==ANCHOR and t['variant']=='primary':
            reference[t['metal']]=float(total[t['metal_index']])
            reference_endpoints[t['metal']]={'manifest':r['manifest'],'task_id':r['task_id'],
                                            'readout':r['native_readout'],'metal_index':t['metal_index']}
    if set(reference)!={'La','Ca'}:raise InvalidArtifact('fixed1H4I reference anchor missing')
    components=[];pairs={}
    for t in detached:
        total=arrays[t['task_id']];i=t['metal_index'];r=saved['rows'][t['task_id']]
        term=float(total[i]);error=(term-reference[t['metal']])*EV_TO_KCAL
        check(t['task_id']+'_constant_metal',error)
        closure=(math.fsum(map(float,total))-r['energy_eV'])*EV_TO_KCAL
        check(t['task_id']+'_readout_closure',closure)
        components.append({'task_id':t['task_id'],'case_id':t['case_id'],'metal':t['metal'],
                           'variant':t['variant'],'isolated_metal_model_eV':term,
                           'nonmetal_sum_model_eV':math.fsum(map(float,np.delete(total,i))),
                           'constant_error_model_kcal':error,'closure_error_model_kcal':closure})
        pairs.setdefault((t['case_id'],t['variant']),{})[t['metal']]=t
    paired=[]
    for (case,variant),pair in pairs.items():
        if set(pair)!={'La','Ca'}:raise InvalidArtifact('unpaired detached variant')
        ca,la=(pair[k] for k in ('Ca','La'));i=ca['metal_index']
        a,b=(xyz(verify(t['xyz'])) for t in (ca,la))
        if i!=la['metal_index'] or a[:i]+a[i+1:]!=b[:i]+b[i+1:]:
            raise InvalidArtifact('factorization requires identical paired nonmetal coordinates/species')
        delta=np.delete(arrays[ca['task_id']]-arrays[la['task_id']],i)*EV_TO_KCAL
        maximum=float(np.max(np.abs(delta)));sum_abs=math.fsum(map(float,np.abs(delta)))
        check(case+'_'+variant+'_nonmetal_max',maximum)
        check(case+'_'+variant+'_nonmetal_sum_abs',sum_abs)
        paired.append({'case_id':case,'variant':variant,'max_abs_nonmetal_change_model_kcal':maximum,
                       'sum_abs_nonmetal_change_model_kcal':sum_abs})
    offset=reference['Ca']-reference['La'];scores={}
    for case,old in saved['scores'].items():
        ca,la=(old['endpoints'][metal]['bound']['energy_eV'] for metal in ('Ca','La'))
        score=(ca-la-offset)*EV_TO_KCAL;error=score-old['R_mask_model_kcal']
        check(case+'_two_vs_four_forward_score',error)
        scores[case]={'two_forward_model_kcal':score,'four_forward_model_kcal':old['R_mask_model_kcal'],
                      'difference_model_kcal':error}
    sodium={t['metal']:saved['rows'][t['task_id']]['energy_eV'] for t in m['tasks']
            if t['ablation_group']=='spectator' and t['position']=='bound'}
    score=(sodium['Ca']-sodium['La']-offset)*EV_TO_KCAL
    error=score-saved['spectator']['R_mask_model_kcal'];check('GGR_sodium_two_vs_four_forward_score',error)
    return {'status':'pass' if all(c['pass'] for c in checks) else 'fail','factorization_id':FACTORIZATION,
            'protocol_id':PROTOCOL,'output_semantics':SEMANTICS,'model':m['model'],'software':m['software'],
            'source_collection':record(collection),'anchor_case':ANCHOR,'reference_endpoints':reference_endpoints,
            'disconnected_atom_model_eV':reference,'Ca_minus_La_disconnected_atom_model_eV':offset,
            'formula':'(T_bound_Ca-T_bound_La-(C_Ca-C_La))*eV_to_model_kcal',
            'conversion':EV_TO_KCAL,'tolerance_model_kcal':TOLERANCE,'components':components,
            'component_summation':'math.fsum_cross_driver_v1',
            'paired_nonmetal_checks':paired,'scores':scores,
            'sodium_score':{'two_forward_model_kcal':score,'difference_model_kcal':error},'checks':checks,
            'reference_interpretation':'model_readout_terms_not_quantum_ion_or_aquo_energies',
            'new_model_forwards':0,'new_DFT_calls':0,'calibration_changed':False,'baseline_changed':False}


@cached_file_checks
def report(collection,agreement,output):
    result=calculate(collection)
    result.update(agreement=record(agreement),implementation=record(__file__))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    lines=['# Disconnected-component factorization','',
           f"Status: {result['status']}; {len(result['checks'])} numerical checks; zero new model forwards.",
           'Modified learned model terms, not quantum ion energies, aquo energies or a fitted calibration.',
           f"Fixed anchor: {ANCHOR}; tolerance {TOLERANCE} model kcal.",'',
           '|Case|Four-call score|Two-call score|Difference, model kcal|','|---|---:|---:|---:|']
    for case,s in result['scores'].items():lines.append(f"|{case}|{s['four_forward_model_kcal']}|{s['two_forward_model_kcal']}|{s['difference_model_kcal']}|")
    lines+=['','The canonical four-call run is unchanged. Extend this equivalence check to its results when complete.',
            'No new accuracy claim, physical energy interpretation or automatic classification follows.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return {'status':result['status'],'result':record(out/'result.json')}


@cached_file_checks
def verified(path):
    saved=read_json(path)
    for key in ('agreement','implementation'):verify(saved[key])
    actual=calculate(verify(saved['source_collection']))
    if saved['status']!='pass' or {k:saved.get(k) for k in actual}!=actual:
        raise InvalidArtifact('factorization reference differs from actual qualified readouts')
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('collection','agreement','output'):p.add_argument('--'+key,required=True)
    print(json.dumps(report(**vars(p.parse_args())),indent=2))
