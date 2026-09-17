"""Panel-wide factorization check and execution-specific research calibration."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL
from mace_omol_ablation_run import PROTOCOL,SEMANTICS,source_task
from mace_omol_ablation_panel import validate as panel_validate,collect as panel_collect
from mace_omol_panel_prepare import validate as prepared
from mace_omol_panel_report import summarize
from mace_omol_backbone_audit import inspect as inspect_backbone
from mace_omol_panel import DECISION
from mace_omol_factorization import verified as verified_factorization,FACTORIZATION,TOLERANCE
from mace_omol_readout import checked_arrays
from mace_global_prepare import FF

SCHEMA='alquemia.mace_mask_calibration.v1'


@cached_file_checks
def calculate(source_report,factorization):
    source=read_json(source_report);saved=read_json(verify(source['collection']))
    mp=verify(saved['manifest']);m=read_json(mp)
    if (source.get('protocol_id')!=PROTOCOL or source.get('output_semantics')!=SEMANTICS
            or m.get('stage')!='ablation_canonical'):
        raise InvalidArtifact('actual canonical masked-descriptor report required')
    verify(source['report_implementation']);verify(source['agreement'])
    panel_validate(mp);actual=panel_collect(mp)
    if saved!=actual or actual['status']!='complete' or not actual['numerical_gate_pass']:
        raise InvalidArtifact('complete actual canonical execution required')
    ref=verified_factorization(factorization)
    if ref['model']!=m['model'] or ref['software']!=m['software']:
        raise InvalidArtifact('factorization model/software differs from canonical run')
    data=prepared(verify(m['preparation']));integrity=inspect_backbone(data)
    endpoints={**actual['rows'],**{k:v['result'] for k,v in actual['reused_rows'].items()}}
    original=summarize(data,endpoints,integrity);recorded={s['case_id']:s for s in source['scores']}
    expected_bands=None
    if original['bands'] is not None:
        b=original['bands'];expected_bands={
            'Ca_max_inclusive_model_kcal':b['Ca_max_inclusive_kcal_mol'],
            'La_min_inclusive_model_kcal':b['La_min_inclusive_kcal_mol'],
            'scale':'R_mask','scope':PROTOCOL+'; canonical PQQ only','evidence_use':b['evidence_use']}
    if (source['bands']!=expected_bands or source['calibration_gap_model_kcal']!=original['calibration_gap_kcal_mol']
            or source['model']!=m['model'] or source['preparation']!=m['preparation']
            or source['agreement']!=m['agreement'] or source['decision_policy']!=DECISION):
        raise InvalidArtifact('released four-call calibration differs from actual fixed-rule algebra')
    checks=[];scores=[];offset=ref['Ca_minus_La_disconnected_atom_model_eV']
    def check(name,error):checks.append({'name':name,'error_model_kcal':float(error),'pass':abs(error)<=TOLERANCE})
    for row in original['scores']:
        name=row['case_id'];old=recorded[name]
        if any(old[k]!=row[k] for k in ('case_id','evaluation_role','expected_class','status','calibrated_class')):
            raise InvalidArtifact('canonical report evidence/state/decision differs from actual records')
        if old['R_mask_model_kcal']!=row['R_coord_kcal_mol']:
            raise InvalidArtifact('canonical report value differs from actual endpoint algebra')
        score={k:old.get(k) for k in ('case_id','evaluation_role','expected_class','sequence_accession_group',
                                     'evidence_stratum','prospectively_blind','preparation_status',
                                     'outside_reported_training_charge_range','outside_reported_training_size_range')}
        score.update(status=row['status'],four_call_model_kcal=row['R_coord_kcal_mol'],
                     two_call_model_kcal=None,calibrated_class=None,expected_direction_reached=False)
        if row['status']=='computed':
            ca,la=(row['endpoints'][metal]['bound']['energy_eV'] for metal in ('Ca','La'))
            value=(ca-la-offset)*EV_TO_KCAL;error=value-row['R_coord_kcal_mol']
            check(name+'_two_vs_four',error)
            score.update(two_call_model_kcal=value,difference_model_kcal=error,
                         bound_model_eV={'Ca':ca,'La':la},nonmetal_cancellation=None)
            nonmetal={};coordinates={};terms={}
            for metal in ('La','Ca'):
                r=row['endpoints'][metal]['detached'];t,_=source_task(r);a=checked_arrays(r,t)
                if r['input_state_check']['metal_neighbor_edge_count']!=0:
                    raise InvalidArtifact('canonical reference metal has graph edges')
                i=t['metal_index'];total=a['node_energy_eV']+a['embedding_energy_eV']
                terms[metal]=float(total[i]);check(name+'_'+metal+'_constant',
                    (terms[metal]-ref['disconnected_atom_model_eV'][metal])*EV_TO_KCAL)
                check(name+'_'+metal+'_closure',(math.fsum(map(float,total))-r['energy_eV'])*EV_TO_KCAL)
                nonmetal[metal]=np.delete(total,i)
                coords=xyz(verify(t['xyz']));coordinates[metal]=coords[:i]+coords[i+1:]
            if coordinates['Ca']!=coordinates['La']:raise InvalidArtifact('canonical paired nonmetal coordinates differ')
            delta=(nonmetal['Ca']-nonmetal['La'])*EV_TO_KCAL
            maximum=float(np.max(np.abs(delta)));sum_abs=math.fsum(map(float,np.abs(delta)))
            check(name+'_nonmetal_max',maximum);check(name+'_nonmetal_sum_abs',sum_abs)
            score.update(disconnected_atom_model_eV=terms,
                         nonmetal_cancellation={'max_abs_model_kcal':maximum,'sum_abs_model_kcal':sum_abs})
        scores.append(score)
    calibration=[s for s in scores if s['evaluation_role']=='calibration'];transfers=[s for s in scores if s['evaluation_role']!='calibration']
    if len(calibration)!=25 or len(transfers)!=3:raise InvalidArtifact('calibration/transfer inventory changed')
    numerical=bool(checks) and all(c['pass'] for c in checks);gap=None;bands=None
    if all(s['status']=='computed' for s in calibration):
        upper=max(s['two_call_model_kcal'] for s in calibration if s['expected_class']=='Ca')
        lower=min(s['two_call_model_kcal'] for s in calibration if s['expected_class']=='La');gap=lower-upper
        # Never let a roundoff crossing rescue a failed original calibration.
        if numerical and original['bands'] is not None and gap>DECISION['minimum_gap_kcal_mol']:
            bands={'Ca_max_inclusive_model_kcal':upper,'La_min_inclusive_model_kcal':lower,
                   'scale':'R_mask','scope':'canonical_PQQ_functional_class_only',
                   'score_evaluation':FACTORIZATION,'evidence_use':'consumed_calibration_not_independent_accuracy'}
    if bands:
        for s in scores:
            value=s['two_call_model_kcal']
            if value is None:continue
            s['calibrated_class']='Ca' if value<=bands['Ca_max_inclusive_model_kcal'] else (
                'La' if value>=bands['La_min_inclusive_model_kcal'] else 'inconclusive')
            s['expected_direction_reached']=s['calibrated_class']==s['expected_class']
    return {'schema_version':SCHEMA,'status':'pass' if numerical and bands else 'unavailable',
            'protocol_id':PROTOCOL,'output_semantics':SEMANTICS,'score_evaluation':FACTORIZATION,
            'source_report':record(source_report),'source_collection':source['collection'],
            'factorization':record(factorization),'model':m['model'],'software':m['software'],
            'preparation':m['preparation'],'decision_policy':DECISION,'scores':scores,'checks':checks,
            'factorization_gate_pass':numerical,'calibration_gap_model_kcal':gap,'bands':bands,
            'original_four_call_bands':source['bands'],
            'calibration_valid_count':sum(s['status']=='computed' for s in calibration),'calibration_total_count':25,
            'transfer_valid_count':sum(s['status']=='computed' for s in transfers),'transfer_total_count':3,
            'transfer_correct_count':sum(s['expected_direction_reached'] for s in transfers),
            'all_three_transfer_gate_pass':all(s['expected_direction_reached'] for s in transfers),
            'compatibility':{'assembly':'deposited_chain_A_with_declared_cofactor_and_site_waters',
                             'cofactor_charge_e':-3,'PQQ_atoms':27,'explicit_water_count':0,
                             'forcefield':record(FF),'spin_multiplicity':1},
            'scientific_scope':'PQQ_functional_metal_association; not generic_affinity_or_occupancy',
            'new_model_forwards':0,'new_DFT_calls':0,'baseline_changed':False,'production_promotion':False}


@cached_file_checks
def release(source_report,factorization,agreement,output):
    result=calculate(source_report,factorization)
    result.update(agreement=record(agreement),implementation=record(__file__))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'reference.json',result)
    lines=['# Execution-specific masked MACE research calibration','',
           f"Reference status: {result['status']}; factorization: {result['factorization_gate_pass']}.",
           f"Calibration: {result['calibration_valid_count']}/25; gap {result['calibration_gap_model_kcal']} model kcal.",
           f"Transfers: {result['transfer_valid_count']}/3 valid, {result['transfer_correct_count']}/3 in expected regions.",
           'Unsupported cases remain in the denominator. No new model forward or independent accuracy observation.','',
           '|Case|Expected|Four-call score|Two-call score|Decision|','|---|---|---:|---:|---|']
    for s in result['scores']:lines.append(f"|{s['case_id']}|{s['expected_class']}|{s['four_call_model_kcal']}|{s['two_call_model_kcal']}|{s['calibrated_class']}|")
    lines+=['','PQQ functional association only. These are learned model units, not binding free energies.',
            'All cases consumed; group related sequences/structures. Non-PQQ affinity has no band from this calibration.',
            'Baseline/default and original four-call results remain unchanged.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return {'status':result['status'],'reference':record(out/'reference.json')}


@cached_file_checks
def verified(path):
    saved=read_json(path)
    for key in ('agreement','implementation'):verify(saved[key])
    actual=calculate(verify(saved['source_report']),verify(saved['factorization']))
    if {k:saved.get(k) for k in actual}!=actual:
        raise InvalidArtifact('calibration differs from its actual source calculations')
    return saved


@cached_file_checks
def decision(score,manifest,calibration):
    """Keep valid raw scores when a calibration is unavailable or inapplicable."""
    ref=verified(calibration);p=read_json(verify(manifest['preparation']));policy=ref['compatibility']
    result={'calibration':record(calibration),'calibrated_class':None,'decision_status':'calibration_unavailable'}
    if ref['status']!='pass' or ref['bands'] is None:return result
    if (manifest['model']!=ref['model'] or manifest['software']!=ref['software']
            or score.get('score_evaluation')!=ref['score_evaluation']
            or manifest.get('factorization')!=ref['factorization']):
        return {**result,'decision_status':'incompatible_model_or_numeric_evaluation'}
    if (p['assembly']!=policy['assembly'] or p['cofactor_charge_e']!=policy['cofactor_charge_e']
            or sum(a['kind']=='frozen_PQQ' for a in p['physical_atoms'])!=policy['PQQ_atoms']
            or len(p['explicit_waters'])!=policy['explicit_water_count'] or p['forcefield']!=policy['forcefield']
            or any(e['spin_multiplicity']!=1 for e in p['endpoints'].values())):
        return {**result,'decision_status':'outside_canonical_PQQ_calibration_scope'}
    value=score['R_mask_model_kcal']
    if value is None or not score['numerical_gate_pass']:return {**result,'decision_status':'calculation_unavailable'}
    bands=ref['bands'];label='Ca' if value<=bands['Ca_max_inclusive_model_kcal'] else (
        'La' if value>=bands['La_min_inclusive_model_kcal'] else 'inconclusive')
    return {**result,'calibrated_class':label,'decision_status':'research_calibrated_PQQ_functional_association',
            'bands':bands,'classification_does_not_establish':'generic_affinity_or_physiological_occupancy'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('source-report','factorization','agreement','output'):p.add_argument('--'+key,required=True)
    print(json.dumps(release(**vars(p.parse_args())),indent=2))
