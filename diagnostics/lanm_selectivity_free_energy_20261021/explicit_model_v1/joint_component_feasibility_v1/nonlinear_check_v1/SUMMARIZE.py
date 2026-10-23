#!/usr/bin/env python3
"""Summarize the one declared native candidate; no engines or parameter selection."""
import argparse, hashlib, json, math, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
BASE=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1'

def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def read(p):
    if isinstance(p,dict):
        assert pin(p['path'])==p, 'Changed pinned input'
        p=p['path']
    return json.loads(Path(p).read_text())

def main(result_path,out):
    wall,cpu=time.monotonic(),time.process_time()
    decl=read(HERE/'ROOT_DECLARATION.json'); cert=read(decl['half_certificate'])
    failed=read(decl['failed_proposal']); native=read(result_path)
    cfg=read(native['configuration']); saved=read(BASE/'actual_v1/RESULT.json')
    admitted=read(BASE/'actual_v1/ROOT_REVIEW.json')
    # The diagnostic's immutable result is admitted independently; locate its exact pin.
    assert admitted['status']=='admitted_saved_local_component_feasibility_diagnostic'
    assert admitted['result']==pin(BASE/'actual_v1/RESULT.json')
    assert cfg['declaration']==pin(HERE/'ROOT_DECLARATION.json')
    assert cfg['component_diagnostic']==pin(BASE/'actual_v1/RESULT.json')
    assert native['implementation']==cfg['implementation']
    assert decl['normalized_direction']==cert['normalized_half_direction']
    assert saved['status']=='complete_saved_joint_component_feasibility'
    assert pin(BASE/'actual_v1/RESULT.json')['sha256']=='ce2fd33896b823c53b43eca4c49813b56fe1ba547bf774adb02602b873cd1dc0'
    assert native['status']=='complete_one_fixed_half_native_pilot', 'Only a completed native point can be compared'
    assert cfg['candidate']['parameters']==decl['physical_parameters']==cert['physical_half_parameters']
    assert cfg['candidate']['key']==cert['half_parameter_token']
    assert cfg['origin']==decl['origin_artifact'] and cfg['half_review']==decl['half_certificate']
    assert failed['status']=='failed_numerical_local_proposal' and failed['solver']['success'] is False
    origin=read(decl['origin_artifact']); ans=native['answer']; comp=native['components']
    derivative=next(r for r in saved['origins'] if r['key']==decl['origin'])
    assert derivative['Jacobian_components_order']==saved['component_order']==['P8','B8','P9','B9']
    assert len(derivative['Jacobian_components'])==4
    assert all(len(row)==12 for row in derivative['Jacobian_components'])
    assert len(decl['normalized_direction'])==12
    for part in ('origin','candidate'):
        for n in ('8','9'):
            assert comp[part]['paired']['W'+n]['reference']==saved['reference_targets'][n]
    errors=saved['vectors'][decl['origin']]['component_errors_kcal_mol']
    predictions={label:errors[label]+math.fsum(j*d for j,d in zip(row,decl['normalized_direction']))
                 for label,row in zip(saved['component_order'],derivative['Jacobian_components'])}
    channels={}
    for name in saved['component_order']:
        term,n=name[0],name[1:]; old=comp['origin']['paired']['W'+n]['error'][term]
        actual=comp['candidate']['paired']['W'+n]['error'][term]
        assert abs(old-errors[name])<1e-9
        channels[name]={'origin_error_kcal_mol':old,'predicted_error_kcal_mol':predictions[name],
                        'actual_error_kcal_mol':actual,'actual_minus_predicted_kcal_mol':actual-predictions[name]}
    blocks={}
    for key,score_key in [('old_EF','old_EF'),('pair_mu','pair'),('C_mu','C'),('chi','chi')]:
        old=origin['score']['block_objectives'][score_key];new=ans['score']['block_objectives'][score_key]
        blocks[key]={'origin_squared_norm_unweighted':old,'predicted_squared_norm_unweighted':failed['half']['blocks'][key]['predicted_score'],
                     'actual_squared_norm_unweighted':new,'actual_nonincrease_exact':new<=old,'change':new-old}
    for key,term in [('P_energy','P'),('B_energy','B')]:
        old=comp['origin']['J_'+term];new=comp['candidate']['J_'+term]
        blocks[key]={'origin_RMS_kcal_mol':math.sqrt(old),'predicted_RMS_kcal_mol':failed['half']['blocks'][key]['predicted_norm'],
                     'actual_RMS_kcal_mol':math.sqrt(new),'actual_squared_norm':new,'actual_strict_improvement':new<old,
                     'fractional_RMS_improvement':1-math.sqrt(new/old)}
    old=ans['score']['old_EF'];gates={}
    for group,entries in old['energy_controls'].items():
        for name,row in entries.items():gates['energy/'+group+'/'+name]=row
    for group,entries in old['per_fixture'].items():
        for name,row in entries.items():
            gates['force/'+group+'/'+name]={k:row[k] for k in ['rigid_RMS_kcal_mol_A','maximum_body_component_error_kcal_mol_A','passed']}
    assert len(gates)==25 and all(r['passed'] is True for r in gates.values())==ans['score']['original_EF_checks_passed']
    both=all(blocks[k]['actual_strict_improvement'] for k in ('P_energy','B_energy'))
    old_preserved=all(blocks[k]['actual_nonincrease_exact'] for k in ('old_EF','pair_mu','C_mu','chi'))
    gates_pass=all(r['passed'] is True for r in gates.values())
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    summary={'status':'complete_consumed_native_prediction_comparison','implementation':pin(__file__),
             'declaration':pin(HERE/'ROOT_DECLARATION.json'),'native_result':pin(result_path),'saved_diagnostic':pin(BASE/'actual_v1/RESULT.json'),
             'half_review':decl['half_certificate'],'failed_optimizer':decl['failed_proposal'],'original_optimizer_passed':False,
             'blocks':blocks,'signed_component_errors':channels,'all25':gates,'all25_pass':gates_pass,
             'both_components_improve':both,'all_old_blocks_nonincrease_exact':old_preserved,
             'declared_joint_native_outcome_passes':both and old_preserved and gates_pass,
             'fresh_validation':False,'physical_model_qualified':False,'native_cost':native['cost'],
             'cost':{'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'new_native_calls':0}}
    (out/'RESULT.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:summary[k] for k in ['status','all25_pass','both_components_improve','all_old_blocks_nonincrease_exact','declared_joint_native_outcome_passes','blocks']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',required=True);p.add_argument('--out',required=True);a=p.parse_args();main(a.result,a.out)
