#!/usr/bin/env python3
"""Small terminal saved-score/Jacobian summary; no native or scoring imports."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time
import numpy as np

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
SOURCE='c3ee89ae069495e7513566f65c8f6f581da88da01b60da58c792020b182119b4'
OBJECTIVE='60274b54f5e1b71fed236818ca03e226ea46538720a5968100d81a973c286e80'
GEOMETRIES=('water_minus','water_central','water_plus','aquo8_central','aquo9_central')
SIZES=(9,9,9,51,57)

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']);return json.loads(Path(p['path']).read_text())
def close(a,b,label):require(np.allclose(a,b,rtol=1e-10,atol=1e-11),label)

def brief(score):
    energy=score['energy_controls'];force=score['per_fixture']
    require(sum(len(v) for v in energy.values())==10 and set(force)=={'La','Dy','paired'} and all(set(v)==set(GEOMETRIES) for v in force.values()),'Control denominator differs')
    for values in energy.values():
        for v in values.values():require(v['passed']==(abs(v['error_kcal_mol'])<=1),'Stored energy gate inconsistent')
    for values in force.values():
        for v in values.values():require(v['passed']==(v['rigid_RMS_kcal_mol_A']<=2 and v['maximum_body_component_error_kcal_mol_A']<=5),'Stored force gate inconsistent')
    require(score['development_checks_passed']==(all(v['passed'] for group in energy.values() for v in group.values()) and all(v['passed'] for group in force.values() for v in group.values())),'Stored full-control flag inconsistent')
    r=np.asarray(score['weighted_residuals']);require(r.shape==(413,) and np.isfinite(r).all(),'Invalid saved residual')
    close(float(r@r),score['objective'],'Stored loss/residual mismatch')
    return {'objective':score['objective'],'development_checks_passed':score['development_checks_passed'],
        'energy_controls':energy,'energy_passed_count':sum(v['passed'] for g in energy.values() for v in g.values()),'energy_total':10,
        'force_controls':{kind:{q:{k:v[k] for k in ('rigid_RMS_kcal_mol_A','maximum_body_component_error_kcal_mol_A','passed')} for q,v in group.items()} for kind,group in force.items()},
        'force_passed_count':sum(v['passed'] for g in force.values() for v in g.values()),'force_total':15,
        'development_partitions':score.get('development_partitions'),'weighted_block_losses':{k:.25*v for k,v in score['blocks'].items()}}

def decompose(r):
    r=np.asarray(r);require(r.shape==(413,) and np.isfinite(r).all(),'Residual shape/finiteness differs')
    result={'weighted_loss':float(r@r),'blocks':{'energy':float(r[:8]@r[:8])},'force_fixtures':{}};start=8
    for kind in ('La','Dy','paired'):
        result['force_fixtures'][kind]={};blockstart=start
        for q,size in zip(GEOMETRIES,SIZES):
            s=r[start:start+size];loss=float(s@s);result['force_fixtures'][kind][q]={'weighted_loss':loss,'linear_residual_rigid_RMS_kcal_mol_A':4*math.sqrt(5*loss)};start+=size
        s=r[blockstart:start];result['blocks'][kind]=float(s@s)
    require(start==413,'Residual layout differs');close(sum(result['blocks'].values()),result['weighted_loss'],'Partition closure differs');return result

def summarize(result):
    require(result['status']=='bounded_joint_coordination_water_complete' and result['error'] is None and result['candidate_frozen'],'Terminal successful17676 result required; no polling')
    require(result['job_id']=='17676','Wrong owned run')
    records=result['vectors'];selected=result['selected'];initial=result['initial'];params=result['family']['parameters'];require(len(params)==6,'Six parameters required')
    lo=np.array([p['lower'] for p in params]);hi=np.array([p['upper'] for p in params]);width=hi-lo
    require(np.all(width>0),'Invalid parameter intervals');values=np.array(selected['parameters'])
    require(values.shape==(6,) and np.isfinite(values).all() and np.all(values>=lo) and np.all(values<=hi),'Selected parameters invalid')
    bounds=[{'name':p['name'],'initial':v0,'selected':v,'lower':p['lower'],'upper':p['upper'],'normalized_move_from_initial':(v-v0)/w,
        'lower_contact':v==p['lower'],'upper_contact':v==p['upper']} for p,v0,v,w in zip(params,initial['parameters'],selected['parameters'],width)]
    steps=[]
    for entry in result['iterations']:
        require(entry['status']=='complete' and len(entry['derivative_keys'])==12 and len(entry['trial_keys'])==4,'Incomplete terminal batch')
        origin=records[entry['origin_key']];picked=records[entry['selected_key']]
        require(all(records[k]['status']=='complete' for k in entry['derivative_keys']+entry['trial_keys']),'Unavailable linked vector')
        best=origin
        for k in entry['trial_keys']:
            if records[k]['score']['objective']<best['score']['objective']:best=records[k]
        require(best['parameters']==picked['parameters'],'Current-first selection differs')
        delta=origin['score']['objective']-picked['score']['objective']
        steps.append({'outer_step':entry['outer_step'],'origin_loss':origin['score']['objective'],'selected_loss':picked['score']['objective'],
            'actual_improvement':delta,'selected_parameters':picked['parameters'],'all_controls_passed':picked['score']['development_checks_passed'],
            'passing_unselected_trials':[slot['id'] for slot,k in zip(entry['trial_slots'],entry['trial_keys']) if records[k]['score']['development_checks_passed'] and records[k]['parameters']!=picked['parameters']]})
    require(len(steps)==result['completed_outer_steps'],'Completed-step count differs')
    if steps:require(records[result['iterations'][-1]['selected_key']]['parameters']==selected['parameters'],'Terminal selection differs from last completed batch')
    out={'status':'saved_terminal_residual_summary_complete','stopping_reason':result['stopping_reason'],'completed_outer_steps':result['completed_outer_steps'],
        'initial':brief(initial['score']),'selected':brief(selected['score']),'parameters':bounds,'steps':steps,
        'total_actual_loss_improvement':initial['score']['objective']-selected['score']['objective'],
        'last_Jacobian':{'status':'unavailable','reason':'No recorded complete Jacobian','selected_point_projection_available':False},
        'recorded_cost':{k:result.get(k) for k in ('new_parameter_vectors','new_native_program_invocations','elapsed_seconds','parent_process_CPU_seconds','aggregate_native_child_user_CPU_seconds','aggregate_native_child_system_CPU_seconds')},
        'new_native_program_invocations':0,'new_model_evaluations':0,'raw_output_replay_performed':False,'objective_reevaluation_performed':False,
        'selected_candidate_changed':False,'physical_accuracy_qualified':False,'global_family_impossibility_claimed':False}
    entries=[x for x in result['iterations'] if 'Jacobian_normalized' in x]
    if not entries:return out
    entry=entries[-1];origin=records[entry['origin_key']];r=np.asarray(origin['score']['weighted_residuals']);J=np.asarray(entry['Jacobian_normalized']);require(J.shape==(413,6) and np.isfinite(J).all(),'Invalid saved Jacobian')
    columns=[]
    for j in range(6):
        a,b=(records[k] for k in entry['derivative_keys'][2*j:2*j+2]);span=b['normalized_parameters'][j]-a['normalized_parameters'][j]
        require(span>0 and span==entry['derivative_actual_widths'][j],'Derivative span differs')
        columns.append((np.asarray(b['score']['weighted_residuals'])-np.asarray(a['score']['weighted_residuals']))/span)
    require(np.column_stack(columns).tolist()==J.tolist(),'Saved Jacobian fails exact derivative reconstruction')
    direction,_,rank,singular=np.linalg.lstsq(J,-r,rcond=1e-6);projection=r+J@direction
    saved=entry['linear_directions'];require(int(rank)==saved['rank'],'Saved retained rank differs');close(singular,saved['singular_values'],'Saved singular values differ')
    close(direction,saved['unrestricted_GN_not_evaluated'],'Saved unconstrained direction differs')
    same=origin['parameters']==selected['parameters'];z=(np.asarray(origin['parameters'])-lo)/width;trials=[]
    for slot,k in zip(entry['trial_slots'],entry['trial_keys']):
        row=records[k];delta=(np.asarray(row['parameters'])-np.asarray(origin['parameters']))/width;close(delta,slot['normalized_delta'],'Actual trial displacement differs')
        require(np.max(np.abs(delta))<=.025+1e-12,'Trial exceeds trust box');pred=r+J@delta;predloss=float(pred@pred);close(predloss,slot['predicted_objective'],'Recorded trial prediction differs')
        trials.append({'id':slot['id'],'predicted_loss':predloss,'actual_loss':row['score']['objective'],'actual_improvement':origin['score']['objective']-row['score']['objective'],
            'actual_controls_passed':row['score']['development_checks_passed'],'linear_residual_error_norm':float(np.linalg.norm(np.asarray(row['score']['weighted_residuals'])-pred))})
    originparts=decompose(r)
    for kind,fixtures in originparts['force_fixtures'].items():
        for q,v in fixtures.items():close(v['linear_residual_rigid_RMS_kcal_mol_A'],origin['score']['per_fixture'][kind][q]['rigid_RMS_kcal_mol_A'],'Fixture scaling differs')
    out['last_Jacobian']={'status':'available','outer_step':entry['outer_step'],'origin_key':entry['origin_key'],'origin_parameters':origin['parameters'],
        'origin_equals_selected':same,'selected_point_projection_available':same,'projection_applies_to':'selected_point' if same else 'last_recorded_origin_only',
        'origin_actual_controls':brief(origin['score']),'rcond':1e-6,'absolute_singular_cutoff':float(singular[0]*1e-6),
        'rank':int(rank),'singular_values':singular.tolist(),'condition':float(singular[0]/singular[-1]) if singular[-1]>0 else None,
        'unrestricted_normalized_direction':direction.tolist(),'unrestricted_exceeds_trust':bool(np.max(np.abs(direction))>.025),
        'unrestricted_exceeds_parameter_box':bool(np.any(z+direction<0) or np.any(z+direction>1)),
        'origin_decomposition':originparts,'unrestricted_linear_projection':decompose(projection),'actual_trials':trials,
        'projection_is_prediction_only':True,'projected_control_pass_claimed':False,'global_family_impossibility_claimed':False}
    return out

def main(path,output):
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread required');started=time.monotonic();cpu=time.process_time();parent=pin(path);result=read(parent)
    require(result['implementation']['sha256']==SOURCE,'Wrong fitter source');cfg=read(result['config']);require(cfg['implementation']==result['implementation'] and cfg['objective']['sha256']==OBJECTIVE,'Scientific source differs')
    for k in ('implementation','objective','joint_math','specification','root_declaration'):require(pin(cfg[k]['path'])==cfg[k],'Changed source: '+k)
    summary=summarize(result);require(pin(path)==parent,'Parent mutated during read')
    summary.update(implementation=pin(__file__),specification=pin(Path(__file__).with_name('SPEC.md')),parent_result=parent,parent_result_unchanged=True,
        wall_seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,threads=1)
    out=Path(output);out.mkdir(parents=True,exist_ok=False);(out/'RESULT.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    selected=summary['selected'];lines=['# Terminal joint-water residual summary','',f"Stop: `{summary['stopping_reason']}` after {summary['completed_outer_steps']} completed steps. Selected checks: {selected['energy_passed_count']}/10 energy, {selected['force_passed_count']}/15 force.",'',f"Weighted loss: {summary['initial']['objective']:.9f} → {selected['objective']:.9f}.",'','| Energy control | Error (kcal/mol) | Pass |','|---|---:|---|']
    for group,rows in selected['energy_controls'].items():
        for q,v in rows.items():lines.append(f"| {group}/{q} | {v['error_kcal_mol']:.6f} | {v['passed']} |")
    lines+=['','| Force control | RMS | Maximum component | Pass |','|---|---:|---:|---|']
    for kind,rows in selected['force_controls'].items():
        for q,v in rows.items():lines.append(f"| {kind}/{q} | {v['rigid_RMS_kcal_mol_A']:.6f} | {v['maximum_body_component_error_kcal_mol_A']:.6f} | {v['passed']} |")
    j=summary['last_Jacobian']
    if j['status']=='available':
        condition='undefined' if j['condition'] is None else f"{j['condition']:.3f}"
        lines+=['',f"Last Jacobian origin equals selected: **{j['origin_equals_selected']}**. Rank {j['rank']}/6, condition {condition}. Its origin loss {j['origin_decomposition']['weighted_loss']:.9f} has unrestricted local linear residual loss {j['unrestricted_linear_projection']['weighted_loss']:.9f}. Projection applies to `{j['projection_applies_to']}` only; it is not an evaluated candidate or a global family limitation."]
    else:lines+=['','No recorded Jacobian is available; no span conclusion.']
    lines+=['','The JSON retains all parameter bounds, complete-step improvements, actual last trials, weighted residual partitions and recorded process/native costs. No raw native admission or objective recalculation was duplicated. No physical transfer, fit extension or candidate change follows.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'result':pin(out/'RESULT.json'),'report':pin(out/'REPORT.md')}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',required=True);p.add_argument('--output',required=True);a=p.parse_args();main(a.result,a.output)
