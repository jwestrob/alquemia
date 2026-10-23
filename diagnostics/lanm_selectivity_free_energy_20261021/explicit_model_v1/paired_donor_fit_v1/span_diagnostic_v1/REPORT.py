#!/usr/bin/env python3
"""Terminal saved-array diagnostic only; no evaluator or native imports."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
import numpy as np

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
SOURCE='9559a2692042781268567933f805a9e2827ba6e48352e00f12e79bddc4ce4c76'
OBJECTIVE='86e85fcb0537bdf5b459874eeda0c8bb68aead9758c712014b18b2db94757974'
BODY='bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819'

def require(ok,message):
    if not ok:raise ValueError(message)

def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def read(record):
    require(pin(record['path'])==record,'Changed authority: '+record['path'])
    return json.loads(Path(record['path']).read_text())

def module(record,name):
    require(pin(record['path'])==record,'Changed source')
    s=importlib.util.spec_from_file_location(name,record['path']);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def close(a,b,label):
    require(np.allclose(a,b,rtol=1e-10,atol=1e-11),label)

def layout(score,data,model,obj,body):
    segments=[];reconstructed=[]
    def add(values,**labels):
        start=len(reconstructed);reconstructed.extend(values)
        segments.append(dict(start=start,stop=len(reconstructed),**labels))
    for group,weight in obj.ENERGY_WEIGHTS.items():
        values=score['energy_errors_kcal_mol'][group]
        for name,value in values.items():
            env=('mixed_wet_minus_dry' if group=='water_addition' else
                 'dry' if name.startswith('matched_dry_') else 'wet')
            add([value*math.sqrt(.25*weight/len(values))],block='energy',subblock=group,fixture=name,environment=env)
    require(len(reconstructed)==27,'Energy residual denominator differs')
    for kind in ('La','Dy','paired'):
        for q in obj.DONOR:
            row=data['rows'][q+'_La']
            if kind=='paired':
                a=body.subtract(model[q+'_Dy']['forces'],model[q+'_La']['forces'])
                b=body.subtract(data['rows'][q+'_Dy']['reference']['forces'],row['reference']['forces'])
            else:a,b=model[q+'_'+kind]['forces'],data['rows'][q+'_'+kind]['reference']['forces']
            weighted,_=body.rigid_residual(row,a,b)
            add([v*math.sqrt(.25/len(obj.DONOR)) for v in weighted],block=kind,subblock='rigid',fixture=q,
                environment='dry' if q.startswith('matched_dry_') else 'wet')
    close(reconstructed,score['weighted_residuals'],'Objective residual ordering/scaling differs')
    return segments

def decompose(vector,segments):
    vector=np.asarray(vector);out={'loss':float(vector@vector),'blocks':{},'energy_subblocks':{},'environments':{},'force_fixtures':{}}
    covered=[]
    for seg in segments:
        a=vector[seg['start']:seg['stop']];loss=float(a@a);block=seg['block'];env=seg['environment']
        covered.extend(range(seg['start'],seg['stop']))
        out['blocks'][block]=out['blocks'].get(block,0.)+loss
        out['environments'][env]=out['environments'].get(env,0.)+loss
        if block=='energy':
            name=seg['subblock'];out['energy_subblocks'][name]=out['energy_subblocks'].get(name,0.)+loss
        else:
            out['force_fixtures'].setdefault(block,{})[seg['fixture']]={'weighted_loss':loss,
                'rigid_RMS_from_linear_residual_kcal_mol_A':4*math.sqrt(12)*math.sqrt(loss),'environment':env}
    require(covered==list(range(len(vector))),'Incomplete or overlapping residual layout')
    close(sum(out['blocks'].values()),out['loss'],'Block loss closure')
    close(sum(out['environments'].values()),out['loss'],'Environment loss closure')
    return out

def brief(score):
    return {'objective':score['objective'],'development_checks_passed':score['development_checks_passed'],
        'energy_checks_passed':score['energy_checks_passed'],
        'per_fixture':{kind:{q:{k:v[k] for k in ('rigid_RMS_kcal_mol_A','maximum_body_component_error_kcal_mol_A','passed')}
                           for q,v in rows.items()} for kind,rows in score['per_fixture'].items()}}

def main(path,output):
    started=time.monotonic();cpu=time.process_time()
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread required')
    source=pin(path);result=read(source)
    require(result['status']=='bounded_paired_donor_complete' and result['error'] is None,'Terminal successful result required; no polling')
    require(result['implementation']['sha256']==SOURCE,'Different donor implementation')
    cfg=read(result['config']);require(cfg['implementation']==result['implementation'],'Implementation/config mismatch')
    require(cfg['objective']['sha256']==OBJECTIVE and cfg['body_helper']['sha256']==BODY,'Objective/helper changed')
    require(cfg['joint_math']['sha256']=='741167a8c5ce3c9e07af6a94baf40894c277b5c47950f31b3c526da125ae11ca','Direction algorithm changed')
    for key in ('implementation','objective','body_helper','joint_math','specification','root_declaration'):
        require(pin(cfg[key]['path'])==cfg[key],'Source pin differs: '+key)
    data=read(cfg['training_result']);obj=module(cfg['objective'],'span_objective');body=module(cfg['body_helper'],'span_body')
    iterations=[e for e in result['iterations'] if 'Jacobian_normalized' in e]
    require(iterations,'No measured Jacobian available')
    entry=iterations[-1];require(entry['status']=='complete','Last Jacobian trial batch incomplete')
    records=result['vectors'];origin=records[entry['origin_key']];selected=result['selected']
    water={gid:result['initial']['model'][gid] for gid in obj.WATER_CELLS}
    keys=list(dict.fromkeys([entry['origin_key']]+entry['derivative_keys']+entry['trial_keys']))
    for k in keys:
        row=records[k];require(row['status']=='complete','Unavailable linked vector')
        require(obj.evaluate(data,row['model'],body,water)==row['score'],'Linked score replay differs: '+k)
    require(obj.evaluate(data,selected['model'],body,water)==selected['score'],'Selected score replay differs')
    params=result['family']['parameters'];require(len(params)==8,'Eight parameters required')
    lo=np.array([p['lower'] for p in params]);hi=np.array([p['upper'] for p in params]);width=hi-lo
    center=np.array(origin['parameters']);z=(center-lo)/width;r=np.array(origin['score']['weighted_residuals'])
    close(z,origin['normalized_parameters'],'Origin normalization differs')
    derivative=[records[k] for k in entry['derivative_keys']];require(len(derivative)==16,'Sixteen derivative slots required')
    columns=[];actual_widths=[]
    for j in range(8):
        a,b=derivative[2*j:2*j+2]
        for row,sign in ((a,-1),(b,1)):
            target=z.copy();target[j]=np.clip(target[j]+sign*.005,0.,1.)
            expected=center.copy();expected[j]=np.clip(center[j]+(target[j]-z[j])*width[j],lo[j],hi[j])
            require(row['parameters']==expected.tolist(),'Probe differs from declared normalized step')
            close(row['normalized_parameters'],(np.asarray(row['parameters'])-lo)/width,'Probe normalization differs')
        span=b['normalized_parameters'][j]-a['normalized_parameters'][j];require(span>0,'Collapsed probe interval');actual_widths.append(span)
        columns.append((np.asarray(b['score']['weighted_residuals'])-np.asarray(a['score']['weighted_residuals']))/span)
    J=np.column_stack(columns);require(J.tolist()==entry['Jacobian_normalized'],'Saved Jacobian does not exactly replay')
    require(actual_widths==entry['derivative_actual_widths'],'Probe widths differ')
    full,_,rank,singular=np.linalg.lstsq(J,-r,rcond=1e-6);saved=entry['linear_directions']
    require(int(rank)==saved['rank'],'Retained rank differs');close(singular,saved['singular_values'],'Singular values differ')
    close(full,saved['unrestricted_GN_not_evaluated'],'Unrestricted direction differs')
    projected=r+J@full;segments=layout(origin['score'],data,origin['model'],obj,body)
    original=decompose(r,segments);projection=decompose(projected,segments)
    for kind,value in origin['score']['blocks'].items():close(original['blocks'][kind],.25*value,'Block weighting differs')
    for kind,fixtures in original['force_fixtures'].items():
        for q,value in fixtures.items():close(value['rigid_RMS_from_linear_residual_kcal_mol_A'],origin['score']['per_fixture'][kind][q]['rigid_RMS_kcal_mol_A'],'Fixture scaling differs')
    same_origin=origin['parameters']==selected['parameters']
    lower,upper=np.maximum(-.025,-z),np.minimum(.025,1-z)
    directions={}
    for name in ('scaled_GN','bounded_LS'):
        d=np.array(saved[name]);require(np.all(d>=lower-1e-12) and np.all(d<=upper+1e-12),'Saved direction violates trust/bounds')
        directions[name]={'normalized_direction':d.tolist(),'linear_prediction':decompose(r+J@d,segments)}
    trials=[];best=origin
    require(len(entry['trial_keys'])==len(entry['trial_slots'])==4,'Complete four-trial batch required')
    for slot,key in zip(entry['trial_slots'],entry['trial_keys']):
        row=records[key];require(row['parameters']==slot['parameters'],'Trial slot mismatch')
        delta=(np.array(row['parameters'])-center)/width;close(delta,slot['normalized_delta'],'Trial displacement differs')
        require(np.max(np.abs(delta))<=.025+1e-12,'Actual trial trust limit exceeded')
        prediction=r+J@delta;actual=np.array(row['score']['weighted_residuals']);pred_loss=float(prediction@prediction)
        close(pred_loss,slot['predicted_objective'],'Saved trial prediction differs')
        trials.append({'id':slot['id'],'key':key,'normalized_delta':delta.tolist(),'predicted_loss':pred_loss,
            'actual_loss':row['score']['objective'],'predicted_improvement':original['loss']-pred_loss,
            'actual_improvement':original['loss']-row['score']['objective'],
            'weighted_residual_linearization_error_norm':float(np.linalg.norm(actual-prediction)),
            'actual_minus_predicted_loss':row['score']['objective']-pred_loss,'actual_controls_passed':row['score']['development_checks_passed'],
            'predicted_decomposition':decompose(prediction,segments),'actual_decomposition':decompose(actual,segments)})
        if row['score']['objective']<best['score']['objective']:best=row
    require(best['parameters']==records[entry['selected_key']]['parameters']==selected['parameters'],'Original complete-batch selection differs')
    probes=[]
    for j,row in enumerate(derivative):
        delta=(np.array(row['parameters'])-center)/width;prediction=r+J@delta;actual=np.array(row['score']['weighted_residuals'])
        probes.append({'parameter':params[j//2]['name'],'side':'minus' if j%2==0 else 'plus',
            'actual_loss':row['score']['objective'],'predicted_loss':float(prediction@prediction),
            'weighted_residual_linearization_error_norm':float(np.linalg.norm(actual-prediction))})
    summary={'status':'saved_span_diagnostic_complete','implementation':pin(__file__),'specification':pin(Path(__file__).with_name('SPEC.md')),
        'parent_result':source,'parent_config':result['config'],'objective':cfg['objective'],'body_helper':cfg['body_helper'],
        'new_native_program_invocations':0,'selected_candidate_changed':False,'physical_accuracy_qualified':False,'reserved_outcomes_opened':False,
        'last_measured_outer_step':entry['outer_step'],'origin_key':entry['origin_key'],'origin_parameters':origin['parameters'],
        'selected_parameters':selected['parameters'],'Jacobian_origin_equals_selected':same_origin,
        'projection_applies_to':'final_selected_point' if same_origin else 'last_measured_origin_only',
        'origin_actual_controls':brief(origin['score']),'selected_actual_controls':brief(selected['score']),
        'weighted_residual_count':len(r),'parameter_names':[p['name'] for p in params],
        'rcond':1e-6,'absolute_singular_cutoff':float(singular[0]*1e-6),'retained_rank':int(rank),'singular_values':singular.tolist(),
        'full_condition':float(singular[0]/singular[-1]) if singular[-1]>0 else None,
        'unrestricted_normalized_direction':full.tolist(),'unrestricted_exceeds_trust':bool(np.any(full<lower) or np.any(full>upper)),
        'unrestricted_exceeds_parameter_box':bool(np.any(z+full<0) or np.any(z+full>1)),
        'unrestricted_direction_evaluated':False,'origin_decomposition':original,'unrestricted_linear_residual_decomposition':projection,
        'locally_projected_loss_fraction':1-projection['loss']/original['loss'] if original['loss'] else 0.,
        'saved_trust_directions':directions,'actual_trials':trials,'derivative_linearity':probes,'residual_layout':segments,
        'parent_stopping_reason':result['stopping_reason'],'parent_completed_outer_steps':result['completed_outer_steps'],
        'parent_cost':{k:result[k] for k in ('new_parameter_vectors','new_native_program_invocations','elapsed_seconds','parent_process_CPU_seconds','aggregate_native_child_user_CPU_seconds','aggregate_native_child_system_CPU_seconds')},
        'limitations':['Local retained-span projection is not a global family impossibility result.',
            'Unrestricted and projected responses are linear predictions, not actual new evaluations or accuracy controls.',
            'Loss partitions are correlated weighted residual contributions, not physical energy attribution.',
            'No selected-point projection is available when the measured origin differs from the final selection.']}
    require(pin(path)==source,'Parent result mutated during read')
    summary.update(parent_result_unchanged=True,elapsed_seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,local_numerical_threads=1)
    out=Path(output);out.mkdir(parents=True,exist_ok=False);(out/'RESULT.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    condition='undefined' if summary['full_condition'] is None else f"{summary['full_condition']:.3f}"
    lines=['# Donor residual-span diagnostic','',f"Actual stop: `{result['stopping_reason']}` after {result['completed_outer_steps']} steps. Selected development controls: **{selected['score']['development_checks_passed']}**.",'',
        f"The measured Jacobian is at iteration {entry['outer_step']}'s origin. It equals the selected point: **{same_origin}**. Rank {rank}/8, condition {condition}, unchanged rcond 1e-6.",'',
        f"At that origin the actual weighted loss is {original['loss']:.9f}; unrestricted local linear residual loss is {projection['loss']:.9f} ({summary['locally_projected_loss_fraction']:.3%} projected). This is an unevaluated local linear prediction. It exceeds the trust box: {summary['unrestricted_exceeds_trust']}.",'',
        '| Block | Origin weighted loss | Linear projected residual loss |','|---|---:|---:|']
    for block in ('energy','La','Dy','paired'):lines.append(f"| {block} | {original['blocks'][block]:.9f} | {projection['blocks'][block]:.9f} |")
    lines+=['','| Trial | Predicted loss | Actual loss | Actual improvement | Controls passed |','|---|---:|---:|---:|---|']
    for row in trials:lines.append(f"| {row['id']} | {row['predicted_loss']:.9f} | {row['actual_loss']:.9f} | {row['actual_improvement']:.9f} | {row['actual_controls_passed']} |")
    lines+=['','The JSON retains every wet/dry force fixture, mixed wet-minus-dry energy partition and derivative linearity check. These describe the weighted development loss; they do not attribute physical energy. Only actual controls determine a pass. No global family, protein-transfer or physical-accuracy conclusion follows.','',
        f"Saved-only analysis: {summary['elapsed_seconds']:.3f} s wall, {summary['CPU_seconds']:.3f} CPU-s, one thread, zero native calls. Parent result is unchanged."]
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'result':pin(out/'RESULT.json'),'report':pin(out/'REPORT.md')}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();main(a.result,a.output)
