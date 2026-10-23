#!/usr/bin/env python3
"""Bounded saved arithmetic review; no solver, producer import, or native calls."""
import collections
import hashlib
import json
import math
import os
from pathlib import Path
import time
import numpy as np

BASE = Path(__file__).resolve().parents[5]
REL = Path('explicit_model_v1/joint_component_feasibility_v1')
W = BASE/'workspaces/lanm_selectivity_free_energy_20261021'/REL
SLICES = dict(old_EF=(0,413), pair_mu=(413,515), C_mu=(515,527), chi=(527,531))
ORDER = ('P8','B8','P9','B9')

def pin(p):
    p=Path(p).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def read(p):
    if isinstance(p,dict):
        assert pin(p['path'])==p
        p=p['path']
    return json.loads(Path(p).read_text())

def equal(a,b):
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    assert a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    err=np.abs(a-b)
    assert np.all(err <= 1e-10*np.maximum(1,np.maximum(abs(a),abs(b))))
    return float(np.max(err,initial=0))

def main():
    tick,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','OMP_THREAD_LIMIT'))
    resultpin=pin(W/'actual_v1/RESULT.json'); result=read(resultpin)
    assert result['status']=='complete_saved_joint_component_feasibility' and result['error'] is None
    src=result['sources']['artifacts'];fit=read(src['fit_result']);manifest=read(src['evaluator_manifest']);domain=read(src['domain'])
    assert len(result['vectors'])==673 and len(result['origins'])==24 and len(result['trials'])==96
    assert result['counts']['LPs']==sum(len(o['LP']['modes']) for o in result['origins'])==72
    assert result['degenerate_origin_steps']==[]
    lower=np.array([r['lower'] for r in domain['parameters']]);width=np.array([r['upper']-r['lower'] for r in domain['parameters']])
    origin=result['origins'][-1];entry=fit['iterations'][-1];key=origin['key']
    assert key==entry['origin_key']==fit['last_Jacobian_origin_key']!=fit['selected_key']
    keys=[key]+entry['derivative_keys']+entry['trial_keys'];saved={k:read(fit['vectors'][k]['artifact']) for k in set(keys)}
    assert len(keys)==len(saved)==29
    errors={};residuals={};params={};cycle_error=0
    for k,v in saved.items():
        assert v['status']=='complete' and len(v['cells'])==46 and set(v['cells'])==set(manifest['rows'])
        assert all(r['status']=='complete' for r in v['cells'].values())
        energy={n:r['energy_kcal_mol'] for n,r in manifest['static'].items()}
        energy.update({n:r['energy'] for n,r in v['cells'].items()})
        parents={}
        for n in (8,9):
            for m in ('La','Dy'):
                name=f'aquo{n}_central_{m}';ion=energy['ion_'+m]
                pair=math.fsum(energy[name+f'__dimer_{i:02d}']-ion-energy[name+f'__water_{i:02d}'] for i in range(1,n+1))
                collective=energy[name]-ion-energy[name+'__shell']
                parents[name]=dict(C=collective,P=pair,B=collective-pair)
                for channel,value in parents[name].items():
                    cycle_error=max(cycle_error,equal(value,result['vectors'][k]['parents_kcal_mol'][name][channel]))
        e=[]
        for n in (8,9):
            for channel in ('P','B'):
                e.append(parents[f'aquo{n}_central_Dy'][channel]-parents[f'aquo{n}_central_La'][channel]-result['reference_targets'][str(n)][channel])
        errors[k]=np.asarray(e);residuals[k]=np.asarray(v['score']['weighted_residuals']);params[k]=np.asarray(v['parameters'])
        equal(errors[k],[result['vectors'][k]['component_errors_kcal_mol'][n] for n in ORDER])
        assert residuals[k].shape==(531,)
        assert v['score']['original_EF_checks_passed']==result['vectors'][k]['original25_pass']
    J=[];K=[];z=(params[key]-lower)/width
    for j in range(12):
        minus,plus=entry['derivative_keys'][2*j:2*j+2]
        for k,s in ((minus,-1),(plus,1)):
            expected=params[key].copy();target=np.clip(z[j]+s*.005,0,1)
            expected[j]=np.clip(params[key][j]+(target-z[j])*width[j],lower[j],lower[j]+width[j])
            assert expected.tolist()==params[k].tolist()
        span=(params[plus][j]-params[minus][j])/width[j]
        assert span==entry['derivative_spans'][j]>0
        J.append((residuals[plus]-residuals[minus])/span);K.append((errors[plus]-errors[minus])/span)
    J=np.asarray(J).T;K=np.asarray(K).T
    jac_error=equal(J,origin['Jacobian_original531']);equal(J,entry['Jacobian_normalized'])
    component_jac_error=equal(K,origin['Jacobian_components'])
    blocks={n:(residuals[key][a:b],J[a:b]) for n,(a,b) in SLICES.items()}
    blocks.update(P_energy=(errors[key][[0,2]]/math.sqrt(2),K[[0,2]]/math.sqrt(2)),B_energy=(errors[key][[1,3]]/math.sqrt(2),K[[1,3]]/math.sqrt(2)))
    gradients={n:(r@j)/(r@r) for n,(r,j) in blocks.items()}
    for n,g in gradients.items():equal(g,origin['LP']['fractional_norm_gradients'][n])
    modes={}
    for name,mode in origin['LP']['modes'].items():
        dx=np.array(mode['normalized_direction']);t=mode['t'];x=np.r_[dx,t]
        names=['P_energy','B_energy']+mode['protected_blocks']
        A=np.array([np.r_[gradients[n],float(n in ('P_energy','B_energy'))] for n in names])
        lo=np.maximum(-.025,-z);hi=np.minimum(.025,1-z);cost=np.r_[np.zeros(12),-1]
        y=np.array(mode['inequality_marginals']);lm=np.array(mode['lower_marginals']);um=np.array(mode['upper_marginals'])
        checks=dict(primal=max(0,float(np.max(A@x)),float(np.max(lo-dx)),float(np.max(dx-hi)),-t),
            dual_sign=max(0,float(np.max(y)),float(np.max(-lm)),float(np.max(um)),abs(float(um[-1]))),
            stationarity=float(np.max(abs(cost-A.T@y-lm-um))),
            gap=abs(float(cost@x-lo@lm[:12]-hi@um[:12])),
            complementarity=max(float(np.max(abs(y*(A@x)))),float(np.max(abs(lm[:12]*(dx-lo)))),float(np.max(abs(um[:12]*(hi-dx)))),abs(float(lm[-1]*t))))
        assert max(checks.values())<=1e-7
        predictions={}
        for n,(r,j) in blocks.items():
            v=j@dx;pred=dict(origin=float(np.linalg.norm(r)),half=float(np.linalg.norm(r+.5*v)),full=float(np.linalg.norm(r+v)),linear_norm_squared_change=float(2*r@v),quadratic_norm_squared_change=float(v@v))
            old=mode['block_predictions'][n];equal([pred['origin'],pred['half'],pred['full']],[old['origin_norm'],old['half_step_linearized_norm'],old['full_step_linearized_norm']])
            predictions[n]=pred
        modes[name]=dict(t=t,certificate=checks,predictions=predictions)
    counts=collections.Counter();by_step=collections.Counter();qualified_step=collections.Counter();final_trials=[]
    for i,(iteration,actual_origin) in enumerate(zip(fit['iterations'],result['origins']),1):
        assert iteration['outer_step']==actual_origin['outer_step']==i
        for j,(tk,point) in enumerate(zip(iteration['trial_keys'],iteration['trial_points'])):
            row=result['trials'][4*(i-1)+j];okey=iteration['origin_key']
            assert (row['outer_step'],row['trial_index'],row['trial_key'],row['origin_key'])==(i,j,tk,okey)
            assert row['trial_original25_pass']==fit['vectors'][tk]['original_EF_checks_passed']
            assert row['origin_original25_pass']==fit['vectors'][okey]['original_EF_checks_passed']
            gate_values=[]
            for label,g in row['actual_original25'].items():
                passes=abs(g['error_kcal_mol'])<=1 if label.startswith('energy/') else g['rigid_RMS_kcal_mol_A']<=2 and g['maximum_body_component_error_kcal_mol_A']<=5
                assert passes==g['passed'];gate_values.append(passes)
            assert len(gate_values)==25 and all(gate_values)==row['trial_original25_pass']
            a=row['actual_scores'];o=result['vectors'][okey]['scores'];eps={n:1e-10*max(1,abs(a[n]),abs(o[n])) for n in a}
            both=all(a[n]<o[n]-eps[n] for n in ('P_energy','B_energy'));protected=all(a[n]<=o[n]+eps[n] for n in SLICES)
            assert both==row['flags']['both_component_scores_decrease'] and protected==row['flags']['all_original_four_nonincrease']
            counts['trials']+=1;counts['from_qualified_origin']+=int(row['origin_original25_pass'])
            if both and protected and all(gate_values):
                counts['compatible_trials']+=1;by_step[i]+=1
                if row['origin_original25_pass']:counts['compatible_from_qualified_origin']+=1;qualified_step[i]+=1
            if i==24:
                dz=(params[tk]-params[key])/width;pr=residuals[key]+J@dz;pe=errors[key]+K@dz
                equal(pr,row['predicted_original_residuals531']);equal(pe,row['predicted_component_errors_kcal_mol']);equal(residuals[tk],row['actual_original_residuals531']);equal(errors[tk],row['actual_component_errors_kcal_mol'])
                final_trials.append(dict(trial_index=j,method=point['method'],fraction=point['fraction'],key=tk,original25_pass=all(gate_values),component_max_linearization_error=float(np.max(abs(errors[tk]-pe))),EF_linearization_L2=float(np.linalg.norm(residuals[tk][:413]-pr[:413])),component_errors=errors[tk].tolist(),original_block_scores=a))
    out=W/'math_review_v1/actual_v1';out.mkdir(parents=True,exist_ok=False)
    review=dict(status='passed_bounded_saved_component_math_review',implementation=pin(__file__),result=resultpin,
        sources={n:src[n] for n in ('fit_result','domain','evaluator_manifest')},raw_vectors_replayed={k:fit['vectors'][k]['artifact'] for k in keys},
        counts=dict(counts),compatible_trials_by_step=dict(by_step),compatible_from_qualified_origin_by_step=dict(qualified_step),
        final_origin=key,final_origin_is_selected=False,cycle_max_abs_difference=cycle_error,original_Jacobian_max_abs_difference=jac_error,component_Jacobian_max_abs_difference=component_jac_error,
        final_origin_LP=modes,final_four_trials=final_trials,
        scope=['All96 trial identities/derived flags/25 stored metric gates; native force generation not repeated','Final29 raw saved vectors, literal cycle arithmetic,24 probe ordering/spans,531 and4 component secants,three existing LP certificates, four trial predictions','No producer/LP import or solve, new point, native call, fitting or selection'],
        limits=['Positive tangent t does not establish a finite residual-norm-preserving step','Last origin differs from terminal selected vector','Compatible historical trials demonstrate small local improvement only; large component errors remain','No global no-go, physical validation or new functional-term identification'],
        cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1,new_native_calls=0,new_parameter_evaluations=0,new_optimizations=0))
    p=out/'REVIEW.json';p.write_text(json.dumps(review,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(review=pin(p),counts=dict(counts),cost=review['cost']),indent=2))

if __name__=='__main__':main()
