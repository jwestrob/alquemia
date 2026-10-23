#!/usr/bin/env python3
"""Check the already stored half point by scalar sums. No optimization or native API."""
import hashlib
import json
import math
import os
from pathlib import Path
import time

BASE=Path(__file__).resolve().parents[6]
REL=Path('explicit_model_v1/joint_component_feasibility_v1/constrained_step_v1')
W=BASE/'workspaces/lanm_selectivity_free_energy_20261021'/REL
SLICES={'old_EF':(0,413),'pair_mu':(413,515),'C_mu':(515,527),'chi':(527,531)}

def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def read(p):
    if isinstance(p,dict):assert pin(p['path'])==p;p=p['path']
    return json.loads(Path(p).read_text())

def dot(a,b):return math.fsum(x*y for x,y in zip(a,b))

def norm(a):return math.sqrt(dot(a,a))

def near(a,b):
    assert math.isfinite(a) and math.isfinite(b)
    assert abs(a-b)<=1e-10*max(1,abs(a),abs(b))
    return abs(a-b)

def main():
    wall,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    rp=pin(W/'actual_v1/RESULT.json');r=read(rp)
    assert r['status']=='failed_numerical_local_proposal' and r['solver']['success'] is False and r['solver']['status']==8
    assert r['solver_calls']==1 and r['new_native_calls']==r['new_parameter_evaluations']==0
    sources=r['sources']['artifacts'];fit=read(sources['fit']);domain=read(sources['domain']);diag=read(sources['diagnostic'])
    origin=r['origin'];entry=fit['iterations'][-1]
    assert entry['origin_key']==origin==fit['last_Jacobian_origin_key']!=fit['selected_key'] and entry['outer_step']==24
    assert entry['derivative_keys']==r['derivative_keys'] and entry['derivative_spans']==r['derivative_spans']
    keys=[origin]+entry['derivative_keys'];vectors={k:read(fit['vectors'][k]['artifact']) for k in keys}
    assert len(vectors)==25 and all(v['status']=='complete' for v in vectors.values())
    center=vectors[origin]['parameters'];width=[p['upper']-p['lower'] for p in domain['parameters']]
    assert center==r['origin_parameters'] and vectors[origin]['score']['original_EF_checks_passed'] is True
    d=[x/2 for x in r['normalized_direction']];t=r['common_fractional_squared_norm_reduction']/2
    half=[x+w*y for x,w,y in zip(center,width,d)];z=[(x-p['lower'])/w for x,p,w in zip(center,domain['parameters'],width)]
    boxmargins=[]
    for j,p in enumerate(domain['parameters']):
        near(r['physical_direction'][j],width[j]*r['normalized_direction'][j]);near(r['unevaluated_parameters'][j],center[j]+r['physical_direction'][j])
        boxmargins += [d[j]-max(-.025,-z[j]),min(.025,1-z[j])-d[j],(half[j]-p['lower'])/width[j],(p['upper']-half[j])/width[j]]
    boxmargins += [t,1-t]
    assert min(boxmargins)>0
    # Build scalar secant rows from the original pinned vectors, not producer functions.
    def response(getter):
        base=getter(vectors[origin],origin);jac=[]
        for row in range(len(base)):
            jac.append([(getter(vectors[entry['derivative_keys'][2*j+1]],entry['derivative_keys'][2*j+1])[row]-getter(vectors[entry['derivative_keys'][2*j]],entry['derivative_keys'][2*j])[row])/entry['derivative_spans'][j] for j in range(12)])
        return base,[x+dot(j,d) for x,j in zip(base,jac)]
    residual,predicted=response(lambda v,k:v['score']['weighted_residuals'])
    component_order=['P8','B8','P9','B9']
    component,predicted_component=response(lambda v,k:[diag['vectors'][k]['component_errors_kcal_mol'][n] for n in component_order])
    blocks={n:(residual[a:b],predicted[a:b]) for n,(a,b) in SLICES.items()}
    for n,ix in (('P_energy',(0,2)),('B_energy',(1,3))):blocks[n]=([component[i]/math.sqrt(2) for i in ix],[predicted_component[i]/math.sqrt(2) for i in ix])
    constraints={};block_values={};gates={};maxdiff=0
    for n,(old,new) in blocks.items():
        c=1-(t if n in ('P_energy','B_energy') else 0)-dot(new,new)/dot(old,old);constraints['block/'+n]=c
        stored=r['half']['blocks'][n]
        maxdiff=max(maxdiff,near(norm(old),stored['origin_norm']),near(norm(new),stored['predicted_norm']))
        block_values[n]=dict(origin_norm=norm(old),half_norm=norm(new),squared_norm_ratio=dot(new,new)/dot(old,old),fractional_squared_norm_reduction=1-dot(new,new)/dot(old,old))
    energy_names=[(g,n) for g,rows in vectors[origin]['score']['old_EF']['energy_controls'].items() for n in rows]
    for group,name in energy_names:
        _,v=response(lambda item,k:[item['score']['old_EF']['energy_controls'][group][name]['error_kcal_mol']]);v=v[0];label=group+'/'+name
        constraints['energy/upper/'+label]=1-v;constraints['energy/lower/'+label]=1+v
        stored=r['half']['original25']['energy/'+label];maxdiff=max(maxdiff,near(v,stored['error_kcal_mol']))
        assert abs(v)<=1 and stored['passed_exact_threshold'] is True
        gates['energy/'+label]=dict(error_kcal_mol=v,absolute_margin=1-abs(v))
    for kind in ('La','Dy','paired'):
        for geo in ('water_minus','water_central','water_plus','aquo8_central','aquo9_central'):
            _,new=response(lambda item,k:item['score']['old_EF']['per_fixture'][kind][geo]['raw_body_components_kcal_mol_A'])
            waters=8 if geo=='aquo8_central' else 9 if geo=='aquo9_central' else 1
            assert len(new)==3+6*waters
            weights=[1/6]*3+[1/(12*waters)]*(6*waters)
            rms2=math.fsum(w*v*v for w,v in zip(weights,new));rms=math.sqrt(rms2);mx=max(abs(v) for v in new);label=kind+'/'+geo
            constraints['force/RMS/'+label]=1-rms2/4
            for j,v in enumerate(new):
                constraints['force/upper/'+label+'/'+str(j)]=1-v/5;constraints['force/lower/'+label+'/'+str(j)]=1+v/5
            stored=r['half']['original25']['force/'+label]
            maxdiff=max(maxdiff,near(rms,stored['rigid_RMS_kcal_mol_A']),near(mx,stored['maximum_body_component_error_kcal_mol_A']))
            for a,b in zip(new,stored['raw_body_components_kcal_mol_A']):maxdiff=max(maxdiff,near(a,b))
            assert rms<=2 and mx<=5 and stored['passed_exact_threshold'] is True
            gates['force/'+label]=dict(rigid_RMS_kcal_mol_A=rms,maximum_body_component_error_kcal_mol_A=mx,RMS_margin=2-rms,maximum_margin=5-mx)
    assert len(constraints)==851 and len(gates)==25 and set(constraints)==set(r['half']['constraints'])
    for n,v in constraints.items():maxdiff=max(maxdiff,near(v,r['half']['constraints'][n]))
    limiting=min(constraints,key=constraints.get)
    assert constraints[limiting]>0
    assert all(v['fractional_squared_norm_reduction']>0 for v in block_values.values())
    out=W/'independent_half_review_v1';out.mkdir(parents=True,exist_ok=False)
    review=dict(status='passed_saved_half_linear_feasibility_review',implementation=pin(__file__),original_result=rp,original_status=r['status'],original_protocol_passed=False,original_solver_status=8,
        original_implementation=r['implementation'],source_inputs=r['inputs'],specification=r['specification'],origin=origin,origin_artifact=r['origin_artifact'],raw_vector_pins={k:fit['vectors'][k]['artifact'] for k in keys},
        selection_rule='The already reported half point: exactly0.5 times the sole returned normalized direction and t; no search or optimization.',
        normalized_half_direction=d,physical_half_parameters=half,common_fractional_squared_norm_reduction_t=t,half_parameter_token=hashlib.sha256(json.dumps(half,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest(),
        minimum_normalized_box_margin=min(boxmargins),constraint_count=851,constraints=constraints,minimum_constraint=dict(name=limiting,value=constraints[limiting]),maximum_scalar_replay_difference=maxdiff,
        original25=gates,all25_exact_threshold_pass=True,blocks=block_values,strict_linear_feasibility_verified=True,native_candidate_evaluated=False,native_candidate_qualified=False,global_optimum_claim=False,
        interpretation=['Original failed solver attempt remains failed under unchanged protocol.','The fixed stored half point is an independently verified interior witness for the saved secant constraints, regardless of solver termination.','A separately declared single-candidate consumed-data native pilot is defensible; it requires its own exact bindings and full original gates, and is not an automatic continuation or original protocol pass.','Predicted pair-dipole improvement is small; nonlinear numerical and physical outcomes remain unknown. No new term, fitted target, threshold, fresh validation or global capability claim follows.'],
        cost=dict(wall_seconds=time.monotonic()-wall,CPU_seconds=time.process_time()-cpu,threads=1,new_native_calls=0,new_parameter_evaluations=0,optimizer_calls=0))
    p=out/'REVIEW.json';p.write_text(json.dumps(review,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(review=pin(p),minimum_constraint=review['minimum_constraint'],maximum_scalar_replay_difference=maxdiff,minimum_box_margin=min(boxmargins),cost=review['cost']),indent=2))

if __name__=='__main__':main()
