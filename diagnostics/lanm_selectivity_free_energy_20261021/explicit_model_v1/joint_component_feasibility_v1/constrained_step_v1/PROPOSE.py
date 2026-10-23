#!/usr/bin/env python3
"""One saved linear-response constrained proposal; never invokes an evaluator."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time
import numpy as np
import scipy
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
THREADS = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
ORIGIN = 'e63e7189dbecdb8e4112aceb5e4515f69644dbfa574b5fe07644d01ae736a037'
GEOMETRIES = ('water_minus','water_central','water_plus','aquo8_central','aquo9_central')
SLICES = {'old_EF':(0,413),'pair_mu':(413,515),'C_mu':(515,527),'chi':(527,531)}
OPTIONS = {'maxiter':200,'ftol':1e-12,'disp':False}


def require(ok,msg):
    if not ok: raise ValueError(msg)


def pin(path):
    path=Path(path).resolve();return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def read(record):
    require(pin(record['path'])==record,'Changed artifact '+record['path'])
    return json.loads(Path(record['path']).read_text())


def close(a,b,label):
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    require(a.shape==b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)),label+' shape/finite')
    require(np.all(abs(a-b)<=1e-10*np.maximum(1.,np.maximum(abs(a),abs(b)))),label+' arithmetic')


def token(v):
    return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def physical(score):
    """Recover exact original physical gates and verify their 413-row weighting."""
    old=score['old_EF']; energy=[]; names=[]
    for group,rows in old['energy_controls'].items():
        for name,row in rows.items():
            energy.append(row['error_kcal_mol']);names.append(group+'/'+name)
            require(row['passed']==(abs(energy[-1])<=1),'Original energy gate')
    require(len(energy)==10,'Ten energy gates')
    force={}; reconstruct=[]
    errs=old['energy_errors_kcal_mol']
    for group in ('balanced_exchange','individual_radial_work'):
        reconstruct.extend(v*math.sqrt(.5/4)*.5 for v in errs[group].values())
    for kind in ('La','Dy','paired'):
        for geo in GEOMETRIES:
            row=old['per_fixture'][kind][geo]
            raw=np.asarray(row['raw_body_components_kcal_mol_A'],dtype=float)
            n=8 if geo=='aquo8_central' else 9 if geo=='aquo9_central' else 1
            require(raw.shape==(3+6*n,),'Exact body-component denominator')
            w=np.r_[np.full(3,1/6),np.full(6*n,1/(12*n))]
            rms=float(np.sqrt(np.dot(w,raw*raw))); maximum=float(np.max(abs(raw)))
            close([rms,maximum],[row['rigid_RMS_kcal_mol_A'],row['maximum_body_component_error_kcal_mol_A']],'Original physical force gate replay')
            require(row['passed']==(rms<=2 and maximum<=5),'Original force gate boolean')
            force[kind+'/'+geo]={'raw':raw,'weights':w,'passed':row['passed']}
            reconstruct.extend((raw*np.sqrt(w)/(4*math.sqrt(5))).tolist())
    close(reconstruct,old['weighted_residuals'],'Complete original413 physical reconstruction')
    close(np.asarray(reconstruct)/2,score['weighted_residuals'][:413],'Joint531 EF slice')
    paired=errs['paired_radial_work'];work=errs['individual_radial_work']
    for sign in ('minus','plus'):close(paired[sign],work['Dy_'+sign]-work['La_'+sign],'Paired radial gate closure')
    return {'energy':np.asarray(energy),'energy_names':names,'force':force}


def prepare(inputs_path):
    src=json.loads(Path(inputs_path).read_text()); records=src['artifacts']
    data={k:read(v) for k,v in records.items() if Path(v['path']).suffix=='.json'}
    for v in records.values():require(pin(v['path'])==v,'Direct source pin')
    diag,fit,domain=data['diagnostic'],data['fit'],data['domain']
    review,root=data['full_fit_review'],data['root_fit_review']
    require(review['status']=='passed_terminal_joint_fit_review_continuation' and root['status']=='admitted_completed_frozen_fit_provenance_and_numerical_review' and root['independent_review']==records['full_fit_review'],'Admitted original fit')
    require(review['fit_result']==root['fit_result']==records['fit'] and root['selected']==fit['selected'],'Exact reviewed fit/selection joins')
    require(diag['status']=='complete_saved_joint_component_feasibility' and diag['counts']['vectors']==673 and diag['counts']['origins']==24 and diag['counts']['saved_trials']==96,'Complete saved component diagnostic')
    require(diag['sources']['artifacts']['fit_result']==records['fit'] and diag['sources']['artifacts']['domain']==records['domain'],'Exact fit/domain diagnostic join')
    em=read(diag['sources']['artifacts']['evaluator_manifest'])
    require(em['sources']['old_objective']==records['old_objective'] and em['sources']['body_helper']==records['body_helper'],'Exact original physical gate definitions')
    require(fit['status']=='bounded_joint_response_fit_complete' and fit['last_Jacobian_origin_key']==ORIGIN and fit['selected_key']!=ORIGIN,'Original last origin, not selected')
    entry=fit['iterations'][-1];drow=diag['origins'][-1]
    require(entry['outer_step']==24 and entry['origin_key']==drow['key']==ORIGIN and drow['original25_pass'],'Fixed qualified origin24')
    require(entry['derivative_keys']==drow['derivative_keys'] and entry['derivative_spans']==drow['derivative_spans'],'Exact probe/span join')
    keys=[ORIGIN]+entry['derivative_keys']; values={}
    lower=np.asarray([p['lower'] for p in domain['parameters']]);upper=np.asarray([p['upper'] for p in domain['parameters']]);width=upper-lower
    require(lower.shape==(12,) and np.all(width>0),'Twelve box variables')
    for key in keys:
        record=fit['vectors'][key];v=read(record['artifact'])
        require(record['status']==v['status']=='complete' and v['parameters']==record['parameters'] and token(v['parameters'])==key,'Exact saved vector identity')
        require(diag['vectors'][key]['artifact']==record['artifact'],'Diagnostic vector pin')
        require(len(v['cells'])==46 and all(r['status']=='complete' for r in v['cells'].values()),'Complete saved46 cells')
        values[key]={'raw':v,'physical':physical(v['score'])}
    origin=values[ORIGIN];center=np.asarray(origin['raw']['parameters']);z=(center-lower)/width
    require(origin['raw']['score']['original_EF_checks_passed'] is True,'Origin passes all25')
    require(np.all(z>=0) and np.all(z<=1),'Origin domain')
    oldcols=[];energycols=[];forcecols={k:[] for k in origin['physical']['force']}
    for j in range(12):
        a,b=[values[k] for k in entry['derivative_keys'][2*j:2*j+2]]
        for item,sign in ((a,-1),(b,1)):
            target=z.copy();target[j]=np.clip(target[j]+sign*.005,0,1)
            expect=center.copy();expect[j]=np.clip(center[j]+(target[j]-z[j])*width[j],lower[j],upper[j])
            require(item['raw']['parameters']==expect.tolist(),'Exact probe ordering and displacement')
            require(item['physical']['energy_names']==origin['physical']['energy_names'],'Physical gate order')
        span=(b['raw']['parameters'][j]-a['raw']['parameters'][j])/width[j]
        require(span==entry['derivative_spans'][j] and span>0,'Actual normalized span')
        oldcols.append((np.asarray(b['raw']['score']['weighted_residuals'])-np.asarray(a['raw']['score']['weighted_residuals']))/span)
        energycols.append((b['physical']['energy']-a['physical']['energy'])/span)
        for name in forcecols:
            require(np.array_equal(a['physical']['force'][name]['weights'],origin['physical']['force'][name]['weights']) and np.array_equal(b['physical']['force'][name]['weights'],origin['physical']['force'][name]['weights']),'Unchanged body weights')
            forcecols[name].append((b['physical']['force'][name]['raw']-a['physical']['force'][name]['raw'])/span)
    J=np.column_stack(oldcols);close(J,entry['Jacobian_normalized'],'Original saved531 Jacobian');close(J,drow['Jacobian_original531'],'Component diagnostic Jacobian')
    r=np.asarray(origin['raw']['score']['weighted_residuals']); blocks={k:{'r':r[a:b],'J':J[a:b]} for k,(a,b) in SLICES.items()}
    require(diag['component_order']==drow['Jacobian_components_order']==['P8','B8','P9','B9'],'Component order')
    error=np.array([diag['vectors'][ORIGIN]['component_errors_kcal_mol'][k] for k in diag['component_order']]);jc=np.asarray(drow['Jacobian_components'])
    # Reconstruct the new columns from the already admitted per-vector literal-cycle values.
    for j in range(12):
        ka,kb=entry['derivative_keys'][2*j:2*j+2]
        ea=np.array([diag['vectors'][ka]['component_errors_kcal_mol'][k] for k in diag['component_order']]);eb=np.array([diag['vectors'][kb]['component_errors_kcal_mol'][k] for k in diag['component_order']])
        close(jc[:,j],(eb-ea)/entry['derivative_spans'][j],'Component secant replay')
    for name,ix in (('P_energy',[0,2]),('B_energy',[1,3])):blocks[name]={'r':error[ix]/math.sqrt(2),'J':jc[ix]/math.sqrt(2)}
    forces={k:{'r':v['raw'],'weights':v['weights'],'J':np.column_stack(forcecols[k])} for k,v in origin['physical']['force'].items()}
    return {'source':src,'origin':ORIGIN,'parameters':center,'z':z,'width':width,'lower':lower,'upper':upper,
            'bounds_lower':np.r_[np.maximum(-.025,-z),0.],'bounds_upper':np.r_[np.minimum(.025,1-z),1.],
            'blocks':blocks,'energy':{'r':origin['physical']['energy'],'J':np.column_stack(energycols),'names':origin['physical']['energy_names']},
            'forces':forces,'derivative_keys':entry['derivative_keys'],'derivative_spans':entry['derivative_spans'],
            'origin_artifact':fit['vectors'][ORIGIN]['artifact']}


def constraints(problem,x):
    """All constraints return >=0; analytic Jacobian has13 columns."""
    d,t=x[:-1],x[-1];values=[];rows=[];names=[]
    def add(name,value,gradient):names.append(name);values.append(float(value));rows.append(gradient)
    for name,b in problem['blocks'].items():
        v=b['r']+b['J']@d;den=float(b['r']@b['r']);require(den>1e-24,'Degenerate block '+name)
        component=name in ('P_energy','B_energy')
        add('block/'+name,1.-(t if component else 0.)-float(v@v)/den,np.r_[-2*v@b['J']/den,-1. if component else 0.])
    e=problem['energy'];v=e['r']+e['J']@d
    for i,name in enumerate(e['names']):
        add('energy/upper/'+name,1-v[i],np.r_[-e['J'][i],0.]);add('energy/lower/'+name,1+v[i],np.r_[e['J'][i],0.])
    for name,b in problem['forces'].items():
        v=b['r']+b['J']@d;w=b['weights']
        add('force/RMS/'+name,1-float(np.dot(w,v*v))/4,np.r_[-(w*v)@b['J']/2,0.])
        for i in range(len(v)):
            add('force/upper/'+name+'/'+str(i),1-v[i]/5,np.r_[-b['J'][i]/5,0.]);add('force/lower/'+name+'/'+str(i),1+v[i]/5,np.r_[b['J'][i]/5,0.])
    return np.asarray(values),np.asarray(rows),names


def predictions(problem,x):
    d=x[:-1];blocks={};gates={}
    for name,b in problem['blocks'].items():
        v=b['r']+b['J']@d;factor=4 if name in SLICES else 1
        blocks[name]={'origin_norm':float(np.linalg.norm(b['r'])),'predicted_norm':float(np.linalg.norm(v)),
                      'origin_score':factor*float(b['r']@b['r']),'predicted_score':factor*float(v@v)}
    e=problem['energy'];v=e['r']+e['J']@d
    for i,name in enumerate(e['names']):gates['energy/'+name]={'error_kcal_mol':float(v[i]),'passed_exact_threshold':bool(abs(v[i])<=1),'excess_kcal_mol':max(0.,float(abs(v[i])-1))}
    for name,b in problem['forces'].items():
        v=b['r']+b['J']@d;rms=float(np.sqrt(np.dot(b['weights'],v*v)));maximum=float(np.max(abs(v)))
        gates['force/'+name]={'rigid_RMS_kcal_mol_A':rms,'maximum_body_component_error_kcal_mol_A':maximum,
                              'raw_body_components_kcal_mol_A':v.tolist(),'passed_exact_threshold':rms<=2 and maximum<=5,
                              'RMS_excess':max(0.,rms-2),'maximum_excess':max(0.,maximum-5)}
    values,_,names=constraints(problem,x)
    bound=max(0.,float(np.max(problem['bounds_lower']-x)),float(np.max(x-problem['bounds_upper'])))
    return {'blocks':blocks,'original25':gates,'all25_pass_exact_threshold':all(g['passed_exact_threshold'] for g in gates.values()),
            'constraints':dict(zip(names,values.tolist())),'constraint_violation':max(0.,float(-np.min(values))),
            'bound_violation':bound,'numerically_feasible':bool(max(bound,float(-np.min(values)))<=1e-8)}


def run(inputs_path,release_path,out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    result={'status':'started','implementation':pin(__file__),'inputs':pin(inputs_path),'new_native_calls':0,'new_reference_calls':0,'new_parameter_evaluations':0,'solver_calls':0,'error':None}
    try:
        require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
        src=json.loads(Path(inputs_path).read_text());require(pin(src['specification']['path'])==src['specification'],'Frozen specification')
        release=json.loads(Path(release_path).read_text())
        require(release['status']=='released_saved_constrained_component_proposal' and release['implementation']==pin(__file__) and release['inputs']==pin(inputs_path) and release['specification']==src['specification'] and release['solver_call_cap']==1 and release['native_program_cap']==0,'Exact root numerical-only release')
        problem=prepare(inputs_path)
        result.update(specification=src['specification'],root_release=pin(release_path),sources=src,origin=ORIGIN,origin_artifact=problem['origin_artifact'],derivative_keys=problem['derivative_keys'],derivative_spans=problem['derivative_spans'],numerical_versions={'numpy':np.__version__,'scipy':scipy.__version__})
        degenerate=[k for k,b in problem['blocks'].items() if np.linalg.norm(b['r'])<=1e-12]
        result['degenerate_blocks']=degenerate
        require(not degenerate,'Degenerate block norms: '+repr(degenerate))
        x0=np.zeros(13);initial=predictions(problem,x0);require(initial['numerically_feasible'] and initial['all25_pass_exact_threshold'],'Feasible original zero step')
        result['initial']=initial;result['solver_calls']=1
        opt=minimize(lambda x:-float(x[-1]),x0,jac=lambda x:np.r_[np.zeros(12),-1.],method='SLSQP',
                     bounds=list(zip(problem['bounds_lower'],problem['bounds_upper'])),
                     constraints={'type':'ineq','fun':lambda x:constraints(problem,x)[0],'jac':lambda x:constraints(problem,x)[1]},options=OPTIONS)
        require(opt.x.shape==(13,) and np.all(np.isfinite(opt.x)),'Finite solver point')
        full=predictions(problem,opt.x);half=predictions(problem,opt.x*.5)
        result.update(status='numerically_admitted_unevaluated_local_proposal' if opt.success and full['numerically_feasible'] else 'failed_numerical_local_proposal',
                      solver={'method':'SLSQP','options':OPTIONS,'success':bool(opt.success),'status':int(opt.status),'message':str(opt.message),'iterations':int(opt.nit),'function_evaluations':int(opt.nfev),'Jacobian_evaluations':int(opt.njev),'multipliers':np.asarray(opt.get('multipliers',[])).tolist()},
                      normalized_direction=opt.x[:-1].tolist(),physical_direction=(problem['width']*opt.x[:-1]).tolist(),origin_parameters=problem['parameters'].tolist(),
                      unevaluated_parameters=(problem['parameters']+problem['width']*opt.x[:-1]).tolist(),common_fractional_squared_norm_reduction=float(opt.x[-1]),
                      trust_lower=problem['bounds_lower'].tolist(),trust_upper=problem['bounds_upper'].tolist(),half=half,full=full,
                      global_optimum_claim=False,nonlinear_native_feasibility_qualified=False,physical_accuracy_qualified=False,fresh_validation_released=False)
    except Exception as exc:result.update(status='failed_numerical_local_proposal',error=repr(exc))
    result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'native_programs':0,'scheduler_allocated_CPU_seconds':0}
    with (out/'RESULT.json').open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',required=True);p.add_argument('--root-release',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    result=run(a.inputs,a.root_release,a.out);print(json.dumps({k:result.get(k) for k in ('status','error','solver_calls','cost')},indent=2))
    raise SystemExit(0 if result['status']=='numerically_admitted_unevaluated_local_proposal' else 1)
