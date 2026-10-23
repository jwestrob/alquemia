#!/usr/bin/env python3
"""Saved-only alpha-response analysis; projections never become candidates."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time
import numpy as np

GEOMETRIES=('water_minus','water_central','water_plus','aquo8_central','aquo9_central')
N_WATERS=(1,1,1,8,9)
PROBES=tuple(m+'_'+s for m in ('La','Dy') for s in ('minus5','minus1','plus1','plus5'))
BASE_ALPHA={'La':1.6597,'Dy':.7280}

def require(value,message):
    if not value:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(record):
    require(pin(record['path'])==record,'Changed input: '+record['path'])
    return json.loads(Path(record['path']).read_text())

def linear_checks(residual):
    """Decode the declared 413-vector into unchanged predicted control units."""
    r=np.asarray(residual);require(r.shape==(413,) and np.isfinite(r).all(),'Invalid residual')
    energy=r[:8]/math.sqrt(.25*.5/4)
    errors={name:float(value) for name,value in zip(
        ('exchange_water_minus','exchange_water_central','exchange_water_plus','exchange_W9',
         'work_La_minus','work_La_plus','work_Dy_minus','work_Dy_plus'),energy)}
    errors.update(paired_work_minus=float(energy[6]-energy[4]),paired_work_plus=float(energy[7]-energy[5]))
    forces={};offset=8
    for kind in ('La','Dy','paired'):
        forces[kind]={}
        for g,nw in zip(GEOMETRIES,N_WATERS):
            n=3+6*nw;part=r[offset:offset+n];offset+=n
            weights=np.array([math.sqrt(.25/5)/(2*math.sqrt(2*3))]*3+
                             [math.sqrt(.25/5)/(2*math.sqrt(2*nw*6))]*(6*nw))
            raw=part/weights;rms=float(2*np.linalg.norm(part)/math.sqrt(.25/5));maximum=float(np.max(np.abs(raw)))
            forces[kind][g]={'RMS':rms,'maximum':maximum,'passed':rms<=2 and maximum<=5}
    require(offset==413,'Residual partition differs')
    count=sum(abs(x)<=1 for x in errors.values())+sum(x['passed'] for group in forces.values() for x in group.values())
    return {'energy_errors':errors,'force_controls':forces,'passed_checks':int(count),'declared_checks':25,'all_checks_passed':count==25}

def projection(matrix,residual):
    step,_,rank,singular=np.linalg.lstsq(matrix,-residual,rcond=1e-6)
    remaining=residual+matrix@step
    return {'rank':int(rank),'singular_values':singular.tolist(),'relative_cutoff':1e-6,
            'absolute_cutoff':float(singular[0]*1e-6),'step':step.tolist(),
            'remaining_weighted_loss':float(remaining@remaining),'predicted_checks':linear_checks(remaining)}

def analyze(data,result):
    records=result['records'];require(set(records)==set(PROBES),'All eight probes required')
    require(result.get('error') is None and result.get('global_stop') is False,'Native run incomplete')
    base=data['baseline_score'];r=np.asarray(base['weighted_residuals']);J=np.asarray(data['six_column_Jacobian'])
    require(J.shape==(413,6) and np.isfinite(J).all(),'Invalid six-column response')
    require(np.isclose(r@r,base['objective'],rtol=1e-12,atol=1e-12),'Baseline loss mismatch')
    baseline_checks=linear_checks(r)
    require(baseline_checks['all_checks_passed']==base['development_checks_passed'],'Baseline gate mismatch')
    for kind,group in baseline_checks['force_controls'].items():
        for g,value in group.items():
            actual=base['per_fixture'][kind][g]
            require(np.allclose([value['RMS'],value['maximum']],[actual['rigid_RMS_kcal_mol_A'],actual['maximum_body_component_error_kcal_mol_A']],rtol=1e-10,atol=1e-10),'Baseline force units mismatch')
    for p in PROBES:
        require(records[p]['status']=='complete','Missing probe: '+p)
        v=np.asarray(records[p]['score']['weighted_residuals']);require(v.shape==(413,) and np.isfinite(v).all(),'Invalid probe residual')
        require(np.isclose(v@v,records[p]['score']['objective'],rtol=1e-12,atol=1e-12),'Probe loss mismatch')
        require(linear_checks(v)['all_checks_passed']==records[p]['score']['development_checks_passed'],'Probe gate mismatch')
    lo=np.array([p['lower'] for p in data['family']['parameters']]);hi=np.array([p['upper'] for p in data['family']['parameters']])
    z=np.array(data['selected_normalized_parameters']);columns={};scales={}
    require(lo.shape==hi.shape==z.shape==(6,) and np.isfinite(lo).all() and np.isfinite(hi).all() and np.isfinite(z).all(),'Six finite parameter bounds and coordinates required')
    require(np.all(hi>lo) and np.all(z>=0) and np.all(z<=1),'Invalid parameter bounds or selected coordinate')
    for suffix,h in (('1',.01),('5',.05)):
        cols=[]
        for metal in ('La','Dy'):
            minus,plus=(records[metal+'_'+sign+suffix] for sign in ('minus','plus'))
            require(minus['fraction']==-h and plus['fraction']==h and minus['metal']==plus['metal']==metal,'Probe coordinate differs')
            for record,sign in ((minus,-1),(plus,1)):
                actual=(record['alpha']-BASE_ALPHA[metal])/BASE_ALPHA[metal]
                require(math.isfinite(actual) and math.isclose(actual,record['eta'],rel_tol=1e-12,abs_tol=1e-14),'Alpha coordinate inconsistent with actual token')
                require(math.isclose(actual,sign*h,rel_tol=0,abs_tol=1e-10),'Alpha token differs from declared probe')
            spacing=plus['eta']-minus['eta'];require(spacing>0,'Nonpositive probe spacing')
            cols.append((np.asarray(plus['score']['weighted_residuals'])-np.asarray(minus['score']['weighted_residuals']))/spacing)
        columns[suffix]=np.column_stack(cols);A=np.column_stack((J,columns[suffix]));fit=projection(A,r);step=np.asarray(fit['step'])
        fit.update(alpha_relative_changes={'La':float(step[6]),'Dy':float(step[7])},
                   implied_alpha={'La':float(1.6597*(1+step[6])),'Dy':float(.7280*(1+step[7]))},
                   implied_six_parameters=(lo+(z+step[:6])*(hi-lo)).tolist(),
                   original_six_bounds_satisfied=bool(np.all(z+step[:6]>=0) and np.all(z+step[:6]<=1)),
                   original_six_trust_satisfied=bool(np.max(np.abs(step[:6]))<=.025),
                   alpha_positive=bool(np.all(1+step[6:]>0)),alpha_step_within_measured_range=bool(np.max(np.abs(step[6:]))<=h),
                   alpha_columns=columns[suffix].tolist(),prediction_only=True,candidate_selected=False)
        scales[suffix+'percent']=fit
    differences={}
    for j,metal in enumerate(('La','Dy')):
        a,b=columns['1'][:,j],columns['5'][:,j];norm=float(np.linalg.norm(a));den=norm*float(np.linalg.norm(b))
        differences[metal]={'one_percent_norm':norm,'five_percent_norm':float(np.linalg.norm(b)),
                            'relative_column_change':float(np.linalg.norm(b-a)/norm) if norm>0 else None,
                            'cosine':float(a@b/den) if den>0 else None}
    return {'status':'saved_alpha_response_analysis_complete','baseline_controls':baseline_checks,
            'original_six_projection':projection(J,r),'augmented_projections':scales,'scale_dependence':differences,
            'actual_probes':{p:{'objective':records[p]['score']['objective'],'all_controls_passed':records[p]['score']['development_checks_passed'],
                                'controls':linear_checks(records[p]['score']['weighted_residuals'])} for p in PROBES},
            'coordinate_convention':'original six normalized coordinates plus relative metal alpha changes',
            'alpha_changes_susceptibility_and_screening':True,'new_native_calls':0,'candidate_selected':False,
            'global_family_impossibility_claimed':False,'physical_accuracy_qualified':False}

def main(result_path,out):
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'One numerical thread required')
    started=time.monotonic();cpu=time.process_time();parent=pin(result_path);result=read(parent);manifest=read(result['manifest']);data=read(manifest['data'])
    require(result.get('status')=='metal_alpha_diagnostic_complete','Parent native diagnostic incomplete')
    require(result.get('native_program_invocations')==80,'Native program count differs from declared 80')
    summary=analyze(data,result);require(pin(result_path)==parent,'Parent changed during analysis')
    summary.update(implementation=pin(__file__),parent_result=parent,manifest=result['manifest'],data=manifest['data'],wall_seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu)
    out=Path(out);out.mkdir(parents=True,exist_ok=False);(out/'RESULT.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    lines=['# Metal-polarizability response diagnostic','','All results use consumed configurations. Projected changes are unevaluated; no candidate is selected.','','| Response matrix | Rank | Predicted residual loss | Predicted W9 energy error |','|---|---:|---:|---:|']
    for name,p in [('Original six',summary['original_six_projection'])]+list(summary['augmented_projections'].items()):
        lines.append(f"| {name} | {p['rank']} | {p['remaining_weighted_loss']:.9f} | {p['predicted_checks']['energy_errors']['exchange_W9']:.6f} |")
    lines+=['','The JSON retains every actual probe, the unchanged25 controls, scale dependence, implied changes, parameter/trust/probe-range flags and separate singular cutoffs. Alpha changes both susceptibility and damping length. These results do not isolate a unique physical cause or qualify a new model.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'result':pin(out/'RESULT.json'),'report':pin(out/'REPORT.md')}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',required=True);p.add_argument('--out',required=True);a=p.parse_args();main(a.result,a.out)
