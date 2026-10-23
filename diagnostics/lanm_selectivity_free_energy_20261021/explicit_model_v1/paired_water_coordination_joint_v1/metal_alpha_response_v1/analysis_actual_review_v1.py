#!/usr/bin/env python3
"""Independent compact actual-column and QR projection review; no raw replay."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_coordination_joint_v1/metal_alpha_response_v1'
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):
    assert pin(p['path'])==p,p['path'];return json.loads(Path(p['path']).read_text())

def main():
    wall,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    ap=pin(W/'analysis_v1/RESULT.json');assert ap['sha256']=='0c0df0b2880d264b69aeab9f9f3162688045c489a31c05667765b43956e00a61'
    a=read(ap);result=read(a['parent_result']);manifest=read(a['manifest']);data=read(a['data']);release_path=W/'analysis_submission_v1/ROOT_RELEASE.json';release=json.loads(release_path.read_text());execution_path=W/'analysis_submission_v1/EXECUTION.json';execution=json.loads(execution_path.read_text())
    control_source=pin(HERE/'analysis_saved_controls_v2.py');assert control_source['sha256']=='386a43e95422495d08d8755ebc103901f0ae3f4e4b0ce5f5b620fd2b24f089d4'
    spec=importlib.util.spec_from_file_location('independent_alpha_units',control_source['path']);unit=importlib.util.module_from_spec(spec);spec.loader.exec_module(unit)
    checks=[]
    def check(name,test):
        assert test,name;checks.append(name)
    check('Frozen analysis/source/parent/release joins',release['implementation']==a['implementation'] and release['parent']==a['parent_result'] and pin(a['implementation']['path'])==a['implementation'] and manifest['data']==a['data'] and result['manifest']==a['manifest'])
    check('Saved analysis executed successfully',execution['returncode']==0 and json.loads(execution['stdout'])['result']==ap)
    check('Native summary complete80 without claiming raw review',result['status']=='metal_alpha_diagnostic_complete' and result['native_program_invocations']==80 and result['error'] is None and not result['global_stop'] and release['independent_raw_review_required_before_final_admission'])
    r=np.array(data['baseline_score']['weighted_residuals']);J=np.array(data['six_column_Jacobian']);lo=np.array([p['lower'] for p in data['family']['parameters']]);hi=np.array([p['upper'] for p in data['family']['parameters']]);z=np.array(data['selected_normalized_parameters'])
    def all_checks(v,reported,label):
        e,f,n=unit.decode(v)
        check(label+' full25 energy/force denominator',len(e)==10 and sum(map(len,f.values()))==15 and reported['declared_checks']==25 and n==reported['passed_checks'] and reported['all_checks_passed']==(n==25))
        check(label+' ten energy units/signs',all(abs(e[k]-reported['energy_errors'][k])<1e-10 for k in e))
        check(label+' fifteen RMS/max gates',all(abs(g['rigid_RMS_kcal_mol_A']-reported['force_controls'][kind][q]['RMS'])<1e-10 and abs(g['maximum_body_component_error_kcal_mol_A']-reported['force_controls'][kind][q]['maximum'])<1e-10 and g['passed']==reported['force_controls'][kind][q]['passed'] for kind,rows in f.items() for q,g in rows.items()))
        return {'energy_errors':e,'force_controls':f,'passed_checks':n,'declared_checks':25,'all_checks_passed':n==25}
    baseline_checks=all_checks(r,a['baseline_controls'],'Actual baseline')
    expected_ids={p['id'] for p in manifest['probes']};check('Exact eight actual probes',set(result['records'])==set(a['actual_probes'])==expected_ids)
    for p,record in result['records'].items():
        check('Complete and consistent actual probe '+p,record['status']=='complete' and abs(np.array(record['score']['weighted_residuals'])@np.array(record['score']['weighted_residuals'])-record['score']['objective'])<1e-10)
        all_checks(record['score']['weighted_residuals'],a['actual_probes'][p]['controls'],'Actual probe '+p)
    projections={};cols={};alpha0={'La':1.6597,'Dy':.7280}
    for name in ('six','1percent','5percent'):
        if name=='six':A=J;reported=a['original_six_projection']
        else:
            suffix=name[0];two=[]
            for metal in ('La','Dy'):
                minus=result['records'][metal+'_minus'+suffix];plus=result['records'][metal+'_plus'+suffix]
                em=(minus['alpha']-alpha0[metal])/alpha0[metal];ep=(plus['alpha']-alpha0[metal])/alpha0[metal]
                check(name+' actual-token spacing '+metal,em==minus['eta'] and ep==plus['eta'] and ep>em)
                two.append((np.array(plus['score']['weighted_residuals'])-np.array(minus['score']['weighted_residuals']))/(ep-em))
            columns=np.column_stack(two);cols[name]=columns;A=np.column_stack((J,columns));reported=a['augmented_projections'][name]
            check(name+' exact alpha-column reconstruction',np.array_equal(columns,np.array(reported['alpha_columns'])))
        singular=np.linalg.svd(A,compute_uv=False);cutoff=singular[0]*1e-6;rank=int(np.count_nonzero(singular>cutoff))
        check(name+' singular spectrum/cutoff/full rank',np.allclose(singular,reported['singular_values'],rtol=1e-12,atol=1e-12) and abs(cutoff-reported['absolute_cutoff'])<1e-14 and reported['relative_cutoff']==1e-6 and rank==reported['rank']==A.shape[1])
        # Full rank was independently checked. Reduced QR gives an independent
        # projection without another lstsq call or a fit/new candidate.
        Q,R=np.linalg.qr(A,mode='reduced');step=np.linalg.solve(R,-Q.T@r);remaining=r+Q@(-Q.T@r)
        check(name+' independent QR step/loss',np.allclose(step,reported['step'],rtol=1e-9,atol=1e-10) and abs(remaining@remaining-reported['remaining_weighted_loss'])<1e-11)
        gates=all_checks(remaining,reported['predicted_checks'],name+' projected')
        item={'rank':rank,'singular_values':singular.tolist(),'absolute_cutoff':float(cutoff),'smallest_retained_to_cutoff_ratio':float(singular[-1]/cutoff),'condition_number_in_declared_coordinates':float(singular[0]/singular[-1]),'column_norms':np.linalg.norm(A,axis=0).tolist(),'QR_step':step.tolist(),'QR_remaining_weighted_loss':float(remaining@remaining),'all25_predicted_checks':gates}
        if name!='six':
            alpha={m:alpha0[m]*(1+step[6+j]) for j,m in enumerate(('La','Dy'))};six=lo+(z+step[:6])*(hi-lo)
            flags={'original_six_bounds_satisfied':bool(np.all((z+step[:6]>=0)&(z+step[:6]<=1))),'original_six_trust_satisfied':bool(np.max(np.abs(step[:6]))<=.025),'alpha_positive':bool(np.all(step[6:]>-1)),'alpha_step_within_measured_range':bool(np.max(np.abs(step[6:]))<=(.01 if name=='1percent' else .05))}
            check(name+' implied physical parameters and extrapolation flags',all(abs(alpha[m]-reported['implied_alpha'][m])<1e-9 for m in alpha) and np.allclose(six,reported['implied_six_parameters'],rtol=1e-10,atol=1e-10) and all(flags[k]==reported[k] for k in flags))
            item.update(implied_alpha_A3=alpha,implied_six_parameters=six.tolist(),flags=flags)
        projections[name]=item
    scale={}
    for j,m in enumerate(('La','Dy')):
        c1,c5=cols['1percent'][:,j],cols['5percent'][:,j];relative=float(np.linalg.norm(c5-c1)/np.linalg.norm(c1));cosine=float(c1@c5/(np.linalg.norm(c1)*np.linalg.norm(c5)))
        check('Actual scale dependence '+m,abs(relative-a['scale_dependence'][m]['relative_column_change'])<1e-12 and abs(cosine-a['scale_dependence'][m]['cosine'])<1e-12)
        scale[m]={'relative_column_change':relative,'cosine':cosine}
    check('No evaluated alpha probe improves baseline loss or passes all25',all(p['objective']>=data['baseline_score']['objective'] and not p['all_controls_passed'] for p in a['actual_probes'].values()))
    check('No candidate/physical claim in analysis',not a['candidate_selected'] and not a['physical_accuracy_qualified'] and not a['global_family_impossibility_claimed'])
    check('All input pins unchanged',pin(ap['path'])==ap and all(pin(a[k]['path'])==a[k] for k in ('parent_result','manifest','data','implementation')))
    output={'status':'independent_actual_projection_review_passed','implementation':pin(__file__),'analysis':ap,'native_result':a['parent_result'],'manifest':a['manifest'],'data':a['data'],'analysis_root_release':pin(release_path),'analysis_execution':pin(execution_path),'unit_decoder_source':control_source,'checks_passed':len(checks),'checks':checks,'baseline_checks':baseline_checks,'independent_projections':projections,'scale_dependence':scale,
        'interpretation':['Two additional linearly independent response directions are retained at both measured scales in the declared six-normalized-plus-two-relative-alpha coordinate convention.','Smallest augmented singular value is over600 times the frozen cutoff; rank8 is not a threshold-adjacent rank decision. Rank/condition under truncated SVD remain coordinate-scaling dependent, not statistical identifiability.','The large unconstrained loss reduction requires negative La and Dy polarizabilities, negative Dy damping and large changes outside original bounds/trust/probe ranges. It is an unevaluated extrapolation, not a usable or physically qualified candidate.','Both augmented projections still fail two of the25 predicted gates: Dy/W9 and paired/W9 maximum body components.','No evaluated one-metal alpha probe improves the baseline loss or passes all25 gates. This does not exclude jointly adjusted constrained directions.','Close1%/5% derivative columns establish local scale consistency only; they do not justify extrapolating to changes of order100–700 percent.','Alpha changes susceptibility and native damping length simultaneously; no unique omitted physical term follows.'],
        'raw_native_admission_performed':False,'independent_raw_review_required_before_final_admission':True,'candidate_selected':False,'new_native_calls':0,'new_LP_or_synthetic_runs':0,'threads':1,'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu}
    out=W/'analysis_actual_review_v1';out.mkdir(parents=True,exist_ok=False);dest=out/'REVIEW.json';dest.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'review':pin(dest),'checks_passed':len(checks),'projection_ranks':{k:v['rank'] for k,v in projections.items()},'wall_seconds':output['wall_seconds'],'CPU_seconds':output['CPU_seconds']}))

if __name__=='__main__':main()
