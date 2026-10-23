#!/usr/bin/env python3
"""Explicit synthetic alpha responses and immutable consumed-data unit checks."""
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
WORK=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_coordination_joint_v1/metal_alpha_response_v1'
GEOMS=('water_minus','water_central','water_plus','aquo8_central','aquo9_central')
ENERGY_NAMES=('exchange_water_minus','exchange_water_central','exchange_water_plus','exchange_W9','work_La_minus','work_La_plus','work_Dy_minus','work_Dy_plus','paired_work_minus','paired_work_plus')

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def dump(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
def decode(v):
    """Independent inverse factors rather than the analyzer's weight division."""
    v=np.array(v,dtype=float);e=v[:8]*np.sqrt(32);energy=dict(zip(ENERGY_NAMES[:8],e.tolist()))
    energy.update(paired_work_minus=float(e[6]-e[4]),paired_work_plus=float(e[7]-e[5]))
    per={};index=8
    for kind in ('La','Dy','paired'):
        per[kind]={}
        for name,nw in zip(GEOMS,(1,1,1,8,9)):
            n=3+6*nw;segment=v[index:index+n];index+=n
            raw=segment*np.r_[np.full(3,np.sqrt(480)),np.full(6*nw,np.sqrt(960*nw))]
            rms=float(np.sqrt(80)*np.linalg.norm(segment));mx=float(np.max(np.abs(raw)))
            per[kind][name]={'rigid_RMS_kcal_mol_A':rms,'maximum_body_component_error_kcal_mol_A':mx,'passed':rms<=2 and mx<=5}
    assert index==413
    passed=sum(abs(x)<=1 for x in energy.values())+sum(x['passed'] for g in per.values() for x in g.values())
    return energy,per,passed
def score(v):
    v=np.array(v,dtype=float);energy,per,passed=decode(v)
    return {'weighted_residuals':v.tolist(),'objective':float(v@v),'development_checks_passed':passed==25,'per_fixture':per}

def main():
    start,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    out=WORK/'analysis_controls_v2';out.mkdir(parents=True,exist_ok=False)
    source=HERE/'ANALYZE.py';source_pin=pin(source)
    spec=importlib.util.spec_from_file_location('reviewed_alpha_analyzer',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    previous=WORK/'analysis_controls_v1/RESULT.json';previous_pin=pin(previous)
    manifest_path=WORK/'preparation_v2/MANIFEST.json';manifest=json.loads(manifest_path.read_text());data_path=Path(manifest['data']['path']);assert pin(data_path)==manifest['data'];actual=json.loads(data_path.read_text())
    checks=[]
    def check(name,condition):
        assert condition,name
        checks.append(name)
    def reject(name,fn):
        try:fn()
        except (ValueError,KeyError,TypeError,np.linalg.LinAlgError):checks.append(name)
        else:raise AssertionError(name)
    try:
        check('Actual preparation is zero-call and exactly eight probes',manifest['new_native_calls']==0 and len(manifest['probes'])==8 and {p['id'] for p in manifest['probes']}==set(m.PROBES))
        actual_check=m.linear_checks(actual['baseline_score']['weighted_residuals']);energy,force,count=decode(actual['baseline_score']['weighted_residuals'])
        check('Actual413 control counts and scientific failure',actual_check['declared_checks']==25 and actual_check['passed_checks']==count==22 and not actual_check['all_checks_passed'])
        for name,value in energy.items():check('Actual energy reconstruction '+name,abs(value-actual_check['energy_errors'][name])<1e-12)
        names=[('balanced_exchange',q) for q in GEOMS[:3]+('aquo9_central',)]+[('individual_radial_work',m+'_'+s) for m in ('La','Dy') for s in ('minus','plus')]+[('paired_radial_work',s) for s in ('minus','plus')]
        for name,(group,key) in zip(ENERGY_NAMES,names):check('Actual signed energy gate '+name,abs(actual_check['energy_errors'][name]-actual['baseline_score']['energy_errors_kcal_mol'][group][key])<1e-12)
        for kind in force:
            for geom,x in force[kind].items():
                got=actual_check['force_controls'][kind][geom];saved=actual['baseline_score']['per_fixture'][kind][geom]
                check('Actual force units/gate '+kind+'/'+geom,abs(got['RMS']-x['rigid_RMS_kcal_mol_A'])<1e-12 and abs(got['maximum']-x['maximum_body_component_error_kcal_mol_A'])<2e-12 and abs(got['RMS']-saved['rigid_RMS_kcal_mol_A'])<1e-12 and got['passed']==saved['passed'])
        probe_meta={}
        for p in manifest['probes']:
            metal=p['metal'];eta=(p['alpha']-m.BASE_ALPHA[metal])/m.BASE_ALPHA[metal]
            check('Actual rounded token eta '+p['id'],eta==p['eta'] and abs(eta-p['fraction'])<1e-10)
            check('Five actual prepared key tokens '+p['id'],len(p['inputs'])==5)
            for gid,row in p['inputs'].items():
                kp=row['files']['key'];check('Pinned prepared key '+p['id']+'/'+gid,pin(kp['path'])==kp)
                lines=[s.split() for s in Path(kp['path']).read_text().splitlines() if s.split() and s.split()[0].lower()=='polarize' and s.split()[1]==('500' if metal=='La' else '501')]
                check('Numeric metal-alpha token '+p['id']+'/'+gid,len(lines)==1 and float(lines[0][2])==p['alpha'])
            probe_meta[p['id']]={k:p[k] for k in ('id','metal','fraction','eta','alpha')}
        # This compatibility fixture has explicit synthetic zero response. It
        # contains no native probe result or measured polarizability response.
        compat={'status':'metal_alpha_diagnostic_complete','native_program_invocations':80,'error':None,'global_stop':False,'records':{p:dict(meta,status='complete',score=copy.deepcopy(actual['baseline_score'])) for p,meta in probe_meta.items()}}
        compatibility=m.analyze(actual,compat)
        check('Actual preparation-v2 schema and eta accepted with synthetic zero response',compatibility['status']=='saved_alpha_response_analysis_complete' and all(np.linalg.norm(v['alpha_columns'])==0 for v in compatibility['augmented_projections'].values()))
        dump(out/'SYNTHETIC_actual_metadata_zero_response.json',compat)

        # Orthogonal exact synthetic linear fixture: old six basis directions,
        # two new alpha directions, and an irreducible ninth-coordinate residual.
        r=np.zeros(413);r[0]=.01;r[6]=.02;r[7]=-.03;r[8]=.005
        J=np.zeros((413,6));J[:6]=np.eye(6)
        cols=np.zeros((413,2));cols[6,0]=1;cols[7,1]=1
        family={'parameters':[{'name':'synthetic_'+str(i),'lower':0.,'upper':1.} for i in range(6)]}
        d={'baseline_score':score(r),'six_column_Jacobian':J.tolist(),'family':family,'selected_normalized_parameters':[.5]*6}
        def response_fixture(linear=cols,cubic=None,quadratic=None):
            result={'status':'metal_alpha_diagnostic_complete','native_program_invocations':80,'error':None,'global_stop':False,'records':{}}
            for p,meta in probe_meta.items():
                j=0 if meta['metal']=='La' else 1;e=meta['eta'];v=r+e*linear[:,j]
                if cubic is not None:v=v+e**3*cubic[:,j]
                if quadratic is not None:v=v+e**2*quadratic[:,j]
                result['records'][p]=dict(meta,status='complete',score=score(v))
            return result
        syn=response_fixture();a=m.analyze(d,syn)
        check('Synthetic rank six to eight',a['original_six_projection']['rank']==6 and all(v['rank']==8 for v in a['augmented_projections'].values()))
        check('Synthetic irreducible loss and original loss',abs(a['original_six_projection']['remaining_weighted_loss']-.001325)<1e-14 and all(abs(v['remaining_weighted_loss']-.000025)<1e-14 for v in a['augmented_projections'].values()))
        for scale,x in a['augmented_projections'].items():
            check('Synthetic exact step and alpha coordinates '+scale,np.allclose(x['step'],[-.01,0,0,0,0,0,-.02,.03],atol=1e-12,rtol=0) and np.allclose(x['implied_six_parameters'],[.49,.5,.5,.5,.5,.5],atol=1e-12,rtol=0))
            check('Synthetic flags retained '+scale,x['prediction_only'] and not x['candidate_selected'] and x['original_six_bounds_satisfied'] and x['original_six_trust_satisfied'] and x['alpha_positive'])
        check('Measured range flags distinguish1 versus5 percent',not a['augmented_projections']['1percent']['alpha_step_within_measured_range'] and a['augmented_projections']['5percent']['alpha_step_within_measured_range'])
        check('Linear scale independence',all(v['relative_column_change']<1e-12 and abs(v['cosine']-1)<1e-12 for v in a['scale_dependence'].values()))
        check('Zero alpha columns preserve six-column projection',m.analyze(d,response_fixture(np.zeros((413,2))))['augmented_projections']['1percent']['rank']==6)
        cubic=cols*np.array([40.,-20.]);quadratic=np.zeros((413,2));quadratic[9]=[.3,-.4]
        nonlinear=m.analyze(d,response_fixture(cols,cubic,quadratic))
        for suffix in ('1','5'):
            expected=[]
            for j,metal in enumerate(('La','Dy')):
                p,q=probe_meta[metal+'_plus'+suffix]['eta'],probe_meta[metal+'_minus'+suffix]['eta'];expected.append(cols[:,j]+(p*p+p*q+q*q)*cubic[:,j]+(p+q)*quadratic[:,j])
            check('Actual-spacing finite central derivative '+suffix,np.allclose(nonlinear['augmented_projections'][suffix+'percent']['alpha_columns'],np.column_stack(expected),atol=1e-12,rtol=1e-12))
        check('Nonlinear scale difference remains visible',all(v['relative_column_change']>.04 for v in nonlinear['scale_dependence'].values()))
        matrix=np.column_stack((J,cols));scaled=matrix.copy();scaled[:,6:]*=10
        unscaled=m.projection(matrix,r);rescaled=m.projection(scaled,r)
        check('Coordinate scaling changes steps inversely without changing span',np.allclose(np.array(rescaled['step'])[6:]*10,np.array(unscaled['step'])[6:],atol=1e-12) and abs(rescaled['remaining_weighted_loss']-unscaled['remaining_weighted_loss'])<1e-14)
        cutoff=matrix.copy();cutoff[:,6]*=.5e-6;cutoff[:,7]*=2e-6
        cut=m.projection(cutoff,r)
        check('Frozen relative singular cutoff drops only subcutoff column',cut['rank']==7 and cut['relative_cutoff']==1e-6 and abs(cut['absolute_cutoff']-1e-6)<1e-15 and abs(cut['remaining_weighted_loss']-(.02**2+.005**2))<1e-12)
        check('Zero-rank matrix remains explicit',m.projection(np.zeros((413,6)),r)['rank']==0)

        # Independent threshold isolation: energy, RMS, and max-component gate.
        for x,passed in ((.9999,True),(1.0001,False)):
            v=np.zeros(413);v[3]=x/np.sqrt(32);check('Isolated energy threshold '+str(x),m.linear_checks(v)['all_checks_passed']==passed)
        for x,passed in ((1.9999,True),(2.0001,False)):
            v=np.zeros(413);v[8]=x/np.sqrt(80);check('Isolated RMS threshold '+str(x),m.linear_checks(v)['force_controls']['La']['water_minus']['passed']==passed)
        for x,passed in ((4.9999,True),(5.0001,False)):
            v=np.zeros(413);v[89]=x/np.sqrt(960*9);check('Isolated W9 body-component threshold '+str(x),m.linear_checks(v)['force_controls']['La']['aquo9_central']['passed']==passed)
        for label,value in [('short',np.zeros(412)),('NaN',np.full(413,np.nan)),('Inf',np.full(413,np.inf))]:reject('Invalid413 residual '+label,lambda value=value:m.linear_checks(value))
        def bad_data(name,edit):
            value=copy.deepcopy(d);edit(value);reject(name,lambda:m.analyze(value,syn))
        def bad_result(name,edit):
            value=copy.deepcopy(syn);edit(value);reject(name,lambda:m.analyze(d,value))
        bad_data('Nonfinite six-column Jacobian',lambda x:x['six_column_Jacobian'][0].__setitem__(0,float('nan')))
        bad_data('Wrong Jacobian dimensions',lambda x:x.__setitem__('six_column_Jacobian',np.zeros((412,6)).tolist()))
        bad_data('Baseline objective mismatch',lambda x:x['baseline_score'].__setitem__('objective',123.))
        bad_data('Baseline force-unit mismatch',lambda x:x['baseline_score']['per_fixture']['La']['aquo9_central'].__setitem__('rigid_RMS_kcal_mol_A',123.))
        for key,value in [('lower',float('nan')),('upper',float('inf')),('upper',0.)]:bad_data('Invalid parameter '+key+str(value),lambda x,key=key,value=value:x['family']['parameters'][0].__setitem__(key,value))
        bad_data('Wrong parameter count',lambda x:x['family']['parameters'].pop())
        for value in (float('nan'),float('inf'),-.01,1.01):bad_data('Invalid normalized origin '+str(value),lambda x,value=value:x['selected_normalized_parameters'].__setitem__(0,value))
        bad_result('Missing probe retained as failure',lambda x:x['records'].pop('La_plus5'))
        bad_result('Incomplete probe rejected',lambda x:x['records']['Dy_minus1'].__setitem__('status','failed'))
        bad_result('Global stop rejected',lambda x:x.__setitem__('global_stop',True))
        bad_result('Native error rejected',lambda x:x.__setitem__('error','explicit synthetic native failure'))
        bad_result('Nonfinite probe rejected',lambda x:x['records']['La_plus1']['score']['weighted_residuals'].__setitem__(0,float('nan')))
        bad_result('Probe objective mismatch rejected',lambda x:x['records']['La_plus1']['score'].__setitem__('objective',123.))
        bad_result('Probe gate mismatch rejected',lambda x:x['records']['La_plus1']['score'].__setitem__('development_checks_passed',False))
        bad_result('Wrong declared fraction rejected',lambda x:x['records']['La_plus1'].__setitem__('fraction',.02))
        bad_result('Wrong alpha token rejected',lambda x:x['records']['La_plus1'].__setitem__('alpha',1.7))
        bad_result('Wrong actual eta rejected',lambda x:x['records']['La_plus1'].__setitem__('eta',.02))
        bad_result('Nonfinite alpha token rejected',lambda x:x['records']['Dy_minus1'].__setitem__('alpha',float('inf')))
        bad_result('Probe metal mismatch rejected',lambda x:x['records']['Dy_minus1'].__setitem__('metal','La'))
        # Exercise the actual CLI entry's complete/80 gates using explicit
        # synthetic files, not a native run or a real result path.
        sd=out/'SYNTHETIC_DATA.json';dump(sd,d);sm=out/'SYNTHETIC_MANIFEST.json';dump(sm,{'data':pin(sd)});full=copy.deepcopy(syn);full['manifest']=pin(sm)
        for name,key,value in [('running','status','running'),('79programs','native_program_invocations',79),('81programs','native_program_invocations',81)]:
            value_result=copy.deepcopy(full);value_result[key]=value;p=out/('SYNTHETIC_'+name+'.json');dump(p,value_result);dest=out/('REJECTED_'+name)
            reject('Main complete80 admission '+name,lambda p=p,dest=dest:m.main(p,dest));check('Rejected main creates no analysis '+name,not dest.exists())
        complete=out/'SYNTHETIC_COMPLETE80.json';dump(complete,full);m.main(complete,out/'SYNTHETIC_ANALYSIS')
        emitted=json.loads((out/'SYNTHETIC_ANALYSIS/RESULT.json').read_text());check('Complete synthetic80 CLI persists rank8 prediction only',emitted['augmented_projections']['1percent']['rank']==8 and not emitted['candidate_selected'] and emitted['new_native_calls']==0)
        check('Analyzer source unchanged through controls',pin(source)==source_pin)
        check('Prior analysis controls remain immutable',pin(previous)==previous_pin)
        result={'status':'passed','checks_passed':len(checks),'checks':checks,'implementation':pin(__file__),'analyzer':source_pin,'actual_prepared_manifest':pin(manifest_path),'actual_prepared_data':pin(data_path),'prior_controls_unchanged':previous_pin,'synthetic_response_fixtures_explicit':True,'actual_alpha_response_outputs_opened':False,'native_calls':0,'parameter_updates':0,'threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    except Exception as exc:
        result={'status':'failed','error':repr(exc),'checks_passed':len(checks),'checks':checks,'implementation':pin(__file__),'analyzer':source_pin,'native_calls':0,'actual_alpha_response_outputs_opened':False,'threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    dump(out/'RESULT.json',result)
    print(json.dumps({'result':pin(out/'RESULT.json'),'status':result['status'],'checks_passed':len(checks),'error':result.get('error'),'wall_seconds':result['wall_seconds'],'CPU_seconds':result['CPU_seconds']}))
    return 0 if result['status']=='passed' else 1

if __name__=='__main__':raise SystemExit(main())
