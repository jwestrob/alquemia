#!/usr/bin/env python3
"""Two declared half-points,112 native maximum; no solver, optimization or selection."""
import argparse,copy,hashlib,importlib.util,json,math,os,resource,sys,time
from pathlib import Path
import NATIVE_ADAPTER as A
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/'local_comparison_v2'
sys.path.insert(0,str(PARENT))
spec=importlib.util.spec_from_file_location('radial_native_parent',PARENT/'RUN.py');P=importlib.util.module_from_spec(spec);spec.loader.exec_module(P)
pin,verify,read,save,digest,module,require=P.pin,P.verify,P.read,P.save,P.digest,P.module,P.require
THREADS=P.THREADS
LIMITS={'new_parameter_vectors':2,'new_native_programs':112,'reused_derivative_vectors':0,'new_center_calls':0,'workers':12,'solver_calls':0,'numerical_threads':1}
CASES=('original12','radial14')
def fingerprint(cfg):return digest({k:v for k,v in cfg.items() if k not in ('root_release','source_review','preparation','fingerprint')})

def metadata(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['adapter']==pin(A.__file__) and cfg['previous_runner']==pin(P.__file__),'Exact policy sources')
    require(cfg['limits']==LIMITS and cfg['root_declaration']==pin(HERE/'ROOT_DECLARATION.json'),'Exact finite scope')
    dec=read(cfg['root_declaration']);oldcfg=read(cfg['comparison_configuration']);prior=read(dec['producer']);affine=read(dec['affine_result'])
    require(dec['status']=='declared_two_fixed_half_direction_native_tests' and tuple(dec['points'])==CASES,'Exact two reported half-points')
    require(prior['status']=='completed_radial_local_comparison' and prior['error'] is None and prior['native_programs']==112 and prior['solver_calls']==2 and prior['configuration']==cfg['comparison_configuration'] and prior['affine_result']==dec['affine_result'],'Completed comparison authority')
    require(affine['status']=='completed_two_matched_radial_affine_comparisons' and affine['solver_calls']==2,'Two completed saved solves')
    # All inherited scientific/native configuration values are literal; runtime paths/policy change only.
    changed={'implementation','adapter','limits','run_id','output_directory','root_release','source_review','preparation','fingerprint'}
    require(all(cfg[k]==v for k,v in oldcfg.items() if k not in changed),'Inherited scientific configuration unchanged')
    points=[]
    for case in CASES:
        row=dec['points'][case];pred=next(x for x in affine['cases'][case]['predictions'] if x['fraction']==.5)
        require(row['prediction']==pred and row['fraction']==.5 and row['parameters']==pred['parameters'] and row['key']==pred['key']==digest(row['parameters']) and row['solver_status']==affine['cases'][case]['solver_status'],'Literal declared point/prediction '+case)
        require(pred['retained_strict_positive'] and pred['all25_predicted'] and pred['signed_W9_safeguard_pass'],'Declared predicted safeguards')
        points.append({'ordinal':len(points)+1,'case':case,'key':row['key'],'parameters':row['parameters'],'prediction':pred,'solver_status':row['solver_status']})
    require(dec['center']==cfg['center'] and points[0]['parameters'][12:]==[0.,0.] and len({x['key'] for x in points})==2,'Matched zero extension and unique points')
    domain=read(read(cfg['parent_configuration'])['domain']);center=read(cfg['center']);base=module(cfg['parent_evaluator'],'radial_metadata_base')
    for point in points:
        require(len(point['parameters'])==14 and all(math.isfinite(x) for x in point['parameters']),'Fourteen finite parameters')
        base.checked(point['parameters'][:12],domain);require(all(-1<=x<=1 for x in point['parameters'][12:]),'Radial box')
        for metal in ('La','Dy'):require(base.effective(point['parameters'][:12],metal,domain)!=base.effective(center['parameters'],metal,domain),'Both metals changed; no center-cell reuse')
    manifest=read(oldcfg['inputs'])['sources'];em=read(read(cfg['parent_configuration'])['manifest'])
    require(len(em['rows'])==46 and sum(1+int(r['needs_gradient']) for r in em['rows'].values())==56,'Exact46/56 denominator')
    return dec,oldcfg,points

def release_gate(cfg):
    dec,oldcfg,points=metadata(cfg);root=read(cfg['comparison_root_admission'])
    require(root['status']=='admitted_radial_local_comparison' and root['result']==dec['producer'] and root['affine_result']==dec['affine_result'],'Actual root comparison admission')
    actual=read(root['independent_review']);require(actual['status']=='passed_radial_local_comparison_actual_review' and actual['result']==dec['producer'] and actual['affine_result']==dec['affine_result'] and actual['native_candidate_admission'] is False,'Independent actual comparison review; no candidate admission')
    review=read(cfg['source_review']);prep=read(cfg['preparation']);release=read(cfg['root_release']);fp=fingerprint(cfg)
    require(review['status']=='passed_radial_native_comparison_source_review' and review['configuration_fingerprint']==fp,'Exact independent source review')
    require(prep['status']=='prepared_radial_native_comparison_zero_calls' and prep['configuration_fingerprint']==fp and prep['native_programs']==prep['solver_calls']==0,'Exact zero-call preparation')
    require(release['status']=='released_two_fixed_radial_native_comparisons' and release['configuration_fingerprint']==fp,'One exact root release')
    for k in ('implementation','adapter','previous_runner','root_declaration','comparison_configuration','comparison_root_admission','report_helper','limits','inputs','output_directory','root_cache_bridge','qualification_review','qualification_root_review','source_review','preparation'):
        require(release[k]==cfg[k],'Exact release '+k)

def summarize(answer,prediction,core,fixed,targets,center,helper):
    r=helper.summarize(answer,prediction,core,fixed,targets);metrics=r['metrics'];original=core.assess(center,fixed)
    electric={k:{'actual_score':metrics['original_block_scores'][k],'center_score':original['original_block_scores'][k],'ratio':metrics['original_block_scores'][k]/original['original_block_scores'][k],'nonincrease':metrics['original_block_scores'][k]<=original['original_block_scores'][k]} for k in ('pair_mu','C_mu','chi')}
    C=center['response']['C_dipoles_Debye'];ref=r['W9_C_Dy_minus_La_Debye']['reference'][2]
    ceiling=abs(math.fsum([C['aquo9_central_Dy'][2],-C['aquo9_central_La'][2],-ref]));actual=math.fsum([answer['response']['C_dipoles_Debye']['aquo9_central_Dy'][2],-answer['response']['C_dipoles_Debye']['aquo9_central_La'][2],-ref])
    require(ceiling==prediction['signed_W9_z_absolute_ceiling_Debye'],'Exact original signed ceiling')
    signed={'actual_model_minus_reference_Debye':actual,'absolute_ceiling_Debye':ceiling,'passed':abs(actual)<=ceiling,'predicted_residual_Debye':prediction['signed_W9_z_residual_Debye'],'actual_minus_predicted_Debye':actual-prediction['signed_W9_z_residual_Debye']}
    r.update(electrical_norm_nonincrease=electric,signed_W9_z_safeguard=signed,all_declared_safeguards_pass=metrics['all25_original_physical_gates'] and all(x['nonincrease'] for x in electric.values()) and signed['passed'],native_candidate_admission=False)
    return r

def run(mode,config_path,out):
    tick,cpu=time.monotonic(),time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN);cfg=json.loads(Path(config_path).read_text());cfg['fingerprint']=fingerprint(cfg)
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False);ev=None
    result={'status':'started','mode':mode,'implementation':pin(__file__),'configuration':pin(config_path),'configuration_fingerprint':cfg['fingerprint'],'limits':LIMITS,'native_programs':0,'solver_calls':0,'vectors':{},'error':None,'selected_candidate':None,'native_candidate_admission':False,'reference_calls':0,'heldout_reads':0,'job_id':None}
    def persist(terminal=False):P.persist_result(out,result,ev,terminal)
    try:
        dec,parent_cfg,points=metadata(cfg);result['root_declaration']=cfg['root_declaration'];result['points']=[{k:v for k,v in row.items() if k!='prediction'} for row in points]
        require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
        if mode=='prepare':result.update(status='prepared_radial_native_comparison_zero_calls',inherited_preparation=parent_cfg['preparation'],new_cache_constructed=False,old46_reparsed=False,old24_replayed=False,expected_native_programs=112)
        else:
            require(str(out)==cfg['output_directory'],'Exact fresh output');release_gate(cfg)
            s,base,old,original,attr,checks,parent,center,inv=P.context(parent_cfg)
            evaluator=A.factory(P.A)(old,base,checks,release_gate)
            ev=evaluator(cfg,out/'cache','run',{'new_probes':points},center,parent)
            declaration=read(s['campaign_declaration']);fixed=original.make_fixed(read(declaration['start']['artifact']),ev.m,read(s['decomposition']),ev.domain,declaration)
            fixed['domain']=P.M.extend_domain(fixed['domain']);core=P.M.core14(original);zero=P.M.view(center,cfg['center'],cfg['root_cache_bridge']);helper=module(cfg['report_helper'],'native_saved_reporting')
            targets=read(cfg['targets']);result['center_metrics']=core.assess(zero,fixed);persist();(out/'vectors').mkdir()
            for row in points:
                rec={'case':row['case'],'status':'started','parameters':row['parameters'],'fraction':.5,'solver_status_retained':row['solver_status']};result['vectors'][row['key']]=rec;persist()
                try:
                    answer=ev.evaluate(row['parameters']);p=out/'vectors'/(row['key']+'.json');save(p,answer);rec.update(status='complete',artifact=pin(p),comparison=summarize(answer,row['prediction'],core,fixed,targets,zero,helper));persist()
                except Exception as exc:rec.update(status='failed',error=repr(exc));raise
            require(ev.ledger['native_program_invocations']==112 and len(result['vectors'])==2,'Exact two-vector112 denominator')
            result['status']='completed_two_fixed_radial_native_comparisons'
    except Exception as exc:
        result.update(status='failed_radial_native_comparison' if mode=='run' else 'failed_preparation',error=repr(exc))
        if ev is not None:
            try:ev._stop(exc)
            except Exception as other:result['stop_reporting_error']=repr(other)
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':ch.ru_utime-ch0.ru_utime,'child_system_CPU_seconds':ch.ru_stime-ch0.ru_stime,'execution_mode':'login','native_workers':12,'numerical_threads':1,'scheduler_allocated_CPU_seconds':0};persist(True)
    print(json.dumps({'status':result['status'],'error':result['error'],'native_programs':result['native_programs']}));return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('prepare','run'));p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();raise SystemExit(bool(run(a.mode,a.config,a.out)['error']))
