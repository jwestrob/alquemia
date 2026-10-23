#!/usr/bin/env python3
"""Four radial derivative vectors and two matched affine comparisons; no fit loop."""
import argparse,copy,hashlib,importlib.util,json,os,resource,socket,sys,time
from pathlib import Path
import ADAPTER as A
import MATH as M
HERE=Path(__file__).resolve().parent
LIMITS={'new_parameter_vectors':4,'new_native_programs':112,'reused_derivative_vectors':24,'new_center_calls':0,'workers':12,'solver_calls':2,'numerical_threads':1}
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')

def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(isinstance(p,dict) and pin(p['path'])==p,'Changed/unbound source');return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def fingerprint(c):return digest({k:v for k,v in c.items() if k not in ('root_release','source_review','preparation','fingerprint')})
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def context(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['adapter']==pin(A.__file__) and cfg['math']==pin(M.__file__) and cfg['solver']==pin(HERE/'SOLVE.py'),'Frozen new code')
    require(cfg['limits']==LIMITS and cfg['threads']==list(THREADS) and all(os.environ.get(k)=='1' for k in THREADS),'Finite limits/numerical threads')
    require(not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID') and cfg['hostname']==socket.gethostname(),'Truthful login context')
    require(cfg['specification']==pin(HERE/'SPEC.md'),'Frozen implementation specification')
    inp=read(cfg['inputs']);s=inp['sources'];[verify(p) for p in s.values()]
    for name in ('parent_evaluator','parent_configuration','core','original_adapter','settings','targets','environment','center','native_source_manifest','build'):
        require(cfg[name]==s[name],'Direct/source binding '+name)
    require(cfg['reuse_inventory']==cfg['inputs'],'Immutable evidence alias')
    bridge=read(s['bridge_assessment']);require(bridge['status']=='conditional_radial_zero_cache_bridge_assessment' and not bridge['cache_bridge_admitted'],'Conditional assessment retained')
    require(bridge['old_center']==cfg['center'] and bridge['new_source_manifest']==cfg['native_source_manifest'] and bridge['new_build']==cfg['build'] and bridge['new_executables']==cfg['executables'],'Zero bridge native/center authority')
    parent=read(cfg['parent_configuration']);base=module(cfg['parent_evaluator'],'radial_original_evaluator');old=module(cfg['original_adapter'],'radial_original_adapter')
    require(parent['implementation']==cfg['parent_evaluator'] and old.J.pin(old.J.__file__)==s['journal'],'Exact evaluator/journal source')
    original=module(cfg['core'],'radial_original_core');attr=module(s['attribution_math'],'radial_original_attribution_math')
    sys.path.insert(0,str(Path(s['qualification_checks']['path']).parent));checks=module(s['qualification_checks'],'radial_admitted_readback')
    center=read(cfg['center']);require(digest(center['parameters'])==bridge['old_center_key'] and center['status']=='complete','Exact finalfeca')
    prior=read(bridge['old24_result']);inv=read(bridge['old24_inventory']);require(prior['status']=='completed_final_local_attribution_panel' and prior['native_programs']==672 and prior['solver_calls']==4,'Complete inherited final24')
    require(cfg['inherited_actual_spans']==[r['actual_normalized_span'] for r in bridge['actual_normalized_spans']],'Literal old secant denominators')
    rows=bridge['probes_in_declared_order'];require(len(rows)==24 and [r['old_parameters'] for r in rows]==original.probe_points(center['parameters'],read(parent['domain'])),'Literal old24 order and clipped spans')
    for r in rows:
        require(prior['vectors'][r['old_key']]['artifact']==r['old_artifact'] and r['extended_parameters_view']==r['old_parameters']+[0.,0.],'Exact inherited artifacts/views');verify(r['old_artifact'])
    new=[]
    for i,metal in enumerate(('La','Dy')):
        for sign in (-1,1):
            v=center['parameters']+[0.,0.];v[12+i]=sign*.01;new.append({'ordinal':len(new)+1,'metal':metal,'sign':sign,'parameters':v,'key':digest(v),'normalized_span':.01,'new_native_programs':28})
    env=read(cfg['environment']);require(str(Path(sys.prefix).resolve())==env['prefix'],'Pinned existing numerical environment')
    import importlib.metadata
    for name,version in env['versions'].items():require(importlib.metadata.version(name)==version,'Exact package '+name)
    for p in env['files']:verify(p)
    build=read(cfg['build']);require(build['status']=='complete' and build['manifest']==cfg['native_source_manifest'] and build['executables']==cfg['executables'],'Actual radial build');[verify(p) for p in cfg['executables'].values()]
    inventory={'old_probes':rows,'old_normalized_spans':bridge['actual_normalized_spans'],'new_probes':new,'bridge_assessment':s['bridge_assessment']}
    bridge_gate(cfg)
    return s,base,old,original,attr,checks,parent,center,inventory

def bridge_gate(cfg):
    bridge=read(cfg['root_cache_bridge']);require(bridge['status']=='admitted_radial_zero_cache_equivalence' and bridge['root_approval'] and bridge['cache_bridge_admitted'],'Explicit root zero-cache bridge')
    sources=read(cfg['inputs'])['sources'];assessment=read(sources['bridge_assessment'])
    require(bridge['conditional_assessment']==sources['bridge_assessment'] and bridge['old_center']==cfg['center'] and bridge['new_source_manifest']==cfg['native_source_manifest'] and bridge['new_build']==cfg['build'] and bridge['new_executables']==cfg['executables'],'Exact bridge scientific/native authority')
    require(bridge['probes_in_declared_order']==assessment['probes_in_declared_order'] and bridge['actual_normalized_spans']==assessment['actual_normalized_spans'],'Exact ordered old24 bridge')
    require(bridge['independent_actual_review']==cfg['qualification_review'] and bridge['root_qualification_admission']==cfg['qualification_root_review'],'Exact complete qualification admissions')
    q=read(cfg['qualification_review']);root=read(cfg['qualification_root_review'])
    require(q['status']=='passed_radial_drive_qualification_actual_review' and q['numerical_checks_passed'] and q['zero_parity_passed'],'Actual radial49/132 numerical qualification')
    require(root['status']=='admitted_radial_drive_numerical_qualification' and root['result']==q['result']==bridge['qualification_result'] and root['independent_actual_review']==cfg['qualification_review'] and root['counts']['native_programs']==49 and root['counts']['internal_FD_energies']==132,'Root exact finite qualification')
    result=read(q['result']);require(result['source_manifest']==cfg['native_source_manifest'] and result['build']==cfg['build'] and result['executables']==cfg['executables'],'Admitted actual native source/executables')

def release_gate(cfg):
    bridge_gate(cfg)
    review=read(cfg['source_review']);prep=read(cfg['preparation']);release=read(cfg['root_release']);fp=fingerprint(cfg)
    require(review['status']=='passed_radial_local_comparison_source_review' and review['configuration_fingerprint']==fp,'Exact independent source review')
    require(prep['status']=='prepared_radial_local_comparison_zero_calls' and prep['configuration_fingerprint']==fp and prep['native_programs']==prep['solver_calls']==0,'Exact zero-call preparation')
    require(release['status']=='released_four_radial_probes_two_affine_comparisons' and release['configuration_fingerprint']==fp,'Exact finite root release')
    for k in ('implementation','adapter','math','solver','inputs','limits','output_directory','root_cache_bridge','qualification_review','qualification_root_review','source_review','preparation'):
        require(release[k]==cfg[k],'Release field '+k)


def run(mode,config_path,out):
    tick,cpu=time.monotonic(),time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN);cfg=json.loads(Path(config_path).read_text());cfg['fingerprint']=fingerprint(cfg)
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False);ev=None
    result={'status':'started','mode':mode,'implementation':pin(__file__),'configuration':pin(config_path),'configuration_fingerprint':cfg['fingerprint'],'limits':LIMITS,'native_programs':0,'solver_calls':0,'vectors':{},'error':None,'selected_candidate':None,'native_trials':0,'reference_calls':0,'heldout_reads':0,'job_id':None}
    def persist():
        if ev is not None:result['evaluator']=ev.snapshot();result['native_programs']=ev.ledger['native_program_invocations']
        save(out/'RESULT.tmp',result);os.replace(out/'RESULT.tmp',out/'RESULT.json')
    try:
        s,base,old,original,attr0,checks,parent,center,inv=context(cfg)
        if mode=='run':require(str(out)==cfg['output_directory'],'Exact fresh run directory');release_gate(cfg)
        ev=A.evaluator_class(old,base,checks,release_gate)(cfg,out/'cache',mode,inv,center,parent)
        save(out/'INVENTORY.json',inv);result['inventory']=pin(out/'INVENTORY.json');result['center']=cfg['center'];result['root_cache_bridge']=cfg.get('root_cache_bridge')
        decl=read(s['campaign_declaration']);fixed=original.make_fixed(read(decl['start']['artifact']),ev.m,read(s['decomposition']),ev.domain,decl)
        fixed['domain']=M.extend_domain(fixed['domain']);core=M.core14(original);attr=M.attribution14(attr0)
        zero=M.view(center,cfg['center'],cfg.get('root_cache_bridge'));result['center_metrics']=core.assess(zero,fixed)
        result['dimension_only_edits']={'core':core.dimension_edits,'attribution':attr.dimension_edits};persist()
        if mode=='prepare':
            result.update(status='prepared_radial_local_comparison_zero_calls',preparation_checks={'old24_metadata_pins_only':True,'old46_center_import_exact':True,'four_new_probes':True,'112_after_exact_unchanged_metal_reuse':True,'new_center_calls':0,'actual_linearization_or_solve':False})
        else:
            probes=[M.view(read(r['old_artifact']),r['old_artifact'],cfg['root_cache_bridge']) for r in inv['old_probes']];(out/'vectors').mkdir()
            for row in inv['new_probes']:
                rec={'status':'started','parameters':row['parameters'],'ordinal':row['ordinal']};result['vectors'][row['key']]=rec;persist()
                try:
                    answer=ev.evaluate(row['parameters']);p=out/'vectors'/(row['key']+'.json');save(p,answer);rec.update(status='complete',artifact=pin(p));probes.append(answer);persist()
                except Exception as e:rec.update(status='failed',error=repr(e));raise
            require(len(probes)==28 and ev.ledger['native_program_invocations']==112,'Exact derivative denominator before two solves')
            solver=module(cfg['solver'],'radial_two_case_solver');result['affine_diagnostics']=solver.run(core,attr,zero,probes,fixed,cfg,out/'affine_cases');result['solver_calls']=result['affine_diagnostics']['solver_calls'];require(result['solver_calls']==2,'Exactly two attempted affine comparisons')
            result['status']='completed_radial_local_comparison'
    except Exception as e:
        if ev is not None:ev._stop(e)
        result.update(status='failed_radial_local_comparison' if mode=='run' else 'failed_preparation',error=repr(e))
    finally:
        if (out/'affine_cases/RESULT.json').exists():result['affine_result']=pin(out/'affine_cases/RESULT.json');result['solver_calls']=json.loads((out/'affine_cases/RESULT.json').read_text())['solver_calls']
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':ch.ru_utime-ch0.ru_utime,'child_system_CPU_seconds':ch.ru_stime-ch0.ru_stime,'native_workers':12,'numerical_threads':1,'execution_mode':'login','scheduler_allocated_CPU_seconds':0,'inherited_costs_not_added':True};persist()
    print(json.dumps({'status':result['status'],'native_programs':result['native_programs'],'solver_calls':result['solver_calls'],'error':result['error']}));return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('prepare','run'));p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();raise SystemExit(bool(run(a.mode,a.config,a.out)['error']))
