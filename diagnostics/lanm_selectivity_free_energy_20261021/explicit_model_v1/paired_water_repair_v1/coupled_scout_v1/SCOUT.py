#!/usr/bin/env python3
"""Four fixed, root-released coupled probes; no optimizer or model changes."""
import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import time

THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
CELLS = ('water_central_La', 'water_central_Dy', 'aquo8_central_La', 'aquo8_central_Dy')


def require(ok, message):
    if not ok: raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path']) == record, 'Changed evidence: '+record['path'])
    return Path(record['path'])


def read(record): return json.loads(verify(record).read_text())


def module(record, name):
    spec = importlib.util.spec_from_file_location(name, verify(record)); value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value); return value


def save(path, obj):
    p = Path(path); temp = p.with_name(p.name+'.tmp')
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n'); os.replace(temp, p)


def fingerprint(cfg):
    return hashlib.sha256(json.dumps({k:v for k,v in cfg.items() if k not in ('root_release','preparation_result')},
        sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def predicted_controls(parent, weighted):
    # Exact original block/body weights, reconstructed without new molecular arrays.
    score = parent['selected_score']; offset = 1; controls = {}
    for kind in ('La', 'Dy', 'paired'):
        begin, raw = offset, []
        for q in ('water_central', 'aquo8_central'):
            bodies = score['per_fixture'][kind][q]['model_body_response']
            cats = sorted(set(b['category'] for b in bodies))
            for category in cats:
                subset = [b for b in bodies if b['category'] == category]
                for body in subset:
                    terms = 1 if category == 'ion' else 2
                    weight = .5*math.sqrt(1/(len(cats)*len(subset)*terms*3))/2/math.sqrt(2)
                    for _ in range(terms*3): raw.append(weighted[offset]/weight); offset += 1
        rms, maximum = 4*math.sqrt(sum(v*v for v in weighted[begin:offset])), max(abs(v) for v in raw)
        controls[kind] = {'rigid_RMS_kcal_mol_A': rms, 'maximum_component_kcal_mol_A': maximum,
                          'predicted_pass': rms <= 2 and maximum <= 5}
    require(offset == len(weighted), 'Predicted residual coverage')
    energy = weighted[0]*2
    return {'objective':sum(v*v for v in weighted), 'balanced_energy_error_kcal_mol':energy,
            'controls':controls, 'all_predicted_controls_pass':abs(energy)<=1 and all(v['predicted_pass'] for v in controls.values())}


def inputs(cfg, out):
    require(cfg['implementation'] == pin(__file__) and all(os.environ.get(k)=='1' for k in THREADS), 'Source/threads differ')
    verify(cfg['specification']); parent = read(cfg['parent_result']); diagnosis = read(cfg['diagnosis'])
    require(cfg['parent_result']['sha256']=='524d093fb920304e254369097a96843408ab21597aa7898f5d1a3947c64fe6b6' and
            cfg['diagnosis']['sha256']=='7943ce93d34ce7434d84e836e8e13639435564ab4c86e18a472eb8a25f505647', 'Wrong frozen parent')
    require(parent['status']=='finite_pilot_complete' and diagnosis['candidate_reselected'] is False and
            diagnosis['selected_candidate_unchanged']==parent['candidate_parameters'], 'Parent candidate differs')
    oldcfg=read(parent['config'])
    for key in ('evaluator','evaluator_config','objective','body_helper','water48'):
        require(cfg[key]==oldcfg[key], 'Original helper/config changed: '+key); verify(cfg[key])
    initial=parent['candidate_parameters']; widths=[p['upper']-p['lower'] for p in parent['family']['parameters']]
    vectors=[]
    for name,label in (('uniformly_scaled_Gauss_Newton','scaled_GN'),('bounded_linear_least_squares','bounded_LS')):
        direction=diagnosis[name]
        require(max(abs(v) for v in direction['normalized_delta'])<=.025, 'Direction outside declared radius')
        for fraction,suffix in ((.5,'half'),(1.,'full')):
            delta=[fraction*v for v in direction['normalized_delta']]
            values=[v+d*w for v,d,w in zip(initial,delta,widths)]
            require(all(p['lower']<=v<=p['upper'] for p,v in zip(parent['family']['parameters'],values)), 'Probe outside original bounds')
            weighted=[(1-fraction)*a+fraction*b for a,b in zip(parent['selected_score']['weighted_residuals'],direction['weighted_residual'])]
            vectors.append({'id':label+'_'+suffix,'direction':name,'fraction':fraction,'normalized_delta':delta,
                            'parameters':values,'linear_prediction':predicted_controls(parent,weighted)})
    require(len({tuple(v['parameters']) for v in vectors})==4, 'Probe duplication')
    # A new ledger is mandatory. Original complete native directories are only read.
    cache=out/'cache'; (cache/'cells').mkdir(parents=True,exist_ok=False)
    native_pins={}
    for evaluation in parent['search']['evaluations'].values():
        require(evaluation['status']=='complete','Parent vector incomplete')
        for row in evaluation['model'].values():
            p=row['provenance']
            if not p['reused_baseline']: native_pins[p['cache_key']]=p['native_authority']
    seeded=[]
    for key,record in sorted(native_pins.items()):
        path=verify(record); native=read(record)
        require(path.parent.name==key and native['status']=='complete' and native['cache_key']==key,'Parent native cache mismatch')
        require(native['cache_identity']['implementation']==cfg['evaluator'],'Parent evaluator changed')
        for p in native['input'].values():verify(p)
        for key2 in ('analyze_output','gradient_output'):verify(native['native'][key2])
        require(len(native['receipts'])==2 and all(v['returncode']==0 for v in native['receipts']),'Parent native receipts')
        target=cache/'cells'/key; target.symlink_to(path.parent,target_is_directory=True)
        seeded.append({'cache_key':key,'authority':record,'symlink':str(target),'read_only_use':True})
    require(len(seeded)==142,'Parent native-cell denominator')
    native=module(cfg['evaluator'],'coupled_scout_evaluator'); ecfg=read(cfg['evaluator_config'])
    ev=native.Evaluator(ecfg,cache,root_release=cfg.get('root_release'))
    original_gate=ev._gate
    def scoped_gate():
        original_gate()
        with ev.lock:
            require(ev._ledger['native_program_invocations']<32,'Four-vector scout call ceiling reached')
    ev._gate=scoped_gate
    objective=module(cfg['objective'],'coupled_scout_objective'); helper=module(cfg['body_helper'],'coupled_scout_body')
    water=read(cfg['water48']); panel=read(water['panel'])
    baseline=ev.evaluate(initial)
    for gid in CELLS:
        require(baseline[gid]['energy']==parent['selected_model'][gid]['energy'] and
                baseline[gid]['forces']==parent['selected_model'][gid]['forces'],'Selected baseline raw replay differs')
    score=objective.evaluate(ev.data,baseline,helper)
    require(score==parent['selected_score'],'Original candidate score replay differs')
    require(ev.snapshot()['native_program_invocations']==0,'Preparation made native calls')
    verify(parent['evaluator']['ledger']); verify(cfg['parent_result'])
    return parent,ev,objective,helper,water,panel,baseline,vectors,seeded


def main(stage,config_path,output):
    out=Path(output);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    children=resource.getrusage(resource.RUSAGE_CHILDREN);cfg=json.loads(Path(config_path).read_text());ev=None
    result={'status':'unavailable','stage':stage,'implementation':pin(__file__),'config':pin(config_path),
        'scientific_config_sha256':fingerprint(cfg),'rows':{k:{'status':'not_attempted','cells':{g:'not_attempted' for g in CELLS}}
            for k in ('scaled_GN_half','scaled_GN_full','bounded_LS_half','bounded_LS_full')},
        'new_native_program_invocations':0,'candidate_reselected':False,'physical_accuracy_qualified':False,
        'validation_outcomes_opened':False,'automatic_extension':False,'job_id':os.environ.get('SLURM_JOB_ID'),'error':None}
    try:
        if stage=='run':
            require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and
                    int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==4,'Exact four-CPU allocation required')
            prep=read(cfg['preparation_result']);release=read(cfg['root_release'])
            require(prep['status']=='prepared' and prep['new_native_program_invocations']==0 and
                    prep['scientific_config_sha256']==fingerprint(cfg),'Actual preparation differs')
            require(release['scout_status']=='approved_four_fixed_coupled_probes' and
                    release['scout_max_vectors']==4 and release['scout_max_native_program_invocations']==32 and release['scout_workers']==4 and
                    release['scout_config_sha256']==fingerprint(cfg),'Scoped root release differs')
            for key in ('implementation','specification','parent_result','diagnosis','preparation_result'):
                require(release['scout_'+key]==cfg[key],'Root scout binding differs: '+key)
        else:
            require(cfg.get('root_release') is None and cfg.get('preparation_result') is None,'Preparation must be unreleased')
        parent,ev,obj,helper,water,panel,baseline,vectors,seeded=inputs(cfg,out)
        result.update(parent_result=cfg['parent_result'],diagnosis=cfg['diagnosis'],vectors=vectors,seeded_cache=seeded,
            baseline_model=baseline,baseline_score=parent['selected_score'],baseline_parameters=parent['candidate_parameters'],
            original_ledger=parent['evaluator']['ledger'],evaluator=ev.snapshot())
        if stage=='prepare':result['status']='prepared'
        else:
            require(vectors==prep['vectors'],'Frozen four vectors differ')
            errors=[]
            def evaluate(vector):
                require(not ev.stop_event.is_set(),'Global stop active')
                model=ev.evaluate(vector['parameters']);score=obj.evaluate(ev.data,model,helper)
                return {'status':'complete','model':model,'score':score,
                    'cells':{k:'complete' for k in CELLS},'method_sensitivity':obj.method_sensitivity(ev.data,model,water,panel)}
            with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
                futures={pool.submit(evaluate,v):v for v in vectors}
                for future in concurrent.futures.as_completed(futures):
                    vector=futures[future];key=vector['id']
                    try:result['rows'][key]=future.result()
                    except Exception as exc:
                        ev.stop_event.set();result['rows'][key].update(status='unavailable',error=repr(exc));errors.append(key)
                        for other in futures:
                            if other is not future:other.cancel()
                    result['evaluator']=ev.snapshot();result['new_native_program_invocations']=result['evaluator']['native_program_invocations']
                    save(out/'RESULT.json',result)
            if errors:
                for vector in vectors:
                    if vector['id'] in errors:
                        requests=[r for r in ev.snapshot()['requests'] if r['vector']==vector['parameters']]
                        if requests:result['rows'][vector['id']]['cells']=requests[-1]['cells']
            require(not errors,'At least one probe unavailable; no selection')
            selected={'id':'original_candidate','parameters':parent['candidate_parameters'],'model':baseline,'score':parent['selected_score']}
            for vector in vectors:
                row=result['rows'][vector['id']]
                row['actual_minus_linear_objective']=row['score']['objective']-vector['linear_prediction']['objective']
                if row['score']['objective']<selected['score']['objective']:
                    selected={'id':vector['id'],'parameters':vector['parameters'],'model':row['model'],'score':row['score']}
            result.update(status='four_probe_scout_complete',selected=selected,candidate_reselected=selected['id']!='original_candidate',
                development_checks_passed=selected['score']['development_checks_passed'],automatic_validation_release=False)
        verify(cfg['parent_result']);verify(parent['evaluator']['ledger'])
        result['original_result_and_ledger_unchanged']=True
    except Exception as exc:result.update(status='failed_or_unavailable',error=repr(exc),candidate_reselected=False)
    finally:
        if ev is not None:result['evaluator']=ev.snapshot();result['new_native_program_invocations']=result['evaluator']['native_program_invocations']
        end=resource.getrusage(resource.RUSAGE_CHILDREN)
        result.update(elapsed_seconds=time.monotonic()-start,parent_process_CPU_seconds=time.process_time()-cpu,
            aggregate_native_child_user_CPU_seconds=end.ru_utime-children.ru_utime,
            aggregate_native_child_system_CPU_seconds=end.ru_stime-children.ru_stime)
        save(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'new_native_program_invocations':result['new_native_program_invocations'],'result':pin(out/'RESULT.json')}))
    return 0 if result['status'] in ('prepared','four_probe_scout_complete') else 1


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('prepare','run'));p.add_argument('--config',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();raise SystemExit(main(a.stage,a.config,a.output))
