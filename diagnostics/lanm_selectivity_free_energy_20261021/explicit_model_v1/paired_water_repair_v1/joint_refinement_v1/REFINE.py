#!/usr/bin/env python3
"""Twelve bounded joint steps of the unchanged four-parameter water objective."""
import argparse
import concurrent.futures
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import resource
import time
import numpy as np

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
CELLS=('water_central_La','water_central_Dy','aquo8_central_La','aquo8_central_Dy')


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path'])==record,'Changed evidence: '+record['path']);return Path(record['path'])


def read(record):return json.loads(verify(record).read_text())


def token(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def fingerprint(cfg):return token({k:v for k,v in cfg.items() if k not in ('root_release','preparation_result')})


def module(record,name):
    s=importlib.util.spec_from_file_location(name,verify(record));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def save(path,value):
    p=Path(path);temp=p.with_name(p.name+'.tmp');temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(temp,p)


def prepare(cfg,out):
    require(cfg['implementation']==pin(__file__) and all(os.environ.get(k)=='1' for k in THREADS),'Source/threads differ')
    verify(cfg['specification']);scout=read(cfg['scout_result'])
    require(cfg['scout_result']['sha256']=='dfb777e38e152f7333960ac0989243c2c911af677903a16073e87ffb2546b13c' and
            scout['status']=='four_probe_scout_complete' and scout['selected']['id']=='bounded_LS_full','Wrong start scout')
    oldcfg=read(scout['config']);require(cfg['seed_helper']==oldcfg['implementation'],'Wrong seed helper')
    for key in ('evaluator','evaluator_config','objective','body_helper','water48'):
        require(cfg[key]==oldcfg[key],'Original model helper changed: '+key);verify(cfg[key])
    helper=module(cfg['seed_helper'],'joint_saved_seed_helper')
    parent,ev,obj,body,water,panel,oldbase,unused,seeded=helper.inputs(oldcfg,out)
    additions={}
    for row in scout['rows'].values():
        require(row['status']=='complete','Scout row unavailable')
        for cell in row['model'].values():
            provenance=cell['provenance']
            if not provenance['reused_baseline']:additions[provenance['cache_key']]=provenance['native_authority']
    added=[]
    for key,authority in sorted(additions.items()):
        path=verify(authority);native=read(authority)
        require(native['status']=='complete' and native['cache_key']==key and
                native['cache_identity']['implementation']==cfg['evaluator'],'Scout native authority differs')
        for p in native['input'].values():verify(p)
        for name in ('analyze_output','gradient_output'):verify(native['native'][name])
        require(len(native['receipts'])==2 and all(r['returncode']==0 for r in native['receipts']),'Scout native receipts')
        link=out/'cache/cells'/key
        if link.exists():require(link.resolve()==path.parent,'Conflicting native cache authority')
        else:link.symlink_to(path.parent,target_is_directory=True);added.append({'cache_key':key,'authority':authority,'symlink':str(link)})
    require(len(added)==16 and len(seeded)==142,'Native parent denominators differ')
    ev.root_release=cfg.get('root_release')
    # Replace only the previous instance's scout32 guard, preserving the frozen
    # evaluator class's source/config/release/native admission and1920 ceiling.
    original_gate=type(ev)._gate.__get__(ev,type(ev))
    def bounded_gate():
        original_gate()
        with ev.lock:require(ev._ledger['native_program_invocations']<1152,'Joint refinement native ceiling reached')
    ev._gate=bounded_gate
    initial=scout['selected']['parameters'];model=ev.evaluate(initial);score=obj.evaluate(ev.data,model,body)
    require(score==scout['selected']['score'],'Selected scout score replay differs')
    for gid in CELLS:
        require(model[gid]['energy']==scout['selected']['model'][gid]['energy'] and
                model[gid]['forces']==scout['selected']['model'][gid]['forces'],'Selected scout raw model differs')
    require(ev.snapshot()['native_program_invocations']==0,'Preparation made native calls')
    protected=[cfg['scout_result'],scout['evaluator']['ledger'],scout['parent_result'],parent['evaluator']['ledger']]
    for p in protected:verify(p)
    return ev,obj,body,water,panel,initial,model,score,seeded+added,protected


def directions(J,r,z):
    full,_,rank,singular=np.linalg.lstsq(J,-r,rcond=1e-6)
    lower,upper=np.maximum(-.025,-z),np.minimum(.025,1-z);factor=1.
    for v,lo,hi in zip(full,lower,upper):
        if v>0:factor=min(factor,hi/v)
        elif v<0:factor=min(factor,lo/v)
    scaled=full*factor;bounded=None;loss=None
    for pattern in itertools.product((-1,0,1),repeat=4):
        fixed=[j for j,v in enumerate(pattern) if v];free=[j for j,v in enumerate(pattern) if not v];d=np.zeros(4)
        for j in fixed:d[j]=lower[j] if pattern[j]<0 else upper[j]
        if free:d[free]=np.linalg.lstsq(J[:,free],-(r+J[:,fixed]@d[fixed]),rcond=1e-6)[0]
        if np.any(d<lower-1e-12) or np.any(d>upper+1e-12):continue
        d=np.clip(d,lower,upper);value=float(np.linalg.norm(r+J@d)**2)
        if loss is None or value<loss:bounded,loss=d,value
    require(bounded is not None,'No feasible bounded direction')
    return scaled,bounded,{'rank':int(rank),'singular_values':singular.tolist(),
        'unrestricted_GN_not_evaluated':full.tolist(),'scaled_GN':scaled.tolist(),'bounded_LS':bounded.tolist()}


def main(stage,config_path,output):
    out=Path(output);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    children=resource.getrusage(resource.RUSAGE_CHILDREN);cfg=json.loads(Path(config_path).read_text());ev=None
    result={'status':'unavailable','stage':stage,'implementation':pin(__file__),'config':pin(config_path),
        'scientific_config_sha256':fingerprint(cfg),'new_native_program_invocations':0,'new_parameter_vectors':0,
        'iterations':[],'vectors':{},'physical_accuracy_qualified':False,'validation_outcomes_opened':False,
        'automatic_extension':False,'global_optimum_claimed':False,'job_id':os.environ.get('SLURM_JOB_ID'),'error':None}
    try:
        if stage=='run':
            require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and
                    int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==8,'Exact eight-CPU allocation required')
            prep=read(cfg['preparation_result']);release=read(cfg['root_release'])
            require(prep['status']=='prepared' and prep['new_native_program_invocations']==0 and
                    prep['scientific_config_sha256']==fingerprint(cfg),'Actual preparation differs')
            require(release['joint_status']=='approved_bounded_joint_water_refinement' and
                    release['joint_max_new_vectors']==144 and release['joint_max_native_program_invocations']==1152 and
                    release['joint_max_outer_steps']==12 and release['joint_workers']==8 and
                    release['joint_config_sha256']==fingerprint(cfg),'Joint root scope differs')
            for key in ('implementation','specification','scout_result','preparation_result'):
                require(release['joint_'+key]==cfg[key],'Joint release binding differs: '+key)
        else:require(cfg.get('root_release') is None and cfg.get('preparation_result') is None,'Preparation must be unreleased')
        ev,obj,body,water,panel,initial,model,score,seeded,protected=prepare(cfg,out)
        parameters=ev.family['parameters'];lower=np.asarray([p['lower'] for p in parameters]);upper=np.asarray([p['upper'] for p in parameters]);width=upper-lower
        current={'parameters':initial,'normalized_parameters':((np.asarray(initial)-lower)/width).tolist(),
                 'status':'complete','role':'initial_scout','model':model,'score':score,'cells':{g:'complete' for g in CELLS}}
        records=result['vectors'];records[token(initial)]=current
        result.update(initial=current,seeded_cache=seeded,protected_parents=protected,evaluator=ev.snapshot(),family=ev.family)
        def persist():
            result.update(evaluator=ev.snapshot(),new_native_program_invocations=ev.snapshot()['native_program_invocations'],new_parameter_vectors=len(records)-1)
            save(out/'RESULT.json',result)
        def batch(points,role,iteration):
            keys=[token(p['parameters']) for p in points]
            novel={k for k in keys if k not in records}
            require(len(records)-1+len(novel)<=144,'New-vector ceiling reached')
            todo=[]
            for point,key in zip(points,keys):
                if key in records:
                    require(records[key]['status']=='complete' or key in todo,'Cannot repeat unavailable point');continue
                records[key]={'parameters':point['parameters'],'normalized_parameters':((np.asarray(point['parameters'])-lower)/width).tolist(),
                    'status':'scheduled','role':role,'outer_step':iteration,'cells':{g:'not_attempted' for g in CELLS}}
                todo.append(key)
            persist();errors=[]
            def call(key):
                require(not ev.stop_event.is_set(),'Global native stop active')
                values=records[key]['parameters'];m=ev.evaluate(values)
                return {'model':m,'score':obj.evaluate(ev.data,m,body),'cells':{g:'complete' for g in CELLS}}
            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
                futures={pool.submit(call,k):k for k in todo}
                for future in concurrent.futures.as_completed(futures):
                    key=futures[future]
                    try:records[key].update(status='complete',**future.result())
                    except Exception as exc:
                        ev.stop_event.set();records[key].update(status='unavailable',error=repr(exc));errors.append(key)
                        for other in futures:
                            if other is not future:other.cancel()
                    persist()
            if errors:
                requests=ev.snapshot()['requests']
                for key in errors:
                    matched=[r for r in requests if r['vector']==records[key]['parameters']]
                    if matched:records[key]['cells']=matched[-1]['cells']
                persist();raise RuntimeError('Native/vector failure; no next batch')
            return [records[k] for k in keys]
        if stage=='prepare':result['status']='prepared'
        else:
            require(initial==prep['initial']['parameters'],'Prepared starting vector differs')
            stop='outer_step_limit'
            for iteration in range(1,13):
                if current['score']['development_checks_passed']:stop='all_existing_controls_passed';break
                center=np.asarray(current['parameters']);z=(center-lower)/width
                entry={'outer_step':iteration,'origin_key':token(current['parameters']),'status':'derivatives',
                       'derivative_slots':[f'{j}_{s}' for j in range(4) for s in ('minus','plus')],'trial_slots':[]}
                result['iterations'].append(entry);points=[]
                for j in range(4):
                    for sign in (-1,1):
                        target=z.copy();target[j]=np.clip(target[j]+sign*.005,0.,1.)
                        values=center.copy();values[j]=np.clip(center[j]+(target[j]-z[j])*width[j],lower[j],upper[j])
                        points.append({'parameters':values.tolist()})
                entry['derivative_points']=points
                derivatives=batch(points,'derivative',iteration);columns=[]
                for j in range(4):
                    a,b=derivatives[2*j:2*j+2];span=b['normalized_parameters'][j]-a['normalized_parameters'][j]
                    require(span>0,'Collapsed derivative interval')
                    columns.append((np.asarray(b['score']['weighted_residuals'])-np.asarray(a['score']['weighted_residuals']))/span)
                J=np.column_stack(columns);r=np.asarray(current['score']['weighted_residuals'])
                scaled,bounded,diagnostic=directions(J,r,z)
                entry.update(status='trials',derivative_keys=[token(v['parameters']) for v in derivatives],
                    derivative_actual_widths=[derivatives[2*j+1]['normalized_parameters'][j]-derivatives[2*j]['normalized_parameters'][j] for j in range(4)],
                    Jacobian_normalized=J.tolist(),linear_directions=diagnostic)
                trials=[]
                for name,d in (('scaled_GN',scaled),('bounded_LS',bounded)):
                    for fraction,suffix in ((.5,'half'),(1.,'full')):
                        values=np.clip(center+fraction*d*width,lower,upper);actual=(values-center)/width
                        require(np.max(np.abs(actual))<=.025+1e-12,'Trial exceeds fixed trust box')
                        trials.append({'id':name+'_'+suffix,'parameters':values.tolist(),
                            'normalized_delta':actual.tolist(),'predicted_objective':float(np.linalg.norm(r+J@actual)**2)})
                entry['trial_slots']=trials;completed=batch(trials,'trial',iteration);selected=current
                for value in completed:
                    if value['score']['objective']<selected['score']['objective']:selected=value
                improved=selected is not current;current=selected
                entry.update(status='complete',trial_keys=[token(v['parameters']) for v in completed],
                    selected_key=token(current['parameters']),actual_improvement=improved,
                    all_controls_passed=current['score']['development_checks_passed'])
                result['current_key']=token(current['parameters']);persist()
                if current['score']['development_checks_passed']:stop='all_existing_controls_passed';break
                if not improved:stop='no_actual_improvement_at_complete_trial_batch';break
            result.update(status='bounded_joint_refinement_complete',stopping_reason=stop,selected=current,
                development_checks_passed=current['score']['development_checks_passed'],
                selected_method_sensitivity=obj.method_sensitivity(ev.data,current['model'],water,panel),
                candidate_frozen=True,automatic_validation_release=False,completed_outer_steps=len(result['iterations']))
        for p in protected:verify(p)
        result['original_results_and_ledgers_unchanged']=True
    except Exception as exc:result.update(status='failed_or_unavailable',error=repr(exc),candidate_frozen=False)
    finally:
        if ev is not None:result['evaluator']=ev.snapshot();result['new_native_program_invocations']=result['evaluator']['native_program_invocations']
        end=resource.getrusage(resource.RUSAGE_CHILDREN)
        result.update(elapsed_seconds=time.monotonic()-start,parent_process_CPU_seconds=time.process_time()-cpu,
            aggregate_native_child_user_CPU_seconds=end.ru_utime-children.ru_utime,aggregate_native_child_system_CPU_seconds=end.ru_stime-children.ru_stime,
            numpy_version=np.__version__)
        save(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'new_native_program_invocations':result['new_native_program_invocations'],'result':pin(out/'RESULT.json')}))
    return 0 if result['status'] in ('prepared','bounded_joint_refinement_complete') else 1


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('prepare','run'));p.add_argument('--config',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();raise SystemExit(main(a.stage,a.config,a.output))
