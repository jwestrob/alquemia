#!/usr/bin/env python3
"""Six-parameter paired water development from exact17593 selection."""
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
CELLS=tuple(q+'_'+m for q in ('water_minus','water_central','water_plus','aquo8_central') for m in ('La','Dy'))


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
    for key in ('specification','evaluator','evaluator_config','objective','body_helper','water48','center_objective','joint_math','training_result','prior_result','epsilon_result','root_declaration'):verify(cfg[key])
    data=read(cfg['training_result']);require(data['status']=='radial_training_complete' and data['complete']==8 and data['saved_controls_passed'],'Data unavailable')
    require(cfg['training_result']['sha256']=='3aefeb74f40aaf51e327fdc5c0ffd85f7f00e6e769b07c080f3ece16c345192d','Wrong consumed data')
    evaluator=module(cfg['evaluator'],'radial_native_evaluator');ecfg=read(cfg['evaluator_config'])
    require(ecfg['training_result']==cfg['training_result'] and ecfg['implementation']==cfg['evaluator'],'Evaluator/data join differs')
    ev=evaluator.Evaluator(ecfg,out/'cache',root_release=cfg.get('root_release'))
    obj=module(cfg['objective'],'radial_objective');body=module(cfg['body_helper'],'radial_body')
    center=module(cfg['center_objective'],'original_four_center_objective');water=read(cfg['water48']);panel=read(water['panel'])
    prior=read(cfg['prior_result'])
    require(cfg['prior_result']['sha256']=='fa2ef37ca71e79b27fdd8c27f8e7b694682ec2416e6ba477c7e64b77c9e79e09' and ecfg['prior_result']==cfg['prior_result'],'Exact continuation parent differs')
    require(prior['implementation']['sha256']=='f6550a635a6f4fd0ee3e801011bb99c869d4ca7352f209c656fee52d59ca96cd' and prior['status']=='bounded_radial_continuation_complete' and prior['error'] is None,'Prior result/source differs');verify(prior['implementation'])
    require(prior['stopping_reason']=='no_actual_improvement_at_complete_trial_batch','Prior local stop differs')
    epsilon=read(cfg['epsilon_result']);require(epsilon['status']=='epsilon_diagnostic_complete' and epsilon['error'] is None and epsilon['job_id']=='17596' and epsilon['native_program_invocations']==32 and epsilon['analysis']['parameter_linearity_passed'] and epsilon['analysis']['retained_increment_rank']==2,'Epsilon diagnostic unavailable')
    require(epsilon['config']['candidate_result']==cfg['prior_result'] and ecfg['epsilon_result']==cfg['epsilon_result'],'Epsilon lineage differs')
    initial=prior['selected']['parameters']+[.24429636467571708,.37332330078584813];model=ev.evaluate(initial);score=obj.evaluate(ev.data,model,body)
    require(initial==ev.continuation_start and score==prior['selected']['score'],'Exact17593 selected score differs')
    for gid in CELLS:
        require(model[gid]['energy']==prior['selected']['model'][gid]['energy'] and model[gid]['forces']==prior['selected']['model'][gid]['forces'],'Initial raw model differs')
    require(ev.snapshot()['native_program_invocations']==0,'Initial admission made calls')
    require(center_diagnostics(center,ev.data,model,body,water,panel)==prior['selected_center_diagnostics'],'Prior selected center diagnostics differ')
    protected=[data['transfer_result'],data['transfer_bound'],data['candidate_result'],cfg['training_result'],cfg['prior_result'],prior['config'],prior['evaluator']['ledger'],cfg['epsilon_result'],epsilon['manifest']]
    for p in protected:verify(p)
    return ev,obj,body,center,water,panel,initial,model,score,protected


def center_diagnostics(center,data,model,body,water,panel):
    ids=center.CELLS;subset={k:model[k] for k in ids};olddata=dict(data,rows={k:dict(data['rows'][k],baseline=data['rows'][k]['original_baseline']) for k in ids})
    return {'original_four_center_score':center.evaluate(olddata,subset,body),'method_sensitivity':center.method_sensitivity(olddata,subset,water,panel)}



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
            require(release['six_status']=='approved_bounded_six_parameter_water' and
                    release['six_max_new_vectors']==384 and release['six_max_native_program_invocations']==4096 and
                    release['six_max_outer_steps']==24 and release['six_workers']==8 and
                    release['six_config_sha256']==fingerprint(cfg),'Radial root scope differs')
            for key in ('implementation','specification','training_result','preparation_result','prior_result','epsilon_result','root_declaration'):
                require(release['six_'+key]==cfg[key],'Radial release binding differs: '+key)
        else:require(cfg.get('root_release') is None and cfg.get('preparation_result') is None,'Preparation must be unreleased')
        ev,obj,body,center_helper,water,panel,initial,model,score,protected=prepare(cfg,out)
        math_helper=module(cfg['joint_math'],'frozen_joint_direction_math')
        require(cfg['joint_math']['sha256']=='741167a8c5ce3c9e07af6a94baf40894c277b5c47950f31b3c526da125ae11ca','Direction mathematics changed')
        parameters=ev.family['parameters'];require(len(parameters)==6,'Exactly six parameters required');lower=np.asarray([p['lower'] for p in parameters]);upper=np.asarray([p['upper'] for p in parameters]);width=upper-lower
        current={'parameters':initial,'normalized_parameters':((np.asarray(initial)-lower)/width).tolist(),
                 'status':'complete','role':'initial_selected17593_original_epsilons','model':model,'score':score,'cells':{g:'complete' for g in CELLS}}
        records=result['vectors'];records[token(initial)]=current
        result.update(initial=current,prior_result=cfg['prior_result'],protected_parents=protected,evaluator=ev.snapshot(),family=ev.family,initial_center_diagnostics=center_diagnostics(center_helper,ev.data,model,body,water,panel))
        def persist():
            result.update(evaluator=ev.snapshot(),new_native_program_invocations=ev.snapshot()['native_program_invocations'],new_parameter_vectors=len(records)-1)
            save(out/'RESULT.json',result)
        def batch(points,role,iteration):
            keys=[token(p['parameters']) for p in points]
            novel={k for k in keys if k not in records}
            require(len(records)-1+len(novel)<=384,'New-vector ceiling reached')
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
            for iteration in range(1,25):
                if current['score']['development_checks_passed']:stop='all_existing_controls_passed';break
                center=np.asarray(current['parameters']);z=(center-lower)/width
                entry={'outer_step':iteration,'origin_key':token(current['parameters']),'status':'derivatives',
                       'derivative_slots':[f'{j}_{s}' for j in range(6) for s in ('minus','plus')],'trial_slots':[]}
                result['iterations'].append(entry);points=[]
                for j in range(6):
                    for sign in (-1,1):
                        target=z.copy();target[j]=np.clip(target[j]+sign*.005,0.,1.)
                        values=center.copy();values[j]=np.clip(center[j]+(target[j]-z[j])*width[j],lower[j],upper[j])
                        points.append({'parameters':values.tolist()})
                entry['derivative_points']=points
                derivatives=batch(points,'derivative',iteration);columns=[]
                for j in range(6):
                    a,b=derivatives[2*j:2*j+2];span=b['normalized_parameters'][j]-a['normalized_parameters'][j]
                    require(span>0,'Collapsed derivative interval')
                    columns.append((np.asarray(b['score']['weighted_residuals'])-np.asarray(a['score']['weighted_residuals']))/span)
                J=np.column_stack(columns);r=np.asarray(current['score']['weighted_residuals'])
                scaled,bounded,diagnostic=math_helper.directions(J,r,z)
                entry.update(status='trials',derivative_keys=[token(v['parameters']) for v in derivatives],
                    derivative_actual_widths=[derivatives[2*j+1]['normalized_parameters'][j]-derivatives[2*j]['normalized_parameters'][j] for j in range(6)],
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
            result.update(status='bounded_six_parameter_water_complete',stopping_reason=stop,selected=current,
                development_checks_passed=current['score']['development_checks_passed'],
                selected_center_diagnostics=center_diagnostics(center_helper,ev.data,current['model'],body,water,panel),
                candidate_frozen=True,automatic_validation_release=False,completed_outer_steps=len(result['iterations']))
        for p in protected:verify(p)
        result['original_results_unchanged']=True
    except Exception as exc:result.update(status='failed_or_unavailable',error=repr(exc),candidate_frozen=False)
    finally:
        if ev is not None:result['evaluator']=ev.snapshot();result['new_native_program_invocations']=result['evaluator']['native_program_invocations']
        end=resource.getrusage(resource.RUSAGE_CHILDREN)
        result.update(elapsed_seconds=time.monotonic()-start,parent_process_CPU_seconds=time.process_time()-cpu,
            aggregate_native_child_user_CPU_seconds=end.ru_utime-children.ru_utime,aggregate_native_child_system_CPU_seconds=end.ru_stime-children.ru_stime,
            numpy_version=np.__version__)
        save(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'new_native_program_invocations':result['new_native_program_invocations'],'result':pin(out/'RESULT.json')}))
    return 0 if result['status'] in ('prepared','bounded_six_parameter_water_complete') else 1


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('prepare','run'));p.add_argument('--config',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();raise SystemExit(main(a.stage,a.config,a.output))
