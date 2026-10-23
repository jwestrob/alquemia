#!/usr/bin/env python3
"""Four epsilon probes at one frozen point; no fit or selected new candidate."""
import argparse
import concurrent.futures
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import threading
import time
import numpy as np

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
EPS={'La':.24429636467571708,'Dy':.37332330078584813}
POINT=[3.832755583060068,.31164088100092946,3.525530001514476,.32135326991694113]
CELLS=tuple(q+'_'+m for q in ('water_minus','water_central','water_plus','aquo8_central')for m in ('La','Dy'))


def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed evidence: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def token(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def fingerprint(c):return token({k:v for k,v in c.items()if k not in ('root_release','preparation_result')})
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def save(p,x):
    p=Path(p);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)


def epsilon_key(base,prm,metal,epsilon,resolver):
    mt=500 if metal=='La'else 501;before=resolver.parse_parameters(prm,base);pair=resolver.resolve_pair(before,mt,36)
    require(pair['epsilon_kcal_mol']==EPS[metal],'Frozen epsilon differs')
    lines=base.splitlines(keepends=True);changed=[]
    for j,line in enumerate(lines):
        f=line.split()
        if f and f[0].lower()=='vdwpr'and sorted(map(int,f[1:3]))==[34,mt]:
            require(len(f)==5,'Unexpected water pair line');f[4]=format(epsilon,'.17g');lines[j]=' '.join(f)+'\n';changed.append(j)
    require(len(changed)==1,'Exactly one explicit water pair required');key=''.join(lines);after=resolver.parse_parameters(prm,key)
    for k in before:
        if k=='vdwpr':require({x:v for x,v in before[k].items()if x!=(34,mt)}=={x:v for x,v in after[k].items()if x!=(34,mt)},'Other vdW pair changed')
        else:require(before[k]==after[k],'Non-epsilon parameter changed: '+k)
    final=resolver.resolve_pair(after,mt,36)
    require(final['Rmin_pair_A']==pair['Rmin_pair_A']and final['thole']==pair['thole']and final['epsilon_kcal_mol']==epsilon,'Rendered epsilon/Rmin/damping differs')
    require(resolver.resolve_pair(before,mt,37)==resolver.resolve_pair(after,mt,37),'Water-H pair changed')
    return key,{'metal':metal,'baseline_pair':pair,'probe_pair':final,'changed_line_index':changed[0],'only_water_pair_epsilon_changed':True}


def prepare(cfg,out):
    require(cfg['implementation']==pin(__file__)and all(os.environ.get(k)=='1'for k in THREADS),'Source/thread controls differ')
    for k,p in cfg.items():
        if isinstance(p,dict)and set(p)=={'path','sha256'}:verify(p)
    require(cfg['candidate_result']['sha256']=='fa2ef37ca71e79b27fdd8c27f8e7b694682ec2416e6ba477c7e64b77c9e79e09','Wrong frozen17593 result')
    r=read(cfg['candidate_result']);require(r['status']=='bounded_radial_continuation_complete'and r['error']is None and r['job_id']=='17593'and r['stopping_reason']=='no_actual_improvement_at_complete_trial_batch','Wrong terminal point')
    selected=r['selected'];require(selected['parameters']==POINT and selected['status']=='complete','Frozen parameters differ')
    c=read(r['config']);ecfg=read(c['evaluator_config'])
    require(c['implementation']==cfg['candidate_implementation']and c['evaluator']==cfg['candidate_evaluator']and c['objective']==cfg['objective']and c['training_result']==cfg['training']and c['body_helper']==cfg['body_helper'],'Consumed source joins differ')
    adapter=module(cfg['authority_helper'],'epsilon_saved_authority');data=adapter.loaded(verify(cfg['original32']));training=read(cfg['training'])
    body=module(cfg['body_helper'],'epsilon_body');objective=module(cfg['objective'],'epsilon_objective');native=module(data['config']['native_helper'],'epsilon_native');resolver=module(data['config']['resolver'],'epsilon_resolver')
    prior=read(adapter.PINS['pins']['continuation_prior_result']);previous=read(adapter.PINS['pins']['previous_transfer_result'])
    admitted={}
    def admit(vector):
        require(vector['status']=='complete'and set(vector['model'])==set(CELLS),'Incomplete derivative/center vector')
        require(objective.evaluate(training,vector['model'],body)==vector['score'],'Saved score replay differs')
        for gid,value in vector['model'].items():
            row=data['rows'][gid];pv=value['provenance'];inp=pv['effective_input'];key,_=adapter.d.render(data['config'],row,vector['parameters'])
            require(verify(inp['key']).read_text()==key and inp['xyz']['sha256']==row['input']['xyz']['sha256']and inp['parameters']['sha256']==row['input']['parameters']['sha256'],'Saved physical input differs')
            keyid=(gid,pv['native_authority']['path'],pv['native_authority']['sha256'])
            if keyid not in admitted:
                branch=adapter.selected_authority(gid,value,row,inp,ecfg,cfg['candidate_result']['path'],previous,prior)
                adapter.replay_native(value,row,native);admitted[keyid]={'cell':gid,'authority':pv['native_authority'],'branch':branch}
    admit(selected);last=r['iterations'][-1];require(last['status']=='complete'and last['origin_key']==token(POINT)==last['selected_key'],'Last Jacobian is not at selected point')
    require(len(last['derivative_keys'])==8,'Exact eight derivative slots required');derivatives=[r['vectors'][k]for k in last['derivative_keys']]
    lower=np.asarray([p['lower']for p in r['family']['parameters']]);upper=np.asarray([p['upper']for p in r['family']['parameters']]);width=upper-lower;z=(np.asarray(POINT)-lower)/width;columns=[];spans=[]
    for j in range(4):
        for offset,sign in enumerate((-1,1)):
            value=derivatives[2*j+offset];target=z.copy();target[j]=np.clip(target[j]+sign*.005,0.,1.);expected=np.asarray(POINT).copy();expected[j]=np.clip(expected[j]+(target[j]-z[j])*width[j],lower[j],upper[j])
            require(value['parameters']==expected.tolist(),'Original FD point/step differs');admit(value)
        a,b=derivatives[2*j:2*j+2];span=b['normalized_parameters'][j]-a['normalized_parameters'][j];spans.append(span)
        columns.append((np.asarray(b['score']['weighted_residuals'])-np.asarray(a['score']['weighted_residuals']))/span)
    J4=np.column_stack(columns);require(J4.tolist()==last['Jacobian_normalized']and spans==last['derivative_actual_widths'],'Original Jacobian replay differs')
    probes=[]
    for metal in ('La','Dy'):
        for sign,label in ((-1,'minus'),(1,'plus')):
            epsilon=EPS[metal]*(1+sign*.05);entry={'id':metal+'_'+label,'metal':metal,'epsilon':epsilon,'eta':(epsilon-EPS[metal])/EPS[metal],'inputs':{},'reused_cells':{}}
            for gid in CELLS:
                if not gid.endswith('_'+metal):entry['reused_cells'][gid]=selected['model'][gid];continue
                value=selected['model'][gid];inp=value['provenance']['effective_input'];key,render=epsilon_key(verify(inp['key']).read_text(),verify(inp['parameters']).read_text(),metal,epsilon,resolver)
                directory=out/'inputs'/entry['id']/gid;directory.mkdir(parents=True)
                shutil.copyfile(verify(inp['xyz']),directory/'system.xyz');shutil.copyfile(verify(inp['parameters']),directory/'amoeba09.prm');(directory/'system.key').write_text(key)
                entry['inputs'][gid]={'files':{k:pin(directory/n)for k,n in(('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))},'render':render,'baseline_input':inp}
            require(len(entry['inputs'])==len(entry['reused_cells'])==4,'Opposite-metal reuse differs');probes.append(entry)
    manifest={'status':'prepared','implementation':pin(__file__),'scientific_config_sha256':fingerprint(cfg),'candidate_result':cfg['candidate_result'],'baseline':selected,'probes':probes,
        'four_column_Jacobian':J4.tolist(),'four_column_actual_spans':spans,'derivative_keys':last['derivative_keys'],'admitted_native_authorities':list(admitted.values()),
        'parameter_coordinate_convention':{'original_four':'(value-lower)/(upper-lower)','epsilon':'(epsilon-epsilon0)/epsilon0'},
        'maximum_new_native_programs':32,'maximum_new_physical_cells':16,'workers':4,'new_native_calls':0,
        'protected_parents':[cfg['candidate_result'],r['config'],r['evaluator']['ledger'],cfg['training'],cfg['original32'],adapter.PINS['pins']['continuation_prior_result']]}
    for p in manifest['protected_parents']:verify(p)
    save(out/'MANIFEST.json',manifest);return manifest


def analyze(manifest,records,data,training,body,objective):
    baseline=manifest['baseline'];r=np.asarray(baseline['score']['weighted_residuals']);J4=np.asarray(manifest['four_column_Jacobian']);Je=[];physical=[];linearity={};responses={}
    for metal in ('La','Dy'):
        minus,plus=(records[metal+'_'+s]for s in('minus','plus'));span=plus['eta']-minus['eta'];epsspan=plus['epsilon']-minus['epsilon']
        delta=np.asarray(plus['score']['weighted_residuals'])-np.asarray(minus['score']['weighted_residuals']);Je.append(delta/span);physical.append(delta/epsspan)
        for gid in CELLS:
            if not gid.endswith('_'+metal):continue
            a,b,c=minus['model'][gid],plus['model'][gid],baseline['model'][gid]
            energy=(a['energy']+b['energy'])/2-c['energy'];force=max(abs((x+y)/2-z)for av,bv,cv in zip(a['forces'],b['forces'],c['forces'])for x,y,z in zip(av,bv,cv))
            components={k:{'minus_minus_center':a['native']['components'][k]['energy_kcal_mol']-v['energy_kcal_mol'],'plus_minus_center':b['native']['components'][k]['energy_kcal_mol']-v['energy_kcal_mol']}for k,v in c['native']['components'].items()}
            linearity[gid]={'midpoint_energy_error_kcal_mol':energy,'maximum_midpoint_force_component_error_kcal_mol_A':force,'passed':abs(energy)<=.0001 and force<=.001,'native_energy_components':components}
        energies={k:{q:(plus['score']['energy_errors_kcal_mol'][k][q]-v)/span for q,v in values.items()}for k,values in minus['score']['energy_errors_kcal_mol'].items()}
        rigid={};slopes={}
        for gid in CELLS:
            row=training['rows'][gid];a=body.response(row,minus['model'][gid]['forces']);b=body.response(row,plus['model'][gid]['forces'])
            rigid[gid]=[{'category':x['category'],'indices':x['indices'],**{k:[(yy-xx)/span for xx,yy in zip(x[k],y[k])]for k in('net_force','torque_over_1A')}}for x,y in zip(a,b)]
            if gid.startswith('water_'):
                m=gid.rsplit('_',1)[1];path=data['paths']['water_'+m];ids=path['moving_atoms'];u=path['direction']
                slopes[gid]=(-sum(plus['model'][gid]['forces'][i][j]*u[j]for i in ids for j in range(3))+sum(minus['model'][gid]['forces'][i][j]*u[j]for i in ids for j in range(3)))/span
        for q in('water_minus','water_central','water_plus'):slopes[q+'_Dy_minus_La']=slopes[q+'_Dy']-slopes[q+'_La']
        responses[metal]={'eta_span':span,'epsilon_span_kcal_mol':epsspan,'energy_error_derivatives_per_eta':energies,'rigid_body_derivatives_per_eta':rigid,'outward_coordinate_slope_derivatives_per_eta':slopes}
    Je=np.column_stack(Je);U,s4,_=np.linalg.svd(J4,full_matrices=False);threshold=float(s4[0]*1e-6);B=U[:,s4>threshold];old_projection=B@(B.T@r);W=Je-B@(B.T@Je)
    Ue,se,_=np.linalg.svd(W,full_matrices=False);C=Ue[:,se>threshold];require(np.max(np.abs(B.T@C),initial=0.)<1e-10,'Increment is not orthogonal');extended_projection=old_projection+C@(C.T@r)
    slices={'energy':(0,7)};index=7
    for kind in('La','Dy','paired'):
        begin=index
        for q in('water_minus','water_central','water_plus','aquo8_central'):
            row=training['rows'][q+'_La'];n=len(body.rigid_residual(row,row['reference']['forces'],row['reference']['forces'])[0]);slices[kind+':'+q]=(index,index+n);index+=n
        slices[kind]=(begin,index)
    require(index==len(r),'Objective block/fixture layout differs')
    remaining={'original':r,'after_original_four_span':r-old_projection,'after_nested_six_span':r-extended_projection}
    return {'parameter_linearity_passed':all(x['passed']for x in linearity.values()),'linearity':linearity,'epsilon_columns_per_eta':Je.tolist(),'epsilon_columns_per_kcal_mol':np.column_stack(physical).tolist(),
        'signed_epsilon_responses':responses,'original_singular_values':s4.tolist(),'absolute_singular_threshold':threshold,'orthogonal_increment_singular_values':se.tolist(),
        'retained_original_rank':B.shape[1],'retained_increment_rank':C.shape[1],'augmented_singular_values_descriptive_only':np.linalg.svd(np.column_stack((J4,Je)),compute_uv=False).tolist(),
        'column_norms':{m:{'epsilon_norm':float(np.linalg.norm(Je[:,j])),'outside_original_span_norm':float(np.linalg.norm(W[:,j])),'outside_fraction':float(np.linalg.norm(W[:,j])/np.linalg.norm(Je[:,j]))if np.linalg.norm(Je[:,j])else None}for j,m in enumerate(('La','Dy'))},
        'weighted_residual_losses':{k:float(v@v)for k,v in remaining.items()},'residual_components':{k:v.tolist()for k,v in remaining.items()},
        'weighted_residual_decomposition':{k:{label:float(v[a:b]@v[a:b])for label,(a,b)in slices.items()}for k,v in remaining.items()},
        'parameter_solution_computed':False,'new_candidate_selected':False,'physical_accuracy_qualified':False}


def run(cfg,out):
    start=time.monotonic();cpu=time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN)
    require(os.environ.get('SLURM_JOB_ID','').isdigit()and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==4,'Exact four-CPU Slurm task required')
    require(cfg['implementation']==pin(__file__)and all(os.environ.get(k)=='1'for k in THREADS),'Source/thread controls differ')
    prep=read(cfg['preparation_result']);manifest=read(prep['manifest']);release=read(cfg['root_release'])
    require(prep['status']=='prepared'and prep['new_native_calls']==0 and manifest['scientific_config_sha256']==fingerprint(cfg),'Actual preparation differs')
    require(release['status']=='approved_frozen_epsilon_diagnostic'and release['implementation']==cfg['implementation']and release['preparation_result']==cfg['preparation_result']and release['scientific_config_sha256']==fingerprint(cfg)and release['maximum_new_native_programs']==32 and release['workers']==4,'Root finite release differs')
    adapter=module(cfg['authority_helper'],'epsilon_run_authority');data=adapter.loaded(verify(cfg['original32']));training=read(cfg['training']);body=module(cfg['body_helper'],'epsilon_run_body');objective=module(cfg['objective'],'epsilon_run_objective');native=module(data['config']['native_helper'],'epsilon_run_native');resolver=module(data['config']['resolver'],'epsilon_run_resolver')
    for p in manifest['protected_parents']:verify(p)
    result={'status':'running','implementation':pin(__file__),'config':cfg,'preparation_result':cfg['preparation_result'],'manifest':prep['manifest'],'job_id':os.environ['SLURM_JOB_ID'],
        'native_program_invocations':0,'records':{p['id']:{'status':'not_attempted','metal':p['metal'],'epsilon':p['epsilon'],'eta':p['eta'],'cells':{g:{'status':'not_attempted'}for g in CELLS}}for p in manifest['probes']},'error':None,'new_candidate_selected':False,'reserved_validation_opened':False}
    lock=threading.RLock();stop=threading.Event()
    def persist():
        with lock:save(out/'RESULT.json',result)
    def worker(probe):
        record=result['records'][probe['id']]
        with lock:
            record['status']='started';record['model']=copy.deepcopy(probe['reused_cells'])
            for gid,value in probe['reused_cells'].items():record['cells'][gid]={'status':'complete','reused':True,'authority':value['provenance']['native_authority']}
        try:
            for gid,inp in probe['inputs'].items():
                require(not stop.is_set(),'Global native stop active');row=training['rows'][gid];directory=out/probe['id']/gid;directory.mkdir(parents=True);cell={'status':'started','input':inp,'receipts':[]}
                with lock:record['cells'][gid]=cell
                key,render=epsilon_key(verify(inp['baseline_input']['key']).read_text(),verify(inp['baseline_input']['parameters']).read_text(),probe['metal'],probe['epsilon'],resolver)
                require(key==verify(inp['files']['key']).read_text()and render==inp['render'],'Prepared epsilon key changed')
                for k,n in(('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(inp['files'][k]),directory/n)
                save(directory/'RESULT.json',cell)
                try:
                    for program,args in(('analyze',['system.xyz','E']),('testgrad',['system.xyz','Y','N'])):
                        with lock:
                            require(not stop.is_set()and result['native_program_invocations']<32,'Native stop/call ceiling');result['native_program_invocations']+=1;save(out/'CALL_COUNT.json',{'native_program_invocations':result['native_program_invocations']})
                        receipts=[]
                        try:
                            native.invoke([data['config']['executables'][program]['path']]+args,directory,directory/(program+'.out'),receipts)
                            adapter.program_admission(directory/(program+'.out'),program,native)
                        except Exception:
                            stop.set();raise
                        finally:
                            with lock:cell['receipts'].extend(receipts);save(directory/'RESULT.json',cell)
                    parsed=native.parse(directory/'analyze.out',directory/'testgrad.out',dict(row['input'],derivative_test=False));value={'energy':parsed['total_energy_kcal_mol'],'forces':parsed['force_kcal_mol_A'],'native':parsed}
                    with lock:cell.update(status='complete',model=value);save(directory/'RESULT.json',cell);record['model'][gid]=value
                except Exception as exc:
                    stop.set()
                    with lock:cell.update(status='failed',error=repr(exc));save(directory/'RESULT.json',cell)
                    raise
            score=objective.evaluate(training,record['model'],body)
            with lock:record['score']=score;record['status']='complete'
        except Exception as exc:
            stop.set()
            with lock:record.update(status='failed_or_unavailable',error=repr(exc));result['error']=repr(exc);persist()
    try:
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=4)as pool:list(pool.map(worker,manifest['probes']))
        require(not stop.is_set()and all(r['status']=='complete'for r in result['records'].values()),'Native probe denominator unavailable')
        result['analysis']=analyze(manifest,result['records'],data,training,body,objective)
        result['status']='epsilon_diagnostic_complete'if result['analysis']['parameter_linearity_passed']else'parameter_linearity_failed'
        for p in manifest['protected_parents']:verify(p)
        result['protected_parent_pins_unchanged']=True
    except Exception as exc:result.update(status='failed_or_unavailable',error=repr(exc))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result.update(elapsed_seconds=time.monotonic()-start,parent_cpu_seconds=time.process_time()-cpu,aggregate_child_user_cpu_seconds=ch.ru_utime-ch0.ru_utime,aggregate_child_system_cpu_seconds=ch.ru_stime-ch0.ru_stime);persist()
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('prepare','run'));p.add_argument('--config',required=True);p.add_argument('--output',required=True);a=p.parse_args();cfg=json.loads(Path(a.config).read_text());out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    try:
        if a.stage=='prepare':
            require(cfg['root_release']is None and cfg['preparation_result']is None,'Preparation must be unreleased');m=prepare(cfg,out);r={'status':'prepared','implementation':pin(__file__),'config':pin(a.config),'manifest':pin(out/'MANIFEST.json'),'new_native_calls':0,'elapsed_seconds':time.monotonic()-start,'parent_cpu_seconds':time.process_time()-cpu};save(out/'RESULT.json',r)
        else:r=run(cfg,out)
    except Exception as exc:r={'status':'precall_or_preparation_failed','error':repr(exc),'implementation':pin(__file__)};save(out/'FAILURE.json',r)
    print(json.dumps({'status':r['status'],'output':str(out)}));raise SystemExit(0 if r['status']in('prepared','epsilon_diagnostic_complete')else 1)
