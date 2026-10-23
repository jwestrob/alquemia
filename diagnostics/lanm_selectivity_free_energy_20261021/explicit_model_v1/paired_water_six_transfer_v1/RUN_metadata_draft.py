#!/usr/bin/env python3
"""Conditional passing six-parameter candidate transfer; unchanged native and score kernels."""
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

HERE=Path(__file__).resolve().parent
PINS=json.loads((HERE/'PINS.json').read_text())
def require(ok,message):
    if not ok:raise ValueError(message)
require(hashlib.sha256((HERE/'PINS.json').read_bytes()).hexdigest()=='264887b029b815ed9a6cfcfe9b992dbc2dacb75ac042b9a44f338c081715d30d','Plan pin inventory changed')
def initial_module(record,name):
    require(hashlib.sha256(Path(record['path']).read_bytes()).hexdigest()==record['sha256'],'Changed source')
    spec=importlib.util.spec_from_file_location(name,record['path']);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
old=initial_module(PINS['pins']['previous_run_source'],'unchanged_consumed_transfer')
d=old.d
pin,verify,read,require,module,save=old.pin,old.verify,old.read,old.require,old.module,old.save
loaded,gradient,program_admission,replay_native=old.loaded,old.gradient,old.program_admission,old.replay_native
renderer=module(PINS['pins']['renderer'],'six_water_transfer_renderer')
def render(cfg,row,values):return renderer.render(d,cfg,row,values)
WATER=d.CENTERS+d.RADIAL
MEMBERSHIPS={k:PINS['inventory'][k]for k in ('original_center_development_4','radial_water_development_8','hydrated_donor_transfer_18','unchanged_dry_controls_6')}


def family_gate(family,training,ecfg):
    evaluator=module(PINS['pins']['six_evaluator'],'six_transfer_family_constants')
    require(family['original_pair_epsilon']==ecfg['original_pair_epsilon']==evaluator.EPSILON,'Original epsilon conventions differ')
    require(family['initial_vector']==ecfg['initial_vector']==list(evaluator.START),'Six-value start differs')
    require(ecfg['baseline_vector']==training['baseline_parameters']==list(evaluator.INITIAL),'Historical baseline differs')
    require(family['original_vector_for_bounds']==training['family']['initial_vector']==list(evaluator.ORIGINAL_INITIAL),'Original bound origins differ')
    require(family['specification']==ecfg['specification']==PINS['pins']['six_specification'],'Six-family specification differs')
    require(len(family['parameters'])==6 and ecfg['parameter_names']==evaluator.NAMES,'Exactly six parameters required')
    for j,p in enumerate(family['parameters']):
        if j<4:
            oldp=training['family']['parameters'][j];require(all(p[k]==oldp[k]for k in ('name','lower','upper')),'Original four parameter bounds differ')
        else:
            epsilon=evaluator.EPSILON['La' if j==4 else 'Dy'];require(p['lower']==.5*epsilon and p['upper']==1.5*epsilon,'Epsilon bounds differ')
        require(p['name']==evaluator.NAMES[j] and p['initial']==evaluator.START[j],'Parameter order/start differs')


def candidate_gate(candidate):
    require(candidate['status']=='bounded_six_parameter_water_complete' and candidate['error'] is None and candidate['candidate_frozen'],'Only a complete error-free frozen six fit may bind')
    require(candidate['implementation']==PINS['pins']['six_source'] and candidate['config']==PINS['pins']['six_run_config'] and candidate['job_id']==PINS['candidate_job_id'],'Wrong actual six-fit producer')
    require(candidate['original_results_unchanged'] and candidate['development_checks_passed'] and candidate['selected']['score']['development_checks_passed'],'All original water controls must pass before binding')
    require(candidate['selected']['status']=='complete' and set(candidate['selected']['model'])==set(WATER),'Eight complete selected water rows required')
    require(candidate['vectors'][d.e.digest(candidate['selected']['parameters'])]==candidate['selected'],'Selected point absent from actual fit vector ledger')


def selected_authority(gid,value,row,inp,ecfg,candidate_path,previous,sources,values):
    pv=value['provenance'];ar=read(pv['native_authority']);metal=row['geometry']['metal'];i=0 if metal=='La' else 2
    effective=values[i:i+2]+[values[4 if metal=='La' else 5]]
    require(pv['requested_vector']==values and pv['effective_parameters']==effective,'Selected six-vector/effective values differ')
    require(pv['source_input']=={k:row['input'][k]for k in ('xyz','key','parameters')},'Selected source input differs')
    for k,p in pv['effective_input'].items():verify(p);require(p['sha256']==inp[k]['sha256'],'Selected effective input differs')
    request={'cell_id':gid,'effective_values':effective,'xyz':row['input']['xyz'],'parameters':row['input']['parameters'],'rendered_key_sha256':inp['key']['sha256'],'implementation':ecfg['implementation'],**{k:ecfg[k] for k in ('resolver','native_helper','executables','training_result')}}
    require(d.e.digest(request)==pv['cache_key'] and pv['request_cache_implementation']==ecfg['implementation'],'Current six-value cache identity differs')
    if pv['reused_baseline']:
        require(pv['native_authority']==PINS['pins']['previous_transfer_result'],'Baseline is not exact17578 authority')
        entry=previous['rows'][gid]
        require(entry['status']=='complete' and all(entry['model'][k]==value[k]for k in ('energy','forces','native')) and pv['effective_input']==entry['input']['files'],'Selected17578 baseline differs')
        require(pv['native_implementation'] is None,'Baseline native source relabeled')
        return 'exact17578_candidate_baseline'
    if pv.get('inherited_source_kind')=='17596_epsilon':
        epsilon=sources['17596'];matches=[]
        for probe_id,probe in epsilon['records'].items():
            cell=probe['cells'][gid]
            if cell.get('reused'):continue
            path=Path(PINS['pins']['epsilon_result']['path']).parent/probe_id/gid/'RESULT.json'
            if pin(path)==pv['native_authority']:matches.append(cell)
        require(len(matches)==1 and ar==matches[0] and ar['status']=='complete','Selected epsilon authority absent from exact17596')
        require(pv['native_cache_key'] is None and pv['native_implementation']==epsilon['implementation'] and pv['inherited_native'],'Epsilon native source relabeled')
        require(all(ar['model'][k]==value[k]for k in ('energy','forces','native')),'Epsilon raw model differs')
        for k,p in ar['input']['files'].items():verify(p);require(p['sha256']==inp[k]['sha256'],'Epsilon physical input differs')
        receipts=ar['receipts'];branch='inherited17596_epsilon_program_cell'
    else:
        require(ar['status']=='complete' and ar['cell_id']==gid and ar['cache_key']==pv['native_cache_key']==d.e.digest(ar['cache_identity']),'Selected cache authority differs')
        identity=ar['cache_identity']
        if pv['inherited_native']:
            kind=pv.get('inherited_source_kind');require(kind in ('17593_native','17589_native'),'Unknown inherited producer')
            label=kind[:5];parent=sources[label];parent_pin=PINS['pins']['six_prior_result' if label=='17593' else 'continuation_prior_result']
            producer=PINS['pins']['continuation_evaluator' if label=='17593' else 'radial_evaluator']
            require(Path(pv['native_authority']['path']).is_relative_to(Path(parent_pin['path']).parent/'cache/cells'),'Inherited native outside actual producer cache')
            require(any(v['provenance']['native_authority']==pv['native_authority']for vector in parent['vectors'].values()for k,v in vector['model'].items()if k==gid),'Inherited authority absent from exact producer inventory')
            require(identity['effective_values']==effective[:2] and effective[2]==d.e.EPSILON[metal],'Inherited fixed-epsilon input differs')
            branch='inherited'+label+'_native_cache'
        else:
            producer=PINS['pins']['six_evaluator'];require(Path(pv['native_authority']['path']).is_relative_to(Path(candidate_path).resolve().parent/'cache/cells'),'New native outside exact six-fit cache')
            require(identity['effective_values']==effective,'New six-cache parameters differ');branch='new_six_parameter_native_cache'
        require(pv['native_implementation']==producer==identity['implementation'],'Actual native producer relabeled')
        for k in ('native_helper','resolver','executables','training_result'):require(identity[k]==ecfg[k],'Selected engine/source differs: '+k)
        require(identity['cell_id']==gid and identity['xyz']==row['input']['xyz'] and identity['parameters']==row['input']['parameters'] and identity['rendered_key_sha256']==inp['key']['sha256'],'Native identity/input differs')
        require(ar['input']==pv['effective_input'] and all(ar[k]==value[k]for k in ('energy','forces','native')),'Selected cache values/input differ')
        receipts=ar['receipts']
    require(len(receipts)==2 and all(r['returncode']==0 for r in receipts),'Selected native receipts incomplete')
    for receipt in receipts:verify(receipt['output'])
    return branch


def previous_reuse(gid,inp,previous,row,native):
    prior=previous['rows'][gid]
    require(prior['status']=='complete','Prior transfer incomplete')
    if not all(prior['input']['files'][k]['sha256']==p['sha256']for k,p in inp.items()):return None
    for p in prior['input']['files'].values():verify(p)
    value=copy.deepcopy(prior['model']);replay_native(value,row,native)
    for r in prior['receipts']:require(r['returncode']==0,'Prior native return failure');verify(r['output'])
    check=None
    if gid in d.FD:
        check=prior.get('cartesian_check');require(check is not None,'Prior Cartesian evidence absent')
        again=gradient(verify(check['output']),row['geometry']['atoms'],True,native)
        require(again==check and check['fd_passed'],'Prior Cartesian gate differs')
        require(abs(check['energy']-value['energy'])<=.0001 and
            max(abs(x-y)for a,b in zip(check['forces'],value['forces'])for x,y in zip(a,b))<=.001,'Prior Cartesian E/F mismatch')
    return {'source':'identical_previous_transfer_input','authority':PINS['pins']['previous_transfer_result'],
        'value':value,'cartesian_check':check,'cartesian_check_authority':PINS['pins']['previous_transfer_result'] if check else None}


def bind(data_path,candidate_path,out):
    require(pin(data_path)==PINS['pins']['prepared32'],'Wrong original32 data')
    data=loaded(data_path);cfg=data['config'];candidate=json.loads(Path(candidate_path).read_text())
    candidate_gate(candidate)
    ccfg=read(candidate['config'])
    for key,name in (('implementation','six_source'),('specification','six_specification'),('objective','radial_objective'),('evaluator','six_evaluator'),('training_result','radial_training'),('evaluator_config','six_evaluator_config'),('preparation_result','six_preparation'),('prior_result','six_prior_result'),('epsilon_result','epsilon_result'),('root_release','six_root_release')):
        require(ccfg[key]==PINS['pins'][name],'Wrong six-fit binding: '+key);verify(ccfg[key])
    require(ccfg['body_helper']==cfg['body_helper'],'Candidate body weights differ')
    fit_source=module(ccfg['implementation'],'six_transfer_actual_source');release=read(ccfg['root_release'])
    require(release['six_status']=='approved_bounded_six_parameter_water' and release['six_config_sha256']==fit_source.fingerprint(ccfg) and release['six_preparation_result']==ccfg['preparation_result'],'Actual fit root release differs')
    ecfg=read(ccfg['evaluator_config']);require(ecfg['implementation']==ccfg['evaluator'] and ecfg['training_result']==ccfg['training_result'],'Six evaluator identity differs')
    for key in ('native_helper','resolver','executables'):require(ecfg[key]==cfg[key],'Candidate native context differs')
    native_sources={'17593':read(PINS['pins']['six_prior_result']),'17589':read(PINS['pins']['continuation_prior_result']),'17596':read(PINS['pins']['epsilon_result'])}
    fit_parent=native_sources['17593']
    require(candidate['prior_result']==PINS['pins']['six_prior_result'],'Wrong six-fit parent')
    require(candidate['initial']['parameters']==fit_parent['selected']['parameters']+[d.e.EPSILON['La'],d.e.EPSILON['Dy']] and candidate['initial']['score']==fit_parent['selected']['score'],'Six fit did not start at exact17593 selection')
    for gid in WATER:require(all(candidate['initial']['model'][gid][k]==fit_parent['selected']['model'][gid][k]for k in ('energy','forces','native')),'Six-fit initial raw model differs')
    verify(fit_parent['config']);verify(fit_parent['evaluator']['ledger'])
    selected=candidate['selected'];values=selected['parameters']
    require(candidate['development_checks_passed'] and selected['score']['development_checks_passed'],'This conditional transfer requires all unchanged water development controls to pass')
    require(selected['status']=='complete' and set(selected['model'])==set(WATER),'Eight complete selected water rows required')
    require(candidate['development_checks_passed']==selected['score']['development_checks_passed'],'Candidate control labels differ')
    training=read(ccfg['training_result']);family_gate(candidate['family'],training,ecfg)
    require(len(values)==6 and all(p['lower']<=v<=p['upper']for v,p in zip(values,candidate['family']['parameters'])),'Candidate parameter bounds differ')
    native=module(cfg['native_helper'],'radial_transfer_bind_native');body=module(cfg['body_helper'],'radial_transfer_bind_body')
    objective=module(ccfg['objective'],'radial_transfer_objective');center=module(ccfg['center_objective'],'radial_transfer_center')
    require(objective.evaluate(training,selected['model'],body)==selected['score'],'Selected radial objective replay differs')
    center_score=center.evaluate({'rows':{k:data['rows'][k]for k in d.CENTERS}}, {k:selected['model'][k]for k in d.CENTERS},body)
    require(center_score==candidate['selected_center_diagnostics']['original_four_center_score'],'Original center diagnostics differ')
    previous=read(PINS['pins']['previous_transfer_result']);previous_bound=read(previous['bound'])
    require(previous['status']=='consumed_transfer_complete' and previous['error'] is None and
        previous['implementation']==PINS['pins']['previous_run_source'] and previous_bound['data']==pin(data_path),'Previous transfer not admitted')
    inputs={};reusable={};authority=[];branches={}
    for gid in d.CELLS:
        row=data['rows'][gid];key,rendered=render(cfg,row,values);directory=out/'inputs'/gid;directory.mkdir(parents=True)
        shutil.copyfile(verify(row['input']['xyz']),directory/'system.xyz');shutil.copyfile(verify(row['input']['parameters']),directory/'amoeba09.prm');(directory/'system.key').write_text(key)
        inp={k:pin(directory/file)for k,file in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
        require(inp['xyz']['sha256']==row['input']['xyz']['sha256'] and inp['parameters']['sha256']==row['input']['parameters']['sha256'],'Input copy differs')
        inputs[gid]={'files':inp,'render':rendered,'input_identity':d.e.digest({'cell':gid,'xyz_sha256':inp['xyz']['sha256'],'key_sha256':inp['key']['sha256'],'prm_sha256':inp['parameters']['sha256'],'native':cfg['native_helper'],'executables':cfg['executables']})}
        if gid in WATER:
            value=copy.deepcopy(selected['model'][gid]);branches[gid]=selected_authority(gid,value,row,inp,ecfg,candidate_path,previous,native_sources,values)
            replay_native(value,row,native);authority.append(value['provenance']['native_authority'])
            reusable[gid]={'source':'selected_six_parameter_development','authority':value['provenance']['native_authority'],'value':value}
            prior_reuse=previous_reuse(gid,inp,previous,row,native)
            if prior_reuse and prior_reuse['cartesian_check'] is not None:
                reusable[gid]['cartesian_check']=prior_reuse['cartesian_check'];reusable[gid]['cartesian_check_authority']=prior_reuse['cartesian_check_authority']
        elif rendered['original_key_reused']:
            value=copy.deepcopy(row['baseline']);replay_native(value,row,native)
            reusable[gid]={'source':'unchanged_original_input','authority':row['baseline_result'],'value':value}
        else:
            prior_reuse=previous_reuse(gid,inp,previous,row,native)
            if prior_reuse:reusable[gid]=prior_reuse
    require(all(k in reusable and not inputs[k]['render']['water_recipient']for k in d.DRY),'Dry control reuse not established')
    require(all(k in reusable for k in WATER),'Eight water rows must be reusable before donor work')
    flags={'development_checks_passed':candidate['development_checks_passed'],'selected_score_development_passed':selected['score']['development_checks_passed'],
        'selected_per_fixture_controls':selected['score']['per_fixture'],'stopping_reason':candidate['stopping_reason'],
        'completed_outer_steps':candidate['completed_outer_steps'],'rank_records_at_iteration_origins':[{'origin_key':x['origin_key'],'linear_directions':x['linear_directions']}for x in candidate['iterations']],
        'final_point_rank_not_inferred':True,'at_parameter_bound':[v in(p['lower'],p['upper'])for v,p in zip(values,candidate['family']['parameters'])],
        'original_four_center_diagnostics':candidate['selected_center_diagnostics']}
    score=module(PINS['pins']['previous_score_source'],'original_transfer_score')
    previous_score=score.evaluate(data,{k:r['model']for k,r in previous['rows'].items()},body)
    require(previous_score==previous['comparison'],'Previous full32 score changed')
    bound={'status':'candidate_bound_no_native_calls','implementation':pin(__file__),'data_source':PINS['pins']['previous_data_source'],
        'score_source':PINS['pins']['previous_score_source'],'specification':PINS['pins']['plan'],'renderer':PINS['pins']['renderer'],'pin_inventory':pin(HERE/'PINS.json'),
        'data':pin(data_path),'candidate_result':pin(candidate_path),'candidate_config':candidate['config'],'parameters':values,
        'fit_parent_result':PINS['pins']['six_prior_result'],'selected_water_model':selected['model'],'candidate_flags':flags,'inputs':inputs,'reuse':reusable,
        'selected_native_authorities':authority,'selected_authority_branches':branches,'memberships':MEMBERSHIPS,
        'previous_transfer_result':PINS['pins']['previous_transfer_result'],'previous_candidate_comparison':previous_score,
        'maximum_native_programs':40,'maximum_internal_FD_energies':1380,'maximum_workers':8,'new_native_calls':0,'reserved_validation_opened':False}
    save(out/'BOUND.json',bound);return bound

def run(bound_path,release_path,out):
    start=time.monotonic();parent=time.process_time();child0=resource.getrusage(resource.RUSAGE_CHILDREN)
    bound=json.loads(Path(bound_path).read_text());data=loaded(verify(bound['data']));cfg=data['config'];release=json.loads(Path(release_path).read_text())
    require(bound['implementation']==pin(__file__) and bound['data_source']==PINS['pins']['previous_data_source'] and bound['score_source']==PINS['pins']['previous_score_source'],'Bound source changed')
    require(release['status']=='approved_consumed_six_parameter_water_transfer' and release['bound']==pin(bound_path) and release['implementation']==pin(__file__),'Exact root release required')
    require(release['maximum_native_programs']==40 and release['maximum_workers']==8 and release['reserved_validation_authorized'] is False,'Release scope differs')
    require(all(bound['candidate_flags'][k] for k in ('development_checks_passed','selected_score_development_passed')),'All original water controls must remain passed')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==8,'One root-owned eight-CPU Slurm task required')
    require(all(os.environ.get(k)=='1' for k in d.e.THREADS),'One numerical thread per native process required')
    verify(bound['candidate_result']);verify(bound['candidate_config']);verify(bound['fit_parent_result']);verify(bound['renderer']);verify(bound['specification']);verify(bound['pin_inventory']);verify(bound['previous_transfer_result'])
    require(bound['memberships']==MEMBERSHIPS and bound['maximum_native_programs']==40 and bound['maximum_internal_FD_energies']==1380,'Bound finite inventory differs')
    require(all(k in bound['reuse'] for k in WATER),'All eight selected water records must be reused')
    for p in bound['selected_native_authorities']:verify(p)
    native=module(cfg['native_helper'],'transfer_run_native');body=module(cfg['body_helper'],'transfer_run_body');score=module(bound['score_source'],'transfer_score')
    for gid,i in bound['inputs'].items():
        for p in i['files'].values():verify(p)
        key,rendered=render(cfg,data['rows'][gid],bound['parameters']);require(key==Path(i['files']['key']['path']).read_text() and rendered==i['render'],'Rendered input changed')
    for gid,reuse in bound['reuse'].items():verify(reuse['authority']);replay_native(reuse['value'],data['rows'][gid],native)
    result={'status':'running','bound':pin(bound_path),'release':pin(release_path),'implementation':pin(__file__),'job_id':os.environ['SLURM_JOB_ID'],'candidate_flags':bound['candidate_flags'],'memberships':MEMBERSHIPS,'previous_candidate_comparison':bound['previous_candidate_comparison'],'rows':{k:{'status':'not_attempted','receipts':[]}for k in d.CELLS},'native_program_invocations':0,'internal_FD_energies_requested':0,'internal_FD_energies_completed':0,'stage_A_passed':False,'error':None,'reserved_validation_opened':False,'physical_accuracy_qualified':False}
    lock=threading.RLock();stop=threading.Event()
    def persist():save(out/'RESULT.json',result)
    def fail(exc):
        stop.set()
        with lock:result['error']=repr(exc);persist()
    def cell(gid):
        row=data['rows'][gid];target=result['rows'][gid];fd=gid in d.FD;center=gid in d.CENTERS
        if stop.is_set():return
        try:
            directory=out/gid;directory.mkdir();meta=dict(row['input'],derivative_test=fd)
            for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(bound['inputs'][gid]['files'][k]),directory/name)
            with lock:target.update(status='started',input=bound['inputs'][gid]);persist()
            reuse=bound['reuse'].get(gid);old=reuse['value']if reuse else None
            if old:
                saved_fd=reuse.get('cartesian_check')
                if fd and saved_fd is not None:
                    verify(reuse['cartesian_check_authority'])
                    g=gradient(verify(saved_fd['output']),row['geometry']['atoms'],True,native)
                    require(g==saved_fd and g['fd_passed'],'Saved Cartesian check differs')
                    require(abs(g['energy']-old['energy'])<=.0001 and max(abs(x-y)for a,b in zip(g['forces'],old['forces'])for x,y in zip(a,b))<=.001,'Saved Cartesian E/F differs')
                    with lock:target.update(status='complete',model=old,reused_from=reuse['authority'],cartesian_check=g,cartesian_check_authority=reuse['cartesian_check_authority']);persist()
                    return
                ap=verify(old['native']['analyze_output']);gp=verify(old['native']['gradient_output'])
                if not fd or old['native']['finite_difference'] is not None:
                    extra={}
                    if fd:
                        g=gradient(gp,row['geometry']['atoms'],True,native);require(g['fd_passed'],'Reused Cartesian force check fails');extra['cartesian_check']=g
                    with lock:target.update(status='complete',model=old,reused_from=reuse['authority'],**extra);persist()
                    return
            else:ap=directory/'analyze.out'
            programs=[]if old else[('analyze',['system.xyz','E'])]
            programs.append(('testgrad',['system.xyz','Y','Y','0.00001']if fd else['system.xyz','Y','N']))
            for program,args in programs:
                with lock:
                    require(not stop.is_set(),'Global native stop active');require(result['native_program_invocations']<40,'40-program ceiling')
                    result['native_program_invocations']+=1
                    if program=='testgrad' and fd:result['internal_FD_energies_requested']+=6*row['geometry']['atoms'];require(result['internal_FD_energies_requested']<=1380,'FD energy ceiling')
                    target['active_program']={'program':program,'args':args,'status':'started'}
                    persist()
                output=directory/(program+'.out');local_receipts=[];updates={}
                try:
                    native.invoke([cfg['executables'][program]['path']]+args,directory,output,local_receipts);program_admission(output,program,native)
                    if program=='testgrad':
                        g=gradient(output,row['geometry']['atoms'],fd,native)
                        if fd:require(g['fd_passed'],'Cartesian force consistency failed')
                        if center:
                            require(abs(g['energy']-old['energy'])<=.0001,'Frozen-center energy replay failed')
                            require(max(abs(x-y)for a,b in zip(g['forces'],old['forces'])for x,y in zip(a,b))<=.001,'Frozen-center force replay failed')
                            updates.update(cartesian_check=g,model=old)
                        else:
                            parsed=native.parse(ap,output,meta);updates['model']={'energy':parsed['total_energy_kcal_mol'],'forces':parsed['force_kcal_mol_A'],'native':parsed}
                            if fd:updates['cartesian_check']=g
                        if fd:
                            with lock:result['internal_FD_energies_completed']+=6*row['geometry']['atoms']
                except Exception as exc:fail(exc);raise
                finally:
                    with lock:
                        target['receipts'].extend(local_receipts);target['active_program']['status']='returned';target.update(updates);persist()
            with lock:target.update(status='complete');persist()
        except Exception as exc:
            with lock:target.update(status='failed',error=repr(exc));fail(exc)
    try:
        for gid in d.DRY:result['rows'][gid].update(status='complete',model=bound['reuse'][gid]['value'],reused_from=bound['reuse'][gid]['authority'])
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=4)as pool:list(pool.map(cell,d.CENTERS))
        require(not stop.is_set() and all(result['rows'][k]['status']=='complete' for k in d.CENTERS),'Fitted-center gate failed')
        for gid in d.RADIAL:cell(gid)
        require(not stop.is_set() and all(result['rows'][k]['status']=='complete' for k in WATER),'Eight water admission failed')
        result['stage_A_passed']=True;persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=8)as pool:list(pool.map(cell,d.HYDRATED))
        require(not stop.is_set() and all(r['status']=='complete'for r in result['rows'].values()),'Transfer native/Cartesian gate failed')
        result['comparison']=score.evaluate(data,{k:r['model']for k,r in result['rows'].items()},body);result['status']='consumed_six_parameter_transfer_complete'
        result['historical_score_partitions_not_independent_validation']=True
        for p in (bound['candidate_result'],bound['candidate_config'],bound['previous_transfer_result'],bound['data']):verify(p)
        result['protected_parent_pins_unchanged']=True
    except Exception as exc:result.update(status='native_or_cartesian_failure',error=repr(exc));stop.set()
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);elapsed=time.monotonic()-start
        result.update(elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*8,parent_cpu_seconds=time.process_time()-parent,aggregate_child_user_cpu_seconds=ch.ru_utime-child0.ru_utime,aggregate_child_system_cpu_seconds=ch.ru_stime-child0.ru_stime,per_invocation_child_cpu_may_overlap_not_summed=True,global_stop=stop.is_set());persist()
    return result

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('bind');a.add_argument('--data',required=True);a.add_argument('--candidate',required=True);a.add_argument('--out',required=True)
    a=sub.add_parser('run');a.add_argument('--bound',required=True);a.add_argument('--release',required=True);a.add_argument('--out',required=True)
    args=ap.parse_args();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:
        if args.action=='bind':r=bind(args.data,args.candidate,out)
        else:r=run(args.bound,args.release,out)
    except Exception as exc:
        r={'status':'precall_admission_failed','error':repr(exc),'implementation':pin(__file__),'rows':{k:{'status':'not_attempted'}for k in d.CELLS},'new_native_calls':0};save(out/'FAILURE.json',r)
        print(json.dumps({'status':r['status'],'result':pin(out/'FAILURE.json'),'error':r['error']}));return 1
    print(json.dumps({'status':r['status'],'result':pin(out/('BOUND.json'if args.action=='bind'else'RESULT.json'))}));return 0 if r['status']in ('candidate_bound_no_native_calls','consumed_six_parameter_transfer_complete')else 1

if __name__=='__main__':raise SystemExit(main())
