#!/usr/bin/env python3
"""Four original water orientations; exact passing historical water model only."""
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

HERE=Path(__file__).resolve().parent;EXPLICIT=HERE.parents[1];ROOT=EXPLICIT.parents[2]
WORK=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
CELLS=tuple('validation_aquo8_orientation_'+s+'_'+m for s in ('minus','plus') for m in ('La','Dy'))
ANCHORS=('aquo8_central_La','aquo8_central_Dy')
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
WATER=(3.906136716492464,.28782327521600637,3.5865909437359376,.2837564269122712,.20765190997435956,.3173248056679709)

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def save(path,value):
    p=Path(path);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def default_config():
    fit=read(pin(WORK/'paired_water_six_parameter_v1/run_v1/execution/RESULT.json'));fcfg=read(fit['config']);ecfg=read(fcfg['evaluator_config'])
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'root_declaration':pin(HERE.parent/'ROOT_DECLARATION.md'),
        'water_result':pin(WORK/'paired_water_six_parameter_v1/run_v1/execution/RESULT.json'),
        'water_review':pin(WORK/'paired_water_six_parameter_v1/result_review_v1/REVIEW.json'),
        'transfer_result':pin(WORK/'paired_water_six_transfer_v1/run_v1/execution/RESULT.json'),
        'transfer_review':pin(WORK/'paired_water_six_transfer_v1/actual_transfer_review_v1/REVIEW.json'),
        'reserved_panel':pin(WORK/'matched_dry_v1/reserved_validation_v1/RESERVED_PANEL.json'),
        'reserved_native_panel':pin(WORK/'matched_dry_v1/reserved_validation_v1/NATIVE_PANEL.json'),
        'typing_helper':pin(EXPLICIT/'paired_donor_validation_stageA_v1/model/RUN.py'),
        'numerical_helper':pin(EXPLICIT/'paired_water_transfer_v1/RUN.py'),'native_helper':ecfg['native_helper'],
        'body_helper':fcfg['body_helper'],'resolver':ecfg['resolver'],'executables':ecfg['executables'],
        'maximum_native_programs':8,'maximum_workers':4,'new_FD_energies':0,'root_release':None}

def context(cfg):
    require(cfg==default_config(),'Frozen source/input configuration differs')
    expected={'water_result':'6027e319a261eb85171fb06e30aa3e519411f6371a7905eeca0e8a79812c48fc',
        'water_review':'a3147adcd7448012437e3ac275872fa58fe65e98a114619a95ec01e17a0df187',
        'transfer_result':'e3eeb870aac22f408c8b38a5f050404cfee5a8f4c0468633c07122fbc5e07254',
        'transfer_review':'3f7b3a9b375cf3cdad0580a75852ce9488a664232aa1da8b1254b547ac2b4ce4',
        'reserved_panel':'eed99b0fb6dbe7b49857d68853fe809481f6b657c928a7e75197d77c61f3432b',
        'reserved_native_panel':'75d3fc69703c8a5acca61f547e3fcf5d614cfe6ac653d2f3ecd3f34a93bcaf5a',
        'typing_helper':'33a085bef1b51dac097cc398b78c298457c84f3d7f8f4e002ea5072c87166d27'}
    for key,value in expected.items():require(cfg[key]['sha256']==value,'Unexpected authority: '+key);verify(cfg[key])
    for key in ('specification','root_declaration','native_helper','numerical_helper','body_helper','resolver'):verify(cfg[key])
    for p in cfg['executables'].values():verify(p)
    water=read(cfg['water_result']);transfer=read(cfg['transfer_result']);wr=read(cfg['water_review']);tr=read(cfg['transfer_review'])
    require(water['status']=='bounded_six_parameter_water_complete' and water['job_id']=='17605' and water['error'] is None and water['candidate_frozen'] and water['development_checks_passed'],'Passing17605 unavailable')
    require(water['selected']['parameters']==list(WATER) and water['selected']['score']['development_checks_passed'],'Frozen water selection differs')
    require(wr['status']=='saved_result_review_passed' and wr['result']==cfg['water_result'] and wr['all_selected_raw_replays_exact'] and wr['exact_evaluator_replay'],'Water independent review differs')
    require(transfer['status']=='consumed_six_parameter_transfer_complete' and transfer['job_id']=='17612' and transfer['error'] is None and transfer['stage_A_passed'] and not transfer['global_stop'],'Actual17612 numerical transfer unavailable')
    require(tr['status']=='passed_actual_raw_numerical_and_comparison_review' and tr['actual_result']==cfg['transfer_result'],'Transfer independent review differs')
    bound=read(transfer['bound']);require(bound['candidate_result']==cfg['water_result'] and bound['parameters']==list(WATER),'Actual water/transfer join differs')
    data=read(bound['data']);fcfg=read(water['config']);ecfg=read(fcfg['evaluator_config']);training=read(fcfg['training_result'])
    require(fcfg['evaluator']['sha256']=='c73e62c1b6eaba44c3a595c2f14f682c8311c721966d5e97e27d612c6b69c3dd','Original water evaluator differs')
    for key in ('implementation','evaluator','objective','body_helper','specification','evaluator_config','root_release'):verify(fcfg[key])
    native=module(cfg['native_helper'],'orientation_native');helper=module(cfg['numerical_helper'],'orientation_numeric');body=module(cfg['body_helper'],'orientation_body');obj=module(fcfg['objective'],'orientation_water_objective')
    require(obj.evaluate(training,water['selected']['model'],body)==water['selected']['score'],'Water selected objective replay differs')
    require(water['vectors'][digest(list(WATER))]==water['selected'],'Selected vector absent from exact ledger')
    for gid,value in water['selected']['model'].items():
        pv=value['provenance'];ar=read(pv['native_authority']);meta=data['rows'][gid];inp=bound['inputs'][gid]['files']
        require(not pv['reused_baseline'] and not pv['inherited_native'],'Unexpected native branch for this exact selection')
        require(pv['requested_vector']==list(WATER) and pv['native_implementation']==fcfg['evaluator']==ar['cache_identity']['implementation'],'Water producer differs')
        require(ar['status']=='complete' and ar['cell_id']==gid and ar['cache_key']==pv['native_cache_key']==pv['cache_key']==digest(ar['cache_identity']),'Selected cache identity differs')
        require(pv['native_authority']==pin(Path(cfg['water_result']['path']).parent/'cache/cells'/pv['cache_key']/'RESULT.json'),'Selected authority path differs')
        require(ar['input']==pv['effective_input'] and all(ar[k]==value[k] for k in ('energy','forces','native')),'Selected raw/input values differ')
        for key,p in inp.items():verify(p);verify(pv['effective_input'][key]);require(p['sha256']==pv['effective_input'][key]['sha256'],'Transferred water input bytes differ')
        require(all(value[k]==transfer['rows'][gid]['model'][k] for k in ('energy','forces','native')),'Transferred water E/F differs')
        require(len(ar['receipts'])==2 and all(r['returncode']==0 for r in ar['receipts']),'Selected native receipts incomplete')
        for receipt in ar['receipts']:verify(receipt['output'])
        helper.replay_native(value,meta,native)
    cartesian={}
    for gid in ANCHORS:
        row=transfer['rows'][gid];require(row['status']=='complete' and row['input']['files']==bound['inputs'][gid]['files'],'Aquo numerical input differs')
        for receipt in row['receipts']:require(receipt['returncode']==0,'Cartesian program failed');verify(receipt['output'])
        g=helper.gradient(verify(row['cartesian_check']['output']),25,True,native)
        require(g==row['cartesian_check'] and g['fd_passed'],'Aquo Cartesian consistency failed')
        value=row['model'];require(abs(g['energy']-value['energy'])<=1e-4 and max(abs(x-y) for a,b in zip(g['forces'],value['forces']) for x,y in zip(a,b))<=1e-3,'Aquo saved energy/force replay failed')
        cartesian[gid]=g
    panel=read(cfg['reserved_panel']);npanel=read(cfg['reserved_native_panel']);require(npanel['source_inventory']==cfg['reserved_panel'] and panel['status']=='reserved_complete','Original geometry preparation unavailable')
    geometries={g['geometry_id']:g for g in panel['geometries'] if g['geometry_id'] in CELLS};checks={x['id']:x for x in panel['checks']}
    require(set(geometries)==set(CELLS),'Four reserved water slots required')
    for gid,g in geometries.items():
        require(checks[gid]['passed'] and g['validation_status']=='reserved' and g['atoms']==25 and g['charge']==3 and g['physical_multiplicity']==(1 if g['metal']=='La' else 6),'Reserved water state/gate differs');verify(g['xyz'])
        require(any(c['geometry']==g and c['recipe']=='ecp_pbe0' for c in npanel['cells']),'Reserved native geometry join differs')
    return data,geometries,water,transfer,bound,cartesian

def bind(config_path,out):
    cfg=json.loads(Path(config_path).read_text());data,geometries,water,transfer,previous,cartesian=context(cfg)
    typing=module(cfg['typing_helper'],'orientation_typing');body=module(cfg['body_helper'],'orientation_body_bind');resolver=module(cfg['resolver'],'orientation_resolver')
    rows={}
    for gid in CELLS:
        g=geometries[gid];metal=g['metal'];parent='aquo8_central_'+metal;source=data['rows'][parent];old=previous['inputs'][parent]['files']
        require(set(source['input']['types'])=={36,37,500 if metal=='La' else 501} and len(source['input']['types'])==25,'Nonwater atom type present')
        directory=out/'inputs'/gid;directory.mkdir(parents=True);(directory/'system.xyz').write_text(typing.typed_xyz(g,source))
        for key,name in (('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(old[key]),directory/name)
        files={k:pin(directory/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
        require(all(files[k]['sha256']==old[k]['sha256'] for k in ('key','parameters')),'Water key/parameter bytes changed')
        parameters=resolver.parse_parameters(verify(files['parameters']).read_text(),verify(files['key']).read_text());pair=resolver.resolve_pair(parameters,500 if metal=='La' else 501,36);i=0 if metal=='La' else 2
        require([pair['Rmin_pair_A'],pair['thole'],pair['epsilon_kcal_mol']]==[WATER[i],WATER[i+1],WATER[4 if metal=='La' else 5]],'Selected water pair differs')
        meta=copy.deepcopy(source['input']);meta.update(id=gid,derivative_test=False,source_geometry=cfg['reserved_panel'],**files)
        rows[gid]={'geometry':g,'input':meta,'bodies':copy.deepcopy(source['bodies']),'files':files,'parent_training_id':parent,'resolved_water_pair':pair}
        require(body.bodies(rows[gid])==rows[gid]['bodies'],'Water body graph differs')
    result={'status':'frozen_water_orientations_bound_no_calls','implementation':pin(__file__),'config':pin(config_path),
        'water_result':cfg['water_result'],'water_review':cfg['water_review'],'transfer_result':cfg['transfer_result'],'transfer_review':cfg['transfer_review'],
        'parameters':list(WATER),'rows':rows,'anchor_model':{gid:copy.deepcopy(transfer['rows'][gid]['model']) for gid in ANCHORS},
        'aquo_cartesian_checks':cartesian,'all_eight_selected_water_raw_replayed':True,'new_native_calls':0,'reserved_outcomes_opened':False,
        'full16_validation_complete':False,'donor_transfer_failures_preserved':not transfer['comparison']['all_development_checks_passed']}
    save(out/'BOUND.json',result);return result

def bound_inputs(cfg,bound):
    require(cfg==default_config() and bound['implementation']==pin(__file__) and bound['status']=='frozen_water_orientations_bound_no_calls','Bound source mismatch')
    for p in cfg['executables'].values():verify(p)
    require(set(bound['rows'])==set(CELLS) and bound['parameters']==list(WATER) and bound['all_eight_selected_water_raw_replayed'],'Bound scope mismatch')
    for key in ('water_result','water_review','transfer_result','transfer_review'):require(bound[key]==cfg[key],'Bound authority mismatch');verify(cfg[key])
    for row in bound['rows'].values():
        for p in row['files'].values():verify(p)

def admit_result(config,bound,result):
    """Pure saved admission. Dict inputs; four availability rows, no engines."""
    bound_inputs(config,bound);require(result['implementation']==pin(__file__) and read(result['bound'])==bound,'Actual result/bound mismatch')
    require(set(result['rows'])==set(CELLS) and 0<=result['native_program_invocations']<=8 and result['new_FD_energies']==0,'Actual finite scope mismatch')
    release=read(result['release']);require(release['implementation']==pin(__file__) and release['bound']==result['bound'] and release['status']=='approved_frozen_water_orientation_model','Actual release differs')
    native=module(config['native_helper'],'orientation_saved_native');helper=module(config['numerical_helper'],'orientation_saved_numeric');admitted={}
    for gid in CELLS:
        row=result['rows'][gid]
        try:
            require(row['status']=='complete' and row['input']==bound['rows'][gid]['files'],'Native cell unavailable/input differs')
            require(len(row['receipts'])==2 and all(x['returncode']==0 for x in row['receipts']),'Native receipt incomplete')
            for record,program,args in zip(row['receipts'],('analyze','testgrad'),(['system.xyz','E'],['system.xyz','Y','N'])):
                require(record['argv']==[config['executables'][program]['path']]+args,'Native command differs');verify(record['output'])
            value=row['model'];helper.replay_native(value,bound['rows'][gid],native)
            require(value['native']['finite_difference'] is None,'Unexpected reserved finite differences')
            admitted[gid]={'status':'complete','model':value,'input':row['input'],'receipts':row['receipts']}
        except Exception as exc:admitted[gid]={'status':'unavailable','model':None,'reason':repr(exc)}
    complete=sum(x['status']=='complete' for x in admitted.values())
    require(complete<4 or result['status']=='water_orientation_model_complete' and result['error'] is None and not result['global_stop'],'Contradictory terminal labels')
    return {'status':'saved_model_admission_complete' if complete==4 else 'saved_model_admission_partial','rows':admitted,'complete':complete,
        'anchor_model':bound['anchor_model'],'new_native_calls':0,'full16_validation_complete':False,'physical_accuracy_qualified':False}

def run(bound_path,release_path,out):
    started=time.monotonic();cpu=time.process_time();children=resource.getrusage(resource.RUSAGE_CHILDREN)
    bp=pin(bound_path);bound=read(bp);cfg=read(bound['config']);bound_inputs(cfg,bound);release=json.loads(Path(release_path).read_text())
    require(release['status']=='approved_frozen_water_orientation_model' and release['bound']==bp and release['implementation']==pin(__file__) and release['water_result']==cfg['water_result'],'Exact root release required')
    require(release['maximum_native_programs']==8 and release['maximum_workers']==4 and release['new_FD_energies']==0,'Release ceiling mismatch')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==4 and all(os.environ.get(k)=='1' for k in THREADS),'Exact4CPU/single-thread runtime required')
    native=module(cfg['native_helper'],'orientation_runtime_native');helper=module(cfg['numerical_helper'],'orientation_runtime_numeric')
    result={'status':'running','implementation':pin(__file__),'bound':bp,'release':pin(release_path),'water_result':cfg['water_result'],'transfer_result':cfg['transfer_result'],
        'job_id':os.environ['SLURM_JOB_ID'],'rows':{g:{'status':'not_attempted','receipts':[]} for g in CELLS},'native_program_invocations':0,'new_FD_energies':0,
        'error':None,'anchor_model':bound['anchor_model'],'full16_validation_complete':False,'physical_accuracy_qualified':False}
    lock=threading.RLock();stop=threading.Event()
    def persist():save(out/'RESULT.json',result)
    def fail(exc):
        stop.set()
        with lock:result['error']=repr(exc);persist()
    def cell(gid):
        if stop.is_set():return
        target=result['rows'][gid];row=bound['rows'][gid]
        try:
            directory=out/gid;directory.mkdir()
            for key,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
                shutil.copyfile(verify(row['files'][key]),directory/name);require(pin(directory/name)['sha256']==row['files'][key]['sha256'],'Runtime copy differs')
            with lock:target.update(status='started',input=row['files']);persist()
            for program,args in (('analyze',['system.xyz','E']),('testgrad',['system.xyz','Y','N'])):
                with lock:require(not stop.is_set(),'Global stop active');require(result['native_program_invocations']<8,'Native ceiling');result['native_program_invocations']+=1;target['active_program']=program;persist()
                receipts=[];output=directory/(program+'.out')
                try:
                    native.invoke([cfg['executables'][program]['path']]+args,directory,output,receipts);helper.program_admission(output,program,native)
                    if program=='testgrad':
                        n=native.parse(directory/'analyze.out',output,dict(row['input'],derivative_test=False))
                        with lock:target['model']={'energy':n['total_energy_kcal_mol'],'forces':n['force_kcal_mol_A'],'native':n,'provenance':{'bound':bp,'input':row['files'],'water_result':cfg['water_result']}}
                except Exception as exc:fail(exc);raise
                finally:
                    with lock:target['receipts'].extend(receipts);persist()
            with lock:target['status']='complete';persist()
        except Exception as exc:
            stop.set()
            with lock:target.update(status='failed',error=repr(exc));result['error']=repr(exc);persist()
    try:
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(cell,CELLS))
        require(not stop.is_set() and all(r['status']=='complete' for r in result['rows'].values()),'Incomplete water orientation cells')
        result['status']='water_orientation_model_complete'
        verify(cfg['water_result']);verify(cfg['transfer_result']);result['original_results_unchanged']=True
    except Exception as exc:stop.set();result.update(status='native_or_admission_failure',error=repr(exc))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result.update(global_stop=stop.is_set(),elapsed_seconds=time.monotonic()-started,parent_CPU_seconds=time.process_time()-cpu,aggregate_child_user_CPU_seconds=ch.ru_utime-children.ru_utime,aggregate_child_system_CPU_seconds=ch.ru_stime-children.ru_stime,per_program_child_CPU_may_overlap_not_summed=True);persist()
    return result

def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('bind');a.add_argument('--config',required=True)
    a=s.add_parser('run');a.add_argument('--bound',required=True);a.add_argument('--release',required=True)
    for a in s.choices.values():a.add_argument('--out',required=True)
    args=p.parse_args();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:r=bind(args.config,out) if args.action=='bind' else run(args.bound,args.release,out)
    except Exception as exc:r={'status':'adapter_failed_or_unavailable','error':repr(exc),'new_native_calls':0 if args.action=='bind' else None,'rows':{g:{'status':'not_attempted'} for g in CELLS}};save(out/'FAILURE.json',r)
    print(json.dumps({'status':r['status'],'out':str(out),'error':r.get('error')}));return 0 if r['status'] in ('frozen_water_orientations_bound_no_calls','water_orientation_model_complete') else 1
if __name__=='__main__':raise SystemExit(main())
