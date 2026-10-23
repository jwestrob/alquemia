#!/usr/bin/env python3
"""One terminal saved replay of40 new cells/40 reused entries/eight scores."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[5]
REL=Path('lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_coordination_joint_v1/metal_alpha_response_v1')
WORK=ROOT/'workspaces'/REL;DIAG=ROOT/'diagnostics'/REL
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,p['path'];return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def save(p,x):
    with p.open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
    return pin(p)

def main():
    wall,cpu=time.monotonic(),time.process_time();assert all(os.environ.get(k)=='1' for k in THREADS)
    out=WORK/'actual_review_v1';out.mkdir(exist_ok=False);checks=[]
    def check(name,ok):
        checks.append({'name':name,'passed':bool(ok)});assert ok,name
    try:
        parent=pin(WORK/'run_v1/execution/RESULT.json');r=read(parent);configpin=pin(WORK/'run_v1/CONFIG_RUN.json');cfg=read(configpin)
        check('terminal actual application',r['status']=='metal_alpha_diagnostic_complete' and r['job_id']=='17695' and r['error'] is None and r['global_stop'] is False)
        check('exact submitted source/config',r['implementation']==cfg['implementation'] and r['config']==cfg and cfg['implementation']['sha256']=='54cb4777d81d48bf7215fe8cf5e357f2a1a4e6b6cf99123b05c6b07a605d9be3')
        completion_source=WORK/'run_v1/completion_event.json';event_bytes=completion_source.read_bytes();(out/'completion_snapshot.json').write_bytes(event_bytes);event=json.loads(event_bytes);scheduler=event['records']['17695']
        check('terminal scheduler completion',scheduler['state']=='COMPLETED' and scheduler['exit_code']=='0:0' and int(scheduler['allocated_cpus'])==8)
        probe=module(cfg['implementation'],'actual_alpha_probe');probe.checked_config(cfg)
        prep=read(cfg['preparation_result']);manifest=read(prep['manifest']);data=read(manifest['data']);release=read(cfg['root_release']);review=read(cfg['selected_review']);training=read(cfg['training'])
        check('exact prepared authority',r['manifest']==prep['manifest'] and r['preparation_result']==cfg['preparation_result'] and prep['status']=='prepared' and prep['implementation']==cfg['implementation'] and manifest['scientific_config_sha256']==probe.fingerprint(cfg))
        check('root release finite scope',release['status']=='approved_fixed_metal_alpha_diagnostic' and release['implementation']==cfg['implementation'] and release['preparation_result']==cfg['preparation_result'] and release['scientific_config_sha256']==probe.fingerprint(cfg) and release['workers']==8 and release['maximum_new_physical_cells']==40 and release['maximum_new_native_programs']==80 and release['candidate_selection_authorized'] is False)
        check('unchanged consumed references/data',training['rows']==data['rows'] and set(data['rows'])==set(probe.CELLS) and data['provenance']['candidate_result']==cfg['candidate_result'] and data['provenance']['selected_review']==cfg['selected_review'])
        check('raw-pass science-fail retained',review['result']==cfg['candidate_result'] and review['status']=='selected_ten_raw_review_passed' and review['all_selected_raw_replays_exact'] and review['exact_selected_score_replayed'] and review['all_development_controls_passed'] is False and data['baseline_score']['development_checks_passed'] is False)
        native=module(cfg['native_helper'],'actual_alpha_native');resolver=module(cfg['resolver'],'actual_alpha_resolver');obj=module(cfg['objective'],'actual_alpha_objective');body=module(cfg['body_helper'],'actual_alpha_body');auth=module(cfg['reviewer'],'actual_alpha_producer')
        ev=module(cfg['candidate_evaluator'],'actual_alpha_admission');guard=object.__new__(ev.Evaluator);guard.native=native
        baseline_evidence={};reference_evidence={}
        for gid,row in data['rows'].items():
            value=data['baseline_model'][gid];producer=auth.native_authority(gid,value)
            parsed=native.parse(verify(value['native']['analyze_output']),verify(value['native']['gradient_output']),dict(row['input'],derivative_test=False))
            check('unique baseline raw replay '+gid,parsed==value['native'] and parsed['total_energy_kcal_mol']==value['energy'] and parsed['force_kcal_mol_A']==value['forces'])
            baseline_evidence[gid]=producer;rr=read(row['reference_result'])
            check('reference authority available '+gid,rr['status']=='complete' and rr['returncode']==0 and all(rr['checks'].values()))
            for key in ('output','engrad','input','geometry','panel','adapter'):verify(rr[key])
            reference_evidence[gid]={'result':row['reference_result'],'raw_pins':{key:rr[key] for key in ('output','engrad','input','geometry','panel','adapter')},'unchanged_exact_training_reference':True}
        check('eight complete full matrices',set(r['records'])=={p['id'] for p in manifest['probes']} and all(v['status']=='complete' and set(v['model'])==set(v['cells'])==set(probe.CELLS) for v in r['records'].values()))
        native_cells=[];reused=[];scores={};receipts=0
        for spec in manifest['probes']:
            name=spec['id'];record=r['records'][name]
            check('probe coordinates '+name,all(record[k]==spec[k] for k in ('metal','alpha','eta','fraction')) and spec['alpha']==probe.ALPHA[spec['metal']]*(1+spec['fraction']))
            for gid in spec['reused_cells']:
                value=record['model'][gid];cell=record['cells'][gid]
                check('unchanged reused entry '+name+'/'+gid,value==data['baseline_model'][gid] and cell=={'status':'complete','reused':True,'authority':value['provenance']['native_authority']})
                reused.append({'probe':name,'cell':gid,'authority':cell['authority']})
            for gid,inp in spec['inputs'].items():
                directory=WORK/'run_v1/execution'/name/gid;authority=pin(directory/'RESULT.json');cell=read(authority);value=record['model'][gid]
                check('native cell authority '+name+'/'+gid,cell==record['cells'][gid] and cell['status']=='complete' and cell['cell_id']==gid and cell['input']==inp and cell['model']==value)
                key,render=probe.alpha_key(verify(inp['baseline_input']['key']).read_text(),verify(inp['baseline_input']['parameters']).read_text(),spec['metal'],spec['alpha'],resolver)
                check('exact active alpha token '+name+'/'+gid,key==verify(inp['files']['key']).read_text()==verify(cell['runtime_input']['key']).read_text() and render==inp['render'])
                for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
                    check('actual input bytes '+name+'/'+gid+'/'+k,pin(directory/n)==cell['runtime_input'][k] and cell['runtime_input'][k]['sha256']==inp['files'][k]['sha256'])
                check('unchanged physical source '+name+'/'+gid,all(cell['runtime_input'][k]['sha256']==inp['baseline_input'][k]['sha256'] for k in ('xyz','parameters')))
                check('two native receipts '+name+'/'+gid,len(cell['receipts'])==2)
                for receipt,program,args in zip(cell['receipts'],('analyze','testgrad'),(['system.xyz','E'],['system.xyz','Y','N'])):
                    check('actual argv/threads '+name+'/'+gid+'/'+program,receipt['argv']==[cfg['executables'][program]['path']]+args and receipt['returncode']==0 and Path(receipt['cwd'])==directory and all(receipt['threads'].get(k)=='1' for k in THREADS))
                    output=verify(receipt['output']);check('actual output '+name+'/'+gid+'/'+program,output==directory/(program+'.out'));guard._program_admission(output,program);receipts+=1
                parsed=native.parse(directory/'analyze.out',directory/'testgrad.out',dict(data['rows'][gid]['input'],derivative_test=False))
                check('native energy force replay '+name+'/'+gid,parsed==value['native'] and parsed['finite_difference'] is None and parsed['total_energy_kcal_mol']==value['energy'] and parsed['force_kcal_mol_A']==value['forces'])
                check('new producer provenance '+name+'/'+gid,value['provenance']=={'probe_id':name,'alpha':spec['alpha'],'baseline':cfg['candidate_result'],'effective_input':cell['runtime_input'],'native_implementation':cfg['implementation']})
                native_cells.append({'probe':name,'cell':gid,'authority':authority,'runtime_input':cell['runtime_input'],'native_output':{k:parsed[k] for k in ('analyze_output','gradient_output')}})
            score=obj.evaluate({'rows':data['rows']},record['model'],body);check('full unchanged score replay '+name,score==record['score'] and len(score['weighted_residuals'])==413 and score['declared_energy_controls']==10 and score['declared_force_controls']==15)
            scores[name]={'objective':score['objective'],'W9_balanced_error_kcal_mol':score['energy_controls']['balanced_exchange']['aquo9_central']['error_kcal_mol'],'energy_checks_passed':sum(v['passed'] for values in score['energy_controls'].values() for v in values.values()),'force_checks_passed':sum(v['passed'] for values in score['per_fixture'].values() for v in values.values()),'all_controls_passed':score['development_checks_passed'],'energy_controls':score['energy_controls'],'force_controls':{kind:{g:{k:v[k] for k in ('passed','rigid_RMS_kcal_mol_A','maximum_body_component_error_kcal_mol_A')} for g,v in values.items()} for kind,values in score['per_fixture'].items()}}
        counter=read(pin(WORK/'run_v1/execution/CALL_COUNT.json'))
        check('complete finite denominator',len(native_cells)==40 and len(reused)==40 and len(scores)==8 and receipts==r['native_program_invocations']==counter['native_program_invocations']==counter['maximum']==80)
        check('no candidate selection or fresh outcome',r['new_candidate_selected'] is False and r['reserved_outcomes_opened'] is False)
        for parent_pin in manifest['protected_parents']:verify(parent_pin)
        check('parent scientific bytes unchanged',pin(parent['path'])==parent and r['protected_parent_pins_unchanged'])
        result={'status':'actual_alpha_raw_score_review_passed','implementation':pin(__file__),'result':parent,'configuration':configpin,'root_release':cfg['root_release'],'preparation':cfg['preparation_result'],'manifest':r['manifest'],'data':manifest['data'],'completion_snapshot':pin(out/'completion_snapshot.json'),'source_review':pin(WORK/'independent_review_v1/REVIEW.json'),'checks':checks,'check_count':len(checks),'new_cells_raw_replayed':40,'new_program_receipts_admitted':receipts,'reused_entries_exact':40,'unique_reused_producers_replayed':10,'full_scores_exact':8,'native_cells':native_cells,'reused_entries':reused,'baseline_producers':baseline_evidence,'reference_authorities':reference_evidence,'baseline_score_summary':{'objective':data['baseline_score']['objective'],'W9_balanced_error_kcal_mol':data['baseline_score']['energy_controls']['balanced_exchange']['aquo9_central']['error_kcal_mol']},'probe_scores':scores,'actual_native_cost':{'scheduler_seconds':int(scheduler['elapsed_seconds']),'allocated_CPUs':8,'allocated_CPU_seconds':8*int(scheduler['elapsed_seconds']),**{k:r[k] for k in ('elapsed_seconds','parent_cpu_seconds','aggregate_child_user_cpu_seconds','aggregate_child_system_cpu_seconds')},'overlapping_per_program_CPU_not_summed':True},'new_native_calls_this_review':0,'fresh_validation_outcomes_read':False,'candidate_selected':False,'physical_accuracy_qualified':False,'alpha_changes_susceptibility_and_damping_length':True,'augmented_projection_reviewed':False,'historical17670_review_loss_preserved':True,'local_threads':1,'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu}
        rp=save(out/'REVIEW.json',result)
        lines=['# Actual metal-alpha diagnostic review','','All40 new energy/force cells,40 reused entries and8 complete scores replay exactly. Eighty native program receipts pass; no extra native call was made.','','| Probe | Objective | W9 paired error (kcal/mol) | Energy gates | Force gates |','|---|---:|---:|---:|---:|']
        for name,v in scores.items():lines.append(f"| {name} | {v['objective']:.9f} | {v['W9_balanced_error_kcal_mol']:.6f} | {v['energy_checks_passed']}/10 | {v['force_checks_passed']}/15 |")
        lines+=['',f"Actual cost: {scheduler['elapsed_seconds']} scheduler seconds ×8CPUs = {8*int(scheduler['elapsed_seconds'])} allocated CPU-seconds; application {r['elapsed_seconds']:.3f}s. Measured parent CPU {r['parent_cpu_seconds']:.3f}s; aggregate child user/system {r['aggregate_child_user_cpu_seconds']:.3f}/{r['aggregate_child_system_cpu_seconds']:.3f}s. Per-program overlapping counters are not summed.",'','These are fixed responses on consumed fixtures, not a fit or selected candidate. Alpha changes both susceptibility and damping length. The unchanged provisional ECP reference, failed parent scientific gate and17670 lost historical review remain explicit. Fresh validation stays held. Root owns the separate two-scale augmented-response analysis.','',f"Review: {rp['path']}\nSHA256: {rp['sha256']}",'']
        (HERE/'REPORT.md').write_text('\n'.join(lines));print(json.dumps({'review':rp,'report':pin(HERE/'REPORT.md'),'checks':len(checks),'scores':{k:{q:v[q] for q in ('objective','W9_balanced_error_kcal_mol','energy_checks_passed','force_checks_passed','all_controls_passed')} for k,v in scores.items()},'wall_seconds':result['wall_seconds'],'CPU_seconds':result['CPU_seconds']}))
    except Exception as exc:
        save(out/'FAILURE.json',{'status':'actual_review_failed','error':repr(exc),'implementation':pin(__file__),'checks':checks,'new_native_calls':0});raise

if __name__=='__main__':main()
