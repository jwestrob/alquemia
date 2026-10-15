#!/usr/bin/env python3
"""Finite24-cell PQQ pose vacuum matrix through the existing warm CPU worker."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import signal
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/pqq_pose'
PROTOCOL='ExaA_intact_PQQ_e3_plus_minus_one_degree_MACE_common9_v1'
BASES=('origin','adaptive_Ca','adaptive_La');SIGNS=('minus','plus');METALS=('Ca','La')
VARIANTS=tuple(f'{q}__PQQ_{sign}1deg' for q in BASES for sign in SIGNS)
ORDER=BASES+VARIANTS;EV=23.06054783061903

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def check(p):
    if pin(p['path'])!=p:raise ValueError('Changed artifact: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def require(value,reason):
    if not value:raise ValueError(reason)

def prepare(preparation,output):
    """Bind only actual admitted preparation; metadata copying/hashing only."""
    pp=pin(preparation);p=read(pp);pm=read(p['manifest']);inv=read(p['inventory'])
    require(p['status']=='prepared' and p['prepared_cells']==24 and len(p['rows'])==2,'Actual complete24-cell preparation required')
    previous=ROOT/'workspaces/exaa_metal_position_20261014/scoring_v1/CONFIG.json';old=read(previous)
    template=read(old['worker_manifest']);cells=[]
    for row in p['rows']:
        require(row['status']=='prepared' and set(row['matrix'])==set(METALS),'Prepared source unavailable')
        src,=[r for r in inv['rows'] if r['case_id']==row['source_case_id']]
        for z in METALS:
            require(set(row['matrix'][z])==set(VARIANTS),'Wrong six pose variants')
            for q in VARIANTS:
                v=row['matrix'][z][q];require(v['status']=='prepared' and v['geometry_checks']['pass'],'Pose geometry unavailable')
                for key in ('xyz','geometry_admission'):check(v[key])
                cells.append({'task_id':f"sample{row['sample']}_{z}_{q}",'sample':row['sample'],'role':row['role'],
                  'source_case_id':row['source_case_id'],'case_id':row['case_id'],'metal':z,'candidate':q,
                  **{k:v[k] for k in ('xyz','charge','multiplicity','atom_count','valence_electrons','base_candidate','sign','geometry_admission')},
                  **{k:src[k] for k in ('maps','state_signature','source_preparation')}})
    require({r['sample'] for r in p['rows']}=={0,3} and len(cells)==24,'Wrong finite matrix')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('SCORE.py','SCORE.sbatch'):shutil.copyfile(HERE/name,out/name)
    worker=copy.deepcopy(template)
    worker.update(protocol_id=PROTOCOL,tasks=cells,cases=[],fixtures=[],origin_task_count=24,preparation=pp,
      finite_scope={'fixed_geometry_calls':24,'qualification_calls':0,'optimizer_starts':0,'SCC_calls':0,'DFT_calls':0},
      authorization='Declared finite intact-PQQ pose diagnostic; actual root-owned Slurm only')
    write(out/'worker_manifest.json',worker)
    c={'protocol':PROTOCOL,'implementation':pin(out/'SCORE.py'),'wrapper':pin(out/'SCORE.sbatch'),
      'plan':pin(HERE.parent/'PQQ_POSE_PLAN.md'),'preparation':pp,'preparation_manifest':p['manifest'],
      'inventory':p['inventory'],'rows':p['rows'],'cells':cells,'worker_manifest':pin(out/'worker_manifest.json'),
      'previous_score_config':pin(previous),'previous_score_adapter':old['implementation'],'worker_template':old['worker_manifest'],
      'python':pin(template['optimizer_python']),
      'original_union_preparation':pm['origin_union_preparation'],'donor_union_preparation':pm['donor_union_preparation'],
      'force_result':pm['force_result'],'original_order':list(BASES),'new_order':list(ORDER),
      'declared_sources':2,'declared_new_cells':24,'reused_original_cells':12,'CPUs':7,
      'energy_conversion_kcal_per_eV':EV,'operational_origin_tolerance_kcal_mol':template['pool_settings']['origin_selection_tolerance_kcal_mol'],
      'solvent_review_absolute_math_pool_change_kcal_mol':1.0,'stage':'MACE_vacuum_only',
      'new_SCC_calls':0,'new_DFT_calls':0,'new_optimizations':0,'output':str(out/'execution')}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def preflight(c,b,s):
    for key in ('plan','wrapper','force_result','original_union_preparation','donor_union_preparation'):check(c[key])
    require(pin(sys.executable)==c['python'],'Parent interpreter differs from qualified numerical environment')
    old=read(c['previous_score_config']);oldpf=s.preflight(old,b)
    worker=read(c['worker_manifest']);template=read(c['worker_template']);p=read(c['preparation']);pm=read(c['preparation_manifest']);inv=read(c['inventory'])
    require(p['manifest']==c['preparation_manifest'] and pm['plan']==c['plan'] and p['status']=='prepared' and p['rows']==c['rows'],'Actual preparation join differs')
    require(c['protocol']==PROTOCOL and c['stage']=='MACE_vacuum_only' and c['declared_new_cells']==24 and c['reused_original_cells']==12,'Protocol scope differs')
    require(c['original_order']==list(BASES) and c['new_order']==list(ORDER) and c['solvent_review_absolute_math_pool_change_kcal_mol']==1.0,'Fixed pool/gate changed')
    require(c['energy_conversion_kcal_per_eV']==EV and c['operational_origin_tolerance_kcal_mol']==0.1,'Selection/units changed')
    for key in ('implementation','model','settings','pool_settings','optimizer_software','source_modules','source_snapshot','calculator_source','model_python','optimizer_python','model_executable','optimizer_executable','gates'):
        require(worker[key]==template[key],'Original warm runtime changed: '+key)
    require(worker['tasks']==c['cells'] and len(c['cells'])==24 and worker['threads_per_model']==7 and worker['lanes']==1 and not worker['fixtures'] and not worker['cases'],'Wrong native call scope')
    require(c['python']==worker['optimizer_executable'],'Qualified parent interpreter pin differs')
    baseline={};checks=[]
    for row in c['rows']:
        sample=row['sample'];src,=[r for r in inv['rows'] if r['case_id']==row['source_case_id']]
        admission=read(row['geometry_admission']);require(admission['pass'] and admission['status']=='prepared' and row['status']=='prepared','Pose admission failed')
        for key in ('maps','state_signature','source_preparation','prepared_admission'):require(row[key]==src[key],'Source provenance differs')
        state=read(row['state_signature']);require(state==src['state'] and not state['water_inventory'],'State changed')
        ids=set(row['PQQ_indices_zero_based'].values());require(len(ids)==27,'Wrong complete cofactor')
        baseline[str(sample)]={z:{} for z in METALS};oldcoords={};paired={}
        for z in METALS:
            for base in BASES:
                v=src['matrix'][z][base];r=read(v['MACE']);request=read(r['request']);check(v['xyz']);check(v['native_frozen_ALPB_result'])
                s.admitted_result(r,request,{**v,'xyz':v['xyz']},worker['model'],oldpf['versions'])
                require(r['energy_eV']==v['MACE_eV'] and request['xyz']['sha256']==v['xyz']['sha256'],'Exact reused energy/coordinate differs')
                baseline[str(sample)][z][base]={'status':'complete','energy_eV':r['energy_eV'],'receipt':v['MACE'],'xyz':v['xyz'],
                  'forces':r['forces'],'request':r['request'],'charge':v['charge'],'multiplicity':v['multiplicity'],
                  'existing_native_ALPB':v['native_frozen_ALPB_result'],'reused':True}
                oldcoords[z,base]=check(v['xyz']).read_text().splitlines(keepends=True)
        for t in (t for t in c['cells'] if t['sample']==sample):
            z,q,base=t['metal'],t['candidate'],t['base_candidate'];v=row['matrix'][z][q]
            require(t['xyz']==v['xyz'] and v['geometry_admission']==row['geometry_admission'] and v['geometry_checks']==admission['variants'][q] and v['geometry_checks']['pass'],'Pose cell gate differs')
            require(t['charge']==state['endpoint_charges'][z] and t['multiplicity']==1 and t['atom_count']==192,'Pose endpoint state differs')
            lines=check(t['xyz']).read_text().splitlines(keepends=True);atoms=b.xyz(check(t['xyz']))
            require(len(atoms)==192 and atoms[0][0]==z and all(math.isfinite(x) for a in atoms for x in a[1:]),'Pose coordinates invalid')
            require(all(lines[i+2]==oldcoords[z,base][i+2] for i in range(192) if i not in ids),'Non-PQQ coordinate line changed')
            require(v['original_xyz']==src['matrix'][z][base]['xyz'],'Wrong original base identity')
            paired[z,q]=atoms
        for q in VARIANTS:
            ca,la=paired['Ca',q],paired['La',q];require(ca[1:]==la[1:] and ca[0][1:]==la[0][1:],'Pose paired coordinates differ')
        checks.append({'sample':sample,'geometry_admission':row['geometry_admission'],'all12_new_cells_admitted':True,'six_original_cells_reused':True})
    return {'status':'passed','versions':oldpf['versions'],'runtime_qualification':oldpf['qualification'],'checks':checks,'baseline':baseline,'new_model_calls':0}

def choose(energies,order,tolerance):
    if any(energies[z].get(q) is None for z in METALS for q in order):return None
    works={z:{q:(energies[z][q]-energies[z]['origin'])*EV for q in order} for z in METALS};out={'works_kcal_mol':works}
    for policy in ('mathematical','operational'):
        selected={}
        for z in METALS:
            best=min(order,key=lambda q:works[z][q])
            selected[z]=best if policy=='mathematical' or works[z][best]<-tolerance else 'origin'
        chosen={z:works[z][selected[z]] for z in METALS}
        out[policy]={'selected_geometries':selected,'per_metal_work_kcal_mol':chosen,'Ca_minus_La_work_kcal_mol':chosen['Ca']-chosen['La']}
    return out

def collect(c,pf,results):
    rows=[]
    for source in c['rows']:
        sample=source['sample'];matrix=copy.deepcopy(pf['baseline'][str(sample)])
        for t in (t for t in c['cells'] if t['sample']==sample):
            matrix[t['metal']][t['candidate']]=results.get(t['task_id'],{'status':'unavailable','energy_eV':None,'reason':'not dispatched','xyz':t['xyz']})
        energy={z:{q:(v['energy_eV'] if v['status']=='complete' else None) for q,v in matrix[z].items()} for z in METALS}
        pools={name:choose(energy,order,c['operational_origin_tolerance_kcal_mol']) for name,order in (('original3',BASES),('expanded9',ORDER))}
        changes={}
        for policy in ('mathematical','operational'):
            old,new=pools['original3'],pools['expanded9'];changes[policy]=None
            if old and new:
                a=old[policy];b=new[policy];delta={z:b['per_metal_work_kcal_mol'][z]-a['per_metal_work_kcal_mol'][z] for z in METALS}
                changes[policy]={'per_metal_change_kcal_mol':delta,'paired_contrast_change_kcal_mol':delta['Ca']-delta['La'],
                  'before':a['selected_geometries'],'after':b['selected_geometries']}
        allworks={z:{q:((v-energy[z]['origin'])*EV if v is not None else None) for q,v in energy[z].items()} for z in METALS}
        finite={}
        for q in VARIANTS:
            base=q.split('__PQQ_',1)[0]
            work={z:(allworks[z][q]-allworks[z][base] if allworks[z][q] is not None else None) for z in METALS}
            finite[q]={'base':base,'per_metal_pose_work_kcal_mol':work,'Ca_minus_La_pose_work_kcal_mol':work['Ca']-work['La'] if all(v is not None for v in work.values()) else None}
        rows.append({'sample':sample,'role':source['role'],'source_case_id':source['source_case_id'],'matrix':matrix,
          'all_origin_relative_works_kcal_mol':allworks,'finite_pose_responses':finite,'pools':pools,'common3_to9':changes,'classification':None})
    valid=[r for r in rows if r['common3_to9']['mathematical'] is not None]
    qualifies=any(abs(r['common3_to9']['mathematical']['paired_contrast_change_kcal_mol'])>=1.0 for r in valid)
    decision='qualifies_for_root_solvent_review' if qualifies else ('close_finite_angle_test' if len(valid)==2 else 'unresolved_incomplete_matrix')
    return {'protocol':PROTOCOL,'terminal':True,'preparation':c['preparation'],'inventory':c['inventory'],
      'declared_sources':2,'declared_new_cells':24,'reused_original_cells':12,'complete_new_cells':sum(r['status']=='complete' for r in results.values()),
      'rows':rows,'pool_orders':{'original3':list(BASES),'expanded9':list(ORDER)},
      'stage_gate':{'criterion':'absolute mathematical common3_to9 paired contrast change >=1kcal/mol in either source; either sign',
        'threshold_kcal_mol':1.0,'decision':decision,'solvent_executed':False,'automatic_solvent_launch':False},
      'conversion_kcal_per_eV':EV,'new_SCC_calls':0,'new_DFT_calls':0,'new_optimizations':0,'production_changed':False,
      'limitations':'MACE vacuum finite-pose response, developmental consumed sources; no composite gradient, calibration, affinity, DFT PQQ validation or class calls'}

def run(config,action):
    require(bool(os.environ.get('SLURM_JOB_ID')),'Actual Slurm required before scientific imports')
    c=read(config);require(c['implementation']==pin(__file__),'Wrong scorer snapshot')
    s=load('_pqq_pose_prior_score',c['previous_score_adapter']);placement=s.placement()
    out=Path(c['output'] if action=='run' else str(Path(c['output']).parent/'preflight_'+os.environ['SLURM_JOB_ID']))
    out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time();transport=None;handlers={};results={};pf=None
    receipt={'status':'failed','config':pin(config),'placement':placement,'action':action}
    with (Path(config).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            m=read(c['worker_manifest']);b=load('_pqq_pose_warm_worker',m['implementation']);pf=preflight(c,b,s);write(out/'PREFLIGHT.json',pf)
            if action=='run':
                attempt=out/'mace';(attempt/'lane0').mkdir(parents=True);handlers=b.handlers();transport=b.WarmCPU(check(c['worker_manifest']),m,attempt,0)
                require(transport.ready['threads']==7 and transport.ready['versions']==pf['versions'],'Loaded runtime differs')
                for t in c['cells']:
                    request={'manifest':c['worker_manifest'],'task_id':t['task_id'],**{k:t[k] for k in ('xyz','charge','multiplicity')}}
                    value,rp=b.cell(attempt/'cells'/t['task_id'],request,transport,m['model'])
                    result={'status':'unavailable','receipt':rp,'xyz':t['xyz'],'energy_eV':None,'reason':value.get('reason')}
                    try:
                        s.admitted_result(value,read(value['request']),t,m['model'],pf['versions'])
                        result.update(status='complete',energy_eV=value['energy_eV'],forces=value['forces'],request=value['request'])
                    except Exception as e:result['reason']=str(e)
                    results[t['task_id']]=result
                transport.close();transport=None
                receipt['status']='complete' if len(results)==24 and all(r['status']=='complete' for r in results.values()) else 'partial'
            else:receipt['status']='complete'
        except BaseException as e:receipt.update(error=str(e),traceback=traceback.format_exc());raise
        finally:
            try:
                if transport:transport.abort()
                for sig,handler in handlers.items():signal.signal(sig,handler)
                write(out/'PARTIAL_RESULTS.json',{'config':pin(config),'rows':results,'declared_cells':24})
                if pf and action=='run':write(out/'FINAL_COLLECTION.json',{'config':pin(config),**collect(c,pf,results)})
            except BaseException as e:
                receipt.update(status='failed',reporting_error={'error':str(e),'traceback':traceback.format_exc()});raise
            finally:
                elapsed=time.monotonic()-start
                receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,parent_process_CPU_seconds=time.process_time()-cpu)
                try:
                    receipt.update(actual_model_forwards=sum(read(p).get('actual_model_forwards',0) for p in out.glob('mace/cells/*/mace_result.json')),
                      worker_resource_receipts=[pin(p) for p in out.glob('mace/lane*/CPU_WORKER_SUMMARY.json')])
                except BaseException as e:receipt.update(status='failed',resource_reporting_error=str(e))
                write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--preparation',required=True);a.add_argument('--output',required=True)
    for action in ('preflight','run'):
        a=sub.add_parser(action);a.add_argument('--config',required=True);a.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.preparation,a.output)
    else:run(check({'path':str(Path(a.config).resolve()),'sha256':a.config_sha}),a.action)
