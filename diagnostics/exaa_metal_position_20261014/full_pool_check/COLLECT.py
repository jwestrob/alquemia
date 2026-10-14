#!/usr/bin/env python3
"""Saved-only historical-three versus common-five geometry pool comparison."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE=ROOT/'diagnostics/exaa_metal_position_20261014/full_pool_check'
BASE=ROOT/'workspaces/exaa_metal_position_20261014'
OLD=('origin','adaptive_Ca','adaptive_La');NEW=OLD+('transplanted_origin','transplanted_adaptive_La')
POLICIES=('operational','mathematical','origin');METALS=('Ca','La');TOL=1e-8
PROTOCOL='ExaA_saved_historical3_to_common5_pool_v1'

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])!=p:raise ValueError('Changed pin: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def prepare(output):
    ip=ROOT/'workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json'
    sp=BASE/'solvent_v1/execution/FINAL_COLLECTION.json';cp=sp.parent.parent/'CONFIG.json'
    bp=ROOT/'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/FINAL_COLLECTION.json'
    inventory=read(ip);solvent=read(sp);cfg=read(cp);baseline=read(bp);rows=[]
    if solvent['config']!=pin(cp) or solvent['available_new_cells']!=8:raise ValueError('Completed reciprocal solvent inputs required')
    for source in inventory['rows']:
        sample=int(source['case_id'].rsplit('sample-',1)[1]);new,=[r for r in solvent['rows'] if r['source_case_id']==source['case_id']]
        historical,=[r for r in baseline['rows'] if r['record_id']==source['record_id']]
        if sample not in (0,3):raise ValueError('Unexpected source')
        cells=[]
        for z in METALS:
            for q in NEW:
                if q in OLD:saved=source['matrix'][z][q]
                else:
                    v=new['matrix']['transplanted'][z][q.removeprefix('transplanted_')]
                    if v['status']!='available':raise ValueError('Missing completed transplanted cell')
                    saved={'xyz':v['xyz'],'charge':source['state']['endpoint_charges'][z],'multiplicity':1,'atom_count':192,
                       'MACE':v['MACE'],'MACE_eV':v['MACE_eV'],'native_frozen_ALPB_result':v['native_result'],'native_frozen_ALPB_hartree':v['G_qvac_hartree']}
                cells.append({'candidate':q,'metal':z,'original':saved,**{k:source[k] for k in ('state_signature','maps','source_preparation')}})
        rows.append({'sample':sample,'role':source['role'],'source_case_id':source['case_id'],'record_id':source['record_id'],
          'cells':cells,'historical_saved_pool':historical['pool']})
    if len(rows)!=2 or sum(len(r['cells']) for r in rows)!=20:raise ValueError('Wrong finite20 saved-cell scope')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'COLLECT.py');shutil.copyfile(HERE/'RUN.sbatch',out/'RUN.sbatch')
    c={'protocol_id':PROTOCOL,'implementation':pin(out/'COLLECT.py'),'wrapper':pin(out/'RUN.sbatch'),'plan':pin(HERE/'PLAN.md'),
       'inventory':pin(ip),'solvent_collection':pin(sp),'solvent_config':pin(cp),'solvent_adapter':cfg['implementation'],
       'historical_collection':pin(bp),'scorer':baseline['scorer'],'rows':rows,'historical_order':OLD,'extended_order':NEW,
       'policies':POLICIES,'algebra_tolerance_kcal_mol':TOL,'new_molecular_calls':0,'classification':None}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def choose(selector,matrix,order):
    previous=selector.ORDER
    try:
        # Private module in this one-process collector; source files stay immutable.
        selector.ORDER=tuple(order)
        return selector.choose({z:{q:matrix[z][q] for q in order} for z in METALS})
    finally:selector.ORDER=previous

def close(a,b,residuals):
    if not math.isfinite(a) or not math.isfinite(b):raise ValueError('Nonfinite saved algebra')
    d=abs(a-b);residuals.append(d)
    if d>TOL:raise ValueError('Saved historical/algebra disagreement: '+str(d))

def reproduce(actual,saved,residuals):
    if actual['status']!='available':return {'status':'unavailable','reason':'required original cell failed validation'}
    if saved['status']!='available':raise ValueError('Saved historical pool is unavailable')
    for policy in POLICIES:
        a,b=actual[policy],saved[policy];close(a['R_model_kcal_mol'],b['R_model_kcal_mol'],residuals)
        for z in METALS:
            if a[z]['candidate']!=b[z]['candidate']:raise ValueError('Historical selected geometry differs')
            close(a[z]['work_kcal_mol'],b[z]['work_kcal_mol'],residuals)
            for q in OLD:close(a[z]['all_works_kcal_mol'][q],b[z]['all_works_kcal_mol'][q],residuals)
            for k in ('MACE_eV','solvent_hartree'):close(a[z]['components'][k],b[z]['components'][k],residuals)
            for k in ('xyz','MACE','native_result'):
                if a[z]['components'][k]!=b[z]['components'][k]:raise ValueError('Historical selected evidence differs: '+k)
    return {'status':'passed','policies':list(POLICIES)}

def execute(path):
    job=os.environ.get('SLURM_JOB_ID')
    if not job:raise ValueError('Actual Slurm required before scientific parsing/collection')
    raw=subprocess.check_output(['scontrol','show','job','-o',job],text=True);fields=dict(x.split('=',1) for x in raw.split() if '=' in x)
    if fields['JobState']!='RUNNING' or fields['NumCPUs']!='1' or fields['NumNodes']!='1':raise ValueError('One actual Slurm CPU required')
    c=read(path)
    if c['implementation']!=pin(__file__) or c['protocol_id']!=PROTOCOL or tuple(c['historical_order'])!=OLD or tuple(c['extended_order'])!=NEW:raise ValueError('Frozen collector/order differs')
    out=Path(path).parent/'execution';out.mkdir(exist_ok=False);started=time.monotonic();receipt={'status':'failed','config':pin(path),'job_id':job,'hostname':socket.gethostname(),'scontrol':raw,'new_molecular_calls':0}
    try:
        for k in ('plan','wrapper','inventory','solvent_collection','historical_collection'):verify(c[k])
        cfg=read(c['solvent_config']);adapter=load('_saved_solvent_validation',c['solvent_adapter'])
        sys.path.insert(0,str(ROOT/'scripts'))
        n=load('_saved_native_component_validation',cfg['native_adapter']);r=load('_saved_component_v2_validation',cfg['component'])
        selector=load('_saved_frozen_distribution_policy',c['scorer'])
        if tuple(selector.ORDER)!=OLD or tuple(selector.POLICIES)!=POLICIES or selector.EV!=adapter.EV or selector.HA!=adapter.HA:raise ValueError('Historical selection/energy policy differs')
        if c['algebra_tolerance_kcal_mol']!=TOL or c['new_molecular_calls']!=0:raise ValueError('Frozen scope/tolerance differs')
        rows=[];residuals=[]
        for source in c['rows']:
            matrix={z:{q:None for q in NEW} for z in METALS};evidence={z:{} for z in METALS};failures=[]
            for cell in source['cells']:
                z,q=cell['metal'],cell['candidate']
                try:
                    adapter.mace(cell['original'],cfg)
                    v=adapter.reusable(cfg,cell,n,r);evidence[z][q]=v
                    matrix[z][q]={'MACE_eV':v['MACE_eV'],'solvent_hartree':v['G_qvac_hartree'],
                       'MACE':v['MACE'],'native_result':v['native_result'],'xyz':v['xyz']}
                except Exception as e:failures.append({'metal':z,'candidate':q,'reason':type(e).__name__+': '+str(e)})
            historical=choose(selector,matrix,OLD);extended=choose(selector,matrix,NEW)
            regression=reproduce(historical,source['historical_saved_pool'],residuals)
            changes={};works={z:{} for z in METALS}
            for z in METALS:
                for q in NEW:works[z][q]=adapter.work(evidence[z].get(q,{'status':'unavailable'}),evidence[z].get('origin',{'status':'unavailable'}))
            for policy in POLICIES:
                a,b=historical[policy],extended[policy]
                if a is None or b is None:changes[policy]=None;continue
                permetal={z:adapter.difference(works[z][b[z]['candidate']],works[z][a[z]['candidate']]) for z in METALS}
                differential=adapter.difference(permetal['Ca'],permetal['La'])
                delta=b['R_model_kcal_mol']-a['R_model_kcal_mol'];close(delta,differential['total'],residuals)
                changes[policy]={'selected':{z:{'historical':a[z]['candidate'],'extended':b[z]['candidate'],
                    'historical_xyz':a[z]['components']['xyz'],'extended_xyz':b[z]['components']['xyz']} for z in METALS},
                    'per_metal_change_kcal_mol':permetal,'differential_change_kcal_mol':differential,'R_model_change_kcal_mol':delta}
            rows.append({'sample':source['sample'],'role':source['role'],'source_case_id':source['source_case_id'],'record_id':source['record_id'],
                'historical3':historical,'extended5':extended,'historical_reproduction':regression,'changes':changes,
                'all_origin_relative_works_kcal_mol':works,'evidence':evidence,'failed_cells':failures,'classification':None})
        result={'protocol_id':PROTOCOL,'config':pin(path),'terminal':True,'declared_sources':2,'declared_saved_cells':20,
          'available_cells':sum(len(r['evidence'][z]) for r in rows for z in METALS),'rows':rows,
          'actual_saved_algebra_checks':{'tolerance_kcal_mol':TOL,'count':len(residuals),'maximum_residual':max(residuals,default=None)},
          'historical_order':OLD,'extended_order':NEW,'new_molecular_calls':0,'classification':None,
          'interpretation':'Within-source incremental pool utility only; reciprocal donor response remains a separate result; no affinity/class claim'}
        write(out/'FINAL_COLLECTION.json',result);receipt.update(status='complete',result=pin(out/'FINAL_COLLECTION.json'))
    except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
    finally:receipt['elapsed_seconds']=time.monotonic()-started;write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('prepare');a.add_argument('--output',type=Path,required=True)
    a=s.add_parser('execute');a.add_argument('--config',type=Path,required=True);a.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.output)
    else:execute(verify({'path':str(a.config.resolve()),'sha256':a.config_sha}))
