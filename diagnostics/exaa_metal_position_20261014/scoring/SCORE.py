#!/usr/bin/env python3
"""Eight reciprocal metal-position calls through the unchanged warm CPU worker."""
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
import socket
import subprocess
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/exaa_metal_position_20261014/scoring'
PACK=ROOT/'workspaces/exaa_source_packing_20261014/packing_score_v1'
INVENTORY=ROOT/'workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json'
PROTOCOL='ExaA_reciprocal_metal_position_MACE_vacuum_v1'
ORDER=('origin','adaptive_La');METALS=('Ca','La');EV=23.06054783061903
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for data in iter(lambda:f.read(1048576),b''):h.update(data)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])!=p:raise ValueError('Changed pin: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,data):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')
def load(p):
    s=importlib.util.spec_from_file_location('_pinned_metal_position_worker',verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def prepare(preparation,output):
    """Only metadata copying/hashing; root runs actual preparation first."""
    p=read(preparation);inventory=read(INVENTORY);old=read(PACK/'worker_manifest.json')
    if p['status']!='prepared' or len(p['rows'])!=2:raise ValueError('Both real prepared source rows required')
    rows=[];tasks=[]
    for row in p['rows']:
        source,=[r for r in inventory['rows'] if r['case_id']==row['source_case_id']]
        sample=int(source['case_id'].rsplit('sample-',1)[1])
        if row['status']!='prepared' or row['sample']!=sample or row['role']!=source['role']:raise ValueError('Prepared source identity differs')
        joined={**row,'original':source};rows.append(joined)
        for z in METALS:
            if set(row['matrix'][z])!=set(ORDER):raise ValueError('Wrong fixed candidate set')
            for q in ORDER:
                cell=row['matrix'][z][q]
                tasks.append({'task_id':f'sample{sample}_{z}_{q}_transplanted','case_id':row['case_id'],
                    'sample':sample,'metal':z,'candidate':q,
                    **{k:cell[k] for k in ('xyz','charge','multiplicity','atom_count')}})
    if {r['sample'] for r in rows}!={0,3} or len(tasks)!=8:raise ValueError('Wrong finite source/cell set')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'SCORE.py');shutil.copyfile(HERE/'RUN.sbatch',out/'RUN.sbatch')
    worker=copy.deepcopy(old);worker.update(protocol_id=PROTOCOL,tasks=tasks,cases=[],fixtures=[],origin_task_count=8,
        preparation=pin(preparation),finite_scope={'fixed_geometry_calls':8,'qualification_calls':0,'optimizer_starts':0,'GFN2_calls':0,'DFT_calls':0},
        authorization='Autonomous ExaA reciprocal metal-position diagnostic; root-owned Slurm execution only')
    write(out/'worker_manifest.json',worker)
    c={'protocol_id':PROTOCOL,'implementation':pin(out/'SCORE.py'),'slurm_script':pin(out/'RUN.sbatch'),
       'plan':pin(HERE.parent/'PLAN.md'),'preparation':pin(preparation),'inventory':pin(INVENTORY),
       'template':pin(PACK/'worker_manifest.json'),'worker_manifest':pin(out/'worker_manifest.json'),
       'qualification':pin(PACK/'execution/mace/QUALIFICATION.json'),
       'qualified_worker_summary':pin(PACK/'execution/mace/lane0/CPU_WORKER_SUMMARY.json'),
       'rows':rows,'cells':tasks,'CPUs':7,'new_calls':8,'reused_original_cells':8,
       'new_optimization':False,'energy_conversion_kcal_per_eV':EV,
       'qualification_policy':'Reuse exact previously qualified seven-thread native worker/model; check runtime identity and all actual call receipts; no new fixture calls',
       'output':str(out/'execution')}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def placement():
    job=os.environ.get('SLURM_JOB_ID')
    if not job or not job.isdigit():raise ValueError('Actual Slurm required')
    raw=subprocess.check_output(['scontrol','show','job','-o',job],text=True)
    fields=dict(x.split('=',1) for x in raw.split() if '=' in x)
    if fields['JobState']!='RUNNING' or fields['NumCPUs']!='7' or fields['NumNodes']!='1':raise ValueError('Actual one-node seven-CPU allocation required')
    if int(os.environ['SLURM_CPUS_ON_NODE'])!=7 or len(os.sched_getaffinity(0))<7:raise ValueError('Insufficient locally allocated CPUs')
    if any(os.environ.get(k)!='1' for k in THREADS):raise ValueError('Parent numerical threads must equal one')
    return {'job_id':job,'hostname':socket.gethostname(),'scontrol':raw,'affinity':sorted(os.sched_getaffinity(0)),'allocated_CPUs':7}

def admitted_result(result,request,task,model,versions):
    if result['status']!='complete' or result['model']!=model:raise ValueError('Native result/model unavailable')
    if (request['charge'],request['multiplicity'])!=(task['charge'],task['multiplicity']):raise ValueError('Native request state differs')
    if request['xyz']['sha256']!=task['xyz']['sha256']:raise ValueError('Native request geometry differs')
    verify(request['xyz']);verify(result['forces'])
    s=result['state']
    if (s['charge'],s['spin_multiplicity'],s['atoms'],s['graph_count'],s['head'],s['selected_metal_index'])!=(task['charge'],1,192,1,'omol',0):raise ValueError('Actual OMOL batch state differs')
    if (result['device'],result['threads'],result['actual_model_forwards'],result['parameter_versions_unchanged'])!=('cpu',7,1,True):raise ValueError('Actual native execution differs')
    if result['versions']!=versions or not math.isfinite(result['energy_eV']):raise ValueError('Native versions/energy differ')

def preflight(c,b):
    m=read(c['worker_manifest']);old=read(c['template']);prep=read(c['preparation']);inv=read(c['inventory']);verify(c['plan']);verify(c['slurm_script'])
    if c['protocol_id']!=PROTOCOL or c['new_calls']!=8 or c['reused_original_cells']!=8 or c['energy_conversion_kcal_per_eV']!=EV:raise ValueError('Frozen finite protocol changed')
    for k in ('implementation','model','settings','pool_settings','optimizer_software','source_modules','source_snapshot','calculator_source','model_python','optimizer_python','model_executable','optimizer_executable','gates'):
        if m[k]!=old[k]:raise ValueError('Original native runtime changed: '+k)
    for p in m['source_modules'].values():verify(p)
    for k in ('implementation','calculator_source','model_executable','optimizer_executable'):verify(m[k])
    verify(m['model']['checkpoint']);b.modules(m)
    if m['tasks']!=c['cells'] or len(m['tasks'])!=8 or m['fixtures'] or m['cases'] or m['threads_per_model']!=7 or m['lanes']!=1:raise ValueError('Worker scope changed')
    qualification=read(c['qualification']);summary=read(c['qualified_worker_summary'])
    if not qualification['pass'] or summary['threads']!=7 or summary['manifest']!=c['template']:raise ValueError('Archived runtime qualification unavailable')
    for q in qualification['checks']:
        r=read(q['receipt']);verify(r['forces']);read(r['request'])
        if not q['pass'] or r['model']!=m['model'] or r['versions']!=summary['versions']:raise ValueError('Archived runtime qualification identity differs')
    checks=[];reused={}
    for row in c['rows']:
        source,=[r for r in inv['rows'] if r['case_id']==row['source_case_id']]
        saved,=[r for r in prep['rows'] if r['case_id']==row['case_id']]
        if {k:v for k,v in row.items() if k!='original'}!=saved or row['original']!=source or saved['status']!='prepared':raise ValueError('Preparation/source join differs')
        for k in ('maps','state_signature','source_preparation'):
            if row[k]!=source[k]:raise ValueError('Physical source provenance changed: '+k)
        state=read(row['state_signature']);verify(row['source_preparation'])
        for p in row['maps'].values():verify(p)
        if state['water_inventory'] or state['endpoint_charges']!={'Ca':-3,'La':-2} or state['multiplicities']!={'Ca':1,'La':1}:raise ValueError('Physical state differs')
        metal_points=[]
        for q in ORDER:
            paired={}
            for z in METALS:
                t,=[t for t in c['cells'] if (t['sample'],t['candidate'],t['metal'])==(row['sample'],q,z)]
                original=source['matrix'][z][q];new=b.xyz(verify(t['xyz']));before=b.xyz(verify(original['xyz']))
                if t['xyz']!=row['matrix'][z][q]['xyz'] or len(new)!=192 or t['atom_count']!=192 or new[0][0]!=z or new[1:]!=before[1:]:raise ValueError('Only the physical metal coordinate may change')
                if t['charge']!=state['endpoint_charges'][z] or t['multiplicity']!=1 or not all(math.isfinite(v) for a in new for v in a[1:]):raise ValueError('Endpoint state/coordinates invalid')
                r=read(original['MACE']);request=read(r['request']);admitted_result(r,request,{**t,'xyz':original['xyz']},m['model'],summary['versions'])
                if r['energy_eV']!=original['MACE_eV']:raise ValueError('Archived original energy differs')
                reused[t['task_id']]={'status':'reused','MACE':original['MACE'],'MACE_request':r['request'],'xyz':original['xyz'],'energy_eV':r['energy_eV'],'charge':t['charge'],'multiplicity':1}
                paired[z]=new
            if paired['Ca'][1:]!=paired['La'][1:] or paired['Ca'][0][1:]!=paired['La'][0][1:]:raise ValueError('Paired physical coordinates differ')
            metal_points.append(paired['Ca'][0][1:])
        if metal_points[0]!=metal_points[1]:raise ValueError('Transplanted metal point changed between host candidates')
        checks.append({'sample':row['sample'],'case_id':row['case_id'],'only_metal_changed':True,'paired_coordinates_exact':True,'original_cells_reused':4})
    return {'status':'passed','rows':checks,'reused_originals':reused,'versions':summary['versions'],'qualification':c['qualification'],'new_model_calls':0}

def collect(c,pf,results):
    rows=[]
    for row in c['rows']:
        cells={};works={}
        for z in METALS:
            cells[z]={}
            for q in ORDER:
                t,=[t for t in c['cells'] if (t['sample'],t['metal'],t['candidate'])==(row['sample'],z,q)]
                cells[z][q]={'original':pf['reused_originals'][t['task_id']],
                    'transplanted':results.get(t['task_id'],{'status':'unavailable','reason':'not dispatched','energy_eV':None})}
            works[z]={}
            for scope in ('original','transplanted'):
                o=cells[z]['origin'][scope];p=cells[z]['adaptive_La'][scope]
                works[z][scope]=(p['energy_eV']-o['energy_eV'])*EV if p['energy_eV'] is not None and o['energy_eV'] is not None else None
        differentials={s:works['Ca'][s]-works['La'][s] if all(works[z][s] is not None for z in METALS) else None for s in ('original','transplanted')}
        rows.append({'sample':row['sample'],'role':row['role'],'case_id':row['case_id'],'source_case_id':row['source_case_id'],
          'matrix':cells,'work_kcal_mol':works,'Ca_minus_La_work_kcal_mol':differentials,
          'delta_D_kcal_mol':differentials['transplanted']-differentials['original'] if all(v is not None for v in differentials.values()) else None})
    contrast={}
    for scope in ('original','transplanted'):
        d={r['sample']:r['Ca_minus_La_work_kcal_mol'][scope] for r in rows}
        contrast[scope]=d[3]-d[0] if all(v is not None for v in d.values()) else None
    return {'protocol_id':PROTOCOL,'terminal':True,'preparation':c['preparation'],'inventory':c['inventory'],
       'declared_new_cells':8,'reused_original_cells':8,'complete_new_cells':sum(r['status']=='complete' for r in results.values()),
       'rows':rows,'target_minus_control_differential_work_kcal_mol':contrast,'conversion_kcal_per_eV':EV,
       'scope':'Fixed vacuum MACE response only; no optimization, solvent, DFT, classifier call or calibration'}

def run(config,action):
    placement_record=placement();c=read(config)
    if c['implementation']!=pin(__file__):raise ValueError('Wrong frozen scoring adapter')
    out=Path(c['output'] if action=='run' else str(Path(c['output']).parent/'preflight_'+os.environ['SLURM_JOB_ID']))
    out.mkdir(parents=True,exist_ok=False);began=time.monotonic();transport=None;handlers={};results={};pf=None
    receipt={'status':'failed','config':pin(config),'placement':placement_record,'action':action}
    with (Path(config).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            m=read(c['worker_manifest']);b=load(m['implementation']);pf=preflight(c,b);write(out/'PREFLIGHT.json',pf)
            if action=='run':
                attempt=out/'mace';(attempt/'lane0').mkdir(parents=True);handlers=b.handlers()
                transport=b.WarmCPU(verify(c['worker_manifest']),m,attempt,0)
                if transport.ready['threads']!=7 or transport.ready['versions']!=pf['versions']:raise ValueError('Loaded worker runtime differs')
                for t in c['cells']:
                    request={'manifest':c['worker_manifest'],'task_id':t['task_id'],**{k:t[k] for k in ('xyz','charge','multiplicity')}}
                    r,rp=b.cell(attempt/'cells'/t['task_id'],request,transport,m['model'])
                    result={'status':'unavailable','receipt':rp,'xyz':t['xyz'],'energy_eV':None,'reason':r.get('reason')}
                    try:
                        admitted_result(r,read(r['request']),t,m['model'],pf['versions'])
                        result.update(status='complete',energy_eV=r['energy_eV'],forces=r['forces'])
                    except Exception as e:result['reason']=str(e)
                    results[t['task_id']]=result
                transport.close();transport=None
                receipt['status']='complete' if all(r['status']=='complete' for r in results.values()) else 'partial'
            else:receipt['status']='complete'
        except BaseException as e:
            receipt.update(error=str(e),traceback=traceback.format_exc());raise
        finally:
            if transport:transport.abort()
            for sig,h in handlers.items():signal.signal(sig,h)
            if pf and action=='run':write(out/'FINAL_COLLECTION.json',{'config':pin(config),**collect(c,pf,results)})
            receipt.update(elapsed_seconds=time.monotonic()-began);write(out/'EXECUTION.json',receipt)

def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('prepare');a.add_argument('--preparation',type=Path,required=True);a.add_argument('--output',type=Path,required=True)
    for action in ('preflight','run'):
        a=s.add_parser(action);a.add_argument('--config',type=Path,required=True);a.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.preparation,a.output)
    else:run(verify({'path':str(a.config.resolve()),'sha256':a.config_sha}),a.action)
if __name__=='__main__':main()
