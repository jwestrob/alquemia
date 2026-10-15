#!/usr/bin/env python3
"""Slurm-only saved six-cell vacuum comparison; no molecule calls or fitted bands."""
import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import traceback

HA_TO_KCAL=627.509474
EV_TO_KCAL=23.06054783061903
ORDER=('origin','adaptive_Ca','adaptive_La')

def pin(path):
    path=Path(path).resolve();h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(path),'sha256':h.hexdigest()}

def check(r):
    if pin(r['path'])['sha256']!=r['sha256']:raise ValueError('Changed artifact: '+r['path'])
    return Path(r['path'])

def read(r):return json.loads((check(r) if isinstance(r,dict) else Path(r)).read_text())
def write(p,data):
    with p.open('x') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')

def paired(work,selection):
    values={z:work[z].get(selection[z]) for z in ('Ca','La')}
    available=all(v is not None for v in values.values())
    return {'status':'available_diagnostic' if available else 'unavailable','selected_geometries':selection,
            'works_kcal_mol':values,'Ca_minus_La_kcal_mol':values['Ca']-values['La'] if available else None}

def pool(work):
    if any(work[z].get(q) is None for z in ('Ca','La') for q in ORDER):
        return {'status':'unavailable','selected_geometries':None,'works_kcal_mol':None,'Ca_minus_La_kcal_mol':None,
                'reason':'Complete six-cell common pool required; no partial minima'}
    selection={z:min(ORDER,key=lambda q:work[z][q]) for z in ('Ca','La')}
    return paired(work,selection)

def compare(mp,m,collections):
    if pin(__file__)!=m['implementation']['COMPARE.py']['snapshot']:raise ValueError('Wrong comparison snapshot')
    inv=read(m['source_collection']);source,=m['sources'];saved,=[r for r in inv['rows'] if r['case_id']==source['case_id']]
    known={t['task_id']:t for t in m['tasks']}
    if len(known)!=6 or source['sample']!=0 or source['role']!='target' or m['candidate_order']!=list(ORDER):raise ValueError('Wrong frozen six-cell target scope')
    observed={};receipts=[]
    for path in collections:
        if not path.exists():receipts.append({'path':str(path.resolve()),'status':'missing_collection'});continue
        c=read(path);receipts.append(pin(path))
        if c['manifest']!=pin(mp) or c['protocol_id']!=m['protocol_id']:raise ValueError('Incompatible native collection')
        for r in c['cells']:
            tid=r['task_id']
            if tid not in known:raise ValueError('Undeclared collected cell')
            if r['status']=='not_completed':continue
            if tid in observed and observed[tid]!=r:raise ValueError('Conflicting collected cell: '+tid)
            if r['status']=='complete_diagnostic':
                if not all(r['checks'].values()):raise ValueError('Admitted cell has failed checks')
                check(r['output']);execution=read(r['execution']);t=known[tid]
                if execution['manifest']!=pin(mp) or execution['task_id']!=tid or execution['returncode']!=0:raise ValueError('Native receipt mismatch')
                if execution['artifacts']['xyz']['sha256']!=t['xyz']['sha256']:raise ValueError('Native coordinate mismatch')
            observed[tid]=r
    cells=[];energies={'DFT':{z:{} for z in ('Ca','La')},'MACE':{z:{} for z in ('Ca','La')}};cost=[]
    for tid,t in known.items():
        row={**observed.get(tid,{'task_id':tid,'status':'not_completed'}),'metal':t['metal'],'candidate':t['candidate']}
        cells.append(row);z=t['metal'];q=t['candidate'];cell=saved['matrix'][z][q]
        mr=read(t['saved_MACE']);request=read(t['saved_MACE_request']);check(t['source_xyz']);check(t['xyz'])
        if cell['MACE']!=t['saved_MACE'] or cell['xyz']!=t['source_xyz'] or t['xyz']['sha256']!=t['source_xyz']['sha256']:raise ValueError('Saved MACE/source binding changed')
        if mr['request']!=t['saved_MACE_request'] or request['xyz']!=t['source_xyz'] or request['charge']!=t['charge'] or request['multiplicity']!=1:raise ValueError('Saved MACE request state mismatch')
        if mr['status']!='complete' or mr['energy_eV']!=cell['MACE_eV'] or mr['energy_eV']!=t['saved_MACE_eV']:raise ValueError('Saved MACE energy mismatch')
        energies['MACE'][z][q]=mr['energy_eV'];energies['DFT'][z][q]=row.get('energy_hartree') if row['status']=='complete_diagnostic' else None
        if energies['DFT'][z][q] is not None and not math.isfinite(energies['DFT'][z][q]):raise ValueError('Nonfinite admitted energy')
        if 'execution' in row:
            receipt=read(row['execution']);start=dt.datetime.fromisoformat(receipt['started_at_utc']);end=dt.datetime.fromisoformat(receipt['finished_at_utc'])
            cost.append({'task_id':tid,'native_receipt':row['execution'],'native_wall_seconds':(end-start).total_seconds(),
                         'MPI_ranks':receipt['parallelism']['nprocs'],'allocation':receipt['allocation'],
                         'scope':'Native receipt interval only; excludes preflight/collection/queue, not end-to-end latency'})
    works={}
    for model,factor in (('DFT',HA_TO_KCAL),('MACE',EV_TO_KCAL)):
        works[model]={z:{q:(energies[model][z][q]-energies[model][z]['origin'])*factor
            if energies[model][z][q] is not None and energies[model][z]['origin'] is not None else None for q in ORDER} for z in ('Ca','La')}
    paths={}
    for q in ORDER:
        paths[q]={model:paired(works[model],{'Ca':q,'La':q}) for model in works}
        paths[q]['DFT_minus_MACE_kcal_mol']=paths[q]['DFT']['Ca_minus_La_kcal_mol']-paths[q]['MACE']['Ca_minus_La_kcal_mol'] if paths[q]['DFT']['status']=='available_diagnostic' else None
    pools={model:pool(work) for model,work in works.items()}
    pools['DFT_minus_MACE_kcal_mol']=pools['DFT']['Ca_minus_La_kcal_mol']-pools['MACE']['Ca_minus_La_kcal_mol'] if pools['DFT']['status']=='available_diagnostic' else None
    chosen=m['historical_operational_selection']
    if chosen!=saved['stored_score']['geometry_choices'] or any(q not in ORDER for q in chosen.values()):raise ValueError('Historical selections changed')
    historical={model:paired(work,chosen) for model,work in works.items()}
    historical['saved_composite_result']=m['historical_score_result'];check(m['historical_score_result'])
    historical['saved_composite_accommodation_kcal_mol']=saved['stored_score']['accommodation_R_shift']
    historical['definition']='Vacuum work evaluated at historical composite-selected geometries; no new operational selection or class call'
    result={'schema':'mmol_common3_vacuum_comparison.v1','manifest':pin(mp),'inventory':m['source_collection'],'collections':receipts,
        'implementation':pin(__file__),'source_case_id':source['source_case_id'],'source1':'inventory only; not evaluated',
        'declared_cells':6,'available_cells':sum(c['status']=='complete_diagnostic' for c in cells),'cells':cells,
        'energies':energies,'per_metal_works_kcal_mol':works,'common_path_responses':paths,'complete_vacuum_pool':pools,
        'historical_composite_selected_geometries':historical,'actual_native_receipt_costs':cost,
        'units':{'DFT_energy':'hartree','MACE_energy':'eV','work':'kcal/mol','HA_TO_KCAL':HA_TO_KCAL,'EV_TO_KCAL':EV_TO_KCAL},
        'formulas':{'work':'W_M(q)=E_M(q)-E_M(origin)','path':'D(q)=W_Ca(q)-W_La(q)',
                    'vacuum_pool':'min_q W_Ca(q)-min_q W_La(q); exact tie order origin,adaptive_Ca,adaptive_La'},
        'limitations':'Fixed saved common pool, not DFT optimization; electronic-state uniqueness/forces unqualified; no solvent, bands, affinity or class decision',
        'slurm_job_id':os.environ['SLURM_JOB_ID']}
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('manifest',type=Path);p.add_argument('--manifest-sha',required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--collections',nargs='+',type=Path,required=True);a=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Scientific saved-artifact comparison requires actual Slurm')
    mp=check({'path':str(a.manifest.resolve()),'sha256':a.manifest_sha});m=read(mp);a.output.mkdir(parents=True,exist_ok=False)
    status={'status':'started','manifest':pin(mp),'started_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'slurm_job_id':os.environ['SLURM_JOB_ID']}
    try:
        result=compare(mp,m,a.collections);write(a.output/'RESULT.json',result)
        status.update(status='complete',result=pin(a.output/'RESULT.json'))
        print(json.dumps({'available_cells':result['available_cells'],'declared_cells':6,'common_paths':result['common_path_responses'],'vacuum_pool':result['complete_vacuum_pool']}))
    except BaseException as e:
        status.update(status='failed',error=f'{type(e).__name__}: {e}',traceback=traceback.format_exc());raise
    finally:
        status['finished_utc']=dt.datetime.now(dt.timezone.utc).isoformat();write(a.output/'EXECUTION.json',status)

if __name__=='__main__':main()
