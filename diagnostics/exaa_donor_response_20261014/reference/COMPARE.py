#!/usr/bin/env python3
"""Join actual completed native collections to the frozen common-path inventory."""
import argparse
import hashlib
import json
import os
from pathlib import Path

HA_TO_KCAL = 627.509474
EV_TO_KCAL = 23.06054783061903

def pin(p):
    p=Path(p).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}

def check(r):
    if pin(r['path'])['sha256']!=r['sha256']:raise ValueError('Hash mismatch: '+r['path'])
    return Path(r['path'])

def main():
    p=argparse.ArgumentParser();p.add_argument('manifest',type=Path);p.add_argument('--manifest-sha',required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--collections',nargs='+',type=Path,required=True)
    a=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required for scientific comparison')
    mp={'path':str(a.manifest.resolve()),'sha256':a.manifest_sha};m=json.loads(check(mp).read_text())
    inventory=json.loads(check(m['source_collection']).read_text());known={t['task_id']:t for t in m['tasks']}
    if len(known)!=8:raise ValueError('Expected eight-cell declared denominator')
    observed={};receipts=[]
    for path in a.collections:
        if not path.exists():
            receipts.append({'path':str(path.resolve()),'status':'missing_collection'});continue
        c=json.loads(path.read_text());receipts.append(pin(path))
        if c['manifest']!=mp:raise ValueError('Incompatible native collection')
        for r in c['cells']:
            tid=r['task_id']
            if tid not in known:raise ValueError('Undeclared collected cell')
            if r['status']=='not_completed':continue
            if tid in observed and observed[tid]!=r:raise ValueError('Conflicting collected cell: '+tid)
            if r['status']=='complete_diagnostic':
                if not all(r['checks'].values()):raise ValueError('Admitted cell has failed check')
                check(r['output']);check(r['execution'])
            observed[tid]=r
    cells=[];lookup={}
    for tid,t in known.items():
        r=observed.get(tid,{'task_id':tid,'status':'not_completed'})
        r={**r,'sample':t['sample'],'metal':t['metal'],'candidate':t['candidate']}
        cells.append(r);lookup[t['sample'],t['metal'],t['candidate']]=r
    results=[]
    for s in m['sources']:
        saved,=[r for r in inventory['rows'] if r['case_id']==s['source_case_id']]
        row={'sample':s['sample'],'role':s['role'],'case_id':s['case_id'],'works':{}}
        for metal in ('Ca','La'):
            origin=lookup[s['sample'],metal,'origin'];moved=lookup[s['sample'],metal,'adaptive_La']
            mace=(saved['matrix'][metal]['adaptive_La']['MACE_eV']-saved['matrix'][metal]['origin']['MACE_eV'])*EV_TO_KCAL
            work={'MACE_vacuum_kcal_mol':mace,'DFT_vacuum_kcal_mol':None,'DFT_minus_MACE_kcal_mol':None,
                  'status':'unavailable','required_cells':[origin['task_id'],moved['task_id']]}
            if origin['status']==moved['status']=='complete_diagnostic':
                dft=(moved['energy_hartree']-origin['energy_hartree'])*HA_TO_KCAL
                work.update(status='available_diagnostic',DFT_vacuum_kcal_mol=dft,DFT_minus_MACE_kcal_mol=dft-mace)
            row['works'][metal]=work
        ca,la=row['works']['Ca'],row['works']['La']
        row['differential_work']={'MACE_vacuum_kcal_mol':ca['MACE_vacuum_kcal_mol']-la['MACE_vacuum_kcal_mol'],
                                  'DFT_vacuum_kcal_mol':None,'DFT_minus_MACE_kcal_mol':None,'status':'unavailable'}
        if ca['status']==la['status']=='available_diagnostic':
            d=ca['DFT_vacuum_kcal_mol']-la['DFT_vacuum_kcal_mol']
            row['differential_work'].update(status='available_diagnostic',DFT_vacuum_kcal_mol=d,
                DFT_minus_MACE_kcal_mol=d-row['differential_work']['MACE_vacuum_kcal_mol'])
        results.append(row)
    control=next(r for r in results if r['role']=='control')['differential_work']
    target=next(r for r in results if r['role']=='target')['differential_work']
    contrast={'definition':'D_target-D_control','status':'unavailable',
              'MACE_vacuum_kcal_mol':target['MACE_vacuum_kcal_mol']-control['MACE_vacuum_kcal_mol'],
              'DFT_vacuum_kcal_mol':None,'DFT_minus_MACE_kcal_mol':None}
    if target['status']==control['status']=='available_diagnostic':
        d=target['DFT_vacuum_kcal_mol']-control['DFT_vacuum_kcal_mol']
        contrast.update(status='available_diagnostic',DFT_vacuum_kcal_mol=d,DFT_minus_MACE_kcal_mol=d-contrast['MACE_vacuum_kcal_mol'])
    result={'schema':'exaa_common_path_response_comparison.v1','manifest':mp,'collections':receipts,
            'implementation':pin(__file__),'cells':cells,'sources':results,'declared_cells':8,
            'target_minus_control_differential':contrast,
            'available_cells':sum(r['status']=='complete_diagnostic' for r in cells),
            'units':{'input_DFT':'hartree','input_MACE':'eV','works':'kcal/mol','HA_TO_KCAL':HA_TO_KCAL,'EV_TO_KCAL':EV_TO_KCAL},
            'formula':'W_M=E_M(adaptive_La)-E_M(origin); D=W_Ca-W_La',
            'limitations':'Common physical path, not operational candidate minima; reference state stability unestablished; no class decision or affinity',
            'slurm_job_id':os.environ['SLURM_JOB_ID']}
    a.output.mkdir(parents=True,exist_ok=False)
    with (a.output/'RESULT.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'available_cells':result['available_cells'],'declared_cells':8,'sources':results}))

if __name__=='__main__':main()
