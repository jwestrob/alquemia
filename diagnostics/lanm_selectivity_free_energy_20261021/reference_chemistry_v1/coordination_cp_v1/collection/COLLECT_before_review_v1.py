#!/usr/bin/env python3
"""Saved eight-fragment coordination CP cycle; no executable is invoked."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import time

HERE=Path(__file__).resolve().parent
GIDS=tuple(f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy'))
ROLES=('ion_ghost_shell','shell_ghost_ion')
IDS=tuple(f'pbe0__aquo{n}__{m}__{role}' for n in (8,9) for m in ('La','Dy') for role in ROLES)
SIGNS={'aquo8_central_La':1,'aquo8_central_Dy':-1,'aquo9_central_La':-1,'aquo9_central_Dy':1}
HA=627.509474

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def same(a,b,label):require(json.loads(json.dumps(a))==json.loads(json.dumps(b)),label)
def finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)
def save(path,x):Path(path).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def delta(values):return math.fsum(SIGNS[g]*values[g] for g in GIDS)

def calculate(complexes,fragments,model,cancellations):
    require(set(complexes)==set(model)==set(GIDS) and set(fragments)==set(IDS),'Incomplete cycle denominator')
    require(cancellations and all(v is True for v in cancellations.values()),'Own-fragment cancellations unavailable')
    for value in complexes.values():
        require(all(finite(value[k]) for k in ('electronic_energy_hartree','D4_energy_hartree','final_energy_hartree')),'Nonfinite original energy')
        require(abs(value['electronic_energy_hartree']+value['D4_energy_hartree']-value['final_energy_hartree'])<1e-8,'Original electronic/D4 closure differs')
    require(all(finite(v) for v in model.values()),'Nonfinite model energy')
    for v in fragments.values():
        require(v.get('status')=='complete' and finite(v.get('electronic_energy_hartree')) and v.get('D4_energy_hartree')==0.,'Fragment unavailable or contains dispersion')
    G={g:math.fsum(fragments[f'pbe0__aquo{g[4]}__{g.rsplit("_",1)[1]}__{role}']['electronic_energy_hartree'] for role in ROLES) for g in GIDS}
    electronic=delta({g:v['electronic_energy_hartree'] for g,v in complexes.items()});d4=delta({g:v['D4_energy_hartree'] for g,v in complexes.items()})
    total=delta({g:v['final_energy_hartree'] for g,v in complexes.items()});ghost=delta(G);cp_total=total-ghost;cp_electronic=electronic-ghost;model_delta=delta(model)
    require(abs(total-electronic-d4)<4e-8 and abs(cp_total-cp_electronic-d4)<4e-8,'Paired dispersion closure differs')
    return {'delta_convention':'Dy9 - La9 - Dy8 + La8','hartree_to_kcal_mol':HA,'G_electronic_hartree':G,
        'hartree':{'original_electronic':electronic,'original_D4_unchanged':d4,'original_total':total,'delta_G_electronic':ghost,'CP_shift':-ghost,'CP_electronic':cp_electronic,'CP_total':cp_total},
        'kcal_mol':{'original_electronic':electronic*HA,'original_D4_unchanged':d4*HA,'original_total':total*HA,'delta_G_electronic':ghost*HA,'CP_shift':-ghost*HA,
            'CP_electronic':cp_electronic*HA,'CP_total':cp_total*HA,'unchanged_model':model_delta,'model_minus_original':model_delta-total*HA,'model_minus_CP':model_delta-cp_total*HA},
        'own_fragment_terms_cancel_algebraically':True,'individual_absolute_CP_complex_energies_computed':False,'fragment_D4_added':False,'fit_targets_changed':False}

def cancellation_check(panel,data):
    require(set(panel['saved_complexes'])==set(panel['geometries'])==set(GIDS),'Original complex inventory differs')
    checks={};inputs={g:verify(panel['saved_complexes'][g]['input']).read_text() for g in GIDS}
    for metal in ('La','Dy'):
        a,b=(panel['geometries'][f'aquo{n}_central_{metal}'] for n in (8,9))
        checks['own_'+metal+'_same_recipe_basis_ECP']=inputs['aquo8_central_'+metal]==inputs['aquo9_central_'+metal]
        checks['own_'+metal+'_same_position_and_parent_charge']=a['positions_A'][0]==b['positions_A'][0]==[0.,0.,0.] and a['charge']==b['charge']==3
        ions=[c for c in panel['cells'] if c['metal']==metal and c['role']==ROLES[0]]
        checks['own_'+metal+'_same_effective_state']=len(ions)==2 and all(c['charge']==3 and c['multiplicity']==1 and c['expected_electrons']==8 and c['expected_ecp_core']==(46 if metal=='La' else 55) for c in ions)
    for n in (8,9):
        a,b=(panel['geometries'][f'aquo{n}_central_{metal}'] for metal in ('La','Dy'))
        checks['own_shell'+str(n)+'_same_coordinates_and_symbols']=a['positions_A'][1:]==b['positions_A'][1:] and a['symbols'][1:]==b['symbols'][1:]==['O','H','H']*n
        shells=[c for c in panel['cells'] if c['water_count']==n and c['role']==ROLES[1]]
        checks['own_shell'+str(n)+'_same_state_and_water_basis']=len(shells)==2 and all(c['charge']==0 and c['multiplicity']==1 and c['expected_electrons']==10*n and c['expected_basis_functions_by_center'][1:]==[31,6,6]*n for c in shells)
    preambles=[text.split('%basis',1)[0] for text in inputs.values()]
    checks['same_global_electronic_recipe']=all(x==preambles[0] for x in preambles) and 'def2-TZVP def2/J RIJCOSX' in preambles[0]
    checks['only_metal_basis_overrides']=all(re.findall(r'(?im)^New(?:ECP|GTO|AuxJGTO)\s+(\w+)',text)==[g.rsplit('_',1)[1]]*3 for g,text in inputs.items())
    checks['exact_consumed_geometry_and_original_authority']=all(panel['geometries'][g]==data['rows'][g]['geometry'] and panel['saved_complexes'][g]['result']==data['rows'][g]['reference_result'] for g in GIDS)
    require(all(checks.values()),'Own-term cancellation source check failed: '+str([k for k,v in checks.items() if not v]))
    return checks

def admit_original(gid,panel,data,adapter):
    saved=panel['saved_complexes'][gid];r=read(saved['result']);old=adapter.b.load_panel(verify(panel['parent_panel']))
    require(r['status']=='complete' and r['returncode']==0 and all(r['checks'].values()),'Original result not admitted')
    same(r['adapter'],panel['base_helper'],'Original collector changed');same(r['panel'],panel['parent_panel'],'Original panel changed')
    c=next(c for c in old['cells'] if c['cell_id']=='ecp_pbe0__'+gid)
    same(c['geometry'],panel['geometries'][gid],'Original geometry changed')
    for key in ('output','input','engrad'):same(r[key],saved[key],'Original '+key+' changed');verify(r[key])
    same(verify(r['input']).read_text(),adapter.b.deck(old,c,4,3000),'Original executed input differs')
    raw=adapter.b.collect_one(old,c,Path(saved['result']['path']).parent,r['returncode'])
    for k,v in raw.items():same(r[k],v,'Original raw replay differs: '+k)
    es,ed,ef=adapter.dft.electronic_parts(verify(r['output']).read_text());require(len(es)==len(ed)==len(ef)==1,'Original energy stages ambiguous')
    terms={'electronic_energy_hartree':es[0],'D4_energy_hartree':ed[0],'final_energy_hartree':ef[0]}
    for k,v in terms.items():require(finite(v) and v==saved[k],'Original energy term differs: '+k)
    require(abs(es[0]+ed[0]-ef[0])<1e-8 and ef[0]==r['energy_hartree'],'Original D4 closure')
    require(r['energy_hartree']*HA==data['rows'][gid]['reference']['energy'] and r['forces_kcal_mol_A']==data['rows'][gid]['reference']['forces'],'Original consumed reference differs')
    return dict(status='complete',**terms,result=saved['result'],source_input=r['input'],output=r['output'],checks=r['checks'],job_id=r['job_id'],allocated_CPU_seconds=r.get('allocated_cpu_seconds'),new_native_calls=0)

def admit_fragment(path,cell,panel,panel_pin,adapter):
    path=Path(path);r=json.loads(path.read_text());require(r['status']=='complete' and r['returncode']==0 and all(r['checks'].values()),'Fragment native admission failed')
    same(r['cell_id'],cell['cell_id'],'Fragment ID differs');same(r['cell'],cell,'Fragment physical definition differs');same(r['panel'],panel_pin,'Fragment panel differs');same(r['adapter'],panel['adapter'],'Fragment adapter differs')
    for k in ('output','stderr','input'):verify(r[k])
    startpath=path.parent/'START.json';start=json.loads(startpath.read_text());same(start['cell'],cell,'Launch cell differs');same(start['panel'],panel_pin,'Launch panel differs');same(start['input'],r['input'],'Launch input differs');same(start['release'],r['release'],'Release lineage differs')
    release=read(r['release']);require(release['status']=='approved_coordination_cp_eight_cells' and release['panel']==panel_pin and release['implementation']==panel['adapter'],'Native release differs')
    require(release['cell_ids']==[c['cell_id'] for c in panel['cells']] and release['maximum_native_calls']==8 and release['ranks']==4,'Native released denominator differs')
    require(path.resolve()==Path(release['execution_root']).resolve()/'tasks'/cell['cell_id']/'execution/RESULT.json','Native attempt path differs')
    verify(release['preparation']);verify(release['independent_review'])
    require(start['ranks']==4 and start['memory_MiB']>0 and start['maxcore_MiB']==math.floor(.75*start['memory_MiB']/4),'Actual rank/memory differs')
    require(start['command']==[panel['orca']['path'],'endpoint.inp','--bind-to none'] and start['maximum_native_calls']==r['native_program_invocations_attempted']==1,'Native command/call count differs')
    same(verify(r['input']).read_text(),adapter.deck(panel,cell,4,start['maxcore_MiB']),'Actual ghost deck differs')
    parsed=adapter.parse(verify(r['output']).read_text(),cell,r['returncode'])
    for k,v in parsed.items():same(r[k],v,'Native fragment replay differs: '+k)
    require(parsed['status']=='complete' and all(parsed['checks'].values()) and finite(parsed['electronic_energy_hartree']) and parsed['D4_energy_hartree']==0.,'Fragment electronic energy unavailable')
    require(abs(parsed['energy_hartree']-parsed['electronic_energy_hartree'])<1e-8,'Ghost energy includes non-electronic contribution')
    return dict(status='complete',result=pin(path),start=pin(startpath),input=r['input'],output=r['output'],release=r['release'],
        electronic_energy_hartree=parsed['electronic_energy_hartree'],D4_energy_hartree=parsed['D4_energy_hartree'],final_energy_hartree=parsed['energy_hartree'],
        checks=parsed['checks'],printed_scf_criteria=parsed['printed_scf_criteria'],native_scf_mode=parsed['native_scf_mode'],native_scf_stop_signals=parsed['native_scf_stop_signals'],
        energy_stage_metadata=parsed['energy_stage_metadata'],final_COSX_correction_hartree=parsed['final_COSX_correction_hartree'],job_id=r['job_id'],allocated_CPU_seconds=r.get('allocated_cpu_seconds'))

def unavailable(path,exc=None):
    p=Path(path);v={'status':'missing' if not p.exists() else 'failed_or_unavailable','path':str(p),'error':None if exc is None else repr(exc)}
    if p.exists():
        v['result']=pin(p)
        try:
            raw=json.loads(p.read_text());v['reported_status']=raw.get('status');v['reported_error']=raw.get('exception',raw.get('error'))
        except Exception:pass
    return v

def report(result):
    lines=['# Coordination counterpoise diagnostic','',f"Status: **{result['status']}**. Admitted fragments {result['fragment_complete']}/8; original complexes {result['original_complete']}/4.",'']
    if result['cycle'] is None:lines+=['The corrected cycle is unavailable. No partial or failed fragment energy was substituted.','']
    else:
        x=result['cycle']['kcal_mol'];lines+=['| Paired W9−W8 quantity | kcal/mol |','|---|---:|']
        for k,label in [('original_electronic','Original electronic'),('original_D4_unchanged','Original D4, unchanged'),('original_total','Original total'),('delta_G_electronic','ΔG electronic'),('CP_shift','CP shift = −ΔG'),('CP_total','Alternative CP total'),('unchanged_model','Unchanged model'),('model_minus_original','Model minus original'),('model_minus_CP','Model minus CP')]:lines.append(f'| {label} | {x[k]:+.9f} |')
        lines+=['','Δ uses Dy9 − La9 − Dy8 + La8. Own-ion terms cancel separately within each metal across sizes; own-shell terms cancel between metals at each size. Newly computed ghost terms do **not** cancel by assumption.','']
    lines+=['| Fragment | Status | Electronic Eh |','|---|---|---:|']
    for cid,v in result['fragments'].items():lines.append(f"| {cid} | {v['status']} | {v['electronic_energy_hartree']:.12f} |" if v['status']=='complete' else f"| {cid} | {v['status']} | unavailable |")
    lines+=['','All original terms, native checks, stopping signals, stage metadata and unavailable-slot details are retained in RESULT.json.',
        '', 'This is metal-versus-whole-shell CP only. Intra-shell BSSE, basis-limit and frozen-core/reference uncertainty remain; W9 also rearranges the other eight waters. No corrected forces, free energies, new fitting targets, parameter selection or sampling release follows.',
        '',f"Recorded new-fragment allocated CPU seconds: {result['cost']['admitted_new_fragment_allocated_CPU_seconds']}; reused original cost is separate. Collector molecular calls: zero."]
    return '\n'.join(lines)+'\n'

def collect(config,out):
    start=time.monotonic();cpu=time.process_time();cfg=json.loads(Path(config).read_text());out=Path(out);out.mkdir(parents=True,exist_ok=False)
    result={'status':'unavailable','implementation':pin(__file__),'config':pin(config),'fragments':{cid:unavailable(cfg.get('result_paths',{}).get(cid,'MISSING_PATH')) for cid in IDS},
        'original_complexes':{g:{'status':'not_admitted'} for g in GIDS},'fragment_complete':0,'original_complete':0,'cycle':None,'cancellation_checks':None,'error':None,'new_native_calls':0,'fit_targets_changed':False,'model_selected':False,'physical_accuracy_qualified':False}
    try:
        require(cfg['implementation']==pin(__file__) and cfg['specification']==pin(HERE/'SPEC.md'),'Collector source/spec differs')
        require(set(cfg['result_paths'])==set(IDS),'Configured fragment denominator differs');verify(cfg['root_inputs']);root=read(cfg['root_inputs']);verify(root['declaration'])
        require(cfg['root_inputs']['sha256']=='77c1ea223d87f99fa25746102f76dc2cb01aae80f48d6293eecb9a4531d247a6','Root input scope differs')
        adapter=module(cfg['adapter'],'coordination_saved_cp');panel=adapter.load(verify(cfg['panel']));same(panel['adapter'],cfg['adapter'],'Panel adapter differs');same(panel['root_inputs'],cfg['root_inputs'],'Panel root lineage differs')
        data=read(root['consumed_input_data']);same(panel['consumed_input_data'],root['consumed_input_data'],'Consumed source differs')
        result.update(panel=cfg['panel'],root_inputs=cfg['root_inputs'],consumed_input_data=root['consumed_input_data'])
        result['cancellation_checks']=cancellation_check(panel,data)
        for gid in GIDS:
            try:result['original_complexes'][gid]=admit_original(gid,panel,data,adapter)
            except Exception as exc:result['original_complexes'][gid]=unavailable(panel['saved_complexes'][gid]['result']['path'],exc)
        for cell in panel['cells']:
            cid=cell['cell_id'];path=cfg['result_paths'][cid]
            try:result['fragments'][cid]=admit_fragment(path,cell,panel,cfg['panel'],adapter)
            except Exception as exc:result['fragments'][cid]=unavailable(path,exc)
        result['fragment_complete']=sum(v['status']=='complete' for v in result['fragments'].values());result['original_complete']=sum(v['status']=='complete' for v in result['original_complexes'].values())
        model={g:data['baseline_model'][g]['energy'] for g in GIDS};same(model,panel['baseline_model_energies'],'Unchanged model lineage differs');result['unchanged_model_energies_kcal_mol']=model
        if result['fragment_complete']==8 and result['original_complete']==4:
            result['cycle']=calculate(result['original_complexes'],result['fragments'],model,result['cancellation_checks']);result['status']='coordination_cp_collection_complete'
        else:result['status']='coordination_cp_collection_incomplete'
    except Exception as exc:result.update(status='coordination_cp_collection_unavailable',error=repr(exc))
    result['cost']={'admitted_new_fragment_allocated_CPU_seconds':sum(v.get('allocated_CPU_seconds') or 0 for v in result['fragments'].values() if v['status']=='complete'),
        'reused_original_allocated_CPU_seconds':sum(v.get('allocated_CPU_seconds') or 0 for v in result['original_complexes'].values() if v['status']=='complete'),
        'costs_for_unadmitted_or_missing_rows_not_included':True,'collector_new_native_calls':0,'collector_wall_seconds':time.monotonic()-start,'collector_CPU_seconds':time.process_time()-cpu,'local_threads':1}
    save(out/'RESULT.json',result);(out/'REPORT.md').write_text(report(result));print(json.dumps({k:result[k] for k in ('status','fragment_complete','original_complete','error')}));return result

def configure(panel_path,task_root,out):
    pp=pin(panel_path);p=read(pp);cfg={'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'root_inputs':pin(HERE.parent/'ROOT_INPUTS.json'),'adapter':p['adapter'],'panel':pp,
        'result_paths':{cid:str(Path(task_root).resolve()/cid/'execution/RESULT.json') for cid in IDS},'new_native_calls':0}
    save(out,cfg);return cfg

if __name__=='__main__':
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    p=sub.add_parser('configure');p.add_argument('--panel',required=True);p.add_argument('--task-root',required=True);p.add_argument('--out',required=True)
    p=sub.add_parser('collect');p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=ap.parse_args()
    if a.action=='configure':configure(a.panel,a.task_root,a.out)
    else:raise SystemExit(0 if collect(a.config,a.out)['status']=='coordination_cp_collection_complete' else 1)
