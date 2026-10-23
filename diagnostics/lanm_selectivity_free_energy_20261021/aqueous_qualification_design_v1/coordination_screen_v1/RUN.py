#!/usr/bin/env python3
"""Two new frozen-candidate W9 cells; all quantum and W8 evidence is saved."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shutil
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
W = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021'
D = ROOT / 'diagnostics/lanm_selectivity_free_energy_20261021'
CELLS = ('aquo9_central_La', 'aquo9_central_Dy')
THREADS = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')

def require(ok, message):
    if not ok: raise ValueError(message)
def pin(path):
    p=Path(path).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']); return Path(p['path'])
def read(p): return json.loads(verify(p).read_text())
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p)); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
def save(path,value):
    p=Path(path); tmp=p.with_name(p.name+'.tmp'); tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n'); os.replace(tmp,p)

def default_config():
    base=D/'explicit_model_v1/paired_water_orientation_validation_v1/model'
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),
        'candidate_helper':pin(base/'RUN.py'),'candidate_config':pin(base/'CONFIG_SOURCE_v2.json'),
        'original_manifest':pin(W/'explicit_model_v1/prepared/MANIFEST.json'),
        'references':{f'aquo{n}_central_{m}':pin(W/f'reference_chemistry_v1/ecp_remaining_v1/tasks/ecp_pbe0__aquo{n}_central_{m}/execution/RESULT.json') for n in (8,9) for m in ('La','Dy')},
        'maximum_native_programs':4,'allocated_cpus':1,'new_FD_energies':0}

def prepare(config_path,out):
    start=time.monotonic();cpu=time.process_time();cfg=json.loads(Path(config_path).read_text());require(cfg==default_config(),'Exact preparation configuration required')
    require(cfg['candidate_helper']['sha256']=='81aaca64d8327a40af48782a61c2ecd26de74c34d7609b3adc33c386c1c1b946','Candidate admission helper differs')
    helper=module(cfg['candidate_helper'],'coord_candidate');ccfg=read(cfg['candidate_config'])
    data,unused,water,transfer,tbound,cartesian=helper.context(ccfg)
    require(not any('aquo9' in g for v in water['vectors'].values() for g in v.get('model',{})),'W9 appeared in candidate fit')
    require(not any('aquo9' in g for g in transfer['rows']),'W9 appeared in candidate transfer')
    original=read(cfg['original_manifest']);old={r['id']:r for r in original['rows']};body=module(ccfg['body_helper'],'coord_body')
    reference={};geometries={};states={};raw_records={}
    for gid,p in cfg['references'].items():
        raw=read(p);require(raw['status']=='complete' and raw['returncode']==0 and all(raw['checks'].values()),'Saved reference failed')
        require(raw['adapter']['sha256']=='eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129','Reference source differs')
        native=module(raw['adapter'],'coord_reference_'+gid);panel=native.load_panel(verify(raw['panel']));cell=next(c for c in panel['cells'] if c['cell_id']=='ecp_pbe0__'+gid)
        replay=native.collect_one(panel,cell,Path(p['path']).parent,0);require(all(raw[k]==v for k,v in replay.items()),'Reference raw replay differs')
        n=9 if gid in CELLS else 8;metal=cell['metal']
        require(cell['recipe']=='ecp_pbe0' and cell['charge']==3 and cell['multiplicity']==1 and cell['physical_multiplicity']==(1 if metal=='La' else 6) and cell['ecp_core']==(46 if metal=='La' else 55) and cell['explicit_electrons']==8+10*n,'Reference state mismatch')
        require(cell['geometry']['atoms']==1+3*n and cell['geometry']['water_count']==n,'Reference inventory mismatch')
        reference[gid]={'energy':raw['energy_hartree']*native.HA,'forces':raw['forces_kcal_mol_A'],'raw_result':p,'physical_state_qualified':False,'method_accuracy_qualified':False}
        geometries[gid]=cell['geometry'];states[gid]={k:cell[k] for k in ('charge','multiplicity','physical_multiplicity','ecp_core','explicit_electrons')};raw_records[gid]=raw
    for metal in ('La','Dy'):
        a,b=(f'aquo{n}_central_{metal}' for n in (8,9));require(raw_records[a]['panel']==raw_records[b]['panel'] and raw_records[a]['adapter']==raw_records[b]['adapter'] and raw_records[a]['input']['sha256']==raw_records[b]['input']['sha256'],'W8/W9 reference recipe differs')
        anchor=copy.deepcopy(data['rows'][a]['geometry']);verify(anchor['xyz']);verify(geometries[a]['xyz']);require(anchor['xyz']['sha256']==geometries[a]['xyz']['sha256'],'W8 geometry alias bytes differ');anchor['xyz']=geometries[a]['xyz'];require(geometries[a]==anchor,'W8 anchor geometry differs')
    rows={}
    for gid in CELLS:
        g=geometries[gid];metal=g['metal'];prior=copy.deepcopy(old[gid]);parent='aquo8_central_'+metal
        for key in ('xyz','key','parameters'):verify(prior[key])
        source=copy.deepcopy(prior['source_geometry']);verify(source['xyz']);require(source['xyz']['sha256']==g['xyz']['sha256'],'Historical geometry alias differs');source['xyz']=g['xyz'];require(source==g,'Original W9 geometry differs')
        require(set(prior['types'])=={36,37,500 if metal=='La' else 501} and len(prior['types'])==28,'Nonwater W9 atom types')
        lines=verify(prior['xyz']).read_text().splitlines();require(int(lines[0].split()[0])==28 and len(lines)==29,'Typed W9 atom count')
        for i,line in enumerate(lines[1:]):
            t=line.split();require(int(t[0])==i+1 and t[1]==g['symbols'][i] and [float(v) for v in t[2:5]]==g['positions_A'][i] and int(t[5])==prior['types'][i] and [int(v)-1 for v in t[6:]]==prior['bonds_zero_based'][i],'Typed W9 coordinate/topology differs')
        directory=out/'inputs'/gid;directory.mkdir(parents=True)
        source_files={'xyz':prior['xyz'],**{k:tbound['inputs'][parent]['files'][k] for k in ('key','parameters')}}
        files={}
        for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
            shutil.copyfile(verify(source_files[k]),directory/name);files[k]=pin(directory/name);require(files[k]['sha256']==source_files[k]['sha256'],'Input copy differs')
        prior.update(**files,source_geometry=g,derivative_test=False)
        row={'geometry':g,'input':prior,'files':files,'source_files':source_files,'reference':reference[gid],'state':states[gid]};row['bodies']=body.bodies(row);require(len(row['bodies'])==10 and sum(b['category']=='water' for b in row['bodies'])==9,'W9 body partition differs');rows[gid]=row
    require(geometries[CELLS[0]]['positions_A']==geometries[CELLS[1]]['positions_A'],'Paired W9 coordinates differ')
    result={'status':'prepared_no_calls','implementation':pin(__file__),'config':pin(config_path),'rows':rows,'candidate_result':ccfg['water_result'],'candidate_parameters':water['selected']['parameters'],
        'anchor_model':{m:copy.deepcopy(transfer['rows']['aquo8_central_'+m]['model']) for m in ('La','Dy')},
        'anchor_reference':{m:reference['aquo8_central_'+m] for m in ('La','Dy')},'candidate_config':cfg['candidate_config'],'executables':ccfg['executables'],'native_helper':ccfg['native_helper'],'numerical_helper':ccfg['numerical_helper'],'body_helper':ccfg['body_helper'],
        'prior_candidate_W9_absent':True,'prior_fit_vectors_inspected':len(water['vectors']),'prior_transfer_cells_inspected':len(transfer['rows']),'original_baseline_W9_preserved':True,'new_native_calls':0,'physical_accuracy_qualified':False,'new_FD_energies':0,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    save(out/'BOUND.json',result);return result

def score(bound,values,body):
    missing=[g for g in CELLS if values.get(g) is None]
    unknown={'status':'unavailable','passed':None,'missing':missing}
    result={'energy':copy.deepcopy(unknown),'rigid':{m:copy.deepcopy(unknown) for m in ('La','Dy','paired')},'declared_comparisons':4,'available_comparisons':0,'passed':False,'physical_accuracy_qualified':False}
    def admitted(v):
        return isinstance(v['energy'],(int,float)) and not isinstance(v['energy'],bool) and math.isfinite(v['energy']) and len(v['forces'])==28 and all(len(a)==3 and all(math.isfinite(x) for x in a) for a in v['forces'])
    for g in CELLS:
        if values.get(g) is not None:require(admitted(values[g]),'Invalid candidate E/F')
    refs={g:bound['rows'][g]['reference'] for g in CELLS}
    if not missing:
        def contrast(v,a):return (v[CELLS[1]]['energy']-v[CELLS[0]]['energy'])-(a['Dy']['energy']-a['La']['energy'])
        mv,rv=contrast(values,bound['anchor_model']),contrast(refs,bound['anchor_reference']);error=mv-rv
        result['energy']={'status':'complete','model':mv,'reference':rv,'error_kcal_mol':error,'passed':abs(error)<=1};result['available_comparisons']+=1
    for metal in ('La','Dy','paired'):
        needed=CELLS if metal=='paired' else ('aquo9_central_'+metal,)
        if any(values.get(g) is None for g in needed):continue
        if metal=='paired':a=body.subtract(values[CELLS[1]]['forces'],values[CELLS[0]]['forces']);b=body.subtract(refs[CELLS[1]]['forces'],refs[CELLS[0]]['forces'])
        else:a=values[needed[0]]['forces'];b=refs[needed[0]]['forces']
        weighted,raw=body.rigid_residual(bound['rows'][CELLS[0]],a,b);rms=2*math.sqrt(sum(x*x for x in weighted));maximum=max(abs(x) for x in raw)
        result['rigid'][metal]={'status':'complete','RMS_kcal_mol_A':rms,'maximum_component_kcal_mol_A':maximum,'passed':rms<=2 and maximum<=5};result['available_comparisons']+=1
    result['passed']=result['available_comparisons']==4 and result['energy']['passed'] and all(v['passed'] for v in result['rigid'].values());return result

def run(bound_path,release_path,out):
    start=time.monotonic();cpu=time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN);bp=pin(bound_path);bound=read(bp);cfg=read(bound['config']);release=json.loads(Path(release_path).read_text())
    require(cfg['implementation']==bound['implementation']==pin(__file__) and bound['status']=='prepared_no_calls','Source/binding changed');verify(cfg['specification']);verify(bound['candidate_result'])
    require(release['status']=='approved_frozen_W9_development_screen' and release['implementation']==pin(__file__) and release['bound']==bp and release['maximum_native_programs']==4 and release['new_FD_energies']==0,'Exact root release required')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1 and int(os.environ.get('SLURM_NTASKS','1'))==1 and all(os.environ.get(k)=='1' for k in THREADS),'One allocated native thread required')
    for p in bound['executables'].values():verify(p)
    native=module(bound['native_helper'],'W9_native');helper=module(bound['numerical_helper'],'W9_numerical');body=module(bound['body_helper'],'W9_body')
    result={'status':'running','implementation':pin(__file__),'bound':bp,'release':pin(release_path),'job_id':os.environ['SLURM_JOB_ID'],'rows':{g:{'status':'not_attempted','model':None,'receipts':[]} for g in CELLS},'native_program_invocations':0,'new_FD_energies':0,'error':None,'physical_accuracy_qualified':False}
    def persist():save(out/'RESULT.json',result)
    try:
        persist()
        for gid in CELLS:
            row=bound['rows'][gid];target=result['rows'][gid];directory=out/gid;directory.mkdir()
            for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['files'][k]),directory/name);require(pin(directory/name)['sha256']==row['files'][k]['sha256'],'Runtime copy mismatch')
            target.update(status='started',input=row['files']);persist()
            for program,args in (('analyze',['system.xyz','E']),('testgrad',['system.xyz','Y','N'])):
                require(result['native_program_invocations']<4,'Native program ceiling');result['native_program_invocations']+=1;persist();output=directory/(program+'.out')
                native.invoke([bound['executables'][program]['path']]+args,directory,output,target['receipts']);helper.program_admission(output,program,native);persist()
            parsed=native.parse(directory/'analyze.out',directory/'testgrad.out',row['input']);value={'energy':parsed['total_energy_kcal_mol'],'forces':parsed['force_kcal_mol_A'],'native':parsed};helper.replay_native(value,row,native);target.update(status='complete',model=value);persist()
        result['status']='complete'
    except Exception as exc:
        result.update(status='native_or_admission_failure',error=repr(exc))
        for v in result['rows'].values():
            if v['status']=='started':v.update(status='failed',error=repr(exc))
    finally:
        result['comparison']=score(bound,{g:v['model'] for g,v in result['rows'].items()},body);now=resource.getrusage(resource.RUSAGE_CHILDREN);result.update(elapsed_seconds=time.monotonic()-start,parent_CPU_seconds=time.process_time()-cpu,aggregate_native_child_user_CPU_seconds=now.ru_utime-child.ru_utime,aggregate_native_child_system_CPU_seconds=now.ru_stime-child.ru_stime);persist()
    return result

def main():
    parser=argparse.ArgumentParser();sub=parser.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--config',required=True)
    a=sub.add_parser('run');a.add_argument('--bound',required=True);a.add_argument('--release',required=True)
    for a in sub.choices.values():a.add_argument('--out',required=True)
    args=parser.parse_args();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:r=prepare(args.config,out) if args.action=='prepare' else run(args.bound,args.release,out)
    except Exception as exc:
        r={'status':'adapter_failed','error':repr(exc),'implementation':pin(__file__),'new_native_calls':0 if args.action=='prepare' else None};save(out/'FAILURE.json',r)
    print(json.dumps({'status':r['status'],'out':str(out)}));return 0 if r['status'] in ('prepared_no_calls','complete') else 1
if __name__=='__main__':raise SystemExit(main())
