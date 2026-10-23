#!/usr/bin/env python3
"""Sixteen consumed dry-mono-acetate SPs, staged4/12; shared strict native parser."""
import argparse
import copy
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('dry_amide_helper',HERE.parent/'hydrated_amide_force_qualification_v1/REFERENCE.py')
h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
w=h.w;b=h.b
PROTOCOL='dry_acetate_force_qualification_v1'
METALS=('La','Dy'); DIRECTIONS=('radial','rotation')
SHIFTS={'minus_010':-.01,'minus_005':-.005,'plus_005':.005,'plus_010':.01}
STEPS=(.01,.005); RESOLUTION=.1; LENGTH_A=1.
WS=b.ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
TRAINING=WS/'explicit_model_v1/matched_dry_v1/training_v1/RESULT.json'
TRAINING_SHA='99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb'
CENTRAL_SHA={'La':'4301f2f00874ae52feb59d3db48edb2d060c36e34f468b92f6c1e473ae2d64b5','Dy':'b4f813d35e4219d606a09a3fc9c96dbec2934c83660fa0d1e31588246799e587'}
HELPER_SHA='7de9c17aa9198ba07c66c3f934390fc2a67d94a77860eca8ba4c0378b25ad677'
WATER_SPEC_SHA='f313f90a541b28d77a356c140944ab7a1b0d80e5971c0efc26c23b7f79461d96'
AMIDE_SPEC_SHA='f5de75cfa8b006a07e93267f6ebcdeef1d02e37a31d544872656698bc92112d4'

def source_gates():
    h.source_gates()
    b.check(b.pin(h.__file__)['sha256']==HELPER_SHA,'Amide helper changed')
    b.check(b.pin(Path(w.__file__).with_name('SPEC.md'))['sha256']==WATER_SPEC_SHA,'Resolution authority changed')
    b.check(b.pin(Path(h.__file__).with_name('SPEC.md'))['sha256']==AMIDE_SPEC_SHA,'Amide resolution authority changed')
    b.check(h.RESOLUTION==w.RESOLUTION==RESOLUTION and h.SHIFTS==SHIFTS,'Steps/resolution changed')

def cid(m,d,s):return f'pbe0__{m}__{d}__{s}__complex'
ALL=tuple(cid(m,d,s) for m in METALS for d in DIRECTIONS for s in SHIFTS)
SCOUT=tuple(cid('La','radial',s) for s in SHIFTS)
REMAINING=tuple(c for c in ALL if c not in SCOUT)

def frame(g):
    b.check(g['atoms']==8 and g['metal_index']==0 and g['ligand_indices']==list(range(1,8)) and g['donor_indices']==[3],'Exact dry mono frame differs')
    b.check(g['symbols'][3:5]==['O','O'] and g['water_count']==0,'Exact donor O indices differ')
    x=g['positions_A'];u=g['radial_direction'];center=[sum(x[i][k] for i in range(1,8))/7 for k in range(3)]
    b.check(abs(b.norm(u)-1)<1e-12 and b.norm(b.sub(u,b.unit(b.sub(x[3],x[0]))))<1e-12,'Radial direction differs')
    cross=b.cross(u,b.sub(x[4],x[3]));b.check(b.norm(cross)>1e-12,'Degenerate geometry-defined rotation axis')
    return {'radial_unit':u,'centroid_A':center,'rotation_unit':b.unit(cross),'rotation_length_A':LENGTH_A,'ligand_indices':list(range(1,8))}

def transform(g,direction,q):
    f=frame(g);x=copy.deepcopy(g['positions_A'])
    for i in f['ligand_indices']:
        if direction=='radial':x[i]=b.add(x[i],b.mul(f['radial_unit'],q))
        else:
            b.check(direction=='rotation','Unknown direction');v=b.sub(x[i],f['centroid_A']);n=f['rotation_unit'];theta=q/LENGTH_A
            rotated=b.add(b.mul(v,math.cos(theta)),b.add(b.mul(b.cross(n,v),math.sin(theta)),b.mul(n,b.dot(n,v)*(1-math.cos(theta)))))
            x[i]=b.add(f['centroid_A'],rotated)
    return x

def projection(forces,g):
    f=frame(g);ids=f['ligand_indices'];net=[sum(forces[i][k] for i in ids) for k in range(3)]
    torques=[b.cross(b.sub(g['positions_A'][i],f['centroid_A']),forces[i]) for i in ids]
    torque=[sum(t[k] for t in torques) for k in range(3)]
    return {'radial':-b.dot(net,f['radial_unit']),'rotation':-b.dot(torque,f['rotation_unit'])/LENGTH_A,
        'ligand_net_force_kcal_mol_A':net,'ligand_torque_kcal_mol':torque,'total_force_kcal_mol_A':[sum(v[k] for v in forces) for k in range(3)],'frame':f}

def deck(panel,cell,ranks,maxcore):return h.deck(panel,cell,ranks,maxcore)

def xc_metadata(text):
    sections=re.findall(r'^DFT GRID GENERATION\s*\n(.*?)^COSX GRID GENERATION',text,re.M|re.S)
    post=re.findall(r'^\s*Gradients option\s+PostSCFGGA\s+\.{2,}\s+(\S+)',text,re.M)
    b.check(len(sections)==1 and len(post)==1,'Native DFT/XC grid or PostSCFGGA metadata missing/ambiguous')
    lines=[s.strip() for s in sections[0].splitlines() if '...' in s or 'Angular grids for' in s]
    count_words=('Total number','Average number')
    return {'settings':[s for s in lines if not s.startswith(count_words)],'counts':[s for s in lines if s.startswith(count_words)],'PostSCFGGA':post[0]}

def parse(text,cell,returncode):
    result=h.parse(text,cell,returncode);result['dft_xc_grid_metadata']=xc_metadata(text);return result

def prepare(out):
    tick=time.monotonic();cpu=time.process_time();source_gates()
    b.check(b.pin(TRAINING)['sha256']==TRAINING_SHA,'Consumed training changed');data=json.loads(TRAINING.read_text())
    b.check(data['status']=='training_complete' and data['complete']==28,'Original training incomplete')
    out.mkdir(parents=True,exist_ok=False)
    for name in ('assets','xyz','inputs'):(out/name).mkdir()
    central={};assets={};cells=[];coverage=[];orca=None
    for name,row in data['rows'].items():
        if not name.startswith('matched_dry_'):continue
        saved=json.loads(b.verify(row['reference_result']).read_text())
        coverage.append({'geometry_id':name,'geometry':row['geometry']['xyz'],'reference_result':row['reference_result'],'native_analytic_admitted':saved['status']=='complete' and all(saved['checks'].values()),'matched_small_step_reference_radial_qualification_found':False,'matched_reference_rotation_qualification_found':False})
    b.check(len(coverage)==6 and all(x['native_analytic_admitted'] for x in coverage),'Six consumed centers unavailable')
    for metal in METALS:
        row=data['rows']['matched_dry_acetate_mono_central_'+metal];g=row['geometry'];saved_pin=row['reference_result']
        b.check(saved_pin['sha256']==CENTRAL_SHA[metal],'Wrong consumed center');saved=json.loads(b.verify(saved_pin).read_text())
        parent=b.load_panel(b.verify(saved['panel']));original=next(c for c in parent['cells'] if c['cell_id']==saved['cell_id'])
        b.check(original['geometry']==g and original['recipe']=='ecp_pbe0' and saved['geometry']==g['xyz'],'Center geometry/recipe differs')
        b.check(original['charge']==2 and original['multiplicity']==1 and original['explicit_electrons']==40 and original['ecp_core']==(46 if metal=='La' else 55),'Center state differs')
        start_path=Path(saved_pin['path']).parent/'START.json';start=json.loads(start_path.read_text())
        b.check(start['cell']==original and start['panel']==saved['panel'] and start['input']==saved['input'],'Center runtime join differs')
        b.check(b.verify(saved['input']).read_text()==b.deck(parent,original,start['ranks'],start['maxcore_MiB']),'Executed center deck not reproduced')
        replay=b.collect_one(parent,original,Path(saved_pin['path']).parent,saved['returncode'])
        b.check(replay['status']=='complete' and replay['forces_kcal_mol_A']==row['reference']['forces']==saved['forces_kcal_mol_A'] and replay['energy_hartree']==row['reference']['energy_hartree']==saved['energy_hartree'],'Central raw E/F replay differs')
        meta=h.metadata(g);raw_output=b.verify(saved['output']).read_text();strict=w.parse(raw_output,meta,saved['returncode']);strict['dft_xc_grid_metadata']=xc_metadata(raw_output)
        b.check(strict['status']=='complete' and strict['energy_hartree']==saved['energy_hartree'] and meta['expected_basis_functions']==189,'Strict central native replay failed')
        if orca is None:orca=parent['orca']
        b.check(orca==parent['orca'],'Mixed native binary');assets[metal]={}
        for name in ('basis','aux'):
            dst=out/'assets'/f'{metal}_{name}.inc';shutil.copyfile(b.verify(parent['assets'][metal][name]),dst);assets[metal][name]=b.pin(dst)
        central[metal]={'original_cell':original,'saved_result':saved_pin,'saved_start':b.pin(start_path),'executed_input':saved['input'],'output':saved['output'],'engrad':saved['engrad'],'original_panel':saved['panel'],'original_adapter':saved['adapter'],'energy_hartree':saved['energy_hartree'],'analytic':projection(saved['forces_kcal_mol_A'],g),'strict_native_replay':strict,'original_resource_values':{'ranks':start['ranks'],'maxcore_MiB':start['maxcore_MiB']}}
        for direction in DIRECTIONS:
            for label,q in SHIFTS.items():
                cell_id=cid(metal,direction,label);xyz=out/'xyz'/(cell_id+'.xyz');b.xyz_write(xyz,g['symbols'],transform(g,direction,q),f'Consumed dry mono-acetate {direction}; q={q} A; metal fixed')
                symbols,x=b.xyz_read(xyz);ng=copy.deepcopy(g);ng.update(positions_A=x,xyz=b.pin(xyz));cell=h.metadata(ng)
                cell.update(cell_id=cell_id,direction=direction,shift_label=label,q_A=q);cells.append(cell)
    p={'protocol':PROTOCOL,'created_utc':b.now(),'adapter':b.pin(__file__),'specification':b.pin(HERE/'SPEC.md'),'training_result':b.pin(TRAINING),'orca':orca,'assets':assets,'central':central,'cells':cells,'helpers':{'amide':b.pin(h.__file__),'water48':b.pin(w.__file__),'original':b.pin(b.__file__),'donor':b.pin(w.helper.__file__),'dft':b.pin(w.dft.__file__),'canonical':b.pin(w.cp.__file__)},'resolution_authorities':[b.pin(Path(w.__file__).with_name('SPEC.md')),b.pin(Path(h.__file__).with_name('SPEC.md'))],'maximum_native_calls':16,'ranks':4,'steps_A':list(STEPS),'rotation_length_A':LENGTH_A,'derivative_resolution_kcal_mol_A':RESOLUTION,'scout_cell_ids':list(SCOUT),'remaining_cell_ids':list(REMAINING),'coverage':coverage,'physical_state_qualified':False,'method_accuracy_qualified':False,'validation_outcomes_opened':False}
    for cell in cells:
        path=out/'inputs'/(cell['cell_id']+'.inp');path.write_text(deck(p,cell,4,3000));cell['template_input']=b.pin(path)
    pp=b.write(out/'PANEL.json',p);val=validate(pp['path']);val.update(central_raw_replays=2,cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu},new_native_calls=0)
    pf=b.write(out/'PREFLIGHT.json',val)
    print(json.dumps({'panel':pp,'preflight':pf,'status':val['status'],'cost':val['cost'],'new_native_calls':0},indent=2))

def load(path):
    source_gates();p=json.loads(Path(path).read_text());b.check(p['protocol']==PROTOCOL and p['adapter']==b.pin(__file__),'Protocol/source differs');w.recursive_verify(p)
    b.check(len(p['cells'])==16 and {c['cell_id'] for c in p['cells']}==set(ALL) and p['maximum_native_calls']==16,'Sixteen-slot scope differs')
    b.check(p['scout_cell_ids']==list(SCOUT) and p['remaining_cell_ids']==list(REMAINING) and p['steps_A']==list(STEPS) and p['ranks']==4 and p['derivative_resolution_kcal_mol_A']==RESOLUTION and p['rotation_length_A']==LENGTH_A,'Staging/settings differ')
    return p

def validate(path):
    p=load(path);records=[]
    for c in p['cells']:
        origin=p['central'][c['metal']];g=origin['original_cell']['geometry'];symbols,x=b.xyz_read(b.verify(c['source_geometry']));expected=transform(g,c['direction'],c['q_A'])
        b.check(symbols==g['symbols']==c['symbols'] and x==c['positions_A'] and x[0]==g['positions_A'][0],'Labels/coordinates/metal differ')
        err=max(abs(x[i][j]-expected[i][j]) for i in range(8) for j in range(3));b.check(err<=1e-12,'Prescribed rigid transform differs')
        ids=g['ligand_indices'];internal=max(abs(b.norm(b.sub(x[i],x[j]))-b.norm(b.sub(g['positions_A'][i],g['positions_A'][j]))) for i in ids for j in ids);b.check(internal<=1e-9,'Rigid internal distances differ')
        r=origin['original_resource_values'];b.check(deck(p,c,r['ranks'],r['maxcore_MiB'])==b.verify(origin['executed_input']).read_text().replace(' EnGrad',''),'Scientific deck delta differs')
        raw=deck(p,c,4,3000);b.check(raw==b.verify(c['template_input']).read_text() and 'EnGrad' not in raw and '* xyzfile 2 1 geometry.xyz' in raw,'Template/state differs')
        b.check(c['expected_electrons']==40 and c['expected_basis_functions']==189 and not any(c['ghost_mask']),'Native state/dimensions differ')
        records.append({'cell_id':c['cell_id'],'maximum_transform_error_A':err,'maximum_internal_distance_change_A':internal,'scientific_deck_only_EnGrad_removed':True})
    for d in DIRECTIONS:
        for s in SHIFTS:
            aa=[next(c for c in p['cells'] if c['cell_id']==cid(m,d,s)) for m in METALS];b.check(aa[0]['positions_A']==aa[1]['positions_A'],'Paired geometry differs')
    return {'status':'passed_static_preflight','panel':b.pin(path),'records':records,'cells':16,'new_native_calls':0,'threads':1}

def release_gate(p,path,cell_id):
    r=json.loads(b.verify(path).read_text()) if isinstance(path,dict) else json.loads(Path(path).read_text())
    pp=b.pin(p['_path']);stage=r['stage'];expected=SCOUT if stage=='scout' else REMAINING
    b.check(stage in ('scout','remaining') and r['status']=='released_dry_acetate_force_'+stage,'Unreleased stage')
    b.check(r['implementation']==b.pin(__file__) and r['panel']==pp and r['ranks']==4 and r['maximum_total_native_calls']==16 and r['maximum_stage_calls']==len(expected),'Release binding/scope differs')
    b.check(r['cell_ids']==list(expected) and cell_id in expected and set(r['task_output_directories'])==set(expected),'Exact stage inventory required')
    review=json.loads(b.verify(r['independent_source_review']).read_text());b.check(review['status']=='passed_dry_acetate_force_source_review' and review['implementation']==b.pin(__file__) and review['panel']==pp,'Source review differs')
    if stage=='remaining':
        scout=json.loads(b.verify(r['scout_result']).read_text());sr=json.loads(b.verify(r['scout_review']).read_text())
        b.check(scout['status']=='dry_acetate_force_scout_complete' and scout['stage']=='scout' and scout['adapter']==b.pin(__file__) and scout['panel']==pp and scout['admitted_cells']==4 and scout['stage_numerical_checks_passed'] is True,'Actual scout qualification required')
        b.check(sr['status']=='passed_dry_acetate_force_scout' and sr['result']==r['scout_result'] and sr['implementation']==b.pin(__file__) and sr['panel']==pp,'Actual root scout review differs')
    return r

def run(path,cell_id,release_path,out):
    p=load(path);p['_path']=str(Path(path).resolve());r=release_gate(p,release_path,cell_id);c=next(c for c in p['cells'] if c['cell_id']==cell_id)
    b.check(out.resolve()==Path(r['task_output_directories'][cell_id]).resolve(),'Released output path differs')
    b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native work requires Slurm');ranks=int(os.environ.get('SLURM_NTASKS','0'));b.check(ranks==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread MPI ranks required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));b.check(mem>0,'Actual memory grant required');maxcore=math.floor(.75*mem/ranks)
    out.mkdir(parents=True,exist_ok=False);shutil.copyfile(b.verify(c['source_geometry']),out/'geometry.xyz');(out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore));env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[key]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','');cmd=[p['orca']['path'],'endpoint.inp','--bind-to none']
    b.write(out/'START.json',{'cell':c,'panel':b.pin(path),'adapter':b.pin(__file__),'release':b.pin(release_path),'input':b.pin(out/'endpoint.inp'),'geometry':b.pin(out/'geometry.xyz'),'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'ranks':ranks,'memory_MiB':mem,'maxcore_MiB':maxcore,'command':cmd,'utc':b.now()});tick=time.monotonic()
    with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:proc=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False)
    elapsed=time.monotonic()-tick
    try:result=parse((out/'endpoint.out').read_text(errors='replace'),c,proc.returncode)
    except Exception as exc:result={'status':'failed_collection','exception':repr(exc)}
    result.update(cell_id=cell_id,cell=c,panel=b.pin(path),adapter=b.pin(__file__),release=b.pin(release_path),returncode=proc.returncode,elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],new_native_calls=1,output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'),input=b.pin(out/'endpoint.inp'),geometry=b.pin(out/'geometry.xyz'),finished_utc=b.now())
    b.write(out/'RESULT.json',result);return 0 if result['status']=='complete' else 1

def score(p,rows):
    directions={};checks={};components=('electronic','D4','PBE0_D4');fields=('electronic_energy_hartree','D4_energy_hartree','energy_hartree')
    for direction in DIRECTIONS:
        dd={}
        for kind in ('La','Dy','paired'):
            needed=METALS if kind=='paired' else (kind,)
            available=all(cid(m,direction,s) in rows for m in needed for s in SHIFTS)
            rec=None
            if available:
                derivatives={};analytic=p['central'][kind]['analytic'][direction] if kind!='paired' else p['central']['Dy']['analytic'][direction]-p['central']['La']['analytic'][direction]
                for step in STEPS:
                    suffix='010' if step==.01 else '005';v={}
                    for component,field in zip(components,fields):
                        per={m:(rows[cid(m,direction,'plus_'+suffix)][field]-rows[cid(m,direction,'minus_'+suffix)][field])*b.HA/(2*step) for m in needed}
                        v[component]=per[kind] if kind!='paired' else per['Dy']-per['La']
                    derivatives[str(step)]=v
                delta={k:derivatives['0.005'][k]-derivatives['0.01'][k] for k in components};errors={str(step):derivatives[str(step)]['PBE0_D4']-analytic for step in STEPS}
                b.check(all(math.isfinite(v) for d in derivatives.values() for v in d.values()) and math.isfinite(analytic),'Nonfinite derivative')
                rec={'analytic_kcal_mol_A':analytic,'derivatives_kcal_mol_A':derivatives,'smaller_minus_larger_step':delta,'component_two_step_resolved':{k:abs(v)<=RESOLUTION for k,v in delta.items()},'total_FD_minus_analytic':errors,'prespecified_smaller_step_summary':{k:derivatives['0.005'][k] if abs(delta[k])<=RESOLUTION else None for k in components}}
            dd[kind]=rec
            for name in ('two_step','analytic_0.01','analytic_0.005'):
                error=None if rec is None else (rec['smaller_minus_larger_step']['PBE0_D4'] if name=='two_step' else rec['total_FD_minus_analytic'][name.split('_',1)[1]])
                checks[direction+'/'+kind+'/'+name]={'error_kcal_mol_A':error,'passed':None if error is None else abs(error)<=RESOLUTION,'resolution_kcal_mol_A':RESOLUTION}
        directions[direction]=dd
    return {'directions':directions,'primary_checks':checks,'declared_checks':18,'available_checks':sum(v['passed'] is not None for v in checks.values()),'passed_checks':sum(v['passed'] is True for v in checks.values()),'all_numerical_checks_passed':all(v['passed'] is True for v in checks.values())}

def admit(p,pp,filename):
    result=json.loads(Path(filename).read_text());cell_id=result['cell_id'];c=next(c for c in p['cells'] if c['cell_id']==cell_id)
    b.check(result['status']=='complete' and all(result['checks'].values()) and result['cell']==c and result['panel']==pp and result['adapter']==b.pin(__file__) and result['new_native_calls']==1,'Native result unavailable or wrong authority')
    release=release_gate(p,result['release'],cell_id);directory=Path(filename).resolve().parent;b.check(directory==Path(release['task_output_directories'][cell_id]).resolve(),'Native directory differs')
    replay=parse(b.verify(result['output']).read_text(),c,result['returncode']);b.check(replay['status']=='complete','Raw native admission failed')
    for k,v in json.loads(json.dumps(replay)).items():b.check(result[k]==v,'Stored/native parse differs: '+k)
    start=json.loads((directory/'START.json').read_text());b.check(start['cell']==c and start['panel']==pp and start['adapter']==b.pin(__file__) and start['release']==result['release'] and start['input']==result['input'] and start['geometry']==result['geometry'] and start['ranks']==4,'Native START joins differ')
    b.check(start['maxcore_MiB']==math.floor(.75*start['memory_MiB']/4) and b.verify(result['input']).read_text()==deck(p,c,4,start['maxcore_MiB']),'Executed deck/resource policy differs')
    b.check(b.verify(result['geometry']).read_bytes()==b.verify(c['source_geometry']).read_bytes(),'Executed XYZ differs')
    return result

def collect(path,inventory_path,stage,out):
    tick=time.monotonic();cpu=time.process_time();p=load(path);p['_path']=str(Path(path).resolve());pp=b.pin(path);inventory=json.loads(Path(inventory_path).read_text());expected=SCOUT if stage=='scout' else ALL
    b.check(stage in ('scout','full') and set(inventory)==set(expected),'Exact collection inventory required')
    rows={};slots={c:{'status':'not_requested' if c not in expected else 'unavailable'} for c in ALL};costs=[];failures={};reviews=[]
    for cell_id in expected:
        filename=inventory[cell_id]
        if filename is None or not Path(filename).is_file():continue
        try:
            saved=json.loads(Path(filename).read_text())
            if saved.get('new_native_calls')==1:
                elapsed=saved['elapsed_seconds'];allocated=saved['allocated_cpu_seconds']
                b.check(all(math.isfinite(v) and v>=0 for v in (elapsed,allocated)),'Invalid native cost')
                costs.append({'cell_id':cell_id,'elapsed_seconds':elapsed,'allocated_cpu_seconds':allocated,'job_id':saved['job_id'],'reported_status':saved['status']})
            result=admit(p,pp,filename);b.check(result['cell_id']==cell_id,'Result inventory ID differs');rows[cell_id]=result;slots[cell_id]={'status':'complete','result':b.pin(filename)}
            old=p['central'][result['cell']['metal']]['strict_native_replay'];xc=result['dft_xc_grid_metadata']
            reviews.append({'cell_id':cell_id,'printed_scf_criteria':result['printed_scf_criteria'],'SCF_mode':result['native_scf_mode'],'SCF_stop_signals':result['native_scf_stop_signals'],'energy_stages':result['native_energy_stages'],'center_stage_comparison':w.compare_stages(old['native_energy_stages'],result['native_energy_stages']),'dft_xc_grid_metadata':xc,'center_xc_comparison':{'settings_equal':old['dft_xc_grid_metadata']['settings']==xc['settings'],'PostSCFGGA_equal':old['dft_xc_grid_metadata']['PostSCFGGA']==xc['PostSCFGGA'],'counts_equal':old['dft_xc_grid_metadata']['counts']==xc['counts']}})
        except Exception as exc:slots[cell_id]={'status':'unavailable','result':b.pin(filename),'error':repr(exc)};failures[cell_id]=repr(exc)
    scored=score(p,rows);relevant=[v for k,v in scored['primary_checks'].items() if stage=='full' or k.startswith('radial/La/')]
    complete=len(rows)==len(expected);passed=complete and all(v['passed'] is True for v in relevant)
    status=('dry_acetate_force_scout_complete' if stage=='scout' else 'dry_acetate_force_complete') if complete else 'dry_acetate_force_incomplete'
    result={'status':status,'adapter':b.pin(__file__),'panel':pp,'stage':stage,'inventory':b.pin(inventory_path),'slots':slots,'admitted_cells':len(rows),'stage_expected_cells':len(expected),'total_declared_cells':16,'stage_numerical_checks_passed':passed,'scoring':scored,'native_reviews':reviews,'failures':failures,'native_costs':costs,'native_allocated_cpu_seconds_sum':sum(r['allocated_cpu_seconds'] for r in costs),'new_native_calls':0,'physical_state_qualified':False,'method_accuracy_qualified':False,'held_outcomes_opened':False,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    out.mkdir(parents=True,exist_ok=False);rp=b.write(out/'RESULT.json',result)
    lines=['# Dry mono-acetate reference derivative diagnostic','',f'Stage {stage}: {len(rows)}/{len(expected)} admitted; numerical checks passed: {passed}. Total inventory16. No center rerun.','', '| Direction | Comparison | Analytic | FD0.010 | FD0.005 |','|---|---|---:|---:|---:|']
    for d,groups in scored['directions'].items():
        for m,v in groups.items():
            if v is not None:lines.append(f"| {d} | {m} | {v['analytic_kcal_mol_A']:.9f} | {v['derivatives_kcal_mol_A']['0.01']['PBE0_D4']:.9f} | {v['derivatives_kcal_mol_A']['0.005']['PBE0_D4']:.9f} |")
    lines+=['','Units kcal/mol/Å; rotation coordinate q=1Å theta. Resolution0.1; both steps and individual/paired checks remain separate. Electronic/D4 derivatives, SCF warnings and energy-stage changes are retained in RESULT. This is sampled numerical consistency, not physical reference accuracy.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'result':rp,'status':status,'stage_numerical_checks_passed':passed,'cost':result['cost']},indent=2))

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--cell',required=True);a.add_argument('--release',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--inventory',required=True);a.add_argument('--stage',choices=('scout','full'),required=True);a.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();preexisting=args.out.exists();tick=time.monotonic();cpu=time.process_time()
    try:
        if args.action=='prepare':prepare(args.out)
        elif args.action=='run':return run(args.panel,args.cell,args.release,args.out)
        else:collect(args.panel,args.inventory,args.stage,args.out)
    except Exception as exc:
        if not preexisting:
            args.out.mkdir(parents=True,exist_ok=True)
            if not (args.out/'FAILURE.json').exists():b.write(args.out/'FAILURE.json',{'status':'failed','implementation':b.pin(__file__),'error':repr(exc),'action':args.action,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu},'native_calls_not_inferred':True})
        raise
    return 0

if __name__=='__main__':sys.exit(main())
