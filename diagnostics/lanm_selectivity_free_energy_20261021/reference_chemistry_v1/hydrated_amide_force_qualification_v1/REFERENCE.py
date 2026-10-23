#!/usr/bin/env python3
"""Eight exact hydrated-amide SPs; preparation/collection are saved-data only."""
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
sp=importlib.util.spec_from_file_location('water48',HERE.parent/'water_force_qualification_v1/REFERENCE.py')
w=importlib.util.module_from_spec(sp);sp.loader.exec_module(w);b=w.b
WS=b.ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
TRAINING=WS/'explicit_model_v1/matched_dry_v1/training_v1/RESULT.json'
TRAINING_SHA='99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb'
WATER_SHA='2871d929c576281fba152de0a7f58069e8c8e88068cd173f970c11fe0952db7c'
BASE_SHA='eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129'
CENTRAL_SHA={'La':'77b64d307a4f6609d43bbd3efc72cfc227668a32965fbc16d9fb518833dee81d',
             'Dy':'459e63b573a78fc68b7fe812ddf834a580057002be78c9fedfd126342e083553'}
PROTOCOL='hydrated_amide_force_qualification_v1';METALS=('La','Dy');SHIFTS=dict(w.SHIFTS);RESOLUTION=.1


def source_gates():
    b.check(b.pin(w.__file__)['sha256']==WATER_SHA and b.pin(b.__file__)['sha256']==BASE_SHA,'Frozen helpers changed')
    b.check(w.RESOLUTION==RESOLUTION and SHIFTS=={'minus_010':-.01,'minus_005':-.005,'plus_005':.005,'plus_010':.01},'Water48 criteria changed')


def metadata(g):
    z=[11 if i==g['metal_index'] else b.Z[s] for i,s in enumerate(g['symbols'])]
    counts=[w.helper.NBF[s] for s in g['symbols']]
    return {'method':'pbe0','d4':True,'role':'complex','metal':g['metal'],'symbols':g['symbols'],'positions_A':g['positions_A'],
        'source_geometry':g['xyz'],'ghost_mask':[False]*g['atoms'],'charge':g['charge'],'multiplicity':1,
        'physical_multiplicity_of_parent':g['physical_multiplicity'],'expected_effective_nuclear_charges':z,
        'expected_electrons':sum(z)-g['charge'],'expected_correlated_electrons':sum(z)-g['charge'],
        'expected_basis_functions':sum(counts),'expected_basis_functions_by_center':counts,'expected_ecp_core':46 if g['metal']=='La' else 55}


def projection(forces,g):
    ids=g['ligand_indices'];u=g['radial_direction']
    ligand=[sum(forces[i][j] for i in ids) for j in range(3)]
    net=[sum(x[j] for x in forces) for j in range(3)]
    return {'moving_indices':ids,'outward_unit_vector':u,'ligand_net_force_kcal_mol_A':ligand,
        'coordinate_derivative_kcal_mol_A':-b.dot(ligand,u),'total_force_kcal_mol_A':net}


def deck(p,c,ranks,maxcore):
    text=b.deck(p,p['central'][c['metal']]['original_cell'],ranks,maxcore)
    b.check(text.count(' EnGrad')==1,'Original EnGrad token differs')
    return text.replace(' EnGrad','')


def parse(text,c,returncode):
    r=w.parse(text,c,returncode)
    r['checks']['single_point_no_gradient_stage']=not r['native_energy_stages']['has_analytic_gradient_stage']
    r['status']='complete' if all(r['checks'].values()) else 'failed_admission'
    return r


def prepare(out):
    tick=time.monotonic();source_gates();b.check(b.pin(TRAINING)['sha256']==TRAINING_SHA,'Training changed')
    data=json.loads(TRAINING.read_text());b.check(data['status']=='training_complete' and data['complete']==28,'Training unavailable')
    out.mkdir(parents=True,exist_ok=False)
    for name in ('assets','xyz','inputs'):(out/name).mkdir()
    assets={};central={};cells=[];controls=[];parent_pin=None;orca=None
    for metal in METALS:
        row=data['rows']['hydrated_acetamide_central_'+metal];g=row['geometry'];saved_pin=row['reference_result']
        b.check(saved_pin['sha256']==CENTRAL_SHA[metal],'Wrong consumed center');saved=json.loads(b.verify(saved_pin).read_text())
        b.check(saved['status']=='complete' and all(saved['checks'].values()),'Central native unavailable')
        parent=b.load_panel(b.verify(saved['panel']));parent_pin=saved['panel'];orca=parent['orca']
        original=next(c for c in parent['cells'] if c['cell_id']==saved['cell_id'])
        bound=original['geometry']
        b.check({k:v for k,v in bound.items() if k not in ('xyz','source_xyz')}=={k:v for k,v in g.items() if k!='xyz'} and
                bound['source_xyz']==g['xyz'] and saved['geometry']==bound['xyz'] and original['recipe']=='ecp_pbe0','Central geometry/recipe differs')
        b.check(b.verify(bound['xyz']).read_bytes()==b.verify(g['xyz']).read_bytes(),'Bound/source geometry bytes differ')
        start_path=Path(saved_pin['path']).parent/'START.json';start=json.loads(start_path.read_text())
        b.check(start['cell']==original and start['panel']==saved['panel'] and start['input']==saved['input'],'Central runtime join differs')
        raw=b.verify(saved['input']).read_text()
        b.check(raw==b.deck(parent,original,start['ranks'],start['maxcore_MiB']),'Original renderer/executed scientific settings mismatch')
        b.check(original['charge']==3 and original['multiplicity']==1 and original['explicit_electrons']==110 and g['atoms']==31,'Center state differs')
        replay=b.collect_one(parent,original,Path(saved_pin['path']).parent,saved['returncode'])
        b.check(replay['status']=='complete' and replay['forces_kcal_mol_A']==saved['forces_kcal_mol_A'] and replay['energy_hartree']==saved['energy_hartree'],'Central E/gradient replay differs')
        b.check(row['reference']['forces']==saved['forces_kcal_mol_A'] and row['reference']['energy_hartree']==saved['energy_hartree'],'Training/native center differs')
        c=metadata(g);strict=w.parse(b.verify(saved['output']).read_text(),c,saved['returncode'])
        b.check(strict['status']=='complete' and strict['energy_hartree']==saved['energy_hartree'],'Strict central native replay failed')
        assets[metal]={}
        for name in ('basis','aux'):
            dst=out/'assets'/f'{metal}_{name}.inc';shutil.copyfile(b.verify(parent['assets'][metal][name]),dst);assets[metal][name]=b.pin(dst)
        central[metal]={'original_cell':original,'saved_result':saved_pin,'saved_start':b.pin(start_path),'executed_input':saved['input'],
            'output':saved['output'],'engrad':saved['engrad'],'original_panel':saved['panel'],'original_adapter':saved['adapter'],
            'energy_hartree':saved['energy_hartree'],'analytic':projection(saved['forces_kcal_mol_A'],g),
            'strict_native_replay':strict,'original_resource_values':{'ranks':start['ranks'],'maxcore_MiB':start['maxcore_MiB']}}
        for label,q in SHIFTS.items():
            x=copy.deepcopy(g['positions_A']);u=g['radial_direction']
            for i in g['ligand_indices']:x[i]=[round(v+q*d,12) for v,d in zip(x[i],u)]
            xyz=out/'xyz'/f'{metal}__{label}.xyz';b.xyz_write(xyz,g['symbols'],x,f'Intact acetamide q={q} A; metal and waters fixed')
            symbols,x=b.xyz_read(xyz);newg=copy.deepcopy(g);newg.update(positions_A=x,xyz=b.pin(xyz))
            c=metadata(newg);c.update(cell_id=f'pbe0__{metal}__{label}__complex',shift_label=label,q_A=q)
            cells.append(c)
        # Saved-output negative controls exercise state/identity/termination admissions.
        text=b.verify(saved['output']).read_text()
        for name,changed,rc in [('nonzero_exit',text,1),('missing_termination',text.replace('ORCA TERMINATED NORMALLY','ORCA TERMINATED ABNORMALLY'),0)]:
            rr=w.parse(changed,metadata(g),rc);b.check(rr['status']=='failed_admission','Failed rejection '+name)
            controls.append({'metal':metal,'case':name,'passed':True})
        wrong=metadata(g);wrong['expected_electrons']+=2
        b.check(w.parse(text,wrong,0)['status']=='failed_admission','Wrong electron target accepted')
        controls.append({'metal':metal,'case':'wrong_expected_electrons','passed':True})
    p={'protocol':PROTOCOL,'created_utc':b.now(),'adapter':b.pin(__file__),'specification':b.pin(HERE/'SPEC.md'),
        'training_result':b.pin(TRAINING),'parent_panel':parent_pin,'orca':orca,'assets':assets,'central':central,'cells':cells,
        'helpers':{'water48':b.pin(w.__file__),'original':b.pin(b.__file__),'donor':b.pin(w.helper.__file__),
            'dft':b.pin(w.dft.__file__),'canonical':b.pin(w.cp.__file__)},
        'maximum_native_calls':8,'ranks':4,'steps_A':[.01,.005],'derivative_resolution_kcal_mol_A':RESOLUTION,
        'coarse_screen_preserved':[v for v in data['radial_derivative_screen'] if v['family']=='acetamide'],
        'physical_state_qualified':False,'method_accuracy_qualified':False,'validation_outcomes_opened':False}
    for c in cells:
        path=out/'inputs'/(c['cell_id']+'.inp');path.write_text(deck(p,c,4,3000));c['template_input']=b.pin(path)
    panel=b.write(out/'PANEL.json',p);val=validate(panel['path']);val.update(saved_failure_controls=controls,
        actual_central_native_replays=2,command=sys.argv,elapsed_seconds=time.monotonic()-tick)
    preflight=b.write(out/'PREFLIGHT.json',val);print(json.dumps({'status':val['status'],'panel':panel,'preflight':preflight,'molecular_calls':0}))


def load(path):
    source_gates();p=json.loads(Path(path).read_text());b.check(p['protocol']==PROTOCOL and p['adapter']==b.pin(__file__),'Protocol/source differs');w.recursive_verify(p)
    expected={f'pbe0__{m}__{s}__complex' for m in METALS for s in SHIFTS}
    b.check(len(p['cells'])==p['maximum_native_calls']==8 and {c['cell_id'] for c in p['cells']}==expected,'Eight-slot denominator differs')
    b.check(p['steps_A']==[.01,.005] and p['ranks']==4 and p['derivative_resolution_kcal_mol_A']==RESOLUTION,'Criteria/resources differ')
    return p


def validate(path):
    p=load(path);records=[]
    for c in p['cells']:
        origin=p['central'][c['metal']];g=origin['original_cell']['geometry'];u=g['radial_direction'];ids=g['ligand_indices']
        b.check(ids==list(range(1,10)) and g['metal_index']==0 and abs(b.norm(u)-1)<1e-14 and b.dot(u,b.sub(g['positions_A'][1],g['positions_A'][0]))>0,'Actual motion/sign differs')
        symbols,x=b.xyz_read(b.verify(c['source_geometry']));b.check(symbols==c['symbols'] and x==c['positions_A'],'XYZ differs')
        b.check(all(x[i]==g['positions_A'][i] for i in range(31) if i not in ids),'Metal/water moved')
        displacement=max(abs(x[i][j]-g['positions_A'][i][j]-c['q_A']*u[j]) for i in ids for j in range(3));b.check(displacement<1e-12,'Ligand translation differs')
        internal=max(abs(b.norm(b.sub(x[i],x[j]))-b.norm(b.sub(g['positions_A'][i],g['positions_A'][j]))) for i in ids for j in ids);b.check(internal<2e-12,'Ligand internal geometry changed')
        original=b.verify(origin['executed_input']).read_text();r=origin['original_resource_values']
        b.check(deck(p,c,r['ranks'],r['maxcore_MiB'])==original.replace(' EnGrad',''),'Scientific deck mismatch')
        raw=deck(p,c,4,3000);b.check(raw==b.verify(c['template_input']).read_text() and 'EnGrad' not in raw and '* xyzfile 3 1 geometry.xyz' in raw,'Prepared deck differs')
        b.check(c['expected_electrons']==110 and c['expected_electrons']==sum(c['expected_effective_nuclear_charges'])-3 and c['expected_basis_functions']==sum(c['expected_basis_functions_by_center']),'State/basis dimensions')
        records.append({'cell_id':c['cell_id'],'q_A':c['q_A'],'atoms':31,'ligand_atoms':len(ids),'electrons':c['expected_electrons'],
            'basis_functions':c['expected_basis_functions'],'maximum_translation_error_A':displacement,'maximum_internal_distance_change_A':internal,'exact_scientific_deck_minus_EnGrad':True})
    for label in SHIFTS:
        a,z=[next(c for c in p['cells'] if c['metal']==m and c['shift_label']==label) for m in METALS]
        b.check(a['positions_A']==z['positions_A'] and a['symbols'][1:]==z['symbols'][1:],'Matched-metal geometry differs')
    return {'status':'passed_static_preflight','panel':b.pin(path),'records':records,'molecular_calls':0,'local_threads':1,'utc':b.now()}


def run(path,cell_id,out):
    p=load(path);b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native work requires Slurm');c=next(c for c in p['cells'] if c['cell_id']==cell_id)
    ranks=int(os.environ.get('SLURM_NTASKS','0'));b.check(ranks==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread MPI ranks required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));b.check(mem>0,'Actual memory grant required');maxcore=math.floor(.75*mem/ranks)
    out.mkdir(parents=True,exist_ok=False);shutil.copyfile(b.verify(c['source_geometry']),out/'geometry.xyz');(out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore));env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[key]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','');cmd=[p['orca']['path'],'endpoint.inp','--bind-to none']
    b.write(out/'START.json',{'cell':c,'panel':b.pin(path),'input':b.pin(out/'endpoint.inp'),'geometry':b.pin(out/'geometry.xyz'),
        'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'ranks':ranks,'memory_MiB':mem,'maxcore_MiB':maxcore,'command':cmd,'utc':b.now()});tick=time.monotonic()
    with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:proc=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False)
    elapsed=time.monotonic()-tick
    try:r=parse((out/'endpoint.out').read_text(errors='replace'),c,proc.returncode)
    except Exception as exc:r={'status':'failed_collection','exception':repr(exc)}
    r.update(cell_id=cell_id,cell=c,panel=b.pin(path),adapter=b.pin(__file__),returncode=proc.returncode,elapsed_seconds=elapsed,
        allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],finished_utc=b.now(),new_native_calls=1,
        output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'),input=b.pin(out/'endpoint.inp'),geometry=b.pin(out/'geometry.xyz'))
    b.write(out/'RESULT.json',r);print(json.dumps({'status':r['status'],'cell_id':cell_id,'seconds':elapsed}));return 0 if r['status']=='complete' else 1


def derivatives(rows,metal,h):
    suffix='010' if h==.01 else '005';a=rows[f'pbe0__{metal}__plus_{suffix}__complex'];z=rows[f'pbe0__{metal}__minus_{suffix}__complex']
    return {key:(a[field]-z[field])*b.HA/(2*h) for key,field in [('electronic','electronic_energy_hartree'),('D4','D4_energy_hartree'),('PBE0_D4','energy_hartree')]}


def collect(path,paths):
    p=load(path);rows={};reviews=[];result={'status':'unavailable','panel':b.pin(path),'adapter':b.pin(__file__),
        'slots':{c['cell_id']:{'status':'unavailable'} for c in p['cells']},'molecular_calls':0,'physical_state_qualified':False,'method_accuracy_qualified':False}
    for filename in paths:
        r=json.loads(Path(filename).read_text());cid=r['cell_id'];b.check(cid in result['slots'] and cid not in rows,'Wrong/duplicate slot');c=next(c for c in p['cells'] if c['cell_id']==cid)
        result['slots'][cid]={'status':r['status'],'result':b.pin(filename)}
        if r['status']!='complete':continue
        b.check(all(r['checks'].values()) and r['panel']==b.pin(path) and r['adapter']==p['adapter'] and r['cell']==c and r['new_native_calls']==1,'Native provenance differs')
        replay=parse(b.verify(r['output']).read_text(),c,r['returncode']);b.check(replay['status']=='complete' and all(replay['checks'].values()),'Native replay failed')
        for key in ('energy_hartree','electronic_energy_hartree','D4_energy_hartree','native_energy_stages'):b.check(replay[key]==r[key],'Stored/native field differs')
        start=json.loads((Path(filename).parent/'START.json').read_text());b.check(start['cell']==c and start['panel']==b.pin(path) and start['input']==r['input'] and start['geometry']==r['geometry'] and start['ranks']==4,'START provenance differs')
        b.check(start['maxcore_MiB']==math.floor(.75*start['memory_MiB']/4) and b.verify(r['input']).read_text()==deck(p,c,4,start['maxcore_MiB']),'Executed deck/resource policy differs')
        b.check(b.verify(r['geometry']).read_bytes()==b.verify(c['source_geometry']).read_bytes(),'Executed geometry differs')
        rows[cid]=r;reviews.append({'result':b.pin(filename),'checks':replay['checks'],'printed_scf_criteria':replay['printed_scf_criteria'],
            'SCF_mode':replay['native_scf_mode'],'SCF_stop_signals':replay['native_scf_stop_signals'],'energy_stages':replay['native_energy_stages']})
    result.update(native_replay_reviews=reviews,admitted_cells=len(rows),expected_cells=8)
    if len(rows)!=8:return result
    per={}
    for metal in METALS:
        center=p['central'][metal];analytic=center['analytic']['coordinate_derivative_kcal_mol_A'];ds={str(h):derivatives(rows,metal,h) for h in p['steps_A']}
        differences={k:ds['0.005'][k]-ds['0.01'][k] for k in ds['0.01']};errors={h:v['PBE0_D4']-analytic for h,v in ds.items()}
        per[metal]={'central_analytic':center['analytic'],'derivatives_kcal_mol_A':ds,'smaller_minus_larger_step':differences,
            'two_step_resolved':{k:abs(v)<=RESOLUTION for k,v in differences.items()},
            'prespecified_smaller_step_summary':{k:ds['0.005'][k] if abs(v)<=RESOLUTION else None for k,v in differences.items()},
            'PBE0_D4_FD_minus_central_analytic':errors,'sampled_analytic_consistency_at_both_steps':all(abs(v)<=RESOLUTION for v in errors.values()),
            'energy_stage_comparisons':{s:w.compare_stages(center['strict_native_replay']['native_energy_stages'],rows[f'pbe0__{metal}__{s}__complex']['native_energy_stages']) for s in SHIFTS}}
    paired={str(h):{k:per['Dy']['derivatives_kcal_mol_A'][str(h)][k]-per['La']['derivatives_kcal_mol_A'][str(h)][k] for k in ('electronic','D4','PBE0_D4')} for h in p['steps_A']}
    analytic=p['central']['Dy']['analytic']['coordinate_derivative_kcal_mol_A']-p['central']['La']['analytic']['coordinate_derivative_kcal_mol_A']
    changes={k:paired['0.005'][k]-paired['0.01'][k] for k in paired['0.01']};errors={h:v['PBE0_D4']-analytic for h,v in paired.items()}
    result.update(status='complete',per_metal=per,Dy_minus_La={'central_analytic_derivative_kcal_mol_A':analytic,'derivatives_kcal_mol_A':paired,
        'smaller_minus_larger_step':changes,'two_step_resolved':{k:abs(v)<=RESOLUTION for k,v in changes.items()},
        'prespecified_smaller_step_summary':{k:paired['0.005'][k] if abs(v)<=RESOLUTION else None for k,v in changes.items()},
        'PBE0_D4_FD_minus_central_analytic':errors,'sampled_analytic_consistency_at_both_steps':all(abs(v)<=RESOLUTION for v in errors.values())},
        coarse_screen_preserved=p['coarse_screen_preserved'],derivative_resolution_kcal_mol_A=RESOLUTION,no_extrapolation=True,no_additional_step=True)
    return result


def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True);a=sub.add_parser('prepare');a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('validate');a.add_argument('--panel',required=True);a.add_argument('--out',required=True)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--cell',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--out',required=True);a.add_argument('results',nargs='+');a=ap.parse_args()
    if a.action=='prepare':prepare(a.out)
    elif a.action=='validate':b.write(a.out,validate(a.panel))
    elif a.action=='run':return run(a.panel,a.cell,a.out)
    else:b.write(a.out,collect(a.panel,a.results))

if __name__=='__main__':sys.exit(main())
