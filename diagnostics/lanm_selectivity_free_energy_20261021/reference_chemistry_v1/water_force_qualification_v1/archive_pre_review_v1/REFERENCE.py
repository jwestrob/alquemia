#!/usr/bin/env python3
"""Finite rigid-water CC/PBE0 derivative panel; native work requires Slurm."""
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
sp=importlib.util.spec_from_file_location('donor_cp',HERE.parent/'donor_cp_v1/REFERENCE.py');helper=importlib.util.module_from_spec(sp);sp.loader.exec_module(helper)
b=helper.b;cp=helper.cp;dft=helper.dft;ROOT=b.ROOT
WS=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
PARENT=WS/'reference_chemistry_v1/dft_water_cp_v1/prepared_v3/PANEL.json'
PARENT_SHA='3efe421275c6e9833eea757b6b0990e04195283a571326f65c40708ece1f1d0e'
BASELINE=WS/'reference_chemistry_v1/dft_water_cp_v1/final_collection_v1/RESULT.json'
BASELINE_SHA='6465f8c9002c31ee9c79a945152476c91ad740da0a5ddca4e2327d98321bda64'
AMOEBA=WS/'explicit_model_v1/controls_v1/RESULT.json'
AMOEBA_SHA='f2d70fc06fad5ae19566c9fd779b14cbcf48f170614ccae24320a482c4c0deaf'
PROTOCOL='matched_water_force_qualification_v1'
METALS=('La','Dy');METHODS=('cc','pbe0');ROLES=('complex','ion_ghost_water','water_ghost_ion')
SHIFTS={'minus_010':-.01,'minus_005':-.005,'plus_005':.005,'plus_010':.01}
RESOLUTION=.1
deck=helper.deck

def norm(x):return math.sqrt(sum(v*v for v in x))
def dot(x,y):return sum(a*z for a,z in zip(x,y))
def subtract(x,y):return [a-z for a,z in zip(x,y)]
def recursive_verify(obj):
    if isinstance(obj,dict):
        if set(obj)=={'path','sha256'}:b.verify(obj)
        else:
            for value in obj.values():recursive_verify(value)
    elif isinstance(obj,list):
        for value in obj:recursive_verify(value)

def energy_stage_metadata(text):
    native=text.split('****END OF INPUT****')[-1];lines=native.splitlines()
    selected=[{'line':i+1,'text':line} for i,line in enumerate(lines) if re.search(r'COSX|gridx[123]|final integration|change in energy|Total Energy\s*:|FINAL SINGLE POINT|Dispersion correction|Energy\+Gradient|SCF GRADIENT|SCF CONVERGED',line,re.I)]
    grid=native.split('COSX GRID GENERATION',1)[1].split('SCF SETTINGS',1)[0][:10000] if 'COSX GRID GENERATION' in native else ''
    return {'selected_native_lines':selected,'COSX_grid_generation_excerpt':grid,
        'recomputed_final_exchange':bool(re.search(r'Recomputing exchange energy using gridx3',native)),
        'final_integration_energy_hartree':[cp.number(v) for v in re.findall(r'Total energy after final integration\s*:\s*('+cp.FLOAT+')',native)],
        'has_analytic_gradient_stage':'ORCA SCF GRADIENT CALCULATION' in native,
        'surface_equivalence_assumed':False}

def parse(text,c,returncode):
    r=helper.parse(text,c,returncode)
    if c['method']=='pbe0':r['native_energy_stages']=energy_stage_metadata(text)
    return r

def analytic_record(forces,u):
    b.check(len(forces)==4 and all(len(x)==3 and all(math.isfinite(v) for v in x) for x in forces),'Invalid central force')
    water=[sum(row[j] for row in forces[1:]) for j in range(3)];net=[sum(row[j] for row in forces) for j in range(3)]
    return {'coordinate_derivative_kcal_mol_A':-dot(water,u),'water_radial_force_kcal_mol_A':dot(water,u),
        'metal_radial_force_kcal_mol_A':dot(forces[0],u),'net_force_kcal_mol_A':net,'net_force_norm_kcal_mol_A':norm(net)}

def prepare(out):
    tick=time.monotonic()
    for path,sha in ((PARENT,PARENT_SHA),(BASELINE,BASELINE_SHA),(AMOEBA,AMOEBA_SHA)):b.check(b.pin(path)['sha256']==sha,'Declared saved source changed')
    old=dft.load(PARENT);canonical=cp.load(b.verify(old['canonical_panel']));baseline=json.loads(BASELINE.read_text());amoeba=json.loads(AMOEBA.read_text())
    out.mkdir(parents=True,exist_ok=False)
    for name in ('assets','inputs','xyz'):(out/name).mkdir()
    assets={};central={};geometries=[];cells=[]
    for metal in METALS:
        assets[metal]={}
        for key in ('basis','aux'):
            dst=out/'assets'/f'{metal}_{key}.inc';shutil.copyfile(b.verify(old['assets'][metal][key]),dst);assets[metal][key]=b.pin(dst)
        origin=next(x for x in canonical['cells'] if x['cell_id']==metal+'__complex');x=origin['positions_A'];v=subtract(x[1],x[0]);length=norm(v);u=[a/length for a in v]
        saved=old['saved_complexes'][metal];r=json.loads(b.verify(saved['result']).read_text());b.check(r['status']=='complete' and all(r['checks'].values()),'Saved PBE0 unavailable')
        b.check(r['geometry']==origin['source_geometry'],'Central PBE0 geometry differs');recursive_verify(saved)
        ar=next(row for row in amoeba['rows'] if row['id']=='water_central_'+metal);b.check(ar['status']=='complete','Saved AMOEBA unavailable');g=ar['input']['source_geometry']
        b.check(g['symbols']==origin['symbols'] and g['positions_A']==x,'AMOEBA/PBE0 central coordinates differ');recursive_verify(ar['input']);recursive_verify(ar['native'])
        b.check(all(abs(f+gr)<1e-12 for fr,grr in zip(ar['native']['force_kcal_mol_A'],ar['native']['gradient_kcal_mol_A']) for f,gr in zip(fr,grr)),'AMOEBA force sign differs')
        central[metal]={'origin_cell':origin,'unit_metal_to_O':u,'metal_O_distance_A':length,'PBE0_saved':saved,'PBE0_analytic':analytic_record(r['forces_kcal_mol_A'],u),
            'PBE0_native_energy_stages':energy_stage_metadata(b.verify(saved['output']).read_text()),'PBE0_engrad':r['engrad'],'AMOEBA_input':ar['input'],
            'AMOEBA_native':ar['native'],'AMOEBA_analytic':analytic_record(ar['native']['force_kcal_mol_A'],u)}
        for label,q in SHIFTS.items():
            positions=[list(x[0])]+[[round(a+q*z,12) for a,z in zip(row,u)] for row in x[1:]]
            xyz=out/'xyz'/f'{metal}__{label}.xyz';b.xyz_write(xyz,origin['symbols'],positions,f'Rigid water q={q} A; metal fixed')
            geometry={'geometry_id':metal+'__'+label,'metal':metal,'shift_label':label,'q_A':q,'unit_metal_to_O':u,'symbols':origin['symbols'],'positions_A':positions,'xyz':b.pin(xyz)};geometries.append(geometry)
            for method in METHODS:
                for role in ROLES:
                    c=copy.deepcopy(next(row for row in canonical['cells'] if row['cell_id']==metal+'__'+role));c.pop('template_input')
                    c.update(cell_id=f'{method}__{metal}__{label}__{role}',method=method,d4=(method=='pbe0' and role=='complex'),positions_A=positions,
                        source_geometry=geometry['xyz'],origin_geometry=origin['source_geometry'],geometry_id=geometry['geometry_id'],q_A=q,shift_label=label)
                    cells.append(c)
    p={'protocol':PROTOCOL,'created_utc':b.now(),'adapter':b.pin(__file__),'specification':b.pin(HERE/'SPEC.md'),
        'helpers':{'donor':b.pin(HERE.parent/'donor_cp_v1/REFERENCE.py'),'canonical':old['canonical_helper'],'dft':old['adapter'],'base':old['base_helper']},
        'parent_panel':b.pin(PARENT),'water_baseline':b.pin(BASELINE),'canonical_collection':baseline['canonical_collection'],'AMOEBA_control':b.pin(AMOEBA),
        'orca':old['orca'],'assets':assets,'central':central,'geometries':geometries,'cells':cells,'maximum_native_calls':48,'ranks':4,
        'scout_cell_id':'cc__La__plus_010__complex','derivative_resolution_kcal_mol_A':RESOLUTION,'steps_A':[.01,.005],
        'free_water_D4_constant_cancels':True,'physical_state_qualified':False,'method_accuracy_qualified':False}
    for c in cells:
        path=out/'inputs'/(c['cell_id']+'.inp');path.write_text(deck(p,c,4,3000));c['template_input']=b.pin(path)
    pin=b.write(out/'PANEL.json',p);v=validate(pin['path']);v.update(command=sys.argv,elapsed_seconds=time.monotonic()-tick);b.write(out/'PREFLIGHT.json',v)
    print(json.dumps({'status':v['status'],'panel':pin,'cells':len(cells),'molecular_calls':0}))

def load(path):
    p=json.loads(Path(path).read_text());b.check(p['protocol']==PROTOCOL and p['adapter']==b.pin(__file__),'Wrong or changed adapter');recursive_verify(p)
    b.check(p['maximum_native_calls']==len(p['cells'])==48 and len(p['geometries'])==8,'Finite denominator changed')
    expected={f'{method}__{metal}__{label}__{role}' for method in METHODS for metal in METALS for label in SHIFTS for role in ROLES}
    b.check({c['cell_id'] for c in p['cells']}==expected,'Inventory differs');b.check(p['steps_A']==[.01,.005] and p['derivative_resolution_kcal_mol_A']==RESOLUTION,'Resolution/steps differ')
    return p

def validate(path):
    p=load(path);records=[]
    for g in p['geometries']:
        origin=p['central'][g['metal']]['origin_cell'];x=origin['positions_A'];u=p['central'][g['metal']]['unit_metal_to_O'];v=subtract(x[1],x[0]);b.check(abs(norm(u)-1)<1e-14 and dot(u,v)>0,'Direction normalization/sign')
        symbols,actual=b.xyz_read(g['xyz']['path']);b.check(symbols==g['symbols'] and actual==g['positions_A'],'XYZ differs');b.check(actual[0]==x[0],'Metal moved')
        b.check(max(abs((actual[i][j]-x[i][j])-g['q_A']*u[j]) for i in (1,2,3) for j in range(3))<1e-12,'Water is not rigid translation')
        b.check(max(abs(norm(subtract(actual[i],actual[j]))-norm(subtract(x[i],x[j]))) for i,j in ((1,2),(1,3),(2,3)))<2e-12,'Water internal geometry changed')
        b.check(abs(norm(subtract(actual[1],actual[0]))-norm(v)-g['q_A'])<1e-12,'Displacement sign/distance differs')
    for label in SHIFTS:
        a,z=[next(g for g in p['geometries'] if g['geometry_id']==m+'__'+label) for m in METALS];b.check(a['positions_A']==z['positions_A'] and a['symbols'][1:]==z['symbols'][1:],'Paired geometry differs')
    for c in p['cells']:
        raw=deck(p,c,4,3000);b.check(raw==b.verify(c['template_input']).read_text(),'Deck differs');g=next(g for g in p['geometries'] if g['geometry_id']==c['geometry_id']);b.check(c['positions_A']==g['positions_A'] and c['source_geometry']==g['xyz'],'Cell geometry differs')
        header,*lines=raw.split('* xyz ',1)[1].splitlines();lines=lines[:-1];labels=[line.split()[0] for line in lines];charge,mult=map(int,header.split());b.check(len(labels)==4 and [s.endswith(':') for s in labels]==c['ghost_mask'],'Center/ghost inventory')
        actualZ=[0 if s.endswith(':') else (11 if i==0 else b.Z[s]) for i,s in enumerate(labels)];ne={'complex':18,'ion_ghost_water':8,'water_ghost_ion':10}[c['role']]
        b.check(actualZ==c['expected_effective_nuclear_charges'] and sum(actualZ)-charge==ne==c['expected_electrons']==c['expected_correlated_electrons'] and mult==1,'Electron/state inventory')
        b.check(c['expected_basis_functions']==90 and c['expected_basis_functions_by_center']==[47,31,6,6],'AO inventory')
        b.check(('NewECP' in raw)==bool(c['expected_ecp_core']),'ECP inventory')
        if c['ghost_mask'][0]:b.check('DelECP '+c['metal'] in raw and 'NewGTO '+c['metal'] in raw,'Ghost ECP/orbital basis')
        if c['method']=='cc':b.check('NoFrozenCore' in raw and 'STol 1e-8' in raw and not any(k in raw for k in ('D4','RIJCOSX','EnGrad')),'CC method')
        else:b.check(('D4' in raw)==(c['role']=='complex') and 'RIJCOSX' in raw and 'DefGrid3' in raw and 'EnGrad' not in raw,'DFT method/D4')
        records.append({'cell_id':c['cell_id'],'charge':charge,'electrons':ne,'basis_functions':90,'q_A':c['q_A'],'ghost_centers':sum(c['ghost_mask'])})
    return {'status':'passed_static_preflight','panel':b.pin(path),'independent_deck_inventory':records,'molecular_calls':0,'local_threads':1,'utc':b.now()}

def run(path,cell_id,out):
    p=load(path);b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native work requires Slurm');c=next(c for c in p['cells'] if c['cell_id']==cell_id)
    ranks=int(os.environ.get('SLURM_NTASKS','0'));b.check(ranks==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread MPI ranks required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));b.check(mem>0,'Actual memory grant required');maxcore=math.floor(.75*mem/ranks)
    out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore));env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[key]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','');cmd=[p['orca']['path'],'endpoint.inp','--bind-to none']
    b.write(out/'START.json',{'cell':c,'panel':b.pin(path),'input':b.pin(out/'endpoint.inp'),'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'ranks':ranks,'memory_MiB':mem,'maxcore_MiB':maxcore,'command':cmd,'utc':b.now()});tick=time.monotonic()
    with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:proc=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False)
    elapsed=time.monotonic()-tick
    try:r=parse((out/'endpoint.out').read_text(errors='replace'),c,proc.returncode)
    except Exception as exc:r={'status':'failed_collection','exception':repr(exc)}
    r.update(cell_id=cell_id,cell=c,panel=b.pin(path),adapter=b.pin(__file__),returncode=proc.returncode,elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],finished_utc=b.now(),output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'),input=b.pin(out/'endpoint.inp'))
    b.write(out/'RESULT.json',r);print(json.dumps({'status':r['status'],'cell_id':cell_id,'seconds':elapsed}));return 0 if r['status']=='complete' else 1

def energy_terms(rows,method,metal,label):
    trio=[rows[f'{method}__{metal}__{label}__{role}'] for role in ROLES];key='energy_hartree' if method=='cc' else 'electronic_energy_hartree'
    result={'complex':trio[0][key],'CP':trio[0][key]-trio[1][key]-trio[2][key]}
    if method=='pbe0':result.update(D4=trio[0]['D4_energy_hartree'],complex_D4=result['complex']+trio[0]['D4_energy_hartree'],CP_D4=result['CP']+trio[0]['D4_energy_hartree'])
    return result

def derivatives(rows,method,metal,step):
    suffix='010' if step==.01 else '005';plus=energy_terms(rows,method,metal,'plus_'+suffix);minus=energy_terms(rows,method,metal,'minus_'+suffix)
    return {key:(plus[key]-minus[key])*b.HA/(2*step) for key in plus}

def collect(path,paths):
    p=load(path);rows={};reviews=[]
    for filename in paths:
        r=json.loads(Path(filename).read_text());c=next(c for c in p['cells'] if c['cell_id']==r['cell_id'])
        b.check(r['status']=='complete' and all(r['checks'].values()) and r['panel']==b.pin(path) and r['adapter']==p['adapter'] and r['cell']==c,'Native provenance/admission differs');b.check(r['cell_id'] not in rows,'Duplicate result')
        replay=parse(b.verify(r['output']).read_text(),c,r['returncode']);b.check(replay['status']=='complete' and all(replay['checks'].values()),'Native replay failed')
        for key in ('energy_hartree','electronic_energy_hartree','D4_energy_hartree'):
            if key in replay:b.check(replay[key]==r[key],'Stored/native energy differs')
        start=json.loads((Path(filename).parent/'START.json').read_text());b.check(start['cell']==c and start['panel']==b.pin(path) and start['ranks']==4,'START provenance differs');b.check(b.verify(r['input']).read_text()==deck(p,c,4,start['maxcore_MiB']),'Executed deck differs')
        reviews.append({'result':b.pin(filename),'checks':replay['checks'],'energy_stages':replay.get('native_energy_stages'),'SCF_mode':replay['native_scf_mode'],'SCF_stop_signals':replay['native_scf_stop_signals'],'printed_scf_criteria':replay['printed_scf_criteria']});rows[r['cell_id']]=r
    b.check(set(rows)=={c['cell_id'] for c in p['cells']},'Complete48 denominator required');results={}
    for metal in METALS:
        central=p['central'][metal];a=central['AMOEBA_analytic']['coordinate_derivative_kcal_mol_A'];native=central['PBE0_analytic']['coordinate_derivative_kcal_mol_A'];methods={}
        for method in METHODS:
            perstep={str(h):derivatives(rows,method,metal,h) for h in p['steps_A']};differences={key:perstep['0.005'][key]-perstep['0.01'][key] for key in perstep['0.01']}
            methods[method]={'derivatives_kcal_mol_A':perstep,'smaller_minus_larger_step':differences,'two_step_resolved':{k:abs(v)<=RESOLUTION for k,v in differences.items()},
                'AMOEBA_minus_derivative':{h:{key:a-value for key,value in ds.items() if key!='D4'} for h,ds in perstep.items()}}
        pbe=methods['pbe0'];pbe['complex_D4_FD_minus_central_analytic']={h:ds['complex_D4']-native for h,ds in pbe['derivatives_kcal_mol_A'].items()}
        pbe['sampled_analytic_consistency_at_both_steps']=all(abs(v)<=RESOLUTION for v in pbe['complex_D4_FD_minus_central_analytic'].values())
        pbe['central_energy_stage_metadata']=central['PBE0_native_energy_stages'];pbe['SP_vs_EnGrad_surface_equivalence_assumed']=False
        results[metal]={'central_AMOEBA_analytic':central['AMOEBA_analytic'],'central_PBE0_D4_analytic':central['PBE0_analytic'],'methods':methods}
    paired={method:{str(h):{key:results['Dy']['methods'][method]['derivatives_kcal_mol_A'][str(h)][key]-results['La']['methods'][method]['derivatives_kcal_mol_A'][str(h)][key] for key in results['La']['methods'][method]['derivatives_kcal_mol_A'][str(h)]} for h in p['steps_A']} for method in METHODS}
    return {'status':'complete','panel':b.pin(path),'native_replay_reviews':reviews,'per_metal':results,'Dy_minus_La_derivatives':paired,'derivative_resolution_kcal_mol_A':RESOLUTION,
        'molecular_calls':0,'physical_state_qualified':False,'method_accuracy_qualified':False,'no_force_field_fit':True}

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
