#!/usr/bin/env python3
"""Eight declared electronic metal/whole-shell CP single points; root submits."""
import argparse
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

HERE = Path(__file__).resolve().parent
sp = importlib.util.spec_from_file_location('coordination_donor_helper', HERE.parent/'donor_cp_v1/REFERENCE.py')
donor = importlib.util.module_from_spec(sp); sp.loader.exec_module(donor)
b, dft, cp = donor.b, donor.dft, donor.cp
PROTOCOL = 'coordination_metal_whole_shell_cp_v1'
GIDS = [f'aquo{n}_central_{m}' for n in (8, 9) for m in ('La', 'Dy')]
ROLES = ('ion_ghost_shell', 'shell_ghost_ion')
HELPERS = {k: HERE.parent/v for k, v in {
    'base_helper':'REFERENCE.py', 'canonical_helper':'canonical_water_cp_v2/REFERENCE.py',
    'dft_helper':'dft_water_cp_v1/REFERENCE.py', 'donor_helper':'donor_cp_v1/REFERENCE.py'}.items()}

def read(pin): return json.loads(b.verify(pin).read_text())

def make_cell(g, role):
    n = g['water_count']; m = g['metal']; count = 1+3*n
    b.check(n in (8, 9) and g['symbols']==[m]+['O','H','H']*n, 'Whole water-shell inventory differs')
    ghost = [False]+[True]*(count-1) if role==ROLES[0] else [True]+[False]*(count-1)
    charge = 3 if role==ROLES[0] else 0
    z = [0 if ghost[i] else (11 if i==0 else b.Z[s]) for i,s in enumerate(g['symbols'])]
    return dict(cell_id=f'pbe0__aquo{n}__{m}__{role}', geometry_id=g['geometry_id'], water_count=n,
        method='pbe0', metal=m, role=role, d4=False, symbols=g['symbols'], positions_A=g['positions_A'],
        source_geometry=g['xyz'], ghost_mask=ghost, charge=charge, multiplicity=1,
        physical_multiplicity_of_parent=g['physical_multiplicity'], expected_effective_nuclear_charges=z,
        expected_electrons=sum(z)-charge, expected_correlated_electrons=sum(z)-charge,
        expected_basis_functions=47+43*n, expected_basis_functions_by_center=[47]+[31,6,6]*n,
        expected_ecp_core=(46 if m=='La' else 55) if not ghost[0] else 0)

def deck(p, c, ranks, maxcore): return donor.deck(p, c, ranks, maxcore)

def parse(text, c, returncode):
    r = donor.parse(text, c, returncode)
    native = text.split('****END OF INPUT****')[-1]
    r['energy_stage_metadata'] = [{'line':i+1, 'text':line.strip()} for i,line in enumerate(text.splitlines())
        if re.search(r'final.*(?:grid|COSX|exchange)|COSX.*(?:correction|energy)|FinalGrid|GridX|Grid[123]|SCF CONVERGED AFTER|signals convergence', line, re.I)]
    r['final_COSX_correction_hartree'] = [cp.number(v) for v in re.findall(r'Final COSX correction\s*[:.]?\s*('+cp.FLOAT+r')', native, re.I)]
    return r

def prepare(config, out):
    tick, cpu = time.monotonic(), time.process_time(); cfg = json.loads(Path(config).read_text())
    b.check(cfg['implementation']==b.pin(__file__), 'Preparation source differs')
    for key,path in HELPERS.items(): b.check(cfg[key]==b.pin(path), 'Helper differs: '+key)
    root = read(cfg['root_inputs']); b.verify(root['declaration']); b.verify(root['prior_alpha_admission'])
    b.check(root['source_cells']==GIDS and root['new_reference_cells']==8 and root['new_native_calls_authorized_by_this_file'] is False, 'Root finite scope differs')
    data = read(root['consumed_input_data']); old = b.load_panel(cp.PARENT)
    b.check(b.pin(cp.PARENT)['sha256']==cp.PARENT_SHA, 'Original reference panel differs')
    out.mkdir(parents=True, exist_ok=False); (out/'assets').mkdir(); (out/'inputs').mkdir()
    assets = {}
    for m in ('La','Dy'):
        assets[m] = {}
        for key in ('basis','aux'):
            dst=out/'assets'/f'{m}_{key}.inc'; shutil.copyfile(b.verify(old['assets'][m][key]),dst); assets[m][key]=b.pin(dst)
    p = dict(protocol=PROTOCOL, adapter=b.pin(__file__), configuration=b.pin(config), root_inputs=cfg['root_inputs'],
        declaration=root['declaration'], consumed_input_data=root['consumed_input_data'], parent_panel=b.pin(cp.PARENT),
        orca=old['orca'], assets=assets, source_assets={m:{k:old['assets'][m][k] for k in ('basis','aux')} for m in ('La','Dy')},
        geometries={}, saved_complexes={}, baseline_model_energies={}, cells=[], maximum_native_calls=8, ranks=4,
        run_wrapper=b.pin(HERE/'RUN.sbatch'), interpretation='Two-fragment metal|whole-shell BSSE only; original targets unchanged; no force/state/basis-limit/free-energy qualification',
        created_utc=b.now(), **{k:cfg[k] for k in HELPERS})
    replay = {}
    for gid in GIDS:
        row=data['rows'][gid]; g=row['geometry']; p['geometries'][gid]=g
        src=row['reference_result']; b.check(src==row['reference'].get('raw_result',src), 'Reference authority differs')
        r=read(src); oc=next(c for c in old['cells'] if c['cell_id']=='ecp_pbe0__'+gid)
        b.check(oc['geometry']==g and r['panel']==p['parent_panel'] and r['adapter']==p['base_helper'], 'Original cell/panel join differs')
        b.check(r['cell_id']==oc['cell_id'] and r['status']=='complete' and all(r['checks'].values()), 'Original complex not admitted')
        for key in ('output','input','engrad','geometry'): b.verify(r[key])
        b.check(b.verify(r['input']).read_text()==b.deck(old,oc,4,3000), 'Original executed recipe differs')
        q=b.collect_one(old,oc,Path(src['path']).parent,r['returncode'])
        b.check(all(q[k]==r[k] for k in q), 'Original raw complex replay differs')
        b.check(r['energy_hartree']*b.HA==row['reference']['energy'] and r['forces_kcal_mol_A']==row['reference']['forces'], 'Consumed reference differs')
        es,ed,ef=dft.electronic_parts(b.verify(r['output']).read_text())
        b.check(len(es)==len(ed)==len(ef)==1 and abs(es[0]+ed[0]-ef[0])<1e-8 and ef[0]==r['energy_hartree'], 'Original electronic/D4 closure differs')
        p['saved_complexes'][gid]=dict(result=src, output=r['output'], input=r['input'], engrad=r['engrad'], source_geometry=r['geometry'],
            electronic_energy_hartree=es[0], D4_energy_hartree=ed[0], final_energy_hartree=ef[0], job_id=r['job_id'])
        p['baseline_model_energies'][gid]=data['baseline_model'][gid]['energy']
        p['cells'].extend(make_cell(g,role) for role in ROLES); replay[gid]=True
    for n in (8,9):
        a,z=(p['geometries'][f'aquo{n}_central_{m}'] for m in ('La','Dy'))
        b.check(a['positions_A']==z['positions_A'] and a['symbols'][1:]==z['symbols'][1:], 'Own-shell cancellation invalid')
    for m in ('La','Dy'):
        a,z=(p['saved_complexes'][f'aquo{n}_central_{m}']['input'] for n in (8,9))
        b.check(a['sha256']==z['sha256'], 'Own-metal recipe cancellation invalid')
        aa,zz=(p['geometries'][f'aquo{n}_central_{m}'] for n in (8,9))
        b.check(aa['positions_A'][0]==zz['positions_A'][0]==[0.,0.,0.] and aa['charge']==zz['charge']==3 and aa['physical_multiplicity']==zz['physical_multiplicity'], 'Own-ion state/position differs')
    p['cancellation_checks']={'paired_coordinates_and_light_basis':True,'same_metal_recipe_across_sizes':True,
        'own_ion':{'same_position':[0.,0.,0.],'same_effective_charge':3,'same_effective_multiplicity':1,'same_per_metal_basis_aux_ECP':True},
        'own_shell':{'same_order_and_coordinates_between_metals':True,'light_basis':'def2-TZVP','auxiliary_basis':'def2/J','charge':0,'multiplicity':1,'whole_interacting_shell':True}}
    for c in p['cells']:
        f=out/'inputs'/(c['cell_id']+'.inp'); f.write_text(deck(p,c,4,3000)); c['template_input']=b.pin(f)
    pin=b.write(out/'PANEL.json',p); val=validate(pin['path'])
    val.update(original_raw_replays=replay, command=sys.argv, wall_seconds=time.monotonic()-tick, parent_cpu_seconds=time.process_time()-cpu)
    b.write(out/'PREFLIGHT.json',val); print(json.dumps({'panel':pin,'status':val['status'],'new_native_calls':0}))

def load(path):
    p=json.loads(Path(path).read_text()); b.check(p['protocol']==PROTOCOL and p['adapter']==b.pin(__file__), 'Protocol/source differs')
    for key in ('adapter','configuration','root_inputs','declaration','consumed_input_data','parent_panel','orca','run_wrapper',*HELPERS): b.verify(p[key])
    for key,src in HELPERS.items(): b.check(p[key]==b.pin(src), 'Executing helper differs')
    b.check(len(p['cells'])==p['maximum_native_calls']==8 and p['ranks']==4 and set(p['geometries'])==set(GIDS), 'Finite inventory differs')
    b.check({c['cell_id'] for c in p['cells']}=={make_cell(g,r)['cell_id'] for g in p['geometries'].values() for r in ROLES}, 'Cell IDs differ')
    for m,a in p['assets'].items():
        for key,pin in a.items(): b.verify(pin); b.verify(p['source_assets'][m][key]); b.check(pin['sha256']==p['source_assets'][m][key]['sha256'],'Copied asset differs')
    for r in p['saved_complexes'].values():
        for key in ('result','output','input','engrad','source_geometry'): b.verify(r[key])
    for c in p['cells']:
        expected=make_cell(p['geometries'][c['geometry_id']],c['role'])
        b.check({k:v for k,v in c.items() if k!='template_input'}==expected,'Fragment definition differs')
        b.verify(c['source_geometry']); b.verify(c['template_input'])
        b.check(b.verify(c['template_input']).read_text()==deck(p,c,4,3000),'Frozen template differs')
    return p

def validate(path):
    p=load(path); rows=[]
    for c in p['cells']:
        text=b.verify(c['template_input']).read_text(); sy,x=b.xyz_read(c['source_geometry']['path'])
        b.check(sy==c['symbols'] and x==c['positions_A'], 'Exact source geometry differs')
        header,*lines=text.split('* xyz ',1)[1].splitlines(); lines=lines[:-1]; q,m=map(int,header.split())
        labels=[line.split()[0] for line in lines]; actual_z=[0 if s.endswith(':') else (11 if i==0 else b.Z[s]) for i,s in enumerate(labels)]
        b.check(q==c['charge'] and m==1 and actual_z==c['expected_effective_nuclear_charges'], 'Independent charge/ghost inventory differs')
        b.check(sum(actual_z)-q==c['expected_electrons']==(8 if c['role']==ROLES[0] else 10*c['water_count']), 'Electron inventory differs')
        b.check([s.endswith(':') for s in labels]==c['ghost_mask'] and len(labels)==1+3*c['water_count'], 'Ghost mask differs')
        b.check(sum(c['expected_basis_functions_by_center'])==c['expected_basis_functions']==47+43*c['water_count'], 'Basis inventory differs')
        b.check('D4' not in text and 'EnGrad' not in text and 'MaxIter 300' in text and 'DefGrid3' in text,'Single-point recipe differs')
        b.check(('NewECP' in text)==bool(c['expected_ecp_core']) and ('DelECP '+c['metal'] in text)==c['ghost_mask'][0], 'ECP/ghost role differs')
        b.check(cp.metal_blocks(b.verify(p['assets'][c['metal']]['basis']).read_text(),c['metal'])[1] in text and b.verify(p['assets'][c['metal']]['aux']).read_text() in text,'Metal orbital/aux basis differs')
        rows.append({k:c[k] for k in ('cell_id','charge','expected_electrons','expected_basis_functions','expected_ecp_core')})
    return dict(status='passed_static_preflight',panel=b.pin(path),independent_deck_inventory=rows,new_native_calls=0,local_threads=1,
        native_ghost_semantics_admitted=False,utc=b.now())

def gate(p, path, release, cell_id, out):
    r=json.loads(Path(release).read_text())
    b.check(r['status']=='approved_coordination_cp_eight_cells' and r['implementation']==p['adapter'] and r['panel']==b.pin(path),'Root release differs')
    b.check(r['cell_ids']==[c['cell_id'] for c in p['cells']] and r['maximum_native_calls']==8 and r['ranks']==4,'Released finite scope differs')
    pre=read(r['preparation']); b.verify(r['independent_review'])
    b.check(pre['status']=='passed_static_preflight' and pre['panel']==b.pin(path),'Preparation not admitted')
    b.check(out.resolve()==Path(r['execution_root']).resolve()/'tasks'/cell_id/'execution','One-attempt cell path differs')
    b.check(not out.exists(), 'Existing attempt must be preserved')
    return r

def run(path, release, cell_id, out):
    p=load(path); matches=[c for c in p['cells'] if c['cell_id']==cell_id]; b.check(len(matches)==1,'Unknown cell'); c=matches[0]
    gate(p,path,release,cell_id,out); b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native calls require Slurm')
    ranks=int(os.environ.get('SLURM_NTASKS','0')); b.check(ranks==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread ranks required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0')); b.check(mem>0,'Actual memory grant required')
    maxcore=math.floor(.75*mem/ranks); b.check(maxcore>0,'Invalid per-rank memory')
    out.mkdir(parents=True,exist_ok=False); (out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore))
    env=os.environ.copy()
    for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'): env[k]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE'); env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','')
    cmd=[p['orca']['path'],'endpoint.inp','--bind-to none']
    b.write(out/'START.json',dict(cell=c,panel=b.pin(path),release=b.pin(release),input=b.pin(out/'endpoint.inp'),job_id=os.environ['SLURM_JOB_ID'],
        host=socket.gethostname(),ranks=ranks,memory_MiB=mem,maxcore_MiB=maxcore,command=cmd,maximum_native_calls=1,utc=b.now()))
    tick=time.monotonic(); code=None
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:
            proc=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False); code=proc.returncode
        r=parse((out/'endpoint.out').read_text(errors='replace'),c,code)
    except Exception as exc: r=dict(status='failed_execution_or_collection',exception=repr(exc))
    elapsed=time.monotonic()-tick
    r.update(cell_id=cell_id,cell=c,panel=b.pin(path),adapter=p['adapter'],release=b.pin(release),returncode=code,elapsed_seconds=elapsed,
        allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],finished_utc=b.now(),input=b.pin(out/'endpoint.inp'),native_program_invocations_attempted=1,
        output=b.pin(out/'endpoint.out') if (out/'endpoint.out').exists() else None,stderr=b.pin(out/'endpoint.err') if (out/'endpoint.err').exists() else None)
    b.write(out/'RESULT.json',r); print(json.dumps({'status':r['status'],'cell_id':cell_id})); return 0 if r['status']=='complete' else 1

def main():
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare'); a.add_argument('--config',required=True); a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('validate'); a.add_argument('--panel',required=True); a.add_argument('--out',required=True)
    a=sub.add_parser('run'); a.add_argument('--panel',required=True); a.add_argument('--release',required=True); a.add_argument('--cell',required=True); a.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.action=='prepare': prepare(a.config,a.out)
    elif a.action=='validate': b.write(a.out,validate(a.panel))
    else: return run(a.panel,a.release,a.cell,a.out)

if __name__=='__main__': sys.exit(main())
