#!/usr/bin/env python3
"""Literal 34 extracted dimers plus two own-basis ions; root releases each SP."""
import argparse
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

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('matched_pair_cp_helper',HERE.parent/'coordination_cp_v1/REFERENCE.py')
helper=importlib.util.module_from_spec(sp);sp.loader.exec_module(helper)
b=helper.b;GIDS=helper.GIDS;PROTOCOL='literal_matched_pair_decomposition_v1'
HELPERS={**helper.HELPERS,'coordination_helper':HERE.parent/'coordination_cp_v1/REFERENCE.py'}
read=helper.read

def geometry(parent, water_index):
    m=parent['metal'];n=parent['water_count'];ion=water_index is None
    indices=[0] if ion else [0,3*water_index-2,3*water_index-1,3*water_index]
    b.check(ion or 1<=water_index<=n,'Water index outside parent')
    gid=f'{m}_bare_ion' if ion else f'aquo{n}_{m}_water_{water_index:02d}'
    return dict(geometry_id=gid,symbols=[parent['symbols'][i] for i in indices],positions_A=[parent['positions_A'][i] for i in indices],
        atoms=len(indices),metal=m,metal_index=0,charge=3,physical_multiplicity=parent['physical_multiplicity'],
        parent_geometry_id=parent['geometry_id'],water_index=water_index,parent_atom_indices=indices,
        kind='literal_own_basis_bare_ion' if ion else 'literal_extracted_metal_water')

def definition(g):
    m=g['metal'];ion=g['atoms']==1;role='bare_ion' if ion else 'dimer';n=g['parent_geometry_id'].split('_')[0]
    cid=f'pbe0__{m}__bare_ion' if ion else f'pbe0__{n}__{m}__water_{g["water_index"]:02d}'
    z=[11] if ion else [11,8,1,1];nbf=[47] if ion else [47,31,6,6]
    return dict(cell_id=cid,geometry_id=g['geometry_id'],geometry=g,parent_geometry_id=g['parent_geometry_id'],water_index=g['water_index'],
        parent_atom_indices=g['parent_atom_indices'],method='pbe0',metal=m,role=role,d4=True,symbols=g['symbols'],positions_A=g['positions_A'],
        source_geometry=g['xyz'],ghost_mask=[False]*g['atoms'],charge=3,multiplicity=1,physical_multiplicity_of_parent=g['physical_multiplicity'],
        expected_effective_nuclear_charges=z,expected_electrons=sum(z)-3,expected_correlated_electrons=sum(z)-3,
        expected_basis_functions=sum(nbf),expected_basis_functions_by_center=nbf,expected_ecp_core=46 if m=='La' else 55)

def deck(p,c,ranks,maxcore):return helper.deck(p,c,ranks,maxcore)

def parse(text,c,returncode):
    r=helper.parse(text,c,returncode)
    if c['role']=='bare_ion':
        r['checks']['isolated_atom_D4_zero']=r['D4_energy_hartree']==0.
        r['status']='complete' if all(r['checks'].values()) else 'failed_admission'
    return r

def prepare(config,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config).read_text())
    b.check(cfg['implementation']==b.pin(__file__),'Executing source differs')
    for k,f in HELPERS.items():b.check(cfg[k]==b.pin(f),'Helper changed: '+k)
    root=read(cfg['root_inputs']);b.verify(root['declaration']);b.verify(root['cp_final_admission']);b.verify(root['cp_collection'])
    b.check(root['source_cells']==GIDS and root['new_reference_cells']==36 and root['new_model_analyze_calls']==34 and root['new_native_calls_authorized_by_this_file'] is False,'Root finite scope differs')
    old=helper.load(b.verify(root['cp_panel']));data=read(root['consumed_input_data'])
    b.check(old['consumed_input_data']==root['consumed_input_data'],'Original consumed DATA differs')
    out.mkdir(parents=True,exist_ok=False)
    for name in ('assets','xyz','inputs'):(out/name).mkdir()
    assets={}
    for m in ('La','Dy'):
        assets[m]={}
        for key in ('basis','aux'):
            target=out/'assets'/f'{m}_{key}.inc';shutil.copyfile(b.verify(old['assets'][m][key]),target);assets[m][key]=b.pin(target)
    p=dict(protocol=PROTOCOL,adapter=b.pin(__file__),configuration=b.pin(config),root_inputs=cfg['root_inputs'],declaration=root['declaration'],
        parent_cp_panel=root['cp_panel'],consumed_input_data=root['consumed_input_data'],cp_collection_context_only=root['cp_collection'],
        parent_reference_panel=old['parent_panel'],orca=old['orca'],assets=assets,source_assets=old['assets'],parent_geometries=old['geometries'],
        saved_complexes=old['saved_complexes'],baseline_model_energies=old['baseline_model_energies'],cells=[],maximum_native_calls=36,ranks=4,
        run_wrapper=b.pin(HERE/'RUN.sbatch'),cancellation_checks=old['cancellation_checks'],created_utc=b.now(),
        interpretation='Literal own-basis paired dimer sum and cooperative remainder; no ghost corrections, forces, fitting or physical truth claim',**{k:cfg[k] for k in HELPERS})
    for gid in GIDS:
        parent=p['parent_geometries'][gid];row=data['rows'][gid]
        b.check(parent==row['geometry'] and p['saved_complexes'][gid]['result']==row['reference_result'],'Consumed parent join differs')
        b.check(p['baseline_model_energies'][gid]==data['baseline_model'][gid]['energy'],'Frozen model energy differs')
        for i in range(1,parent['water_count']+1):
            g=geometry(parent,i);g['xyz']=b.xyz_write(out/'xyz'/(g['geometry_id']+'.xyz'),g['symbols'],g['positions_A'],g['geometry_id']);p['cells'].append(definition(g))
    for m in ('La','Dy'):
        parent=p['parent_geometries']['aquo8_central_'+m];other=p['parent_geometries']['aquo9_central_'+m]
        b.check(parent['positions_A'][0]==other['positions_A'][0]==[0.,0.,0.] and parent['physical_multiplicity']==other['physical_multiplicity'],'Own-ion position/state differs')
        g=geometry(parent,None);g['xyz']=b.xyz_write(out/'xyz'/(g['geometry_id']+'.xyz'),g['symbols'],g['positions_A'],g['geometry_id']);p['cells'].append(definition(g))
    for c in p['cells']:
        f=out/'inputs'/(c['cell_id']+'.inp');f.write_text(deck(p,c,4,3000));c['template_input']=b.pin(f)
    pin=b.write(out/'PANEL.json',p);v=validate(pin['path']);v.update(command=sys.argv,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu)
    b.write(out/'PREFLIGHT.json',v);print(json.dumps({'panel':pin,'status':v['status'],'new_native_calls':0}))

def load(path):
    p=json.loads(Path(path).read_text());b.check(p['protocol']==PROTOCOL and p['adapter']==b.pin(__file__),'Protocol/source differs')
    for key in ('adapter','configuration','root_inputs','declaration','parent_cp_panel','consumed_input_data','cp_collection_context_only','parent_reference_panel','orca','run_wrapper',*HELPERS):b.verify(p[key])
    for key,f in HELPERS.items():b.check(p[key]==b.pin(f),'Executing helper changed')
    b.check(p['maximum_native_calls']==len(p['cells'])==36 and p['ranks']==4 and set(p['parent_geometries'])==set(GIDS),'Finite denominator differs')
    expected={f'pbe0__aquo{n}__{m}__water_{i:02d}' for n in (8,9) for m in ('La','Dy') for i in range(1,n+1)}|{f'pbe0__{m}__bare_ion' for m in ('La','Dy')}
    b.check({c['cell_id'] for c in p['cells']}==expected,'Literal inventory differs')
    for m,a in p['assets'].items():
        for key,pin in a.items():b.verify(pin);b.verify(p['source_assets'][m][key]);b.check(pin['sha256']==p['source_assets'][m][key]['sha256'],'Copied asset changed')
    for r in p['saved_complexes'].values():
        for k in ('result','output','input','engrad','source_geometry'):b.verify(r[k])
    for g in p['parent_geometries'].values():b.verify(g['xyz'])
    for c in p['cells']:
        g=c['geometry'];parent=p['parent_geometries'][c['parent_geometry_id']]
        b.check({k:v for k,v in g.items() if k!='xyz'}==geometry(parent,c['water_index']),'Extracted parent/atom mapping differs')
        b.check({k:v for k,v in c.items() if k!='template_input'}==definition(g),'Native cell definition differs')
        b.verify(c['source_geometry']);b.verify(c['template_input']);b.check(b.verify(c['template_input']).read_text()==deck(p,c,4,3000),'Exact template changed')
    return p

def validate(path):
    p=load(path);rows=[]
    for c in p['cells']:
        text=b.verify(c['template_input']).read_text();sy,x=b.xyz_read(c['source_geometry']['path']);g=c['geometry'];parent=p['parent_geometries'][c['parent_geometry_id']]
        b.check(sy==c['symbols'] and x==c['positions_A']==[parent['positions_A'][i] for i in c['parent_atom_indices']],'Literal subset coordinates differ')
        header,*lines=text.split('* xyz ',1)[1].splitlines();lines=lines[:-1];q,m=map(int,header.split());labels=[s.split()[0] for s in lines]
        ion=c['role']=='bare_ion';expected=[c['metal']] if ion else [c['metal'],'O','H','H']
        b.check(labels==expected and q==3 and m==1 and not any(c['ghost_mask']),'Own-basis charge/atom/ghost inventory differs')
        b.check(c['expected_electrons']==(8 if ion else 18) and c['expected_basis_functions']==(47 if ion else 90),'Electron/AO inventory differs')
        b.check('D4' in text and 'EnGrad' not in text and 'DelECP' not in text and 'NewECP '+c['metal'] in text,'Real-atom SP recipe differs')
        b.check('N_core '+str(c['expected_ecp_core']) in text and 'MaxIter 300' in text and 'VeryTightSCF' in text and 'DefGrid3' in text and 'Guess PModel' in text,'Frozen core/SCF/grid differs')
        b.check(helper.cp.metal_blocks(b.verify(p['assets'][c['metal']]['basis']).read_text(),c['metal'])[1] in text and b.verify(p['assets'][c['metal']]['aux']).read_text() in text,'Metal basis/aux differs')
        rows.append({k:c[k] for k in ('cell_id','parent_geometry_id','water_index','parent_atom_indices','expected_electrons','expected_basis_functions','expected_ecp_core')})
    return dict(status='passed_static_preflight',panel=b.pin(path),declared=36,dimers=34,own_ions=2,independent_deck_inventory=rows,new_native_calls=0,local_threads=1,native_admission_pending=True,utc=b.now())

def gate(p,path,release,cell_id,out):
    r=json.loads(Path(release).read_text())
    b.check(r['status']=='approved_matched_pair_36_cells' and r['implementation']==p['adapter'] and r['panel']==b.pin(path),'Root release differs')
    b.check(r['cell_ids']==[c['cell_id'] for c in p['cells']] and r['maximum_native_calls']==36 and r['ranks']==4,'Released scope differs')
    pre=read(r['preparation']);b.verify(r['independent_review']);b.check(pre['status']=='passed_static_preflight' and pre['panel']==b.pin(path),'Preparation not admitted')
    b.check(out.resolve()==Path(r['execution_root']).resolve()/'tasks'/cell_id/'execution','One-attempt path differs');b.check(not out.exists(),'Existing attempt must be preserved')

def run(path,release,cell_id,out):
    p=load(path);matches=[c for c in p['cells'] if c['cell_id']==cell_id];b.check(len(matches)==1,'Unknown cell');c=matches[0]
    gate(p,path,release,cell_id,out);b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native calls require Slurm')
    ranks=int(os.environ.get('SLURM_NTASKS','0'));b.check(ranks==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four one-thread ranks required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));b.check(mem>0,'Actual memory grant required')
    maxcore=math.floor(.75*mem/ranks);b.check(maxcore>0,'Invalid memory per rank');out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore))
    env=os.environ.copy()
    for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[k]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','');cmd=[p['orca']['path'],'endpoint.inp','--bind-to none']
    start=b.write(out/'START.json',dict(cell=c,panel=b.pin(path),release=b.pin(release),input=b.pin(out/'endpoint.inp'),job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),ranks=ranks,memory_MiB=mem,maxcore_MiB=maxcore,command=cmd,maximum_native_calls=1,utc=b.now()))
    tick=time.monotonic();code=None
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:
            proc=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False);code=proc.returncode
        r=parse((out/'endpoint.out').read_text(errors='replace'),c,code)
    except Exception as exc:r=dict(status='failed_execution_or_collection',exception=repr(exc))
    elapsed=time.monotonic()-tick;r.update(cell_id=cell_id,cell=c,panel=b.pin(path),adapter=p['adapter'],release=b.pin(release),start=start,returncode=code,elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],finished_utc=b.now(),input=b.pin(out/'endpoint.inp'),native_program_invocations_attempted=1,output=b.pin(out/'endpoint.out') if (out/'endpoint.out').exists() else None,stderr=b.pin(out/'endpoint.err') if (out/'endpoint.err').exists() else None)
    b.write(out/'RESULT.json',r);print(json.dumps({'status':r['status'],'cell_id':cell_id}));return 0 if r['status']=='complete' else 1

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--config',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('validate');a.add_argument('--panel',required=True);a.add_argument('--out',required=True)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--release',required=True);a.add_argument('--cell',required=True);a.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.action=='prepare':prepare(a.config,a.out)
    elif a.action=='validate':b.write(a.out,validate(a.panel))
    else:return run(a.panel,a.release,a.cell,a.out)

if __name__=='__main__':sys.exit(main())
