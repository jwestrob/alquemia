#!/usr/bin/env python3
"""68 literal fragments in each consumed parent's complete orbital/auxiliary basis."""
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

HERE=Path(__file__).resolve().parent

def module(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
helper=module(HERE.parent/'coordination_cp_v1/REFERENCE.py','common_basis_cp_helper')
prior=module(HERE.parent/'coordination_cp_v1/collection/COLLECT.py','common_basis_saved_cp')
b=helper.b;read=helper.read;GIDS=helper.GIDS
PROTOCOL='common_parent_basis_decomposition_v1'
ROLES=('metal_water','water_only')
SCOUTS=[f'pbe0__aquo9__Dy__water_01__{r}' for r in ROLES]
HELPERS={**helper.HELPERS,'coordination_helper':HERE.parent/'coordination_cp_v1/REFERENCE.py',
    'coordination_collector':HERE.parent/'coordination_cp_v1/collection/COLLECT.py',
    'matched_pair_helper':HERE.parent/'matched_pair_decomposition_v1/REFERENCE.py'}
AUX={f'aquo{n}_central_{m}':v for n,m,v in [(8,'La',919),(8,'Dy',937),(9,'La',990),(9,'Dy',1008)]}

def definition(g,i,role):
    n=g['water_count'];m=g['metal'];b.check(role in ROLES and 1<=i<=n,'Invalid fragment')
    water=[3*i-2,3*i-1,3*i];real=([0]+water) if role=='metal_water' else water
    ghost=[j not in real for j in range(len(g['symbols']))];charge=3 if role=='metal_water' else 0
    z=[0 if ghost[j] else (11 if j==0 else b.Z[s]) for j,s in enumerate(g['symbols'])]
    return dict(cell_id=f'pbe0__aquo{n}__{m}__water_{i:02d}__{role}',geometry_id=g['geometry_id'],parent_geometry_id=g['geometry_id'],
        water_count=n,water_index=i,real_parent_atom_indices=real,water_parent_atom_indices=water,method='pbe0',metal=m,role=role,d4=False,
        symbols=g['symbols'],positions_A=g['positions_A'],source_geometry=g['xyz'],ghost_mask=ghost,charge=charge,multiplicity=1,
        physical_multiplicity_of_parent=g['physical_multiplicity'],expected_effective_nuclear_charges=z,
        expected_electrons=sum(z)-charge,expected_correlated_electrons=sum(z)-charge,expected_basis_functions=47+43*n,
        expected_basis_functions_by_center=[47]+[31,6,6]*n,expected_aux_J_basis_functions=AUX[g['geometry_id']],
        expected_ecp_core=(46 if m=='La' else 55) if role=='metal_water' else 0)

def deck(p,c,ranks,maxcore):return helper.deck(p,c,ranks,maxcore)

def auxiliary_counts(text):
    native=text.split('****END OF INPUT****')[-1]
    return {'basis_header':[int(v) for v in re.findall(r'# of basis functions in Aux-J\s*\.{3,}\s*(\d+)',native)],
        'integral_header':[int(v) for v in re.findall(r'Number of AuxJ basis functions\s*\.{3,}\s*(\d+)',native)]}

def auxiliary_admission(text,expected):
    counts=auxiliary_counts(text)
    return counts,all(values and all(v==expected for v in values) for values in counts.values())

def parse(text,c,returncode):
    r=helper.parse(text,c,returncode);counts,passed=auxiliary_admission(text,c['expected_aux_J_basis_functions'])
    r['native_aux_J_basis_function_counts']=counts;r['checks']['full_parent_aux_J_basis']=passed
    r['status']='complete' if all(r['checks'].values()) else 'failed_admission'
    return r

def reused_admission(cp_panel,cp_pin,cp_result,data):
    b.check(cp_result['status']=='coordination_cp_collection_complete','Prior CP cycle unavailable')
    rows={}
    for gid in GIDS:
        q=prior.admit_original(gid,cp_panel,data,helper)
        counts,passed=auxiliary_admission(b.verify(q['output']).read_text(),AUX[gid]);b.check(passed,'Original auxiliary basis differs: '+gid)
        rows[gid]=dict(q,expected_aux_J_basis_functions=AUX[gid],native_aux_J_basis_function_counts=counts)
    for c in cp_panel['cells']:
        src=cp_result['fragments'][c['cell_id']]['result'];q=prior.admit_fragment(b.verify(src),c,cp_panel,cp_pin,helper)
        counts,passed=auxiliary_admission(b.verify(q['output']).read_text(),AUX[c['geometry_id']]);b.check(passed,'Reused fragment auxiliary basis differs')
        rows[c['cell_id']]=dict(q,expected_aux_J_basis_functions=AUX[c['geometry_id']],native_aux_J_basis_function_counts=counts)
    b.check(len(rows)==12,'Reused denominator differs');return rows

def prepare(config,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config).read_text())
    b.check(cfg['implementation']==b.pin(__file__),'Preparation source differs')
    for key,path in HELPERS.items():b.check(cfg[key]==b.pin(path),'Helper differs: '+key)
    for key in ('declaration','proposal','cp_final_admission','pair_final_admission'):b.verify(cfg[key])
    old=helper.load(b.verify(cfg['cp_panel']));pairs=module(HELPERS['matched_pair_helper'],'common_basis_pair_helper').load(b.verify(cfg['pair_panel']))
    data=read(old['consumed_input_data']);cp_result=read(cfg['cp_collection']);pair_result=read(cfg['pair_collection'])
    b.check(pair_result['status']=='matched_pair_collection_complete' and pair_result['reference_complete']==36 and pair_result['model_complete']==34,'Prior literal pair result unavailable')
    b.check(pairs['parent_cp_panel']==cfg['cp_panel'] and pairs['parent_geometries']==old['geometries'],'Prior parent geometries differ')
    reused=reused_admission(old,cfg['cp_panel'],cp_result,data)
    out.mkdir(parents=True,exist_ok=False)
    for name in ('assets','inputs'):(out/name).mkdir()
    assets={}
    for m in ('La','Dy'):
        assets[m]={}
        for key in ('basis','aux'):
            dst=out/'assets'/f'{m}_{key}.inc';shutil.copyfile(b.verify(old['assets'][m][key]),dst);assets[m][key]=b.pin(dst)
    p=dict(protocol=PROTOCOL,adapter=b.pin(__file__),configuration=b.pin(config),declaration=cfg['declaration'],proposal=cfg['proposal'],
        cp_panel=cfg['cp_panel'],pair_panel=cfg['pair_panel'],cp_collection=cfg['cp_collection'],pair_collection=cfg['pair_collection'],
        cp_final_admission=cfg['cp_final_admission'],pair_final_admission=cfg['pair_final_admission'],consumed_input_data=old['consumed_input_data'],
        parent_reference_panel=old['parent_panel'],geometries=old['geometries'],orca=old['orca'],assets=assets,source_assets=old['assets'],
        saved_reuse=reused,cells=[],maximum_native_calls=68,reused_cells=12,ranks=4,scout_cell_ids=SCOUTS,run_wrapper=b.pin(HERE/'RUN.sbatch'),
        created_utc=b.now(),interpretation='Common-parent-basis electronic C/P/B; N-fold ion subtraction; saved D4 added once separately; original targets unchanged',**{k:cfg[k] for k in HELPERS})
    for gid in GIDS:
        g=p['geometries'][gid]
        for i in range(1,g['water_count']+1):
            own=next(c for c in pairs['cells'] if c['cell_id']==f'pbe0__aquo{g["water_count"]}__{g["metal"]}__water_{i:02d}')
            indices=[0,3*i-2,3*i-1,3*i]
            b.check(own['parent_atom_indices']==indices and own['positions_A']==[g['positions_A'][j] for j in indices],'Literal extracted dimer differs')
            p['cells'].extend(definition(g,i,r) for r in ROLES)
    for c in p['cells']:
        f=out/'inputs'/(c['cell_id']+'.inp');f.write_text(deck(p,c,4,3000));c['template_input']=b.pin(f)
    pin=b.write(out/'PANEL.json',p);v=validate(pin['path']);v.update(reused_raw_admissions=12,reused_auxiliary_counts={k:v['native_aux_J_basis_function_counts'] for k,v in reused.items()},
        command=sys.argv,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu)
    b.write(out/'PREFLIGHT.json',v);print(json.dumps({'panel':pin,'status':v['status'],'new_native_calls':0}))

def load(path):
    p=json.loads(Path(path).read_text());b.check(p['protocol']==PROTOCOL and p['adapter']==b.pin(__file__),'Protocol/source differs')
    for key in ('adapter','configuration','declaration','proposal','cp_panel','pair_panel','cp_collection','pair_collection','cp_final_admission','pair_final_admission','consumed_input_data','parent_reference_panel','orca','run_wrapper',*HELPERS):b.verify(p[key])
    for key,src in HELPERS.items():b.check(p[key]==b.pin(src),'Executing helper differs')
    b.check(len(p['cells'])==p['maximum_native_calls']==68 and p['ranks']==4 and p['reused_cells']==len(p['saved_reuse'])==12 and p['scout_cell_ids']==SCOUTS and set(p['geometries'])==set(GIDS),'Finite inventory differs')
    expected=[definition(p['geometries'][gid],i,r) for gid in GIDS for i in range(1,p['geometries'][gid]['water_count']+1) for r in ROLES]
    b.check([c['cell_id'] for c in p['cells']]==[c['cell_id'] for c in expected],'Ordered literal inventory differs')
    for m,a in p['assets'].items():
        for key,pin in a.items():b.verify(pin);b.verify(p['source_assets'][m][key]);b.check(pin['sha256']==p['source_assets'][m][key]['sha256'],'Copied asset differs')
    for row in p['saved_reuse'].values():
        for key in ('result','output'):b.verify(row[key])
    for c,e in zip(p['cells'],expected):
        b.check({k:v for k,v in c.items() if k!='template_input'}==e,'Fragment definition differs');b.verify(c['source_geometry'])
        b.check(b.verify(c['template_input']).read_text()==deck(p,c,4,3000),'Frozen template differs')
    return p

def validate(path):
    p=load(path);rows=[]
    for c in p['cells']:
        text=b.verify(c['template_input']).read_text();sy,x=b.xyz_read(c['source_geometry']['path'])
        b.check(sy==c['symbols'] and x==c['positions_A'],'Exact parent geometry differs')
        head,*lines=text.split('* xyz ',1)[1].splitlines();lines=lines[:-1];q,m=map(int,head.split());tokens=[v.split() for v in lines]
        labels=[t[0] for t in tokens];xyz=[[float(v) for v in t[1:]] for t in tokens]
        z=[0 if s.endswith(':') else (11 if j==0 else b.Z[s]) for j,s in enumerate(labels)]
        b.check(q==c['charge'] and m==1 and z==c['expected_effective_nuclear_charges'] and sum(z)-q==c['expected_electrons']==(18 if c['role']=='metal_water' else 10),'Charge/electron/ghost inventory differs')
        b.check([s.endswith(':') for s in labels]==c['ghost_mask'] and [s.rstrip(':') for s in labels]==sy and xyz==x,'Full centers or literal coordinates differ')
        b.check(sum(c['expected_basis_functions_by_center'])==c['expected_basis_functions']==47+43*c['water_count'],'Orbital inventory differs')
        b.check('D4' not in text and 'EnGrad' not in text and 'MaxIter 300' in text and 'DefGrid3' in text and 'VeryTightSCF' in text and 'Guess PModel' in text,'Electronic SP recipe differs')
        b.check(('NewECP' in text)==bool(c['expected_ecp_core']) and ('DelECP '+c['metal'] in text)==c['ghost_mask'][0],'Ghost ECP differs')
        b.check(helper.cp.metal_blocks(b.verify(p['assets'][c['metal']]['basis']).read_text(),c['metal'])[1] in text and b.verify(p['assets'][c['metal']]['aux']).read_text() in text,'Orbital/auxiliary metal basis differs')
        rows.append({k:c[k] for k in ('cell_id','real_parent_atom_indices','charge','expected_electrons','expected_basis_functions','expected_aux_J_basis_functions','expected_ecp_core')})
    return dict(status='passed_static_preflight',panel=b.pin(path),declared=68,reused=12,independent_deck_inventory=rows,new_native_calls=0,local_threads=1,native_admission_pending=True,utc=b.now())

def release_contract(p,path,r):
    b.check(r['status']=='approved_common_basis_decomposition_cells' and r['implementation']==p['adapter'] and r['panel']==b.pin(path),'Root release differs')
    ids=[c['cell_id'] for c in p['cells']]
    b.check(r['cell_ids']==ids and r['maximum_native_calls']==68 and r['ranks']==4,'Released full scope differs')
    pre=read(r['preparation']);b.verify(r['independent_review']);b.check(pre['status']=='passed_static_preflight' and pre['panel']==b.pin(path),'Preparation not admitted')
    expected=SCOUTS if r['stage']=='scouts' else [v for v in ids if v not in SCOUTS]
    b.check(r['stage'] in ('scouts','remaining') and r['released_cell_ids']==expected,'Released subset differs')
    if r['stage']=='remaining':
        review=read(r['scout_review']);b.check(review['status']=='passed_scout_review' and review['panel']==b.pin(path) and review['independent_review']==r['independent_review'],'Actual scout review missing or different')
        b.check(set(r['scout_results'])==set(SCOUTS) and review['scout_results']==r['scout_results'],'Scout authority denominator differs')
        for cid,pin in r['scout_results'].items():
            c=next(c for c in p['cells'] if c['cell_id']==cid);admit_result(b.verify(pin),c,p,path,require_scout=True)
    return r

def gate(p,path,release,cell_id,out):
    r=release_contract(p,path,json.loads(Path(release).read_text()));b.check(cell_id in r['released_cell_ids'],'Cell not released')
    b.check(out.resolve()==Path(r['execution_root']).resolve()/'tasks'/cell_id/'execution','One-attempt path differs');b.check(not out.exists(),'Existing attempt must be preserved');return r

def admit_result(path,c,p,panel_path,require_scout=False):
    path=Path(path);r=json.loads(path.read_text());b.check(r['status']=='complete' and r['returncode']==0 and all(r['checks'].values()),'Native admission failed')
    b.check(r['cell_id']==c['cell_id'] and r['cell']==c and r['panel']==b.pin(panel_path) and r['adapter']==p['adapter'],'Native source/cell/panel differs')
    start=read(r['start']);b.check(r['start']==b.pin(path.parent/'START.json') and start['cell']==c and start['panel']==r['panel'] and start['input']==r['input'] and start['release']==r['release'],'Native START lineage differs')
    release=read(r['release']);b.check(not require_scout or release['stage']=='scouts','Scout must come from original scout release');release_contract(p,panel_path,release)
    b.check(c['cell_id'] in release['released_cell_ids'] and path.resolve()==Path(release['execution_root']).resolve()/'tasks'/c['cell_id']/'execution/RESULT.json','Native attempt path/subset differs')
    b.check(start['ranks']==4 and start['memory_MiB']>0 and start['maxcore_MiB']==math.floor(.75*start['memory_MiB']/4),'Actual resources differ')
    b.check(start['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'] and start['maximum_native_calls']==r['native_program_invocations_attempted']==1,'Native invocation differs')
    b.verify(r['stderr']);b.check(b.verify(r['input']).read_text()==deck(p,c,4,start['maxcore_MiB']),'Executed deck differs')
    q=parse(b.verify(r['output']).read_text(),c,r['returncode']);b.check(all(r[k]==v for k,v in q.items()) and q['status']=='complete','Raw native replay differs')
    b.check(q['D4_energy_hartree']==0. and math.isfinite(q['electronic_energy_hartree']),'Electronic fragment invalid')
    return dict(q,result=b.pin(path),start=r['start'],input=r['input'],output=r['output'],release=r['release'],job_id=r['job_id'],allocated_CPU_seconds=r.get('allocated_cpu_seconds'))

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
