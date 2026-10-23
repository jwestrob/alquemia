#!/usr/bin/env python3
"""Five finite PBE0 counterpoise/dispersion cells; native work requires Slurm."""
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
spec=importlib.util.spec_from_file_location('canonical_cp',HERE.parent/'canonical_water_cp_v2/REFERENCE.py')
cp=importlib.util.module_from_spec(spec);spec.loader.exec_module(cp)
b=cp.base; ROOT=b.ROOT; F=cp.FLOAT
WS=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1'
PARENT=WS/'canonical_water_cp_v2/prepared_v1/PANEL.json'
PARENT_SHA='00f7ab059686b5b8c475c99d7993f75d70bea8075c1b7baf9aa4f8cdc1c894dd'
PROTOCOL='lady_dft_water_cp_v1'

def electronic_parts(text):
    native=text.split('****END OF INPUT****')[-1]
    es=[cp.number(x) for x in re.findall(r'^Total Energy\s*:\s*('+F+r')\s+Eh',native,re.M)]
    ed=[cp.number(x) for x in re.findall(r'^Dispersion correction\s+('+F+r')',native,re.M)]
    ef=[cp.number(x) for x in re.findall(r'FINAL SINGLE POINT ENERGY\s+('+F+r')',native)]
    return es,ed,ef

def parse(text,c,returncode):
    native=text.split('****END OF INPUT****')[-1]; p=cp.parse(text,c,returncode); n=len(c['symbols']); es,ed,ef=electronic_parts(text)
    checks={key:p['checks'][key] for key in ('process_success','normal_termination','scf_converged','charge','multiplicity','total_electrons','basis_function_count','ecp_assignment')}
    rows=p['atomic_table']; coords=len(rows)==n and [x['index'] for x in rows]==list(range(n)) and [x['symbol'] for x in rows]==c['symbols']
    err=max((abs(a-z) for r,xx in zip(rows,c['positions_A']) for a,z in zip(r['positions_A'],xx)),default=math.inf) if coords else None
    checks.update(per_center_basis=p['basis_functions_by_center'][:n]==c['expected_basis_functions_by_center'],coordinates=coords and err<1e-6,
        ghost_nuclear_charges=coords and [r['effective_Z'] for r in rows]==c['expected_effective_nuclear_charges'],
        restricted_reference=bool(re.search(r'Hartree-Fock type\s+HFTyp\s*\.{2,}\s*RHF',native)),
        PBE_exchange=bool(re.search(r'Exchange Functional\s+Exchange\s*\.{2,}\s*PBE\b',native)),
        PBE_correlation=bool(re.search(r'Correlation Functional\s+Correlation\s*\.{2,}\s*PBE\b',native)),
        exact_exchange_25_percent=bool(re.search(r'Fraction HF Exchange\s+ScalHFX\s*\.{2,}\s*0\.250000\b',native)),
        RIJCOSX='RIJ-COSX (HFX calculated with COS-X)).... on' in native,
        dispersion_scope=(len(ed)==1 and 'DFTD4' in native) if c['d4'] else (not ed and 'DFT DISPERSION CORRECTION' not in native and 'DFTD4' not in native),
        electronic_final_closure=bool(es) and bool(ef) and all(math.isfinite(x) for x in es+ed+ef) and abs(es[-1]+(ed[-1] if ed else 0)-ef[-1])<1e-8)
    return {'status':'complete' if all(checks.values()) else 'failed_admission','checks':checks,'energy_hartree':ef[-1] if ef else None,
        'electronic_energy_hartree':es[-1] if es else None,'D4_energy_hartree':ed[-1] if ed else 0.,'printed_scf_criteria':p['printed_scf_criteria'],
        'atomic_table':rows,'basis_functions_by_center':p['basis_functions_by_center'][:n],'ecp_assigned_centers':p['ecp_assigned_centers'],
        'ecp_core_counts':p['ecp_core_counts'],'coordinate_max_error_A':err,'physical_state_qualified':False,'method_accuracy_qualified':False}

def deck(p,c,ranks,maxcore):
    basis=''
    if c['metal']:
        metal=c['metal']; full=b.verify(p['assets'][metal]['basis']).read_text(); ecp,gto=cp.metal_blocks(full,metal)
        basis=(full if c['role']=='ion_ghost_water' else 'DelECP '+metal+'\n'+gto+'\n')+b.verify(p['assets'][metal]['aux']).read_text()
    coords=[s+(':' if ghost else '')+' '+' '.join('%.12f'%z for z in x) for s,x,ghost in zip(c['symbols'],c['positions_A'],c['ghost_mask'])]
    return ('! PBE0 '+('D4 ' if c['d4'] else '')+'def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF\n'
        +f'%maxcore {maxcore}\n%pal nprocs {ranks} end\n%scf\n Guess PModel\n MaxIter 300\nend\n'
        +('%basis\n'+basis+'\nend\n' if basis else '')+f'* xyz {c["charge"]} 1\n'+'\n'.join(coords)+'\n*\n')

def prepare(out):
    tick=time.monotonic(); b.check(b.pin(PARENT)['sha256']==PARENT_SHA,'CP parent changed'); old=cp.load(PARENT); original=b.load_panel(cp.PARENT)
    out.mkdir(parents=True,exist_ok=False); (out/'assets').mkdir(); (out/'inputs').mkdir(); assets={}
    for metal in ('La','Dy'):
        assets[metal]={}
        for key in ('basis','aux'):
            dest=out/'assets'/f'{metal}_{key}.inc'; shutil.copyfile(b.verify(original['assets'][metal][key]),dest); assets[metal][key]=b.pin(dest)
    cells=[]
    for source in old['cells']:
        if source['role']=='complex':continue
        c=copy.deepcopy(source); c.pop('template_input'); c['d4']=False; cells.append(c)
    la=next(c for c in cells if c['cell_id']=='La__water_ghost_ion'); dy=next(c for c in cells if c['cell_id']=='Dy__water_ghost_ion')
    b.check(la['positions_A'][1:]==dy['positions_A'][1:] and la['symbols'][1:]==dy['symbols'][1:],'Free water differs between metals')
    cells.append({'cell_id':'real_water_D4','role':'real_water','metal':None,'symbols':la['symbols'][1:],'positions_A':la['positions_A'][1:],
        'source_geometry':la['source_geometry'],'ghost_mask':[False]*3,'charge':0,'multiplicity':1,'physical_multiplicity_of_parent':1,
        'expected_effective_nuclear_charges':[8,1,1],'expected_electrons':10,'expected_correlated_electrons':10,'expected_basis_functions':43,
        'expected_basis_functions_by_center':[31,6,6],'expected_ecp_core':0,'d4':True})
    saved={}
    for metal in ('La','Dy'):
        path=WS/'pilot_v1/tasks'/f'ecp_pbe0__water_central_{metal}'/'execution/RESULT.json'; r=json.loads(path.read_text()); b.check(r['status']=='complete' and all(r['checks'].values()),'Saved complex unavailable')
        output=b.verify(r['output']); es,ed,ef=electronic_parts(output.read_text()); b.check(len(es)==len(ed)==len(ef)==1 and abs(es[0]+ed[0]-ef[0])<1e-8,'Saved electronic/dispersion closure')
        g=next(x for x in old['cells'] if x['cell_id']==metal+'__complex'); b.check(r['geometry']==g['source_geometry'],'Saved geometry pin mismatch')
        saved[metal]={'result':b.pin(path),'output':b.pin(output),'source_geometry':g['source_geometry'],'electronic_energy_hartree':es[0],'D4_energy_hartree':ed[0],'final_energy_hartree':ef[0]}
    panel={'protocol':PROTOCOL,'adapter':b.pin(__file__),'canonical_helper':old['adapter'],'base_helper':old['helper'],'specification':b.pin(HERE/'SPEC.md'),
        'canonical_panel':b.pin(PARENT),'parent_panel':b.pin(cp.PARENT),'orca':old['orca'],'assets':assets,'cells':cells,'saved_complexes':saved,
        'isolated_ion_D4_hartree':{'La':0.,'Dy':0.},'ion_D4_evidence':'Analytical empty pair/triple sum, not a native calculation',
        'maximum_native_calls':5,'ranks':4,'created_utc':b.now()}
    for c in cells:
        dest=out/'inputs'/(c['cell_id']+'.inp'); dest.write_text(deck(panel,c,4,3000)); c['template_input']=b.pin(dest)
    p=b.write(out/'PANEL.json',panel); val=validate(p['path']); val.update(command=sys.argv,elapsed_seconds=time.monotonic()-tick)
    b.write(out/'PREFLIGHT.json',val); print(json.dumps({'panel':p,'status':val['status'],'native_calls':0}))

def load(path):
    p=json.loads(Path(path).read_text()); b.check(p['protocol']==PROTOCOL,'Wrong protocol')
    for key in ('adapter','canonical_helper','base_helper','specification','canonical_panel','parent_panel','orca'):b.verify(p[key])
    b.check(p['adapter']==b.pin(__file__),'Changed executing source'); b.check(len(p['cells'])==p['maximum_native_calls']==5,'Denominator differs')
    b.check({c['cell_id'] for c in p['cells']}=={m+'__'+r for m in ('La','Dy') for r in ('ion_ghost_water','water_ghost_ion')}|{'real_water_D4'},'Inventory differs')
    for aset in p['assets'].values():
        for a in aset.values():b.verify(a)
    for c in p['cells']:b.verify(c['source_geometry']);b.verify(c['template_input'])
    for r in p['saved_complexes'].values():
        for key in ('result','output','source_geometry'):b.verify(r[key])
    return p

def validate(path):
    p=load(path)
    for c in p['cells']:
        text=deck(p,c,4,3000); b.check(text==b.verify(c['template_input']).read_text(),'Input differs'); n=len(c['symbols'])
        b.check(c['expected_electrons']==sum(c['expected_effective_nuclear_charges'])-c['charge'],'Electrons differ')
        b.check(c['expected_basis_functions']==sum(c['expected_basis_functions_by_center']),'Basis inventory differs')
        b.check(('D4' in text)==c['d4'] and not any(x in text for x in ('EnGrad','CCSD','Opt ')),'Method changed')
        if c['role']=='water_ghost_ion':b.check('NewECP' not in text and 'DelECP '+c['metal'] in text,'Ghost ECP')
        if c['role']=='ion_ghost_water':b.check('N_core '+str(c['expected_ecp_core']) in text,'Real core missing')
        if c['role']=='real_water':b.check(n==3 and not any(c['ghost_mask']) and 'NewECP' not in text,'Real water invalid')
        parent_symbols,parent_x=b.xyz_read(c['source_geometry']['path']); start=1 if c['role']=='real_water' else 0
        b.check(c['symbols']==parent_symbols[start:] and c['positions_A']==parent_x[start:],'Geometry differs')
    return {'status':'passed_static_preflight','panel':b.pin(path),'cells':len(p['cells']),'molecular_calls':0,'local_threads':1,'utc':b.now()}

def run(path,cell_id,out):
    p=load(path); b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required'); c=next(x for x in p['cells'] if x['cell_id']==cell_id)
    ranks=int(os.environ.get('SLURM_NTASKS','0'));b.check(ranks==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread ranks required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));b.check(mem>0,'Actual memory required'); maxcore=math.floor(.75*mem/ranks)
    out.mkdir(parents=True,exist_ok=False); (out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore));env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[key]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH',''); cmd=[p['orca']['path'],'endpoint.inp','--bind-to none']
    b.write(out/'START.json',{'panel':b.pin(path),'cell':c,'input':b.pin(out/'endpoint.inp'),'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'ranks':ranks,'memory_MiB':mem,'maxcore_MiB':maxcore,'command':cmd,'utc':b.now()}); tick=time.monotonic()
    with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:proc=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False)
    elapsed=time.monotonic()-tick
    try:r=parse((out/'endpoint.out').read_text(errors='replace'),c,proc.returncode)
    except Exception as e:r={'status':'failed_collection','exception':repr(e)}
    r.update(cell_id=cell_id,cell=c,panel=b.pin(path),adapter=b.pin(__file__),returncode=proc.returncode,elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'),input=b.pin(out/'endpoint.inp'),finished_utc=b.now())
    b.write(out/'RESULT.json',r); print(json.dumps({'cell':cell_id,'status':r['status']}));return 0 if r['status']=='complete' else 1

def collect(path,paths,canonical):
    p=load(path); rows={}
    for filename in paths:
        r=json.loads(Path(filename).read_text());b.check(r['status']=='complete' and all(r['checks'].values()),'Incomplete native result');b.check(r['panel']==b.pin(path),'Panel differs');b.check(r['cell_id'] not in rows,'Duplicate cell');b.verify(r['output']);rows[r['cell_id']]=r
    b.check(set(rows)=={c['cell_id'] for c in p['cells']},'Five-cell denominator required')
    corr=json.loads(Path(canonical).read_text());b.check(corr['status']=='complete' and corr['panel']==p['canonical_panel'],'Canonical collection differs')
    for x in corr['result_pins']:b.verify(x)
    metals={}; water=rows['real_water_D4']['D4_energy_hartree']
    for m in ('La','Dy'):
        e=(p['saved_complexes'][m]['electronic_energy_hartree']-rows[m+'__ion_ghost_water']['electronic_energy_hartree']-rows[m+'__water_ghost_ion']['electronic_energy_hartree'])*b.HA
        disp=(p['saved_complexes'][m]['D4_energy_hartree']-p['isolated_ion_D4_hartree'][m]-water)*b.HA
        metals[m]={'PBE0_CP_kcal_mol':e,'D4_interaction_kcal_mol':disp,'PBE0_D4_CP_kcal_mol':e+disp,'CCSD_T_CP_kcal_mol':corr['CP_interaction_kcal_mol'][m]}
    contrast={key:metals['Dy'][key]-metals['La'][key] for key in metals['La']}
    return {'status':'complete','panel':b.pin(path),'new_result_pins':[b.pin(x) for x in paths],'canonical_collection':b.pin(canonical),'per_metal':metals,'Dy_minus_La':contrast,'molecular_calls':0,'physical_state_qualified':False,'method_accuracy_qualified':False}

def main():
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('validate');a.add_argument('--panel',required=True);a.add_argument('--out',required=True)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--cell',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--canonical',required=True);a.add_argument('--out',required=True);a.add_argument('results',nargs='+')
    a=ap.parse_args()
    if a.action=='prepare':prepare(a.out)
    elif a.action=='validate':b.write(a.out,validate(a.panel))
    elif a.action=='run':return run(a.panel,a.cell,a.out)
    else:b.write(a.out,collect(a.panel,a.results,a.canonical))

if __name__=='__main__':sys.exit(main())
