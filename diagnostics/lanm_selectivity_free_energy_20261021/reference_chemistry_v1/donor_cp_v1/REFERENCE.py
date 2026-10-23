#!/usr/bin/env python3
"""Finite matched-donor CP references; no native calls outside Slurm."""
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
sp=importlib.util.spec_from_file_location('dft_water',HERE.parent/'dft_water_cp_v1/REFERENCE.py');dft=importlib.util.module_from_spec(sp);sp.loader.exec_module(dft)
cp=dft.cp;b=dft.b;ROOT=b.ROOT;F=cp.FLOAT
PARENT=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/matched_water_removal_v1/prepared_v1/NATIVE_PANEL.json'
PARENT_SHA='304da68b3ee13ca2604cf498c107bd43f934f8ca849e101564e4f37c6ff8a821'
PARENT_TASKS=PARENT.parent.parent/'run_v1/reference/tasks'
WATER=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dft_water_cp_v1/final_collection_v1/RESULT.json'
WATER_SHA='6465f8c9002c31ee9c79a945152476c91ad740da0a5ddca4e2327d98321bda64'
PROTOCOL='matched_donor_cp_v1';DONORS=('acetamide','acetate_bi');METALS=('La','Dy');ROLES=('complex','ion_ghost_ligand','ligand_ghost_ion')
NBF={'La':47,'Dy':47,'C':31,'N':31,'O':31,'H':6}

def deck(p,c,ranks,maxcore):
    basis=''
    if c['metal']:
        metal=c['metal'];full=b.verify(p['assets'][metal]['basis']).read_text();ecp,gto=cp.metal_blocks(full,metal)
        basis=full if not c['ghost_mask'][0] else 'DelECP '+metal+'\n'+gto+'\n'
        if c['method']=='pbe0':basis+=b.verify(p['assets'][metal]['aux']).read_text()
    keyword='RHF CCSD(T) def2-TZVP NoFrozenCore' if c['method']=='cc' else 'PBE0 '+('D4 ' if c['d4'] else '')+'def2-TZVP def2/J RIJCOSX DefGrid3'
    text=f'! {keyword} VeryTightSCF NoAutostart\n%maxcore {maxcore}\n%pal nprocs {ranks} end\n%scf\n Guess PModel\n MaxIter 300\nend\n'
    if c['method']=='cc':text+='%mdci\n MaxIter 100\n STol 1e-8\n PrintLevel 3\nend\n'
    if basis:text+='%basis\n'+basis+'\nend\n'
    text+=f'* xyz {c["charge"]} 1\n'
    text+='\n'.join(s+(':' if ghost else '')+' '+' '.join('%.12f'%z for z in x) for s,x,ghost in zip(c['symbols'],c['positions_A'],c['ghost_mask']))+'\n*\n'
    return text

def parse(text,c,returncode):
    p=(cp.parse if c['method']=='cc' else dft.parse)(text,c,returncode);native=text.split('****END OF INPUT****')[-1];n=len(c['symbols'])
    section=native.split('BASIS SET INFORMATION',1)[-1].split('AUXILIARY',1)[0].split('ECP PARAMETER INFORMATION',1)[0].split('SHARK',1)[0]
    groups={int(i):cp.shell_count(contraction) for i,contraction in re.findall(r'Group\s+(\d+)\s+Type\s+[A-Za-z]+:?\s*:\s*[^\n]*?contracted to\s+([0-9spdfghi]+)',section,re.I)}
    assignment={int(i):int(g) for i,g in re.findall(r'Atom\s+(\d+)\s*[A-Za-z]+:?\s+basis set group\s*=>\s*(\d+)',section)}
    counts=[groups.get(assignment.get(i)) for i in range(n)];rows=p['atomic_table'];same=len(rows)==n and [x['index'] for x in rows]==list(range(n)) and [x['symbol'] for x in rows]==c['symbols']
    error=max((abs(x-y) for row,pos in zip(rows,c['positions_A']) for x,y in zip(row['positions_A'],pos)),default=math.inf) if same else None
    checks=p['checks'];checks.pop('four_centers_and_coordinates',None);checks.pop('coordinates',None)
    checks.update(all_centers_and_coordinates=same and error<1e-6,per_center_basis=counts==c['expected_basis_functions_by_center'],
        ghost_nuclear_charges=same and [x['effective_Z'] for x in rows]==c['expected_effective_nuclear_charges'])
    p['basis_functions_by_center']=counts;p['coordinate_max_error_A']=error
    if c['method']=='cc':
        tolerances=[cp.number(x) for x in re.findall(r'Convergence tolerance \(max\. residuum\)\s*\.{2,}\s*('+F+r')',native)]
        chunks=re.split(r'Iter\s+E\(tot\)',native,maxsplit=1);trace=chunks[1].split('COUPLED CLUSTER ENERGY',1)[0] if len(chunks)==2 else ''
        entries=[list(map(cp.number,x)) for x in re.findall(r'^\s*('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s*$',trace,re.M)]
        checks['declared_CC_residual']=len(tolerances)==1 and abs(tolerances[0]-1e-8)<1e-16 and bool(entries) and abs(entries[-1][3])<=1e-8
        p['CC_residual_trace']=entries;p['CC_printed_tolerance']=tolerances
        p['T1_diagnostic']=[cp.number(x) for x in re.findall(r'T1 diagnostic\s*\.{2,}\s*('+F+r')',native)]
        p['triples_hartree']=[cp.number(x) for x in re.findall(r'^Triples Correction \(T\)\s*\.{2,}\s*('+F+r')',native,re.M)]
        amps=native.split('LARGEST AMPLITUDES',1)[1].split('!!!',1)[0] if 'LARGEST AMPLITUDES' in native else ''
        p['largest_printed_amplitudes']=[cp.number(x) for x in re.findall(r'^\s*-?\d+->\s*-?\d+\s+-?\d+->\s*-?\d+\s+('+F+r')\s*$',amps,re.M)]
    p['native_scf_mode']=re.findall(r'Convergence Check Mode\s+ConvCheckMode\s*\.{2,}\s*([^\n]+)',native)
    p['native_scf_stop_signals']=re.findall(r'\*+\s*(.*?)\s+signals convergence\s*\*+',native)
    p['status']='complete' if all(checks.values()) else 'failed_admission';return p

def prepare(out):
    tick=time.monotonic();b.check(b.pin(PARENT)['sha256']==PARENT_SHA,'Matched geometry parent changed');old=b.load_panel(PARENT)
    b.check(b.pin(WATER)['sha256']==WATER_SHA,'Completed water baseline changed')
    out.mkdir(parents=True,exist_ok=False);(out/'assets').mkdir();(out/'inputs').mkdir();assets={};saved={};cells=[]
    for metal in METALS:
        assets[metal]={}
        for key in ('basis','aux'):
            dst=out/'assets'/f'{metal}_{key}.inc';shutil.copyfile(b.verify(old['assets'][metal][key]),dst);assets[metal][key]=b.pin(dst)
    for donor in DONORS:
        qlig=0 if donor=='acetamide' else -1;pair=[];saved[donor]={}
        for metal in METALS:
            oldcell=next(x for x in old['cells'] if x['geometry_id']==f'matched_dry_{donor}_central_{metal}');g=oldcell['geometry'];pair.append(g)
            for method in ('cc','pbe0'):
                for role in ROLES:
                    if method=='pbe0' and role=='complex':continue
                    n=len(g['symbols']);ghost={'complex':[False]*n,'ion_ghost_ligand':[False]+[True]*(n-1),'ligand_ghost_ion':[True]+[False]*(n-1)}[role]
                    q={'complex':3+qlig,'ion_ghost_ligand':3,'ligand_ghost_ion':qlig}[role];core=46 if metal=='La' else 55
                    za=[0 if isghost else (11 if i==0 else b.Z[s]) for i,(s,isghost) in enumerate(zip(g['symbols'],ghost))];ne=sum(za)-q
                    cells.append({'cell_id':f'{method}__{donor}__{metal}__{role}','method':method,'donor':donor,'metal':metal,'role':role,'d4':False,
                        'symbols':g['symbols'],'positions_A':g['positions_A'],'source_geometry':g['xyz'],'ghost_mask':ghost,'charge':q,'multiplicity':1,
                        'physical_multiplicity_of_parent':g['physical_multiplicity'],'expected_effective_nuclear_charges':za,'expected_electrons':ne,
                        'expected_correlated_electrons':ne,'expected_basis_functions':sum(NBF[x] for x in g['symbols']),
                        'expected_basis_functions_by_center':[NBF[x] for x in g['symbols']],'expected_ecp_core':0 if ghost[0] else core})
            path=PARENT_TASKS/oldcell['cell_id']/'execution/RESULT.json';r=json.loads(path.read_text());b.check(r['status']=='complete' and all(r['checks'].values()),'Saved complex unavailable')
            b.check(r['panel']==b.pin(PARENT) and r['geometry']==g['xyz'] and r['adapter']==old['adapter'],'Saved complex provenance differs')
            text=b.verify(r['output']).read_text();es,ed,ef=dft.electronic_parts(text);b.check(len(es)==len(ed)==len(ef)==1 and abs(es[0]+ed[0]-ef[0])<1e-8,'Saved electronic/D4 closure')
            saved[donor][metal]={'result':b.pin(path),'output':r['output'],'source_geometry':g['xyz'],'electronic_energy_hartree':es[0],'D4_energy_hartree':ed[0],'final_energy_hartree':ef[0]}
        b.check(pair[0]['positions_A']==pair[1]['positions_A'] and pair[0]['symbols'][1:]==pair[1]['symbols'][1:],'Paired geometry differs')
        g=pair[0];symbols=g['symbols'][1:];pos=g['positions_A'][1:];za=[b.Z[s] for s in symbols]
        cells.append({'cell_id':f'pbe0__{donor}__real_ligand_D4','method':'pbe0','donor':donor,'metal':None,'role':'real_ligand','d4':True,
            'symbols':symbols,'positions_A':pos,'source_geometry':g['xyz'],'ghost_mask':[False]*len(symbols),'charge':qlig,'multiplicity':1,
            'physical_multiplicity_of_parent':1,'expected_effective_nuclear_charges':za,'expected_electrons':sum(za)-qlig,'expected_correlated_electrons':sum(za)-qlig,
            'expected_basis_functions':sum(NBF[s] for s in symbols),'expected_basis_functions_by_center':[NBF[s] for s in symbols],'expected_ecp_core':0})
    panel={'protocol':PROTOCOL,'created_utc':b.now(),'adapter':b.pin(__file__),'specification':b.pin(HERE/'SPEC.md'),'parent_panel':b.pin(PARENT),
        'base_helper':old['adapter'],'canonical_helper':b.pin(HERE.parent/'canonical_water_cp_v2/REFERENCE.py'),'dft_helper':b.pin(HERE.parent/'dft_water_cp_v1/REFERENCE.py'),
        'water_baseline':b.pin(WATER),'orca':old['orca'],'assets':assets,'saved_complexes':saved,'cells':cells,'maximum_native_calls':22,'ranks':4,
        'scout_cell_id':'cc__acetamide__La__complex','isolated_ion_D4_hartree':0.,'ion_D4_evidence':'Exact empty-pair/triple sum; not calculated',
        'primary_exchange_resolution_kcal_mol':.25,'interpretation':'Same-basis frozen-f fixed interaction exchange; physical accuracy, forces and free energy unqualified'}
    for c in cells:
        path=out/'inputs'/(c['cell_id']+'.inp');path.write_text(deck(panel,c,4,3000));c['template_input']=b.pin(path)
    pin=b.write(out/'PANEL.json',panel);val=validate(pin['path']);val.update(command=sys.argv,elapsed_seconds=time.monotonic()-tick);b.write(out/'PREFLIGHT.json',val)
    print(json.dumps({'panel':pin,'status':val['status'],'cells':len(cells),'scout':panel['scout_cell_id'],'molecular_calls':0}))

def load(path):
    p=json.loads(Path(path).read_text());b.check(p['protocol']==PROTOCOL,'Wrong protocol')
    for key in ('adapter','specification','parent_panel','base_helper','canonical_helper','dft_helper','water_baseline','orca'):b.verify(p[key])
    b.check(p['adapter']==b.pin(__file__),'Executing adapter differs');b.check(len(p['cells'])==p['maximum_native_calls']==22,'Twenty-two-cell denominator differs')
    expected={f'cc__{d}__{m}__{r}' for d in DONORS for m in METALS for r in ROLES}|{f'pbe0__{d}__{m}__{r}' for d in DONORS for m in METALS for r in ROLES[1:]}|{f'pbe0__{d}__real_ligand_D4' for d in DONORS}
    b.check({c['cell_id'] for c in p['cells']}==expected,'Inventory differs')
    for a in p['assets'].values():
        for pin in a.values():b.verify(pin)
    for pair in p['saved_complexes'].values():
        for r in pair.values():
            for key in ('result','output','source_geometry'):b.verify(r[key])
    for c in p['cells']:b.verify(c['source_geometry']);b.verify(c['template_input'])
    return p

def validate(path):
    p=load(path);records=[]
    for c in p['cells']:
        raw=deck(p,c,4,3000);b.check(raw==b.verify(c['template_input']).read_text(),'Deck changed');symbols,x=b.xyz_read(c['source_geometry']['path']);offset=1 if c['role']=='real_ligand' else 0
        b.check(symbols[offset:]==c['symbols'] and x[offset:]==c['positions_A'],'Coordinates changed');header,*lines=raw.split('* xyz ',1)[1].splitlines();lines=lines[:-1];labels=[line.split()[0] for line in lines];q,mult=map(int,header.split())
        actualZ=sum((11 if s.rstrip(':') in METALS else b.Z[s]) for s in labels if not s.endswith(':'))
        expected_ne={'complex':40,'ion_ghost_ligand':8,'ligand_ghost_ion':32,'real_ligand':32}[c['role']]
        b.check(actualZ-q==expected_ne==c['expected_electrons'] and mult==1,'Independent electron inventory')
        b.check([s.endswith(':') for s in labels]==c['ghost_mask'],'Ghost mask differs');b.check(sum(c['expected_basis_functions_by_center'])==c['expected_basis_functions'],'AO bookkeeping')
        if c['method']=='cc':b.check('NoFrozenCore' in raw and 'STol 1e-8' in raw and not any(x in raw for x in ('D4','RIJCOSX','EnGrad')),'CC method differs')
        else:b.check(('D4' in raw)==c['d4'] and 'RIJCOSX' in raw and 'EnGrad' not in raw,'DFT method differs')
        b.check(('NewECP' in raw)==bool(c['expected_ecp_core']),'ECP inventory differs')
        if c['metal'] and c['ghost_mask'][0]:b.check('DelECP '+c['metal'] in raw and 'NewGTO '+c['metal'] in raw,'Ghost orbital basis/ECP differs')
        records.append({'cell_id':c['cell_id'],'atoms':len(labels),'charge':q,'electrons':expected_ne,'basis_functions':c['expected_basis_functions'],'ghost_centers':sum(c['ghost_mask'])})
    return {'status':'passed_static_preflight','panel':b.pin(path),'independent_deck_inventory':records,'molecular_calls':0,'local_threads':1,'utc':b.now()}

def run(path,cell_id,out):
    p=load(path);b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native work requires Slurm');cells=[c for c in p['cells'] if c['cell_id']==cell_id];b.check(len(cells)==1,'Unknown cell');c=cells[0]
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

def collect(path,paths):
    p=load(path);rows={}
    for filename in paths:
        r=json.loads(Path(filename).read_text());b.check(r['status']=='complete' and all(r['checks'].values()) and r['panel']==b.pin(path) and r['adapter']==p['adapter'],'Invalid native result');b.check(r['cell_id'] not in rows,'Duplicate result');b.verify(r['output']);b.verify(r['input']);rows[r['cell_id']]=r
    b.check(set(rows)=={x['cell_id'] for x in p['cells']},'Complete22 denominator required');water=json.loads(b.verify(p['water_baseline']).read_text());donors={}
    for donor in DONORS:
        pair={};ld4=rows[f'pbe0__{donor}__real_ligand_D4']['D4_energy_hartree']
        for metal in METALS:
            stem=f'cc__{donor}__{metal}__';cc=(rows[stem+'complex']['energy_hartree']-rows[stem+'ion_ghost_ligand']['energy_hartree']-rows[stem+'ligand_ghost_ion']['energy_hartree'])*b.HA
            stem=f'pbe0__{donor}__{metal}__';saved=p['saved_complexes'][donor][metal];pe=(saved['electronic_energy_hartree']-rows[stem+'ion_ghost_ligand']['electronic_energy_hartree']-rows[stem+'ligand_ghost_ion']['electronic_energy_hartree'])*b.HA;disp=(saved['D4_energy_hartree']-p['isolated_ion_D4_hartree']-ld4)*b.HA
            interaction={'CCSD_T_CP_kcal_mol':cc,'PBE0_CP_kcal_mol':pe,'PBE0_D4_CP_kcal_mol':pe+disp};exchange={key:value-water['per_metal'][metal][key] for key,value in interaction.items()}
            pair[metal]={'interaction':interaction,'water_to_donor_exchange':exchange,'D4_interaction_kcal_mol':disp}
        contrast={key:pair['Dy']['water_to_donor_exchange'][key]-pair['La']['water_to_donor_exchange'][key] for key in pair['La']['interaction']};error={key:contrast[key]-contrast['CCSD_T_CP_kcal_mol'] for key in ('PBE0_CP_kcal_mol','PBE0_D4_CP_kcal_mol')}
        donors[donor]={'per_metal':pair,'paired_exchange':contrast,'method_minus_CC_paired_exchange':error,'PBE0_D4_within_declared_resolution':abs(error['PBE0_D4_CP_kcal_mol'])<=p['primary_exchange_resolution_kcal_mol']}
    return {'status':'complete','panel':b.pin(path),'native_result_pins':[b.pin(x) for x in paths],'water_baseline':p['water_baseline'],'donors':donors,'molecular_calls':0,'physical_state_qualified':False,'method_accuracy_qualified':False}

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
