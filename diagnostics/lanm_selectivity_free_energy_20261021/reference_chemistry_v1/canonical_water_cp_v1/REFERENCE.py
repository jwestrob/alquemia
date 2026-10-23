#!/usr/bin/env python3
"""Six canonical water CP cells. Preparation is local; native calls require Slurm."""
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
spec=importlib.util.spec_from_file_location('frozen_reference',HERE.parent/'REFERENCE.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
ROOT=base.ROOT
PARENT=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/prepared_v2/PANEL.json'
PARENT_SHA='d5929526cf70ac957c847a5936ed229e9edf1322764b0cc287a59292e01204f9'
PROTOCOL='canonical_lady_water_cp_v1'
ROLES=('complex','ion_ghost_water','water_ghost_ion')
FLOAT=r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'

def number(text):return float(text.replace('D','E').replace('d','e'))

def metal_blocks(text,metal):
    blocks=re.findall(r'^New(?:ECP|GTO)\s+'+metal+r'\s*\n.*?^end\s*$',text,re.M|re.S)
    base.check(len(blocks)==2 and blocks[0].startswith('NewECP') and blocks[1].startswith('NewGTO'),'Unexpected archived basis structure')
    return blocks

def deck(p,c,ranks,maxcore):
    full=base.verify(p['assets'][c['metal']]).read_text();ecp,gto=metal_blocks(full,c['metal'])
    basis=full if c['role']!='water_ghost_ion' else 'DelECP '+c['metal']+'\n'+gto+'\n'
    coords=[]
    for i,(symbol,x) in enumerate(zip(c['symbols'],c['positions_A'])):
        coords.append(symbol+(':' if c['ghost_mask'][i] else '')+' '+' '.join('%.12f'%v for v in x))
    return (f'! RHF CCSD(T) def2-TZVP NoFrozenCore VeryTightSCF NoAutostart\n%maxcore {maxcore}\n'
        f'%pal nprocs {ranks} end\n%scf\n Guess PModel\n MaxIter 300\nend\n'
        '%mdci\n MaxIter 100\n STol 1e-8\n PrintLevel 3\nend\n'
        '%basis\n'+basis+'\nend\n'+f'* xyz {c["charge"]} 1\n'+'\n'.join(coords)+'\n*\n')

def prepare(out):
    tick=time.monotonic();base.check(base.pin(PARENT)['sha256']==PARENT_SHA,'Parent panel changed');old=base.load_panel(PARENT)
    out.mkdir(parents=True,exist_ok=False);(out/'assets').mkdir();(out/'inputs').mkdir()
    assets={};cells=[]
    for metal in ('La','Dy'):
        src=base.verify(old['assets'][metal]['basis']);dst=out/'assets'/(metal+'_basis.inc');shutil.copyfile(src,dst);assets[metal]=base.pin(dst)
        g=next(g for g in old['geometries'] if g['geometry_id']=='water_central_'+metal)
        core=46 if metal=='La' else 55
        for role in ROLES:
            ghost={'complex':[False]*4,'ion_ghost_water':[False,True,True,True],'water_ghost_ion':[True,False,False,False]}[role]
            charge=0 if role=='water_ghost_ion' else 3
            za=[0 if ghost[i] else base.Z[s]-(core if i==0 else 0) for i,s in enumerate(g['symbols'])]
            ne=sum(za)-charge
            cells.append({'cell_id':metal+'__'+role,'metal':metal,'role':role,'symbols':g['symbols'],'positions_A':g['positions_A'],
                'source_geometry':g['xyz'],'ghost_mask':ghost,'charge':charge,'multiplicity':1,
                'physical_multiplicity_of_parent':g['physical_multiplicity'],'expected_effective_nuclear_charges':za,
                'expected_electrons':ne,'expected_correlated_electrons':ne,'expected_basis_functions':90,
                'expected_basis_functions_by_center':[47,31,6,6],'expected_ecp_core':0 if ghost[0] else core})
    panel={'protocol':PROTOCOL,'created_utc':base.now(),'adapter':base.pin(__file__),'helper':old['adapter'],
        'specification':base.pin(HERE/'SPEC.md'),'parent_panel':base.pin(PARENT),'orca':old['orca'],
        'assets':assets,'cells':cells,'maximum_native_calls':6,'ranks':4,'local_preparation_threads':1,
        'resources':'standard_expand,4MPI,default memory/time,MaxCore75percent actual grant',
        'interpretation':'Frozen-f canonical fixed-geometry water interaction; no state/basis-limit/force/affinity qualification'}
    for c in cells:
        f=out/'inputs'/(c['cell_id']+'.inp');f.write_text(deck(panel,c,4,3000));c['template_input']=base.pin(f)
    pin=base.write(out/'PANEL.json',panel);validation=validate(pin['path'])
    validation.update(command=sys.argv,elapsed_seconds=time.monotonic()-tick)
    base.write(out/'PREFLIGHT.json',validation);print(json.dumps({'panel':pin,'preflight':validation['status'],'cells':6,'molecular_calls':0}))

def load(path):
    p=json.loads(Path(path).read_text());base.check(p['protocol']==PROTOCOL,'Protocol differs')
    for key in ('adapter','helper','specification','parent_panel','orca'):base.verify(p[key])
    base.check(p['adapter']==base.pin(__file__),'Executing source differs')
    for a in p['assets'].values():base.verify(a)
    base.check(len(p['cells'])==p['maximum_native_calls']==6,'Six-cell denominator differs')
    base.check({(c['metal'],c['role']) for c in p['cells']}=={(m,r) for m in ('La','Dy') for r in ROLES},'Inventory differs')
    for c in p['cells']:base.verify(c['source_geometry']);base.verify(c['template_input'])
    return p

def validate(path):
    p=load(path);checks=[]
    for c in p['cells']:
        text=deck(p,c,4,3000);base.check(text==base.verify(c['template_input']).read_text(),'Input template differs')
        symbols,x=base.xyz_read(c['source_geometry']['path']);base.check(symbols==c['symbols'] and x==c['positions_A'],'Input geometry differs')
        lines=text.split('* xyz ')[1].splitlines()[1:-1];base.check(len(lines)==4,'Four basis centers required')
        mask=[line.split()[0].endswith(':') for line in lines];base.check(mask==c['ghost_mask'],'Ghost notation differs')
        base.check(c['expected_electrons'] in (8,10,18) and c['expected_electrons']==sum(c['expected_effective_nuclear_charges'])-c['charge'],'Electron bookkeeping')
        base.check(c['expected_correlated_electrons']==c['expected_electrons'],'Frozen-core mismatch')
        base.check('NoFrozenCore' in text and 'CCSD(T)' in text and not any(k in text for k in ('D4','RIJCOSX','EnGrad','DLPNO','Opt ')),'Method mismatch')
        if c['ghost_mask'][0]:base.check('NewECP' not in text and 'DelECP '+c['metal'] in text,'Ghost metal retains ECP')
        else:base.check('N_core '+str(c['expected_ecp_core']) in text,'Real metal core differs')
        gto=metal_blocks(base.verify(p['assets'][c['metal']]).read_text(),c['metal'])[1]
        base.check(gto in text,'Metal orbital basis changed');checks.append(c['cell_id'])
    return {'status':'passed_static_preflight','panel':base.pin(path),'cells':checks,'molecular_calls':0,'local_threads':1,
        'native_syntax_and_ghost_semantics_validated':False,'utc':base.now()}

def shell_count(pattern):
    angular={s:i for i,s in enumerate('spdfghi')}
    return sum(int(n)*(2*angular[s]+1) for n,s in re.findall(r'(\d+)([spdfghi])',pattern.lower()))

def parse(text,c,returncode):
    # Exclude the echoed input from all native evidence checks.
    native=text.split('****END OF INPUT****')[-1]
    au=native.split('CARTESIAN COORDINATES (A.U.)',1)[-1].split('INTERNAL COORDINATES',1)[0]
    rows=[]
    for line in au.splitlines():
        match=re.match(r'\s*(\d+)\s+([A-Za-z]+):?\s+('+FLOAT+r')\*?\s+(\d+)\s+('+FLOAT+r')\s+('+FLOAT+r')\s+('+FLOAT+r')\s+('+FLOAT+r')\s*$',line)
        if match:
            i,s,z,frag,mass,x,y,zz=match.groups();rows.append({'index':int(i),'symbol':s,'effective_Z':number(z),'positions_A':[number(v)*base.BOHR for v in (x,y,zz)]})
    basis_section=native.split('BASIS SET INFORMATION',1)[-1].split('ECP PARAMETER INFORMATION',1)[0].split('SHARK',1)[0]
    groups={int(i):{'symbol':s,'spherical_functions':shell_count(sh),'contraction':sh} for i,s,sh in re.findall(r'Group\s+(\d+)\s+Type\s+([A-Za-z]+):?\s*:\s*[^\n]*?contracted to\s+([0-9spdfghi]+)',basis_section,re.I)}
    assigned={int(i):int(g) for i,g in re.findall(r'Atom\s+(\d+)\s*[A-Za-z]+:?\s+basis set group\s*=>\s*(\d+)',basis_section)}
    counts=[groups[assigned[i]]['spherical_functions'] if i in assigned and assigned[i] in groups else None for i in range(4)]
    basis_numbers=[int(v) for v in re.findall(r'Number of basis functions\s*\.{2,}\s*(\d+)',native,re.I)]
    total=[int(v) for v in re.findall(r'Number of [Ee]lectrons\s+(?:NEL\s+)?\.{2,}\s*(\d+)',native)]
    corr=[int(v) for v in re.findall(r'Number of correlated electrons\s*\.{2,}\s*(\d+)',native,re.I)]
    mult=[int(v) for v in re.findall(r'Multiplicity\s+(?:Mult\s+)?\.{2,}\s*(\d+)',native)]
    charge=[int(v) for v in re.findall(r'Total Charge\s+(?:Charge\s+)?\.{2,}\s*([-+]?\d+)',native)]
    ecp_assign=[int(v) for v in re.findall(r'Atom\s+(\d+)\s*[A-Za-z]+\s+ECP group\s*=>',native)]
    cores=[int(v) for v in re.findall(r'replacing\s+(\d+)\s+core electrons',native)]
    finals=[number(v) for v in re.findall(r'FINAL SINGLE POINT ENERGY\s+('+FLOAT+r')',native)]
    cc=[number(v) for v in re.findall(r'E\(CCSD\(T\)\)\s*(?:\.{2,}|=|:)\s*('+FLOAT+r')',native)]
    coordinate_ok=len(rows)==4 and [r['index'] for r in rows]==list(range(4)) and [r['symbol'] for r in rows]==c['symbols']
    coord_error=max((abs(a-b) for r,xx in zip(rows,c['positions_A']) for a,b in zip(r['positions_A'],xx)),default=math.inf) if coordinate_ok else None
    checks={'process_success':returncode==0,'normal_termination':'ORCA TERMINATED NORMALLY' in native,
        'scf_converged':bool(re.search(r'SCF CONVERGED AFTER',native)),
        'cc_converged':bool(re.search(r'(?:COUPLED[ -]CLUSTER|CCSD).*?(?:HAS\s+)?CONVERGED',native,re.I)) and not bool(re.search(r'COUPLED[ -]CLUSTER[^\n]*(?:NOT CONVERGED|NOT converged)',native,re.I)),
        'RHF_reference':bool(re.search(r'Reference Wavefunction\s*\.{2,}\s*RHF',native,re.I)),
        'charge':bool(charge) and all(q==c['charge'] for q in charge),
        'multiplicity':bool(mult) and all(v==1 for v in mult),
        'total_electrons':bool(total) and all(v==c['expected_electrons'] for v in total),
        'correlated_electrons':bool(corr) and all(v==c['expected_correlated_electrons'] for v in corr),
        'basis_function_count':bool(basis_numbers) and all(v==90 for v in basis_numbers),
        'per_center_basis':counts==c['expected_basis_functions_by_center'],
        'four_centers_and_coordinates':coordinate_ok and coord_error<1e-6,
        'ghost_nuclear_charges':coordinate_ok and [r['effective_Z'] for r in rows]==c['expected_effective_nuclear_charges'],
        'ecp_assignment':(ecp_assign==[0] and cores==[c['expected_ecp_core']]) if c['expected_ecp_core'] else (not ecp_assign and not cores),
        'final_CC_energy':bool(finals) and bool(cc) and math.isfinite(cc[-1]) and abs(finals[-1]-cc[-1])<1e-8}
    criteria=[]
    for name,value,tolerance in re.findall(r'Last\s+(.+?)\s+\.{3}\s+('+FLOAT+r')\s+Tolerance\s*:\s*('+FLOAT+r')',native):
        criteria.append({'name':name.strip(),'value':number(value),'tolerance':number(tolerance),'passes':abs(number(value))<=number(tolerance)})
    diag=[{'line':i+1,'text':line} for i,line in enumerate(text.splitlines()) if re.search(r'diagnostic|largest.*amplit|triples|E\(CCSD|residual|STol|Frozen core treatment|Internal Orbitals|Number of correlated electrons|Reference Wavefunction',line,re.I)]
    return {'status':'complete' if all(checks.values()) else 'failed_admission','checks':checks,
        'energy_hartree':cc[-1] if cc else None,'printed_final_energies_hartree':finals,'printed_CCSD_T_energies_hartree':cc,
        'atomic_table':rows,'coordinate_max_error_A':coord_error,'basis_groups':groups,'basis_functions_by_center':counts,
        'total_electrons':total,'correlated_electrons':corr,'ecp_assigned_centers':ecp_assign,'ecp_core_counts':cores,
        'printed_scf_criteria':criteria,'native_diagnostic_lines':diag,'physical_state_qualified':False,'method_accuracy_qualified':False}

def run(path,cell_id,out):
    base.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native calls require Slurm');p=load(path)
    chosen=[c for c in p['cells'] if c['cell_id']==cell_id];base.check(len(chosen)==1,'Unknown cell');c=chosen[0]
    ranks=int(os.environ.get('SLURM_NTASKS','0'));base.check(ranks==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread MPI ranks required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));base.check(mem>0,'Actual memory grant required')
    maxcore=math.floor(.75*mem/ranks);out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore))
    env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[key]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(base.verify(p['orca']).parent)+os.pathsep+env.get('PATH','')
    cmd=[p['orca']['path'],'endpoint.inp','--bind-to none'];tick=time.monotonic()
    base.write(out/'START.json',{'utc':base.now(),'cell':c,'panel':base.pin(path),'input':base.pin(out/'endpoint.inp'),
        'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'ranks':ranks,'memory_MiB':mem,'maxcore_MiB':maxcore,'command':cmd})
    with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:
        proc=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False)
    elapsed=time.monotonic()-tick
    try:r=parse((out/'endpoint.out').read_text(errors='replace'),c,proc.returncode)
    except Exception as e:r={'status':'failed_collection','exception':repr(e)}
    r.update(cell_id=cell_id,cell=c,panel=base.pin(path),adapter=base.pin(__file__),returncode=proc.returncode,
        elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,job_id=os.environ['SLURM_JOB_ID'],finished_utc=base.now(),
        output=base.pin(out/'endpoint.out'),stderr=base.pin(out/'endpoint.err'),input=base.pin(out/'endpoint.inp'),
        artifact_sizes={x.name:x.stat().st_size for x in out.iterdir() if x.is_file()})
    base.write(out/'RESULT.json',r);print(json.dumps({'status':r['status'],'cell_id':cell_id,'seconds':elapsed}))
    return 0 if r['status']=='complete' else 1

def collect(path,paths):
    p=load(path);rows={}
    for filename in paths:
        r=json.loads(Path(filename).read_text());base.check(r['panel']==base.pin(path),'Different panel');base.check(r['cell_id'] not in rows,'Duplicate result')
        base.check(r['status']=='complete' and all(r['checks'].values()),'Native cell unavailable');rows[r['cell_id']]=r
    base.check(set(rows)=={c['cell_id'] for c in p['cells']},'Full six-cell denominator required')
    energies={m:(rows[m+'__complex']['energy_hartree']-rows[m+'__ion_ghost_water']['energy_hartree']-rows[m+'__water_ghost_ion']['energy_hartree'])*base.HA for m in ('La','Dy')}
    return {'status':'complete','panel':base.pin(path),'result_pins':[base.pin(x) for x in paths],
        'CP_interaction_kcal_mol':energies,'Dy_minus_La_CP_water_interaction_kcal_mol':energies['Dy']-energies['La'],
        'molecular_calls':0,'affinity_claim':False,'physical_state_qualified':False,'method_accuracy_qualified':False}

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('validate');a.add_argument('--panel',required=True);a.add_argument('--out',required=True)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--cell',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--out',required=True);a.add_argument('results',nargs='+')
    a=ap.parse_args()
    if a.action=='prepare':prepare(a.out)
    elif a.action=='validate':base.write(a.out,validate(a.panel))
    elif a.action=='run':return run(a.panel,a.cell,a.out)
    elif a.action=='collect':base.write(a.out,collect(a.panel,a.results))

if __name__=='__main__':sys.exit(main())
