#!/usr/bin/env python3
"""Finite small-cluster La/Dy panel. Native calculations are Slurm-only."""
import argparse
import datetime as dt
import hashlib
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
ROOT = HERE.parents[2]
PROTOCOL = 'lady_reference_chemistry_v1'
Z = {'H':1,'C':6,'N':7,'O':8,'La':57,'Dy':66}
HA = 627.509474
BOHR = 0.529177210903
ORCA = Path('/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca')
ARCHIVE = ROOT/'workspaces/lanm_coupled_response_20261016/lbt_local_reference_remaining_v2/7CCO_Dy/execution/MANIFEST.json'


def check(ok, why):
    if not ok:
        raise ValueError(why)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as h:
        sha = hashlib.file_digest(h, 'sha256').hexdigest()
    return {'path':str(p), 'sha256':sha}


def verify(p):
    check(pin(p['path']) == p, 'Changed input: '+p['path'])
    return Path(p['path'])


def write(path, value):
    with Path(path).open('x') as h:
        json.dump(value, h, indent=2, allow_nan=False)
        h.write('\n')
    return pin(path)


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def add(a,b): return [x+y for x,y in zip(a,b)]
def sub(a,b): return [x-y for x,y in zip(a,b)]
def mul(a,b): return [x*b for x in a]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def norm(a): return math.sqrt(dot(a,a))
def unit(a): return mul(a,1/norm(a))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]


def xyz_write(path, symbols, positions, comment):
    with Path(path).open('x') as h:
        h.write(str(len(symbols))+'\n'+comment+'\n')
        for s,x in zip(symbols,positions): h.write(s+' '+' '.join('%.12f'%v for v in x)+'\n')
    return pin(path)


def xyz_read(path):
    a=Path(path).read_text().splitlines()
    rows=[r.split() for r in a[2:] if r.strip()]
    check(len(rows)==int(a[0]), 'XYZ atom count differs')
    return [r[0] for r in rows], [[float(v) for v in r[1:4]] for r in rows]


def published_structures():
    text=(HERE/'assets/jp2c07237_si_001.txt').read_text()
    names={'La3_H2O.xyz':'water','La3_CH3COO.xyz':'acetate_bi',
           'La3_CH3COO_mono_7cco.xyz':'acetate_mono','La3_CH3CONH2.xyz':'acetamide'}
    lines=text.splitlines(); out={}
    for i,line in enumerate(lines):
        t=line.split()
        if len(t)==2 and t[0].isdigit() and t[1] in names:
            n=int(t[0]); rows=[r.split() for r in lines[i+1:i+1+n]]
            check(len(rows)==n and [int(r[0]) for r in rows]==list(range(1,n+1)), 'Published atom ordering')
            out[names[t[1]]]={'symbols':[r[1] for r in rows], 'positions':[[float(v) for v in r[2:5]] for r in rows],
                               'tinker_types':[int(r[5]) for r in rows], 'bonds':[[int(v)-1 for v in r[6:]] for r in rows],
                               'source_lines':[i+2,i+1+n], 'source_filename':t[1]}
    check(set(out)==set(names.values()), 'Missing published dimer')
    return out


def aquo(n):
    z=0.5 if n==8 else 0.35
    directions=[]
    for ring in range(2):
        for j in range(4):
            ang=j*math.pi/2+ring*math.pi/4
            directions.append([math.sqrt(1-z*z)*math.cos(ang), math.sqrt(1-z*z)*math.sin(ang), z*(1 if ring==0 else -1)])
    if n==9: directions.append([0.,0.,1.])
    symbols=['La']; x=[[0.,0.,0.]]; groups=[0]; oh=.9572; half=math.radians(104.52/2)
    for i,u in enumerate(directions):
        tangent=unit(cross(u,[0.,1.,0.] if abs(u[2])>.9 else [0.,0.,1.]))
        oxygen=mul(u,2.5)
        symbols+=['O','H','H']; groups += [i+1]*3
        x += [oxygen,add(oxygen,mul(add(mul(u,math.cos(half)),mul(tangent,math.sin(half))),oh)),
                     add(oxygen,mul(add(mul(u,math.cos(half)),mul(tangent,-math.sin(half))),oh))]
    mindist=min(norm(sub(x[i],x[j])) for i in range(1,len(x)) for j in range(i+1,len(x)) if groups[i]!=groups[j])
    check(mindist>1., 'Generated intermolecular overlap')
    return symbols,x,mindist


def prepare(out):
    out.mkdir(parents=True,exist_ok=False)
    (out/'xyz').mkdir(); (out/'assets').mkdir(); (out/'cells').mkdir()
    old=json.loads(ARCHIVE.read_text()); assets={}
    for metal in ('La','Dy'):
        assets[metal]={}
        for kind,record in old['assets'][metal].items():
            src=verify(record); target=out/'assets'/f'{metal}_{kind}.inc'; shutil.copyfile(src,target)
            assets[metal][kind]=pin(target)
    check(pin(ORCA)==old['orca'], 'ORCA changed from native reference lineage')
    source=pin(HERE/'assets/jp2c07237_si_001.txt'); geometries=[]; fragments=[]
    for family,s in published_structures().items():
        pos=s['positions']; donors=[i for i,a in enumerate(s['symbols']) if a=='O']
        if family!='acetate_bi': donors=[min(donors,key=lambda i:norm(sub(pos[i],pos[0])))]
        center=[sum(pos[i][k] for i in donors)/len(donors) for k in range(3)]
        direction=unit(sub(pos[0],center))
        for metal in ('La','Dy'):
            for label,shift in [('minus',-.15),('central',0.),('plus',.15)]:
                symbols=[metal]+s['symbols'][1:]; xx=[add(pos[0],mul(direction,shift))]+pos[1:]
                gid=f'{family}_{label}_{metal}'; xp=xyz_write(out/'xyz'/(gid+'.xyz'),symbols,xx,gid)
                symbols,xx=xyz_read(xp['path'])
                q=2 if family.startswith('acetate') else 3
                geometries.append({'geometry_id':gid,'case_id':family,'distortion':label,'metal':metal,
                    'metal_index':0,'charge':q,'physical_multiplicity':1 if metal=='La' else 6,
                    'symbols':symbols,'positions_A':xx,'xyz':xp,'atoms':len(symbols),'source':source,
                    'source_tinker_filename':s['source_filename'],'radial_direction':direction,
                    'radial_shift_A':shift,'donor_indices':donors,'source_tinker_types':s['tinker_types'],
                    'source_bonds_zero_based':s['bonds'],'kind':'published_dimer'})
        fp=xyz_write(out/'xyz'/(family+'_free.xyz'),s['symbols'][1:],pos[1:],family+' fixed free fragment')
        fragments.append({'geometry_id':family+'_free','case_id':family,'charge':-1 if family.startswith('acetate') else 0,
            'physical_multiplicity':1,'metal_index':None,'xyz':fp,'atoms':len(pos)-1,
            'geometry_convention':'Exact internal geometry of corresponding published complex; no fragment optimization'})
    for n in (8,9):
        symbols,x,mindist=aquo(n)
        for metal in ('La','Dy'):
            ss=[metal]+symbols[1:]; gid=f'aquo{n}_central_{metal}'; xp=xyz_write(out/'xyz'/(gid+'.xyz'),ss,x,gid)
            ss,xx=xyz_read(xp['path'])
            geometries.append({'geometry_id':gid,'case_id':f'aquo{n}','distortion':'central','metal':metal,
                'metal_index':0,'charge':3,'physical_multiplicity':1 if metal=='La' else 6,'symbols':ss,
                'positions_A':xx,'xyz':xp,'atoms':len(ss),'water_count':n,'kind':'constructed_aquo',
                'minimum_interwater_distance_A':mindist,'recipe':'SPEC.md deterministic square-antiprism/rings, radius2.50A'})
    cells=[]
    for recipe in ('ecp_pbe0','ae_pbe0'):
        for g in geometries:
            tid=recipe+'__'+g['geometry_id']; core=46 if g['metal']=='La' else 55
            mult=1 if recipe=='ecp_pbe0' else g['physical_multiplicity']
            ne=sum(Z[s] for s in g['symbols'])-g['charge']-(core if recipe=='ecp_pbe0' else 0)
            check((ne-(mult-1))%2==0 and ne>=mult-1,'Electron/multiplicity parity')
            cells.append({'cell_id':tid,'recipe':recipe,'geometry_id':g['geometry_id'],'charge':g['charge'],
                'multiplicity':mult,'physical_multiplicity':g['physical_multiplicity'],'metal':g['metal'],
                'ecp_core':core if recipe=='ecp_pbe0' else 0,'explicit_electrons':ne,'geometry':g,
                'pilot':g['case_id'] in ('water','acetate_bi') and g['distortion']=='central'})
    panel={'protocol':PROTOCOL,'created_utc':now(),'adapter':pin(__file__),'specification':pin(HERE/'SPEC.md'),
        'published_coordinate_source':source,'archived_recipe_source':pin(ARCHIVE),'orca':pin(ORCA),'assets':assets,
        'geometries':geometries,'free_fragments':fragments,'cells':cells,'pilot_cell_ids':[c['cell_id'] for c in cells if c['pilot']],
        'interpretation':'Fixed-state chemistry; not equilibrium populations or binding free energy',
        'local_preparation':{'hostname':socket.gethostname(),'threads':1,'molecular_calls':0}}
    check(len(geometries)==28 and len(cells)==56 and len(panel['pilot_cell_ids'])==8,'Panel denominator')
    rec=write(out/'PANEL.json',panel)
    write(out/'PREPARATION.json',{'status':'complete','panel':rec,'geometries':28,'cells':56,'pilot_cells':8,'molecular_calls':0})
    print(json.dumps(rec))


def load_panel(path):
    p=json.loads(Path(path).read_text())
    check(p['protocol']==PROTOCOL,'Protocol differs'); verify(p['adapter']); verify(p['specification']); verify(p['orca'])
    check(p['adapter']==pin(__file__),'Executing adapter differs')
    for g in p['geometries']+p['free_fragments']: verify(g['xyz'])
    for a in p['assets'].values():
        for q in a.values(): verify(q)
    return p


def deck(p,c,ranks,maxcore):
    metal=c['metal']
    if c['recipe']=='ecp_pbe0':
        kw='PBE0 D4 def2-TZVP def2/J RIJCOSX'
        basis='\n'.join(verify(p['assets'][metal][k]).read_text() for k in ('basis','aux'))
    else:
        kw='PBE0 D4 ZORA ZORA-def2-TZVP SARC/J RIJCOSX NoECP'
        basis=f' NewGTO {metal} "SARC-ZORA-TZVP" end\n DelECP {metal}\n'
    return (f'! {kw} NoAutostart DefGrid3 VeryTightSCF EnGrad\n%maxcore {maxcore}\n'
            f'%pal nprocs {ranks} end\n%scf\n Guess PModel\n MaxIter 300\nend\n%basis\n{basis}\nend\n'
            f'* xyzfile {c["charge"]} {c["multiplicity"]} geometry.xyz\n')


def parse_engrad(path):
    lines=[r.strip() for r in Path(path).read_text().splitlines() if r.strip() and not r.lstrip().startswith('#')]
    n=int(lines[0]); energy=float(lines[1]); values=[float(r) for r in lines[2:2+3*n]]
    coords=[r.split() for r in lines[2+3*n:]]
    check(len(values)==3*n and len(coords)==n,'Engrad dimensions')
    check(all(math.isfinite(v) for v in [energy]+values),'Nonfinite energy/gradient')
    return n,energy,[values[3*i:3*i+3] for i in range(n)],[int(r[0]) for r in coords],[[float(v)*BOHR for v in r[1:4]] for r in coords]


def collect_one(p,c,out,returncode):
    text=(out/'endpoint.out').read_text(errors='replace'); g=c['geometry']; checks={}
    checks['process_success']=returncode==0
    checks['normal_termination']='ORCA TERMINATED NORMALLY' in text
    checks['scf_converged']=bool(re.search(r'SCF CONVERGED AFTER',text))
    numbers=re.findall(r'Number of Electrons\s+(?:NEL\s+)?\.{2,}\s*([0-9]+)',text)
    multiples=re.findall(r'Multiplicity\s+(?:Mult\s+)?\.{2,}\s*([0-9]+)',text)
    checks['electron_count']=bool(numbers) and int(numbers[-1])==c['explicit_electrons']
    checks['multiplicity']=bool(multiples) and int(multiples[-1])==c['multiplicity']
    if c['ecp_core']:
        checks['core_representation']=bool(re.search(r'replacing\s+'+str(c['ecp_core'])+r'\s+core electrons',text))
    else:
        checks['core_representation']=('ZORA' in text and not re.search(r'replacing\s+\d+\s+core electrons',text))
    n,e,grad,z,x=parse_engrad(out/'endpoint.engrad')
    checks['atom_identity']=n==g['atoms'] and z==[Z[s] for s in g['symbols']]
    checks['coordinates']=max(abs(a-b) for aa,bb in zip(x,g['positions_A']) for a,b in zip(aa,bb))<1e-7
    energies=re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+0-9.Ee]+)',text)
    checks['energy_matches_gradient']=bool(energies) and abs(float(energies[-1])-e)<1e-8
    spin=re.findall(r'Expectation value of <S\*\*2>\s*[:.]?\s*([-+0-9.Ee]+)',text)
    criteria=[]
    for line in text.splitlines():
        if re.search(r'(Energy change|MAX-Density change|RMS-Density change|DIIS Error|Orbital Gradient)',line) and ('YES' in line or 'NO' in line): criteria.append(line.strip())
    force=[[-v*HA/BOHR for v in row] for row in grad]
    result={'status':'complete' if all(checks.values()) else 'failed_admission','cell_id':c['cell_id'],'checks':checks,
        'energy_hartree':e,'gradient_hartree_bohr':grad,'forces_kcal_mol_A':force,
        'radial_metal_force_kcal_mol_A':dot(force[0],g['radial_direction']) if 'radial_direction' in g else None,
        'spin_squared_values':[float(v) for v in spin], 'printed_scf_criteria':criteria,
        'physical_state_qualified':False,'method_accuracy_qualified':False,'output':pin(out/'endpoint.out'),
        'engrad':pin(out/'endpoint.engrad'),'input':pin(out/'endpoint.inp'),'geometry':g['xyz']}
    return result


def run(panel_path,cell_id,out):
    check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native molecular calls require Slurm')
    p=load_panel(panel_path); matches=[c for c in p['cells'] if c['cell_id']==cell_id]; check(len(matches)==1,'Unknown/duplicate cell')
    c=matches[0]; ranks=int(os.environ.get('SLURM_NTASKS','1')); cpus=int(os.environ.get('SLURM_CPUS_PER_TASK','1'))
    check(1<=ranks<=24 and cpus==1,'Expected one numerical thread per MPI rank')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0'))
    if not mem: mem=int(os.environ.get('SLURM_MEM_PER_CPU','0'))*ranks
    check(mem>0,'Actual granted memory required'); maxcore=math.floor(.75*mem/ranks)
    out.mkdir(parents=True,exist_ok=False); shutil.copyfile(verify(c['geometry']['xyz']),out/'geometry.xyz')
    (out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore))
    env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'): env[key]='1'
    env['PATH']=str(ORCA.parent)+os.pathsep+env.get('PATH',''); env['OMP_DYNAMIC']='FALSE'; env['MKL_DYNAMIC']='FALSE'
    start=time.monotonic(); cmd=[str(verify(p['orca'])),'endpoint.inp','--bind-to none']
    write(out/'START.json',{'cell':c,'panel':pin(panel_path),'utc':now(),'host':socket.gethostname(),'job_id':os.environ['SLURM_JOB_ID'],
        'ranks':ranks,'allocated_cpus':ranks,'memory_MiB':mem,'maxcore_MiB':maxcore,'command':cmd,
        'environment':{k:env[k] for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','PATH')},'input':pin(out/'endpoint.inp')})
    with (out/'endpoint.out').open('x') as h, (out/'endpoint.err').open('x') as err:
        proc=subprocess.run(cmd,cwd=out,env=env,stdout=h,stderr=err,check=False)
    elapsed=time.monotonic()-start
    try: result=collect_one(p,c,out,proc.returncode)
    except Exception as e: result={'status':'failed_collection','cell_id':cell_id,'exception':repr(e)}
    result.update({'returncode':proc.returncode,'elapsed_seconds':elapsed,'allocated_cpu_seconds':elapsed*ranks,
        'job_id':os.environ['SLURM_JOB_ID'],'finished_utc':now(),'panel':pin(panel_path),'adapter':pin(__file__)})
    write(out/'RESULT.json',result); print(json.dumps({'status':result['status'],'cell_id':cell_id,'seconds':elapsed}))
    return 0 if result['status']=='complete' else 1


def validate(panel_path):
    p=load_panel(panel_path); checks=[]
    for g in p['geometries']:
        s,x=xyz_read(verify(g['xyz'])); check(s==g['symbols'] and x==g['positions_A'],'Coordinates differ')
        peer=[h for h in p['geometries'] if h['case_id']==g['case_id'] and h['distortion']==g['distortion'] and h['metal']!=g['metal']]
        check(len(peer)==1 and peer[0]['positions_A']==x,'La/Dy geometry mismatch')
        check(len(set((h['charge'], h['atoms']) for h in [g,peer[0]]))==1,'Matched charge/atom count')
        checks.append(g['geometry_id'])
    for c in p['cells']:
        t=deck(p,c,4,3000); check('EnGrad' in t and 'MaxIter 300' in t and '* xyzfile' in t,'Incomplete deck')
        check(('N_core '+str(c['ecp_core']) in t) if c['ecp_core'] else ('NoECP' in t and 'ZORA' in t),'Core deck mismatch')
    return {'status':'complete','geometries_checked':len(checks),'cells_rendered':len(p['cells']),'pilot_cells':p['pilot_cell_ids'],
        'molecular_calls':0,'panel':pin(panel_path),'adapter':pin(__file__),'utc':now(),'local_threads':1}


def collect(panel_path,result_paths):
    p=load_panel(panel_path); rows={}; unavailable=[]
    for path in result_paths:
        r=json.loads(Path(path).read_text()); check(r['panel']==pin(panel_path),'Result panel differs')
        check(r['cell_id'] not in rows,'Duplicate result'); rows[r['cell_id']]=r
        if r['status']!='complete': unavailable.append(r['cell_id'])
    contrasts=[]
    for recipe in ('ecp_pbe0','ae_pbe0'):
        for donor in ('acetate_bi','acetate_mono','acetamide'):
            for dist in ('minus','central','plus'):
                ids=[f'{recipe}__{donor}_{dist}_Dy',f'{recipe}__water_{dist}_Dy',f'{recipe}__{donor}_{dist}_La',f'{recipe}__water_{dist}_La']
                if all(k in rows and rows[k]['status']=='complete' for k in ids):
                    e=[rows[k]['energy_hartree'] for k in ids]
                    contrasts.append({'recipe':recipe,'donor':donor,'distortion':dist,'cells':ids,
                        'Dy_minus_La_water_to_donor_kcal_mol':((e[0]-e[1])-(e[2]-e[3]))*HA})
    coordination=[]
    for recipe in ('ecp_pbe0','ae_pbe0'):
        ids=[f'{recipe}__aquo9_central_Dy',f'{recipe}__aquo8_central_Dy',f'{recipe}__aquo9_central_La',f'{recipe}__aquo8_central_La']
        if all(k in rows and rows[k]['status']=='complete' for k in ids):
            e=[rows[k]['energy_hartree'] for k in ids]
            coordination.append({'recipe':recipe,'cells':ids,'Dy_minus_La_8_to_9_coordination_kcal_mol':((e[0]-e[1])-(e[2]-e[3]))*HA})
    core_sensitivity=[]
    for c in contrasts:
        if c['recipe']!='ecp_pbe0': continue
        peer=[r for r in contrasts if r['recipe']=='ae_pbe0' and r['donor']==c['donor'] and r['distortion']==c['distortion']]
        if peer:
            core_sensitivity.append({'donor':c['donor'],'distortion':c['distortion'],
                'AE_minus_ECP_exchange_contrast_kcal_mol':peer[0]['Dy_minus_La_water_to_donor_kcal_mol']-c['Dy_minus_La_water_to_donor_kcal_mol']})
    return {'status':'partial' if len(rows)<56 or unavailable else 'complete','declared':56,'available':sum(r['status']=='complete' for r in rows.values()),
        'failed':unavailable,'contrasts':contrasts,'coordination_contrasts':coordination,'core_sensitivity':core_sensitivity,
        'result_pins':[pin(q) for q in result_paths],'panel':pin(panel_path),
        'native_calls_this_collection':0,'physical_state_qualified':False,'binding_free_energy_claimed':False,'utc':now()}


def main():
    a=argparse.ArgumentParser(); sub=a.add_subparsers(dest='action',required=True)
    b=sub.add_parser('prepare'); b.add_argument('--out',type=Path,required=True)
    b=sub.add_parser('validate'); b.add_argument('--panel',required=True); b.add_argument('--out',type=Path,required=True)
    b=sub.add_parser('run'); b.add_argument('--panel',required=True); b.add_argument('--cell',required=True); b.add_argument('--out',type=Path,required=True)
    b=sub.add_parser('collect'); b.add_argument('--panel',required=True); b.add_argument('--out',type=Path,required=True); b.add_argument('results',nargs='+')
    c=a.parse_args()
    if c.action=='prepare': prepare(c.out)
    elif c.action=='validate': write(c.out,validate(c.panel)); print('Validation passed; zero molecular calls')
    elif c.action=='run': return run(c.panel,c.cell,c.out)
    elif c.action=='collect': write(c.out,collect(c.panel,c.results))
    return 0


if __name__=='__main__': sys.exit(main())
