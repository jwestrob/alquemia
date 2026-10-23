#!/usr/bin/env python3
"""Finite 18 canonical-overlap controls + 8 local cluster energies; root submits."""
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
s=importlib.util.spec_from_file_location('local_cc_parent',HERE.parent/'coordination_cp_v1/REFERENCE.py')
h=importlib.util.module_from_spec(s);s.loader.exec_module(h)
b,cp,donor=h.b,h.cp,h.donor
read=h.read;F=cp.FLOAT;PROTOCOL='coordination_local_cc_26_v1'
TIERS={'tight':{'TCutPairs':1e-5,'TCutPNO':1e-7,'TCutMKN':1e-3,'TCutDO':5e-3},
       'tighter':{'TCutPairs':1e-6,'TCutPNO':1e-8,'TCutMKN':1e-4,'TCutDO':5e-3}}
SCOUT='cluster__tight__aquo9_central_Dy'
HELPERS={**h.HELPERS,'coordination_helper':HERE.parent/'coordination_cp_v1/REFERENCE.py'}

def deck(p,c,ranks,maxcore):
    m=c['metal'];full=b.verify(p['assets'][m]['basis']).read_text();ecp,gto=cp.metal_blocks(full,m)
    basis=('DelECP '+m+'\n'+gto+'\n') if c['ghost_mask'][0] else full
    basis+='\n'+b.verify(p['assets'][m]['aux_c']).read_text()
    opts='\n'.join(f' {k} {v:.12g}' for k,v in TIERS[c['tier']].items())
    text=(f'! RHF DLPNO-CCSD(T1) TightPNO def2-TZVP def2-TZVP/C NoFrozenCore VeryTightSCF NoAutostart PrintBasis\n'
          f'%maxcore {maxcore}\n%pal nprocs {ranks} end\n%scf\n Guess PModel\n MaxIter 300\nend\n'
          '%mdci\n MaxIter 100\n STol 1e-8\n PrintLevel 3\n TriTolE 1e-8\n'+opts+'\nend\n%basis\n'+basis+'\nend\n')
    text+=f'* xyz {c["charge"]} 1\n'
    text+='\n'.join(s+(':' if g else '')+' '+' '.join('%.12f'%v for v in x)
        for s,x,g in zip(c['symbols'],c['positions_A'],c['ghost_mask']))+'\n*\n'
    return text

def auxiliary(text,n):
    native=text.split('****END OF INPUT****')[-1]
    match=re.search(r'AUXILIARY/C BASIS SET INFORMATION\s*[-]+(.*?)(?=\n[-]{5,}\n|\nECP PARAMETER|\nSHARK|\Z)',native,re.S)
    section=match.group(1) if match else ''
    groups={int(i):sh for i,sh in re.findall(r'Group\s+(\d+)\s+Type\s+[A-Za-z]+:?\s*:\s*[^\n]*?contracted to\s+([0-9spdfghi]+)',section,re.I)}
    assigned={int(i):int(g) for i,g in re.findall(r'Atom\s+(\d+)\s*[A-Za-z]+:?\s+basis set group\s*=>\s*(\d+)',section)}
    contractions=[groups.get(assigned.get(i)) for i in range(n)]
    counts=[cp.shell_count(x) if x else None for x in contractions]
    totals=[int(x) for x in re.findall(r'(?:# of basis functions in Aux-C|Number of AuxC basis functions)\s*\.{2,}\s*(\d+)',native)]
    return dict(contractions=contractions,counts=counts,totals=totals)

def iteration_rows(text,header,columns):
    parts=re.split(header,text)
    out=[]
    for part in parts[1:]:
        rows=[];started=False
        for line in part.splitlines()[1:]:
            tokens=line.split()
            if len(tokens) in columns and all(re.fullmatch(F,x) for x in tokens):
                rows.append([cp.number(x) for x in tokens]);started=True
            elif started and (re.search(r'converged|ENERGY|TIMINGS|TRIPLES|correction|^-{5}',line,re.I)):
                break
        out.append(rows)
    return out

def parse(text,c,returncode):
    # Preserve canonical physical inventory parser; replace only method-specific traces.
    q=donor.parse(text,dict(c,method='cc'),returncode);native=text.split('****END OF INPUT****')[-1]
    checks=q['checks'];checks.pop('declared_CC_residual',None)
    traces=iteration_rows(native,r'Iter\s+E\(tot\)[^\n]*', (5,6))
    valid=[x for x in traces if x]
    ccres=[x[-1][-2] for x in valid]
    tol=[cp.number(x) for x in re.findall(r'Convergence tolerance \(max\. residuum\)\s*\.{2,}\s*('+F+r')',native)]
    triples=iteration_rows(native,r'Iter\s+E\(Triple\)[^\n]*',(4,))
    tri_tail=[x[-1][2] for x in triples if x]
    tri_section=native.split('Parameters used in iterative triples',1)[-1]
    etol=[cp.number(x) for x in re.findall(r'^ETol\s*\.{2,}\s*('+F+r')',tri_section,re.M)]
    energies=[cp.number(x) for x in re.findall(r'^E\(CCSD\)\s*\.{2,}\s*('+F+r')',native,re.M)]
    correction=[cp.number(x) for x in re.findall(r'^Triples Correction \(T\)\s*\.{2,}\s*('+F+r')',native,re.M)]
    finals=q['printed_final_energies_hartree'];local={}
    for key,value in TIERS[c['tier']].items():
        vals=[cp.number(x) for x in re.findall(r'^\s*'+key+r'\s*(?:=|:)\s*('+F+r')',native,re.M)]
        local[key]=vals
        checks['local_'+key]=bool(vals) and all(abs(v-value)<=abs(value)*1e-9 for v in vals)
    aux=auxiliary(native,len(c['symbols']));auxcounts=aux['counts']
    checks.update(declared_CC_residual=bool(tol) and all(abs(v-1e-8)<1e-16 for v in tol) and bool(ccres) and all(abs(v)<=1e-8 for v in ccres),
        iterative_triples=bool(triples) and all(triples) and bool(etol) and all(abs(v-1e-8)<1e-16 for v in etol) and all(abs(v)<=1e-8 for v in tri_tail),
        CCSD_plus_triples_closure=len(energies)==len(correction)==len(finals)==1 and abs(energies[0]+correction[0]-finals[0])<=1e-8,
        local_MP2_converged=bool(re.search(r'LMP2\s+iterations\s+CONVERGED',native)),
        localization_evidence=bool(re.search(r'(?:localiz[^\n]*(?:converged|completed|finished)|(?:converged|completed|finished)[^\n]*localiz)',native,re.I)) and not bool(re.search(r'localiz[^\n]*(?:not converged|failed)',native,re.I)),
        complete_aux_C=bool(aux['totals']) and all(isinstance(x,int) and x>0 for x in auxcounts) and all(x==sum(auxcounts) for x in aux['totals']),
        exact_metal_aux_C=auxcounts[0]=={'La':351,'Dy':369}[c['metal']],
        no_D4='DFT DISPERSION CORRECTION' not in native,
        no_native_failure=not bool(re.search(r'\b(?:NaN|Infinity)\b|error termination|(?:iterations.*NOT.*CONVERGED)|^\s*Error\s*\(',native,re.I|re.M)))
    hf=[cp.number(x) for x in re.findall(r'^E\(0\)\s*\.{2,}\s*('+F+r')',native,re.M)]
    if not hf: hf=[cp.number(x) for x in re.findall(r'^Reference energy\s*\.{2,}\s*('+F+r')',native,re.M)]
    checks['finite_RHF_energy']=bool(hf) and all(math.isfinite(x) and abs(x-hf[0])<=1e-9 for x in hf)
    q.update(CC_residual_traces=traces,CC_terminal_residuals=ccres,triples_iteration_traces=triples,triples_terminal_deltas=tri_tail,
        triples_printed_energy_tolerance=etol,CCSD_energies_hartree=energies,triples_hartree=correction,RHF_energies_hartree=hf,
        native_local_thresholds=local,aux_C=aux,
        local_diagnostics=[{'line':i+1,'text':v} for i,v in enumerate(native.splitlines()) if re.search(r'PNO|TNO|localiz|pairs|triples|core|memory|converg|WARNING',v,re.I)])
    q['status']='complete' if all(checks.values()) else 'failed_admission';return q

def saved_canonical(pin):
    r=read(pin);b.check(r['status']=='complete' and all(r['checks'].values()),'Canonical control unavailable')
    for k in ('input','output','stderr'):b.verify(r[k])
    q=(cp.parse if r['panel']['sha256']=='00f7ab059686b5b8c475c99d7993f75d70bea8075c1b7baf9aa4f8cdc1c894dd' else donor.parse)(b.verify(r['output']).read_text(),r['cell'],r['returncode'])
    b.check(q['status']=='complete' and q['energy_hartree']==r['energy_hartree'],'Canonical raw replay differs')
    native=b.verify(r['output']).read_text().split('****END OF INPUT****')[-1]
    hf=[cp.number(x) for x in re.findall(r'^E\(0\)\s*\.{2,}\s*('+F+r')',native,re.M)]
    b.check(bool(hf) and all(math.isfinite(x) and abs(x-hf[0])<=1e-9 for x in hf),'Canonical RHF energy unavailable/inconsistent')
    return dict(result=pin,output=r['output'],input=r['input'],cell=r['cell'],energy_hartree=r['energy_hartree'],RHF_energy_hartree=hf[0],elapsed_seconds=r['elapsed_seconds'])

def definitions(controls,old):
    cells=[]
    for key,saved in controls.items():
        c=copy.deepcopy(saved['cell']);c.pop('template_input',None);c.update(cell_id='control__'+key,method='dlpno',tier='tight',stage='controls',canonical_id=key)
        if 'donor' not in c:c['donor']='water'
        cells.append(c)
    for tier in TIERS:
        for gid in h.GIDS:
            g=old['geometries'][gid];n=g['water_count'];m=g['metal']
            cells.append(dict(cell_id=f'cluster__{tier}__{gid}',geometry_id=gid,method='dlpno',tier=tier,stage='cluster',metal=m,role='complex',d4=False,
                symbols=g['symbols'],positions_A=g['positions_A'],source_geometry=g['xyz'],ghost_mask=[False]*len(g['symbols']),charge=3,multiplicity=1,
                physical_multiplicity_of_parent=g['physical_multiplicity'],expected_effective_nuclear_charges=[11]+[8,1,1]*n,
                expected_electrons=8+10*n,expected_correlated_electrons=8+10*n,expected_basis_functions=47+43*n,
                expected_basis_functions_by_center=[47]+[31,6,6]*n,expected_ecp_core=46 if m=='La' else 55))
    return cells

def prepare(config,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config).read_text())
    b.check(cfg['implementation']==b.pin(__file__),'Source differs')
    for key,p in HELPERS.items():b.check(cfg[key]==b.pin(p),'Helper differs')
    for key in ('specification','feasibility','evidence'):b.verify(cfg[key])
    old=h.load(b.verify(cfg['coordination_panel']));water=read(cfg['water_collection']);donors=read(cfg['donor_collection'])
    controls={}
    for pin in water['result_pins']+donors['native_result_pins']:
        r=read(pin)
        if r['cell'].get('method','cc')!='cc':continue
        key=('water__'+r['cell_id']) if r['cell_id'].startswith(('La__','Dy__')) else r['cell_id']
        controls[key]=saved_canonical(pin)
    b.check(len(controls)==18,'Exact eighteen canonical controls required')
    out.mkdir(parents=True,exist_ok=False);(out/'inputs').mkdir();(out/'assets').mkdir();assets={}
    for m in ('La','Dy'):
        assets[m]={}
        dst=out/'assets'/f'{m}_basis.inc';shutil.copyfile(b.verify(old['assets'][m]['basis']),dst);assets[m]['basis']=b.pin(dst)
        raw=b.verify(old['assets'][m]['aux']).read_text();b.check(raw.count('NewAuxJGTO')==1,'Unexpected metal auxiliary definition')
        dst=out/'assets'/f'{m}_aux_c.inc';dst.write_text(raw.replace('NewAuxJGTO','NewAuxCGTO'));assets[m]['aux_c']=b.pin(dst)
    p=dict(protocol=PROTOCOL,adapter=b.pin(__file__),configuration=b.pin(config),specification=cfg['specification'],feasibility=cfg['feasibility'],evidence=cfg['evidence'],
        coordination_panel=cfg['coordination_panel'],water_collection=cfg['water_collection'],donor_collection=cfg['donor_collection'],
        orca=old['orca'],native_modules=cfg['native_modules'],assets=assets,source_assets=old['assets'],controls=controls,saved_complexes=old['saved_complexes'],
        geometries=old['geometries'],cells=definitions(controls,old),tiers=TIERS,maximum_native_calls=26,ranks=4,scout_cell_id=SCOUT,
        run_wrapper=b.pin(HERE/'RUN.sbatch'),created_utc=b.now(),**{k:cfg[k] for k in HELPERS})
    for c in p['cells']:
        f=out/'inputs'/(c['cell_id']+'.inp');f.write_text(deck(p,c,4,3000));c['template_input']=b.pin(f)
    pin=b.write(out/'PANEL.json',p);v=validate(pin['path']);v.update(canonical_raw_replays=18,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu)
    b.write(out/'PREFLIGHT.json',v);print(json.dumps(dict(panel=pin,status=v['status'],new_native_calls=0)))

def load(path):
    p=json.loads(Path(path).read_text());b.check(p['protocol']==PROTOCOL and p['adapter']==b.pin(__file__),'Protocol/source differs')
    for key in ('adapter','configuration','specification','feasibility','evidence','coordination_panel','water_collection','donor_collection','orca','run_wrapper',*HELPERS):b.verify(p[key])
    for key,v in HELPERS.items():b.check(p[key]==b.pin(v),'Helper source differs')
    for v in p['native_modules'].values():b.verify(v)
    old=read(p['coordination_panel']);expected=definitions(p['controls'],old)
    b.check(len(expected)==len(p['cells'])==p['maximum_native_calls']==26 and p['tiers']==TIERS and p['ranks']==4 and p['scout_cell_id']==SCOUT,'Finite scope differs')
    for m,a in p['assets'].items():
        b.check(b.verify(a['basis']).read_bytes()==b.verify(p['source_assets'][m]['basis']).read_bytes(),'Orbital/ECP asset differs')
        b.check(b.verify(a['aux_c']).read_text()==b.verify(p['source_assets'][m]['aux']).read_text().replace('NewAuxJGTO','NewAuxCGTO'),'AuxC coefficient change')
    for c,e in zip(p['cells'],expected):
        b.check({k:v for k,v in c.items() if k!='template_input'}==e,'Cell changed');b.verify(c['source_geometry'])
        b.check(b.verify(c['template_input']).read_text()==deck(p,c,4,3000),'Input changed')
    return p

def validate(path):
    p=load(path);rows=[]
    for c in p['cells']:
        raw=b.verify(c['template_input']).read_text();header,*lines=raw.split('* xyz ',1)[1].splitlines();lines=lines[:-1]
        q,m=map(int,header.split());tokens=[x.split() for x in lines];labels=[x[0] for x in tokens];xyz=[[float(z) for z in x[1:]] for x in tokens]
        z=[0 if x.endswith(':') else (11 if i==0 else b.Z[x]) for i,x in enumerate(labels)]
        b.check(q==c['charge'] and m==1 and z==c['expected_effective_nuclear_charges'] and sum(z)-q==c['expected_electrons'],'Electron/state inventory')
        b.check([x.rstrip(':') for x in labels]==c['symbols'] and xyz==c['positions_A'] and [x.endswith(':') for x in labels]==c['ghost_mask'],'Coordinates/ghosts differ')
        b.check(b.xyz_read(c['source_geometry']['path'])==(c['symbols'],c['positions_A']),'Source XYZ differs')
        b.check(('NewECP' in raw)==bool(c['expected_ecp_core']) and ('DelECP '+c['metal'] in raw)==c['ghost_mask'][0],'ECP inventory')
        b.check(not any(x in raw for x in ('D4','RIJCOSX','RIJK','EnGrad','AutoAux')) and 'DLPNO-CCSD(T1)' in raw and 'NoFrozenCore' in raw,'Method differs')
        rows.append({k:c[k] for k in ('cell_id','charge','expected_electrons','expected_basis_functions','expected_ecp_core','tier')})
    return dict(status='passed_static_preflight',panel=b.pin(path),declared_cells=26,independent_deck_inventory=rows,new_native_calls=0,local_threads=1,
        native_DLPNO_execution_and_output_grammar_admitted=False,utc=b.now())

def release_contract(p,path,r):
    b.check(r['status']=='approved_coordination_local_cc_subset' and r['implementation']==p['adapter'] and r['panel']==b.pin(path),'Root release differs')
    ids=[c['cell_id'] for c in p['cells']];b.check(r['cell_ids']==ids and r['maximum_native_calls']==26 and r['ranks']==4,'Released denominator differs')
    pre=read(r['preparation']);b.verify(r['independent_review']);b.check(pre['status']=='passed_static_preflight' and pre['panel']==b.pin(path),'Preparation not admitted')
    stage=r['stage'];allowed=[c['cell_id'] for c in p['cells'] if c['stage']=='controls'] if stage=='controls' else ([SCOUT] if stage=='cluster_scout' else [x for x in ids if x.startswith('cluster__') and x!=SCOUT])
    b.check(stage in ('controls','cluster_scout','cluster_remaining') and r['released_cell_ids'] and len(set(r['released_cell_ids']))==len(r['released_cell_ids']) and set(r['released_cell_ids'])<=set(allowed),'Released subset differs')
    if stage!='controls':
        qual=read(r['control_qualification']);b.check(qual['status']=='local_controls_passed' and qual['panel']==b.pin(path) and qual['complete']==18 and qual['all_controls_passed'],'Canonical-overlap controls unavailable')
        b.check(qual['implementation']==p['adapter'],'Qualification source differs')
    if stage=='cluster_remaining':
        scout=read(r['cluster_scout_review']);b.check(scout['status']=='passed_cluster_scout_review' and scout['panel']==b.pin(path) and scout['cost_and_memory_accepted'] is True,'Actual scout cost/native admission missing')
        c=next(c for c in p['cells'] if c['cell_id']==SCOUT);admit_result(b.verify(scout['result']),c,p,path)
    return r

def admit_result(path,c,p,panel_path):
    path=Path(path);r=json.loads(path.read_text());b.check(r['status']=='complete' and r['returncode']==0 and all(r['checks'].values()),'Native result not admitted')
    b.check(r['cell_id']==c['cell_id'] and r['cell']==c and r['panel']==b.pin(panel_path) and r['adapter']==p['adapter'],'Cell/source differs')
    st=read(r['start']);rel=read(r['release']);release_contract(p,panel_path,rel)
    b.check(r['start']==b.pin(path.parent/'START.json') and st['cell']==c and st['panel']==r['panel'] and st['release']==r['release'] and st['input']==r['input'],'START lineage differs')
    b.check(c['cell_id'] in rel['released_cell_ids'] and path.resolve()==Path(rel['execution_root']).resolve()/'tasks'/c['cell_id']/'execution/RESULT.json','One-attempt path/subset differs')
    b.check(st['ranks']==4 and st['memory_MiB']>0 and st['maxcore_MiB']==math.floor(.75*st['memory_MiB']/4),'Resources differ')
    b.check(st['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'] and r['native_program_invocations_attempted']==1,'Invocation differs')
    b.check(b.verify(r['input']).read_text()==deck(p,c,4,st['maxcore_MiB']),'Executed input differs');b.verify(r['stderr'])
    q=parse(b.verify(r['output']).read_text(),c,r['returncode']);b.check(q['status']=='complete' and all(r[k]==v for k,v in q.items()),'Raw replay differs')
    return dict(q,result=b.pin(path),output=r['output'],input=r['input'],start=r['start'],release=r['release'],elapsed_seconds=r['elapsed_seconds'],allocated_cpu_seconds=r['allocated_cpu_seconds'])

def run(path,release,cell_id,out):
    p=load(path);r=release_contract(p,path,json.loads(Path(release).read_text()));b.check(cell_id in r['released_cell_ids'],'Cell not released')
    b.check(out.resolve()==Path(r['execution_root']).resolve()/'tasks'/cell_id/'execution' and not out.exists(),'One-attempt path differs/existing')
    c=next(c for c in p['cells'] if c['cell_id']==cell_id);b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native calls require Slurm')
    ranks=int(os.environ.get('SLURM_NTASKS','0'));b.check(ranks==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread ranks required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));b.check(mem>0,'Actual memory grant required');maxcore=math.floor(.75*mem/ranks)
    out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(deck(p,c,4,maxcore));env=os.environ.copy()
    for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[k]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','');cmd=[p['orca']['path'],'endpoint.inp','--bind-to none']
    st=b.write(out/'START.json',dict(cell=c,panel=b.pin(path),release=b.pin(release),input=b.pin(out/'endpoint.inp'),job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),ranks=4,memory_MiB=mem,maxcore_MiB=maxcore,command=cmd,maximum_native_calls=1,utc=b.now()))
    tick=time.monotonic();code=None
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:code=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False).returncode
        q=parse((out/'endpoint.out').read_text(errors='replace'),c,code)
    except Exception as e:q=dict(status='failed_collection',exception=repr(e))
    elapsed=time.monotonic()-tick;q.update(cell_id=cell_id,cell=c,panel=b.pin(path),adapter=p['adapter'],release=b.pin(release),start=st,returncode=code,
        native_program_invocations_attempted=1,elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*4,job_id=os.environ['SLURM_JOB_ID'],finished_utc=b.now(),
        output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'),input=b.pin(out/'endpoint.inp'),artifact_sizes={x.name:x.stat().st_size for x in out.iterdir() if x.is_file()})
    b.write(out/'RESULT.json',q);print(json.dumps(dict(status=q['status'],cell_id=cell_id,seconds=elapsed)));return 0 if q['status']=='complete' else 1

def control_score(p,rows):
    per={};hf={};aux={};checks={}
    for c in p['cells']:
        if c['stage']!='controls':continue
        r=rows[c['cell_id']];old=p['controls'][c['canonical_id']];key=(c['donor'],c['metal']);per.setdefault(key,{})[c['role']]=(r['energy_hartree'],old['energy_hartree'])
        hf[c['cell_id']]=(r['RHF_energies_hartree'][0]-old['RHF_energy_hartree']);checks['RHF__'+c['cell_id']]=abs(hf[c['cell_id']])<=1e-6
        aux.setdefault(key,[]).append(r['aux_C'])
    cpvals={};errors={}
    for (d,m),roles in per.items():
        full=roles['complex'];frag=[v for k,v in roles.items() if k!='complex'];b.check(len(frag)==2,'CP denominator')
        value=[(full[i]-sum(x[i] for x in frag))*b.HA for i in range(2)];key=d+'__'+m;cpvals[key]=value;errors[key]=value[0]-value[1];checks['CP__'+key]=abs(errors[key])<=.25
        checks['AuxC__'+key]=all(v==aux[(d,m)][0] for v in aux[(d,m)])
    paired={d:[cpvals[d+'__Dy'][i]-cpvals[d+'__La'][i] for i in range(2)] for d in ('water','acetate_bi','acetamide')}
    for d,v in paired.items():checks['paired_CP__'+d]=abs(v[0]-v[1])<=.25
    exchange={d:[paired[d][i]-paired['water'][i] for i in range(2)] for d in ('acetate_bi','acetamide')}
    for d,v in exchange.items():checks['paired_exchange__'+d]=abs(v[0]-v[1])<=.25
    return dict(checks=checks,all_controls_passed=all(checks.values()),per_metal_CP_local_canonical_kcal_mol=cpvals,CP_errors_kcal_mol=errors,paired_CP_local_canonical=paired,paired_exchange_local_canonical=exchange,RHF_errors_hartree=hf)

def collect(path,inventory,out,stage):
    tick,cpu=time.monotonic(),time.process_time();p=load(path);inv=json.loads(Path(inventory).read_text());cells=[c for c in p['cells'] if stage=='full' or c['stage']=='controls']
    b.check(set(inv)=={c['cell_id'] for c in cells},'Exact result inventory including unavailable slots required')
    out.mkdir(parents=True,exist_ok=False);rows={c['cell_id']:dict(status='unavailable',error=None) for c in cells};good={}
    for c in cells:
        try:
            b.check(inv[c['cell_id']] is not None,'Missing native result');admitted=admit_result(inv[c['cell_id']],c,p,path);rows[c['cell_id']]=dict(admitted,status='complete');good[c['cell_id']]=admitted
        except Exception as e:rows[c['cell_id']]=dict(status='unavailable',error=repr(e),path=inv[c['cell_id']])
    q=dict(status='unavailable',implementation=b.pin(__file__),panel=b.pin(path),inventory=b.pin(inventory),stage=stage,declared=len(cells),complete=len(good),rows=rows,all_controls_passed=False,new_native_calls=0)
    if all(c['cell_id'] in good for c in p['cells'] if c['stage']=='controls'):
        sc=control_score(p,good);q.update(sc);q['status']='local_controls_passed' if stage=='controls' and sc['all_controls_passed'] else 'local_controls_failed'
    if stage=='full' and len(good)==26:
        def contrast(e):return e['aquo9_central_Dy']-e['aquo9_central_La']-e['aquo8_central_Dy']+e['aquo8_central_La']
        values={t:{g:good[f'cluster__{t}__{g}']['energy_hartree'] for g in h.GIDS} for t in TIERS}
        cc={t:contrast(v)*b.HA for t,v in values.items()};old=p['saved_complexes'];pe=contrast({g:old[g]['electronic_energy_hartree'] for g in h.GIDS})*b.HA;d4=contrast({g:old[g]['D4_energy_hartree'] for g in h.GIDS})*b.HA
        diff=cc['tighter']-cc['tight'];q.update(status='local_coordination_method_sensitivity_complete',raw_CC_C_kcal_mol=cc,PBE0_electronic_C_kcal_mol=pe,original_D4_C_kcal_mol=d4,PBE0_D4_C_kcal_mol=pe+d4,
            tighter_minus_tight_C_kcal_mol=diff,target_threshold_sensitivity_passed=abs(diff)<=.25,
            endpoint_tighter_minus_tight_kcal_mol={g:(values['tighter'][g]-values['tight'][g])*b.HA for g in h.GIDS},
            CC_minus_PBE0_D4_C_kcal_mol={t:v-pe-d4 for t,v in cc.items()},method_sensitivity_controls_passed=q['all_controls_passed'] and abs(diff)<=.25)
    q.update(physical_accuracy_qualified=False,cluster_CP_or_B_qualified=False,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu)
    b.write(out/'RESULT.json',q);print(json.dumps({k:q[k] for k in ('status','complete','declared','new_native_calls')}));return q

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--config',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('validate');a.add_argument('--panel',required=True);a.add_argument('--out',required=True)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--release',required=True);a.add_argument('--cell',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--inventory',required=True);a.add_argument('--stage',choices=['controls','full'],required=True);a.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.action=='prepare':prepare(a.config,a.out)
    elif a.action=='validate':b.write(a.out,validate(a.panel))
    elif a.action=='run':return run(a.panel,a.release,a.cell,a.out)
    else:collect(a.panel,a.inventory,a.out,a.stage)
if __name__=='__main__':sys.exit(main())
