#!/usr/bin/env python3
"""Ten fixed whole256 CPCM-water singlepoint energies; every stage requires Slurm."""
import argparse
import datetime as dt
from decimal import Decimal, localcontext
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import socket
import sys
import time
import traceback

PROTOCOL='lbt_solvent_response_v1'
IDS=tuple(s+'_'+t for s in ('7CCO','7CCN') for t in ('before_La','after_La','before_Dy','after_Dy'))+('7CCN_step92_La','7CCN_step725_La')
INVENTORY_SHA='08e9b43b8e32a2b7b1919ce5d1fed1dc5f29468f4e112b6e92edcf15fecbfe32'
HA=627.509474
CAVEATS={
 'scope':'Ten fixed whole256 energy-only CPCM-water endpoints; zero new vacuum/ALPB/MACE/GFN/MM calls or optimizer steps.',
 'solvent':'CPCM changes self-consistent density, cavity and solvent treatment relative to fixed-FF ALPB; differences do not uniquely isolate polarization or establish accuracy.',
 'state':'Conditional termini/protonation and charge -2; La physical/effective singlet, Dy physical6/effective1 frozen-f ECP55; no SOC or physical-state qualification.',
 'interpretation':'Conditional deformation works and matched metal contrasts, not binding affinities, equilibrium, force qualification or production promotion.',
 'failure':'No automatic retry, changed SCF controls or expanded selection. Missing native evidence remains unavailable.'}


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve()
    with p.open('rb') as h: sha=hashlib.file_digest(h,'sha256').hexdigest()
    return {'path':str(p),'sha256':sha}


def verify(record):
    require(pin(record['path'])['sha256']==record['sha256'],'Changed pin: '+record['path'])
    return Path(record['path'])


def read(record): return json.loads(verify(record).read_text())


def write(path,value):
    with Path(path).open('x') as h: json.dump(value,h,indent=2,allow_nan=False); h.write('\n')
    return pin(path)


def module(name,record):
    spec=importlib.util.spec_from_file_location(name,verify(record)); obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj); return obj


def binding(cfg):
    return {k:(v['sha256'] if k=='implementation' else v) for k,v in cfg.items()
            if k not in ('preparation_result','execution_decision','endpoint_results','_config_pin')}


def printed_energy(text):
    values=re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+0-9.Ee]+)',text)
    require(len(values)==1 and math.isfinite(float(values[0])),'One finite final printed energy required')
    return values[0]


def inputs(cfg):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    require(cfg['protocol']==PROTOCOL and cfg['implementation']['sha256']==pin(__file__)['sha256']
            and cfg['inventory']['sha256']==INVENTORY_SHA and cfg['mpi_ranks']==12
            and cfg['permitted_tasks']==list(IDS) and cfg['maximum_new_DFT_calls']==10
            and cfg['reporting_resolution_kcal_mol']==1.0,'Frozen ten-cell protocol differs')
    for value in cfg.values():
        if isinstance(value,dict) and set(value)=={'path','sha256'}: verify(value)
    require(pin(sys.executable)['sha256']==cfg['python']['sha256'] and sys.prefix==cfg['python_prefix'],
            'Original reference Python required')
    inv=read(cfg['inventory']); require(inv['protocol']=='lbt_solvent_response_inventory_v1'
            and [r['cell_id'] for r in inv['cells']]==list(IDS),'Ten selected cells differ')
    source=read(inv['source_inventory']); require(inv['source_inventory']['sha256']=='e2c196d1b5961cc74d50647e0565780928c3c1c103b282d563833f4b7e003957','Original eight-state authority differs')
    old=module('solvent_whole_reference_primitives',inv['reference_helper']); native=inv['native_reference']
    compact=read(native['compact_assets_manifest'])
    require(native['compact_assets_manifest']['sha256']==old.COMPACT_SHA and native['method']==old.METHOD
            and native['metal_assets']==compact['assets'] and native['orca']==compact['orca'],'Original native recipe differs')
    for v in compact['implementation'].values(): verify(v)
    for assets in compact['assets'].values():
        for v in assets.values(): verify(v)
    verify(compact['orca']); sys.path.insert(0,str(verify(compact['implementation']['metal_environment_frozen_f.py']).parent))
    names=('metal_environment_frozen_f','run_orca_task_manifest','ggr_sensitivity','mace_omol_vacuum',
           'render_orca_runtime_input','affordable_response','affordable_common','metal_environment_electronic_state')
    mods={n:importlib.import_module(n) for n in names}
    for n,obj in mods.items(): require(pin(obj.__file__)['sha256']==compact['implementation'][n+'.py']['sha256'],'Imported helper differs: '+n)
    runtime=module('solvent_native_runtime',native['runtime_reference_helper']); model=mods['metal_environment_frozen_f']
    require(model.METHOD==old.METHOD and model.HA_TO_KCAL==HA,'Original method/conversion differs')
    import numpy as np
    collectors={}
    for key,p in inv['parents'].items():
        c=read(p['collector']); cc=read(p['config']); verify(p['implementation'])
        require(c['status']=='complete' and c['config']==p['config'] and c['implementation']==p['implementation']
                and cc['implementation']['sha256']==p['implementation']['sha256']
                and c['physical_state_qualified'] is False and c['force_consistency_qualified'] is False,
                'Inherited reference collector admission differs')
        require(c['full_two_available'] if key=='late_7CCN_La' else c['full_four_available'],'Original collector denominator incomplete')
        collectors[key]=c
    require(inv['parents']['late_7CCN_La']['collector']['sha256']=='407bdf88e42c43cb2ae48ed1493cbbc008f5c0eb957628787a6aa1c1b795bb0d','Exact late endpoint pair differs')
    cells={}; evidence={}
    for r in inv['cells']:
        tid=r['cell_id']; metal=r['metal']; c=collectors[r['parent']]; key=r['collector_cell_id']
        nr=c['rows'][key]; ref=r['reference']; t=ref['task']; atoms=read(r['physical']); ids=[a['id'] for a in atoms]
        x=np.asarray([a['xyz_A'] for a in atoms],dtype=float)
        require(len(ids)==len(set(ids))==256 and ids==r['physical_ids'] and ids[-1]=='metal'
                and x.shape==(256,3) and np.isfinite(x).all(),'Whole physical inventory differs')
        pp=verify(r['positions']); px=np.load(pp,allow_pickle=False) if pp.suffix=='.npy' else json.loads(pp.read_text())['positions_A']
        require(np.array_equal(x,px),'Physical/position coordinates differ')
        symbols=[metal if a['element']=='M' else a['element'] for a in atoms]; st=old.state(symbols,metal,model.NUMBERS)
        require(nr['status']=='available' and nr['task']==t and nr['receipt']==ref['native_receipt']
                and nr['output']==ref['output'] and nr['engrad']==ref['engrad'] and nr['energy_hartree']==ref['energy_hartree']
                and nr['electronic_state']==ref['state']==st and nr['physical_multiplicity']==st['physical_multiplicity']
                and t['physical_parent']==r['physical'] and t['positions']==r['positions'] and t['physical_ids']==ids
                and t['charge']==-2 and t['multiplicity']==1 and t['physical_multiplicity']==st['physical_multiplicity']
                and t['source_id']==r['source_id'] and t['metal']==metal and t['atom_count']==256 and t['target_index']==255
                and nr['physical_state_qualified'] is False and nr['force_consistency_qualified'] is False,
                'Inherited native state/geometry/admission differs: '+tid)
        if tid in IDS[:8]:
            sr=next(v for v in source['cells'] if v['cell_id']==tid)
            require(all(sr[k]==r[k] for k in ('source_id','geometry','metal','physical','positions','physical_ids','ALPB'))
                    and sr['reference']['collector_row']==nr,'Original eight-cell inventory join differs')
        else:
            require(r['geometry']==('step92' if key=='old_La' else 'step725') and metal=='La','Late endpoint remapping differs')
        xyz=model.xyz(verify(t['xyz']))
        require([a[0] for a in xyz]==symbols and np.array_equal([a[1:] for a in xyz],x),'Exact original native XYZ differs')
        text=verify(t['input']).read_text(); template=old.asset_checks(compact,metal)
        require(text==template and text.count(' EnGrad\n')==1,'Original complete scientific deck differs')
        deck=text.replace(' EnGrad\n',' CPCM(water)\n',1)
        require('EnGrad' not in deck and deck.count('CPCM(water)')==1,'Only declared input changes permitted')
        er=read(ref['native_receipt']); verify(ref['engrad']); vacuum=printed_energy(verify(ref['output']).read_text())
        require(float(vacuum)==nr['energy_hartree'] and er['normal_termination'] and er['scf_converged'] and er['returncode']==0
                and er['orca_executable']==compact['orca'] and er['artifacts']['template_input']==t['input']
                and er['artifacts']['xyz']==t['xyz'] and er['artifacts']['output']==ref['output']
                and er['artifacts']['engrad']==ref['engrad'],'Original vacuum native receipt/energy differs')
        sol=c['solvent_admission'][key]; ar=read(r['ALPB']['receipt'])
        require(sol['status']=='available' and sol['receipt']==r['ALPB']['receipt']
                and sol['energy_hartree']==r['ALPB']['energy_hartree']==ar['total_hartree']
                and sol['components_hartree']==ar['components_hartree'] and ar['artifacts']==r['ALPB']['artifacts']
                and r['ALPB']['admission_authority']==inv['parents'][r['parent']]['collector'],'Saved ALPB authority differs')
        for v in ar['artifacts'].values(): verify(v)
        start=read(ar['artifacts']['START.json']); ex=read(ar['artifacts']['EXECUTION.json'])
        require(read(start['physical_source'])==atoms and start['xyz']==ar['artifacts']['physical.xyz']
                and start['charges']==ar['artifacts']['charges'] and start['command']==ex['command'] and ex['returncode']==0,
                'Saved ALPB physical/execution join differs')
        ax=model.xyz(verify(start['xyz']))
        require([a[0] for a in ax]==symbols and np.array_equal([a[1:] for a in ax],x)
                and math.isfinite(ar['total_hartree']) and abs(ar['charge_e']+2)<=1e-5,'Saved ALPB coordinates/state differ')
        cells[tid]=dict(row=r,x=x,state=st,deck=deck)
        evidence[tid]=dict(source_id=r['source_id'],geometry=r['geometry'],metal=metal,physical=r['physical'],positions=r['positions'],
            physical_ids=ids,state=st,original_vacuum_task=t,original_vacuum_receipt=ref['native_receipt'],original_vacuum_output=ref['output'],
            vacuum_printed_energy_hartree=vacuum,ALPB=r['ALPB']['receipt'],ALPB_energy_hartree=ar['total_hartree'],
            ALPB_parts_hartree=ar['components_hartree'],vacuum_SCF_evidence=nr['SCF_evidence'],
            native_admission_authority=inv['parents'][r['parent']]['collector'],inherited_force_consistency_qualified=False)
    for s in ('7CCO','7CCN'):
        for g in ('before','after'):
            a,b=(cells[s+'_'+g+'_'+m] for m in ('La','Dy'))
            require(a['row']['physical_ids']==b['row']['physical_ids'] and np.array_equal(a['x'],b['x']),'Matched-metal coordinates differ')
    return dict(cfg=cfg,inv=inv,old=old,compact=compact,modules=mods,runtime=runtime,cells=cells,evidence=evidence,
                method=old.METHOD.replace(' EnGrad',' CPCM(water)'))


def task(w,tid,xyz,deck):
    c=w['cells'][tid]; r=c['row']; directory=verify(deck).parent
    return dict(task_id=tid,source_id=r['source_id'],geometry=r['geometry'],metal=r['metal'],xyz=xyz,input=deck,
        output_path=str(directory/'endpoint.out'),engrad_path=None,task_type='single_point',charge=-2,multiplicity=1,
        physical_multiplicity=c['state']['physical_multiplicity'],electronic_state=c['state'],physical_parent=r['physical'],
        positions=r['positions'],physical_ids=r['physical_ids'],atom_count=256,target_index=255,
        original_vacuum_task=r['reference']['task'],original_vacuum_receipt=r['reference']['native_receipt'],
        cached_ALPB=r['ALPB']['receipt'],caveats=CAVEATS)


def make_manifest(w,out,tids):
    tasks=[]; compact=w['compact']; cfg=w['cfg']
    for tid in tids:
        directory=out/tid; directory.mkdir(); c=w['cells'][tid]
        xyz=directory/'physical.xyz'; xyz.write_bytes(verify(c['row']['reference']['task']['xyz']).read_bytes())
        deck=directory/'endpoint.inp'; deck.write_text(c['deck'])
        t=task(w,tid,pin(xyz),pin(deck)); write(directory/'TASK.json',t); tasks.append(t)
    m=dict(protocol_id=PROTOCOL,agreement=cfg['specification'],adapter=pin(__file__),inputs_config=cfg['_config_pin'],
        scientific_binding=binding(cfg),method=w['method'],assets=compact['assets'],orca=compact['orca'],
        implementation=compact['implementation'],execution_policy=compact['execution_policy'],
        execution_resources={'mpi_ranks':12,'concurrent_tasks':1},tasks=tasks,declared_tasks=list(IDS),caveats=CAVEATS)
    rec=write(out/'MANIFEST.json',m); validate_manifest(w,rec,tids); return rec


def validate_manifest(w,record,tids):
    m=read(record); cfg=w['cfg']; compact=w['compact']
    require((list(tids)==list(IDS) or len(tids)==1 and tids[0] in IDS) and [t['task_id'] for t in m['tasks']]==list(tids)
            and m['protocol_id']==PROTOCOL and m['agreement']==cfg['specification']
            and m['adapter']['sha256']==cfg['implementation']['sha256'] and m['scientific_binding']==binding(cfg)==binding(read(m['inputs_config']))
            and m['method']==w['method'] and m['declared_tasks']==list(IDS) and m['caveats']==CAVEATS
            and m['execution_resources']=={'mpi_ranks':12,'concurrent_tasks':1}
            and all(m[k]==compact[k] for k in ('assets','orca','implementation','execution_policy')),'Frozen native manifest differs')
    verify(m['adapter'])
    for t in m['tasks']:
        c=w['cells'][t['task_id']]
        require(t==task(w,t['task_id'],t['xyz'],t['input']) and verify(t['input']).read_text()==c['deck']
                and verify(t['input']).parent==verify(t['xyz']).parent
                and t['xyz']['sha256']==c['row']['reference']['task']['xyz']['sha256'],'Exact frozen input differs')
    return m


def signature(t):
    return {k:(v['sha256'] if k in ('xyz','input') else v) for k,v in t.items() if k!='output_path'}


def prepared(w):
    cfg=w['cfg']; p=read(cfg['preparation_result']); pcfg=read(p['config']); verify(p['implementation'])
    require(p['protocol']==PROTOCOL and p['stage']=='prepare' and p['status']=='prepared'
            and p['implementation']['sha256']==cfg['implementation']['sha256'] and p['scientific_binding']==binding(cfg)==binding(pcfg)
            and read(p['input_evidence'])==w['evidence'] and p['selected_tasks']==list(IDS)
            and p['new_DFT_calls_attempted']==p['native_runner_invocations']==0,'Actual unchanged ten-cell zero-call preparation required')
    validate_manifest(w,p['manifest'],IDS); return p


def released(w,tid):
    p=prepared(w); cfg=w['cfg']; d=read(cfg['execution_decision'])
    require(d['protocol']==PROTOCOL and d['status']=='approved_ten_fixed_CPCM_singlepoints'
            and d['scientific_binding']==binding(cfg) and d['preparation_result']==cfg['preparation_result']
            and d['permitted_execution_tasks']==list(IDS) and d['mpi_ranks']==12 and d['maximum_new_DFT_calls']==10
            and tid in IDS,'Root ten-singlepoint release required')
    return p


def admit(w,record,tid):
    manifest=validate_manifest(w,record,[tid]); mods=w['modules']; t=manifest['tasks'][0]
    _,rows=mods['ggr_sensitivity'].executed(verify(record)); row=rows[tid]; row.update(task=t,task_id=tid)
    try:
        if Path(t['output_path']).exists(): row['SCF_evidence']=w['runtime'].observations(t['output_path'],mods['mace_omol_vacuum'])
        if row['status']!='complete': return w['old'].canonical(row)
        text=verify(row['output']).read_text(); st=t['electronic_state']
        for pattern,wanted in [(r'Program Version\s+(6\.1\.1)\b','6.1.1'),
            (r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)',str(st['charge'])),
            (r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)','1'),
            (r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)',str(st['explicit_electrons'])),
            (r'Hartree-Fock type\s+HFTyp\s+\.{2,}\s+(\w+)','RHF')]:
            require(re.findall(pattern,text)==[wanted],'Executed state/version differs: '+pattern)
        ecp=re.findall(r'Type\s+(\w+)\s+ECP(?:\s+(\S+))?\s+\(replacing\s+(\d+)\s+core electrons',text)
        require(len(ecp)==1 and ecp[0][0]==t['metal'] and ecp[0][2]==str(st['ecp_core_electrons'])
                and 'DFTD4' in text and 'INITIAL GUESS: MODEL POTENTIAL' in text
                and not re.search(r'Reading point charge file|gCP correction\s+[-+0-9.]|^\s*(?:SMD|COSMO) SOLVATION MODEL',text,re.M),
                'Executed ECP/dispersion/guess/external-field differs')
        solvent={}
        require(len(re.findall(r'^\s*CPCM SOLVATION MODEL\s*$',text,re.M))==1,'Unique native CPCM model required')
        for key,pattern,wanted in [
            ('solvent',r'^\s*Solvent:\s*\.{2,}\s*(\S+)\s*$','WATER'),
            ('epsilon',r'^\s*Epsilon\s+\.{2,}\s*(\S+)\s*$','80.1510'),
            ('refractive_index',r'^\s*Refrac\s+\.{2,}\s*(\S+)\s*$','1.3328'),
            ('solvent_radius_A',r'^\s*Rsolv\s+\.{2,}\s*(\S+)\s*$','1.3000'),
            ('surface',r'^\s*Surface type\s+\.{2,}\s*([^\n]+)$','GAUSSIAN VDW')]:
            found=[v.strip() for v in re.findall(pattern,text,re.M)]; require(found==[wanted],'Native CPCM setting differs: '+key); solvent[key]=wanted
        radii=re.findall(r'Radius for\s+(\w+)\s+used is\s+([-+0-9.]+)\s+Bohr\s+\(=\s+([-+0-9.]+)\s+Ang\.\)',text)
        require(len(radii)==5 and {r[0] for r in radii}=={'C','H','N','O',t['metal']}
                and all(math.isfinite(float(v)) and float(v)>0 for r in radii for v in r[1:]),'Native element-dependent CPCM radii unavailable')
        solvent['element_radii']=radii
        solvent['native_cavity_and_energy_lines']=[line for line in text.splitlines() if re.search(r'Cavity|CPCM|Solvation|Dielectric|Gpol|Gcds',line,re.I)]
        row['CPCM_evidence']=solvent
        require(row['SCF_evidence']['ConvCheckMode']=='Total+1el-Energy' and row['SCF_evidence']['ConvForced']==0,'Original SCF convergence semantics differ')
        energy=printed_energy(text); require(float(energy)==row['energy_hartree'],'Final printed native energy differs')
        receipt=read(row['receipt']); sidecar=read(receipt['artifacts']['runtime_sidecar']); policy=sidecar['memory_policy']
        require(receipt['orca_executable']==w['compact']['orca'] and receipt['returncode']==0 and receipt['normal_termination']
                and receipt['scf_converged'] and 'engrad' not in receipt['artifacts'] and receipt['parallelism']['nprocs']==12,
                'Native energy-only execution failed')
        require(policy['workers']==1 and policy['ranks_per_worker']==12 and policy['memory_source']=='SLURM_MEM_PER_NODE',
                'Original allocation-derived memory policy differs')
        row.update(status='available',printed_energy_hartree=energy,executed_ecp=ecp,electronic_state=st,
            physical_multiplicity=st['physical_multiplicity'],runtime_memory_policy=policy,
            molecular_wall_seconds=(dt.datetime.fromisoformat(receipt['finished_at_utc'])-dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds(),
            physical_state_qualified=False,force_consistency_qualified=False,caveats=CAVEATS)
    except Exception as exc: row.update(status='unavailable',admission_failure=repr(exc),energy_hartree=None)
    return w['old'].canonical(row)


def saved_run(w,record,tid):
    p=read(record); cfg=read(p['config']); verify(p['implementation'])
    require(p['protocol']==PROTOCOL and p['stage']=='run-one' and p['requested_task']==tid
            and p['implementation']['sha256']==w['cfg']['implementation']['sha256'] and binding(cfg)==binding(w['cfg'])
            and cfg['preparation_result']==w['cfg']['preparation_result'] and cfg['execution_decision']==w['cfg']['execution_decision'],
            'Saved singleton identity differs')
    released(dict(w,cfg=cfg),tid)
    if p['status']!='complete':
        require(p['selected_tasks'] in ([],[tid]) and p['new_DFT_calls_attempted'] in (0,1) and p['native_runner_invocations'] in (0,1)
                and (p['selected_tasks']==[tid] or p['new_DFT_calls_attempted']==p['native_runner_invocations']==0),'Failed singleton ledger differs')
        return p,dict(status='unavailable',reason='Failed singleton retained',parent_result=record,native_evidence=p['rows'][tid],failure=p.get('failure'),runner_error=p.get('runner_error'))
    require(p['selected_tasks']==[tid] and p['scientific_binding']==binding(cfg) and p['preparation_result']==cfg['preparation_result']
            and p['execution_decision']==cfg['execution_decision'],'Complete singleton release differs')
    row=admit(w,p['manifest'],tid)
    require(row['status']=='available' and row==p['rows'][tid] and p['new_DFT_calls_attempted']==p['native_runner_invocations']==1,'Saved singlepoint replay differs')
    return p,row


def comparison(w,rows):
    result=dict(status='unavailable',units='kcal/mol',cell_corrections={},deformation_works={},matched_Dy_minus_La={},Dy_minus_La_work_difference={},
                late_La_92_to_725={'status':'unavailable'},reporting_resolution_kcal_mol=1.0,caveats=CAVEATS)
    with localcontext() as ctx:
        ctx.prec=50; factor=Decimal(str(HA)); values={}
        for tid in IDS:
            if rows[tid]['status']!='available': result['cell_corrections'][tid]={'status':'unavailable'}; continue
            e=w['evidence'][tid]; vac=Decimal(e['vacuum_printed_energy_hartree']); alpb=Decimal(str(e['ALPB_energy_hartree'])); cp=Decimal(rows[tid]['printed_energy_hartree'])
            values[tid]=dict(vacuum=vac,ALPB=alpb,CPCM=cp,vacuum_plus_ALPB=vac+alpb,CPCM_minus_vacuum=cp-vac,CPCM_minus_vacuum_minus_ALPB=cp-vac-alpb)
            result['cell_corrections'][tid]=dict(status='available',printed_vacuum_energy_hartree=str(vac),printed_CPCM_energy_hartree=str(cp),
                ALPB_energy_hartree=str(alpb),CPCM_minus_vacuum_kcal_mol=str((cp-vac)*factor),
                ALPB_kcal_mol=str(alpb*factor),CPCM_minus_vacuum_minus_ALPB_kcal_mol=str((cp-vac-alpb)*factor))
        def contrast(weights):
            if any(k not in values for k in weights): return dict(status='unavailable',required_cells=list(weights))
            ds={name:sum((Decimal(v)*values[k][name] for k,v in weights.items()),Decimal(0))*factor for name in next(iter(values.values()))}
            return dict(status='available',weights=weights,kcal_mol={k:float(v) for k,v in ds.items()},decimal_kcal_mol={k:str(v) for k,v in ds.items()})
        for s in ('7CCO','7CCN'):
            for m in ('La','Dy'):
                result['deformation_works'][s+'_'+m]=contrast({s+'_after_'+m:1,s+'_before_'+m:-1})
            for g in ('before','after'):
                result['matched_Dy_minus_La'][s+'_'+g]=contrast({s+'_'+g+'_Dy':1,s+'_'+g+'_La':-1})
            result['Dy_minus_La_work_difference'][s]=contrast({s+'_after_Dy':1,s+'_before_Dy':-1,s+'_after_La':-1,s+'_before_La':1})
        result['late_La_92_to_725']=contrast({'7CCN_step725_La':1,'7CCN_step92_La':-1})
    result['status']='complete' if len(values)==len(IDS) else 'incomplete'; return result


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--stage',choices=('prepare','run-one','collect'),required=True)
    p.add_argument('--task',choices=IDS); p.add_argument('--config',required=True); p.add_argument('--output',required=True); args=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=False); start=time.monotonic()
    result=dict(protocol=PROTOCOL,stage=args.stage,status='unavailable',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),declared_tasks=list(IDS),selected_tasks=[],requested_task=args.task,
        rows={t:dict(status='unavailable',reason='not_attempted') for t in IDS},comparison={'status':'unavailable','required_cells':list(IDS)},
        new_DFT_calls_attempted=0,native_runner_invocations=0,new_vacuum_calls=0,new_MACE_calls=0,new_ALPB_calls=0,new_GFN_calls=0,new_MM_calls=0,
        optimizer_steps=0,new_geometry_count=0,caveats=CAVEATS)
    try:
        cfg=read(result['config']); cfg['_config_pin']=result['config']; w=inputs(cfg)
        result.update(scientific_binding=binding(cfg),input_evidence=write(out/'INPUT_EVIDENCE.json',w['evidence']))
        if args.stage=='prepare':
            require(args.task is None and cfg['preparation_result'] is None and cfg['execution_decision'] is None
                    and cfg['endpoint_results']==dict.fromkeys(IDS),'Preparation cannot bind future outcomes')
            result.update(manifest=make_manifest(w,out,IDS),selected_tasks=list(IDS),status='prepared')
        elif args.stage=='run-one':
            require(args.task in IDS and cfg['endpoint_results']==dict.fromkeys(IDS),'One fresh singleton required')
            tid=args.task; prep=released(w,tid); result.update(preparation_result=cfg['preparation_result'],execution_decision=cfg['execution_decision'])
            result['allocation']=w['runtime'].allocation(w['modules']['render_orca_runtime_input'],out)
            record=make_manifest(w,out,[tid]); t=read(record)['tasks'][0]; prior=next(t for t in read(prep['manifest'])['tasks'] if t['task_id']==tid)
            require(signature(t)==signature(prior),'Singleton differs from prepared task')
            result.update(manifest=record,selected_tasks=[tid]); write(out/'START.json',result)
            try:
                result['native_runner_invocations']+=1
                w['modules']['run_orca_task_manifest'].run_manifest(verify(record),orca_path=verify(w['compact']['orca']),
                    workers=1,nprocs=12,selected_tasks={tid},expected_runner_sha256=w['compact']['execution_policy']['task_runner']['sha256'])
            except Exception as exc: result['runner_error']=repr(exc)
            finally: result['new_DFT_calls_attempted']=int(Path(t['output_path']).exists())
            result['rows'][tid]=admit(w,record,tid); result['status']='complete' if result['rows'][tid]['status']=='available' else 'unavailable'
        else:
            require(args.task is None and set(cfg['endpoint_results'])==set(IDS),'Ten-slot saved collection required')
            prepared(w); result.update(preparation_result=cfg['preparation_result'],execution_decision=cfg['execution_decision'],parent_costs={})
            for tid in IDS:
                if cfg['endpoint_results'][tid] is None: continue
                try:
                    parent,row=saved_run(w,cfg['endpoint_results'][tid],tid); result['rows'][tid]=row
                    result['parent_costs'][tid]=dict(result=cfg['endpoint_results'][tid],slurm_job_id=parent['slurm_job_id'],
                        elapsed_seconds=parent['elapsed_seconds'],application_allocated_CPU_seconds=parent['application_allocated_CPU_seconds'])
                except Exception as exc: result['rows'][tid]=dict(status='unavailable',reason=repr(exc),parent_result=cfg['endpoint_results'][tid])
            result['comparison']=comparison(w,result['rows']); result['status']=result['comparison']['status']
        verify(result['config'])
    except Exception as exc: result.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-start,peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        physical_state_qualified=False,force_accuracy_qualified=False,no_automatic_retry=True)
    result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ.get('SLURM_CPUS_ON_NODE','0'))
    result['scheduler_costs']='Actual scheduler accounting is authoritative; allocated capacity is not utilization; parent RSS is not MPI peak.'
    write(out/'RESULT.json',result); return 0 if result['status'] in ('prepared','complete') else 1


if __name__=='__main__': raise SystemExit(main())
