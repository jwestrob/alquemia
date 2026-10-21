#!/usr/bin/env python3
"""One fixed central PModel memory recovery; numerical work is Slurm-only."""
import argparse
import ast
import datetime as dt
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import socket
import subprocess
import sys
import time
import traceback

PROTOCOL='lbt_omol_memory_recovery_v1'
TID='central_La'
RANKS=12
ALLOCATED_CPUS=24
CPUS_PER_TASK=2
MEMORY_MiB=96000
MAXCORE_MB=6200
METHOD='wB97M-V/def2-TZVPD UKS SCNL, La ECP46; PModel derived-recipe variant on ORCA6.1.1'
HA,BOHR=627.509474,0.529177210903
FIXED={
 'direction_helper':'cd4396c02cf50905307acbf90772d062ae4a6a12aaf9914e01562167dca1ab32',
 'direction_config':'dd4875fabd6d17a12d8c73ad7254068ba83a7070b33d079c62b3b7c2dcd74ff7',
 'external_basis':'7389c56f0cff9e76098b0757a34d1abdd79c48e9249597cea074a0106a3462a6',
 'primary_calc':'6cbf04e71579e9021ba31628f04463a59b28bc063834dc395581f579133c9843',
 'primary_recipe':'15d2e4d2c3e204998f29a07433aa56ad942d273d88e95998b4e1258d18a9d842',
 'failed_result':'881852372e543e2f2c721e1b162997cb348c2bc3ddbbbb7d6e910560f7b5c55f',
 'failed_input':'d482826aa2384c0e00f8c8fb0cfa5dc965e22a3c2a99f17a4e608602a45200c6',
 'failed_config':'1c9f22df1c9ac096568dccac7689a5e6c9f8af4d8f8ae2bd0a3d589aa9ca0c83',
 'failed_implementation':'aa62db4f37de3c50765ad81115dd75cb4aef2a99d92b4465d9bbe59df9a4a74f',
 'oom_result':'0fc5e6a1ddd98de6bfeb233d0f385f315e6f199a9e608e31009b9d3371119d3c',
 'oom_input':'4359537c4a569f2a68793f60a509bda509c854a15a7aaf1060790a9407086896',
 'oom_config':'b654d846a9e000cb7adb0c66718b38df7ab64451f5b242283ff843c176df5a0d',
 'oom_implementation':'27ac4e80b5fb53c049adc6f67822f1629cd3bd7cf372edf1d4708cedcc829856',
 'oom_xyz':'55e07d8146adc13c09203243849b65e80a6afef89b181cfe4c22b265a994b72d',
}
NUMBERS={'H':1,'C':6,'N':7,'O':8,'La':57}
CAVEATS={
 'method':'PModel initial-guess variant of published OMol25 recipe on ORCA6.1.1 with native SCNL; not exact historical6.0.0 reproduction. The SCF root may differ.',
 'state':'Conditional vacuum dianion and UKS singlet; rotation is an initial guess, not physical-state/stability/spin-purity qualification.',
 'force':'Central vacuum comparison only; no independent derivative or force-accuracy qualification.',
 'scope':'One requested central EnGrad/runner invocation; native PModel atomic-density subcalculations are allowed and included in cost. No new MACE, solvent, GFN, MM, optimizer, displaced or other-metal endpoints.',
 'failure':'Unknown native evidence remains pending saved-output review; no molecular retry for parser repair.',
 'resources':'12 MPI processes with one numerical thread each; 24 allocated/charged CPUs provide default96000MiB. Allocated capacity is not measured utilization.',
 'interpretation':'No cross-Hamiltonian absolute-energy work, unique energetic cause, equilibrium, affinity or preference claim.'}


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve()
    with p.open('rb') as h: sha=hashlib.file_digest(h,'sha256').hexdigest()
    return dict(path=str(p),sha256=sha)


def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'],'Changed pin: '+p['path']); return Path(p['path'])


def read(p): return json.loads(verify(p).read_text())


def write(path,value):
    with Path(path).open('x') as h: json.dump(value,h,indent=2,allow_nan=False); h.write('\n')
    return pin(path)


def module(p):
    spec=importlib.util.spec_from_file_location('target_saved_direction',verify(p))
    obj=importlib.util.module_from_spec(spec); spec.loader.exec_module(obj); return obj


def binding(cfg):
    return {k:(v['sha256'] if k=='implementation' else v) for k,v in cfg.items()
            if k not in ('preparation_result','execution_decision','_config_pin')}


def recipe(cfg):
    # Parse constants only; importing the upstream calculator would import ASE/Sella.
    wanted={'ORCA_FUNCTIONAL','ORCA_BASIS','ORCA_SIMPLE_INPUT','ORCA_BLOCKS'}; values={}
    for node in ast.parse(verify(cfg['primary_calc']).read_text()).body:
        if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id in wanted:
            values[node.targets[0].id]=ast.literal_eval(node.value)
    require(set(values)==wanted and values['ORCA_FUNCTIONAL']=='wB97M-V' and values['ORCA_BASIS']=='def2-TZVPD'
            and values['ORCA_SIMPLE_INPUT']==['EnGrad','RIJCOSX','def2/J','NoUseSym','DIIS','NOSOSCF','NormalConv','DEFGRID3','ALLPOP'],
            'Published primary constants differ')
    simple=[values['ORCA_FUNCTIONAL'],values['ORCA_BASIS']]+values['ORCA_SIMPLE_INPUT']+['UKS','PModel','NoAutostart','NONBO','NONPA']
    deck='! '+' '.join(simple)+'\n'+'\n'.join(values['ORCA_BLOCKS'])+'\n%scf rotate {514, 515, 20, 1, 1} end end\n* xyzfile -2 1 physical.xyz\n'
    require('%scf Convergence Tight maxiter 300 end' in deck and '%scf THRESH 1e-12 TCUT 1e-13 end' in deck
            and '%basis GTOName "def2-tzvpd.bas" end' in deck and not any(x in simple for x in ('D4','SCNL','PAtom')),
            'Published recipe scope differs')
    failed_deck=verify(cfg['failed_input']).read_text()
    require(failed_deck.count('PAtom')==1 and deck==failed_deck.replace('PAtom','PModel'),
            'Only PAtom-to-PModel input change is permitted')
    require(deck==verify(cfg['oom_input']).read_text(),'Scientific template must exactly equal failed16611')
    return deck


def basis_inventory(path):
    text=path.read_text(); names={'HYDROGEN':'H','CARBON':'C','NITROGEN':'N','OXYGEN':'O','LANTHANUM':'La'}
    counts={}; current=None; indata=False
    for line in text.splitlines():
        s=line.strip()
        if s=='$DATA': indata=True; continue
        if s=='$END': indata=False; current=None
        if not indata: continue
        if s.isalpha() and s.isupper():
            current=names.get(s)
            if current: require(current not in counts,'Duplicate element basis'); counts[current]=0
        hit=re.fullmatch(r'([SPDFGHI])\s+(\d+)',s)
        if current and hit: counts[current]+=2*'SPDFGHI'.index(hit[1])+1
    require(counts=={'H':9,'C':37,'N':37,'O':40,'La':43} and re.search(r'^LA-ECP GEN\s+46\s+3\s*$',text,re.M),
            'Full CHNOLa orbital basis/ECP metadata differs')
    return dict(functions_per_atom=counts,La_core_electrons=46,orbital_file=pin(path),
                auxiliary='built-in def2/J',ECP='external def2-tzvpd.bas La ECP46; actual native coefficients retained separately')


def inputs(cfg):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    require(cfg['protocol']==PROTOCOL and cfg['implementation']['sha256']==pin(__file__)['sha256'] and cfg['mpi_ranks']==RANKS
            and cfg['allocated_cpus']==ALLOCATED_CPUS and cfg['cpus_per_task']==CPUS_PER_TASK,
            'Frozen central scope differs')
    for k,sha in FIXED.items(): require(cfg[k]['sha256']==sha,'Frozen input differs: '+k)
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}: verify(v)
    require(pin(sys.executable)['sha256']==cfg['python']['sha256'] and sys.prefix==cfg['python_prefix'],'Original reference Python required')
    failed=read(cfg['failed_result']); old_cfg=read(cfg['failed_config'])
    require(failed['protocol']=='lbt_omol_target_reference_v1' and failed['stage']=='run' and failed['status']=='invalid'
            and str(failed['slurm_job_id'])=='16602' and failed['native_runner_invocations']==failed['new_DFT_calls_attempted']==1
            and failed['implementation']==cfg['failed_implementation'] and failed['config']==cfg['failed_config']
            and failed['rows'][TID]['task']['input']==cfg['failed_input']
            and all(cfg[k]==old_cfg[k] for k in ('direction_helper','direction_config','external_basis','primary_calc',
                'primary_recipe','python','python_prefix')) and old_cfg['mpi_ranks']==24,'Exact failed parent/settings required')
    oom=read(cfg['oom_result']); oom_cfg=read(cfg['oom_config']); oom_task=oom['rows'][TID]['task']
    require(oom['protocol']=='lbt_omol_pmodel_reference_v1' and oom['stage']=='run' and oom['status']=='invalid'
            and str(oom['slurm_job_id'])=='16611' and oom['native_runner_invocations']==oom['new_DFT_calls_attempted']==1
            and oom['implementation']==cfg['oom_implementation'] and oom['config']==cfg['oom_config']
            and oom_task['input']==cfg['oom_input'] and oom_task['xyz']==cfg['oom_xyz']
            and oom_task['external_basis']['sha256']==cfg['external_basis']['sha256'] and oom_cfg['mpi_ranks']==24
            and all(cfg[k]==oom_cfg[k] for k in ('direction_helper','direction_config','external_basis','primary_calc',
                'primary_recipe','python','python_prefix','failed_result','failed_input','failed_config','failed_implementation')),
            'Only exact16611 resource recovery permitted')
    d=module(cfg['direction_helper']); dc=read(cfg['direction_config']); dc['_config_pin']=cfg['direction_config']
    require(dc['python']==cfg['python'] and dc['python_prefix']==cfg['python_prefix'],'Saved interpreter binding differs')
    saved=d.inputs(dc); c=saved['cell']; np=saved['np']; x=c['x']; state=c['state']; b=basis_inventory(verify(cfg['external_basis']))
    require(len(c['ids'])==len(set(c['ids']))==256 and c['symbols'][-1]=='La' and c['symbols'].count('La')==1
            and set(c['symbols'])==set(NUMBERS) and x.shape==(256,3) and np.isfinite(x).all()
            and sum(NUMBERS[s] for s in c['symbols'])+2==1076
            and all(state[k]==v for k,v in {'charge':-2,'physical_multiplicity':1,'explicit_electrons':1030,
                'alpha_electrons':515,'beta_electrons':515,'ecp_core_electrons':46,'target_index':255}.items()),
            'Exact central physical/electron identity differs')
    prior_task=failed['rows'][TID]['task']
    require(prior_task['physical_parent']==c['physical_parent'] and prior_task['positions']==c['row']['positions']
            and prior_task['physical_ids']==c['ids'],'Failed parent central geometry authority differs')
    state=dict(state,backend_formalism='UKS',initial_guess='PModel',initial_beta_rotation=[514,515,20,1,1],
               representation='external_orbital_def2_ECP46_UKS_singlet',physical_state_qualified=False)
    require(oom_task['electronic_state']==state and oom_task['physical_parent']==c['physical_parent']
            and oom_task['positions']==c['row']['positions'] and oom_task['physical_ids']==c['ids'],
            'Exact16611 physical/electronic state differs')
    evidence=dict(direction_helper=cfg['direction_helper'],direction_config=cfg['direction_config'],
        central_result=dc['central_result'],central_forces=dc['central_forces'],physical=c['physical_parent'],positions=c['row']['positions'],
        physical_ids=c['ids'],symbols=c['symbols'],state=state,basis=b,
        expected_basis_dimension=sum(b['functions_per_atom'][s] for s in c['symbols']),
        primary_calc=cfg['primary_calc'],primary_recipe=cfg['primary_recipe'],external_basis=cfg['external_basis'],
        saved_direction=saved['direction'],method_port_version='6.1.1',training_backend_version='6.0.0',
        failed_result=cfg['failed_result'],failed_input=cfg['failed_input'],failed_config=cfg['failed_config'],
        failed_implementation=cfg['failed_implementation'],deck_change='Only PAtom replaced by PModel',
        native_NL_treatment='SCNL selected automatically by unchanged EnGrad recipe; verify actual output',
        oom_result=cfg['oom_result'],oom_input=cfg['oom_input'],oom_config=cfg['oom_config'],
        oom_implementation=cfg['oom_implementation'],oom_xyz=cfg['oom_xyz'],
        recovery_change='Resource only: 24MPI/24CPUs to12MPI/24CPUs; exact16611 scientific template/XYZ/basis retained')
    return dict(cfg=cfg,saved=saved,world=saved['world'],np=np,cell=c,state=state,deck=recipe(cfg),evidence=evidence)


def task(w,xyz,deck,basis):
    c=w['cell']; directory=verify(deck).parent
    return dict(task_id=TID,source_id='7CCN',geometry='exact_step725_central',metal='La',xyz=xyz,input=deck,
        external_basis=basis,output_path=str(directory/'endpoint.out'),engrad_path=str(directory/'endpoint.engrad'),
        task_type='analytic_gradient',charge=-2,multiplicity=1,physical_multiplicity=1,electronic_state=w['state'],
        physical_parent=c['physical_parent'],positions=c['row']['positions'],physical_ids=c['ids'],atom_count=256,target_index=255,
        saved_central=w['evidence']['central_result'],caveats=CAVEATS)


def make_manifest(w,out):
    directory=out/TID; directory.mkdir(); cfg=w['cfg']; c=w['cell']; compact=w['world']['compact']
    xyz=directory/'physical.xyz'; xyz.write_text('256\nExact step725 central; published OMol25 recipe port\n'+''.join('%s %.17g %.17g %.17g\n'%(s,*v) for s,v in zip(c['symbols'],c['x'])))
    deck=directory/'endpoint.inp'; deck.write_text(w['deck']); basis=directory/'def2-tzvpd.bas'; shutil.copyfile(verify(cfg['external_basis']),basis)
    require(pin(xyz)['sha256']==cfg['oom_xyz']['sha256'] and pin(deck)['sha256']==cfg['oom_input']['sha256']
            and pin(basis)['sha256']==cfg['external_basis']['sha256'],'Exact16611 input bytes differ')
    t=task(w,pin(xyz),pin(deck),pin(basis)); write(directory/'TASK.json',t)
    m=dict(protocol_id=PROTOCOL,agreement=cfg['specification'],adapter=pin(__file__),inputs_config=cfg['_config_pin'],
        scientific_binding=binding(cfg),method=METHOD,
        external_basis=cfg['external_basis'],orca=compact['orca'],implementation=compact['implementation'],execution_policy=compact['execution_policy'],
        execution_resources={'mpi_ranks':RANKS,'allocated_cpus':ALLOCATED_CPUS,'cpus_per_task':CPUS_PER_TASK,'concurrent_tasks':1},tasks=[t],declared_tasks=[TID],caveats=CAVEATS)
    p=write(out/'MANIFEST.json',m); validate_manifest(w,p); return p


def validate_manifest(w,record):
    m=read(record); cfg=w['cfg']; compact=w['world']['compact']; c=w['cell']; np=w['np']
    require(m['protocol_id']==PROTOCOL and m['agreement']==cfg['specification'] and m['adapter']['sha256']==cfg['implementation']['sha256']
            and m['scientific_binding']==binding(cfg)==binding(read(m['inputs_config'])) and m['external_basis']==cfg['external_basis']
            and m['execution_resources']=={'mpi_ranks':RANKS,'allocated_cpus':ALLOCATED_CPUS,'cpus_per_task':CPUS_PER_TASK,'concurrent_tasks':1} and m['declared_tasks']==[TID]
            and m['method']==METHOD and m['caveats']==CAVEATS
            and [t['task_id'] for t in m['tasks']]==[TID] and all(m[k]==compact[k] for k in ('orca','implementation','execution_policy')),
            'Central native manifest differs')
    verify(m['adapter']); t=m['tasks'][0]; xyz=w['world']['model'].xyz(verify(t['xyz']))
    require(t==task(w,t['xyz'],t['input'],t['external_basis']) and verify(t['input']).read_text()==w['deck']
            and verify(t['external_basis']).name=='def2-tzvpd.bas' and t['external_basis']['sha256']==cfg['external_basis']['sha256']
            and t['xyz']['sha256']==cfg['oom_xyz']['sha256'] and t['input']['sha256']==cfg['oom_input']['sha256']
            and verify(t['xyz']).parent==verify(t['input']).parent==verify(t['external_basis']).parent
            and [a[0] for a in xyz]==c['symbols'] and np.array_equal([a[1:] for a in xyz],c['x']),
            'Exact unchanged geometry, recipe or adjacent external basis differs')
    return m


def signature(t):
    return {k:(v['sha256'] if k in ('input','xyz','external_basis') else v) for k,v in t.items() if k not in ('output_path','engrad_path')}


def released(w):
    cfg=w['cfg']; p=read(cfg['preparation_result']); pc=read(p['config']); verify(p['implementation'])
    require(p['protocol']==PROTOCOL and p['stage']=='prepare' and p['status']=='prepared'
            and p['implementation']['sha256']==cfg['implementation']['sha256'] and binding(pc)==binding(cfg)==p['scientific_binding']
            and read(p['input_evidence'])==w['evidence'] and p['selected_tasks']==[TID]
            and all(p[k]==0 for k in ('native_runner_invocations','new_DFT_calls_attempted','new_MACE_calls','new_ALPB_calls','new_GFN_calls','new_MM_calls')),
            'Actual unchanged zero-call preparation required')
    validate_manifest(w,p['manifest']); decision=read(cfg['execution_decision'])
    require(decision['protocol']==PROTOCOL and decision['status']=='approved_single_memory_recovery_central'
            and decision['scientific_binding']==binding(cfg) and decision['preparation_result']==cfg['preparation_result']
            and decision['permitted_execution_tasks']==[TID] and decision['mpi_ranks']==RANKS
            and decision['allocated_cpus']==ALLOCATED_CPUS and decision['cpus_per_task']==CPUS_PER_TASK
            and decision['maximum_new_DFT_calls']==1,
            'Root central-only execution release required')
    return p


def allocation(w,out):
    text=subprocess.check_output(['scontrol','show','job',os.environ['SLURM_JOB_ID'],'-o'],text=True)
    path=out/'SCONTROL_JOB.txt'; path.write_text(text)
    hit=re.search(r'(?:^|\s)AllocTRES=(\S+)',text); require(hit is not None,'Observed allocation required')
    tres=dict(x.split('=',1) for x in hit[1].split(',')); mem=re.fullmatch(r'(\d+)([KMGTP])',tres.get('mem',''))
    require(mem is not None and int(tres['cpu'])==ALLOCATED_CPUS
            and re.search(r'(?:^|\s)NumNodes=1(?:\s|$)',text)
            and re.search(r'(?:^|\s)NumTasks=12(?:\s|$)',text)
            and re.search(r'(?:^|\s)CPUs/Task=2(?:\s|$)',text)
            and os.environ.get('SLURM_NTASKS')=='12' and os.environ.get('SLURM_CPUS_PER_TASK')=='2',
            'Actual12MPI/24CPU/2CPUs-per-task single-node allocation required')
    memory=int(mem[1])*{'K':1/1024,'M':1,'G':1024,'T':1024**2,'P':1024**3}[mem[2]]
    require(memory==MEMORY_MiB,'Declared default96000MiB allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS',
            'NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')),'One numerical thread per rank required')
    os.environ['SLURM_MEM_PER_NODE']=str(int(memory)); os.environ['METAL_ENV_WORKERS']='1'
    policy=w['world']['modules']['render_orca_runtime_input'].allocation_policy(workers=1,ranks=RANKS)
    require(policy['workers']==1 and policy['ranks_per_worker']==RANKS and policy['allocated_cpus']==ALLOCATED_CPUS
            and policy['local_mpi_slots']==RANKS and policy['memory_source']=='SLURM_MEM_PER_NODE'
            and policy['allocated_memory_MiB']==MEMORY_MiB and policy['maxcore_MB']==MAXCORE_MB
            and policy['reserve_fraction']==.25,'Actual worker/rank/memory allocation differs')
    return dict(policy,scheduler_receipt=pin(path),allocated_TRES=tres)


def native_evidence(text):
    # Do not mistake echoed input for executed method diagnostics.
    lines=[(i+1,s) for i,s in enumerate(text.splitlines()) if not re.match(r'^\s*\|\s*\d+>',s)]
    patterns={
        'method':r'wB97M|libXC|Functional|ScalHFX|ScalDF|range.sep|RIJ.COSX|AuxJ',
        'basis_ECP':r'basis|ECP|core electrons|Group\s+\d|^\s*Atom\s+\d|^\s*[spdfg]\s*[-:]|^\s*\d+\s+[-+0-9.]+\s+[-+0-9.]+\s+\d\s*$',
        'SCF':r'AutoTRAH|SCF CONVERGED|SCF NOT CONVERGED|CONVERGENCE|ConvCheck|ConvForced|Tol[E GX]|TolErr|MaxIter|Thresh|TCut|CNVDIIS|CNVSOSCF|Final (energy|density)|Energy change|DIIS error|Orbital (gradient|rotation)',
        'guess_rotation':r'INITIAL GUESS|PModel|MODEL POTENTIAL|Atomic calculation|Returned from the atomic|PAtom|H.ckel|rotat|\b514\b.*\b515\b|READING.*ORBITAL|MOREAD',
        'grids_NL':r'SCNL|GRID|IntAcc|AngularGrid|GridPruning|WeightScheme|BFCut|WCut|VV10|DFT.NL|NL Energy|NL    Energy|SC\+NL|NL short|[Vv]d[Ww]|non.local',
        'spin_density_gap':r'S\*\*2|S\^2|S\(S\+1\)|<S|spin contamination|Expectation value|HOMO|LUMO|GAP|N\(Alpha\)|N\(Beta\)|N\(Total\)|Number of (Alpha|Beta)|[Nn]umber of (alpha|beta)|NAlpha|NBeta',
    }
    evidence={k:[dict(line=i,text=s) for i,s in lines if re.search(p,s,re.I)] for k,p in patterns.items()}
    evidence['note']='Raw line evidence; absence or unrecognized meaning is unavailable pending saved-output review, never inferred from an input echo.'
    return evidence,'\n'.join(s for _,s in lines)


def admit(w,record,out):
    m=validate_manifest(w,record); t=m['tasks'][0]; world=w['world']; np=w['np']
    _,rows=world['modules']['ggr_sensitivity'].executed(verify(record)); row=dict(rows[TID],task=t,task_id=TID)
    try:
        if Path(t['output_path']).exists():
            text=Path(t['output_path']).read_text(); evidence,body=native_evidence(text)
            row['native_evidence']=write(out/'NATIVE_EVIDENCE.json',evidence)
            try: row['SCF_evidence']=world['runtime'].observations(t['output_path'],world['modules']['mace_omol_vacuum'])
            except Exception as exc: row['SCF_evidence']=dict(status='unavailable',reason=repr(exc))
        if row['status']!='complete': return row
        checks=[(r'Program Version\s+(6\.1\.1)\b','6.1.1'),(r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)','-2'),
            (r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)','1'),(r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)','1030'),
            (r'Hartree-Fock type\s+HFTyp\s+\.{2,}\s+(\w+)','UHF'),
            (r'Basis Dimension\s+Dim\s+\.{2,}\s+(\d+)',str(w['evidence']['expected_basis_dimension'])),
            (r'Maximum # iterations\s+MaxIter\s+\.{2,}\s+(\d+)','300'),
            (r'AO-DIIS\s+CNVDIIS\s+\.{2,}\s+(\w+)','on'),(r'SOSCF\s+CNVSOSCF\s+\.{2,}\s+(\w+)','off')]
        for pattern,wanted in checks: require(re.findall(pattern,body)==[wanted],'Native state/settings unavailable or different: '+pattern)
        ecp=re.findall(r'Type\s+(\w+)\s+ECP(?:\s+(\S+))?\s+\(replacing\s+(\d+)\s+core electrons',body)
        require(len(ecp)==1 and ecp[0][0]=='La' and ecp[0][2]=='46','Actual La ECP46 required')
        require(not re.search(r'DFTD4|DFTD3|Reading point charge file|^\s*(?:CPCM|SMD|COSMO) SOLVATION MODEL|INITIAL GUESS: POLARIZED ATOMIC DENSITIES',body,re.M),
                'Undeclared correction/solvent/old guess')
        energies=re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+0-9.Ee]+)',body)
        require(len(energies)==1 and math.isfinite(float(energies[0])) and float(energies[0])==row['energy_hartree'],'One finite final total energy required')
        receipt=read(row['receipt']); sidecar=read(receipt['artifacts']['runtime_sidecar']); policy=sidecar['memory_policy']
        require(receipt['orca_executable']==world['compact']['orca'] and receipt['returncode']==0 and receipt['normal_termination']
                and receipt['scf_converged'] and receipt['expected_engrad_complete'] and receipt['parallelism']['nprocs']==RANKS
                and policy['workers']==1 and policy['ranks_per_worker']==RANKS and policy['memory_source']=='SLURM_MEM_PER_NODE'
                and policy['allocated_cpus']==ALLOCATED_CPUS and policy['local_mpi_slots']==RANKS
                and policy['allocated_memory_MiB']==MEMORY_MiB and policy['maxcore_MB']==MAXCORE_MB,
                'Native completion/resources differ')
        grad=world['modules']['affordable_response'].read_engrad(t['engrad_path'])
        require(grad['atom_count']==256 and grad['atomic_numbers'].tolist()==[NUMBERS[s] for s in w['cell']['symbols']]
                and abs(grad['energy_Ha']-row['energy_hartree'])<=1e-8
                and np.allclose(grad['coordinates_bohr']*BOHR,w['cell']['x'],atol=1e-6,rtol=0)
                and grad['gradient_Ha_per_bohr'].shape==(256,3) and np.isfinite(grad['gradient_Ha_per_bohr']).all(),
                'Full gradient identity/final-energy/coordinates differ')
        gradient=grad['gradient_Ha_per_bohr']*HA/BOHR
        method_flags=dict(functional_label=bool(re.search(r'wB97M.?V',body,re.I)),
            self_consistent_VV10=bool(re.search(r'Switching to SCNL|self.consistent DFT.NL',body,re.I)),
            post_SCF_VV10_heading=bool(re.search(r'post.SCF DFT.NL dispersion correction',body,re.I)),
            explicit_rotation=bool(re.search(r'(?:rotat[^\n]*514[^\n]*515|514[^\n]*515[^\n]*20)',body,re.I)),
            PModel_guess=bool(re.search(r'INITIAL GUESS:\s*MODEL POTENTIAL',body,re.I)))
        row.update(status='native_core_admitted_pending_method_review',printed_energy_hartree=energies[0],engrad=pin(t['engrad_path']),
            gradient_kcal_mol_per_A=gradient.tolist(),quantity='gradient_not_force',executed_ecp=ecp,electronic_state=w['state'],
            runtime_memory_policy=policy,method_evidence_flags=method_flags,
            method_review=dict(status='pending_saved_output_review',reason='First native instance: inspect actual functional, basis/ECP/auxiliary coefficients, rotation, grids/NL, spin/density/gap before releasing any further calls.'),
            physical_state_qualified=False,force_accuracy_qualified=False,force_consistency_qualified=False,caveats=CAVEATS,
            molecular_wall_seconds=(dt.datetime.fromisoformat(receipt['finished_at_utc'])-dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds())
    except Exception as exc: row.update(status='unavailable',admission_failure=repr(exc))
    return world['old'].canonical(row)


def comparison(w,row):
    if row['status']!='native_core_admitted_pending_method_review': return dict(status='unavailable',reason='Finite native core admission required')
    np=w['np']; d=w['evidence']['saved_direction']; u=np.asarray(d['vector']); ft=-np.asarray(row['gradient_kcal_mol_per_A'])
    fm=np.asarray(d['central_MACE_force']); fp=np.asarray(d['central_PBE0_force'])
    def describe(v):
        norms=np.linalg.norm(v,axis=1)
        return dict(RMS_Cartesian_kcal_mol_A=float(np.sqrt(np.mean(v*v))),max_atom_norm_kcal_mol_A=float(np.max(norms)),
                    projection_kcal_mol_A=float(np.sum(v*u)))
    return dict(status='descriptive_pending_method_review',quantity='force=-gradient',units='kcal/mol/Angstrom',
        source_direction=w['evidence']['saved_direction'],physical=w['cell']['physical_parent'],positions=w['cell']['row']['positions'],
        full256_forces=dict(pmodel_recipe_variant=ft.tolist(),MACE=fm.tolist(),PBE0=fp.tolist()),
        forces={k:describe(v) for k,v in [('pmodel_recipe_variant',ft),('MACE',fm),('PBE0',fp)]},
        differences={k:describe(v) for k,v in [('MACE_minus_port',fm-ft),('MACE_minus_PBE0',fm-fp),('port_minus_PBE0',ft-fp)]},
        physical_ids=w['cell']['ids'],force_accuracy_qualified=False,physical_state_qualified=False,caveats=CAVEATS)


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--stage',choices=('prepare','run'),required=True)
    p.add_argument('--config',required=True); p.add_argument('--output',required=True); args=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=False); tick=time.monotonic()
    result=dict(protocol=PROTOCOL,stage=args.stage,status='unavailable',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),declared_tasks=[TID],selected_tasks=[],
        rows={TID:dict(status='unavailable',reason='not_attempted')},comparison=dict(status='unavailable'),
        new_DFT_calls_attempted=0,native_runner_invocations=0,new_MACE_calls=0,new_ALPB_calls=0,new_GFN_calls=0,new_MM_calls=0,
        optimizer_steps=0,new_geometry_count=0,other_metal_calls=0,displaced_calls=0,caveats=CAVEATS,
        internal_atomic_guess_artifacts=None,requested_whole_molecule_endpoint_limit=1)
    try:
        cfg=read(result['config']); cfg['_config_pin']=result['config']; w=inputs(cfg)
        result.update(scientific_binding=binding(cfg),input_evidence=write(out/'INPUT_EVIDENCE.json',w['evidence']))
        if args.stage=='prepare':
            require(cfg.get('preparation_result') is None and cfg.get('execution_decision') is None,'Preparation cannot bind future outcomes')
            result.update(manifest=make_manifest(w,out),selected_tasks=[TID],status='prepared')
        else:
            prep=released(w); result.update(preparation_result=cfg['preparation_result'],execution_decision=cfg['execution_decision'])
            result['allocation']=allocation(w,out); record=make_manifest(w,out); t=read(record)['tasks'][0]
            require(signature(t)==signature(read(prep['manifest'])['tasks'][0]),'Native task differs from actual preparation')
            result.update(manifest=record,selected_tasks=[TID]); write(out/'START.json',result)
            try:
                result['native_runner_invocations']+=1
                w['world']['modules']['run_orca_task_manifest'].run_manifest(verify(record),orca_path=verify(w['world']['compact']['orca']),
                    workers=1,nprocs=RANKS,selected_tasks={TID},expected_runner_sha256=w['world']['compact']['execution_policy']['task_runner']['sha256'])
            except Exception as exc: result['runner_error']=repr(exc)
            finally:
                result['new_DFT_calls_attempted']=int(Path(t['output_path']).exists())
                atomic_files=[pin(p) for p in sorted(Path(t['output_path']).parent.glob('*_atom*')) if p.is_file()]
                result['internal_atomic_guess_artifacts']=write(out/'ATOMIC_GUESS_ARTIFACTS.json',dict(
                    files=atomic_files,scope='Native PModel initialization, not separately requested molecular endpoints',
                    internal_solver_invocations='Not inferred from artifact count; inspect native output',
                    costs='Included in endpoint and allocation totals'))
            row=admit(w,record,out); result['rows'][TID]=row; result['status']=row['status']
            result['comparison']=write(out/'FORCE_COMPARISON.json',comparison(w,row))
        verify(result['config'])
    except Exception as exc: result.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick,peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        physical_state_qualified=False,force_accuracy_qualified=False,force_consistency_qualified=False,no_automatic_retry=True)
    result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ.get('SLURM_CPUS_ON_NODE','0'))
    result['scheduler_costs']='Actual scheduler accounting is authoritative; allocated capacity is not utilization; parent RSS is not MPI peak.'
    write(out/'RESULT.json',result)
    return 0 if result['status'] in ('prepared','native_core_admitted_pending_method_review') else 1


if __name__=='__main__': raise SystemExit(main())
