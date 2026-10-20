#!/usr/bin/env python3
"""Two fixed whole-La singlepoints along one saved force-discrepancy direction."""
import argparse
import datetime as dt
from decimal import Decimal, localcontext
import hashlib
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

PROTOCOL='lbt_7ccn_reference_direction_v1'
IDS=('minus','plus')
H,HA,GEOM_ATOL=0.01,627.509474,1e-12
FIXED={
 'reference_helper':'f251ca1e9514fb0c76659faf603b03ce2daf3684204fe4e98ef44e5d448bb240',
 'reference_config':'e995bda312a070725daa866f9964c0df42d2bb4cd865f4ba787a44908ffe5f92',
 'central_result':'a53125b98e597ea0473e06d58316ac2b32c914e9647da595f163cf0aa76a8c19',
 'collection_result':'407bdf88e42c43cb2ae48ed1493cbbc008f5c0eb957628787a6aa1c1b795bb0d',
 'central_forces':'e6c0d039500c9deb87a53f31f974dfdc741c1cd801e481b51dc3316d100b7a7c',
}
CAVEATS={'scope':'One direction and one full-Cartesian step; no all-component derivative or physical-accuracy qualification.',
 'reference':'Unchanged PBE0-D4/frozen-f ECP46 La singlet, vacuum dianion, conditional termini; native SCF residual/state limitations retained.',
 'model':'MACE training Hamiltonian differs; no attribution of a failed check uniquely to model, SCF, grid, finite-step or electronic state.',
 'interpretation':'No equilibrium, solvent, preference, transfer or production claim; no retries or follow-on sweeps.'}


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


def module(record):
    spec=importlib.util.spec_from_file_location('direction_central_reference',verify(record))
    obj=importlib.util.module_from_spec(spec); spec.loader.exec_module(obj); return obj


def binding(cfg):
    return {k:(v['sha256'] if k=='implementation' else v) for k,v in cfg.items()
            if k not in ('preparation_result','execution_decision','endpoint_results','_config_pin')}


def inputs(cfg):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    require(cfg['protocol']==PROTOCOL and cfg['implementation']['sha256']==pin(__file__)['sha256']
            and cfg['step_A']==H and cfg['mpi_ranks']==12 and cfg['permitted_tasks']==list(IDS)
            and cfg['force_tolerance_absolute']==0.005 and cfg['force_tolerance_relative']==0.0001,
            'Frozen two-singlepoint scope differs')
    for k,sha in FIXED.items(): require(cfg[k]['sha256']==sha,'Frozen selection differs: '+k)
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}: verify(v)
    require(pin(sys.executable)['sha256']==cfg['python']['sha256'] and sys.prefix==cfg['python_prefix'],
            'Original reference Python required')
    ref=module(cfg['reference_helper']); rcfg=read(cfg['reference_config']); rcfg['_config_pin']=cfg['reference_config']
    require(rcfg['python']==cfg['python'],'Reference interpreter binding differs')
    world=ref.inputs(rcfg); parent,central=ref.saved_run(world,cfg['central_result']); collection=read(cfg['collection_result'])
    ccfg=read(collection['config']); verify(collection['implementation'])
    require(parent['status']=='complete' and parent['slurm_job_id']=='15545' and central['status']=='available'
            and collection['status']=='complete' and collection['stage']=='collect' and collection['slurm_job_id']=='15838'
            and collection['implementation']==cfg['reference_helper'] and ccfg['endpoint_result']==cfg['central_result']
            and ref.binding(ccfg)==ref.binding(rcfg) and collection['rows']['new_La']==central
            and collection['force_reports']['new_La']==cfg['central_forces'],'Actual central reference/collection differs')
    c=world['cells']['new_La']; np=world['np']; report=read(cfg['central_forces'])
    require(report['physical']==c['physical_parent'] and report['positions']==c['row']['positions']
            and [a['physical_id'] for a in report['force_table']]==c['ids'] and len(c['ids'])==256
            and report['quantity']=='force=-gradient' and report['units']=='kcal/mol/Angstrom'
            and report['numerical_force_accuracy_qualified'] is False,'Central force-table identity/units differ')
    force=lambda key: np.asarray([a['force_kcal_mol_A'][key] for a in report['force_table']],dtype=float)
    fr=force('whole_PBE0_vacuum'); fm=force('whole_MACE_vacuum')
    require(fr.shape==fm.shape==(256,3) and np.isfinite(fr).all() and np.isfinite(fm).all()
            and np.array_equal(fr,-np.asarray(central['gradient_kcal_mol_per_A']))
            and np.array_equal(fm,c['MACE_force']),'Central full force replay differs')
    delta=fm-fr; norm=float(np.linalg.norm(delta)); require(norm>0 and math.isfinite(norm),'Nonzero discrepancy direction required')
    u=delta/norm
    require(abs(float(np.linalg.norm(u))-1)<=GEOM_ATOL and float(np.sum(u*delta))>0,'Direction norm/sign differs')
    x=c['x']; displaced={tid:x+sign*H*u for tid,sign in [('minus',-1),('plus',1)]}
    for tid,sign in [('minus',-1),('plus',1)]:
        actual=displaced[tid]-x
        require(np.isfinite(displaced[tid]).all() and np.allclose(actual,sign*H*u,atol=GEOM_ATOL,rtol=0)
                and abs(float(np.linalg.norm(actual))-H)<=GEOM_ATOL
                and float(np.max(np.linalg.norm(actual,axis=1)))<=H+GEOM_ATOL,'Exact declared displacement differs')
    text=verify(central['task']['input']).read_text()
    require(text.count(' EnGrad\n')==1,'Unique original EnGrad token required')
    deck=text.replace(' EnGrad\n','\n',1)
    require('EnGrad' not in deck and deck.endswith('* xyzfile -2 1 physical.xyz\n'),'Only singlepoint task change permitted')
    direction=dict(central_result=cfg['central_result'],collection_result=cfg['collection_result'],central_forces=cfg['central_forces'],
        physical=c['physical_parent'],positions=c['row']['positions'],physical_ids=c['ids'],state=c['state'],step_A=H,
        selection='u=(F_MACE_vac-F_PBE0_vac)/Euclidean_full256x3_norm; no alignment or projection',
        vector=u.tolist(),discrepancy_norm_kcal_mol_A=norm,central_PBE0_force=fr.tolist(),central_MACE_force=fm.tolist(),
        PBE0_projection_kcal_mol_A=float(np.sum(fr*u)),MACE_projection_kcal_mol_A=float(np.sum(fm*u)),
        original_native_gradient=central['engrad'],central_native_receipt=central['receipt'],central_SCF_evidence=central['SCF_evidence'],
        geometric_bookkeeping_absolute_tolerance=GEOM_ATOL)
    return dict(cfg=cfg,ref=ref,world=world,np=np,cell=c,central=central,direction=direction,displaced=displaced,deck=deck)


def task(w,tid,xyz,deck,physical,positions):
    c=w['cell']; directory=verify(deck).parent
    return dict(task_id=tid,source_id='7CCN',geometry='step725_discrepancy_'+tid,metal='La',xyz=xyz,input=deck,
        output_path=str(directory/'endpoint.out'),engrad_path=None,task_type='single_point',charge=-2,multiplicity=1,physical_multiplicity=1,
        electronic_state=c['state'],physical_parent=physical,positions=positions,physical_ids=c['ids'],atom_count=256,target_index=255,
        central_result=w['cfg']['central_result'],central_physical=c['physical_parent'],step_A=H,
        direction_sign=-1 if tid=='minus' else 1,caveats=CAVEATS)


def make_manifest(w,out,tids):
    tasks=[]; c=w['cell']; compact=w['world']['compact']; cfg=w['cfg']
    for tid in tids:
        directory=out/tid; directory.mkdir(); x=w['displaced'][tid]
        xyz=directory/'physical.xyz'; xyz.write_text('256\nDeclared whole-La discrepancy '+tid+'\n'+''.join('%s %.17g %.17g %.17g\n'%(s,*v) for s,v in zip(c['symbols'],x)))
        deck=directory/'endpoint.inp'; deck.write_text(w['deck'])
        physical=write(directory/'PHYSICAL.json',[dict(a,xyz_A=v.tolist()) for a,v in zip(c['atoms'],x)])
        positions=write(directory/'POSITIONS.json',dict(source_ids=c['ids'],positions_A=x.tolist()))
        t=task(w,tid,pin(xyz),pin(deck),physical,positions); write(directory/'TASK.json',t); tasks.append(t)
    m=dict(protocol_id=PROTOCOL,agreement=cfg['specification'],adapter=pin(__file__),inputs_config=cfg['_config_pin'],
        scientific_binding=binding(cfg),method=w['world']['old'].METHOD.removesuffix(' EnGrad'),
        assets=compact['assets'],orca=compact['orca'],implementation=compact['implementation'],execution_policy=compact['execution_policy'],
        execution_resources={'mpi_ranks':12,'concurrent_tasks':1},tasks=tasks,declared_tasks=list(IDS),caveats=CAVEATS)
    rec=write(out/'MANIFEST.json',m); validate_manifest(w,rec,tids); return rec


def validate_manifest(w,record,tids):
    m=read(record); cfg=w['cfg']; compact=w['world']['compact']; c=w['cell']; np=w['np']
    require(list(tids) in (list(IDS),['minus'],['plus']) and [t['task_id'] for t in m['tasks']]==list(tids)
            and m['protocol_id']==PROTOCOL and m['agreement']==cfg['specification']
            and m['adapter']['sha256']==cfg['implementation']['sha256'] and m['scientific_binding']==binding(cfg)==binding(read(m['inputs_config']))
            and m['method']==w['world']['old'].METHOD.removesuffix(' EnGrad') and m['declared_tasks']==list(IDS)
            and m['execution_resources']=={'mpi_ranks':12,'concurrent_tasks':1} and m['caveats']==CAVEATS
            and all(m[k]==compact[k] for k in ('assets','orca','implementation','execution_policy')),'Frozen native manifest differs')
    verify(m['adapter'])
    for t in m['tasks']:
        tid=t['task_id']; x=w['displaced'][tid]; xyz=w['world']['model'].xyz(verify(t['xyz']))
        require(t==task(w,tid,t['xyz'],t['input'],t['physical_parent'],t['positions'])
                and verify(t['input']).read_text()==w['deck'] and verify(t['input']).parent==verify(t['xyz']).parent
                and [a[0] for a in xyz]==c['symbols'] and np.array_equal([a[1:] for a in xyz],x)
                and read(t['physical_parent'])==[dict(a,xyz_A=v.tolist()) for a,v in zip(c['atoms'],x)]
                and read(t['positions'])==dict(source_ids=c['ids'],positions_A=x.tolist()),'Exact displaced physical input/deck differs')
    return m


def signature(t):
    return {k:(v['sha256'] if k in ('xyz','input','physical_parent','positions') else v)
            for k,v in t.items() if k!='output_path'}


def prepared(w):
    cfg=w['cfg']; p=read(cfg['preparation_result']); pcfg=read(p['config']); verify(p['implementation'])
    require(p['protocol']==PROTOCOL and p['stage']=='prepare' and p['status']=='prepared'
            and p['implementation']['sha256']==cfg['implementation']['sha256'] and p['scientific_binding']==binding(cfg)==binding(pcfg)
            and read(p['direction'])==w['direction'] and p['selected_tasks']==list(IDS)
            and p['new_DFT_calls_attempted']==p['native_runner_invocations']==0,'Actual unchanged zero-call preparation required')
    validate_manifest(w,p['manifest'],IDS); return p


def released(w,tid):
    p=prepared(w); cfg=w['cfg']; d=read(cfg['execution_decision'])
    require(d['protocol']==PROTOCOL and d['status']=='approved_two_fixed_singlepoints'
            and d['scientific_binding']==binding(cfg) and d['preparation_result']==cfg['preparation_result']
            and d['permitted_execution_tasks']==list(IDS) and d['mpi_ranks']==12
            and d['central_result']==cfg['central_result'] and d['maximum_new_DFT_calls']==2 and tid in IDS,
            'Root fixed singleton release required')
    return p


def admit(w,record,tid):
    manifest=validate_manifest(w,record,[tid]); world=w['world']; mods=world['modules']; t=manifest['tasks'][0]
    _,rows=mods['ggr_sensitivity'].executed(verify(record)); row=rows[tid]; row.update(task=t,task_id=tid)
    try:
        if Path(t['output_path']).exists(): row['SCF_evidence']=world['runtime'].observations(t['output_path'],mods['mace_omol_vacuum'])
        if row['status']!='complete': return world['old'].canonical(row)
        text=verify(row['output']).read_text(); st=t['electronic_state']
        for pattern,wanted in [(r'Program Version\s+(6\.1\.1)\b','6.1.1'),
            (r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)',str(st['charge'])),
            (r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)','1'),
            (r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)',str(st['explicit_electrons'])),
            (r'Hartree-Fock type\s+HFTyp\s+\.{2,}\s+(\w+)','RHF')]:
            require(re.findall(pattern,text)==[wanted],'Executed state/version differs: '+pattern)
        ecp=re.findall(r'Type\s+(\w+)\s+ECP(?:\s+(\S+))?\s+\(replacing\s+(\d+)\s+core electrons',text)
        require(len(ecp)==1 and ecp[0][0]=='La' and ecp[0][2]==str(st['ecp_core_electrons'])
                and 'DFTD4' in text and 'INITIAL GUESS: MODEL POTENTIAL' in text
                and not re.search(r'Reading point charge file|gCP correction\s+[-+0-9.]|^\s*(?:CPCM|SMD|COSMO) SOLVATION MODEL',text,re.M),
                'Executed vacuum/ECP/dispersion/guess differs')
        require(row['SCF_evidence']['ConvCheckMode']=='Total+1el-Energy' and row['SCF_evidence']['ConvForced']==0,
                'Native SCF convergence semantics changed')
        energies=re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+0-9.Ee]+)',text)
        require(len(energies)==1 and math.isfinite(float(energies[0])) and float(energies[0])==row['energy_hartree'],
                'Final printed native energy differs')
        receipt=read(row['receipt']); sidecar=read(receipt['artifacts']['runtime_sidecar']); policy=sidecar['memory_policy']
        require(receipt['orca_executable']==world['compact']['orca'] and receipt['returncode']==0
                and receipt['normal_termination'] and receipt['scf_converged']
                and 'engrad' not in receipt['artifacts'] and receipt['parallelism']['nprocs']==12,
                'Native singlepoint execution failed')
        require(policy['workers']==1 and policy['ranks_per_worker']==12 and policy['memory_source']=='SLURM_MEM_PER_NODE',
                'Original allocation-derived memory policy differs')
        row.update(status='available',printed_energy_hartree=energies[0],executed_ecp=ecp,electronic_state=st,
            physical_multiplicity=1,runtime_memory_policy=policy,
            molecular_wall_seconds=(dt.datetime.fromisoformat(receipt['finished_at_utc'])-dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds(),
            physical_state_qualified=False,force_consistency_qualified=False,caveats=CAVEATS)
    except Exception as exc: row.update(status='unavailable',admission_failure=repr(exc),energy_hartree=None)
    return world['old'].canonical(row)


def saved_run(w,record,tid):
    p=read(record); cfg=read(p['config']); verify(p['implementation'])
    require(p['protocol']==PROTOCOL and p['stage']=='run-one' and p['requested_task']==tid
            and p['implementation']['sha256']==w['cfg']['implementation']['sha256']
            and binding(cfg)==binding(w['cfg'])
            and cfg['preparation_result']==w['cfg']['preparation_result']
            and cfg['execution_decision']==w['cfg']['execution_decision'], 'Saved singleton identity differs')
    released(dict(w,cfg=cfg),tid)
    if p['status']!='complete':
        require(p['selected_tasks'] in ([],[tid]) and p['new_DFT_calls_attempted'] in (0,1)
                and p['native_runner_invocations'] in (0,1)
                and (p['selected_tasks']==[tid] or p['new_DFT_calls_attempted']==p['native_runner_invocations']==0),
                'Failed singleton call ledger differs')
        return p,dict(status='unavailable',reason='Failed singleton retained',parent_result=record,
                      native_evidence=p['rows'][tid],failure=p.get('failure'),runner_error=p.get('runner_error'))
    require(p['selected_tasks']==[tid] and p['scientific_binding']==binding(cfg)
            and p['preparation_result']==cfg['preparation_result'] and p['execution_decision']==cfg['execution_decision'],
            'Complete singleton release differs')
    row=admit(w,p['manifest'],tid)
    require(row['status']=='available' and row==p['rows'][tid] and p['new_DFT_calls_attempted']==p['native_runner_invocations']==1,
            'Saved singlepoint replay differs')
    return p,row


def compare(direction,rows):
    if any(rows[t]['status']!='available' for t in IDS):
        return dict(status='unavailable',required_tasks=list(IDS),reason='Both declared native energies required')
    with localcontext() as ctx:
        ctx.prec=50
        minus,plus=(Decimal(rows[t]['printed_energy_hartree']) for t in IDS)
        fd=-(plus-minus)*Decimal(str(HA))/(2*Decimal(str(H)))
        analytic=Decimal(str(direction['PBE0_projection_kcal_mol_A'])); delta=fd-analytic
        tol=Decimal('0.005')+Decimal('0.0001')*max(abs(fd),abs(analytic)); passed=abs(delta)<=tol
        quanta={t:Decimal(1).scaleb(Decimal(rows[t]['printed_energy_hartree']).as_tuple().exponent) for t in IDS}
        bound=(quanta['minus']+quanta['plus'])/2*Decimal(str(HA))/(2*Decimal(str(H)))
        return dict(status='passed' if passed else 'failed_directional_consistency',passed=passed,
            definition='force_FD=-(E_plus-E_minus)/(2h); signed residual=FD-minus-central-analytic',units='kcal/mol/Angstrom',
            printed_energy_hartree={t:rows[t]['printed_energy_hartree'] for t in IDS},
            energy_difference_hartree_decimal=str(plus-minus),FD_force_decimal=str(fd),analytic_force_decimal=str(analytic),
            signed_residual_decimal=str(delta),tolerance_decimal=str(tol),FD_force_kcal_mol_A=float(fd),
            analytic_force_kcal_mol_A=float(analytic),signed_residual_kcal_mol_A=float(delta),tolerance_kcal_mol_A=float(tol),
            MACE_force_projection_kcal_mol_A=direction['MACE_projection_kcal_mol_A'],
            printed_energy_quantum_hartree={k:str(v) for k,v in quanta.items()},printed_rounding_force_bound_kcal_mol_A=str(bound),
            force_accuracy_qualified=False,physical_state_qualified=False,caveats=CAVEATS)


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--stage',choices=('prepare','run-one','collect'),required=True)
    p.add_argument('--task',choices=IDS); p.add_argument('--config',required=True); p.add_argument('--output',required=True); args=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=False); start=time.monotonic()
    result=dict(protocol=PROTOCOL,stage=args.stage,status='unavailable',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),declared_tasks=list(IDS),selected_tasks=[],requested_task=args.task,
        rows={t:dict(status='unavailable',reason='not_attempted') for t in IDS},comparison={'status':'unavailable'},
        new_DFT_calls_attempted=0,native_runner_invocations=0,new_MACE_calls=0,new_ALPB_calls=0,new_GFN_calls=0,new_MM_calls=0,
        optimizer_steps=0,central_calls=0,other_metal_calls=0,caveats=CAVEATS)
    try:
        cfg=read(result['config']); cfg['_config_pin']=result['config']; w=inputs(cfg)
        result.update(scientific_binding=binding(cfg),direction=write(out/'DIRECTION.json',w['direction']))
        if args.stage=='prepare':
            require(args.task is None and cfg['preparation_result'] is None and cfg['execution_decision'] is None
                    and cfg['endpoint_results']==dict.fromkeys(IDS),'Preparation cannot bind future outcomes')
            result.update(manifest=make_manifest(w,out,IDS),selected_tasks=list(IDS),status='prepared')
        elif args.stage=='run-one':
            require(args.task in IDS and cfg['endpoint_results']==dict.fromkeys(IDS),'One fresh singleton required')
            tid=args.task; prep=released(w,tid); result.update(preparation_result=cfg['preparation_result'],execution_decision=cfg['execution_decision'])
            world=w['world']; result['allocation']=world['runtime'].allocation(world['modules']['render_orca_runtime_input'],out)
            record=make_manifest(w,out,[tid]); t=read(record)['tasks'][0]
            old=next(t for t in read(prep['manifest'])['tasks'] if t['task_id']==tid)
            require(signature(t)==signature(old),'Singleton differs from prepared task')
            result.update(manifest=record,selected_tasks=[tid]); write(out/'START.json',result)
            try:
                result['native_runner_invocations']+=1
                world['modules']['run_orca_task_manifest'].run_manifest(verify(record),orca_path=verify(world['compact']['orca']),
                    workers=1,nprocs=12,selected_tasks={tid},expected_runner_sha256=world['compact']['execution_policy']['task_runner']['sha256'])
            except Exception as exc: result['runner_error']=repr(exc)
            finally: result['new_DFT_calls_attempted']=int(Path(t['output_path']).exists())
            result['rows'][tid]=admit(w,record,tid); result['status']='complete' if result['rows'][tid]['status']=='available' else 'unavailable'
        else:
            require(args.task is None and set(cfg['endpoint_results'])==set(IDS),'Fixed two-slot collection required')
            prepared(w); result.update(preparation_result=cfg['preparation_result'],execution_decision=cfg['execution_decision'],parent_costs={})
            for tid in IDS:
                if cfg['endpoint_results'][tid] is None: continue
                parent,row=saved_run(w,cfg['endpoint_results'][tid],tid); result['rows'][tid]=row
                result['parent_costs'][tid]=dict(result=cfg['endpoint_results'][tid],slurm_job_id=parent['slurm_job_id'],
                    elapsed_seconds=parent['elapsed_seconds'],application_allocated_CPU_seconds=parent['application_allocated_CPU_seconds'])
            result['comparison']=compare(w['direction'],result['rows'])
            result['status']='complete' if result['comparison']['status']=='passed' else result['comparison']['status']
        verify(result['config'])
    except Exception as exc: result.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-start,peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        physical_state_qualified=False,force_accuracy_qualified=False,no_automatic_retry=True)
    result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ.get('SLURM_CPUS_ON_NODE','0'))
    result['scheduler_costs']='Actual scheduler accounting is authoritative; allocated capacity is not utilization; parent RSS is not MPI peak.'
    write(out/'RESULT.json',result); return 0 if result['status'] in ('prepared','complete') else 1


if __name__=='__main__': raise SystemExit(main())
