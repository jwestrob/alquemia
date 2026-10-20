#!/usr/bin/env python3
"""One fixed La endpoint tests step92-to725 work and residual reference forces."""
import argparse
import datetime as dt
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

PROTOCOL = 'lbt_7ccn_converged_reference_v1'
IDS = ('old_La', 'new_La')
EV, HA, BOHR = 23.06054783061903, 627.509474, 0.529177210903
FIXED = {
    'reference_helper':'bf75f502a2e3c53e51329cfd59b9ff34f661d85acc0ed24483f51868bbf194cd',
    'reference_collection':'5908cb2d3b3fcb07b8ac764e3ddfb78fbc56beb29196790388cf79213af4ad72',
    'reference_config':'2fef75872b58eb11bb6afe5ae6cc38126b3cb46050868e3060d76074dab95aa1',
    'old_reference':'3a173d81072ab269922d17db1596c130e8132bd2f803614b00d2fe59455c10de',
    'continuation_result':'95d870ebbfdedede5166c7499e4b8181dd8c24289503bdb2e35eaba02c4a5068',
    'new_point':'6605f43cc3289509cee65c9be97a195abaad43d2d247ab488f5de889aa6f5419',
    'science_decision':'27377db9c2064650094e9f6106a106b57509a3701fa233cc7dfda5065e1be83b',
}
CAVEATS = {
    'selection':'Exact step725 model-force-converged terminal of the declared 7CCN-La continuation; old step92 native reference reused. Both endpoints are physically admitted; model convergence does not establish reference stationarity.',
    'source':'Lu-conditioned 7CCN peptide; consumed developmental fixture, not a blind protein benchmark.',
    'terminal_chemistry':'Conditional free NH3+/COO- and five deprotonated acidic side chains; experimental terminal specification unresolved.',
    'vacuum_dianion':'Post-added fixed-FF ALPB does not polarize the vacuum density; native convergence does not establish electron binding or a stable ground state.',
    'frozen_f':'Whole charge -2, physical La singlet and effective singlet, archived ECP46; no SOC or physical-state qualification.',
    'force':'Native analytic admission only; no independent finite-difference or force-accuracy qualification. Force differences and projections are descriptive.',
    'interpretation':'Incremental deformation work only; 1 kcal/mol is reporting resolution, not uncertainty or native admission. No minimum, affinity, preference or general accuracy claim.',
}


def require(ok, message):
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
    spec=importlib.util.spec_from_file_location(name,verify(record))
    obj=importlib.util.module_from_spec(spec); spec.loader.exec_module(obj); return obj


def binding(cfg):
    return {k:(v['sha256'] if k=='implementation' else v) for k,v in cfg.items()
            if k not in ('preflight_result','execution_decision','endpoint_result','_config_pin')}


def inputs(cfg):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    require(cfg['protocol']==PROTOCOL and pin(__file__)['sha256']==cfg['implementation']['sha256']
            and cfg['permitted_execution_tasks']==['new_La'] and cfg['mpi_ranks']==12,'One fixed 12-rank endpoint required')
    for k,sha in FIXED.items(): require(cfg[k]['sha256']==sha,'Frozen selection differs: '+k)
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}: verify(v)
    require(pin(sys.executable)['sha256']==cfg['python']['sha256'],'Original Python executable required')
    ref=module('converged_previous_reference',cfg['reference_helper'])
    collection=read(cfg['reference_collection']); ocfg=read(cfg['reference_config'])
    require(collection['status']=='complete' and collection['stage']=='collect' and collection['slurm_job_id']=='15456'
            and collection['full_two_available'] and collection['implementation']==cfg['reference_helper']
            and collection['config']==cfg['reference_config'] and ocfg['endpoint_result']==cfg['old_reference']
            and ocfg['python']==cfg['python'] and ocfg['plan']==cfg['plan'],'Previous incremental reference binding differs')
    ocfg['_config_pin']=cfg['reference_config']; w=ref.inputs(ocfg)
    parent,oldrow=ref.saved_run(w,cfg['old_reference'])
    require(parent['slurm_job_id']=='15321' and oldrow==collection['rows']['new_La']
            and oldrow['status']=='available','Exact step92 native reference replay differs')
    np=w['np']; oldcell=w['cells']['new_La']; oldsol=w['solvent']['new_La']; original_reference=w['reference']
    run=read(cfg['continuation_result']); rcfg=read(run['config']); verify(run['implementation'])
    require(run['status']=='complete' and run['stage']=='run' and run['slurm_job_id']=='15525'
            and run['declared_states']==[['7CCN','La']] and len(run['rows'])==1
            and run['implementation']['sha256']=='1fc6b168c0807a50223b2e5bb2b77837e19efd9d210d4c98f3b3a9e8785a400a'
            and run['config']['sha256']=='1bfaebc040879107ac72d5f7a3777591670c76857bd4bd9af57314fdd173b454',
            'Exact completed generous continuation required')
    rr=run['rows'][0]; search=read(rr['search']); cp=read(rr['checkpoint']); startcp=read(rcfg['checkpoint'])
    previous=read(rcfg['previous_result']); pcfg=read(rcfg['previous_config'])
    require(rr['status']=='complete' and rr['physical_state_admitted'] and rr['converged']
            and rr['termination']==search['termination']=='full_force_converged'
            and rr['terminal_point']==search['terminal_point']==cp['point']==cfg['new_point']
            and rr['additional_accepted_steps']==553 and rr['additional_distinct_evaluations']==567
            and search['accepted_steps']==cp['accepted_steps']==725 and search['distinct_evaluations']==cp['distinct_evaluations']==739
            and search['continuation_checkpoint']==rcfg['checkpoint']
            and previous['config']==rcfg['previous_config'] and previous['rows'][0]['checkpoint']==rcfg['checkpoint']
            and startcp['point']==previous['rows'][0]['terminal_point']
            and pcfg['point']==read(pcfg['checkpoint'])['point']==oldcell['row']['point']
            and cp['initial_point']==startcp['initial_point']
            and search['original_source']==w['inv']['original_source'] and not search['geometry_source_rebased']
            and rr['final_geometry']==search['final_geometry'] and rr['final_geometry']['pass']
            and rr['final_geometry']['state_compatible'] and cfg['new_point'] in run['points'],
            'Fixed converged terminal/step92 lineage/source admission differs')
    decision=read(cfg['science_decision'])
    require(decision['protocol']==PROTOCOL and decision['status']=='prepare_single_fixed_converged_endpoint_reference'
            and decision['continuation_result']==cfg['continuation_result'] and decision['new_terminal']==cfg['new_point']
            and decision['checkpoint']==rr['checkpoint'] and decision['reused_reference_result']==cfg['old_reference']
            and decision['reused_reference_collection']==cfg['reference_collection']
            and decision['previous_reference_adapter']==cfg['reference_helper']
            and decision['previous_reference_collection_config']==cfg['reference_config'] and decision['declared_new_cells']==1,
            'Declared endpoint selection differs')
    verify(decision['continuation_review'])
    p=read(cfg['new_point']); atoms=read(p['physical_source']); x=np.load(verify(p['positions']),allow_pickle=False)
    require(p['status']=='complete' and p['source_id']=='7CCN' and p['metal']=='La' and p['compute_force'] is True
            and p['recipe']==oldcell['point']['recipe']
            and [{k:v for k,v in a.items() if k!='xyz_A'} for a in atoms]
                ==[{k:v for k,v in a.items() if k!='xyz_A'} for a in oldcell['atoms']]
            and x.shape==(256,3) and np.isfinite(x).all() and np.array_equal(x,[a['xyz_A'] for a in atoms]),
            'New physical inventory/recipe differs')
    task=read(p['task']); symbols=oldcell['symbols']; ids=oldcell['ids']; st=w['old'].state(symbols,'La',w['model'].NUMBERS)
    xyz=w['model'].xyz(verify(p['xyz']))
    require(st==oldcell['state'] and task['physical_source']==p['physical_source'] and task['xyz']==p['xyz']
            and task['source_id']=='7CCN' and task['metal']=='La' and task['charge']==-2 and task['multiplicity']==1
            and task['metal_index']==255 and task['atom_count']==256 and task['representation']=='whole'
            and task['source_cap_map'] is None and task['all_electron_count']==st['all_electron_count']
            and [a[0] for a in xyz]==symbols and np.array_equal([a[1:] for a in xyz],x),'New native state/XYZ differs')
    key=hashlib.sha256(json.dumps({'recipe':p['recipe'],'positions_A':x.tolist()},sort_keys=True,
                               separators=(',',':'),allow_nan=False).encode()).hexdigest()
    require(key==p['cache_key'],'New cache identity differs')
    raw=read(p['MACE']); sc=raw['input_state_check']; f=np.asarray(raw['forces_eV_A'])
    require(sc['atoms']==256 and sc['charge']==-2 and sc['spin_multiplicity']==1 and sc['selected_metal_index']==255
            and sc['graph_count']==1 and sc['head']=='omol' and sc['charge_prediction_claimed'] is False
            and sc['metal_neighbor_edge_count']>0 and math.isfinite(raw['energy_eV'])
            and raw['force_convention']=='negative derivative of the same native vacuum energy'
            and f.shape==(256,3) and np.isfinite(f).all(),'Saved MACE energy/force/state differs')
    row=dict(cell_id='new_La',geometry='declared_step725_converged_terminal',metal='La',point=cfg['new_point'],
             physical=p['physical_source'],positions=p['positions'],MACE=p['MACE'],ALPB=p['ALPB'],task=p['task'],xyz=p['xyz'])
    c=dict(row=row,atoms=atoms,ids=ids,x=x,symbols=symbols,state=st,charges=oldcell['charges'],
           physical_parent=p['physical_source'],point=p,raw=raw)
    sol=w['old'].admit_alpb(c,{'ALPB_binary':p['recipe']['native_ALPB_binary']},np,w['model'])
    native=dict(MACE_eV=raw['energy_eV'],field_Eh=0.,MM_kcal_mol=0.,ALPB_Eh=sol['energy_hartree'],ALPB_parts_Eh=sol['components_hartree'])
    require(sol['status']=='available' and p['native_energies']==native
            and p['components']=={'MACE_eV':raw['energy_eV'],'fixed_ALPB_hartree':sol['energy_hartree']},'New native energy replay differs')
    g=np.load(verify(p['gradient_kcal_mol_A']),allow_pickle=False)
    require(g.shape==(256,3) and np.isfinite(g).all()
            and np.array_equal(-f*EV+np.asarray(sol['gradient_hartree_per_bohr'])*(HA/BOHR),g),'New raw composite gradient replay differs')
    c.update(MACE_force=f*EV,composite_force=-g)
    dm=(raw['energy_eV']-oldcell['raw']['energy_eV'])*EV; ds=(sol['energy_hartree']-oldsol['energy_hartree'])*HA
    model_work=dict(MACE_kcal_mol=dm,fixed_ALPB_kcal_mol=ds,composite_kcal_mol=math.fsum((dm,ds)))
    # This spans step92-to725; the run's work_from_checkpoint spans only step172-to725.
    chord=x-oldcell['x']; length=float(np.linalg.norm(chord)); require(length>0 and math.isfinite(length),'Nonzero fixed deformation required')
    directions=dict(physical_ids=ids,old_point=oldcell['row']['point'],new_point=cfg['new_point'],
                    normalization='Full256x3 Euclidean norm; no alignment',chord_length_A=length,
                    directions={'incremental_deformation_chord':(chord/length).tolist()})
    evidence=dict(old_point=oldcell['row']['point'],new_point=cfg['new_point'],old_native_result=cfg['old_reference'],
                  new_search=rr['search'],new_checkpoint=rr['checkpoint'],source=w['inv']['original_source'],
                  physical_guard=rr['final_geometry'],model_work=model_work,state=st,recipe=p['recipe'],
                  old_native_receipt=oldrow['receipt'],old_native_row=oldrow,source_rebased=False)
    w.update(cfg=cfg,reference=original_reference,cells={'old_La':oldcell,'new_La':c},solvent={'old_La':oldsol,'new_La':sol},
             directions=directions,evidence=evidence,oldrow=oldrow,old_parent=parent)
    return w


def native_task(w,xyz,deck):
    c=w['cells']['new_La']; d=verify(deck).parent
    return dict(task_id='new_La',source_id='7CCN',geometry=c['row']['geometry'],metal='La',xyz=xyz,input=deck,
        output_path=str(d/'endpoint.out'),engrad_path=str(d/'endpoint.engrad'),task_type='analytic_gradient',charge=-2,
        multiplicity=1,physical_multiplicity=1,electronic_state=c['state'],physical_parent=c['physical_parent'],
        positions=c['row']['positions'],physical_ids=c['ids'],atom_count=256,target_index=255,ALPB=c['row']['ALPB'],caveats=CAVEATS)


def prepare(w,out):
    d=out/'new_La'; d.mkdir(); c=w['cells']['new_La']; cfg=w['cfg']; compact=w['compact']
    xyz=d/'physical.xyz'; deck=d/'endpoint.inp'
    xyz.write_text('256\nFixed model-converged step725 7CCN-La terminal; no caps or field\n'+''.join('%s %.17g %.17g %.17g\n'%(s,*v) for s,v in zip(c['symbols'],c['x'])))
    deck.write_text(w['old'].asset_checks(compact,'La')); t=native_task(w,pin(xyz),pin(deck)); write(d/'TASK.json',t)
    m=dict(protocol_id=PROTOCOL,agreement=cfg['specification'],plan=cfg['plan'],method=w['old'].METHOD,adapter=pin(__file__),
        inputs_config=cfg['_config_pin'],scientific_binding=binding(cfg),compact_manifest=w['inv']['native_reference']['compact_assets_manifest'],
        assets=compact['assets'],orca=compact['orca'],implementation=compact['implementation'],execution_policy=compact['execution_policy'],
        execution_resources={'mpi_ranks':12,'concurrent_tasks':1},tasks=[t],caveats=CAVEATS,declared_cells=1,
        selection={'continuation_result':cfg['continuation_result'],'new_point':cfg['new_point'],'old_reference':cfg['old_reference']})
    rec=write(out/'MANIFEST.json',m); validate_manifest(w,rec,['new_La']); return rec


def validate_manifest(w,record,tids):
    m=read(record); cfg=w['cfg']; compact=w['compact']; c=w['cells']['new_La']
    require(list(tids)==['new_La'] and [t['task_id'] for t in m['tasks']]==['new_La'] and m['protocol_id']==PROTOCOL
            and m['scientific_binding']==binding(cfg)==binding(read(m['inputs_config']))
            and m['adapter']['sha256']==cfg['implementation']['sha256'] and m['agreement']==cfg['specification']
            and m['plan']==cfg['plan'] and m['method']==w['old'].METHOD
            and m['compact_manifest']==w['inv']['native_reference']['compact_assets_manifest']
            and all(m[k]==compact[k] for k in ('assets','orca','implementation','execution_policy'))
            and m['execution_resources']=={'mpi_ranks':12,'concurrent_tasks':1} and m['declared_cells']==1
            and m['caveats']==CAVEATS and m['selection']=={'continuation_result':cfg['continuation_result'],
                'new_point':cfg['new_point'],'old_reference':cfg['old_reference']},'One-cell native manifest differs')
    verify(m['adapter']); t=m['tasks'][0]; xyz=w['model'].xyz(verify(t['xyz']))
    require(t==native_task(w,t['xyz'],t['input']) and verify(t['xyz']).parent==verify(t['input']).parent
            and verify(t['input']).read_text()==w['old'].asset_checks(compact,'La')
            and verify(t['input']).read_text()==verify(w['oldrow']['task']['input']).read_text()
            and [a[0] for a in xyz]==c['symbols'] and w['np'].array_equal([a[1:] for a in xyz],c['x']),
            'Unchanged vacuum deck/physical coordinates required')
    return m


def preflight_gate(w):
    cfg=w['cfg']; p=read(cfg['preflight_result']); pcfg=read(p['config'])
    require(p['protocol']==PROTOCOL and p['stage']=='preflight' and p['status']=='preflight_passed'
            and p['implementation']['sha256']==cfg['implementation']['sha256'] and binding(pcfg)==binding(cfg)==p['scientific_binding']
            and p['input_evidence']==w['evidence'] and p['solvent_admission']==w['solvent']
            and read(p['directions'])==w['directions']
            and all(p[k]==0 for k in ('new_DFT_calls_attempted','native_runner_invocations','new_ALPB_calls','new_MACE_calls','new_GFN_calls','new_MM_calls')),
            'Actual exact zero-call preflight required')
    validate_manifest(w,p['manifest'],['new_La']); return p


def execution_gate(w):
    p=preflight_gate(w); cfg=w['cfg']; d=read(cfg['execution_decision'])
    require(d['protocol']==PROTOCOL and d['status']=='approved_single_fixed_reference'
            and d['scientific_binding']==binding(cfg) and d['preflight_result']==cfg['preflight_result']
            and d['permitted_execution_tasks']==['new_La'] and d['new_point']==cfg['new_point']
            and d['old_reference']==cfg['old_reference'] and d['mpi_ranks']==12,'Root one-cell execution review required')
    return p


# Native admission below is copied unchanged from pinned c0cca6dd::admit.
def admit(w, manifest_pin, tid):
    validate_manifest(w, manifest_pin, [tid])
    compact, modules, runtime = w['compact'], w['modules'], w['runtime']
    model = modules['metal_environment_frozen_f']; collector = modules['ggr_sensitivity']
    residual = modules['mace_omol_vacuum']
    manifest, rows = collector.executed(verify(manifest_pin))
    require([t['task_id'] for t in manifest['tasks']] == [tid], 'Scout execution denominator differs')
    task = manifest['tasks'][0]; row = rows[tid]
    row.update(task=task, task_id=tid, physical_multiplicity=task['physical_multiplicity'])
    try:
        if Path(task['output_path']).exists():
            row['SCF_evidence'] = runtime.observations(task['output_path'], residual)
        if row['status'] != 'complete':
            return w['old'].canonical(row)
        parsed = model.parse(task, {'electronic_state': task['electronic_state']})
        require(math.isfinite(parsed['energy_hartree']) and w['np'].asarray(parsed['gradient_kcal_mol_per_A']).shape == (256,3)
                and w['np'].isfinite(parsed['gradient_kcal_mol_per_A']).all(),
                'Full finite native gradient unavailable')
        text = Path(task['output_path']).read_text()
        require('INITIAL GUESS: MODEL POTENTIAL' in text and not re.search(r'Reading point charge file|^\s*(?:CPCM|SMD|COSMO) SOLVATION MODEL', text, re.M),
                'Vacuum/guess scope changed')
        require(row['SCF_evidence']['ConvCheckMode'] == 'Total+1el-Energy' and row['SCF_evidence']['ConvForced'] == 0,
                'Native SCF convergence semantics changed')
        receipt = read(row['receipt']); sidecar = read(receipt['artifacts']['runtime_sidecar'])
        require(receipt['orca_executable'] == compact['orca'] and receipt['returncode'] == 0
                and receipt['normal_termination'] and receipt['scf_converged'] and receipt['expected_engrad_complete'], 'Native execution failed')
        policy = sidecar['memory_policy']
        require(policy['workers'] == 1 and policy['ranks_per_worker'] == 12
                and policy['memory_source'] == 'SLURM_MEM_PER_NODE', 'Native allocation-derived memory policy differs')
        row.update(parsed, status='available', runtime_memory_policy=policy,
                   molecular_wall_seconds=(dt.datetime.fromisoformat(receipt['finished_at_utc']) -
                                           dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds(), caveats=CAVEATS)
    except Exception as exc:
        row.update(status='unavailable', admission_failure=repr(exc), energy_hartree=None)
    return w['old'].canonical(row)


def saved_run(w,record):
    p=read(record); pcfg=read(p['config']); verify(p['implementation'])
    require(p['protocol']==PROTOCOL and p['stage']=='run' and p['implementation']['sha256']==w['cfg']['implementation']['sha256']
            and binding(pcfg)==binding(w['cfg']) and p['selected_tasks'] in ([],['new_La'])
            and p['new_DFT_calls_attempted'] in (0,1) and p['native_runner_invocations'] in (0,1),'Saved singleton identity differs')
    execution_gate(dict(w,cfg=pcfg))
    require(p['preflight_result']==pcfg['preflight_result']==w['cfg']['preflight_result']
            and p['execution_decision']==pcfg['execution_decision']==w['cfg']['execution_decision'],
            'Saved singleton preflight/decision differs')
    if p['status']!='complete':
        return p,dict(task_id='new_La',status='unavailable',parent_result=record,reason='Failed singleton retained',
                      native_evidence=p['rows']['new_La'],failure=p.get('failure'),runner_error=p.get('runner_error'))
    require(p['selected_tasks']==['new_La'] and p['new_DFT_calls_attempted']==p['native_runner_invocations']==1,
            'Complete singleton must have exactly one native call')
    row=admit(w,p['manifest'],'new_La')
    require(row['status']=='available' and row==p['rows']['new_La'],'Saved singleton native admission differs')
    return p,row


def comparison(w,rows):
    if any(rows[k]['status']!='available' for k in IDS):
        return dict(status='unavailable',required_cells=list(IDS),reason='Both fixed endpoints required')
    dft=(rows['new_La']['energy_hartree']-rows['old_La']['energy_hartree'])*HA
    model=w['evidence']['model_work']; mace=model['MACE_kcal_mol']; sol=model['fixed_ALPB_kcal_mol']; error=mace-dft
    return dict(status='complete',definition='step725 minus step92; native differences before unit conversion',units='kcal/mol',
        PBE0_vacuum=dft,MACE_vacuum=mace,fixed_FF_ALPB=sol,PBE0_plus_fixed_ALPB=math.fsum((dft,sol)),
        MACE_plus_fixed_ALPB=model['composite_kcal_mol'],MACE_minus_PBE0_work_error=error,
        reporting_resolution_kcal_mol=1.,within_1_kcal_reporting_resolution=abs(error)<=1.,
        reference_composite_lowering=math.fsum((dft,sol))<0.,resolution_is_not_native_admission_or_uncertainty=True,
        force_accuracy_qualified=False,caveats=CAVEATS)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('preflight','run','collect'),required=True)
    parser.add_argument('--config',required=True); parser.add_argument('--output',required=True); args=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=False); tick=time.monotonic()
    result=dict(protocol=PROTOCOL,stage=args.stage,status='execution_failure',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),declared_new_cells=1,declared_comparison_cells=list(IDS),
        selected_tasks=[],rows={k:{'status':'unavailable','reason':'not_attempted'} for k in IDS},
        new_DFT_calls_attempted=0,native_runner_invocations=0,new_ALPB_calls=0,new_MACE_calls=0,new_GFN_calls=0,new_MM_calls=0,
        optimizer_steps=0,new_geometry_count=0,other_metal_calls=0,comparison={'status':'unavailable'},force_reports={},caveats=CAVEATS)
    try:
        cfg=read(result['config']); cfg['_config_pin']=result['config']; w=inputs(cfg)
        result.update(scientific_binding=binding(cfg),input_evidence=w['evidence'],solvent_admission=w['solvent'],
                      directions=write(out/'DIRECTIONS.json',w['directions']))
        result['rows']['old_La']=w['oldrow']; result['old_reference']=cfg['old_reference']
        if args.stage=='preflight':
            require(all(cfg.get(k) is None for k in ('preflight_result','execution_decision','endpoint_result')),
                    'Preflight cannot contain future success bindings')
            result['manifest']=prepare(w,out); result['status']='preflight_passed'
        elif args.stage=='run':
            require(cfg.get('endpoint_result') is None,'Existing endpoint binding forbids repeat execution')
            pre=execution_gate(w); result.update(preflight_result=cfg['preflight_result'],execution_decision=cfg['execution_decision'])
            result['allocation']=w['runtime'].allocation(w['modules']['render_orca_runtime_input'],out)
            record=prepare(w,out); t=read(record)['tasks'][0]
            require(w['reference'].task_signature(t)==w['reference'].task_signature(read(pre['manifest'])['tasks'][0]),
                    'Native task differs from actual preflight')
            result.update(manifest=record,selected_tasks=['new_La']); write(out/'START.json',result)
            try:
                result['native_runner_invocations']+=1
                w['modules']['run_orca_task_manifest'].run_manifest(verify(record),orca_path=verify(w['compact']['orca']),
                    workers=1,nprocs=12,selected_tasks={'new_La'},expected_runner_sha256=w['compact']['execution_policy']['task_runner']['sha256'])
            except Exception as exc: result['runner_error']=repr(exc)
            finally: result['new_DFT_calls_attempted']=int(Path(t['output_path']).exists())
            result['rows']['new_La']=admit(w,record,'new_La')
            result['status']='complete' if result['rows']['new_La']['status']=='available' else 'unavailable'
        else:
            execution_gate(w); result.update(preflight_result=cfg['preflight_result'],execution_decision=cfg['execution_decision'])
            if cfg.get('endpoint_result'):
                parent,row=saved_run(w,cfg['endpoint_result']); result['rows']['new_La']=row
                result['parent_cost']=dict(result=cfg['endpoint_result'],slurm_job_id=parent['slurm_job_id'],
                    elapsed_seconds=parent['elapsed_seconds'],application_allocated_CPU_seconds=parent['application_allocated_CPU_seconds'],
                    native_runner_invocations=parent['native_runner_invocations'],new_DFT_calls_attempted=parent['new_DFT_calls_attempted'])
            else: result['rows']['new_La']={'status':'unavailable','reason':'No actual endpoint result bound'}
            result['comparison']=comparison(w,result['rows'])
            result['status']='complete' if result['comparison']['status']=='complete' else 'unavailable'
        for tid,row in result['rows'].items():
            if row['status']=='available':
                report=w['reference'].force_report(w,tid,row); report['caveats']=CAVEATS
                result['force_reports'][tid]=write(out/(tid+'_FORCES.json'),report)
        result['full_two_available']=all(result['rows'][k]['status']=='available' for k in IDS)
        verify(result['config'])
    except Exception as exc:
        result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
        result['comparison']={'status':'unavailable','reason':'Execution/admission failure retained'}
    result.update(elapsed_seconds=time.monotonic()-tick,peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        parent_RSS_is_not_MPI_peak=True,physical_state_qualified=False,force_consistency_qualified=False)
    result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ.get('SLURM_CPUS_ON_NODE','0'))
    result['scheduler_costs']='Actual scheduler accounting remains authoritative; application timing is not utilization.'
    write(out/'RESULT.json',result)
    return 0 if result['status'] in ('complete','preflight_passed') else 1


if __name__=='__main__': raise SystemExit(main())
