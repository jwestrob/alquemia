#!/usr/bin/env python3
"""Three supplemental Dy forces; eight inherited reference cells; Slurm only."""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import socket
import sys
import time
import traceback

PROTOCOL='lbt_dy_force_completion_v1'
IDS=tuple(s+'_'+t for s in ('7CCO','7CCN') for t in ('before_La','after_La','before_Dy','after_Dy'))
MISSING=('7CCO_before_Dy','7CCN_before_Dy','7CCN_after_Dy')
EV,HA,BOHR=23.06054783061903,627.509474,0.529177210903
ENERGY_TOL=1e-7
NATIVE_SHA='0e48b8c6c2f271512f2c445f826036d8a68f7dca0aabc5923697ab3ecb7d9896'
MCFG_SHA='d3c28300c8c7ef88ef747b733603a344790b932db56aa53b0ccd4e999a7a4023'
SCIENCE_SHA='7c58e45bab117f93d5182ab48ea2d7501671476dc2f6f05ea3a236f65ed4ef21'
INVENTORY_SHA='e2c196d1b5961cc74d50647e0565780928c3c1c103b282d563833f4b7e003957'


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
    spec=importlib.util.spec_from_file_location('dy_force_native',verify(record))
    obj=importlib.util.module_from_spec(spec); spec.loader.exec_module(obj); return obj


def vectors(np,value):
    a=np.asarray(value,dtype=float)
    require(a.shape==(256,3) and np.isfinite(a).all(),'Finite full256 vector array required')
    return a


def state(raw,metal):
    s=raw['input_state_check']
    require(s['atoms']==256 and s['charge']==-2 and s['spin_multiplicity']==(1 if metal=='La' else 6)
            and s['selected_metal_index']==255 and s['graph_count']==1 and s['head']=='omol'
            and s['charge_prediction_claimed'] is False and s['metal_neighbor_edge_count']>0
            and math.isfinite(raw['energy_eV']),'Native MACE state/energy differs')


def prepare(cfg):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    require(cfg['protocol']==PROTOCOL and pin(__file__)['sha256']==cfg['implementation']['sha256']
            and cfg['energy_replay_tolerance_eV']==ENERGY_TOL and cfg['maximum_new_MACE_calls']==3,
            'Frozen force-completion scope differs')
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}: verify(v)
    require(cfg['science_decision']['sha256']==SCIENCE_SHA and socket.gethostname()==cfg['qualified_host']
            and sys.prefix==cfg['python_prefix'] and pin(sys.executable)['sha256']==cfg['python']['sha256']
            and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==12
            and os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles'),'Original allocated GPU runtime required')
    decision=read(cfg['science_decision']); inv=read(cfg['inventory'])
    require(decision['protocol']==PROTOCOL and decision['status']=='declared_three_fixed_force_supplements'
            and decision['specification']==cfg['specification'] and decision['energy_replay_tolerance_eV']==ENERGY_TOL
            and decision['maximum_new_MACE_force_calls']==3
            and [s+'_'+t for s,t in decision['selected_missing_cells']]==list(MISSING)
            and set(s+'_'+t for s,t in decision['full_comparison_cells'])==set(IDS),'Root declaration differs')
    require(cfg['inventory']['sha256']==INVENTORY_SHA
            and inv['protocol']=='lbt_dy_force_completion_inventory_v1' and inv['scientific_execution_performed'] is False
            and [r['cell_id'] for r in inv['cells']]==list(IDS),'Eight frozen cells required')
    parents=inv['parents']; mcfg=read(parents['MACE_config'])
    require(parents['MACE_config']['sha256']==MCFG_SHA and parents['native_MACE_helper']['sha256']==NATIVE_SHA
            and mcfg['energy_replay_tolerance_eV']==ENERGY_TOL,'Original native force precedent differs')
    verify(parents['native_MACE_helper']); verify(mcfg['model']['checkpoint'])
    for v in mcfg['helper_pins'].values(): verify(v)
    require({k:importlib.metadata.version(k) for k in ('torch','mace-torch','ase')}==mcfg['expected_versions'],
            'Original MACE package versions required')
    import numpy as np
    collectors={}; directions={}
    for source in ('7CCO','7CCN'):
        p=parents[source]; c=read(p['collector']); cc=read(p['config']); verify(p['implementation']); read(p['inventory'])
        require(c['status']=='complete' and c['full_four_available'] and c['config']==p['config']
                and c['implementation']==p['implementation'] and cc['implementation']['sha256']==p['implementation']['sha256']
                and c['directions']==p['directions'] and c['physical_state_qualified'] is False
                and c['force_consistency_qualified'] is False,'Admitted reference collector differs')
        collectors[source]=c; directions[source]=read(p['directions'])
    mp=parents['whole_model_completion']; completed=read(mp['result']); mcc=read(mp['config']); verify(mp['implementation'])
    require(completed['status']=='complete' and completed['config']==mp['config']
            and completed['implementation']==mp['implementation'] and mcc['MACE_config']==parents['MACE_config']
            and mcc['native_MACE_helper']==parents['native_MACE_helper'],'Completed whole-model authority differs')
    cells={}
    for r in inv['cells']:
        cid,sid,tid,metal=r['cell_id'],r['source_id'],r['local_cell_id'],r['metal']
        require(cid==sid+'_'+tid and metal==tid.split('_')[1],'Cell/source/metal identity differs')
        c=collectors[sid]; nr=c['rows'][tid]; ref=r['reference']; raw=read(r['MACE']['receipt']); task=read(r['MACE']['task'])
        atoms=read(r['physical']); ids=[a['id'] for a in atoms]; x=vectors(np,[a['xyz_A'] for a in atoms])
        require(len(ids)==len(set(ids))==256 and ids==r['physical_ids'] and ids[-1]=='metal','Physical ID inventory differs')
        pp=verify(r['positions']); pdata=np.load(pp,allow_pickle=False) if pp.suffix=='.npy' else json.loads(pp.read_text())['positions_A']
        require(np.array_equal(x,vectors(np,pdata)),'Saved physical/position coordinates differ')
        symbols=[metal if a['element']=='M' else a['element'] for a in atoms]
        lines=verify(r['MACE']['xyz']).read_text().splitlines(); xyz=[a.split() for a in lines[2:] if a.strip()]
        require(int(lines[0])==len(xyz)==256 and [a[0] for a in xyz]==symbols
                and np.array_equal(x,vectors(np,[[float(v) for v in a[1:]] for a in xyz])),'Saved MACE XYZ differs')
        require(task['physical_source']==r['MACE']['physical_source'] and task['xyz']==r['MACE']['xyz']
                and task['source_id']==sid and task['metal']==metal and task['charge']==-2
                and task['multiplicity']==(1 if metal=='La' else 6) and task['metal_index']==255
                and task['atom_count']==256 and task['representation']=='whole' and task['source_cap_map'] is None,
                'Original whole MACE task differs')
        ma=read(r['MACE']['physical_source'])
        require([a['id'] for a in ma]==ids and np.array_equal(vectors(np,[a['xyz_A'] for a in ma]),x),
                'MACE/reference physical coordinates differ')
        recipe=r['MACE']['recipe']
        require(recipe['model']==mcfg['model'] and recipe['MACE_versions']==mcfg['expected_versions']
                and recipe['configuration']==parents['MACE_config']
                and all(task.get(k)==v for k,v in r['MACE']['state'].items()),'Original model/physical state differs')
        state(raw,metal)
        require(raw['energy_eV']==r['MACE']['energy_eV']
                and r['MACE']['force_available']==(cid not in MISSING),'Saved force availability/state/energy differs')
        if r['MACE']['point']:
            point=read(r['MACE']['point'])
            require(point['status']=='complete' and point['MACE']==r['MACE']['receipt'] and point['task']==r['MACE']['task']
                    and recipe['source_POINT_recipe']==r['MACE']['point']
                    and point['compute_force']==r['MACE']['original_compute_force']
                    and point['physical_source']==r['MACE']['physical_source'] and point['xyz']==r['MACE']['xyz']
                    and point['ALPB']==r['ALPB']['receipt'] and point['native_energies']['MACE_eV']==raw['energy_eV'],
                    'Original POINT components differ')
            pr=point['recipe']
            require(pr['model']==mcfg['model'] and pr['MACE_versions']==mcfg['expected_versions']
                    and pr['physical_ids']==ids and pr['symbols']==symbols and pr['physical_charge']==-2
                    and pr['physical_multiplicity']==task['multiplicity'] and pr['target_index']==255,
                    'Original POINT recipe differs')
            if r['MACE']['stored_composite_gradient']:
                require(point['gradient_kcal_mol_A']==r['MACE']['stored_composite_gradient'],'Archived gradient pin differs')
                verify(r['MACE']['stored_composite_gradient'])
        else:
            mr=next(a for a in completed['rows'] if a['cell_id']==tid)
            require(sid=='7CCO' and mr['MACE']==r['MACE']['receipt'] and mr['ALPB']==r['ALPB']['receipt']
                    and mr['force_arrays']==r['MACE']['force_arrays'] and mr['task']==r['MACE']['task']
                    and mr['physical_source']==r['physical'] and mr['positions']==r['positions']
                    and mr['physical_ids']==ids and mr['physical_multiplicity']==task['multiplicity']
                    and recipe['source_POINT_recipe'] is None and recipe['same_recipe_admitted_by']==mp['result'],
                    'Whole-model supplement provenance differs')
        require(nr==ref['collector_row'] and nr['status']=='available' and nr['receipt']==ref['native_receipt']
                and nr['output']==ref['output'] and nr['engrad']==ref['engrad'] and nr['electronic_state']==ref['state']
                and nr['task']==ref['task'] and nr['task']['physical_parent']==r['physical'] and nr['task']['physical_ids']==ids
                and nr['electronic_state']['charge']==-2 and nr['physical_multiplicity']==task['multiplicity']
                and nr['electronic_state']['effective_multiplicity']==1
                and nr['physical_state_qualified'] is False and nr['force_consistency_qualified'] is False,
                'Inherited native reference admission/state differs')
        for k in ('native_receipt','output','engrad'): verify(ref[k])
        sol=c['solvent_admission'][tid]; ar=read(r['ALPB']['receipt'])
        require(sol['status']=='available' and sol['receipt']==r['ALPB']['receipt']
                and sol['energy_hartree']==r['ALPB']['energy_hartree']==ar['total_hartree']
                and sol['components_hartree']==ar['components_hartree']
                and ar['artifacts']==r['ALPB']['artifacts'] and ar['gradient_units']=='Hartree/Bohr; archived only'
                and r['ALPB']['admission_authority']==parents[sid]['collector'],'Admitted solvent differs')
        for v in r['ALPB']['artifacts'].values(): verify(v)
        sg=vectors(np,sol['gradient_hartree_per_bohr']); sk=vectors(np,sol['gradient_kcal_mol_A'])
        require(np.array_equal(sg,vectors(np,np.loadtxt(r['ALPB']['artifacts']['gradient.dat']['path'])))
                and np.array_equal(sg*HA/BOHR,sk),'Original reference solvent conversion differs')
        report=read(ref['force_report'])
        require(c['force_reports'][tid]==ref['force_report'] and report['numerical_force_accuracy_qualified'] is False
                and [a['physical_id'] for a in report['force_table']]==ids,'Inherited reference force table differs')
        keys=('reference_PBE0_vacuum','reference_PBE0_plus_fixed_ALPB') if sid=='7CCO' else ('whole_PBE0_vacuum','whole_PBE0_plus_fixed_ALPB')
        rf={kind:vectors(np,[a['force_kcal_mol_A'][key] for a in report['force_table']])
            for kind,key in zip(('vacuum','composite'),keys)}
        ng=vectors(np,nr['gradient_kcal_mol_per_A'])
        require(np.array_equal(-ng,rf['vacuum']) and np.array_equal(-(ng+sk),rf['composite']),
                'Exact admitted reference-force representation differs')
        require(directions[sid]['physical_ids']==ids,'Frozen direction atom IDs differ')
        ds={k:vectors(np,v) for k,v in directions[sid]['directions'].items()}
        cell=dict(meta=r,raw=raw,task=task,x=x,symbols=symbols,ids=ids,atoms=atoms,solvent=sol,
                  solvent_gradient=sg,reference=rf,directions=ds,native_row=nr,archived_report=report,
                  force_receipt=r['MACE']['receipt'],force_task=r['MACE']['task'])
        if cid in MISSING:
            require(raw['forces_eV_A'] is None and r['MACE']['original_compute_force'] is False,'Force already available')
        else:
            require(raw['force_convention']=='negative derivative of the same native vacuum energy','Cached force sign differs')
            force=vectors(np,raw['forces_eV_A']); cell['mace']=model_forces(np,cell,force)
            key='whole_MACE_composite' if sid=='7CCO' else 'whole_MACE_plus_fixed_ALPB'
            require(np.array_equal(cell['mace']['composite'],vectors(np,[a['force_kcal_mol_A'][key] for a in report['force_table']])),
                    'Original cached composite-force representation differs')
            if r['MACE']['force_arrays']:
                fa=read(r['MACE']['force_arrays'])
                require(fa['physical_ids']==ids and np.array_equal(cell['mace']['vacuum'],vectors(np,fa['components']['MACE_vacuum']))
                        and np.array_equal(cell['mace']['composite'],vectors(np,fa['components']['composite'])),
                        'Archived whole-model force arrays differ')
        cells[cid]=cell
    for sid in ('7CCO','7CCN'):
        for geom in ('before','after'):
            a,b=(cells[sid+'_'+geom+'_'+m] for m in ('La','Dy'))
            require(a['ids']==b['ids'] and np.array_equal(a['x'],b['x']) and a['symbols'][:-1]==b['symbols'][:-1],
                    'Matched metal coordinates differ')
    return dict(np=np,cells=cells,collectors=collectors,mcfg=mcfg,native=parents['native_MACE_helper'],inventory=inv)


def model_forces(np,cell,force):
    vacuum=force*EV; sg=cell['solvent_gradient']
    # Preserve the two archived producers' associations, without equating their bit patterns.
    composite=vacuum+(-sg*(HA/BOHR)) if cell['meta']['source_id']=='7CCO' else -(-force*EV+sg*(HA/BOHR))
    return {'vacuum':vacuum,'composite':composite}


def describe(np,force,directions):
    force=vectors(np,force); norms=np.linalg.norm(force,axis=1)
    return dict(maximum_atom_kcal_mol_A=float(norms.max()),RMS_atom_vector_kcal_mol_A=float(np.sqrt(np.mean(norms**2))),
                projections_kcal_mol_A={k:float(np.sum(force*d)) for k,d in directions.items()})


def comparison(w,cell):
    np=w['np']; model=cell['mace']; ref=cell['reference']; errors={k:model[k]-ref[k] for k in model}
    return dict(status='available',cell_id=cell['meta']['cell_id'],units='kcal/mol/Angstrom',quantity='force=-gradient',
        physical_ids=cell['ids'],source=cell['meta'],reference_native_evidence=cell['native_row'],
        MACE_force_receipt=cell['force_receipt'],MACE_force_task=cell['force_task'],
        force_table=[dict(index=i,physical_id=a['id'],resname=a.get('resname'),
            MACE={k:v[i].tolist() for k,v in model.items()},PBE0={k:v[i].tolist() for k,v in ref.items()},
            MACE_minus_PBE0={k:v[i].tolist() for k,v in errors.items()}) for i,a in enumerate(cell['atoms'])],
        metrics={family:{k:describe(np,v,cell['directions']) for k,v in values.items()}
                 for family,values in [('MACE',model),('PBE0',ref),('MACE_minus_PBE0',errors)]},
        force_accuracy_qualified=False,shared_solvent_independently_validated=False)


def matched(w,sid,geom,available):
    ids=[sid+'_'+geom+'_'+m for m in ('La','Dy')]
    if not all(k in available for k in ids): return dict(source_id=sid,geometry=geom,status='unavailable',required_cells=ids)
    a,b=(w['cells'][k] for k in ids); np=w['np']
    m={k:b['mace'][k]-a['mace'][k] for k in ('vacuum','composite')}
    r={k:b['reference'][k]-a['reference'][k] for k in m}; e={k:m[k]-r[k] for k in m}
    return dict(source_id=sid,geometry=geom,status='available',definition='Dy-minus-La forces; error=MACE-minus-PBE0',
        units='kcal/mol/Angstrom',physical_ids=a['ids'],
        force_vectors={'MACE':{k:v.tolist() for k,v in m.items()},'PBE0':{k:v.tolist() for k,v in r.items()},
                       'MACE_minus_PBE0':{k:v.tolist() for k,v in e.items()}},
        metrics={family:{k:describe(np,v,a['directions']) for k,v in values.items()}
                 for family,values in [('MACE',m),('PBE0',r),('MACE_minus_PBE0',e)]},
        force_accuracy_qualified=False,affinity_or_preference_claim=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--config',required=True); parser.add_argument('--output',required=True)
    args=parser.parse_args(); require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=False); tick=time.monotonic(); native=None
    result=dict(protocol=PROTOCOL,status='unavailable',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),declared_cells=list(IDS),declared_missing=list(MISSING),
        rows=[dict(cell_id=k,status='unavailable',reason='not_admitted') for k in IDS],
        matched_metal_differences=[dict(source_id=s,geometry=g,status='unavailable',reason='not_admitted',
            required_cells=[s+'_'+g+'_'+m for m in ('La','Dy')]) for s in ('7CCO','7CCN') for g in ('before','after')],
        new_MACE_calls_attempted=0,model_loads_attempted=0,model_initialized=False,
        new_DFT_calls=0,new_ALPB_calls=0,new_GFN_calls=0,new_MM_calls=0,optimizer_steps=0,new_geometry_calls=0)
    try:
        cfg=read(result['config']); w=prepare(cfg); available=set()
        result.update(inventory=cfg['inventory'],science_decision=cfg['science_decision'],
            saved_input_admission='all_eight_passed_before_model_load',
            inherited_energy_comparisons={sid:c['comparison'] for sid,c in w['collectors'].items()},
            inherited_reference_caveats={sid:c['caveats'] for sid,c in w['collectors'].items()})
        write(out/'START.json',result); slots={r['cell_id']:r for r in result['rows']}
        for cid in IDS:
            if cid not in MISSING:
                record=write(out/(cid+'_FORCES.json'),comparison(w,w['cells'][cid]))
                slots[cid].update(status='available',reason='reused_saved_forces',force_report=record); available.add(cid)
        stop=None
        for cid in MISSING:
            slot=slots[cid]; cell=w['cells'][cid]
            if stop:
                slot.update(reason='not_attempted_after_failure',preceding_failure=stop); continue
            directory=out/cid; directory.mkdir(); slot.update(reason='attempting',output_directory=str(directory))
            try:
                if native is None:
                    result['model_loads_attempted']+=1; t=time.monotonic()
                    native=module(w['native']).NativeMACE(w['mcfg'])
                    result.update(model_initialized=True,model_load_seconds=time.monotonic()-t,GPU_device=native.device)
                require(native.calls<3,'Three-call ceiling')
                task=dict(cell['task'],task_id=PROTOCOL+'__'+cid,original_task=cell['meta']['MACE']['task'])
                slot['supplement_task']=write(directory/'TASK.json',task)
                before=native.calls
                try: raw=native.evaluate(task,cell['symbols'],cell['x'].tolist(),compute_force=True)
                finally: result['new_MACE_calls_attempted']=native.calls
                require(native.calls-before==1,'One native call required')
                slot['MACE']=write(directory/'MACE.json',raw); state(raw,cell['meta']['metal'])
                delta=raw['energy_eV']-cell['raw']['energy_eV']; slot['energy_replay_delta_eV']=delta
                require(raw['input_state_check']==cell['raw']['input_state_check'] and abs(delta)<=ENERGY_TOL
                        and raw['force_convention']=='negative derivative of the same native vacuum energy','Original energy/state/force replay differs')
                cell['mace']=model_forces(w['np'],cell,vectors(w['np'],raw['forces_eV_A']))
                cell.update(force_receipt=slot['MACE'],force_task=slot['supplement_task'])
                slot.update(status='available',reason='supplement_complete',
                    force_report=write(out/(cid+'_FORCES.json'),comparison(w,cell)))
                available.add(cid)
            except Exception as exc:
                stop=repr(exc); slot.update(status='unavailable',reason='supplement_failed',failure=stop,traceback=traceback.format_exc())
        result['matched_metal_differences']=[matched(w,s,g,available) for s in ('7CCO','7CCN') for g in ('before','after')]
        result['status']='complete' if len(available)==8 else 'unavailable'; verify(result['config'])
    except Exception as exc: result.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick,force_accuracy_qualified=False,no_automatic_retry=True,
        scientific_accuracy_claim=False,energy_replay_tolerance_eV=ENERGY_TOL,
        interpretation='Consumed fixed deformations; physical Dy6 MACE versus frozen-f effective1 reference. Shared ALPB is not independently validated.')
    write(out/'RESULT.json',result); return 0 if result['status']=='complete' else 1


if __name__=='__main__': raise SystemExit(main())
