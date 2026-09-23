"""Dy-conditioned Hans transfer, reusing the exact previously computed Mex sites."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import time
import math

import lanm_series_followup as original
from affordable_common import InvalidArtifact, HA_TO_KCAL, read_json, record, verify, write_new, xyz
from compact_solvation import completed, diagnostics
from mace_hybrid import EV_TO_KCAL

PROTOCOL = 'LanM_Hans_Dy_crystal_transfer_native_OMOL_fcore_GFN2_v1'


def atom_id(s):
    return [s[k] for k in ('chain','resnum','insertion_code','resname','atom','element')]


def physical_map(repair):
    result = [['selected_metal']]
    for atom in sorted(repair['atom_graph']['source_to_qm'], key=lambda x:x['qm_index']):
        if atom['kind'] == 'source':
            result.append(['source', *atom_id(atom['source'])])
        else:
            result.append(['sigma_link_H', atom_id(atom['retained']), atom_id(atom['omitted']), atom['length_A']])
    return result


def prepare(root, agreement, original_manifest, original_collection, output):
    from protonate_cif import protonate
    from affordable_peptide import repair
    root = Path(root).resolve(); out = Path(output).resolve(); out.mkdir(parents=True,exist_ok=False)
    prior = read_json(original_manifest); old_results = read_json(original_collection)
    if old_results['manifest'] != record(original_manifest) or old_results['complete_endpoints'] != 12:
        raise InvalidArtifact('complete compatible original pilot required')
    source = root/'hans_lanm_protenix_benchmark/references/8FNR.pdb'
    selected = out/'hans_Dy_selected_source.pdb'
    sel = original.select_source_conformers(source,selected)
    pp = out/'hans_Dy_protonated.pdb'
    protonation = protonate(selected,pp,ph=5.,add_missing_residues=False)
    if original.heavy_atoms(source) != original.heavy_atoms(pp):
        raise InvalidArtifact('deposited Hans Dy heavy atom inventory/coordinates changed')
    topology = verify(read_json(verify(prior['sources']))['topology'])
    carver = original.private_carver(); carver.METAL_ELEMENTS |= {'DY'}
    cases=[]
    for index,site in enumerate(original.SITES):
        row={'case_id':'Hans_'+site,'site':site,'status':'unavailable','pdb':'8FNR',
             'chain':'A','metal_resnum':201+index,'conditioning_metal':'Dy'}
        try:
            stem='hans_Dy_'+site;d=out/'sources'/site
            carver.carve(pp,d/'selection_intermediate',stem,site_chain='A',site_resnum=201+index,
                site_icode='',site_atom='DY',qm_inclusion_cut=3.3)
            r=repair(d/'selection_intermediate'/(stem+'_carve_manifest.json'),d/'amide_v3',topology)
            oldcase=next(c for c in prior['cases'] if c['case_id']==row['case_id'])
            oldr=read_json(verify(oldcase['repair_manifest']))
            mapping=physical_map(r); oldmapping=physical_map(oldr)
            match=(mapping==oldmapping and r['outputs']['La']['charge']==oldr['outputs']['La']['charge']
                   and [a[0] for a in xyz(verify(r['outputs']['La']['xyz']))]==[a[0] for a in xyz(verify(oldr['outputs']['La']['xyz']))])
            row.update(status='prepared',repair_manifest=record(d/'amide_v3/repair_manifest.json'),
                original_repair=oldcase['repair_manifest'],physical_mapping=mapping,
                original_physical_mapping=oldmapping,mapping_and_composition_match=match,
                explicit_waters=r['explicit_water_inventory'],atom_count=r['paired_invariants']['atom_count'])
        except Exception as exc:
            row['reason']=f'{type(exc).__name__}: {exc}'
        cases.append(row)
    pins={};impl=out/'implementation';impl.mkdir()
    for p in Path(__file__).parent.glob('*.py'):
        target=impl/p.name;shutil.copyfile(p,target);pins[p.name]=record(target)
    initial=out/'native_initial';initial.mkdir();tasks=[];endpoints=[]
    params=read_json(verify(prior['parameter_export']))
    for case in cases:
        for metal in original.METALS:
            e={'endpoint_id':case['case_id']+'__'+metal,'case_id':case['case_id'],'site':case['site'],
               'metal':metal,'status':'unavailable'}
            endpoints.append(e)
            if case['status']!='prepared':e['reason']=case['reason'];continue
            r=read_json(verify(case['repair_manifest']));old=r['outputs']['La'];d=out/'endpoints'/e['endpoint_id'];d.mkdir(parents=True)
            lines=verify(old['xyz']).read_text().splitlines();lines[2]=metal+lines[2][2:]
            (d/'core.xyz').write_text('\n'.join(lines)+'\n');atoms=xyz(d/'core.xyz')
            st=original.state(atoms,old['charge'],metal,params)
            e.update(status='prepared',xyz=record(d/'core.xyz'),charge=old['charge'],state=st,
                     spin_multiplicity=original.PHYSICAL_MULTIPLICITY[metal])
            for medium in ('vacuum','alpb'):
                tid=e['endpoint_id']+'__'+medium;dd=initial/'tasks'/tid;dd.mkdir(parents=True)
                shutil.copyfile(d/'core.xyz',dd/'core.xyz');(dd/'endpoint.inp').write_text(original.recipe(e['charge'],medium))
                tasks.append({'task_id':tid,'endpoint_id':e['endpoint_id'],'case_id':case['case_id'],
                    'metal':metal,'medium':medium,'charge':e['charge'],'multiplicity':1,'state':st,
                    'xyz':record(dd/'core.xyz'),'input':record(dd/'endpoint.inp'),'output_path':str(dd/'endpoint.out')})
    low_prior=read_json(verify(prior['native_initial_manifest']))
    low={'protocol_id':original.PROTOCOL,'experiment_protocol_id':PROTOCOL,'stage':'initial','tasks':tasks,
         'orca':low_prior['orca'],'parameter_export':prior['parameter_export'],'implementation':pins,
         'execution_resources':{'mpi_ranks':8,'concurrent_tasks':8},
         'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']}}
    write_new(initial/'manifest.json',low)
    m={'protocol_id':PROTOCOL,'agreement':record(agreement),'original_manifest':record(original_manifest),
       'original_collection':record(original_collection),'source':record(source),'selected_source':sel,
       'protonation':protonation,'cases':cases,'endpoints':endpoints,'implementation':pins,
       'model':prior['model'],'parameter_export':prior['parameter_export'],
       'native_initial_manifest':record(initial/'manifest.json'),
       'finite_new_calls':{'MACE':6,'GFN_initial':12,'GFN_self_continuation':12},
       'reused_Mex_endpoints':6,'baseline_changed':False}
    write_new(out/'manifest.json',m);result=validate(out/'manifest.json');write_new(out/'PREFLIGHT.json',result);return result


def validate(manifest):
    from run_orca_task_manifest import load_manifest_tasks
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or len(m['endpoints'])!=6 or [c['site'] for c in m['cases']]!=list(original.SITES):
        raise InvalidArtifact('fixed three-site transfer scope differs')
    for k in ('agreement','original_manifest','original_collection','source','parameter_export'):verify(m[k])
    for p in m['implementation'].values():verify(p)
    prior=read_json(verify(m['original_manifest']))
    if m['model']!=prior['model'] or m['parameter_export']!=prior['parameter_export']:
        raise InvalidArtifact('model differs from reused Mex states')
    eps={e['endpoint_id']:e for e in m['endpoints']}
    for c in m['cases']:
        a,b=[eps[c['case_id']+'__'+z] for z in original.METALS]
        if c['status']!='prepared':continue
        aa,bb=xyz(verify(a['xyz'])),xyz(verify(b['xyz']))
        if aa[1:]!=bb[1:] or aa[0][1:]!=bb[0][1:] or a['charge']!=b['charge']:
            raise InvalidArtifact('La/Dy paired geometry/charge differs')
        r=read_json(verify(c['repair_manifest']))
        if physical_map(r)!=c['physical_mapping']:raise InvalidArtifact('physical mapping differs')
    low=read_json(verify(m['native_initial_manifest']))
    if low['tasks']:load_manifest_tasks(verify(m['native_initial_manifest']))
    for t in low['tasks']:
        e=eps[t['endpoint_id']]
        if verify(t['xyz']).read_bytes()!=verify(e['xyz']).read_bytes() or verify(t['input']).read_text()!=original.recipe(e['charge'],t['medium']):
            raise InvalidArtifact('native geometry/recipe differs')
    return {'status':'prepared','manifest':record(manifest),'supported_sites':sum(c['status']=='prepared' for c in m['cases']),
        'site_denominator':3,'MACE_tasks':sum(e['status']=='prepared' for e in m['endpoints']),
        'GFN_initial_tasks':len(low['tasks']),'mapping_matches':[c.get('mapping_and_composition_match') for c in m['cases']],
        'new_molecular_calls':0}


def mace(manifest,output):
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from mace_omol import input_batch
    validate(manifest);m=read_json(manifest)
    if not os.environ.get('SLURM_JOB_ID') or not torch.cuda.is_available():raise InvalidArtifact('allocated GPU required')
    if record(__file__)['sha256']!=m['implementation']['lanm_structural_transfer.py']['sha256']:raise InvalidArtifact('use pinned adapter')
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK']));torch.set_default_dtype(torch.float64);torch.cuda.reset_peak_memory_stats()
    tick=time.monotonic();calc=mace_omol(model=str(verify(m['model']['checkpoint'])),device='cuda',default_dtype='float64');model=calc.models[0]
    if type(model).__name__!='ScaleShiftMACE' or dict(model.embedding_specs)!=m['model']['embedding_specs'] or calc.energy_units_to_eV!=1.:
        raise InvalidArtifact('model/units differ')
    for p in model.parameters():p.requires_grad_(False)
    load=time.monotonic()-tick;rows=[]
    for e in m['endpoints']:
        row={'endpoint_id':e['endpoint_id'],'status':'unavailable','energy_eV':None,'model_call_attempted':False}
        if e['status']=='prepared':
            start=time.monotonic()
            try:
                atoms=xyz(verify(e['xyz']));a=Atoms([x[0] for x in atoms],positions=[x[1:] for x in atoms],pbc=False)
                a.info.update(charge=e['charge'],spin=e['spin_multiplicity']);batch=calc._atoms_to_batch(a)
                row['input_state_check']=input_batch(calc,a,e['charge'],e['spin_multiplicity'],batch,0)
                with torch.no_grad():
                    row['model_call_attempted']=True
                    value=model(batch.to_dict(),training=False,compute_force=False,compute_virials=False,compute_stress=False,
                        compute_displacement=False,compute_hessian=False,compute_edge_forces=False,compute_atomic_stresses=False)
                energy=float(value['energy'].detach().cpu().item())
                if not math.isfinite(energy) or value['forces'] is not None:raise InvalidArtifact('invalid native scalar output')
                row.update(status='complete',energy_eV=energy,xyz=e['xyz'],state=e['state'])
            except Exception as exc:row.update(status='failed',reason=f'{type(exc).__name__}: {exc}')
            row['wall_seconds']=time.monotonic()-start
        rows.append(row);write_new(out/(e['endpoint_id']+'.json'),row)
    result={'manifest':record(manifest),'rows':rows,'model_load_seconds':load,'wall_seconds':time.monotonic()-tick,
        'job_id':os.environ['SLURM_JOB_ID'],'allocated_CPUs':os.environ['SLURM_CPUS_PER_TASK'],'device':torch.cuda.get_device_name(),
        'peak_host_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'peak_cuda_allocated_bytes':torch.cuda.max_memory_allocated(),
        'actual_MACE_calls':sum(r['model_call_attempted'] for r in rows)}
    write_new(out/'result.json',result);return result


def prepare_continuation(manifest,output):
    validate(manifest);m=read_json(manifest);initial=verify(m['native_initial_manifest']);low=read_json(initial)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];missing=[]
    for t in low['tasks']:
        pin=completed(initial,t['task_id'])
        if pin is None:missing.append({'task_id':t['task_id'],'reason':'initial failed/unrun'});continue
        audit=diagnostics(pin,t)
        if audit['charge_sanity_status']!='pass' or audit['native_valence_electrons']!=t['state']['GFN_valence_electron_count']:
            missing.append({'task_id':t['task_id'],'reason':'initial native state/charge unsupported'});continue
        src=verify(pin['output']).parent
        if not all((src/('endpoint.runtime.'+s)).is_file() for s in ('gbw','xtbw')):
            missing.append({'task_id':t['task_id'],'reason':'matching native seed files unavailable'});continue
        d=out/'tasks'/t['task_id'];d.mkdir(parents=True);shutil.copyfile(verify(t['xyz']),d/'core.xyz')
        (d/'endpoint.inp').write_text(original.recipe(t['charge'],t['medium'],True))
        for s in ('gbw','xtbw'):shutil.copyfile(src/('endpoint.runtime.'+s),d/('seed.immutable.'+s))
        tasks.append({**t,'xyz':record(d/'core.xyz'),'input':record(d/'endpoint.inp'),'output_path':str(d/'endpoint.out'),
            'original':pin,'initial_audit':audit,'immutable_seed':record(d/'seed.immutable.xtbw'),'immutable_gbw':record(d/'seed.immutable.gbw')})
    cm={**low,'parent':record(manifest),'stage':'self_continuation','tasks':tasks,'unavailable':missing}
    write_new(out/'manifest.json',cm)
    from run_orca_task_manifest import load_manifest_tasks
    if tasks:load_manifest_tasks(out/'manifest.json')
    return {'manifest':record(out/'manifest.json'),'prepared':len(tasks),'unavailable':missing,'new_calls':0}


def collect(manifest,mace_collection,continuation_manifest,output):
    validate(manifest);m=read_json(manifest);mr=read_json(mace_collection);cm=read_json(continuation_manifest)
    if mr['manifest']!=record(manifest) or cm['parent']!=record(manifest):raise InvalidArtifact('collection parent differs')
    mc={r['endpoint_id']:r for r in mr['rows']};ct={t['task_id']:t for t in cm['tasks']};rows=[]
    old=read_json(verify(m['original_collection']));oldrows={r['endpoint_id']:r for r in old['endpoints']}
    for e in m['endpoints']:
        r={'endpoint_id':e['endpoint_id'],'status':'unavailable','native_MACE_eV':None,
           'solvent_transfer_kcal_mol':None,'composite_kcal_mol':None,'cells':{}}
        native=mc[e['endpoint_id']];r['MACE_receipt']=native
        if native['status']=='complete':
            if native['xyz']!=e['xyz'] or native['state']!=e['state']:raise InvalidArtifact('MACE state/source changed')
            r['native_MACE_eV']=native['energy_eV']
        for medium in ('vacuum','alpb'):
            tid=e['endpoint_id']+'__'+medium;cell={'status':'unavailable','continuation_hartree':None};r['cells'][medium]=cell
            if tid not in ct:continue
            t=ct[tid];pin=completed(Path(continuation_manifest),tid)
            if pin is None:continue
            try:
                audit=diagnostics(pin,t);first=t['original'];text=verify(pin['output']).read_text()
                if 'INITIAL GUESS: XTBRESTART' not in text:raise InvalidArtifact('native restart unconfirmed')
                if audit['parameter_export']['sha256']!=t['initial_audit']['parameter_export']['sha256'] or audit['charge_sanity_status']!='pass':
                    raise InvalidArtifact('continued parameters/charge failed')
                if audit['native_valence_electrons']!=e['state']['GFN_valence_electron_count']:raise InvalidArtifact('effective electrons changed')
                cell.update(status='complete',continuation=pin,initial=first,continuation_hartree=pin['energy_hartree'],
                    initial_hartree=first['energy_hartree'],audit=audit,
                    continuation_minus_initial_kcal_mol=(pin['energy_hartree']-first['energy_hartree'])*HA_TO_KCAL)
            except (InvalidArtifact,KeyError,OSError,ValueError) as exc:cell['reason']=str(exc)
        if r['native_MACE_eV'] is not None and all(c['status']=='complete' for c in r['cells'].values()):
            transfer=(r['cells']['alpb']['continuation_hartree']-r['cells']['vacuum']['continuation_hartree'])*HA_TO_KCAL
            r.update(status='complete',solvent_transfer_kcal_mol=transfer,composite_kcal_mol=r['native_MACE_eV']*EV_TO_KCAL+transfer)
        rows.append(r)
    newrows={r['endpoint_id']:r for r in rows};vector=[]
    for c in m['cases']:
        site=c['site'];h=[newrows['Hans_'+site+'__'+z] for z in original.METALS];mx=[oldrows['Mex_'+site+'__'+z] for z in original.METALS]
        values=h+mx
        comp=original.balanced_difference(*(x['composite_kcal_mol'] for x in values))
        native=original.balanced_difference(*(None if x['native_MACE_eV'] is None else x['native_MACE_eV']*EV_TO_KCAL for x in values))
        solv=original.balanced_difference(*(x['solvent_transfer_kcal_mol'] for x in values))
        prior=next(x for x in old['ordered_vector'] if x['site']==site)
        work={z:None for z in original.METALS}
        if c.get('mapping_and_composition_match'):
            for z in original.METALS:
                a=newrows['Hans_'+site+'__'+z];b=oldrows['Hans_'+site+'__'+z]
                if a['status']=='complete':work[z]={'composite_kcal_mol':a['composite_kcal_mol']-b['composite_kcal_mol'],
                    'MACE_kcal_mol':(a['native_MACE_eV']-b['native_MACE_eV'])*EV_TO_KCAL,
                    'solvent_kcal_mol':a['solvent_transfer_kcal_mol']-b['solvent_transfer_kcal_mol']}
        vector.append({'site':site,'D_Dy_crystal_kcal_mol':comp,'D_native_MACE_kcal_mol':native,
            'D_solvent_transfer_kcal_mol':solv,'D_La_crystal_kcal_mol':prior['D_composite_kcal_mol'],
            'change_in_D_kcal_mol':None if comp is None else comp-prior['D_composite_kcal_mol'],
            'mapping_and_composition_match':c.get('mapping_and_composition_match',False),'Dy_source_minus_La_source_work':work})
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'fresh_endpoints':rows,'reused_Mex':{k:v for k,v in oldrows.items() if k.startswith('Mex_')},
        'reused_collection':m['original_collection'],'ordered_vector':vector,'fresh_complete':sum(r['status']=='complete' for r in rows),
        'fresh_denominator':6,'Kd':None,'populations':None,'pooling_performed':False,'baseline_changed':False}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='operation',required=True)
    q=sub.add_parser('prepare')
    for k in ('root','agreement','original-manifest','original-collection','output'):q.add_argument('--'+k,required=True)
    for op in ('validate','mace','prepare-continuation'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op!='validate':q.add_argument('--output',required=True)
    q=sub.add_parser('execute-native');q.add_argument('--manifest',required=True);q.add_argument('--authorization',required=True)
    q=sub.add_parser('collect')
    for k in ('manifest','mace-collection','continuation-manifest','output'):q.add_argument('--'+k,required=True)
    a=vars(p.parse_args());op=a.pop('operation').replace('-','_')
    if op=='execute_native':
        if read_json(a['manifest'])['experiment_protocol_id']!=PROTOCOL:raise InvalidArtifact('wrong transfer experiment')
        r=original.execute_native(**a)
    else:r=globals()[op](**a)
    print(json.dumps({k:v for k,v in r.items() if k not in ('fresh_endpoints','rows','reused_Mex')},indent=2))


if __name__=='__main__':main()
