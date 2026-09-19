"""Independent source-water preparation check; original-core DFT and site vector."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_hybrid import accepted_attempt,check_atoms
from hydration_square import checked_parent,endpoint

STAGE='water_reference_orientation'
PROTOCOL='independent_reference_contextual_water_H_transfer_v1'


def prepare(config,output):
    import mace_omol as omol
    from hydration_network import discover,materialize,radial_seed,write_xyz
    cfg=read_json(config);qualified=read_json(verify(cfg['qualified_preparation_manifest']))
    old=read_json(verify(cfg['baseline_manifest']));parents=[]
    for site in cfg['ordered_sites']:
        pin=cfg['core_preparations'][site];p=read_json(verify(pin));energies={}
        for metal in ('Ca','La'):
            e=p['outputs'][metal]
            candidates=[t for t in old['tasks'] if t['xyz']['sha256']==e['xyz']['sha256'] and t['input']['sha256']==e['input']['sha256']]
            if len(candidates)!=1:raise InvalidArtifact('unique compatible archived endpoint required')
            t=candidates[0];energies[metal]=endpoint(record(t['output_path']),record(t['output_path']+'.execution.json'),e['xyz'],e['input'])
        item={'case':cfg['case_id']+'_'+site,'site':site,'preparation':pin,'endpoints':energies}
        checked_parent(item);parents.append(item)
    _,out,m=omol.common(verify(qualified['inventory']),verify(qualified['software']),verify(cfg['agreement']),output,STAGE)
    wet=[p for p in parents if read_json(verify(p['preparation']))['explicit_water_inventory']]
    source={'parents':wet,'water_reference_geometry':qualified['settings']['water_reference_geometry']}
    write_new(out/'source_config.json',source)
    context_config={'source_config':record(out/'source_config.json'),'contact_cutoff_A':qualified['settings']['contact_cutoff_A'],
        'orientation_seeds':qualified['settings']['starts'],'agreement':cfg['agreement']}
    write_new(out/'context_config.json',context_config)
    _,_,data,union=discover(out/'context_config.json');tasks=[];contexts={}
    for d in data:
        atoms,mapping,ledger,waters=materialize(d,union,xyz(verify(source['water_reference_geometry'])))
        name=d['item']['case'];context={'case':name,'parent':d['item']['preparation'],'atom_graph':mapping,
            'water_groups':waters,'charge_ledger':ledger,'contacts':d['contacts'],'protein_fragments_union':union}
        cp=out/'contexts'/name/'context.json';write_new(cp,context);contexts[name]=record(cp)
        q=sum(a['formal_charge'] for a in ledger)
        for metal in ('Ca','La'):
            starts={}
            for seed,rows in [('source',atoms),('radial_away',radial_seed(atoms,waters,mapping))]:
                xp=cp.parent/(metal+'_'+seed+'.xyz');write_xyz(xp,[(metal,*rows[0][1:])]+rows[1:],PROTOCOL);starts[seed]=record(xp)
            charge=q+(2 if metal=='Ca' else 3)
            tasks.append({'task_id':name+'__'+metal,'case_id':name,'metal':metal,'metal_index':0,'kind':'core','variant':'primary',
                'energy_component':omol.COMPONENT,'energy_only':False,'charge':charge,'spin_multiplicity':1,'xyz':starts['source'],
                'state':check_atoms(xyz(verify(starts['source'])),charge),'water_groups':waters,
                'mobile_indices':sorted(i for w in waters if w['role']=='variable' for i in w['hydrogen_indices']),
                'starts':starts,'water_orientation_optimization':True,'context':record(cp)})
    # Copy the demonstrated optimizer exactly; only the case/runner adapter is new.
    for name in ('hydration_proposal_opt.py','hydration_mace.py','contextual_water_context.py'):
        dest=out/'implementation'/name;shutil.copyfile(verify(qualified['implementation'][name]),dest);m['implementation'][name]=record(dest)
    for name in ('water_reference_transfer.py','hydration_square.py','hydration_network.py','hydration_core_transfer.py','affordable_peptide.py','carve_generic.py','coordination_policy.py'):
        dest=out/'implementation'/name;shutil.copyfile(Path(__file__).with_name(name),dest);m['implementation'][name]=record(dest)
    dispatch=out/'implementation/mace_omol.py';text=dispatch.read_text()
    for op in ('validate','collect'):
        header=f'def {op}(manifest):\n'
        if text.count(header)!=1:raise InvalidArtifact('unsupported existing dispatch')
        text=text.replace(header,header+f"    if read_json(manifest).get('stage')=='{STAGE}':\n        from water_reference_transfer import {op} as operation\n        return operation(manifest)\n")
    header='def worker(manifest, task_id, output, memory_mode):\n'
    text=text.replace(header,header+f"    if read_json(manifest).get('stage')=='{STAGE}':\n        from hydration_proposal_opt import worker as operation\n        return operation(manifest,task_id,output,memory_mode)\n")
    dispatch.write_text(text);m['implementation']['mace_omol.py']=record(dispatch)
    m.update(scientific_protocol_id=PROTOCOL,config=record(config),parents=parents,contexts=contexts,
        optimization=qualified['optimization'],tasks=tasks,baseline_changed=False)
    for t in tasks:t['cache_key']=cache_key({'task':t,'model':m['model'],'software':m['software'],'implementation':m['implementation']})
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


def validate(manifest):
    from mace_omol import model
    m=read_json(manifest);cfg=read_json(verify(m['config']));qualified=read_json(verify(cfg['qualified_preparation_manifest']))
    if m['stage']!=STAGE or m['model']!=model(verify(m['software'])) or m['optimization']!=qualified['optimization']:
        raise InvalidArtifact('proposal method changed')
    verify(m['agreement'])
    for pin in m['implementation'].values():verify(pin)
    for t in m['tasks']:
        original=xyz(verify(t['xyz']));fixed=set(range(len(original)))-set(t['mobile_indices'])
        for pin in t['starts'].values():
            rows=xyz(verify(pin))
            if check_atoms(rows,t['charge'])!=t['state'] or any(rows[i]!=original[i] for i in fixed):raise InvalidArtifact('water proposal altered fixed atoms/state')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):raise InvalidArtifact('cache identity changed')
    expected={p['case']+'__'+metal for p in m['parents'] if read_json(verify(p['preparation']))['explicit_water_inventory'] for metal in ('Ca','La')}
    if len(m['tasks'])!=len(expected) or {t['task_id'] for t in m['tasks']}!=expected:raise InvalidArtifact('wet-site paired denominator changed')
    return {'status':'pass','tasks':len(m['tasks']),'ordered_sites':cfg['ordered_sites'],'manifest':record(manifest)}


def collect(manifest):
    from contextual_water_context import verify_proposal_geometry
    validate(manifest);m=read_json(manifest);mp=Path(manifest);rows=[]
    for t in m['tasks']:
        valid=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:valid.append((a,r))
        row={'task_id':t['task_id'],'case_id':t['case_id'],'metal':t['metal'],'status':'unavailable'}
        if valid:
            a,r=valid[-1]
            try:
                verify_proposal_geometry(t,r)
                row.update(status='prepared',result=record(a/'result.json'),proposed_xyz=r['proposed_xyz'],selected_seed=r['selected_seed'])
            except InvalidArtifact as e:row['reason']=str(e)
        rows.append(row)
    return {'status':'complete' if all(r['status']=='prepared' for r in rows) else 'incomplete','manifest':record(manifest),'rows':rows}


def prepare_dft(collection,output):
    from hydration_core_transfer import transfer_water_hydrogens
    from hydration_network import write_xyz
    saved=read_json(collection);source=verify(saved['manifest']);m=read_json(source)
    if collect(source)!=saved or saved['status']!='complete':raise InvalidArtifact('actual qualified proposals required')
    cfg=read_json(verify(m['config']));baseline=read_json(verify(cfg['baseline_manifest']));out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copytree(source.parent/'implementation',out/'implementation');tasks=[];identities=[]
    for parent in m['parents']:
        p=read_json(verify(parent['preparation']))
        if not p['explicit_water_inventory']:identities.append(parent);continue
        context=read_json(verify(m['contexts'][parent['case']]))
        for metal in ('Ca','La'):
            r=next(r for r in saved['rows'] if r['case_id']==parent['case'] and r['metal']==metal)
            original=xyz(verify(p['outputs'][metal]['xyz']));rows,mobile=transfer_water_hydrogens(p,original,context,xyz(verify(r['proposed_xyz'])))
            td=out/'tasks'/(parent['case']+'__'+metal);td.mkdir(parents=True);xp=td/'core.xyz';write_xyz(xp,rows,PROTOCOL)
            old_input=verify(p['outputs'][metal]['input']).read_text();body=old_input.split('* xyzfile')[0]
            ip=td/'endpoint.inp';ip.write_text(body+f"* xyzfile {p['outputs'][metal]['charge']} 1 core.xyz\n")
            t={'task_id':parent['case']+'__'+metal,'case':parent['case'],'metal':metal,'charge':p['outputs'][metal]['charge'],'multiplicity':1,
               'input':record(ip),'xyz':record(xp),'output_path':str(td/'endpoint.out'),'baseline':parent['endpoints'][metal],
               'parent':parent['preparation'],'mobile_H_indices':mobile,'proposal':r['result']}
            t['cache_key']=cache_key(t|{'scientific_protocol_id':PROTOCOL});tasks.append(t)
    result={'protocol_id':PROTOCOL,'tasks':tasks,'identities':identities,'ordered_sites':cfg['ordered_sites'],'parents':m['parents'],
        'source_collection':record(collection),'agreement':m['agreement'],'orca':baseline['orca'],'execution_policy':baseline['execution_policy'],
        'implementation':{p.name:record(p) for p in (out/'implementation').glob('*.py')},
        'execution_resources':{'mpi_ranks':16,'concurrent_tasks':2},'baseline_changed':False}
    write_new(out/'manifest.json',result)
    from affordable_workflow import dry_run
    return dry_run(out/'manifest.json')


def collect_dft(manifest,output):
    m=read_json(manifest);rows=[]
    for parent in m['parents']:
        p=read_json(verify(parent['preparation']));energies={};failures=[]
        for metal in ('Ca','La'):
            if not p['explicit_water_inventory']:energies[metal]=parent['endpoints'][metal];continue
            t=next(t for t in m['tasks'] if t['case']==parent['case'] and t['metal']==metal)
            try:
                original=xyz(verify(p['outputs'][metal]['xyz']));new=xyz(verify(t['xyz']))
                if len(original)!=len(new) or any(original[i]!=new[i] for i in range(len(new)) if i not in t['mobile_H_indices']):raise InvalidArtifact('nonwater-H coordinates changed')
                energies[metal]=endpoint(record(t['output_path']),record(t['output_path']+'.execution.json'),t['xyz'],t['input'])
            except (OSError,ValueError,KeyError) as exc:failures.append({'metal':metal,'reason':str(exc)})
        oldR=(parent['endpoints']['Ca']['energy_hartree']-parent['endpoints']['La']['energy_hartree'])*HA_TO_KCAL
        newR=(energies['Ca']['energy_hartree']-energies['La']['energy_hartree'])*HA_TO_KCAL if len(energies)==2 else None
        rows.append({'site':parent['site'],'case':parent['case'],'status':'complete' if not failures else 'unavailable','failures':failures,
            'water_count':len(p['explicit_water_inventory']),'operation':'context_water_H' if p['explicit_water_inventory'] else 'exact_dry_identity',
            'baseline_R_kcal_mol':oldR,'prepared_R_kcal_mol':newR,'delta_R_kcal_mol':newR-oldR if newR is not None else None,
            'endpoints':energies,'baseline_endpoints':parent['endpoints'],'calibrated_decision':None,'occupancy_probability':None})
    result={'manifest':record(manifest),'ordered_sites':m['ordered_sites'],'rows':rows,'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete',
            'evidence_stratum':'supporting_cross_study_La_Ca_comparison','baseline_changed':False}
    write_new(output,result);return {'status':result['status'],'sites':len(rows)}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='op',required=True)
    for op,fields in [('prepare',('config','output')),('validate',('manifest',)),('collect',('manifest','output')),('prepare-dft',('collection','output')),('collect-dft',('manifest','output'))]:
        q=sub.add_parser(op)
        for name in fields:q.add_argument('--'+name,required=True)
    args=vars(parser.parse_args());op=args.pop('op')
    if op=='collect':dest=args.pop('output');r=collect(**args);write_new(dest,r)
    else:r=globals()[op.replace('-','_')](**args)
    print(json.dumps(r,indent=2))
