"""Matched single-water perturbations using archived real cores and ORCA runner.

This is an electronic hydration diagnostic, not an occupancy/free-energy model.
The preparer never launches calculations or changes the production carver.
"""
from __future__ import annotations

import argparse
from collections import Counter
import datetime as dt
import json
from pathlib import Path
import re
import shutil

import numpy as np

from affordable_common import (HA_TO_KCAL, InvalidArtifact, cache_key, energy,
                               paired, read_json, record, verify, write_new, xyz)

PROTOCOL = 'native_r2scan3c_cpcm_single_water_square_v1'
PARENT_PROTOCOL = 'generic_peptide_amide_vertical_native_r2scan3c_v3'
METHOD = '! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3'
KEYS = ('chain', 'resnum', 'insertion_code', 'resname')


def key(atom):
    return tuple(atom.get(k) for k in KEYS)


def water_groups(parent, atoms):
    mapping = parent['atom_graph']['source_to_qm']
    if {a['qm_index'] for a in mapping} != set(range(1, len(atoms))):
        raise InvalidArtifact('incomplete/duplicate parent atom mapping')
    if len(mapping) != len(atoms)-1:
        raise InvalidArtifact('duplicate parent atom mapping')
    for item in mapping:
        if not np.allclose(item['xyz_A'], atoms[item['qm_index']][1:], rtol=0, atol=1e-8):
            raise InvalidArtifact('parent mapping coordinates disagree')
    groups = []
    for water in parent['explicit_water_inventory']:
        entries = [a for a in mapping if a['kind']=='source' and key(a['source'])==key(water)]
        indices = [a['qm_index'] for a in entries]
        if Counter(atoms[i][0] for i in indices) != Counter({'O':1, 'H':2}):
            raise InvalidArtifact('water is not exactly one complete H2O')
        ledger = [a for a in parent['charge_ledger'] if a['kind']=='water' and key(a['source'])==key(water)]
        if len(ledger)!=1 or ledger[0]['formal_charge']!=0:
            raise InvalidArtifact('water charge ledger is not neutral/unique')
        oxygen = next(i for i in indices if atoms[i][0]=='O')
        hydrogens = [i for i in indices if atoms[i][0]=='H']
        vectors = np.array([atoms[i][1:] for i in hydrogens])-np.array(atoms[oxygen][1:])
        lengths = np.linalg.norm(vectors, axis=1)
        if not np.isfinite(vectors).all() or np.any(lengths < .1):
            raise InvalidArtifact('invalid water coordinates')
        directions = vectors/lengths[:,None]
        if np.linalg.norm(directions.sum(axis=0))<1e-8 or np.linalg.norm(directions[0]-directions[1])<1e-8:
            raise InvalidArtifact('degenerate water plane')
        angle = float(np.degrees(np.arccos(np.clip(directions[0]@directions[1],-1,1))))
        groups.append({'source':water, 'indices':indices, 'oxygen_index':oxygen,
                       'hydrogen_indices':hydrogens, 'OH_lengths_A':lengths.tolist(), 'HOH_angle_deg':angle})
    flat = [i for w in groups for i in w['indices']]
    if len(set(flat))!=len(flat):raise InvalidArtifact('overlapping water groups')
    return groups


def reference_geometry(atoms, groups, reference):
    """Change only water internal geometry; retain each plane and bisector."""
    if Counter(a[0] for a in reference)!=Counter({'O':1,'H':2}):
        raise InvalidArtifact('invalid water geometry reference')
    o = np.array(next(a[1:] for a in reference if a[0]=='O'))
    v = np.array([a[1:] for a in reference if a[0]=='H'])-o
    lengths = np.linalg.norm(v,axis=1)
    if abs(lengths[0]-lengths[1])>1e-6:
        raise InvalidArtifact('reference requires equal OH lengths')
    half_angle = .5*np.arccos(np.clip(v[0]@v[1]/np.prod(lengths),-1,1))
    radial, axial = float(lengths.mean()*np.sin(half_angle)), float(lengths.mean()*np.cos(half_angle))
    result = list(atoms)
    for water in groups:
        origin = np.array(atoms[water['oxygen_index']][1:])
        v = np.array([atoms[i][1:] for i in water['hydrogen_indices']])-origin
        u = v/np.linalg.norm(v,axis=1)[:,None]
        bisector = u[0]+u[1]; bisector/=np.linalg.norm(bisector)
        axis = u[0]-u[1]; axis/=np.linalg.norm(axis)
        for index,sign in zip(water['hydrogen_indices'],(1,-1)):
            result[index] = ('H',*(origin+axial*bisector+sign*radial*axis).tolist())
    return result


def checked_parent(item):
    parent = read_json(verify(item['preparation']))
    if parent['protocol_id']!=PARENT_PROTOCOL:raise InvalidArtifact('unsupported parent protocol')
    outputs = parent['outputs']
    paired(verify(outputs['La']['xyz']),verify(outputs['Ca']['xyz']),outputs['La']['charge'],outputs['Ca']['charge'])
    for metal in ('La','Ca'):
        endpoint = outputs[metal]
        body = verify(endpoint['input']).read_text()
        simple = [line.strip().lower().split()[1:] for line in body.splitlines() if line.strip().startswith('!')]
        if simple != [METHOD.lower().split()[1:]] or re.search(r'%basis|%pointcharges|\b(?:Opt|EnGrad|NumGrad|Freq|TightSCF)\b',body,re.I):
            raise InvalidArtifact('parent method differs from declared recipe')
        if endpoint['multiplicity']!=1 or endpoint['all_electron_count_for_parity']%2:
            raise InvalidArtifact('parent electronic state unsupported')
        archived = item['endpoints'][metal]
        receipt = read_json(verify(archived['receipt']))
        if not receipt['normal_termination'] or not receipt['scf_converged'] or receipt['returncode']!=0:
            raise InvalidArtifact('failed archived endpoint')
        if receipt['artifacts']['xyz']['sha256']!=endpoint['xyz']['sha256']:
            raise InvalidArtifact('archived endpoint has different coordinates')
        if receipt['artifacts']['template_input']['sha256']!=endpoint['input']['sha256'] or receipt['orca_version']!='6.1.1':
            raise InvalidArtifact('archived endpoint has different method/software')
        verify(receipt['artifacts']['xyz']);verify(receipt['artifacts']['template_input'])
        if receipt['artifacts']['output']['sha256']!=archived['output']['sha256']:
            raise InvalidArtifact('archived output receipt mismatch')
        if energy(verify(archived['output']))!=archived['energy_hartree']:
            raise InvalidArtifact('archived energy mismatch')
    atoms = xyz(verify(outputs['La']['xyz']))
    return parent, atoms, water_groups(parent,atoms)


def audit(config, output):
    cfg = read_json(config); rows=[]
    for item in cfg['parents']:
        parent,atoms,groups=checked_parent(item)
        rows.append({'case':item['case'],'preparation':item['preparation'],'water_geometry':groups,
                     'atom_count':len(atoms),'archived_endpoints':item['endpoints']})
    implementation_path=Path(output).with_suffix('.implementation.py')
    with implementation_path.open('x') as handle:handle.write(Path(__file__).read_text())
    result={'config':record(config),'rows':rows,'new_energy_evaluations':0,'implementation':record(implementation_path)}
    write_new(output,result)
    return {'status':'audited','parents':len(rows),'waters':sum(len(r['water_geometry']) for r in rows)}


def prepare(config, output):
    cfg=read_json(config); verify(cfg['agreement']); out=Path(output).resolve()
    out.mkdir(parents=True,exist_ok=False)
    mode=cfg['water_geometry']
    if mode not in ('archived','reference_internal_geometry'):raise InvalidArtifact('unknown geometry policy')
    reference=xyz(verify(cfg['water_reference_geometry']))
    tasks=[]; states=[]; parents=[]
    for item in cfg['parents']:
        parent,atoms,groups=checked_parent(item)
        base=reference_geometry(atoms,groups,reference) if mode=='reference_internal_geometry' else atoms
        parent_record={'case':item['case'],'preparation':item['preparation'],'endpoints':item['endpoints'],
                       'water_geometry':groups,'biological_group':'bovine_alpha_lactalbumin'}
        parents.append(parent_record)
        variants=[('full',None),*[(f"minus_{w['source']['chain']}_{w['source']['resnum']}",w) for w in groups]]
        for name,removed in variants:
            state_id=item['case']+'__'+name
            if removed is None and mode=='archived':
                states.append({'state_id':state_id,'case':item['case'],'removed_water':None,
                               'water_count':len(groups),'reused_endpoints':item['endpoints'],'tasks':[]})
                continue
            remove=set(removed['indices']) if removed else set()
            keep=[i for i in range(len(atoms)) if i not in remove]
            mapping=[{'parent_index':i,'qm_index':j} for j,i in enumerate(keep)]
            d=out/state_id; ids=[]; files={}
            for metal in ('La','Ca'):
                ep=parent['outputs'][metal]; original=verify(ep['xyz']).read_text().splitlines()[2:]
                changed={i for w in groups for i in w['hydrogen_indices']} if mode=='reference_internal_geometry' else set()
                lines=[(f"H {base[i][1]:.10f} {base[i][2]:.10f} {base[i][3]:.10f}" if i in changed else original[i]) for i in keep]
                md=d/metal;md.mkdir(parents=True)
                xp=md/'core.xyz'; xp.write_text(f'{len(lines)}\n{state_id} {metal}; geometry={mode}\n'+'\n'.join(lines)+'\n')
                ip=md/'endpoint.inp'; ip.write_text(f'{METHOD}\n* xyzfile {ep["charge"]} 1 core.xyz\n')
                tid=state_id+'_'+metal;ids.append(tid);files[metal]=xp
                electrons=ep['all_electron_count_for_parity']-(10 if removed else 0)
                if electrons%2:raise InvalidArtifact('odd electron count after neutral water deletion')
                tasks.append({'task_id':tid,'case':state_id,'metal':metal,'charge':ep['charge'],'multiplicity':1,
                              'input':record(ip),'xyz':record(xp),'output_path':str(md/'endpoint.out'),
                              'protocol_id':PROTOCOL,'task_type':'single_point',
                              'all_electron_count_for_parity':electrons})
            inv=paired(files['La'],files['Ca'],parent['outputs']['La']['charge'],parent['outputs']['Ca']['charge'])
            state={'state_id':state_id,'case':item['case'],'removed_water':removed,
                   'water_count':len(groups)-(removed is not None),'tasks':ids,'source_mapping':mapping,
                   'paired_invariants':inv,'water_geometry_policy':mode,'parent_preparation':item['preparation'],
                   'changed_parent_indices':sorted(changed-remove),'retained_parent_indices':keep,
                   'protein_and_heavy_coordinates':'byte_identical_to_parent'}
            state['cache_key']=cache_key(state|{'inputs':[t['xyz'] for t in tasks[-2:]],'method':METHOD})
            write_new(d/'preparation.json',state);states.append(state)
    shutil.copyfile(__file__,out/'implementation.py')
    result={'schema_version':'alquemia.hydration_square.v1','protocol_id':PROTOCOL,
            'configuration':record(config),'agreement':cfg['agreement'],'parents':parents,'states':states,'tasks':tasks,
            'implementation':record(out/'implementation.py'),'execution_policy':cfg['execution_policy'],'orca':cfg['orca'],
            'water_reference_geometry':cfg['water_reference_geometry'],'water_geometry_policy':mode,
            'cost_tracking':{'compute_budget':None,'wall_time_limit':None,'expected_endpoint_count':len(tasks)},
            'interpretation':'vertical_water_metal_coupling; no occupancy, binding free energy or calibrated classification',
            'state_selection':'all single-water deletions from each recorded parent; no outcome selection',
            'bulk_water_chemical_potential':None,'absolute_reference':None,'baseline_changed':False}
    write_new(out/'manifest.json',result)
    return {'manifest':record(out/'manifest.json'),'new_endpoint_count':len(tasks),'states':len(states)}


def endpoint(output, receipt, expected_xyz=None, expected_input=None):
    op=verify(output); rp=verify(receipt); rec=read_json(rp)
    if rec['returncode']!=0 or not rec['normal_termination'] or not rec['scf_converged']:
        raise InvalidArtifact('endpoint did not complete')
    if rec['artifacts']['output']['sha256']!=output['sha256']:
        raise InvalidArtifact('output receipt hash mismatch')
    if expected_xyz and rec['artifacts']['xyz']['sha256']!=expected_xyz['sha256']:
        raise InvalidArtifact('output coordinates differ')
    if expected_input and rec['artifacts']['template_input']['sha256']!=expected_input['sha256']:
        raise InvalidArtifact('output method differs')
    if rec['orca_version']!='6.1.1':raise InvalidArtifact('incompatible ORCA version')
    for name in ('xyz','template_input','runtime_input'):
        verify(rec['artifacts'][name])
    text=op.read_text(); value=energy(op)
    terms={}
    number=r'([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)'
    for name,pattern in {'CPCM_dielectric':r'CPCM Dielectric\s*:\s*'+number,
                         'dispersion':r'Dispersion correction\s+'+number,
                         'gCP':r'gCP correction\s+'+number}.items():
        found=re.findall(pattern,text)
        terms[name]=float(found[-1]) if found else None
    if terms['dispersion'] is not None and terms['gCP'] is not None:
        terms['SCF']=value-terms['dispersion']-terms['gCP']
    timestamps=[dt.datetime.fromisoformat(rec[k]) for k in ('started_at_utc','finished_at_utc')]
    return {'energy_hartree':value,'components_hartree':terms,'output':output,'receipt':receipt,
            'wall_seconds':(timestamps[1]-timestamps[0]).total_seconds(),'orca_version':rec['orca_version']}


def collect(manifest, output):
    mp=Path(manifest); m=read_json(mp); tasks={t['task_id']:t for t in m['tasks']}; rows=[]
    for state in m['states']:
        eps={};failure=[]
        for metal in ('La','Ca'):
            try:
                if state.get('reused_endpoints'):
                    e=state['reused_endpoints'][metal];eps[metal]=endpoint(e['output'],e['receipt'])
                else:
                    task=tasks[state['state_id']+'_'+metal];p=Path(task['output_path']);r=Path(str(p)+'.execution.json')
                    if not p.exists() or not r.exists():raise InvalidArtifact('not_run_or_partial')
                    eps[metal]=endpoint(record(p),record(r),task['xyz'],task['input'])
            except (InvalidArtifact,KeyError,ValueError) as error:failure.append({'metal':metal,'reason':str(error)})
        complete=not failure
        rows.append({'state_id':state['state_id'],'case':state['case'],'removed_water':state['removed_water'],
                     'status':'complete' if complete else 'unavailable','failures':failure,'endpoints':eps,
                     'R_hartree':eps['Ca']['energy_hartree']-eps['La']['energy_hartree'] if complete else None})
    squares=[];repair=[]
    for parent in m['parents']:
        full=next(r for r in rows if r['case']==parent['case'] and r['removed_water'] is None)
        archived={metal:endpoint(e['output'],e['receipt']) for metal,e in parent['endpoints'].items()}
        old=archived['Ca']['energy_hartree']-archived['La']['energy_hartree']
        repair.append({'case':parent['case'],'delta_R_kcal_mol':(full['R_hartree']-old)*HA_TO_KCAL if full['status']=='complete' else None})
        for deleted in [r for r in rows if r['case']==parent['case'] and r['removed_water'] is not None]:
            valid=full['status']==deleted['status']=='complete'
            differences={metal:full['endpoints'][metal]['energy_hartree']-deleted['endpoints'][metal]['energy_hartree'] for metal in ('Ca','La')} if valid else None
            delta=(differences['Ca']-differences['La'])*HA_TO_KCAL if valid else None
            components={}
            if valid:
                for term in ('CPCM_dielectric','dispersion','gCP','SCF'):
                    terms=[r['endpoints'][metal]['components_hartree'].get(term) for r,metal in [(full,'Ca'),(deleted,'Ca'),(full,'La'),(deleted,'La')]]
                    components[term]=(terms[0]-terms[1]-terms[2]+terms[3])*HA_TO_KCAL if None not in terms else None
                if abs(delta-(full['R_hartree']-deleted['R_hartree'])*HA_TO_KCAL)>1e-7:
                    raise InvalidArtifact('square algebra does not close')
            squares.append({'case':parent['case'],'water':deleted['removed_water']['source'],
                            'status':'complete' if valid else 'unavailable','full_state':full['state_id'],'deleted_state':deleted['state_id'],
                            'addition_energy_differences_hartree_without_water_reference':differences,
                            'delta_S_water_addition_kcal_mol':delta,'components_kcal_mol':components,
                            'interpretation':'positive: retaining this water shifts contrast toward La; no absolute occupancy claim'})
    result={'manifest':record(mp),'protocol_id':PROTOCOL,'rows':rows,'squares':squares,'water_geometry_repair_effect':repair,
            'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete',
            'classifier_improvement':None,'water_occupancy':None,'bulk_water_free_energy':None,
            'baseline_changed':False,'implementation':record(__file__)}
    write_new(output,result)
    return {'status':result['status'],'squares':squares,'repair_effect':repair}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='operation',required=True)
    for name in ('audit','prepare'):
        q=s.add_parser(name);q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    args=vars(p.parse_args());op=args.pop('operation');print(json.dumps(globals()[op](**args),indent=2))
