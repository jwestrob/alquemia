#!/usr/bin/env python3
"""Frozen archive-only donor gradients and fixed grouped classifier comparison."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import resource
import time

import gemmi
import numpy as np

from affordable_common import (InvalidArtifact, HA_TO_KCAL, BOHR_TO_A, record,
                               verify, read_json, write_new, xyz)
from affordable_response import cap_jacobians, extract, map_to_source, source_key
from coordination_policy import donor_type, selected_residue_atoms, residue_insertion_code
from site_classifier import fit, predict

PROTOCOL = 'mace_polar_archived_donor_response_probe_v1'
EV_TO_KCAL = 23.06054783061903
CUTOFF_A = 3.2
PRIMARY = ['mean_radial_delta_g_kcal_mol_A', 'tangential_fraction']
ARMS = {'DFT': ['DFT_R_kcal_mol'], 'radial': [PRIMARY[0]],
        'DFT_radial': ['DFT_R_kcal_mol', PRIMARY[0]],
        'tangential': [PRIMARY[1]], 'DFT_tangential': ['DFT_R_kcal_mol', PRIMARY[1]]}
SOURCES = {
    'pqq_manifest': 'workspaces/mace_canonical_20260916/mace_v2/medium/manifest.json',
    'pqq_collection': 'workspaces/mace_canonical_20260916/mace_v2/medium/collection_job_1200776.json',
    'pqq_inventory': 'workspaces/mace_canonical_20260916/audit_v2/inventory.json',
    'direct': 'workspaces/mace_metal_response_20260918/prepared_v2/preparation.json',
    'old_features': 'workspaces/site_classifier_20260918/features_v2/features.json',
    'old_evaluation': 'workspaces/site_classifier_20260918/evaluation_v1/result.json'}


def paired_coordinates(states):
    atoms = {m: xyz(verify(s['xyz'])) for m,s in states.items()}
    a,b = atoms['Ca'],atoms['La']
    if (a[0][0] != 'Ca' or b[0][0] != 'La' or a[0][1:] != b[0][1:] or a[1:] != b[1:]
            or states['La']['charge'] - states['Ca']['charge'] != 1
            or any(s['spin_multiplicity'] != 1 for s in states.values())):
        raise InvalidArtifact('paired coordinates/order/charge/multiplicity mismatch')
    return atoms, np.asarray([a[1:] for a in a], dtype=float)


def native_gradient(result, state, expected_model, expected_software):
    manifest = read_json(verify(result['manifest']))
    model = {k:v for k,v in manifest['model'].items() if k != 'preparation_policy'}
    if model != expected_model or manifest['software'] != expected_software:
        raise InvalidArtifact('native model or software mismatch')
    if result['status'] != 'computed' or result.get('charge_check') is not True:
        raise InvalidArtifact('native endpoint not successfully computed')
    task = next(t for t in manifest['tasks'] if t['task_id'] == result['task_id'])
    if (xyz(verify(task['xyz'])) != xyz(verify(state['xyz'])) or task['charge'] != state['charge']
            or task['spin_multiplicity'] != state['spin_multiplicity']):
        raise InvalidArtifact('force archive does not match prepared endpoint')
    forces = np.load(verify(result['forces']), allow_pickle=False)
    if forces.shape != (len(xyz(verify(state['xyz']))),3) or not np.isfinite(forces).all():
        raise InvalidArtifact('force array shape or finite-value failure')
    return -EV_TO_KCAL*forces


def source_atoms(path):
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise InvalidArtifact('source must contain one selected model')
    atoms = {}
    for ci,chain in enumerate(structure[0]):
        for ri,residue in enumerate(chain):
            for atom in selected_residue_atoms(residue):
                key = f'{ci}/{ri}/{atom.name}'
                atoms[key] = {'chain':chain.name,'chain_index':ci,'residue_index':ri,
                    'resnum':residue.seqid.num,'insertion_code':residue_insertion_code(residue),
                    'resname':residue.name,'atom':atom.name,'element':atom.element.name,
                    'xyz_A':list(atom.pos)}
    return atoms


def canonical_mapping(preparation, coordinates):
    """Expand recorded fragment order, validating source atoms and actual cap geometry."""
    src = preparation.get('source_structure')
    if src is None:
        src = read_json(verify(preparation['protonation_manifest']))['output']
    physical = source_atoms(verify(src)); mapping = []; cofactor = []
    qi = 1
    for fragment in preparation['qm_fragments']:
        first = next(a for a in fragment['atom_records'] if a['origin'] == 'source_heavy_atom')
        chain = fragment['id'].split(':',1)[0]
        candidates = [s for s in physical.values() if s['chain'] == chain and s['atom'] == first['name']
            and np.allclose(s['xyz_A'],first['xyz_A'],atol=1e-6,rtol=0)]
        if len(candidates) != 1:
            raise InvalidArtifact('ambiguous/missing canonical source fragment '+fragment['id'])
        residue = candidates[0]
        for atom in fragment['atom_records']:
            if not np.allclose(coordinates[qi],atom['xyz_A'],atol=1e-6,rtol=0):
                raise InvalidArtifact('canonical fragment atom ordering/coordinates mismatch')
            if atom['origin'] == 'generated_Cbeta_link_cap':
                kr = f"{residue['chain_index']}/{residue['residue_index']}/CB"
                ko = f"{residue['chain_index']}/{residue['residue_index']}/CA"
                xr, xo = (np.asarray(physical[k]['xyz_A']) for k in (kr,ko))
                length = 1.09  # fixed archived Cbeta-H link policy
                target = xr + length*(xo-xr)/np.linalg.norm(xo-xr)
                if not np.allclose(target,coordinates[qi],atol=1e-6,rtol=0):
                    raise InvalidArtifact('canonical cap differs from recorded physical link policy')
                ja,jb = cap_jacobians(xr,xo,length)
                mapping.append({'qm_index':qi,'kind':'cap','retained':kr,'omitted':ko,
                                'J_retained':ja.tolist(),'J_omitted':jb.tolist()})
            elif atom['origin'] in ('source_heavy_atom','source_protonation_hydrogen','generated_opposite_bisector'):
                key = f"{residue['chain_index']}/{residue['residue_index']}/{atom['name']}"
                if atom['element'] != 'H':
                    if key not in physical or not np.allclose(physical[key]['xyz_A'],coordinates[qi],atol=1e-6,rtol=0):
                        raise InvalidArtifact('canonical source-heavy coordinate mismatch')
                elif key not in physical:
                    if atom['origin'] != 'generated_opposite_bisector':
                        raise InvalidArtifact('source protonation hydrogen missing')
                    physical[key] = {**residue,'atom':atom['name'],'element':'H','xyz_A':coordinates[qi].tolist(),
                                     'origin':'physical_prepared_PQQ_proton'}
                mapping.append({'qm_index':qi,'kind':'source','key':key})
            else:
                raise InvalidArtifact('unsupported canonical atom origin '+atom['origin'])
            if fragment['kind'] == 'fixed_core_pqq': cofactor.append(qi)
            qi += 1
    if qi != len(coordinates) or len(cofactor) != 27:
        raise InvalidArtifact('canonical core or complete PQQ atom-count mismatch')
    return mapping,physical,cofactor,src


def generic_mapping(preparation, coordinates):
    physical = source_atoms(verify(preparation['source_structure'])); mapping = []
    for row in preparation['atom_graph']['source_to_qm']:
        qi = row['qm_index']
        if not np.allclose(coordinates[qi],row['xyz_A'],atol=1e-6,rtol=0):
            raise InvalidArtifact('generic source mapping differs from core coordinates')
        if row['kind'] == 'source':
            key = source_key(row['source'])
            if key not in physical or not np.allclose(physical[key]['xyz_A'],coordinates[qi],atol=1e-6,rtol=0):
                raise InvalidArtifact('generic source-heavy or proton coordinates mismatch')
            mapping.append({'qm_index':qi,'kind':'source','key':key})
        elif row['kind'] == 'sigma_link_H':
            kr,ko = (source_key(row[k]) for k in ('retained','omitted'))
            xr,xo = (np.asarray(physical[k]['xyz_A']) for k in (kr,ko))
            length = row['length_A']; target = xr + length*(xo-xr)/np.linalg.norm(xo-xr)
            if not np.allclose(target,coordinates[qi],atol=1e-6,rtol=0):
                raise InvalidArtifact('generic cap geometry mismatch')
            ja,jb = cap_jacobians(xr,xo,length)
            mapping.append({'qm_index':qi,'kind':'cap','retained':kr,'omitted':ko,
                            'J_retained':ja.tolist(),'J_omitted':jb.tolist()})
        else: raise InvalidArtifact('unsupported generic mapping')
    if sorted(r['qm_index'] for r in mapping) != list(range(1,len(coordinates))):
        raise InvalidArtifact('incomplete/duplicated atom mapping')
    return mapping,physical,[],preparation['source_structure']


def map_gradient(gradient, mapping):
    result = {'metal':np.asarray(gradient[0]).copy()}
    for row in mapping:
        g = gradient[row['qm_index']]
        updates = [(row['key'],g)] if row['kind'] == 'source' else [
            (row['retained'],np.asarray(row['J_retained']).T@g),
            (row['omitted'],np.asarray(row['J_omitted']).T@g)]
        for key,value in updates: result[key] = result.get(key,np.zeros(3)) + value
    return result


def selected_donors(mapping, physical, coordinates, preparation):
    candidates = []
    for row in mapping:
        if row['kind'] != 'source': continue
        atom = physical[row['key']]
        kind = donor_type(atom['resname'],atom['atom'],atom['element'])
        if atom['resname'] == 'PQQ':
            pqq = preparation['pqq']
            if 'residues' in pqq:
                allowed = pqq['residues'][0]['allowed_donor_atoms']
            elif pqq.get('schema_id') == 'pdb_ccd_pqq_v1':
                allowed = ('O5','N6','O7A','O7B')
            else: raise InvalidArtifact('unknown PQQ donor schema')
            if atom['atom'] not in allowed: kind = None
        distance = float(np.linalg.norm(np.asarray(atom['xyz_A'])-coordinates[0]))
        if kind is not None and distance <= CUTOFF_A:
            candidates.append({'key':row['key'],'qm_index':row['qm_index'],'donor_type':kind,
                               'distance_A':distance,**atom})
    if not candidates: raise InvalidArtifact('no retained typed donors within fixed cutoff')
    return sorted(candidates,key=lambda r:(r['distance_A'],r['key']))


def observables(delta_gradient, donor_positions, metal_position):
    g = np.asarray(delta_gradient,dtype=float); d = np.asarray(donor_positions)-metal_position
    if g.shape != d.shape or not np.isfinite(g).all(): raise InvalidArtifact('invalid donor gradients')
    length = np.linalg.norm(d,axis=1)
    if np.any(length == 0): raise InvalidArtifact('donor overlaps metal')
    unit = d/length[:,None]; radial = np.einsum('ij,ij->i',g,unit)
    tangent = g-radial[:,None]*unit; total = float(np.sum(g*g))
    if total <= 0: raise InvalidArtifact('undefined zero differential donor load')
    return {PRIMARY[0]:float(np.mean(radial)), PRIMARY[1]:float(np.sum(tangent*tangent)/total),
            'donor_load_coherence':float(np.linalg.norm(g.sum(axis=0))/np.linalg.norm(g,axis=1).sum()),
            'radial_components_kcal_mol_A':radial.tolist(),
            'tangential_components_kcal_mol_A':tangent.tolist()}


def describe_pair(gradients, mapping, physical, coords, donors, cofactor):
    mapped = {m:map_gradient(g,mapping) for m,g in gradients.items()}
    delta = np.asarray([mapped['Ca'][d['key']]-mapped['La'][d['key']] for d in donors])
    positions = np.asarray([d['xyz_A'] for d in donors])
    values = observables(delta,positions,coords[0])
    closure = {}
    for metal,g in gradients.items():
        sources = mapped[metal]
        net = sum(sources.values())
        torque = sum((np.cross(np.asarray(physical[k]['xyz_A'])-coords[0],v)
                      for k,v in sources.items() if k != 'metal'),np.zeros(3))
        closure[metal] = {'net_gradient_error_kcal_mol_A':float(np.linalg.norm(net-g.sum(axis=0))),
                          'torque_error_kcal_mol':float(np.linalg.norm(torque-np.cross(coords-coords[0],g).sum(axis=0)))}
    if max(v for row in closure.values() for v in row.values()) > 1e-4:
        raise InvalidArtifact('physical cap force/torque closure failed')
    # Algebraic invariance on the same observed vectors; not a new model evaluation.
    rotation = np.array([[0.,-1.,0.],[0.,0.,-1.],[1.,0.,0.]])
    transformed = observables(delta@rotation.T,positions@rotation.T+[1.5,-3.,7.],coords[0]@rotation.T+[1.5,-3.,7.])
    rigid_error = max(abs(values[k]-transformed[k]) for k in PRIMARY+['donor_load_coherence'])
    if rigid_error > 1e-9: raise InvalidArtifact('observable rigid transformation algebra failed')
    cofactor_result = {'status':'not_present','force_norm_kcal_mol_A':None,'torque_norm_kcal_mol':None}
    if cofactor:
        force = -(gradients['Ca']-gradients['La'])[cofactor]
        offset = coords[cofactor]-coords[cofactor].mean(axis=0)
        cofactor_result = {'status':'descriptive','atom_count':len(cofactor),
            'origin':'geometric_centroid','force_vector_kcal_mol_A':force.sum(axis=0).tolist(),
            'torque_vector_kcal_mol':np.cross(offset,force).sum(axis=0).tolist(),
            'force_norm_kcal_mol_A':float(np.linalg.norm(force.sum(axis=0))),
            'torque_norm_kcal_mol':float(np.linalg.norm(np.cross(offset,force).sum(axis=0)))}
    return {'features':{k:values[k] for k in PRIMARY},'observables':values,
            'donor_delta_gradient_kcal_mol_A':delta.tolist(),'cofactor':cofactor_result,
            'mapping_checks':closure,'rigid_algebra_max_error':rigid_error}


def prepare(root, agreement, output):
    start,cpu = time.monotonic(),time.process_time(); root = Path(root).resolve()
    sources = {k:record(root/v) for k,v in SOURCES.items()}
    data = {k:read_json(verify(v)) for k,v in sources.items()}
    model = {k:v for k,v in data['pqq_manifest']['model'].items() if k != 'preparation_policy'}
    software = data['pqq_manifest']['software']; rows=[]; failures=[]; receipts=[]
    oldrows = {r['case_id']:r for r in data['old_features']['rows']}
    native = data['pqq_collection']['rows'] | {k:v['result'] for k,v in data['pqq_collection']['reused_rows'].items()}
    for case in data['pqq_inventory']['rows']:
        name=case['case_id']
        try:
            prep = read_json(verify(case['source_manifest'])); states = case['endpoints']
            atoms,coords = paired_coordinates(states)
            gradients = {m:native_gradient(native[name+'_'+m],states[m],model,software) for m in ('Ca','La')}
            mapping,physical,cofactor,source = canonical_mapping(prep,coords)
            donors = selected_donors(mapping,physical,coords,prep)
            description = describe_pair(gradients,mapping,physical,coords,donors,cofactor)
            old=oldrows.get(name)
            row={'case_id':name,'target':'PQQ_functional_class','label':case['expected_class'],
                'role':case['evaluation_role'],'biological_group':old['biological_group'] if old else case['sequence_accession_group'],
                'validation_group':data['old_evaluation']['group_assignments'].get(name),
                'evidence_stratum':case['evidence_stratum'],'DFT_protocol':data['pqq_inventory']['source_protocol_id'],
                'MACE_protocol':PROTOCOL,'preparation':case['source_manifest'],'source_structure':source,
                'states':states,'donors':donors,'physical_mapping':mapping,**description,
                'force_artifacts':{m:native[name+'_'+m]['forces'] for m in ('Ca','La')},
                'prospectively_blind':False,'status':'complete'}
            row['features']['DFT_R_kcal_mol']=case['baseline']['published_R_kcal_mol']
            rows.append(row)
            receipts += [native[name+'_'+m] for m in ('Ca','La')]
        except (InvalidArtifact,KeyError,ValueError,StopIteration) as error:
            failures.append({'case_id':name,'target':'PQQ_functional_class','status':'unavailable','reason':str(error)})
    for name,case in data['direct']['cases'].items():
        try:
            centers=data['direct']['centers'][name]
            prep=read_json(verify(case['source_preparation']))
            states={m:centers[m]['core_state'] for m in ('Ca','La')}; atoms,coords=paired_coordinates(states)
            gradients={m:native_gradient(centers[m]['core'],states[m],model,software) for m in ('Ca','La')}
            mapping,physical,cofactor,source=generic_mapping(prep,coords)
            donors=selected_donors(mapping,physical,coords,prep)
            description=describe_pair(gradients,mapping,physical,coords,donors,cofactor)
            dft={}
            for m in ('Ca','La'):
                art=centers[m]['DFT_artifacts']
                for a in art.values():verify(a)
                if xyz(verify(art['xyz'])) != atoms[m]:raise InvalidArtifact('DFT/MACE coordinate mismatch')
                parsed=extract(art['engrad']['path'],art['output']['path'],art['input']['path'],art['xyz']['path'])
                dft[m]=np.asarray(parsed['gradient_kcal_mol_per_A'])
                # Existing independent mapping adapter must agree with this compact mapper.
                reference=map_to_source(dft[m],prep); mapped=map_gradient(dft[m],mapping)
                if set(reference)!=set(mapped) or any(not np.allclose(reference[k],mapped[k],atol=1e-10,rtol=0) for k in mapped):
                    raise InvalidArtifact('generic mapping differs from existing physical gradient adapter')
            dft_description=describe_pair(dft,mapping,physical,coords,donors,cofactor)
            a=np.asarray(description['donor_delta_gradient_kcal_mol_A']).ravel()
            b=np.asarray(dft_description['donor_delta_gradient_kcal_mol_A']).ravel()
            correspondence={'donor_delta_gradient_cosine':float(a@b/(np.linalg.norm(a)*np.linalg.norm(b))),
                'MACE_vacuum_vs_DFT_CPCM':'different_Hamiltonians_not_force_validation',
                'radial_sign_agreement':int(np.sum(np.sign(description['observables']['radial_components_kcal_mol_A']) ==
                                                   np.sign(dft_description['observables']['radial_components_kcal_mol_A']))),
                'donor_count':len(donors)}
            evidence_case='GGR_1GLG' if name.startswith('GGR') else name
            evidence=oldrows[evidence_case]
            rows.append({'case_id':name,'target':'direct_site_affinity_direction','label':evidence['label'],
                'evidence_stratum':evidence['evidence_stratum'],'biological_group':evidence['biological_group'],
                'metadata_evidence':{'source':sources['old_features'],'case_id':evidence_case},
                'role':'consumed_development','status':'complete','preparation':case['source_preparation'],
                'source_structure':source,'states':states,'physical_mapping':mapping,'donors':donors,**description,
                'DFT_response':dft_description,'DFT_correspondence':correspondence,
                'DFT_artifacts':{m:centers[m]['DFT_artifacts'] for m in ('Ca','La')},
                'force_artifacts':{m:centers[m]['core']['forces'] for m in ('Ca','La')},
                'DFT_R_kcal_mol':(centers['Ca']['DFT_energy_Ha']-centers['La']['DFT_energy_Ha'])*HA_TO_KCAL,
                'prospectively_blind':False,'calibrated_class':None})
            receipts += [centers[m]['core'] for m in ('Ca','La')]
        except (InvalidArtifact,KeyError,ValueError,StopIteration) as error:
            failures.append({'case_id':name,'target':'direct_site_affinity_direction','status':'unavailable','reason':str(error)})
    result={'protocol':PROTOCOL,'agreement':record(agreement),'sources':sources,'model':model,'software':software,
        'arms':ARMS,'rows':rows,'failures':failures,'expected_cases':{'PQQ':28,'direct':4},
        'new_molecular_calls':0,'archived_MACE_endpoints_reused':len(receipts),
        'archived_MACE_model_evaluation_seconds':sum(r['evaluation_seconds'] for r in receipts),
        'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,
        'max_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'implementation':{p.name:record(p) for p in [Path(__file__),Path(__file__).with_name('affordable_response.py'),
            Path(__file__).with_name('coordination_policy.py'),Path(__file__).with_name('site_classifier.py')]},
        'baseline_changed':False,'relaxation_correction_kcal_mol':None,'entropy_correction_kcal_mol':None}
    write_new(output,result)
    return {'output':record(output),'complete_cases':len(rows),'failures':failures,
            'CPU_seconds':result['CPU_seconds'],'wall_seconds':result['wall_seconds']}


def evaluate(features, output):
    start,cpu=time.monotonic(),time.process_time(); data=read_json(features); verify(data['agreement'])
    if data['protocol']!=PROTOCOL or data['arms']!=ARMS:raise InvalidArtifact('changed feature protocol or arms')
    old=read_json(verify(data['sources']['old_evaluation']))
    selected=[r for r in data['rows'] if r['target']=='PQQ_functional_class']
    train=[r for r in selected if r['role']=='calibration']; transfers=[r for r in selected if r['role']!='calibration']
    if len(train)!=25:raise InvalidArtifact('complete frozen canonical population required; do not fit reduced subset')
    if any(r['validation_group'] is None for r in train):raise InvalidArtifact('missing fixed training group')
    results={}; all_models={}; groups=sorted({r['validation_group'] for r in train})
    for arm,names in ARMS.items():
        predictions=[];folds=[]
        for group in groups:
            held=[r for r in train if r['validation_group']==group]
            used=[r for r in train if r['validation_group']!=group]
            model=fit(used,names);folds.append({'held_group':group,'model':model})
            predictions += [{'case_id':r['case_id'],'label':r['label'],'group':group,**predict(model,r)} for r in held]
        model=fit(train,names);all_models[arm]=model
        results[arm]={'correct':sum(r['correct'] for r in predictions),'total':25,'predictions':predictions,'folds':folds,
            'consumed_transfers':[{'case_id':r['case_id'],'label':r['label'],
                'same_group_in_training':None if r['validation_group'] is None else r['validation_group'] in groups,
                **predict(model,r)} for r in transfers]}
    # Standardizing R and S differs only by common offset. Do not quietly change the comparator.
    old_dft={r['case_id']:r for r in old['results']['PQQ_functional_class']['arms']['DFT']['out_of_group']}
    max_dft_replay=max(abs(r['logit']-old_dft[r['case_id']]['logit']) for r in results['DFT']['predictions'])
    if max_dft_replay > 1e-7 or any(r['class']!=old_dft[r['case_id']]['class'] for r in results['DFT']['predictions']):
        raise InvalidArtifact('DFT-only grouped comparator failed archived S-only replay')
    direct=[r for r in data['rows'] if r['target']=='direct_site_affinity_direction']
    direct_order={k:{'MACE_ascending':sorted([{'case_id':r['case_id'],'value':r['features'][k]} for r in direct],key=lambda r:r['value']),
        'DFT_ascending':sorted([{'case_id':r['case_id'],'value':r['DFT_response']['features'][k]} for r in direct],key=lambda r:r['value'])}
        for k in PRIMARY}
    gg={r['case_id']:r for r in direct if r['case_id'].startswith('GGR')}
    partition={k:{'MACE_connected_minus_extended':gg['GGR_connected']['features'][k]-gg['GGR_extended']['features'][k],
                  'DFT_connected_minus_extended':gg['GGR_connected']['DFT_response']['features'][k]-gg['GGR_extended']['DFT_response']['features'][k]}
        for k in PRIMARY} if len(gg)==2 else {'status':'unavailable'}
    result={'protocol':PROTOCOL,'features':record(features),'agreement':data['agreement'],'results':results,'models':all_models,
        'fixed_group_count':len(groups),'DFT_S_to_R_max_logit_replay_error':max_dft_replay,
        'earlier_models':{k:{n:v[n] for n in ('correct','predicted','total')} for k,v in old['results']['PQQ_functional_class']['arms'].items()},
        'direct_orderings':direct_order,'GGR_partition_feature_differences':partition,
        'direct_cases':len(direct),'direct_biological_groups':len({r['biological_group'] for r in direct}),
        'failures':data['failures'],'new_molecular_calls':0,'baseline_changed':False,
        'prospectively_blind':False,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,
        'max_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'implementation':record(__file__)}
    write_new(output,result)
    return {'output':record(output),'summary':{a:{k:r[k] for k in ('correct','total','consumed_transfers')} for a,r in results.items()},
            'direct_orderings':direct_order,'CPU_seconds':result['CPU_seconds'],'wall_seconds':result['wall_seconds']}


def transfer(features, evaluation, agreement, output):
    """Cross-target differences only: reuse saved PQQ weights; never fit or classify."""
    start,cpu=time.monotonic(),time.process_time()
    data=read_json(features); saved=read_json(evaluation)
    if data['protocol']!=PROTOCOL or saved['protocol']!=PROTOCOL:
        raise InvalidArtifact('unexpected response protocol')
    rows={r['case_id']:r for r in data['rows'] if r['target']=='direct_site_affinity_direction'}
    if set(rows)!={'ALPHA_1F6S','ALPHA_6IP9','GGR_extended','GGR_connected'}:
        raise InvalidArtifact('complete declared direct population required')
    old=read_json(verify(data['sources']['old_features']))
    labels={r['case_id']:r for r in old['rows']}
    for name,row in rows.items():
        truth=labels['GGR_1GLG' if name.startswith('GGR') else name]
        if any(row[k]!=truth[k] for k in ('label','biological_group','evidence_stratum')):
            raise InvalidArtifact('direct reporting metadata differs from authoritative evidence ledger')
    results=[]
    for alpha in ('ALPHA_1F6S','ALPHA_6IP9'):
        for ggr in ('GGR_extended','GGR_connected'):
            a,b=rows[alpha],rows[ggr]
            delta={'DFT_R_kcal_mol':a['DFT_R_kcal_mol']-b['DFT_R_kcal_mol'],
                   PRIMARY[0]:a['features'][PRIMARY[0]]-b['features'][PRIMARY[0]]}
            arm_results={}
            for arm in ('DFT','DFT_radial'):
                model=saved['models'][arm]
                if model['features']!=ARMS[arm]:raise InvalidArtifact('changed fitted feature arm')
                if len(model['training_cases'])!=25:raise InvalidArtifact('unexpected fitted PQQ population')
                components={k:float(coef*delta[k]/scale) if active else 0.
                    for k,coef,scale,active in zip(model['features'],model['coefficients'],model['scale'],model['active'])}
                value=sum(components.values())
                arm_results[arm]={'alpha_minus_GGR_logit':value,'components':components,
                                 'expected_order':value>0,'absolute_class':None,'probability':None}
            results.append({'alpha':alpha,'GGR':ggr,'feature_differences':delta,'arms':arm_results})
    result={'protocol':PROTOCOL,'features':record(features),'saved_fits':record(evaluation),'agreement':record(agreement),
        'new_fits':0,'new_molecular_calls':0,'rows':results,
        'expected_order_counts':{a:sum(r['arms'][a]['expected_order'] for r in results) for a in ('DFT','DFT_radial')},
        'comparison_count':4,'biological_groups':2,'independent_biological_comparisons':1,
        'interpretation':'exploratory_cross_target_transfer_ordering_on_consumed_preparations',
        'absolute_affinity_calibration':None,'baseline_changed':False,
        'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'implementation':record(__file__)}
    write_new(output,result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='operation',required=True)
    p=sub.add_parser('prepare')
    for name in ('root','agreement','output'):p.add_argument('--'+name,required=True)
    p=sub.add_parser('evaluate')
    for name in ('features','output'):p.add_argument('--'+name,required=True)
    p=sub.add_parser('transfer')
    for name in ('features','evaluation','agreement','output'):p.add_argument('--'+name,required=True)
    args=vars(parser.parse_args()); operation=args.pop('operation')
    print(json.dumps(globals()[operation](**args),indent=2))
