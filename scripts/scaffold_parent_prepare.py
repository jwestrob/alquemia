"""Exact ff19SB protein parents and source-connected targets; no energy Context."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import resource
import time
import numpy as np
import openmm as mm
from openmm import app, unit
from affordable_common import InvalidArtifact, paired, read_json, record, verify, write_new, xyz
from mace_site_kinematics import Kinematics
from mace_hybrid import write_xyz
import scaffold_environment_inventory as old

PROTOCOL = 'common8_exact_source_ff19SB_parent_and_matched_targets_v1'
TOL = 1e-12
TARGETS = ('origin', 'Ca_adaptive', 'La_adaptive')
DONOR_NAMES = {'anchor_glutamate': ('OE1', 'OE2'), 'anchor_asparagine': ('OD1',),
               'extra_acidic_ligand_homolog': ('OD1', 'OD2')}


def save(path, value):
    write_new(path, value)
    return record(path)


def exact(a, b, what):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise InvalidArtifact(what + ': invalid shape/coordinates')
    delta = float(np.max(np.abs(a-b))) if a.size else 0.
    if delta > TOL: raise InvalidArtifact(what + ': source coordinate mismatch ' + str(delta))
    return delta


def topology_record(topology, positions):
    atoms = [{'parent_index': a.index, 'source_id': old.atom_id(a), 'name': a.name,
              'element': a.element.symbol, 'residue_index': a.residue.index,
              'chain': a.residue.chain.id, 'resnum': a.residue.id,
              'insertion_code': a.residue.insertionCode.strip(), 'resname': a.residue.name}
             for a in topology.atoms()]
    residues = [{'residue_index': r.index, 'chain': r.chain.id, 'resnum': r.id,
                 'insertion_code': r.insertionCode.strip(), 'resname': r.name,
                 'parent_indices': [a.index for a in r.atoms()]} for r in topology.residues()]
    bonds = [{'indices': [a.index, b.index], 'source_length_A': float(np.linalg.norm(positions[a.index]-positions[b.index]))}
             for a, b in topology.bonds()]
    return {'atoms': atoms, 'residues': residues, 'bonds': bonds, 'units': 'angstrom'}


def context_map(g, top, positions):
    atoms = top['atoms']; byid = {a['source_id']: a['parent_index'] for a in atoms}
    pmap = {}; nodes = []; local = set(); links = []; caps = []; other = []; max_delta = 0.
    for i, meta in enumerate(g['source_atom_metadata']):
        if meta is None or meta['resname'] not in old.STANDARD: continue
        aid = old.metadata_id(meta)
        if aid not in byid: raise InvalidArtifact('missing protein source atom ' + aid)
        j = byid[aid]
        if atoms[j]['element'] != meta['element']: raise InvalidArtifact('mapped element mismatch')
        max_delta = max(max_delta, exact(positions[j], g['positions_A'][i], aid)); pmap[i] = j
        nodes.append({'physical_index': i, 'parent_index': j, 'source_id': aid})
    for l in g['core_links']:
        ci, kind = l[:2]
        if kind == 'source' and l[2] in pmap:
            j = pmap[l[2]]; local.add(j)
            links.append({'context_index': ci, 'parent_index': j, 'source_id': atoms[j]['source_id']})
        elif kind == 'cap':
            if l[2] not in pmap or l[3] not in pmap: raise InvalidArtifact('unsupported nonprotein cap')
            caps.append({'context_index': ci, 'retained_parent_index': pmap[l[2]],
                         'omitted_parent_index': pmap[l[3]], 'length_A': l[4], 'q0_xyz_A': l[5]})
        else: other.append(l)
    bonds = {tuple(sorted(b['indices'])) for b in top['bonds']}
    crossing = [list(b) for b in sorted(bonds) if old.partition(b, local) == 'mixed']
    cap_bonds = sorted(tuple(sorted((c['retained_parent_index'], c['omitted_parent_index']))) for c in caps)
    if cap_bonds != [tuple(b) for b in crossing]: raise InvalidArtifact('actual cut bond/cap mismatch')
    return {'physical_nodes': nodes, 'source_context_atoms': links, 'caps': caps,
            'other_source_links': other, 'actual_cut_bonds': crossing,
            'unmapped_fixed_context_indices': sorted(set(range(len(g['core_positions_A'])))-{l[0] for l in g['core_links']}),
            'local_parent_indices': sorted(local), 'max_source_coordinate_difference_A': max_delta}


def reconstruct_context(parent_positions, mapping, g):
    """Keep all fixed source atoms and reproduce the existing cap-offset recipe."""
    p = np.asarray(parent_positions); out = np.asarray(g['core_positions_A'], dtype=float).copy()
    q0 = np.asarray(g['positions_A']); physical = {r['parent_index']: r['physical_index'] for r in mapping['physical_nodes']}
    for r in mapping['source_context_atoms']: out[r['context_index']] = p[r['parent_index']]
    for c in mapping['caps']:
        a, b = c['retained_parent_index'], c['omitted_parent_index']; v = p[b]-p[a]; v0 = q0[physical[b]]-q0[physical[a]]
        if np.linalg.norm(v) < 1e-10: raise InvalidArtifact('collapsed boundary bond')
        # Archived caps may carry rounding offsets. The zero displacement limit is exact.
        out[c['context_index']] += p[a] + c['length_A']*v/np.linalg.norm(v) - np.asarray(c['q0_xyz_A'])
    return out


def mobile_set(top, positions, core_g, core_atoms):
    core_points = np.asarray(core_g['core_positions_A'])[[i for i,a in enumerate(core_atoms) if a[0] != 'H']]
    atoms = top['atoms']; heavy = [a['parent_index'] for a in atoms if a['element'] != 'H']
    nearest = np.linalg.norm(positions[heavy,None,:]-core_points[None,:,:], axis=2).min(axis=1)
    distances = dict(zip(heavy, map(float, nearest)))
    geometric = {atoms[i]['residue_index'] for i in heavy if distances[i] <= 8.}
    byid = {a['source_id']: a['parent_index'] for a in atoms}
    core_residues = set(); core_parent_indices = set()
    for l in core_g['core_links']:
        if l[1] != 'source': continue
        meta = core_g['source_atom_metadata'][l[2]]
        if meta and meta['resname'] in old.STANDARD:
            j = byid[old.metadata_id(meta)]; core_parent_indices.add(j)
            core_residues.add(atoms[j]['residue_index'])
    seed = geometric | core_residues; neighbours = set(); peptide_bonds = []
    for b in top['bonds']:
        a, c = [atoms[i] for i in b['indices']]
        if a['residue_index'] == c['residue_index'] or {a['name'], c['name']} != {'C','N'}: continue
        peptide_bonds.append(b['indices'])
        if a['residue_index'] in seed: neighbours.add(c['residue_index'])
        if c['residue_index'] in seed: neighbours.add(a['residue_index'])
    selected = seed | neighbours
    ids = [a['parent_index'] for a in atoms if a['residue_index'] in selected]
    return {'rule': 'q0_complete_residues_with_heavy_within_8A_of_canonical_core_heavy_plus_core_contributors_plus_one_actual_peptide_neighbor_layer',
            'cutoff_A': 8., 'canonical_core_heavy_positions_A':core_points.tolist(),
            'core_parent_indices':sorted(core_parent_indices), 'geometric_residue_indices': sorted(geometric), 'core_residue_indices': sorted(core_residues),
            'seed_residue_indices': sorted(seed), 'peptide_neighbor_residue_indices': sorted(neighbours-seed),
            'mobile_residue_indices': sorted(selected), 'mobile_parent_indices': ids,
            'frozen_parent_indices': sorted(set(range(len(atoms)))-set(ids)), 'actual_peptide_bonds': peptide_bonds,
            'residue_min_core_heavy_distance_A': {str(r['residue_index']): min(distances[i] for i in r['parent_indices'] if i in distances) for r in top['residues']}}


def donor_targets(carve, top, positions, metal_position):
    byid = {a['source_id']: a for a in top['atoms']}; rows = []
    for role, names in DONOR_NAMES.items():
        r = carve['fixed_core']['requested_roles'][role]
        if isinstance(r, str):
            m = re.fullmatch(r'([^:]+):([A-Z]+)(-?\d+)', r)
            if not m: raise InvalidArtifact('unsupported legacy residue selector ' + r)
            r = {'chain': m[1], 'resname': m[2], 'resnum': int(m[3]), 'icode': ''}
        if role == 'extra_acidic_ligand_homolog' and r['resname'] != 'ASP':
            if r['resname'] in ('ALA', 'ASN', 'SER', 'THR', 'GLY'): continue
            raise InvalidArtifact('unsupported extra-acid donor identity ' + r['resname'])
        if r['resname'] != {'anchor_glutamate':'GLU','anchor_asparagine':'ASN','extra_acidic_ligand_homolog':'ASP'}[role]:
            raise InvalidArtifact('donor role residue mismatch')
        for name in names:
            aid = f"{r['chain']}/{r['resnum']}/{r.get('icode', r.get('insertion_code',''))}/{name}"
            if aid not in byid: raise InvalidArtifact('missing source donor ' + aid)
            a = byid[aid]; j = a['parent_index']
            rows.append({'role': role, 'source_id': aid, 'parent_index': j,
                         'element': a['element'], 'q0_metal_distance_A': float(np.linalg.norm(positions[j]-metal_position))})
    return {'rule':'uniform_role_atom_ids_without_metal_distance_filter', 'atoms': rows,
            'archived_coordination': carve.get('coordination'), 'fixed_roles': carve['fixed_core']['requested_roles']}


def prepare_case(source_row, config, out):
    cid = source_row['case_id']; collection = read_json(verify(source_row['proposal_collection']))
    manifest = read_json(verify(collection['manifest']))
    tasks = {z: next(t for t in manifest['tasks'] if (t['case_id'],t['metal']) == (cid,z)) for z in ('Ca','La')}
    case = next(c for c in manifest['cases'] if c['case_id']==cid); case = case.get('source_case', case)
    paired(verify(tasks['La']['xyz']), verify(tasks['Ca']['xyz']), tasks['La']['charge'], tasks['Ca']['charge'])
    maps = {z: read_json(verify(t['mapping'])) for z,t in tasks.items()}
    if maps['Ca'] != maps['La'] or tasks['Ca']['source_preparation'] != tasks['La']['source_preparation']:
        raise InvalidArtifact('paired physical preparation differs')
    prep = read_json(verify(tasks['Ca']['source_preparation'])); source = verify(prep['source'])
    carve_pin = case['source']['parent']; carve = read_json(verify(carve_pin)); proton = carve['protonation_manifest']
    if read_json(verify(proton))['output'] != prep['source']: raise InvalidArtifact('source protonation mismatch')
    model, positions, system, ff, state = old.parent_model(source, verify(config['forcefield']), config['terminal_repairs'].get(cid))
    top = topology_record(model.topology, positions); g = maps['Ca']['context']
    mapping = context_map(g, top, positions); d = out/cid; d.mkdir(parents=True)
    archive = next((r for r in config['archived_rows'] if r['case_id']==cid), None)
    reuse = None
    if archive:
        if archive['source'] != prep['source']: raise InvalidArtifact('archived parent source differs')
        prior = read_json(verify(archive['artifacts']['parent_atoms.json']))
        if [a['source_id'] for a in prior] != [a['source_id'] for a in top['atoms']]: raise InvalidArtifact('archived parent order differs')
        exact([a['xyz_A'] for a in prior], positions, 'archived parent')
        if read_json(verify(archive['artifacts']['context_mapping.json'])) != mapping: raise InvalidArtifact('archived boundary mapping differs')
        old_system = verify(archive['artifacts']['parent_system.xml'])
        if mm.XmlSerializer.serialize(system) != old_system.read_text(): raise InvalidArtifact('archived protein Hamiltonian differs')
        system_pin = archive['artifacts']['parent_system.xml']; reuse = {'system': system_pin, 'mapping': archive['artifacts']['context_mapping.json']}
    else:
        (d/'parent_system.xml').write_text(mm.XmlSerializer.serialize(system)); system_pin = record(d/'parent_system.xml')
    parent = {'system': system_pin, 'topology': save(d/'topology.json',top), 'atoms': save(d/'atoms.json',top['atoms']),
              'positions_A':save(d/'positions_A.json',positions.tolist()),
              'source_atom_ids':save(d/'source_atom_ids.json',[a['source_id'] for a in top['atoms']]),
              'state':save(d/'source_state.json',state), 'forcefield':config['forcefield'], 'archived_reuse':reuse}
    cp = save(d/'context_parent_mapping.json',mapping)
    core_atoms = xyz(verify(case['source']['endpoints']['Ca']['xyz']))
    exact([a[1:] for a in core_atoms], maps['Ca']['core']['core_positions_A'], 'canonical core')
    mobile = mobile_set(top, positions, maps['Ca']['core'], core_atoms)
    donors = donor_targets(carve, top, positions, np.asarray(g['positions_A'])[g['metal_index']])
    if not {a['parent_index'] for a in donors['atoms']} <= set(mapping['local_parent_indices']): raise InvalidArtifact('donor absent from scoring context')
    pdb = app.PDBFile(str(source)); pp = np.asarray(pdb.positions.value_in_unit(unit.angstrom))
    hetero = [{'source_id':old.atom_id(a),'name':a.name,'element':a.element.symbol,'resname':a.residue.name,'xyz_A':pp[a.index].tolist()}
              for a in pdb.topology.atoms() if a.residue.name not in old.STANDARD]
    nonprotein = {'source_atoms':hetero,'source_water_count':sum(a['resname'] in ('HOH','WAT') and a['element']=='O' for a in hetero),
                  'fixed_context_atoms': [{'context_index':i,'element':a[0],'xyz_A':list(a[1:])} for i,a in enumerate(xyz(verify(tasks['Ca']['xyz'])))
                    if i not in {a['context_index'] for a in mapping['source_context_atoms']} | {a['context_index'] for a in mapping['caps']}],
                  'metal_is_endpoint_dependent_element_only':True,'FF_parameters_assigned':False}
    kin = Kinematics(g); targets = {}; base = next(c for c in read_json(verify(source_row['pool_collection']))['cases'] if c['case_id']==cid)
    for name in TARGETS:
        if name == 'origin': q = np.zeros(len(kin.modes)); archive_coord = tasks['Ca']['xyz']; candidate = None
        else:
            z = name.split('_')[0]; e = next(e for e in collection['endpoints'] if (e['case_id'],e['metal'])==(cid,z))
            if e['status'] != 'candidate_available': raise InvalidArtifact('required adaptive target unavailable')
            candidate = e['candidate']; q = np.asarray(candidate['full_q']); archive_coord = candidate['coordinate']
        physical, context = kin.evaluate(q)[:2]; p = positions.copy()
        for r in mapping['physical_nodes']: p[r['parent_index']] = physical[r['physical_index']]
        if name == 'origin': p = positions.copy()  # preserve exact parent source serialization
        reconstructed = reconstruct_context(p,mapping,g)
        delta = exact(reconstructed, context, name+' reconstructed context')
        exact(context,[a[1:] for a in xyz(verify(archive_coord))], name+' archived context')
        alias = 'origin' if name=='origin' else 'adaptive_'+name.split('_')[0]
        exact(context,[a[1:] for a in xyz(verify(base['aliases'][alias]['source_xyz']))], name+' pool context')
        frozen = mobile['frozen_parent_indices']; exact(p[frozen],positions[frozen],name+' frozen protein')
        td = d/'targets'/name; td.mkdir(parents=True); coords = {}
        for z,t in tasks.items():
            atoms = xyz(verify(t['xyz'])); xp = td/(z+'_context.xyz')
            write_xyz(xp,[(a[0],*map(float,v)) for a,v in zip(atoms,context)]); coords[z]=record(xp)
        targets[name] = {'full_q':q.tolist(), 'parent_positions_A':save(td/'parent_positions_A.json',p.tolist()),
                         'context_coordinates':coords,'actual_archived_coordinate':archive_coord,
                         'proposal_collection':source_row['proposal_collection'],'maximum_reconstruction_difference_A':delta,
                         'donor_positions_A':[p[a['parent_index']].tolist() for a in donors['atoms']],
                         'status':'prepared_geometry_only'}
    row = {'case_id':cid,'status':'prepared','parent':parent,'source':prep['source'],'source_carve':carve_pin,
           'source_protonation':proton,'source_preparation':tasks['Ca']['source_preparation'],
           'context_maps':{z:t['mapping'] for z,t in tasks.items()},
           'origins':{z:{k:t[k] for k in ('xyz','charge','multiplicity')} for z,t in tasks.items()},
           'context_parent_mapping':cp,'mobile_set':save(d/'mobile_set.json',mobile),'donors':save(d/'donors.json',donors),
           'nonprotein_inventory':save(d/'nonprotein_inventory.json',nonprotein),'targets':targets,
           'source_registry_row':source_row, 'parent_atom_count':len(positions), 'mobile_atom_count':len(mobile['mobile_parent_indices']),
           'mobile_residue_count':len(mobile['mobile_residue_indices']), 'donor_atom_count':len(donors['atoms']),
           'cap_count':len(mapping['caps']),'energy_evaluations':0,'force_evaluations':0}
    save(d/'CASE.json',row); return row


def prepare(inputs_path, agreement, archive_path, output):
    start=time.monotonic(); inputs=read_json(inputs_path); archive=read_json(archive_path)
    config=read_json(verify(archive['config'])); config['archived_rows']=[r for r in archive['rows'] if r['case_id'] in ('1H4I','4MAE')]
    if len(inputs['cases']) != 8: raise InvalidArtifact('expected declared common8')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);rows=[]
    for r in inputs['cases']:
        try: row=prepare_case(r,config,out)
        except Exception as exc: row={'case_id':r['case_id'],'status':'unavailable','reason':str(exc),'exception':type(exc).__name__,'source_registry_row':r}
        rows.append(row);print(row['case_id'],row['status'],row.get('reason',''),flush=True)
    result={'protocol_id':PROTOCOL,'inputs':record(inputs_path),'agreement':record(agreement),'archive':record(archive_path),
            'target_mapping_agreement':record(Path(agreement).parent/'TARGET_MAPPING.md'),
            'implementation':record(__file__),'parent_implementation':record(old.__file__),'openmm_version':mm.__version__,
            'forcefield':config['forcefield'],'cases':rows,'case_denominator':8,'prepared_cases':sum(r['status']=='prepared' for r in rows),
            'expected_target_denominator':24,'prepared_targets':sum(len(r.get('targets',{})) for r in rows),
            'geometry_precision_tolerance_A':TOL,'source_H_normalization':False,'new_terminal_repairs':0,
            'energy_evaluations':0,'force_evaluations':0,'molecular_searches':0,'production_changed':False,
            'wall_seconds':time.monotonic()-start,'process_CPU_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime}
    save(out/'manifest.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',required=True);p.add_argument('--agreement',required=True)
    p.add_argument('--archive',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    r=prepare(a.inputs,a.agreement,a.archive,a.output)
    print(json.dumps({k:r[k] for k in ('protocol_id','prepared_cases','prepared_targets','wall_seconds','energy_evaluations')}))
if __name__=='__main__':main()
