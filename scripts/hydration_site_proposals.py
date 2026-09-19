"""Propose water sites by aligning actual protein atoms; no energy or occupancy."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import gemmi
from affordable_common import InvalidArtifact, read_json, record, verify, write_new


def key(source):
    return tuple(source[k] for k in ('chain', 'resnum', 'insertion_code', 'canonical_resname', 'atom'))


def protein_atoms(prep):
    return {key(a['source']): np.array(a['xyz_A']) for a in prep['atom_graph']['source_to_qm']
            if a['kind']=='source' and a['source']['element'] not in ('H', 'D')
            and a['source']['canonical_resname']!='HOH'}


def fit(source, target):
    a, b = np.asarray(source), np.asarray(target)
    if a.shape != b.shape or len(a)<3 or np.linalg.matrix_rank(a-a.mean(axis=0))<2:
        raise InvalidArtifact('insufficient noncollinear source correspondences')
    ac, bc = a.mean(axis=0), b.mean(axis=0)
    u, _, vt = np.linalg.svd((a-ac).T@(b-bc))
    correction = np.eye(3); correction[-1,-1] = np.linalg.det(u@vt)
    rotation = u@correction@vt; translation = bc-ac@rotation
    rmsd = float(np.sqrt(np.mean(np.sum((a@rotation+translation-b)**2,axis=1))))
    return rotation, translation, rmsd


def actual_target_atoms(path):
    structure = gemmi.read_structure(str(path)); waters, other = [], []
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                if atom.element.name in ('H','D'):
                    continue
                data = {'chain':chain.name, 'resnum':residue.seqid.num,
                        'insertion_code':residue.seqid.icode.strip(), 'resname':residue.name,
                        'atom':atom.name, 'altloc':atom.altloc.replace('\x00',''),
                        'xyz_A':list(atom.pos), 'occupancy':atom.occ, 'B_iso':atom.b_iso}
                if residue.is_water():
                    if atom.element.name=='O':waters.append(data)
                elif atom.element.name not in ('Ca','La'):
                    other.append(data)
    return waters, other


def propose(states, agreement, output):
    verify(record(agreement))
    m = read_json(states); preparations = [read_json(verify(pin)) for pin in m['preparations']]
    if {p['case'] for p in preparations}!={'1F6S','6IP9'}:
        raise InvalidArtifact('only approved pair of real structures supported')
    transfers = []
    for source, target in (preparations, preparations[::-1]):
        a, b = protein_atoms(source), protein_atoms(target); keys = sorted(a.keys()&b.keys())
        rotation, translation, rmsd = fit([a[k] for k in keys], [b[k] for k in keys])
        waters, other = actual_target_atoms(verify(target['source_structure']))
        sites = []
        for w in source['water_groups']:
            if w['role']!='variable':continue
            coord = np.array(source['atoms'][w['oxygen_index']][1:])@rotation+translation
            nearest_water = min(waters,key=lambda a:np.linalg.norm(coord-a['xyz_A']))
            nearest_atom = min(other,key=lambda a:np.linalg.norm(coord-a['xyz_A']))
            wd = float(np.linalg.norm(coord-nearest_water['xyz_A']))
            ad = float(np.linalg.norm(coord-nearest_atom['xyz_A']))
            md = float(np.linalg.norm(coord-np.array(target['atoms'][0][1:])))
            reasons = (['near_target_heavy_atom'] if ad<2.0 else [])+(['outside_site_distance'] if md>3.2 else [])
            status = 'represented_by_target_water' if wd<=1.0 else ('flagged_candidate' if reasons else 'missing_site_candidate')
            sites.append({'source_water':w['source'], 'proposed_O_xyz_A':coord.tolist(), 'status':status,
                          'flags':reasons, 'nearest_target_water':nearest_water, 'nearest_water_distance_A':wd,
                          'nearest_target_nonmetal_heavy_atom':nearest_atom, 'nearest_heavy_atom_distance_A':ad,
                          'metal_distance_A':md})
        transfers.append({'source':source['case'], 'target':target['case'],
                          'fit_atom_keys':[list(k) for k in keys], 'fit_RMSD_A':rmsd,
                          'row_vector_rotation':rotation.tolist(), 'translation_A':translation.tolist(), 'sites':sites})
    result = {'states':record(states), 'agreement':record(agreement), 'implementation':record(__file__),
              'transfers':transfers, 'new_energy_evaluations':0, 'inserted_waters':0,
              'thresholds_A':{'represented_water':1.0,'heavy_atom_flag':2.0,'metal_site_distance':3.2},
              'occupancy_probabilities':None, 'baseline_changed':False}
    write_new(output,result)
    return {'transfers':len(transfers),'statuses':[s['status'] for t in transfers for s in t['sites']]}


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('states','agreement','output'):p.add_argument('--'+name,required=True)
    print(json.dumps(propose(**vars(p.parse_args())),indent=2))
