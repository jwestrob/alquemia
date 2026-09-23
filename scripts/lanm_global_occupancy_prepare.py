"""Prepare declared full-chain LanM occupancy hypotheses; no molecular energies."""
from __future__ import annotations
import argparse
import copy
from collections import Counter
from pathlib import Path
import time
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_global_prepare import protein, FF, WATER_FF, STANDARD
from mace_hybrid import write_xyz
from lanm_series_followup import heavy_atoms

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = 'LanM_whole_chain_A_pH5_fixed_crystal_water_occupancy_LaDy_v1'
OCCUPANCIES = {'EF12': ('EF1', 'EF2'), 'EF23': ('EF2', 'EF3'),
               'EF1234': ('EF1', 'EF2', 'EF3', 'EF4')}
SOURCES = {
    'Hans_8DQ2': ('Hans', 110, 'workspaces/benchmark_set_20260915/prepared/hans_lanm_v1/hans_protonation_manifest.json',
                  'workspaces/benchmark_set_20260915/sources/pdb/8DQ2.cif'),
    'Hans_8FNR': ('Hans', 110, 'workspaces/lanm_series_followup_20260923/dy_transfer_v1/hans_Dy_protonation_manifest.json',
                  'hans_lanm_protenix_benchmark/references/8FNR.pdb'),
    'Mex_8FNS': ('Mex', 105, 'workspaces/lanm_series_followup_20260923/sources_v2/mex_protonation_manifest.json',
                 'workspaces/benchmark_set_20260915/evidence/lanthanide/8FNS.cif')}


def source_identity(residue, atom):
    return {'chain': 'A', 'resid': residue.seqid.num,
            'insertion_code': residue.seqid.icode.strip(), 'resname': residue.name,
            'name': atom.name, 'element': atom.element.name, 'xyz_A': list(atom.pos)}


def mapped_protein(atoms, chain):
    """Link OpenMM canonical atom names back to the actual protonated source."""
    residues = {(str(r.seqid.num), r.seqid.icode.strip()): r for r in chain
                if r.name in STANDARD}
    result = []
    for i, a in enumerate(atoms):
        chain_id, resid, icode, name = a['id'].split('/')
        r = residues[(resid, icode)]
        original = a['source_xyz_A']
        matches = [s for s in r if s.element.name == a['element']
                   and np.max(np.abs(np.array(list(s.pos)) - original)) < 1e-8]
        if len(matches) != 1 or chain_id != 'A':
            raise InvalidArtifact('ambiguous protonated-source atom mapping: ' + a['id'])
        b = copy.deepcopy(a)
        b.update(index=i, source=source_identity(r, matches[0]))
        result.append(b)
    return result


def prepare_source(source_id, out):
    import gemmi
    protein_name, count, receipt_path, crystal_path = SOURCES[source_id]
    receipt_pin = record(ROOT / receipt_path)
    receipt = read_json(verify(receipt_pin))
    if (receipt['ph'] != 5.0 or receipt['add_missing_residues'] or
            receipt['repaired_missing_atom_count'] or receipt['repaired_missing_terminal_atom_count']):
        raise InvalidArtifact('source pH or missing-heavy-atom policy differs')
    pp = verify(receipt['output'])
    verify(receipt['source'])
    crystal = ROOT / crystal_path
    before, after = heavy_atoms(crystal), heavy_atoms(pp)
    selected_before = {k: v for k, v in before.items() if k[0] == 'A'}
    selected_after = {k: v for k, v in after.items() if k[0] == 'A'}
    if selected_before != selected_after:
        raise InvalidArtifact('chain A experimental heavy identities/coordinates changed')
    st = gemmi.read_structure(str(pp))
    chain = st[0]['A']
    residues = [r for r in chain if r.name in STANDARD]
    if len(residues) != count or any(r.seqid.icode.strip() for r in residues):
        raise InvalidArtifact('full source construct length/numbering differs')
    if [r.seqid.num for r in residues] != list(range(residues[0].seqid.num,
                                                    residues[0].seqid.num + count)):
        raise InvalidArtifact('missing residues inside selected protein chain')
    row = {'case_id': source_id, 'source_structure': receipt['output'],
           'selected_site': {'xyz_A': [0., 0., 0.]}}
    atoms, charge, audit = protein(row, allow_terminal_completion=False,
                                   check_peptide_connectivity=True)
    if audit['terminal_additions']:
        raise InvalidArtifact('undeclared heavy-atom completion')
    for a, old in zip(atoms, audit['original_protein_atoms']):
        a['source_xyz_A'] = old['xyz_A']
    atoms = mapped_protein(atoms, chain)
    indices = {a['id']: a['index'] for a in atoms}
    bonds = []
    for b in audit['bonds']:
        i, j = indices[b['atom_a_id']], indices[b['atom_b_id']]
        bonds.append({**b, 'indices': [i, j], 'kind': 'protein_covalent',
                      'source_distance_A': float(np.linalg.norm(np.array(atoms[i]['source_xyz_A']) -
                                                                atoms[j]['source_xyz_A'])),
                      'distance_A': float(np.linalg.norm(np.array(atoms[i]['xyz_A']) - atoms[j]['xyz_A']))})
    peptide = [b for b in bonds if atoms[b['indices'][0]]['source']['resid'] !=
               atoms[b['indices'][1]]['source']['resid']]
    if len(peptide) != count - 1:
        raise InvalidArtifact('protein does not have exactly the complete peptide graph')
    waters, sites, removed = [], {}, []
    for r in chain:
        if r.name in STANDARD:
            continue
        if r.name == 'HOH':
            if sorted(a.element.name for a in r) != ['H', 'H', 'O']:
                raise InvalidArtifact('source water lacks exactly O/H/H')
            o = next(a for a in r if a.element.name == 'O')
            ox = np.array(list(o.pos))
            wi = {'chain': 'A', 'resid': r.seqid.num, 'insertion_code': r.seqid.icode.strip(),
                  'resname': 'HOH', 'atom_indices': [], 'hydrogen_moves': []}
            for a in [o, *[a for a in r if a.element.name == 'H']]:
                original = np.array(list(a.pos)); position = original.copy()
                if a.element.name == 'H':
                    vector = original - ox; length = float(np.linalg.norm(vector))
                    if not .5 < length < 1.5:
                        raise InvalidArtifact('unsupported source water H bond')
                    position = ox + vector * (.9572 / length)
                    wi['hydrogen_moves'].append({'name': a.name, 'before_A': original.tolist(),
                                                 'after_A': position.tolist(), 'old_length_A': length})
                idx = len(atoms)
                atoms.append({'index': idx, 'id': f'A/{r.seqid.num}/{r.seqid.icode.strip()}/{a.name}',
                              'element': a.element.name, 'xyz_A': position.tolist(),
                              'source_xyz_A': original.tolist(), 'kind': 'retained_crystal_water',
                              'source': source_identity(r, a)})
                wi['atom_indices'].append(idx)
            oi, hi1, hi2 = wi['atom_indices']
            for hi in (hi1, hi2):
                bonds.append({'indices': [oi, hi], 'atom_a_id': atoms[oi]['id'],
                              'atom_b_id': atoms[hi]['id'], 'kind': 'water_covalent',
                              'forcefield_equilibrium_A': .9572, 'distance_A': .9572,
                              'source_distance_A': float(np.linalg.norm(np.array(atoms[hi]['source_xyz_A']) - ox))})
            vecs = [np.array(atoms[i]['xyz_A']) - ox for i in (hi1, hi2)]
            wi['angle_degrees'] = float(np.degrees(np.arccos(np.clip(np.dot(*vecs)/.9572**2, -1, 1))))
            waters.append(wi)
        elif r.seqid.num in (201, 202, 203, 204) and len(r) == 1 and r[0].element.name in ('La', 'Dy', 'Nd', 'Na'):
            a = r[0]; site = 'EF' + str(r.seqid.num - 200)
            sites[site] = {'site': site, 'source': source_identity(r, a), 'xyz_A': list(a.pos)}
        else:
            raise InvalidArtifact('unsupported chain A heterochemistry: ' + r.name + str(r.seqid))
    if set(sites) != {'EF1', 'EF2', 'EF3', 'EF4'}:
        raise InvalidArtifact('incomplete source site inventory')
    xyzs = np.array([a['xyz_A'] for a in atoms])
    minimum = float('inf'); closest = None
    for i in range(len(atoms)-1):
        distances = np.linalg.norm(xyzs[i+1:] - xyzs[i], axis=1)
        j = int(np.argmin(distances)) + i + 1
        if distances[j-i-1] < minimum:
            minimum = float(distances[j-i-1]); closest = [i, j]
    if minimum < .5:
        raise InvalidArtifact('atom overlap below 0.5 Angstrom after H preparation')
    # These are topology/template calls only, not an energy or minimization.
    from openmm import app
    pdb = app.PDBFile(str(pp)); mod = app.Modeller(pdb.topology, pdb.positions)
    mod.delete([a for a in mod.topology.atoms() if a.residue.chain.id != 'A' or a.residue.name not in STANDARD])
    templates = app.ForceField(str(FF)).getMatchingTemplates(mod.topology)
    out.mkdir(parents=True, exist_ok=False)
    write_new(out/'atoms.json', {'atoms': atoms, 'bonds': bonds})
    write_new(out/'waters.json', waters)
    write_new(out/'audit.json', audit)
    src = {'source_id': source_id, 'protein': protein_name, 'chain': 'A',
           'status': 'prepared', 'protein_residues': count, 'protein_atoms': audit['original_protein_atoms'].__len__(),
           'protein_charge': charge, 'base_atoms': len(atoms), 'water_count': len(waters),
           'crystal': record(crystal), 'protonation_receipt': receipt_pin, 'protonated_source': receipt['output'],
           'atoms': record(out/'atoms.json'), 'water_inventory': record(out/'waters.json'),
           'protein_audit': record(out/'audit.json'), 'sites': sites,
           'residue_templates': [{'resid': r.id, 'resname': r.name, 'template': t.name}
                                 for r, t in zip(mod.topology.residues(), templates)],
           'heavy_atoms_preserved': True, 'missing_heavy_added': 0, 'caps_added': 0,
           'peptide_bonds': len(peptide), 'minimum_all_atom_distance_A': minimum,
           'closest_atom_indices': closest,
           'water_HOH_angle_range_degrees': [min(w['angle_degrees'] for w in waters), max(w['angle_degrees'] for w in waters)],
           'assembly_limitation': 'chain A conditional monomer; no dimer/folding equilibrium',
           'occupancy_interpretation': 'fixed source geometry hypotheses; no natural population assertion'}
    write_new(out/'source.json', src)
    return src


def prepare_state(source, occupancy, parameters, output):
    data = read_json(verify(source['atoms']))
    atoms, bonds = copy.deepcopy(data['atoms']), copy.deepcopy(data['bonds'])
    selected = OCCUPANCIES[occupancy]
    metal_indices = []
    for site in selected:
        s = source['sites'][site]; idx = len(atoms); metal_indices.append(idx)
        atoms.append({'index': idx, 'id': 'metal/' + site, 'kind': 'occupied_lanthanide_site',
                      'element': 'M', 'site': site, 'xyz_A': s['xyz_A'], 'source': s['source'],
                      'source_xyz_A': s['xyz_A']})
    n = len(selected); charge = source['protein_charge'] + 3*n
    state_id = source['source_id'] + '__' + occupancy
    output.mkdir(parents=True, exist_ok=False)
    mapping = {'atoms': atoms, 'bonds': bonds, 'metal_indices': metal_indices,
               'coordination_bonds_included': False}
    write_new(output/'mapping.json', mapping)
    endpoints = {}
    import gemmi
    for metal in ('La', 'Dy'):
        rows = [(metal if a['element'] == 'M' else a['element'], *a['xyz_A']) for a in atoms]
        mult = 1 if metal == 'La' else 1 + 5*n
        electrons = sum(gemmi.Element(a[0]).atomic_number for a in rows) - charge
        valence = sum(sum(parameters['element'][a[0]]['refocc']) for a in rows) - charge
        p = parameters['element'][metal]
        if p['shells'] != ['5d', '6s', '6p'] or p['refocc'] != [1., 1., 1.]:
            raise InvalidArtifact('actual native f-in-core parameter identity differs')
        if electrons % 2 != (mult-1) % 2 or abs(valence-round(valence)) > 1e-9 or round(valence) % 2:
            raise InvalidArtifact('physical/native electron parity differs')
        path = output/(metal + '.xyz'); write_xyz(path, rows)
        endpoints[metal] = {'xyz': record(path), 'charge': charge,
                            'physical_multiplicity': mult, 'native_effective_multiplicity': 1,
                            'spin_multiplicity': mult, 'all_electron_count': electrons,
                            'native_valence_electron_count': round(valence),
                            'native_metal_shells': p['shells'], 'native_metal_refocc': p['refocc']}
    state = {'state_id': state_id, 'protein': source['protein'], 'source_id': source['source_id'],
             'status': 'prepared', 'n': n, 'occupied_sites': list(selected),
             'n_atoms': len(atoms), 'metal_indices': metal_indices, 'mapping': record(output/'mapping.json'),
             'protein_charge': source['protein_charge'], 'charge': charge,
             'source_water_inventory': source['water_inventory'], 'water_count': source['water_count'],
             'source_preparation': record(output.parent.parent/'sources'/source['source_id']/'source.json'),
             'endpoints': endpoints,
             'removed_source_ions': [s for k, s in source['sites'].items() if k not in selected],
             'EF4_Na_replacement_hypothesis': 'EF4' in selected and source['sites']['EF4']['source']['element'] == 'Na',
             'charge_pair_equal': True, 'coordinates_pair_equal': True,
             'Dy_spin_policy': 'maximum_spin_1_plus_5n; other couplings untested'}
    write_new(output/'state.json', state)
    return state


def prepare(output, agreement, sources=None):
    start = time.monotonic(); output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    prior = read_json(ROOT/'workspaces/lanm_series_followup_20260923/prepared_v2/manifest.json')
    parameter_pin = prior['parameter_export']; parameters = read_json(verify(parameter_pin))
    manifest = {'protocol_id': PROTOCOL, 'agreement': record(agreement),
                'preparation_agreement': record(ROOT/'diagnostics/lanm_global_occupancy_20260923/preparation/PLAN.md'),
                'implementation': {name: record(Path(__file__).with_name(name)) for name in
                    ('lanm_global_occupancy_prepare.py', 'mace_global_prepare.py', 'mace_hydrogen.py', 'lanm_series_followup.py')},
                'forcefields': {'protein': record(FF), 'water': record(WATER_FF)},
                'parameter_export': parameter_pin, 'sources': [], 'states': [],
                'molecular_calls': 0, 'new_protonations': 0, 'geometry_optimizations': 0}
    for source_id in sources or SOURCES:
        try:
            s = prepare_source(source_id, output/'sources'/source_id)
            manifest['sources'].append(s)
            for occ in OCCUPANCIES:
                state = prepare_state(s, occ, parameters, output/'states'/(source_id+'__'+occ))
                manifest['states'].append(state)
                print('PREPARED', state['state_id'], state['n_atoms'], state['charge'], flush=True)
        except Exception as exc:
            manifest['sources'].append({'source_id': source_id, 'status': 'unavailable', 'reason': str(exc)})
            for occ, sites in OCCUPANCIES.items():
                if not any(s['state_id'] == source_id+'__'+occ for s in manifest['states']):
                    manifest['states'].append({'state_id': source_id+'__'+occ, 'status': 'unavailable',
                                               'occupied_sites': list(sites), 'reason': str(exc)})
            print('UNAVAILABLE', source_id, str(exc), flush=True)
    manifest['elapsed_seconds'] = time.monotonic()-start
    write_new(output/'manifest.json', manifest)
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--agreement', required=True, type=Path)
    p.add_argument('--sources', nargs='+', choices=tuple(SOURCES))
    a = p.parse_args(); prepare(a.output, a.agreement, a.sources)


if __name__ == '__main__':
    main()
