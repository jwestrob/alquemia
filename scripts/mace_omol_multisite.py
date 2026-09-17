"""Source-backed whole-chain site substitutions with fixed background calcium."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_global_prepare import protein, atom_id, STANDARD, FF, WATER_FF
from mace_gb import RADII
from mace_hybrid import check_atoms, write_xyz

POLICY = 'omol_intact_multisite_fixed_background_Ca_source_acetyl_ff19sb_v1'
ASSEMBLY = 'deposited_chain_A_with_declared_cofactor_and_site_waters'


def build(recipe, site_key):
    import gemmi
    from openmm import app, unit
    cfg = read_json(recipe); panel = read_json(verify(cfg['panel']))
    choices = [r for r in panel['proteins'] if r['panel_key'] == cfg['panel_key']]
    if len(choices) != 1: raise InvalidArtifact('unique frozen source panel entry required')
    source = choices[0]; sites = source['sites']
    if [s['site_key'] for s in sites] != cfg['site_order'] or site_key not in cfg['site_order']:
        raise InvalidArtifact('source site inventory/order differs')
    pm = read_json(verify(cfg['protonation']))
    if (pm['source'] != source['source'] or pm['ph'] != 7.0 or pm['add_missing_residues']
            or pm.get('repaired_missing_atom_count', 0)):
        raise InvalidArtifact('requires archived pH7 preparation without rebuilt residues/heavy atoms')
    raw = gemmi.read_structure(str(verify(source['source']))); raw.setup_entities()
    residues = [r for c in raw[0] if c.name == 'A' for r in c]
    supported = STANDARD | ({'ACE'} if cfg['retain_source_acetyl'] else set())
    inventory = []; raw_atoms = {}; raw_metals = []; raw_waters = []
    for r in residues:
        if r.name not in supported | {'CA','HOH'}:
            raise InvalidArtifact('unsupported raw multisite source chemistry: '+r.name+str(r.seqid))
        if r.name in supported and r.entity_type != gemmi.EntityType.Polymer:
            raise InvalidArtifact('expected source polymer, including actual acetyl when requested')
        inventory.append({'resname': r.name, 'resid': str(r.seqid), 'atoms': len(r)})
        for a in r:
            if a.element.is_hydrogen: continue
            key = f'A/{r.seqid.num}/{r.seqid.icode.strip()}/{a.name}'
            raw_atoms.setdefault(key, []).append({'element': a.element.name, 'xyz_A': list(a.pos),
                                                  'altloc': a.altloc.strip('\x00 '), 'occupancy': a.occ})
            if r.name == 'CA':
                if len(r) != 1 or a.element.name != 'Ca': raise InvalidArtifact('malformed source calcium')
                raw_metals.append(key)
            elif r.name == 'HOH': raw_waters.append(key)
    acetyls = [r for r in residues if r.name == 'ACE']
    if len(acetyls) != int(cfg['retain_source_acetyl']):
        raise InvalidArtifact('actual source acetyl state differs from recipe')
    pdb = app.PDBFile(str(verify(pm['output'])))
    positions = np.array(pdb.positions.value_in_unit(unit.angstrom))
    atoms = [a for a in pdb.topology.atoms() if a.residue.chain.id == 'A']
    byid = {atom_id(a): a for a in atoms}
    if len(byid) != len(atoms): raise InvalidArtifact('duplicate protonated source atom identity')
    metals = {}
    for s in sites:
        sel = s['metal_selector']
        if sel['chain'] != 'A' or sel['resname'] != 'CA': raise InvalidArtifact('background must be actual Ca')
        key = f"A/{sel['resnum']}/{sel['icode']}/{sel['atom']}"
        if key not in byid or byid[key].element.symbol != 'Ca': raise InvalidArtifact('missing source calcium')
        pos = positions[byid[key].index]
        if not np.allclose(pos, s['source_metal_xyz_A'], atol=1e-6, rtol=0):
            raise InvalidArtifact('selected source metal coordinates differ')
        metals[s['site_key']] = {'id': key, 'element': 'Ca', 'xyz_A': pos.tolist(),
                                  'radius_A': 1.8, 'kind': 'background_metal', 'formal_charge_e': 2}
    if sorted(raw_metals) != sorted(a['id'] for a in metals.values()):
        raise InvalidArtifact('source has omitted/extra metal sites')
    water_ids = sorted({f"{w['chain']}/{w['resnum']}//{w['atom']}"
                        for s in sites for w in s['first_shell_waters']})
    if water_ids != cfg['retained_water_oxygen_ids']:
        raise InvalidArtifact('fixed union of archived site waters differs')
    waters = []; water_atoms = []; moves = []
    for key in water_ids:
        if key not in byid or key not in raw_waters: raise InvalidArtifact('missing actual retained water')
        oxygen = byid[key]; members = list(oxygen.residue.atoms())
        if sorted(a.element.symbol for a in members) != ['H','H','O']:
            raise InvalidArtifact('retained water lacks complete archived O/H/H chemistry')
        op = positions[oxygen.index]
        waters.append({'oxygen_id': key, 'source_oxygen': raw_atoms[key]})
        for a in members:
            coord = positions[a.index].copy()
            if a.element.symbol == 'H':
                length = float(np.linalg.norm(coord-op))
                if not .5 < length < 1.5: raise InvalidArtifact('unsupported water O-H geometry')
                after = op+(coord-op)*(.9572/length)
                moves.append({'id': atom_id(a), 'before_A': coord.tolist(), 'after_A': after.tolist(),
                              'old_length_A': length, 'new_length_A': .9572})
                coord = after
            water_atoms.append({'id': atom_id(a), 'element': a.element.symbol, 'xyz_A': coord.tolist(),
                                'radius_A': RADII[a.element.symbol], 'kind': 'retained_site_water'})
    row = {'case_id': cfg['case_prefix']+'_'+site_key, 'source_structure': pm['output'],
           'selected_site': {'xyz_A': metals[site_key]['xyz_A']}}
    pol, q, details = protein(row, allow_terminal_completion=True, check_peptide_connectivity=True,
                              retain_source_acetyl=cfg['retain_source_acetyl'])
    background = [metals[s['site_key']] for s in sites if s['site_key'] != site_key]
    selected = {**metals[site_key], 'element': 'M', 'kind': 'selected_metal', 'formal_charge_e': None}
    physical = [*pol, *water_atoms, *background, selected]
    if len({a['id'] for a in physical}) != len(physical): raise InvalidArtifact('overlapping multisite atoms')
    matched = {}; seen_heavy = set()
    for a in physical:
        if a['element'] == 'H' or a['kind'] == 'terminal_completion': continue
        elem = 'Ca' if a['element'] == 'M' else a['element']
        matches = [r for r in raw_atoms.get(a['id'], []) if r['element'] == elem
                   and np.allclose(r['xyz_A'], a['xyz_A'], atol=1e-6, rtol=0)]
        if len(matches) != 1: raise InvalidArtifact('missing/ambiguous raw source heavy mapping: '+a['id'])
        matched[a['id']] = matches[0]; seen_heavy.add(a['id'])
    expected_heavy = {f'A/{r.seqid.num}/{r.seqid.icode.strip()}/{a.name}' for r in residues
                      if r.name in supported | {'CA'} for a in r if not a.element.is_hydrogen} | set(water_ids)
    if seen_heavy != expected_heavy: raise InvalidArtifact('raw protein/metal heavy atoms were omitted or added')
    links = []
    for b in raw.connections:
        if b.type != gemmi.ConnectionType.Covale or not (b.partner1.chain_name == 'A' or b.partner2.chain_name == 'A'):
            continue
        addresses = []
        for a in (b.partner1, b.partner2):
            addresses.append(f'{a.chain_name}/{a.res_id.seqid.num}/{a.res_id.seqid.icode.strip()}/{a.atom_name}')
        retained = {frozenset((b['atom_a_id'], b['atom_b_id'])) for b in details['bonds']}
        if frozenset(addresses) not in retained: raise InvalidArtifact('raw covalent connection was lost')
        links.append(addresses)
    evidence = cfg['evidence']
    result = {'status': 'prepared', 'case_id': row['case_id'], 'policy_id': POLICY,
              'source': pm['output'], 'source_preparation': cfg['protonation'], 'recipe': record(recipe),
              'selected_site_key': site_key, 'physical_atoms': physical, 'assembly': ASSEMBLY,
              'explicit_waters': waters, 'water_H_moves': moves, 'protein_charge_e': q,
              'cofactor_charge_e': 0, 'background_metal_charge_e': 2*len(background),
              'background_metals': background, 'forcefield': record(FF), 'water_forcefield': record(WATER_FF),
              'preparation_details': details, 'raw_source': source['source'],
              'source_heavy_mapping': matched, 'source_covalent_connections': links,
              'raw_residue_inventory': inventory, 'excluded_water_oxygen_ids': sorted(set(raw_waters)-set(water_ids)),
              'reported_missing_residues': pm.get('detected_missing_residues'),
              'microstate': 'archived_explicit_protonation; actual_source_acetyl_if_present; background_Ca2plus; fixed_neutral_site_water_union',
              'evidence': evidence, 'evidence_use': 'consumed_method_development', 'reference': None,
              'baseline_coordinates_modified': False}
    return result


@cached_file_checks
def audit(path):
    p = read_json(path)
    for pin in p['implementation'].values(): verify(pin)
    verify(p['agreement'])
    expected = build(verify(p['recipe']), p['selected_site_key'])
    if {k:v for k,v in p.items() if k not in ('implementation','agreement','endpoints')} != expected:
        raise InvalidArtifact('multisite preparation differs from full source replay')
    for metal in ('La','Ca'):
        e = p['endpoints'][metal]
        coords = [(metal if a['element']=='M' else a['element'], *a['xyz_A']) for a in p['physical_atoms']]
        q = p['protein_charge_e']+p['background_metal_charge_e']+(3 if metal=='La' else 2)
        background_indices = tuple(i for i,a in enumerate(p['physical_atoms']) if a['kind']=='background_metal')
        if xyz(verify(e['xyz'])) != coords or e['charge'] != q or e['spin_multiplicity'] != 1 or e['state'] != check_atoms(coords,q,background_calcium_indices=background_indices):
            raise InvalidArtifact('multisite paired coordinates/charge/electron state differ')
    return {'status': 'pass', 'preparation': record(path), 'atoms': len(p['physical_atoms']),
            'policy_id': POLICY, 'assembly': ASSEMBLY, 'protein_charge_e': p['protein_charge_e'],
            'background_metal_charge_e': p['background_metal_charge_e'],
            'explicit_water_count': len(p['explicit_waters']), 'background_metals': p['background_metals'],
            'scope': 'actual source replay; selected substitution with fixed other calcium ions and site-water union',
            'new_model_forwards': 0}


@cached_file_checks
def prepare(recipe, site_key, agreement, output):
    p = build(recipe, site_key); out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    from mace_global_benchmark import snapshot
    names = tuple(f.name for f in Path(__file__).parent.glob('mace_omol*.py')) + (
        'mace_canonical_run.py','mace_canonical.py','mace_curvature.py','affordable_response.py',
        'mace_mechanics_run.py','mace_short_engine.py','mace_file_checks.py')
    p['implementation'] = snapshot(out, names); p['agreement'] = record(agreement); p['endpoints'] = {}
    for metal in ('La','Ca'):
        coords = [(metal if a['element']=='M' else a['element'], *a['xyz_A']) for a in p['physical_atoms']]
        q = p['protein_charge_e']+p['background_metal_charge_e']+(3 if metal=='La' else 2)
        background_indices = tuple(i for i,a in enumerate(p['physical_atoms']) if a['kind']=='background_metal')
        state = check_atoms(coords,q,background_calcium_indices=background_indices)
        dest = out/(metal+'.xyz'); write_xyz(dest,coords)
        p['endpoints'][metal] = {'xyz': record(dest), 'charge': q, 'spin_multiplicity': 1, 'state': state}
    final = out/'preparation.json'; write_new(final,p)
    checked = audit(final); write_new(out/'audit.json',checked)
    return {'status': 'prepared', 'preparation': record(final), 'audit': checked}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('recipe','site-key','agreement','output'): parser.add_argument('--'+key,required=True)
    print(json.dumps(prepare(**vars(parser.parse_args())), indent=2))
