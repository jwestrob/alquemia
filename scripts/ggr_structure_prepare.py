#!/usr/bin/env python3
"""Pinned source integrity and preparation for the approved GGR stage B only.

No energies are executed here. Source selection is fixed before scoring; the
existing residue-level alternate-conformer and donor policies are reused.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from pathlib import Path
import shutil
import resource
import time
import xml.etree.ElementTree as ET

import gemmi

import carve_generic as generic
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from affordable_peptide import SourceGraph
from coordination_policy import atom_altloc, selected_residue_atoms

PROTOCOL = 'ggr_source_conformation_audit_v1'
SOURCES = ('2FW0', '2FVY')
CHAIN = 'A'


def residue_id(residue):
    return (residue.seqid.num, residue.seqid.icode.strip(), residue.name)


def locator(residue):
    return {'chain': CHAIN, 'resnum': residue.seqid.num,
            'insertion_code': residue.seqid.icode.strip(), 'resname': residue.name,
            'label_seq_id': residue.label_seq}


def declared_sequence(structure):
    chain = structure[0][CHAIN]
    subchains = {r.subchain for r in chain if r.label_seq is not None}
    matches = [e for e in structure.entities if e.full_sequence and subchains.intersection(e.subchains)]
    if len(matches) != 1:
        raise InvalidArtifact('chain A must map to exactly one declared polymer sequence')
    return list(matches[0].full_sequence)


def observed_sequence_map(structure):
    result = {}
    for residue in structure[0][CHAIN]:
        if residue.label_seq is None:
            continue
        if residue.label_seq in result:
            raise InvalidArtifact('multiple residue identities at one polymer sequence position')
        result[residue.label_seq] = residue
    return result


def implementation_snapshot(output):
    from affordable_common import snapshot_implementation
    inventory = snapshot_implementation(output)
    additional = {}
    for name in ('ggr_structure_prepare.py', 'protonate_cif.py', 'ggr_preparation.py'):
        source = Path(__file__).parent / name
        if source.exists():
            destination = output / name
            shutil.copyfile(source, destination)
            additional[name] = {'live_source': record(source), 'preserved_copy': record(destination)}
    write_new(output / 'stage_b_additions.json', additional)
    return {'shared': inventory, 'stage_b': record(output / 'stage_b_additions.json')}


def inspect_source(source, reference, reference_repair, topology):
    """Audit deposited atoms; do not add atoms, place H, or build QM inputs."""
    structure = gemmi.read_structure(str(source))
    ref = gemmi.read_structure(str(reference))
    if len(structure) != 1 or [c.name for c in structure[0]].count(CHAIN) != 1:
        raise InvalidArtifact('approved single-model chain A source is absent or ambiguous')
    sequence, ref_sequence = declared_sequence(structure), declared_sequence(ref)
    if sequence != ref_sequence:
        raise InvalidArtifact('declared chain A sequence does not exactly match 1GLG')
    observed, ref_observed = observed_sequence_map(structure), observed_sequence_map(ref)
    for position, residue in observed.items():
        if not 1 <= position <= len(sequence) or residue.name != sequence[position - 1]:
            raise InvalidArtifact('observed residue disagrees with declared sequence')
    mapping = [{'sequence_position': pos, 'reference': locator(ref_observed[pos]),
                'source': locator(observed[pos]) if pos in observed else None}
               for pos in sorted(ref_observed)]
    ref_positions = {residue_id(r): pos for pos, r in ref_observed.items()}
    base = read_json(verify(read_json(reference_repair)['baseline']))
    reference_contacts = base['qm_inclusion']['typed_contacts']
    expected_positions = {ref_positions[(c['resnum'], c['insertion_code'], c['resname'])]
                          for c in reference_contacts if c['resname'] not in generic.WATER_NAMES}
    if len(expected_positions) != 6 or any(c['resname'] in generic.WATER_NAMES for c in reference_contacts):
        raise InvalidArtifact('reference no longer describes frozen six-residue, zero-water GGR')
    indexed, _ = generic.index_residues(structure[0])
    metals = [(r, a) for r in indexed if r.chain == CHAIN
              for a in selected_residue_atoms(r.residue) if a.element.name == 'Ca']
    if len(metals) != 1:
        raise InvalidArtifact('chain A must contain exactly one deposited Ca')
    metal_residue, metal = metals[0]
    contacts = generic.collect_contacts(indexed, metal, qm_inclusion_cut_A=3.3)
    qualified = [c for c in contacts if c.qualifies]
    included = [c for c in contacts if c.included_in_qm]
    position_by_id = {residue_id(r): pos for pos, r in observed.items()}
    actual_positions = {position_by_id.get((c.resnum, c.insertion_code, c.resname)) for c in included
                        if c.resname not in generic.WATER_NAMES}
    reasons = []
    if len(qualified) < generic.MIN_ORCA_COORDINATION_NUMBER or sum(c.element == 'N' for c in qualified) > generic.MAX_DIRECT_N:
        reasons.append('frozen typed coordination qualification failed')
    if actual_positions != expected_positions:
        reasons.append('included donor-residue membership changed')
    water_contacts = [c.to_dict() for c in included if c.resname in generic.WATER_NAMES]
    if water_contacts:
        reasons.append('frozen zero-water composition changed')
    unsupported = [c.to_dict() for c in contacts if c.distance_A <= 3.3
                   and generic.canonical_resname(c.resname) not in generic.STANDARD_AA
                   and c.resname not in generic.WATER_NAMES]
    if unsupported:
        reasons.append('unsupported species intersects frozen inclusion shell')

    graph = SourceGraph(source, topology)
    selected = set()
    for contact in included:
        key = contact.residue_key
        if contact.donor_type in ('sidechain_O', 'sidechain_N'):
            selected.update(graph.key(key, n) for n in generic.SIDECHAIN_QM_HEAVY_ATOMS[contact.resname])
        elif contact.donor_type == 'backbone_O':
            selected.update(graph.peptide(key))
            nitrogen = graph.amide_next.get(key)
            if nitrogen is None:
                raise InvalidArtifact('selected backbone carbonyl lacks a connected amide N')
            selected.update((graph.key(key, 'CA'), graph.key(nitrogen[:2], 'CA')))
    selected = graph.close_overlaps(selected)
    needed = selected | {n for k in selected for n in graph.required[k]}
    missing_local = [list(k) for k in needed if k not in graph.atoms]
    if missing_local:
        reasons.append('missing required local heavy atom')
    invalid_bonds = []
    for a, b in sorted(graph.edges):
        if a not in selected and b not in selected:
            continue
        aa, bb = graph.atoms[a], graph.atoms[b]
        distance = aa.pos.dist(bb.pos)
        if not .7 <= distance <= aa.element.covalent_r + bb.element.covalent_r + .35:
            invalid_bonds.append({'a': graph.meta[a], 'b': graph.meta[b], 'distance_A': distance})
    if invalid_bonds:
        reasons.append('invalid local source connectivity')
    local_keys = {k[:2] for k in needed}
    template_atoms = {r.attrib['name']: {a.attrib['name'] for a in r.findall('Atom')
                                      if not a.attrib['name'].startswith('H')}
                      for r in ET.parse(topology).findall('.//Residues/Residue')}
    missing_standard = []
    altlocs, heterogens, waters = [], [], []
    selected_core = [graph.atoms[k] for k in selected]
    for indexed_residue in indexed:
        residue = indexed_residue.residue
        atoms = selected_residue_atoms(residue)
        heavy = [a for a in atoms if not generic.is_hydrogen(a)]
        if not heavy:
            continue
        minimum = min(a.pos.dist(metal.pos) for a in heavy)
        info = {**locator(residue), 'nearest_metal_distance_A': minimum,
                'nearest_selected_core_heavy_distance_A': min(a.pos.dist(b.pos) for a in heavy for b in selected_core)}
        labels = sorted({atom_altloc(a) for a in residue} - {''})
        if labels:
            altlocs.append({**info, 'in_local_representation_or_cap_neighbours': indexed_residue.key in local_keys,
                           'available_labels': labels,
                           'selected_atoms': [{'atom': a.name, 'altloc': atom_altloc(a), 'occupancy': a.occ,
                                               'xyz_A': list(a.pos)} for a in atoms]})
        expected = template_atoms.get(residue.name)
        if expected is None and generic.canonical_resname(residue.name) == 'HIS':
            expected = template_atoms.get('HID')
        if expected is not None:
            missing = sorted(expected - {a.name.strip() for a in heavy})
            if missing:
                missing_standard.append({**info, 'missing_atoms': missing,
                                         'in_local_representation_or_cap_neighbours': indexed_residue.key in local_keys})
                if indexed_residue.key in local_keys:
                    reasons.append('incomplete local standard residue')
        elif residue.name in generic.WATER_NAMES:
            waters.append(info)
        elif residue.name != 'CA':
            heterogens.append({**info, 'atom_names': [a.name for a in heavy]})
    block = gemmi.cif.read(str(source)).sole_block()
    modifications = [dict(zip(list(block.find_mmcif_category('_pdbx_struct_mod_residue.').tags), row))
                     for row in block.find_mmcif_category('_pdbx_struct_mod_residue.')]
    sequence_text = ','.join(sequence)
    return {'status': 'source_integrity_pass' if not reasons else 'unsupported',
            'unsupported_reasons': sorted(set(reasons)), 'source': record(source),
            'reference': record(reference), 'reference_repair': record(reference_repair),
            'chain': CHAIN, 'assembly': 'sequence-matched author chain A monomer, no symmetry mates',
            'sequence_identity': 'exact declared 309-residue identity to 1GLG',
            'sequence_length': len(sequence),
            'sequence_sha256': hashlib.sha256(sequence_text.encode()).hexdigest(),
            'sequence_mapping': mapping,
            'missing_sequence_positions': sorted(set(range(1, len(sequence)+1)) - observed.keys()),
            'site_selector': {**metal_residue.source_dict(), 'atom': metal.name.strip(),
                              'xyz_A': list(metal.pos), 'occupancy': metal.occ},
            'qualification_cutoff_A': 3.1, 'inclusion_cutoff_A': 3.3,
            'diagnostic_cutoff_A': generic.DIAGNOSTIC_CUTOFF_A,
            'altloc_policy': generic.ALTLOC_POLICY_ID,
            'typed_qualification': [c.to_dict() for c in qualified],
            'typed_inclusion': [c.to_dict() for c in included],
            'diagnostic_contacts': [c.to_dict() for c in contacts],
            'expected_donor_sequence_positions': sorted(expected_positions),
            'actual_donor_sequence_positions': sorted(actual_positions, key=lambda x: -1 if x is None else x),
            'explicit_water_inventory': water_contacts, 'unsupported_shell_species': unsupported,
            'local_source_atoms_and_cap_neighbours': [graph.meta[k] for k in sorted(needed) if k in graph.meta],
            'local_missing_heavy_atoms': missing_local, 'invalid_local_bonds': invalid_bonds,
            'missing_standard_heavy_atoms': missing_standard,
            'alternate_conformers': altlocs, 'nonstandard_species': heterogens,
            'deposited_water_inventory': waters, 'deposited_modifications': modifications,
            'radiation_damage_caveat': ('Reported partial Glu149 decarboxylation; retained deposited atoms and frozen altloc policy'
                                        if source.stem.upper() == '2FVY' else None),
            'interpretation': 'structural robustness; one GGR observation, no isolated sugar-effect attribution',
            'new_high_level_evaluations': 0}


def audit(sources, reference, reference_repair, agreement, output):
    if output.exists():
        raise InvalidArtifact('refusing existing audit output')
    output.mkdir(parents=True)
    start, cpu = time.monotonic(), time.process_time()
    repair = read_json(reference_repair)
    topology = verify(repair['topology_definition'])
    report = {'protocol_id': PROTOCOL, 'status': 'audited', 'agreement': record(agreement),
              'reference': record(reference), 'reference_repair': record(reference_repair),
              'topology': record(topology), 'cases': [],
              'implementation_snapshot': implementation_snapshot(output / 'implementation')}
    for sid, source in zip(SOURCES, sources):
        try:
            if gemmi.cif.read(str(source)).sole_block().name.upper() != sid:
                raise InvalidArtifact('deposited source ID does not match the preselected case')
            row = inspect_source(source, reference, reference_repair, topology)
        except Exception as exc:
            row = {'status': 'unsupported', 'source': record(source),
                   'unsupported_reasons': [f'{type(exc).__name__}: {exc}']}
        row['pdb_id'] = sid
        report['cases'].append(row)
    report['accounting'] = {'wall_seconds': time.monotonic()-start,
                            'process_cpu_seconds': time.process_time()-cpu,
                            'high_level_evaluations': 0}
    write_new(output / 'source_audit.json', report)
    return report


def selected_monomer(source, output):
    """Collapse altlocs by the frozen policy before PDBFixer reads the monomer."""
    structure = gemmi.read_structure(str(source))
    if len(structure) != 1:
        raise InvalidArtifact('one deposited model is required')
    for i in reversed(range(len(structure[0]))):
        if structure[0][i].name != CHAIN:
            del structure[0][i]
    selection = []
    for residue in structure[0][CHAIN]:
        chosen = [a.clone() for a in selected_residue_atoms(residue)]
        selection.extend({**locator(residue), 'atom': a.name.strip(),
                          'source_altloc': atom_altloc(a), 'element': a.element.name,
                          'occupancy': a.occ, 'xyz_A': list(a.pos)} for a in chosen)
        for i in reversed(range(len(residue))):
            del residue[i]
        for atom in chosen:
            # A single physical conformer is now present; retain its original
            # alternate identity in the explicit source-selection map above.
            atom.altloc = '\x00'
            residue.add_atom(atom)
    structure.make_mmcif_document().write_file(str(output))
    return selection


def heavy_coordinate_check(selection, protonated):
    expected = {(a['chain'], a['resnum'], a['insertion_code'], a['atom']): a
                for a in selection if a['element'] not in ('H', 'D')}
    actual = {}
    structure = gemmi.read_structure(str(protonated))
    for chain in structure[0]:
        for residue in chain:
            for atom in selected_residue_atoms(residue):
                if generic.is_hydrogen(atom):
                    continue
                key = (chain.name, residue.seqid.num, residue.seqid.icode.strip(), atom.name.strip())
                if key in actual:
                    raise InvalidArtifact('duplicate prepared heavy-atom identifier')
                actual[key] = list(atom.pos)
    missing = [list(k) for k in expected if k not in actual]
    if missing:
        raise InvalidArtifact(f'protonation lost source heavy atoms: {missing}')
    maximum = max((sum((x-y)**2 for x,y in zip(a['xyz_A'], actual[k]))**.5
                   for k,a in expected.items()), default=0.)
    if maximum > .002:
        raise InvalidArtifact('source heavy coordinates moved during protonation')
    return {'source_heavy_atom_count': len(expected), 'maximum_displacement_A': maximum,
            'serialization_tolerance_A': .002,
            'added_heavy_atom_identifiers': [list(k) for k in actual if k not in expected]}


def prepare(source_audit, stage_a_ready, output):
    """One pH-7 source per case, shared by repaired v3 and alpha-capped inputs."""
    from protonate_cif import protonate
    from affordable_peptide import repair
    from ggr_preparation import prepare as prepare_boundary
    if output.exists():
        raise InvalidArtifact('refusing existing preparation directory')
    root = Path(__file__).resolve().parents[1]
    if not output.resolve().is_relative_to((root/'workspaces').resolve()):
        raise InvalidArtifact('scientific products must be under workspaces/')
    prior = read_json(stage_a_ready)
    if prior.get('status') != 'prepared' or len(prior.get('preparations', [])) != 5:
        raise InvalidArtifact('Stage A preparation readiness is not established')
    for item in prior['preparations']:
        verify(item['manifest'])
    audit_record = read_json(source_audit)
    if audit_record['protocol_id'] != PROTOCOL:
        raise InvalidArtifact('unsupported source audit')
    verify(audit_record['agreement'])
    topology = verify(audit_record['topology'])
    output.mkdir(parents=True)
    report = {'schema_version': 'alquemia.ggr_stage_b_preparation.v1',
              'status': 'prepared', 'source_audit': record(source_audit),
              'agreement': audit_record['agreement'], 'stage_a_ready': record(stage_a_ready),
              'implementation_snapshot': implementation_snapshot(output/'implementation'),
              'cases': [], 'preparations': [], 'high_level_evaluations_executed': 0}
    for audited in audit_record['cases']:
        sid = audited['pdb_id']
        row = {'pdb_id': sid, 'status': 'unsupported', 'source': audited['source']}
        report['cases'].append(row)
        if audited['status'] != 'source_integrity_pass':
            row['reason'] = audited['unsupported_reasons']
            continue
        source = verify(audited['source'])
        directory = output/sid.lower()
        directory.mkdir()
        start, cpu = time.monotonic(), time.process_time()
        try:
            selected = directory/'source_selected_chain_A.cif'
            selection = selected_monomer(source, selected)
            selection_path = directory/'source_selection.json'
            write_new(selection_path, {'original_source': record(source),
                      'selected_monomer': record(selected), 'atoms': selection,
                      'altloc_policy': generic.ALTLOC_POLICY_ID,
                      'note': 'Original selected heavy coordinates and all chain-A heterogens/waters retained.'})
            row['source_selection'] = record(selection_path)
            protonated = directory/'source_protonated.pdb'
            protonation_start, protonation_cpu = time.monotonic(), time.process_time()
            row['protonation'] = protonate(selected, protonated, ph=7.0, add_missing_residues=False)
            row['protonation_accounting'] = {'wall_seconds': time.monotonic()-protonation_start,
                                            'process_cpu_seconds': time.process_time()-protonation_cpu}
            row['heavy_coordinate_audit'] = heavy_coordinate_check(selection, protonated)
            # The only permitted added heavy atoms are remote terminal atoms,
            # never a reconstruction inside either selected representation.
            local = {(x['chain'], x['resnum'], x['insertion_code'])
                     for x in audited['local_source_atoms_and_cap_neighbours']}
            if any(tuple(x[:3]) in local for x in row['heavy_coordinate_audit']['added_heavy_atom_identifiers']):
                raise InvalidArtifact('protonation reconstructed a local heavy atom')
            site = audited['site_selector']
            stem = 'ggr_'+sid.lower()+'_GGR'
            base_dir = directory/'historical_v2_intermediate'
            baseline = generic.carve(protonated, base_dir, stem,
                        site_chain=CHAIN, site_resnum=site['resnum'],
                        site_icode=site['insertion_code'], site_atom=site['atom'],
                        qm_inclusion_cut=3.3)
            baseline_path = base_dir/(stem+'_carve_manifest.json')
            row['baseline_manifest'] = record(baseline_path)
            if baseline['status'] != 'ready_for_orca':
                raise InvalidArtifact('frozen protonated source fails generic qualification')
            expected = {(x['chain'], x['resnum'], x['insertion_code'], x['resname'], x['donor_type'])
                        for x in audited['typed_inclusion']}
            observed = {(x['chain'], x['resnum'], x['insertion_code'], x['resname'], x['donor_type'])
                        for x in baseline['qm_inclusion']['typed_contacts']}
            if observed != expected:
                raise InvalidArtifact('protonation changed frozen donor membership or donor type')
            if baseline['charge_ledger']['fragment_formal_charge_sum'] != -3:
                raise InvalidArtifact('incompatible GGR protonation/charge state')
            repaired_dir = directory/'amide_v3'
            repaired = repair(baseline_path, repaired_dir, topology)
            if repaired['explicit_water_inventory'] or repaired['outputs']['La']['charge'] != 0 or repaired['outputs']['Ca']['charge'] != -1:
                raise InvalidArtifact('repaired GGR water/charge invariant failed')
            repair_path = repaired_dir/'repair_manifest.json'
            extended_dir = directory/'alpha_caps'
            extended = prepare_boundary(repair_path, extended_dir, 'alpha_caps')
            manifests = [('amide_v3', repair_path, repaired),
                         ('alpha_caps', extended_dir/'preparation_manifest.json', extended)]
            for representation, manifest, prepared in manifests:
                report['preparations'].append({'name': 'ggr_'+sid.lower()+'_'+representation,
                         'pdb_id': sid, 'representation': representation,
                         'manifest': record(manifest), 'protocol_id': prepared['protocol_id'],
                         'atom_count': prepared['paired_invariants']['atom_count'],
                         'source_coordinate_max_displacement_A': (prepared['heavy_coordinate_max_displacement_A']
                           if representation == 'amide_v3' else prepared['source_coordinate_max_displacement_A']),
                         'ligand_formal_charge': -3, 'new_endpoint_count': 2,
                         'source_protonation': record(protonated),
                         'source_audit': record(source_audit),
                         'biological_group': 'GGR_1GLG', 'evaluation_role': 'method_development'})
            row['status'] = 'prepared'
            row['prepared_manifests'] = [record(m) for _,m,_ in manifests]
        except Exception as exc:
            row.update(status='unsupported', reason=f'{type(exc).__name__}: {exc}')
        row['accounting'] = {'wall_seconds': time.monotonic()-start,
                             'process_cpu_seconds': time.process_time()-cpu,
                             'process_peak_rss_KiB_so_far': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                             'high_level_evaluations': 0}
        write_new(directory/'preparation_report.json', row)
    if any(row['status'] != 'prepared' for row in report['cases']):
        report['status'] = 'partial_or_unsupported'
    write_new(output/'inventory.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='operation', required=True)
    a = sub.add_parser('audit')
    a.add_argument('--source-2fw0', type=Path, required=True)
    a.add_argument('--source-2fvy', type=Path, required=True)
    a.add_argument('--reference', type=Path, required=True)
    a.add_argument('--reference-repair', type=Path, required=True)
    a.add_argument('--agreement', type=Path, required=True)
    a.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--source-audit', type=Path, required=True)
    p.add_argument('--stage-a-ready', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.operation == 'audit':
        result = audit([args.source_2fw0, args.source_2fvy], args.reference,
                       args.reference_repair, args.agreement, args.output)
    else:
        result = prepare(args.source_audit, args.stage_a_ready, args.output)
    for row in result['cases']:
        print(row['pdb_id'], row['status'], row.get('reason', '; '.join(row.get('unsupported_reasons', []))))


if __name__ == '__main__':
    main()
