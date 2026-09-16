#!/usr/bin/env python3
"""Approved GGR boundary diagnostics from pinned, already protonated v3 cores.

No donor reselection, protonation, coordinate optimization, or job submission.
The generic production carver and historical inputs are never modified.
"""
from __future__ import annotations

import argparse
import copy
import importlib.metadata
import platform
import shutil
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

import affordable_peptide as peptide
import carve_generic as generic
from affordable_common import (InvalidArtifact, cache_key, paired, read_json,
                               record, verify, write_new, xyz)

PROTOCOLS = {
    'alpha_caps': 'generic_peptide_alpha_caps_native_r2scan3c_dev_v1',
    'connected_segment': 'ggr_connected_segment_native_r2scan3c_dev_v1',
}
POLICIES = {
    'alpha_caps': 'source_graph_amide_alpha_sigma_caps_v1',
    'connected_segment': 'source_graph_ggr_nearest_donor_native_segment_v1',
}
RECIPE = '! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3'
DEPENDENCIES = ('ggr_preparation.py', 'affordable_peptide.py', 'affordable_common.py',
                'carve_generic.py', 'coordination_policy.py')


def _source_key(item):
    s = item['source']
    return s['chain_index'], s['residue_index'], s['atom']


def _source_selection(graph, parent):
    waters = {graph.locate(s).key for s in parent['explicit_water_inventory']}
    selected = set()
    water_items = []
    keys = set()
    for item in parent['atom_graph']['source_to_qm']:
        if item['kind'] != 'source':
            continue
        key = _source_key(item)
        if key in keys or key not in graph.atoms:
            raise InvalidArtifact(f'duplicate/missing parent source mapping: {key}')
        keys.add(key)
        if list(graph.atoms[key].pos) != item['xyz_A']:
            raise InvalidArtifact(f'parent mapping/source coordinate mismatch: {key}')
        if key[:2] in waters:
            water_items.append(copy.deepcopy(item))
        elif item['source']['element'] not in ('H', 'D'):
            selected.add(key)
    for entry in parent['charge_ledger']:
        if entry['kind'] == 'sidechain':
            r = graph.locate(entry['source'])
            _, check = generic.sidechain_fragment(r, [], 0)
            if check['formal_charge'] != entry['formal_charge']:
                raise InvalidArtifact(f'parent side-chain microstate changed: {r.source_dict()}')
    return selected, water_items


def _materialize(graph, selected, water_items):
    atoms, mapping = graph.materialize(selected)
    for item in water_items:
        atom = generic.atom_to_qm(graph.atoms[_source_key(item)])
        atoms.append(atom)
        item = copy.deepcopy(item)
        item['qm_index'] = len(atoms)
        mapping['source_to_qm'].append(item)
    source_keys = [_source_key(a) for a in mapping['source_to_qm'] if a['kind'] == 'source']
    if len(source_keys) != len(set(source_keys)):
        raise InvalidArtifact('duplicate source mapping after fragment merge')
    return atoms, mapping


def expand_selection(graph, selected, ledger, policy):
    """Return a new heavy-atom set and a traceable, topology-derived extension."""
    if policy not in PROTOCOLS:
        raise InvalidArtifact(f'unknown preparation policy: {policy}')
    amides = [graph.locate(e['source']).key for e in ledger
              if e['kind'] == 'backbone_carbonyl']
    if not amides:
        raise InvalidArtifact('no backbone carbonyl units to extend')
    extended = set(selected)
    description = {'amide_units': [], 'new_formal_charge': 0}
    for rkey in amides:
        n = graph.amide_next.get(rkey)
        if n is None or graph.residues[n[:2]].canonical_resname == 'PRO':
            raise InvalidArtifact('alpha-cap policy requires a supported nonproline peptide N')
        extended.update((graph.key(rkey, 'CA'), graph.key(n[:2], 'CA')))
        description['amide_units'].append({
            'carbonyl': graph.residues[rkey].source_dict(),
            'nitrogen': graph.meta[n],
            'alpha_anchors': [graph.meta[graph.key(rkey, 'CA')],
                              graph.meta[graph.key(n[:2], 'CA')]],
        })
    if policy == 'alpha_caps':
        return extended, description
    if len(amides) != 1:
        raise InvalidArtifact('connected GGR policy requires exactly one selected peptide carbonyl')
    donors = {graph.locate(e['source']).key for e in ledger if e['kind'] == 'sidechain'}
    center = amides[0]
    left = center
    visited = set()
    while left not in donors:
        if left in visited or left not in graph.amide_prev:
            raise InvalidArtifact('no connected preceding selected side-chain donor')
        visited.add(left)
        left = graph.amide_prev[left][:2]
    right = graph.amide_next[center][:2]
    visited = set()
    while right not in donors:
        if right in visited or right not in graph.amide_next:
            raise InvalidArtifact('no connected following selected side-chain donor')
        visited.add(right)
        right = graph.amide_next[right][:2]
    segment = [left]
    while segment[-1] != right:
        n = graph.amide_next.get(segment[-1])
        if n is None or n[:2] in segment:
            raise InvalidArtifact('broken/ambiguous connected segment')
        segment.append(n[:2])
    # This is the approved GGR sequence-defined diagnostic, not an automatic
    # whole-panel region-expansion or charge-changing policy.
    sequence = [graph.residues[k].canonical_resname for k in segment]
    if sequence != ['ASP', 'GLY', 'GLN', 'ILE', 'GLN']:
        raise InvalidArtifact(f'connected segment is outside approved GGR chemistry: {sequence}')
    if graph.residues[center].canonical_resname != 'GLN' or center != segment[2]:
        raise InvalidArtifact('connected segment has a different backbone donor position')
    if any(k in donors for k in segment[1:-1]):
        raise InvalidArtifact('unexpected selected side-chain donor within connected segment')
    templates = {r.attrib['name']: r for r in ET.parse(graph.topology).findall('.//Residues/Residue')}
    for rkey in segment:
        residue = graph.residues[rkey]
        expected = {a.attrib['name'] for a in templates[residue.resname].findall('Atom')}
        actual = set(graph.by_residue[rkey])
        if actual != expected:
            raise InvalidArtifact(f'incomplete/nonstandard full-residue microstate at {residue.source_dict()}: '
                                  f'missing={sorted(expected-actual)}, extra={sorted(actual-expected)}')
        extended.update(k for k in graph.by_residue[rkey].values()
                        if graph.meta[k]['element'] not in ('H', 'D'))
    c_prev = graph.amide_prev.get(left)
    n_next = graph.amide_next.get(right)
    if c_prev is None or n_next is None:
        raise InvalidArtifact('missing actual peptide neighbor for connected terminal cap')
    if graph.residues[n_next[:2]].canonical_resname == 'PRO':
        raise InvalidArtifact('connected terminal cap does not support proline N')
    extended.update(graph.key(c_prev[:2], n) for n in ('CA', 'C', 'O'))
    extended.update(graph.key(n_next[:2], n) for n in ('N', 'CA'))
    description.update({
        'complete_native_residues': [graph.residues[k].source_dict() for k in segment],
        'left_terminal_source': graph.residues[c_prev[:2]].source_dict(),
        'right_terminal_source': graph.residues[n_next[:2]].source_dict(),
        'charge_rationale': 'Existing Asp side-chain charge counted once; newly completed Gly/Gln/Ile '
                           'residues and terminal peptide caps add no formal charge.',
    })
    return extended, description


def _snapshot(output):
    directory = output / 'implementation'
    directory.mkdir()
    result = {}
    for name in DEPENDENCIES:
        source = Path(__file__).resolve().parent / name
        destination = directory / name
        before = record(source)
        shutil.copyfile(source, destination)
        saved = record(destination)
        if before['sha256'] != saved['sha256'] or record(source) != before:
            raise InvalidArtifact(f'implementation changed while copying: {source}')
        result[name] = {'historical_live_source_at_preparation': before, 'preserved_copy': saved}
    write_new(directory / 'inventory.json', result)
    return record(directory / 'inventory.json'), result


def prepare(repair_manifest, output, policy='alpha_caps'):
    started = time.monotonic()
    repair_manifest, output = Path(repair_manifest), Path(output)
    if output.exists():
        raise InvalidArtifact(f'refusing to overwrite preparation: {output}')
    parent = read_json(repair_manifest)
    if parent.get('protocol_id') != peptide.PROTOCOL or parent.get('status') != 'prepared':
        raise InvalidArtifact('only prepared generic peptide-amide v3 cores are accepted; PQQ stays unchanged')
    if policy not in PROTOCOLS:
        raise InvalidArtifact(f'unknown preparation policy: {policy}')
    verify(parent['baseline'])
    source = verify(parent['source_structure'])
    topology = verify(parent['topology_definition'])
    if parent.get('protonation_manifest') is not None:
        verify(parent['protonation_manifest'])
    for metal in ('La', 'Ca'):
        for key in ('xyz', 'input'):
            verify(parent['outputs'][metal][key])
        contents = Path(parent['outputs'][metal]['input']['path']).read_text().splitlines()
        if not contents or contents[0].strip() != RECIPE:
            raise InvalidArtifact('parent electronic recipe differs from the approved frozen recipe')
    paired(parent['outputs']['La']['xyz']['path'], parent['outputs']['Ca']['xyz']['path'],
           parent['outputs']['La']['charge'], parent['outputs']['Ca']['charge'])
    graph = peptide.SourceGraph(source, topology)
    selected, waters = _source_selection(graph, parent)
    original_atoms, _ = _materialize(graph, selected, waters)
    parent_rows = xyz(parent['outputs']['La']['xyz']['path'])
    reconstructed = [(a[0], *(float(f'{v:.10f}') for v in a[1:])) for a in original_atoms]
    if reconstructed != parent_rows[1:]:
        raise InvalidArtifact('source-graph reconstruction does not reproduce the parent v3 coordinates/order')
    if tuple(parent['selected_site']['xyz_A']) != parent_rows[0][1:]:
        raise InvalidArtifact('parent source metal coordinate changed')
    extended, extension = expand_selection(graph, selected, parent['charge_ledger'], policy)
    atoms, mapping = _materialize(graph, extended, waters)
    ligand_charge = sum(e['formal_charge'] for e in parent['charge_ledger'])
    if (parent['outputs']['La']['charge'], parent['outputs']['Ca']['charge']) != (ligand_charge+3, ligand_charge+2):
        raise InvalidArtifact('parent formal charge ledger does not close')
    output.mkdir(parents=True)
    snapshot, code = _snapshot(output)
    result = {
        'schema_version': 'alquemia.ggr_preparation.v1', 'status': 'prepared',
        'protocol_id': PROTOCOLS[policy], 'policy_id': POLICIES[policy],
        'parent_preparation': record(repair_manifest), 'baseline': parent['baseline'],
        'source_structure': record(source), 'topology_definition': record(topology),
        'implementation_snapshot': snapshot,
        'implementation': {k: v['preserved_copy'] for k, v in code.items()},
        'software_versions': {'python': platform.python_version(),
                              'numpy': importlib.metadata.version('numpy'),
                              'gemmi': importlib.metadata.version('gemmi')},
        'expected_endpoint_software': 'ORCA 6.1.1; runtime executable/version recorded by runner',
        'selected_site': parent['selected_site'], 'coordination': parent['coordination'],
        'protonation_manifest': parent.get('protonation_manifest'),
        'protonation_policy': 'reuse_exact_pinned_prepared_source_without_reprotonation',
        'assembly': {'policy': 'same_source_selected_chain_as_parent',
                     'chain': parent['selected_site']['chain']},
        'explicit_water_inventory': parent['explicit_water_inventory'],
        'charge_ledger': parent['charge_ledger'], 'ligand_formal_charge': ligand_charge,
        'extension': extension, 'atom_graph': mapping,
        'atom_graph_index_convention': 'qm_index is 1-based within ligand atoms; equal to zero-based '
                                       'full XYZ index after the metal is prepended',
        'electronic_structure': {k: parent['electronic_structure'][k]
                                 for k in ('basis_policy', 'grid', 'method', 'method_id', 'solvation')},
        'orca_recipe': RECIPE, 'decision': 'uncalibrated_protocol',
        'absolute_reference_status': 'not_registered', 'evaluation_role': 'method_development',
        'state_invariants': {'parent_source_reconstruction': 'exact_at_serialized_precision',
                             'donor_selection': 'unchanged', 'water_inventory': 'unchanged',
                             'microstate': 'unchanged', 'formal_charge': 'unchanged'},
        'outputs': {}, 'tasks': [],
    }
    for metal, charge in (('La', ligand_charge+3), ('Ca', ligand_charge+2)):
        rows = [(metal, *parent['selected_site']['xyz_A'])] + atoms
        ne = generic.require_closed_shell(metal, rows, charge)
        stem = f'{output.name}_{metal}'
        xp, inp = output / f'{stem}.xyz', output / f'{stem}.inp'
        xp.write_text('\n'.join([str(len(rows)), PROTOCOLS[policy]] +
                              [f'{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}' for a in rows]) + '\n')
        inp.write_text(f'{RECIPE}\n* xyzfile {charge} 1 {xp.name}\n')
        result['outputs'][metal] = {'xyz': record(xp), 'input': record(inp), 'charge': charge,
                                    'multiplicity': 1, 'all_electron_count_for_parity': ne}
        result['tasks'].append({'task_id': metal, 'input': record(inp), 'xyz': record(xp),
                                'output_path': str((output / f'{stem}.out').resolve())})
    result['paired_invariants'] = paired(result['outputs']['La']['xyz']['path'],
                                         result['outputs']['Ca']['xyz']['path'], ligand_charge+3, ligand_charge+2)
    serialized = xyz(result['outputs']['La']['xyz']['path'])
    shifts = [float(np.linalg.norm(np.array(serialized[a['qm_index']][1:]) -
                                  np.array(tuple(graph.atoms[_source_key(a)].pos))))
              for a in mapping['source_to_qm'] if a['kind'] == 'source']
    result['source_coordinate_max_displacement_A'] = max(shifts, default=0.)
    if result['source_coordinate_max_displacement_A'] > 1e-8:
        raise InvalidArtifact('serialized source coordinates moved')
    result['cache_key'] = cache_key({
        'protocol': result['protocol_id'], 'policy': result['policy_id'],
        'parent_sha256': record(repair_manifest)['sha256'], 'source_sha256': record(source)['sha256'],
        'topology_sha256': record(topology)['sha256'],
        'implementation_sha256': {k: v['preserved_copy']['sha256'] for k, v in code.items()},
        'assembly': result['assembly'], 'electronic_structure': result['electronic_structure'],
        'protonation_policy': result['protonation_policy'], 'water_inventory': result['explicit_water_inventory'],
        'endpoints': {m: {'xyz_sha256': o['xyz']['sha256'], 'charge': o['charge'], 'multiplicity': o['multiplicity']}
                      for m, o in result['outputs'].items()},
    })
    result['preparation_wall_seconds'] = time.monotonic()-started
    write_new(output / 'preparation_manifest.json', result)
    return result


def batch(config, output):
    config, output = Path(config), Path(output)
    if output.exists():
        raise InvalidArtifact(f'refusing to overwrite batch: {output}')
    entries = read_json(config)['preparations']
    names = [e['name'] for e in entries]
    if len(set(names)) != len(names) or any(Path(n).name != n or n in ('.', '..') for n in names):
        raise InvalidArtifact('invalid/duplicate preparation names')
    # Verify every pinned input before creating any candidate directory.
    for entry in entries:
        verify(entry['repair_manifest'])
    output.mkdir(parents=True)
    inventory = {'schema_version': 'alquemia.ggr_preparation_batch.v1', 'config': record(config),
                 'status': 'prepared', 'preparations': [], 'high_level_evaluations_executed': 0}
    for entry in entries:
        directory = output / entry['name']
        result = prepare(verify(entry['repair_manifest']), directory, entry['policy'])
        if result['paired_invariants']['atom_count'] != entry['expected_atom_count']:
            raise InvalidArtifact(f'prepared atom count differs from frozen plan: {entry["name"]}')
        inventory['preparations'].append({
            'name': entry['name'], 'manifest': record(directory / 'preparation_manifest.json'),
            'protocol_id': result['protocol_id'], 'atom_count': result['paired_invariants']['atom_count'],
            'source_coordinate_max_displacement_A': result['source_coordinate_max_displacement_A'],
            'ligand_formal_charge': result['ligand_formal_charge'],
            'new_endpoint_count': len(result['tasks']),
        })
    write_new(output / 'inventory.json', inventory)
    return inventory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('prepare')
    p.add_argument('--repair-manifest', required=True, type=Path)
    p.add_argument('--policy', choices=tuple(PROTOCOLS), required=True)
    p.add_argument('--output', required=True, type=Path)
    b = commands.add_parser('batch')
    b.add_argument('--config', required=True, type=Path)
    b.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.command == 'prepare':
        result = prepare(args.repair_manifest, args.output, args.policy)
        print(f'{result["protocol_id"]}: {result["paired_invariants"]["atom_count"]} atoms; prepared only')
    else:
        result = batch(args.config, args.output)
        print(f'{len(result["preparations"])} pairs prepared; no electronic-structure calculations run')


if __name__ == '__main__':
    main()
