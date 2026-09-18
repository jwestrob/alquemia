"""Carve-free, source-graph charge groups; no energy or coordinate modification."""
from collections import defaultdict
import ast
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact
from mace_charge_group_prepare import PROTEIN, residue

POLICY = 'typed31_complete_residue_amide_multimetal_charge_groups_v1'


def typed_policy():
    """Load exact pure donor definitions without importing unrelated gemmi code.

    The complete coordination_policy.py source is pinned in the manifest.
    No donor table/function is rewritten; only selected AST definitions from
    that local source are compiled into an isolated namespace.
    """
    wanted = {'COORDINATION_CUTOFF_A', 'STANDARD_AA', 'WATER_NAMES',
              'SIDECHAIN_O_DONORS', 'SIDECHAIN_N_DONORS',
              'COFACTOR_O_DONORS', 'COFACTOR_N_DONORS', 'donor_type'}
    tree = ast.parse(Path(__file__).with_name('coordination_policy.py').read_text())
    nodes = []; found = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
        elif isinstance(node, ast.FunctionDef):
            name = node.name
        else:
            continue
        if name in wanted:
            nodes.append(node); found.append(name)
    if set(found) != wanted or len(found) != len(wanted):
        raise InvalidArtifact('typed donor source definitions changed')
    namespace = {}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), 'pinned_coordination_definitions', 'exec'), namespace)
    return namespace['donor_type'], namespace['COORDINATION_CUTOFF_A']


donor_type, COORDINATION_CUTOFF_A = typed_policy()


def build(prep, ledger, bonds, extra_residues=()):
    atoms = sorted(prep['physical_atoms'], key=lambda a: a['id'])
    byid = {a['id']: a for a in atoms}
    if len(byid) != len(atoms) or COORDINATION_CUTOFF_A != 3.1:
        raise InvalidArtifact('duplicate source identity or changed contact policy')
    metals = [a for a in atoms if a['kind'] in ('selected_metal', 'background_metal')]
    selected = [a for a in metals if a['kind'] == 'selected_metal']
    if len(selected) != 1 or any(a['element'] != 'Ca' for a in metals if a['kind'] == 'background_metal'):
        raise InvalidArtifact('one selected endpoint and explicit background Ca required')
    adjacency = defaultdict(set)
    for a, b in bonds:
        if a not in byid or b not in byid:
            raise InvalidArtifact('source bond atom absent')
        adjacency[a].add(b); adjacency[b].add(a)
    units = {}; base_charges = {}; unit_atoms = defaultdict(list)
    for a in atoms:
        pid, kind = a['id'], a['kind']
        if kind in PROTEIN:
            key = 'protein:' + residue(pid)
            q = ledger[residue(pid)]['formal_charge_e']
        elif kind in ('selected_metal', 'background_metal'):
            key, q = 'metal:' + pid, 2
        elif kind == 'retained_site_water':
            key, q = 'water:' + residue(pid), 0
        elif kind == 'frozen_PQQ':
            key, q = 'cofactor:PQQ', prep['cofactor_charge_e']
            if q != -3:
                raise InvalidArtifact('only the existing complete PQQ -3 state is supported')
        else:
            raise InvalidArtifact('unsupported source group kind: ' + kind)
        units[pid] = key; base_charges[key] = q; unit_atoms[key].append(pid)
    for key, members in unit_atoms.items():
        if key.startswith('water:') and sorted(byid[pid]['element'] for pid in members) != ['H', 'H', 'O']:
            raise InvalidArtifact('incomplete neutral source water')
    parents = {key: key for key in base_charges}

    def root(k):
        while k != parents[k]:
            k = parents[k]
        return k

    joins = []

    def join(a, b, reason):
        ra, rb = root(a), root(b)
        parents[max(ra, rb)] = min(ra, rb)
        joins.append({'unit_a': a, 'unit_b': b, 'reason': reason})

    for a, b in bonds:
        if a.endswith('/SG') and b.endswith('/SG'):
            join(units[a], units[b], 'source_disulfide')
    contacts = []
    for metal in metals:
        for a in atoms:
            if a['kind'] in PROTEIN:
                rn, name = a['resname'], a['id'].rsplit('/', 1)[1]
                dtype = donor_type(rn, name, a['element'])
                if dtype == 'sidechain_N' and any(byid[p]['element'] == 'H' for p in adjacency[a['id']]):
                    continue
            elif a['kind'] == 'retained_site_water':
                dtype = donor_type('HOH', a['id'].rsplit('/', 1)[1], a['element'])
            elif a['kind'] == 'frozen_PQQ':
                dtype = donor_type('PQQ', a['source_record']['name'], a['element'])
            else:
                continue
            distance = float(np.linalg.norm(np.asarray(a['xyz_A']) - metal['xyz_A']))
            if dtype is None or distance > 3.1:
                continue
            contacts.append({'metal_id': metal['id'], 'donor_id': a['id'],
                             'donor_type': dtype, 'distance_A': distance})
            join(units[metal['id']], units[a['id']], 'typed31_contact')
            if dtype == 'backbone_O':
                carbonyl = [p for p in adjacency[a['id']] if p.endswith('/C') and byid[p]['element'] == 'C']
                if len(carbonyl) != 1:
                    raise InvalidArtifact('coordinating carbonyl lacks a unique source-bonded C')
                for p in adjacency[carbonyl[0]]:
                    if p.endswith('/N') and byid[p]['element'] == 'N':
                        join(units[a['id']], units[p], 'actual_carbonyl_amide_CN')
    for metal in metals:
        if not any(c['metal_id'] == metal['id'] for c in contacts):
            raise InvalidArtifact('metal has no supported typed donor in the declared contact radius')
    for r in sorted(extra_residues):
        if r not in ledger:
            raise InvalidArtifact('alternate group residue absent')
        join(units[selected[0]['id']], 'protein:' + r, 'declared_connected_residue_union')
    names = sorted({root(k) for k in parents}); index = {k: i for i, k in enumerate(names)}
    ids = [index[root(units[a['id']])] for a in atoms]
    charges = [0] * len(names)
    for k, q in base_charges.items():
        charges[index[root(k)]] += q
    site = index[root(units[selected[0]['id']])]
    endpoint_charges = {'Ca': charges, 'La': list(charges)}
    endpoint_charges['La'][site] += 1
    for metal, qs in endpoint_charges.items():
        if sum(qs) != prep['endpoints'][metal]['charge']:
            raise InvalidArtifact('formal multimetal charge closure failed')
    return {'policy_id': POLICY, 'atom_ids': [a['id'] for a in atoms],
            'all_Ca_atoms': [('Ca' if a['kind'] in ('selected_metal', 'background_metal') else a['element'], *a['xyz_A']) for a in atoms],
            'group_names': names, 'atom_group_indices': ids, 'endpoint_group_charges_e': endpoint_charges,
            'selected_metal_id': selected[0]['id'], 'selected_metal_index': next(i for i, a in enumerate(atoms) if a['id'] == selected[0]['id']),
            'all_metal_indices': [i for i, a in enumerate(atoms) if a in metals],
            'selected_group_index': site, 'contacts': contacts, 'joins': joins,
            'group_units': [[k for k in sorted(parents) if root(k) == name] for name in names]}
