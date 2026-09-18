"""Prepare complete chemical groups on existing intact, source-mapped proteins."""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import time
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz

POLICY = 'intact_POLAR_formal_chemical_group_constraint_v1'
CASES = {'GGR_1GLG', 'GGR_2FW0', 'GGR_2FVY', 'ALPHA_1F6S', 'ALPHA_6IP9',
         'PQQ_1H4I', 'PQQ_4MAE'}
PROTEIN = {'protein_source', 'terminal_completion'}


def residue(pid):
    fields = pid.split('/')
    if len(fields) != 4:
        raise InvalidArtifact('protein source identity must include chain/residue/insertion/atom')
    return '/'.join(fields[:3])


def protein_ledger(prep):
    """Replay existing source bonds and templates; evaluate no force-field energy."""
    import openmm as mm
    from openmm import app, unit
    physical = prep['physical_atoms']
    rows = [a for a in physical if a['kind'] in PROTEIN]
    top = app.Topology()
    chains, residues, atoms = {}, {}, {}
    for a in rows:
        chain, rid, icode, name = a['id'].split('/')
        key = residue(a['id'])
        if chain not in chains:
            chains[chain] = top.addChain(chain)
        if key not in residues:
            residues[key] = top.addResidue(a['resname'], chains[chain], rid, icode)
        if residues[key].name != a['resname'] or a['id'] in atoms:
            raise InvalidArtifact('inconsistent or duplicate source residue identity')
        atoms[a['id']] = top.addAtom(name, app.element.get_by_symbol(a['element']), residues[key])
    seen = set()
    for b in prep['preparation_details']['bonds']:
        pair = tuple(sorted((b['atom_a_id'], b['atom_b_id'])))
        if pair in seen or any(pid not in atoms for pid in pair):
            raise InvalidArtifact('duplicate or missing source protein bond')
        seen.add(pair)
        top.addBond(*(atoms[pid] for pid in pair))
    ff = app.ForceField(str(verify(prep['forcefield'])))
    system = ff.createSystem(top, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    partial = defaultdict(list)
    for a in rows:
        q = nb.getParticleParameters(atoms[a['id']].index)[0].value_in_unit(unit.elementary_charge)
        partial[residue(a['id'])].append(float(q))
    ledger = {}
    for key, qs in partial.items():
        q = math.fsum(qs)
        if abs(q - round(q)) > 1e-6:
            raise InvalidArtifact('nonintegral template residue charge: ' + key)
        ledger[key] = {'formal_charge_e': round(q), 'template_charge_sum_e': q}
    if sum(x['formal_charge_e'] for x in ledger.values()) != prep['protein_charge_e']:
        raise InvalidArtifact('protein formal charge ledger does not close')
    return ledger, sorted(seen)


def selected_residues(prep, core_pin, core_format):
    core = read_json(verify(core_pin))
    physical = {a['id']: a for a in prep['physical_atoms']}
    wanted = set()
    if core_format == 'physical_mapping':
        for row in core['mapping']:
            if row['kind'] != 'source':
                continue
            pid = row['physical_id']
            if pid not in physical:
                raise InvalidArtifact('core source atom is absent from intact protein')
            a = physical[pid]
            if a['element'] != 'H' and a['kind'] in PROTEIN:
                wanted.add(residue(pid))
    elif core_format == 'canonical_PQQ':
        for fragment in core['qm_fragments']:
            if fragment['kind'] == 'fixed_core_pqq':
                continue
            if fragment['kind'] not in ('fixed_core_protein_sidechain', 'fixed_core_cationic_sidechain'):
                raise InvalidArtifact('unsupported canonical fragment')
            for atom in fragment['atom_records']:
                if atom['origin'] != 'source_heavy_atom':
                    continue
                candidates = [a for a in physical.values() if a['kind'] in PROTEIN
                              and a['id'].rsplit('/', 1)[1] == atom['name']
                              and a['element'] == atom['element']
                              and np.allclose(a['xyz_A'], atom['xyz_A'], atol=1e-6, rtol=0)]
                if len(candidates) != 1:
                    raise InvalidArtifact('ambiguous canonical heavy-atom mapping')
                wanted.add(residue(candidates[0]['id']))
    else:
        raise InvalidArtifact('unknown core mapping format')
    if not wanted:
        raise InvalidArtifact('empty coordination block')
    return wanted


def groups(prep, ledger, bonds, selected):
    physical = prep['physical_atoms']
    ids = [a['id'] for a in physical]
    if len(set(ids)) != len(ids) or not selected <= ledger.keys():
        raise InvalidArtifact('inconsistent atom/group mapping')
    parents = {key: key for key in ledger}
    def root(key):
        while key != parents[key]:
            key = parents[key]
        return key
    def join(a, b):
        ra, rb = root(a), root(b)
        parents[max(ra, rb)] = min(ra, rb)
    ordered = sorted(selected)
    for key in ordered[1:]:
        join(ordered[0], key)
    for a, b in bonds:
        if a.endswith('/SG') and b.endswith('/SG'):
            join(residue(a), residue(b))
    site = root(ordered[0])
    names = sorted({root(key) for key in ledger})
    index = {key: i for i, key in enumerate(names)}
    members = [[] for _ in names]
    charges = [0 for _ in names]
    for key, item in ledger.items():
        charges[index[root(key)]] += item['formal_charge_e']
    assignment = []
    for i, a in enumerate(physical):
        if a['kind'] in PROTEIN:
            key = root(residue(a['id']))
        elif a['kind'] in ('frozen_PQQ', 'retained_site_water', 'selected_metal'):
            key = site
        else:
            raise InvalidArtifact('unsupported atom kind in charge groups: ' + a['kind'])
        assignment.append(index[key])
        members[index[key]].append(i)
    charges[index[site]] += prep['cofactor_charge_e']
    endpoints = {}
    for metal, formal in (('Ca', 2), ('La', 3)):
        q = list(charges)
        q[index[site]] += formal
        if sum(q) != prep['endpoints'][metal]['charge']:
            raise InvalidArtifact('group endpoint charge mismatch')
        endpoints[metal] = q
    if any(not x for x in members) or sorted(sum(members, [])) != list(range(len(physical))):
        raise InvalidArtifact('groups do not partition the full physical system')
    return {'atom_ids': ids, 'atom_group_indices': assignment,
            'group_members': members, 'group_names': names, 'site_group_index': index[site],
            'group_residues': [[key for key in sorted(ledger) if root(key) == name] for name in names],
            'selected_residues': sorted(selected), 'endpoint_group_charges_e': endpoints}


def prepare(config, output):
    start, cpu = time.monotonic(), time.process_time()
    cfg = read_json(config)
    if set(cfg['cases']) != CASES:
        raise InvalidArtifact('declared seven-case inventory differs')
    for key in ('agreement', 'source_inventory'):
        verify(cfg[key])
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    rows = {}
    for name, item in cfg['cases'].items():
        prep = read_json(verify(item['physical']))
        if prep['case_id'] != name:
            raise InvalidArtifact('physical source identity differs')
        physical = prep['physical_atoms']
        for metal in ('Ca', 'La'):
            expected = [(metal if a['kind'] == 'selected_metal' else a['element'], *a['xyz_A']) for a in physical]
            if xyz(verify(prep['endpoints'][metal]['xyz'])) != expected:
                raise InvalidArtifact('paired physical source differs')
        ledger, bonds = protein_ledger(prep)
        primary = selected_residues(prep, item['primary_core'], item['core_format'])
        row = {'sources': item, 'protein_ledger': ledger,
               'primary': groups(prep, ledger, bonds, primary)}
        if 'alternate_core' in item:
            alt = selected_residues(prep, item['alternate_core'], item['core_format'])
            row['connected'] = groups(prep, ledger, bonds, primary | alt)
        path = out / (name + '.json')
        write_new(path, row)
        rows[name] = record(path)
    result = {'policy_id': POLICY, 'status': 'prepared', 'config': record(config), 'cases': rows,
              'implementation': record(__file__), 'new_model_calls': 0, 'new_forcefield_energies': 0,
              'wall_seconds': time.monotonic() - start, 'CPU_seconds': time.process_time() - cpu}
    write_new(out / 'preparation.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.config, args.output), indent=2))
