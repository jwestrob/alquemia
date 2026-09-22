"""Read-only protein parameter/term inventory; never creates an energy Context."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import resource
import time

import numpy as np
import openmm as mm
from openmm import app, unit

from affordable_common import InvalidArtifact, paired, read_json, record, verify, write_new, xyz

CASES = ('1H4I', '4MAE', 'PQQSEQ_83440678cbbd658047c9')
STANDARD = set('ALA ARG ASN ASP CYS GLN GLU GLY HID HIE HIP HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
PROTOCOL = 'source_protein_ff19SB_term_inventory_no_energy_v1'


def atom_id(a):
    return f'{a.residue.chain.id}/{a.residue.id}/{a.residue.insertionCode.strip()}/{a.name}'


def metadata_id(m):
    return f"{m['chain']}/{m['resnum']}/{m['insertion_code']}/{m['atom']}"


def partition(indices, local):
    n = sum(int(i in local) for i in indices)
    return 'local' if n == len(indices) else 'exterior' if n == 0 else 'mixed'


def parent_model(source, ffpath, archived_repair):
    pdb = app.PDBFile(str(source))
    residues = list(pdb.topology.residues())
    unknown = [r for r in residues if r.name not in STANDARD | {'LA', 'CA', 'PQQ'}]
    if unknown:
        raise InvalidArtifact('unsupported source residues: '+str([(r.chain.id, r.id, r.name) for r in unknown]))
    if any(r.chain.id != 'A' for r in residues if r.name in STANDARD):
        raise InvalidArtifact('additional actual protein chain outside fixed chain-A scope')
    excluded = [{'chain': r.chain.id, 'resid': r.id, 'name': r.name,
                 'atoms': len(list(r.atoms()))} for r in residues if r.name not in STANDARD]
    if sorted(r['name'] for r in excluded) not in (['LA', 'PQQ'], ['CA', 'PQQ']):
        raise InvalidArtifact('expected exactly source metal and PQQ; no silent heterogen omission')
    model = app.Modeller(pdb.topology, pdb.positions)
    model.delete([a for a in model.topology.atoms() if a.residue.name not in STANDARD])
    pos = np.asarray(model.positions.value_in_unit(unit.angstrom))
    atoms = list(model.topology.atoms()); original_count = len(atoms)
    ids = [atom_id(a) for a in atoms]
    if len(set(ids)) != len(ids):
        raise InvalidArtifact('nonunique source protein identity')
    ff = app.ForceField(str(ffpath)); exact_error = None
    try:
        system = ff.createSystem(model.topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    except ValueError as exc:
        exact_error = str(exc); system = None
    additions = []
    if archived_repair:
        repair = read_json(verify(archived_repair))
        if repair['source'] != record(source) or repair['forcefield'] != record(ffpath):
            raise InvalidArtifact('old terminal repair has different source/force field')
        additions = repair['preparation_details']['terminal_additions']
        if not exact_error or len(additions) != 1:
            raise InvalidArtifact('unexpected terminal-repair scope')
        addition = additions[0]; byid = {atom_id(a): a for a in atoms}
        if addition['id'] in byid or addition['recipe'] != 'reflect_carbonyl_about_C_CA_v1':
            raise InvalidArtifact('terminal completion mismatch')
        carbon = byid[addition['bond_to']]
        if addition['id'] != atom_id(carbon).rsplit('/', 1)[0]+'/OXT':
            raise InvalidArtifact('terminal repair changes residue')
        neighbours = [b if a == carbon else a for a, b in model.topology.bonds() if carbon in (a, b)]
        if {a.name for a in neighbours} != {'O', 'CA'} or any(a.residue != carbon.residue for a in neighbours):
            raise InvalidArtifact('terminal source connectivity differs')
        a = model.topology.addAtom('OXT', app.element.oxygen, carbon.residue)
        model.topology.addBond(carbon, a)
        pos = np.concatenate((pos, np.asarray(addition['xyz_A'])[None]))
        system = ff.createSystem(model.topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    if system is None:
        raise InvalidArtifact('unparameterized exact parent: '+str(exact_error))
    templates = ff.getMatchingTemplates(model.topology)
    template_rows = [{'chain': r.chain.id, 'resid': r.id, 'icode': r.insertionCode.strip(),
                      'source_name': r.name, 'matched_template': t.name,
                      'actual_H_names': [a.name for a in r.atoms() if a.element.symbol == 'H']}
                     for r, t in zip(model.topology.residues(), templates)]
    return model, pos, system, ff, {'exact_source_template_error': exact_error,
        'original_protein_atom_count': original_count, 'archived_terminal_additions': additions,
        'excluded_nonprotein_inventory': excluded, 'residue_proton_templates': template_rows,
        'source_H_moved': False, 'global_archived_H_normalization_reused': False}


def term_inventory(system, local, positions, atoms):
    force_rows = []; details = []; hbonds = []
    for fi, f in enumerate(system.getForces()):
        name = type(f).__name__; supports = []
        if isinstance(f, mm.HarmonicBondForce):
            for i in range(f.getNumBonds()):
                a, b, length, _ = f.getBondParameters(i); supports.append([int(a), int(b)])
                if 'H' in (atoms[a].element.symbol, atoms[b].element.symbol):
                    actual = float(np.linalg.norm(positions[a]-positions[b]))
                    equilibrium = length.value_in_unit(unit.angstrom)
                    hbonds.append({'indices': [a, b], 'partition': partition([a, b], local),
                                   'source_length_A': actual, 'ff_length_A': equilibrium,
                                   'difference_A': actual-equilibrium})
        elif isinstance(f, mm.HarmonicAngleForce):
            supports = [list(map(int, f.getAngleParameters(i)[:3])) for i in range(f.getNumAngles())]
        elif isinstance(f, mm.PeriodicTorsionForce):
            supports = [list(map(int, f.getTorsionParameters(i)[:4])) for i in range(f.getNumTorsions())]
        elif isinstance(f, mm.CMAPTorsionForce):
            supports = [list(map(int, f.getTorsionParameters(i)[1:])) for i in range(f.getNumTorsions())]
        elif isinstance(f, mm.NonbondedForce):
            supports = [list(map(int, f.getExceptionParameters(i)[:2])) for i in range(f.getNumExceptions())]
        elif isinstance(f, mm.CMMotionRemover):
            force_rows.append({'index': fi, 'type': name, 'energy_term': False}); continue
        else:
            raise InvalidArtifact('unhandled parameterized force: '+name)
        counts = Counter(partition(s, local) for s in supports)
        row = {'index': fi, 'type': name, 'support_counts': dict(counts), 'support_count': len(supports)}
        if isinstance(f, mm.NonbondedForce):
            n = system.getNumParticles(); c = len(local); e = n-c
            all_pairs = {'local': c*(c-1)//2, 'exterior': e*(e-1)//2, 'mixed': c*e}
            row.update({'support_kind': 'explicit_exceptions_including_zero_exclusions',
                        'ordinary_pair_counts': {k: v-counts[k] for k, v in all_pairs.items()},
                        'nonbonded_method': int(f.getNonbondedMethod()),
                        'charge_e': float(sum(f.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(n))),
                        'local_subset_charge_e': float(sum(f.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in local))})
        details.append({'index': fi, 'type': name, 'supports': supports,
                        'partitions': [partition(s, local) for s in supports]})
        force_rows.append(row)
    return force_rows, details, hbonds


def local_template_probe(topology, local, caps, ff):
    """Check real capped protein pieces; never invent charges or parameters."""
    top = app.Topology(); chain = top.addChain('A'); amap = {}
    for r in topology.residues():
        retained_atoms = [a for a in r.atoms() if a.index in local]
        if not retained_atoms: continue
        residue = top.addResidue(r.name, chain, r.id, r.insertionCode)
        for a in retained_atoms:
            amap[a.index] = top.addAtom(a.name, a.element, residue)
        for cap in caps:
            if cap['retained_parent_index'] not in {a.index for a in retained_atoms}: continue
            retained = amap[cap['retained_parent_index']]
            h = top.addAtom('LINK'+str(cap['context_index']), app.element.hydrogen, residue)
            top.addBond(retained, h)
    for a, b in topology.bonds():
        if a.index in local and b.index in local: top.addBond(amap[a.index], amap[b.index])
    return [{'chain': r.chain.id, 'resid': r.id, 'resname': r.name, 'atom_count': len(list(r.atoms()))}
            for r in ff.getUnmatchedResidues(top)]


def inspect_case(cid, source_manifest, config, output):
    tasks = [next(t for t in source_manifest['tasks'] if t['case_id'] == cid and t['metal'] == z) for z in ('Ca', 'La')]
    paired(verify(tasks[1]['xyz']), verify(tasks[0]['xyz']), tasks[1]['charge'], tasks[0]['charge'])
    geometries = [read_json(verify(t['mapping'])) for t in tasks]
    if geometries[0] != geometries[1] or tasks[0]['source_preparation'] != tasks[1]['source_preparation']:
        raise InvalidArtifact('paired physical maps/preparations differ')
    prep = read_json(verify(tasks[0]['source_preparation'])); source = verify(prep['source'])
    case = next(c for c in source_manifest['cases'] if c['case_id'] == cid)
    carve = read_json(verify(case['source']['parent'])); proton = carve['protonation_manifest']
    if read_json(verify(proton))['output'] != prep['source']:
        raise InvalidArtifact('source protonation record mismatch')
    model, positions, system, ff, state = parent_model(source, verify(config['forcefield']), config['terminal_repairs'].get(cid))
    atoms = list(model.topology.atoms()); byid = {atom_id(a): a.index for a in atoms}
    g = geometries[0]['context']; links = g['core_links']; local = set(); physical = []
    pmap = {}; max_delta = 0.
    for i, meta in enumerate(g['source_atom_metadata']):
        if meta is None or meta['resname'] not in STANDARD: continue
        aid = metadata_id(meta)
        if aid not in byid: raise InvalidArtifact('missing exact protein source atom '+aid)
        j = byid[aid]
        if meta['element'] != atoms[j].element.symbol: raise InvalidArtifact('mapped element differs')
        delta = float(np.max(np.abs(positions[j]-g['positions_A'][i])))
        if delta > 1e-12: raise InvalidArtifact('physical map moved source protein coordinates: '+aid)
        max_delta = max(max_delta, delta); pmap[i] = j
        physical.append({'physical_index': i, 'parent_index': j, 'source_id': aid})
    context_links = []; caps = []; other = []
    for link in links:
        ci, kind = link[:2]
        if kind == 'source' and link[2] in pmap:
            j = pmap[link[2]]; local.add(j)
            context_links.append({'context_index': ci, 'parent_index': j, 'source_id': atom_id(atoms[j])})
        elif kind == 'cap':
            if link[2] not in pmap or link[3] not in pmap:
                raise InvalidArtifact('unsupported nonprotein covalent boundary')
            caps.append({'context_index': ci, 'retained_parent_index': pmap[link[2]],
                         'omitted_parent_index': pmap[link[3]], 'length_A': link[4], 'q0_xyz_A': link[5]})
        else:
            other.append(link)
    source_bonds = {tuple(sorted((a.index, b.index))) for a, b in model.topology.bonds()}
    for cap in caps:
        a, b = cap['retained_parent_index'], cap['omitted_parent_index']
        if a not in local or b in local or tuple(sorted((a, b))) not in source_bonds:
            raise InvalidArtifact('cap does not replace exactly one real local/exterior bond')
    crossing = [list(pair) for pair in sorted(source_bonds) if partition(pair, local) == 'mixed']
    cap_pairs = sorted(tuple(sorted((c['retained_parent_index'], c['omitted_parent_index']))) for c in caps)
    if [tuple(p) for p in crossing] != cap_pairs:
        raise InvalidArtifact('real cut bonds differ from cap coverage')
    forces, terms, hbonds = term_inventory(system, local, positions, atoms)
    unmatched = local_template_probe(model.topology, local, caps, ff)
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    atom_rows = [{'parent_index': a.index, 'source_id': atom_id(a), 'resname': a.residue.name,
                  'element': a.element.symbol, 'xyz_A': positions[a.index].tolist(),
                  'in_actual_MACE_context': a.index in local,
                  'source_or_archived_completion': 'source' if a.index < state['original_protein_atom_count'] else 'archived_terminal_completion',
                  'charge_e': nb.getParticleParameters(a.index)[0].value_in_unit(unit.elementary_charge)} for a in atoms]
    out = output/cid; out.mkdir(parents=True)
    (out/'parent_system.xml').write_text(mm.XmlSerializer.serialize(system))
    write_new(out/'parent_atoms.json', atom_rows)
    write_new(out/'term_supports.json', terms)
    write_new(out/'source_proton_inventory.json', state)
    write_new(out/'context_mapping.json', {'physical_nodes': physical, 'source_context_atoms': context_links,
              'caps': caps, 'other_source_links': other, 'actual_cut_bonds': crossing,
              'unmapped_fixed_context_indices': sorted(set(range(len(g['core_positions_A'])))-{l[0] for l in links}),
              'local_parent_indices': sorted(local), 'max_source_coordinate_difference_A': max_delta})
    write_new(out/'source_H_bond_geometry.json', hbonds)
    hdiff = np.asarray([a['difference_A'] for a in hbonds])
    return {'case_id': cid, 'status': 'parameter_inventory_complete_score_coupling_unresolved',
            'source': prep['source'], 'source_preparation': tasks[0]['source_preparation'],
            'source_carve': case['source']['parent'], 'source_protonation': proton,
            'endpoints': {t['metal']: {k: t[k] for k in ('xyz', 'mapping', 'charge', 'multiplicity')} for t in tasks},
            'terminal_repair': config['terminal_repairs'].get(cid), 'source_atom_count': state['original_protein_atom_count'],
            'parent_atom_count': len(atoms), 'parent_H_count': sum(a.element.symbol == 'H' for a in atoms),
            'local_protein_atom_count': len(local), 'exterior_protein_atom_count': len(atoms)-len(local),
            'context_atom_count': len(g['core_positions_A']), 'cap_count': len(caps),
            'parameterized_force_supports': forces, 'local_capped_unmatched_residues': unmatched,
            'source_H_bond_geometry': {'count': len(hbonds), 'max_abs_length_difference_A': float(np.max(np.abs(hdiff))),
                'median_abs_length_difference_A': float(np.median(np.abs(hdiff))), 'abs_difference_over_0p1_A': int(sum(np.abs(hdiff) > .1))},
            'artifacts': {p.name: record(p) for p in sorted(out.iterdir())},
            'constraints': system.getNumConstraints(), 'energy_evaluations': 0, 'force_evaluations': 0,
            'missing_parameters': ['matched capped-local protein reference', 'protein_environment_to_PQQ_or_metal'],
            'absolute_or_relaxation_score_available': False}


def inventory(config_path, output):
    start = time.monotonic(); config = read_json(config_path)
    if tuple(config['cases']) != CASES: raise InvalidArtifact('fixed three-source scope changed')
    manifest = read_json(verify(config['source_manifest'])); verify(config['agreement'])
    output = Path(output).resolve(); output.mkdir(parents=True, exist_ok=False)
    rows = []
    for cid in CASES:
        try: rows.append(inspect_case(cid, manifest, config, output))
        except Exception as exc: rows.append({'case_id': cid, 'status': 'inventory_failed', 'error': str(exc), 'exception': type(exc).__name__})
    result = {'protocol_id': PROTOCOL, 'config': record(config_path), 'implementation': record(__file__),
              'source_manifest': config['source_manifest'], 'forcefield': config['forcefield'],
              'openmm_version': mm.__version__, 'rows': rows, 'case_denominator': 3,
              'energy_evaluations': 0, 'force_evaluations': 0, 'model_or_default_changed': False,
              'wall_seconds': time.monotonic()-start,
              'process_CPU_seconds': resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime}
    write_new(output/'INVENTORY.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', required=True); p.add_argument('--output', required=True)
    a = p.parse_args(); result = inventory(a.config, a.output)
    print(json.dumps({'output': str(Path(a.output).resolve()/'INVENTORY.json'),
                      'rows': [{k: r[k] for k in ('case_id', 'status')} for r in result['rows']],
                      'energy_evaluations': 0, 'force_evaluations': 0}, indent=2))
    return 0 if all(r['status'].startswith('parameter_inventory_complete') for r in result['rows']) else 1


if __name__ == '__main__':
    raise SystemExit(main())
