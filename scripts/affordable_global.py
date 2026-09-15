"""Audit saved boundary terms and export a cap-free, whole-chain development pair.

This prepares a proposed new Hamiltonian; it neither runs nor calibrates it.
Only the already-consumed 1H4I boundary-test physical system is supported.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
import time

import gemmi
import numpy as np

from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new
from affordable_environment import COULOMB_KCAL_A, physical_boundary_key, validate_state
from affordable_state import source_protein
from pqq_microstates import CCD_SCHEMA

PROTOCOL = 'whole_chain_native_gfn2_alpb_water_vertical_proposed_v1'
Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'S': 16, 'Ca': 20, 'La': 57}


def load_primary(result_path):
    result = read_json(result_path)
    states, components, provenance = {}, {}, {}
    for row in result['numerical_checks']:
        if not row['label'].endswith('/primary'):
            continue
        key = row['label'].split('/')[0]
        m = read_json(verify(row['result']['source_manifest']))
        state_path = verify(m['state'])
        s = read_json(state_path)
        validate_state(s)
        states[key] = s
        components[key] = row['result']['components']
        provenance[key] = record(state_path)
    if set(states) != {f'1h4i_qm{n}_{metal}' for n in (33, 36) for metal in ('La', 'Ca')}:
        raise InvalidArtifact('requires four complete, pinned boundary-development states')
    return states, components, provenance


def residue_id(atom_id):
    return '/'.join(atom_id.split('/')[:3])


def direct_decomposition(la, ca):
    if la['environment_atoms'] != ca['environment_atoms']:
        raise InvalidArtifact('permanent environment changed between metals')
    if [(a['id'], a['xyz_A']) for a in la['core_atoms']] != [(a['id'], a['xyz_A']) for a in ca['core_atoms']]:
        raise InvalidArtifact('core coordinates/order changed between metals')
    e = la['environment_atoms']
    pos = np.array([a['xyz_A'] for a in e])
    charge = np.array([a['charge_e'] for a in e])
    contributions = np.zeros(len(e))
    core_terms = []
    for a, b in zip(la['core_atoms'], ca['core_atoms']):
        dq = b['charge_e'] - a['charge_e']
        term = COULOMB_KCAL_A * dq * charge / np.linalg.norm(pos - a['xyz_A'], axis=1)
        contributions += term
        core_terms.append({'id': a['id'], 'Ca_minus_La_charge_e': dq,
                           'direct_score_kcal_mol': float(term.sum())})
    residues = defaultdict(float)
    atom_terms = []
    for a, value in zip(e, contributions):
        residues[residue_id(a['id'])] += float(value)
        atom_terms.append({'id': a['id'], 'permanent_charge_e': a['charge_e'],
                           'direct_score_kcal_mol': float(value)})
    return {'direct_score_kcal_mol': float(contributions.sum()),
            'per_environment_residue_kcal_mol': dict(sorted(residues.items())),
            'per_environment_atom': atom_terms, 'per_core_atom': core_terms}


def audit(result_path, output):
    start = time.monotonic()
    states, components, provenance = load_primary(result_path)
    rows = {}
    for n in (33, 36):
        key = f'1h4i_qm{n}'
        d = direct_decomposition(states[key+'_La'], states[key+'_Ca'])
        contrasts = {k: components[key+'_Ca'][k] - components[key+'_La'][k]
                     for k in components[key+'_La']}
        if abs(d['direct_score_kcal_mol'] - contrasts['direct_core_environment_kcal_mol']) > 1e-9:
            raise InvalidArtifact('saved direct energy does not reproduce')
        rows[key] = dict(d, saved_component_contrasts_kcal_mol=contrasts)
    report = {'protocol': 'saved_boundary_component_audit_v1', 'implementation': record(__file__),
              'source_result': record(result_path), 'source_states': provenance, 'partitions': rows,
              'new_solver_or_quantum_evaluations': 0,
              'interpretation': 'Exact decomposition of the saved descriptor, not causal isolation of a physical effect.',
              'wall_seconds': time.monotonic()-start}
    write_new(output, report)
    return report


def prepare(result_path, carve_path, output, orca):
    start = time.monotonic()
    states, _, provenance = load_primary(result_path)
    if len({physical_boundary_key(s['physical_atoms']) for s in states.values()}) != 1:
        raise InvalidArtifact('the old partitions do not describe one physical atom set')
    s = states['1h4i_qm33_La']
    carve = read_json(carve_path)
    if carve['source_structure'] != s['source'] or carve['pqq']['microstate_id'] != 'pqq_ox_3minus_v1':
        raise InvalidArtifact('source or PQQ state mismatch')
    source = verify(s['source'])
    atoms, pos, ffcharges, _, excluded = source_protein(source)
    physical = [{k: a[k] for k in ('id', 'element', 'xyz_A')} for a in s['physical_atoms']]
    index = {a['id']: i for i, a in enumerate(physical)}
    if len(index) != len(physical) or any(a['id'].startswith('cap/') for a in physical):
        raise InvalidArtifact('duplicate physical atom or synthetic cap')
    # OpenMM source connectivity is already validated against complete ff19SB
    # templates by source_protein; do not manufacture peptide links by numbering.
    bonds = {tuple(sorted((a.index, b.index))) for a, b in atoms[0].residue.chain.topology.bonds()}
    max_source_error = float(np.max(np.abs(pos - np.array([a['xyz_A'] for a in physical[:len(atoms)]]))))
    if max_source_error > 1e-9:
        raise InvalidArtifact('whole-chain export would alter source coordinates')
    pqq = carve['pqq']['residues'][0]
    if pqq['atom_naming_schema'] != CCD_SCHEMA.schema_id:
        raise InvalidArtifact('unsupported PQQ graph schema')
    resid = pqq['residue']
    structure = gemmi.read_structure(str(source))
    candidates = [r for c in structure[0] if c.name == resid['chain'] for r in c
                  if r.name == 'PQQ' and r.seqid.num == resid['resnum'] and r.seqid.icode.strip() == resid['icode']]
    if len(candidates) != 1:
        raise InvalidArtifact('ambiguous source PQQ')
    named = {a.name: (a.element.name, [a.pos.x, a.pos.y, a.pos.z]) for a in candidates[0] if a.element.name != 'H'}
    named.update({h['name']: ('H', h['xyz_A']) for h in pqq['hydrogens']})
    pqq_map = {}
    for name, (element, coordinates) in named.items():
        matches = [i for i, a in enumerate(physical) if a['id'].startswith('pqq/')
                   and a['element'] == element and np.max(np.abs(np.array(a['xyz_A'])-coordinates)) < 1e-9]
        if len(matches) != 1:
            raise InvalidArtifact('PQQ source mapping missing or ambiguous')
        pqq_map[name] = matches[0]
    if set(pqq_map.values()) != {i for i, a in enumerate(physical) if a['id'].startswith('pqq/')}:
        raise InvalidArtifact('unaccounted PQQ atom')
    bonds.update(tuple(sorted((pqq_map[a], pqq_map[b]))) for a, b in CCD_SCHEMA.bonds)
    bonds.update(tuple(sorted((pqq_map[h['name']], pqq_map[h['parent_atom']]))) for h in pqq['hydrogens'])
    chain_charge = float(ffcharges.sum())
    if abs(chain_charge-round(chain_charge)) > 1e-6:
        raise InvalidArtifact('nonintegral source protein charge')
    residue_charges = defaultdict(float)
    for a, q in zip(atoms, ffcharges):
        residue_charges[residue_id(physical[a.index]['id'])] += float(q)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    tasks = []
    for metal, formal in [('La', 3), ('Ca', 2)]:
        full = [dict(a, element=metal if a['id']=='metal' else a['element']) for a in physical]
        charge = round(chain_charge) + pqq['formal_charge'] + formal
        electrons = sum(Z[a['element']] for a in full)-charge
        if electrons % 2:
            raise InvalidArtifact('paired closed-shell all-electron parity failed')
        directory = output/metal
        directory.mkdir()
        coordinates = directory/'whole.xyz'
        coordinates.write_text(f'{len(full)}\n1H4I deposited chain A + frozen oxidized PQQ3- + {metal}; no synthetic caps\n' +
                               ''.join(a['element']+' '+' '.join(f'{v:.10f}' for v in a['xyz_A'])+'\n' for a in full))
        inp = directory/'endpoint.inp'
        inp.write_text('! Native-GFN2-xTB ALPB(water) TightSCF\n%scf\n SmearTemp 300\nend\n'
                       '%method\n WriteXTBParam true\nend\n'+f'* xyzfile {charge} 1 whole.xyz\n')
        tasks.append({'task_id': f'1h4i_global_{metal}', 'metal': metal, 'charge': charge, 'multiplicity': 1,
                      'all_electron_count_for_parity_only': electrons,
                      'input': record(inp), 'xyz': record(coordinates),
                      'output_path': str(directory/'endpoint.out')})
    report = {'protocol_id': PROTOCOL, 'status': 'prepared_proposed_not_executed',
              'source': s['source'], 'source_carve': record(carve_path), 'source_states': provenance,
              'implementation': record(__file__), 'orca': record(orca), 'assembly': s['assembly'],
              'execution_policy': {'task_runner': record(Path(__file__).with_name('run_orca_task_manifest.py')),
                  'runtime_renderer': record(Path(__file__).with_name('render_orca_runtime_input.py'))},
              'preparation_dependencies': {name: record(Path(__file__).with_name(name)) for name in
                  ('affordable_state.py', 'affordable_environment.py', 'affordable_common.py', 'pqq_microstates.py')},
              'excluded_deposited_chains': excluded, 'use': 'development_already_consumed_1H4I',
              'microstate': s['microstate'], 'explicit_waters': [],
              'physical_atoms': physical, 'physical_boundary_hash': physical_boundary_key(s['physical_atoms']),
              'bonds_zero_based': sorted(bonds), 'bond_orders': 'not assigned; source connectivity only',
              'pqq_source_name_to_global_index': pqq_map, 'synthetic_caps': [], 'cut_bonds': [],
              'protein_charge_e': round(chain_charge), 'pqq_charge_e': pqq['formal_charge'],
              'source_residue_formal_charge_ledger_e': {k: round(v,8) for k,v in residue_charges.items()},
              'ff_charge_use': 'source protonation/formal-charge audit only; no permanent FF charges in the proposed electronic model',
              'source_coordinate_max_difference_A': max_source_error,
              'tasks': tasks, 'settings': {'method': 'ORCA native GFN2-xTB', 'solvent': 'ALPB water',
                  'scf': 'TightSCF, default native xTB mixer', 'electronic_temperature_K': 300,
                  'geometry': 'frozen', 'parameters': 'installed ORCA defaults; export requested',
                  'reference': None, 'classification_threshold': None},
              'score': None, 'raw_contrast_expression': 'E_Ca(full)-E_La(full); no inherited aquo offset or bands',
              'new_evaluations': 0, 'wall_seconds': time.monotonic()-start}
    report['scientific_cache_key'] = cache_key({k: report[k] for k in (
        'protocol_id','source','assembly','microstate','explicit_waters','physical_atoms','tasks','settings','implementation','orca','preparation_dependencies')})
    write_new(output/'global_manifest.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['audit', 'prepare'])
    parser.add_argument('--solver-result', type=Path, required=True)
    parser.add_argument('--carve-manifest', type=Path)
    parser.add_argument('--orca', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.operation == 'audit':
        report = audit(args.solver_result, args.output)
        print({k: v['direct_score_kcal_mol'] for k,v in report['partitions'].items()})
    else:
        if args.carve_manifest is None or args.orca is None:
            parser.error('prepare requires --carve-manifest and --orca')
        report = prepare(args.solver_result, args.carve_manifest, args.output, args.orca)
        print({'status': report['status'], 'atoms': len(report['physical_atoms']),
               'bonds': len(report['bonds_zero_based']), 'charges': [t['charge'] for t in report['tasks']]})


if __name__ == '__main__':
    main()
