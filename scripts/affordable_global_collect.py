"""Collect the manifested whole-protein pilot, keeping failed endpoints visible."""
import argparse
from collections import defaultdict
import math
from pathlib import Path
import re

from affordable_common import InvalidArtifact, energy, read_json, record, verify, write_new, xyz, HA_TO_KCAL


def mulliken(output, coordinates, charge):
    text = Path(output).read_text()
    tables = re.findall(r'MULLIKEN ATOMIC CHARGES\s*\n-+\n(.*?)Sum of atomic charges:\s*([-+\d.]+)', text, re.S)
    if len(tables) != 1:
        raise InvalidArtifact('missing/ambiguous closed-shell Mulliken table')
    table, total = tables[0]
    rows = re.findall(r'^\s*(\d+)\s+([A-Z][a-z]?)\s*:\s*([-+\d.]+)\s*$', table, re.M)
    atoms = xyz(coordinates)
    if len(rows) != len(atoms) or [int(r[0]) for r in rows] != list(range(len(atoms))):
        raise InvalidArtifact('incomplete Mulliken atom indexing')
    if [r[1] for r in rows] != [a[0] for a in atoms]:
        raise InvalidArtifact('Mulliken/XYZ ordering mismatch')
    charges = [float(r[2]) for r in rows]
    # The table has six decimal places; account only for its maximum rounding
    # error, not an adjustable scientific charge-closure tolerance.
    if not all(math.isfinite(q) for q in charges) or abs(float(total)-charge)>1e-6 or abs(sum(charges)-charge) > len(charges)*0.5e-6+1e-6:
        raise InvalidArtifact('printed Mulliken charge closure failed')
    return {'charge_e': charges, 'printed_sum_e': float(total), 'sum_of_printed_atoms_e': sum(charges),
            'representation': 'native-model Mulliken population; not ORCA MBIS',
            'printed_atom_precision_decimal_places': 6}


def collect(manifest_path, destination):
    m = read_json(manifest_path)
    verify(m['orca'])
    rows = []
    for task in m['tasks']:
        row = {'metal': task['metal'], 'task_id': task['task_id'], 'energy_hartree': None,
               'status': 'unavailable', 'population': None}
        try:
            coordinates = verify(task['xyz']); verify(task['input'])
            path = Path(task['output_path'])
            receipt_path = path.with_name(path.name+'.execution.json')
            receipt = read_json(receipt_path)
            verify(receipt['manifest'])
            if receipt['manifest'] != record(manifest_path) or receipt['orca_executable'] != m['orca']:
                raise InvalidArtifact('execution manifest or executable mismatch')
            for artifact in receipt['artifacts'].values():
                if artifact is not None: verify(artifact)
            row.update(execution_receipt=record(receipt_path), output=record(path),
                       parallelism=receipt['parallelism'], allocation=receipt['allocation'],
                       started_at_utc=receipt['started_at_utc'], finished_at_utc=receipt['finished_at_utc'])
            row['energy_hartree'] = energy(path)
            row['status'] = 'energy_complete_population_unavailable'
            population = mulliken(path, coordinates, task['charge'])
            grouped = defaultdict(float)
            population['atoms'] = []
            for atom, q in zip(m['physical_atoms'], population['charge_e']):
                ident = atom['id']
                group = 'pqq' if ident.startswith('pqq/') else '/'.join(ident.split('/')[:3])
                grouped[group] += q
                population['atoms'].append({'id': ident, 'charge_e': q})
            population['residue_charge_e'] = dict(grouped)
            row.update(status='complete', population=population)
        except (InvalidArtifact, FileNotFoundError, KeyError, ValueError) as exc:
            row['reason'] = str(exc)
        rows.append(row)
    result = {'protocol_id': m['protocol_id'], 'manifest': record(manifest_path), 'implementation': record(__file__),
              'rows': rows, 'R_hartree': None, 'R_kcal_mol': None, 'score': None,
              'decision': 'uncalibrated_global_protocol', 'scientific_validation': 'not_established',
              'relative_population_Ca_minus_La_e': None}
    by_metal = {r['metal']:r for r in rows}
    if all(r['energy_hartree'] is not None for r in rows):
        result['R_hartree'] = by_metal['Ca']['energy_hartree']-by_metal['La']['energy_hartree']
        result['R_kcal_mol'] = result['R_hartree']*HA_TO_KCAL
    if all(r['status']=='complete' for r in rows):
        ca, la = (by_metal[metal]['population']['residue_charge_e'] for metal in ('Ca','La'))
        if ca.keys() != la.keys(): raise InvalidArtifact('population source identities differ')
        result['relative_population_Ca_minus_La_e'] = {k: ca[k]-la[k] for k in la}
    write_new(destination, result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = collect(args.manifest, args.output)
    print({r['metal']:r['status'] for r in result['rows']})
