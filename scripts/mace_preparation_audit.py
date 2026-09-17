"""Read-only bonded-geometry audit on a pinned physical protein state."""
import argparse
import csv
import importlib.metadata
from pathlib import Path
import time
import numpy as np
import openmm
from openmm import app, unit
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from affordable_state import FF, source_protein


def audit(state_path, plan, output):
    start = time.monotonic()
    s = read_json(state_path)
    atoms, pos, _, physical, _ = source_protein(verify(s['source']))
    selected = [a for a in s['physical_atoms'] if a.get('source_index') is not None]
    if [a['id'] for a in physical] != [a['id'] for a in selected]:
        raise InvalidArtifact('physical protein atom mapping differs')
    source_error = float(np.max(np.abs(pos - np.array([a['xyz_A'] for a in selected]))))
    if source_error > .002:
        raise InvalidArtifact('physical/source protein geometry mismatch')
    pos = np.array([a['xyz_A'] for a in selected])
    pdb = app.PDBFile(str(verify(s['source'])))
    model = app.Modeller(pdb.topology, pdb.positions)
    allowed_ids = set(a['id'] for a in selected)
    model.delete([a for a in model.topology.atoms() if
                  f'{a.residue.chain.id}/{a.residue.id}/{a.residue.insertionCode}/{a.name}' not in allowed_ids])
    ids = [f'{a.residue.chain.id}/{a.residue.id}/{a.residue.insertionCode}/{a.name}' for a in model.topology.atoms()]
    if ids != [a['id'] for a in selected]:
        raise InvalidArtifact('force-field topology order mismatch')
    system = app.ForceField(str(FF)).createSystem(model.topology, nonbondedMethod=app.NoCutoff,
                                               constraints=None, rigidWater=False)
    bonds = next(f for f in system.getForces() if isinstance(f, openmm.HarmonicBondForce))
    rows = []
    for i in range(bonds.getNumBonds()):
        a, b, length, stiffness = bonds.getBondParameters(i)
        measured = float(np.linalg.norm(pos[a] - pos[b]))
        equilibrium = length.value_in_unit(unit.angstrom)
        rows.append({'atom_a_index': a, 'atom_a_id': ids[a], 'atom_b_index': b, 'atom_b_id': ids[b],
                     'elements': '-'.join(sorted((selected[a]['element'], selected[b]['element']))),
                     'measured_A': measured, 'forcefield_equilibrium_A': equilibrium,
                     'deviation_A': measured - equilibrium,
                     'forcefield_k_kJ_mol_nm2': stiffness.value_in_unit(unit.kilojoule_per_mole / unit.nanometer**2)})
    groups = {}
    for name in sorted({r['elements'] for r in rows}):
        values = np.array([r['deviation_A'] for r in rows if r['elements'] == name])
        groups[name] = {'bonds': len(values), 'mean_deviation_A': float(values.mean()),
                        'RMS_deviation_A': float(np.sqrt(np.square(values).mean())),
                        'max_absolute_deviation_A': float(np.abs(values).max())}
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    table = out / 'protein_bonds.tsv'
    with table.open('x') as f:
        w = csv.DictWriter(f, fieldnames=rows[0], delimiter='\t')
        w.writeheader(); w.writerows(rows)
    result = {'status': 'complete', 'source_state': record(state_path), 'plan': record(plan),
              'implementation': record(__file__), 'protein_helper': record(Path(__file__).with_name('affordable_state.py')),
              'forcefield': record(FF), 'openmm_version': importlib.metadata.version('openmm'),
              'protein_atoms': len(atoms), 'source_coordinate_max_error_A': source_error,
              'bond_groups': groups, 'largest_deviations': sorted(rows, key=lambda r: abs(r['deviation_A']), reverse=True)[:30],
              'all_bonds': record(table), 'wall_seconds': time.monotonic() - start,
              'new_MACE_calls': 0, 'new_DFT_calls': 0,
              'excluded_from_bond_comparison': 'metal and PQQ; no missing charges or parameters invented'}
    write_new(out / 'result.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('state', 'plan', 'output'):
        p.add_argument('--' + key, required=True)
    a = p.parse_args()
    print(audit(a.state, a.plan, a.output)['status'])
