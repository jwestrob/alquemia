#!/usr/bin/env python3
"""Collect the ten completed DFT cells of the frozen22-cell donor panel.

Saved-only subset of the original declared comparison: replay all native
admissions, report both donor CP interaction/exchange contrasts, and explicitly
leave the correlated-method comparison unavailable. No molecular calls.
"""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('donor_cp', HERE / 'REFERENCE.py')
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
b = ref.b
BASE = b.ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/donor_cp_v1'


def main():
    path = BASE / 'prepared_v1/PANEL.json'
    panel = ref.load(path)
    cells = {c['cell_id']: c for c in panel['cells'] if c['method'] == 'pbe0'}
    out = BASE / 'dft_subset_collection_v1'
    out.mkdir(exist_ok=False)
    rows = {}
    pins = []
    for cell_id, cell in cells.items():
        filename = BASE / 'dft_fragments_v1/tasks' / cell_id / 'execution/RESULT.json'
        result = json.loads(filename.read_text())
        b.check(result['status'] == 'complete' and result['returncode'] == 0, 'Unsuccessful result')
        b.check(result['panel'] == b.pin(path) and result['adapter'] == panel['adapter'], 'Binding mismatch')
        b.check(result['cell'] == cell, 'Cell mismatch')
        b.verify(result['input'])
        replay = ref.parse(b.verify(result['output']).read_text(), cell, result['returncode'])
        b.check(replay['status'] == 'complete' and replay['checks'] == result['checks'], 'Native admission mismatch')
        for key in ('energy_hartree', 'electronic_energy_hartree', 'D4_energy_hartree'):
            b.check(replay[key] == result[key], 'Native energy replay mismatch')
        rows[cell_id] = replay
        pins.append(b.pin(filename))
    b.check(len(rows) == 10, 'Incomplete DFT subset')
    water = json.loads(b.verify(panel['water_baseline']).read_text())
    donors = {}
    for donor in ref.DONORS:
        pair = {}
        ligand_d4 = rows[f'pbe0__{donor}__real_ligand_D4']['D4_energy_hartree']
        for metal in ref.METALS:
            prefix = f'pbe0__{donor}__{metal}__'
            saved = panel['saved_complexes'][donor][metal]
            electronic = (saved['electronic_energy_hartree'] - rows[prefix + 'ion_ghost_ligand']['electronic_energy_hartree']
                          - rows[prefix + 'ligand_ghost_ion']['electronic_energy_hartree']) * b.HA
            dispersion = (saved['D4_energy_hartree'] - panel['isolated_ion_D4_hartree'] - ligand_d4) * b.HA
            interaction = {'PBE0_CP_kcal_mol': electronic, 'PBE0_D4_CP_kcal_mol': electronic + dispersion}
            pair[metal] = {'interaction': interaction,
                           'water_to_donor_exchange': {key: val - water['per_metal'][metal][key] for key, val in interaction.items()}}
        donors[donor] = {'per_metal': pair, 'paired_exchange': {
            key: pair['Dy']['water_to_donor_exchange'][key] - pair['La']['water_to_donor_exchange'][key]
            for key in pair['La']['interaction']}}
    result = {'status': 'DFT_subset_complete_correlated_comparison_pending', 'panel': b.pin(path),
              'implementation': b.pin(__file__), 'native_result_pins': pins,
              'water_baseline': panel['water_baseline'], 'donors': donors,
              'new_molecular_calls': 0, 'physical_accuracy_qualified': False,
              'force_targets_changed': False, 'complete_full22_comparison': False, 'utc': b.now()}
    b.write(out / 'RESULT.json', result)
    print(json.dumps({'result': str(out / 'RESULT.json'), 'donors': donors}, indent=2))


if __name__ == '__main__':
    main()
