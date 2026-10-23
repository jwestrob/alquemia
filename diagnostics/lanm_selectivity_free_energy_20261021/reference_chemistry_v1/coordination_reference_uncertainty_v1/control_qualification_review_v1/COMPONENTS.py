"""Saved-only component accounting; no native execution or admission changes."""
import hashlib
import json
import re
import time
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_reference_uncertainty_v1'
OUT = BASE / 'control_qualification_review_v1/actual_v1/COMPONENTS.json'
KCAL = Decimal('627.509474063')


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verified_read(binding):
    raw = Path(binding['path']).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == binding['sha256'], binding['path']
    return raw.decode()


def energies(binding):
    text = verified_read(binding)
    patterns = {
        'CCSD': r'^\s*E\(CCSD\)\s+\.\.\.\s+([-+\d.]+)\s*$',
        'triples': r'^\s*Triples Correction \(T\)\s+\.\.\.\s+([-+\d.]+)\s*$',
        'total': r'^\s*E\(CCSD\(T\)\)\s+\.\.\.\s+([-+\d.]+)\s*$',
        'RHF': r'^\s*E\(0\)\s+\.\.\.\s+([-+\d.]+)\s*$',
        'T1_diagnostic': r'^\s*T1 diagnostic\s+\.\.\.\s+([-+\d.]+)\s*$',
    }
    result = {'output': binding, 'evidence': {}}
    values = {}
    for name, pattern in patterns.items():
        hits = [(n, m.group(1)) for n, line in enumerate(text.splitlines(), 1)
                if (m := re.match(pattern, line))]
        assert hits and len({v for _, v in hits}) == 1, (binding, name, hits)
        values[name] = Decimal(hits[-1][1])
        result['evidence'][name] = [{'line': n, 'token': v} for n, v in hits]
    closure = values['total'] - values['CCSD'] - values['triples']
    assert abs(closure) <= Decimal('0.000000002'), (binding, closure)
    result['printed_component_closure_hartree'] = str(closure)
    result['values'] = {k: str(v) for k, v in values.items()}
    return result, values


def main():
    wall, cpu = time.monotonic(), time.process_time()
    assert not OUT.exists(), 'Preserve completed component accounting'
    review_path = BASE / 'control_qualification_review_v1/actual_v1/REVIEW.json'
    review = json.loads(review_path.read_text())
    result = json.loads(verified_read(review['collector_result']))
    panel = json.loads(verified_read(review['panel']))
    cells = [c for c in panel['cells'] if c['stage'] == 'controls']
    assert len(cells) == len(result['rows']) == 18
    rows, groups = {}, {}
    for c in cells:
        cid = c['cell_id']
        old = panel['controls'][c['canonical_id']]
        new = result['rows'][cid]
        canonical, eold = energies(old['output'])
        local, enew = energies(new['output'])
        assert abs(enew['total'] - Decimal(str(new['energy_hartree']))) <= Decimal('0.000000001')
        assert abs(eold['total'] - Decimal(str(old['energy_hartree']))) <= Decimal('0.000000001')
        shifts = {k: (enew[k] - eold[k]) * KCAL for k in ('RHF', 'CCSD', 'triples', 'total')}
        assert shifts['RHF'] == 0
        rows[cid] = {'metal': c['metal'], 'donor': c['donor'], 'role': c['role'],
                     'canonical': canonical, 'local': local,
                     'local_minus_canonical_kcal_mol': {k: float(v) for k, v in shifts.items()}}
        group = groups.setdefault(c['donor'] + '__' + c['metal'], {})
        sign = 1 if c['role'] == 'complex' else -1
        for k, v in shifts.items():
            group[k] = group.get(k, Decimal(0)) + sign * v
    cp = {}
    for key, group in groups.items():
        target = Decimal(str(result['CP_errors_kcal_mol'][key]))
        assert abs(group['total'] - target) < Decimal('0.0000001'), (key, group, target)
        cp[key] = {k: float(v) for k, v in group.items()}
        cp[key]['component_closure_kcal_mol'] = float(group['total'] - group['CCSD'] - group['triples'])
    paired = {donor: {k: float(groups[donor + '__Dy'][k] - groups[donor + '__La'][k])
                      for k in ('RHF', 'CCSD', 'triples', 'total')}
              for donor in ('water', 'acetamide', 'acetate_bi')}
    receipt = {
        'status': 'saved_component_accounting_complete', 'implementation': pin(__file__),
        'review': pin(review_path), 'collector_result': review['collector_result'],
        'panel': review['panel'], 'rows': rows, 'complete_local': 18, 'complete_canonical': 18,
        'formula': 'CP = complex - ion_in_full_basis - ligand_in_full_basis; error = local - canonical',
        'hartree_to_kcal_mol': str(KCAL), 'CP_errors_kcal_mol': cp,
        'paired_CP_errors_kcal_mol': paired, 'all_RHF_differences_exactly_zero_at_printed_precision': True,
        'gate_changes': False, 'native_calls': 0,
        'interpretation_limit': 'CCSD/RI/locality and triples differences are accounting components, not individually identified causes; no tighter-threshold outcome is inferred.',
        'cost': {'wall_seconds': time.monotonic() - wall, 'cpu_seconds': time.process_time() - cpu, 'threads': 1},
    }
    OUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'output': pin(OUT), 'CP_errors_kcal_mol': cp, 'paired': paired, 'cost': receipt['cost']}, indent=2))


if __name__ == '__main__':
    main()
