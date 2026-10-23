"""Read the already printed T0/T1 accounting; do not run or re-score a method."""
import hashlib
import json
import re
import time
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_reference_uncertainty_v1/control_qualification_review_v1/actual_v1'


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    wall, cpu = time.monotonic(), time.process_time()
    out = BASE / 'TRIPLES_DETAIL.json'
    assert not out.exists()
    components = json.loads((BASE / 'COMPONENTS.json').read_text())
    k = Decimal(components['hartree_to_kcal_mol'])
    rows, groups = {}, {}
    for cid, row in components['rows'].items():
        binding = row['local']['output']
        raw = Path(binding['path']).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == binding['sha256']
        text = raw.decode()
        fields = {}
        for name, label in [('conventional_T0', 'Conventional (T0)'), ('truncated_T0', 'Truncated    (T0)'), ('corrected_T', 'corrected     (T)')]:
            hits = [(n, m.group(1)) for n, line in enumerate(text.splitlines(), 1)
                    if (m := re.match(r'\s*The\s+' + re.escape(label) + r'\s+is\s+([-+\d.]+)\s*$', line))]
            assert len(hits) == 1, (cid, name, hits)
            fields[name] = {'line': hits[0][0], 'hartree_token': hits[0][1]}
        values = {name: Decimal(q['hartree_token']) for name, q in fields.items()}
        assert abs(values['corrected_T'] - Decimal(row['local']['values']['triples'])) <= Decimal('0.00000000051')
        canonical = Decimal(row['canonical']['values']['triples'])
        shifts = {name: (v - canonical) * k for name, v in values.items()}
        rows[cid] = {'output': binding, 'printed_triples_components': fields,
                     'difference_from_canonical_triples_kcal_mol': {name: float(v) for name, v in shifts.items()}}
        group = groups.setdefault(row['donor'] + '__' + row['metal'], {})
        sign = 1 if row['role'] == 'complex' else -1
        for name, v in shifts.items():
            group[name] = group.get(name, Decimal(0)) + sign * v
    result = {'status': 'saved_triples_detail_complete', 'implementation': pin(__file__),
              'components': pin(BASE / 'COMPONENTS.json'), 'rows': rows,
              'CP_triples_difference_from_canonical_kcal_mol': {key: {name: float(v) for name, v in q.items()} for key, q in groups.items()},
              'interpretation': 'Native Conventional T0 is a local semicanonical intermediate, not the archived canonical CCSD(T) result. Corrected T is the actual iterative T1 contribution. No alternative method result or new gate is admitted.',
              'native_calls': 0, 'cost': {'wall_seconds': time.monotonic()-wall, 'cpu_seconds': time.process_time()-cpu, 'threads': 1}}
    out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'output': pin(out), 'CP_triples_difference_from_canonical_kcal_mol': result['CP_triples_difference_from_canonical_kcal_mol'], 'cost': result['cost']}, indent=2))


if __name__ == '__main__':
    main()
