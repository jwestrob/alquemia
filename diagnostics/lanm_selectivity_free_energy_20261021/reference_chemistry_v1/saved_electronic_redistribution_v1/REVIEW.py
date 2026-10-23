#!/usr/bin/env python3
"""Saved population/dipole differences; never invokes a molecular program."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time

NUMBER = r'[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?'
GIDS = ('aquo8_central_La', 'aquo8_central_Dy', 'aquo9_central_La', 'aquo9_central_Dy')


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def readpin(value):
    assert pin(value['path']) == value, 'Changed input: ' + value['path']
    return Path(value['path']).read_text()


def main(config_path, output):
    start, cpu = time.monotonic(), time.process_time()
    assert all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'))
    cfg = json.loads(config_path.read_text())
    assert cfg['implementation'] == pin(__file__)
    collection = json.loads(readpin(cfg['collection']))
    readpin(cfg['independent_review'])
    assert collection['status'] == 'common_basis_collection_complete'
    assert (collection['new_complete'], collection['reused_complete']) == (68, 12)
    assert collection['error'] is None
    assert not output.exists(), 'Preserve previous execution'
    output.mkdir(parents=True)
    checks = []

    def check(name, condition):
        assert condition, name
        checks.append(name)

    rows = {}
    sources = dict(collection['reused_terms'], **collection['new_fragments'])
    for cid, source in sources.items():
        text = readpin(source['output'])
        n = 8 if 'aquo8' in cid else 9
        metal = 'La' if 'La' in cid else 'Dy'
        labels = [metal] + ['O', 'H', 'H'] * n
        expected_q = 0 if 'water_only' in cid or 'shell_ghost_ion' in cid else 3
        row = {'output': source['output'], 'expected_total_charge_e': expected_q, 'labels': labels}
        for scheme in ('MULLIKEN', 'LOEWDIN'):
            headers = list(re.finditer(r'^' + scheme + r' ATOMIC CHARGES\s*$', text, re.M))
            check(cid + ':' + scheme + ':one_table', len(headers) == 1)
            tail = text[headers[0].end():].lstrip('\n- ')
            table = []
            for line in tail.splitlines():
                match = re.match(r'^\s*(\d+)\s+([A-Za-z]+)\s*:\s*(' + NUMBER + r')\s*$', line)
                if not match:
                    break
                table.append((int(match[1]), match[2], float(match[3])))
            check(cid + ':' + scheme + ':centers', [(i, s) for i, s, q in table] == list(enumerate(labels)))
            values = [q for i, s, q in table]
            check(cid + ':' + scheme + ':finite_sum', all(math.isfinite(q) for q in values) and abs(sum(values) - expected_q) <= 5e-5)
            row[scheme.lower()] = values
        for key, label in (('electronic_dipole_au', 'Electronic contribution'), ('nuclear_dipole_au', 'Nuclear contribution'), ('total_dipole_au', 'Total Dipole Moment')):
            hits = re.findall(r'^\s*' + label + r'\s*:\s*(' + NUMBER + r')\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s*$', text, re.M)
            check(cid + ':' + key + ':one_vector', len(hits) == 1)
            row[key] = list(map(float, hits[0]))
            check(cid + ':' + key + ':finite', all(math.isfinite(x) for x in row[key]))
        check(cid + ':dipole_components', max(abs(t - e - z) for t, e, z in zip(row['total_dipole_au'], row['electronic_dipole_au'], row['nuclear_dipole_au'])) <= 5e-9)
        rows[cid] = row

    def combine(terms, key):
        size = len(rows[terms[0][1]][key])
        return [sum(coef * rows[cid][key][i] for coef, cid in terms) for i in range(size)]

    parents = {}
    for gid in GIDS:
        n = 8 if 'aquo8' in gid else 9
        metal = gid.rsplit('_', 1)[1]
        prefix = f'pbe0__aquo{n}__{metal}'
        ion, shell = prefix + '__ion_ghost_shell', prefix + '__shell_ghost_ion'
        dimers = [prefix + f'__water_{i:02d}__metal_water' for i in range(1, n + 1)]
        waters = [prefix + f'__water_{i:02d}__water_only' for i in range(1, n + 1)]
        c_terms = [(1, gid), (-1, ion), (-1, shell)]
        p_terms = [(1, x) for x in dimers] + [(-n, ion)] + [(-1, x) for x in waters]
        b_terms = [(1, gid), (n - 1, ion), (-1, shell)] + [(-1, x) for x in dimers] + [(1, x) for x in waters]
        parent = {'water_count': n, 'combinations': {'C': c_terms, 'P': p_terms, 'B': b_terms}, 'vectors': {}}
        for key in ('mulliken', 'loewdin', 'electronic_dipole_au', 'nuclear_dipole_au', 'total_dipole_au'):
            c, p, b = [combine(terms, key) for terms in (c_terms, p_terms, b_terms)]
            check(gid + ':' + key + ':expanded_closure', max(abs(x - y - z) for x, y, z in zip(c, p, b)) <= 1e-10)
            parent['vectors'][key] = {'C': c, 'P': p, 'B': b}
            if key in ('mulliken', 'loewdin'):
                check(gid + ':' + key + ':neutral_differences', max(abs(sum(v)) for v in (c, p, b)) <= 5e-4)
            if key == 'nuclear_dipole_au':
                check(gid + ':nuclear_cancellation', max(abs(x) for v in (c, p, b) for x in v) <= 1e-6)
        parent['metal_population_charge_differences_e'] = {s: {k: v[0] for k, v in parent['vectors'][s].items()} for s in ('mulliken', 'loewdin')}
        parents[gid] = parent

    paired = {}
    for scheme in ('mulliken', 'loewdin'):
        paired[scheme] = {q: sum(sign * parents[gid]['metal_population_charge_differences_e'][scheme][q] for sign, gid in zip((1, -1, -1, 1), GIDS)) for q in ('C', 'P', 'B')}
    result = {'status': 'complete_saved_redistribution_diagnostic', 'implementation': pin(__file__), 'config': pin(config_path), 'rows': rows,
              'parents': parents, 'paired_metal_population_charge_differences_e': paired, 'checks': checks,
              'new_native_calls': 0, 'unique_charge_transfer_identified': False, 'model_change_released': False,
              'cost': {'wall_seconds': time.monotonic() - start, 'CPU_seconds': time.process_time() - cpu, 'local_threads': 1}}
    (output / 'RESULT.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = ['# Saved electronic redistribution', '', f'All {len(rows)} raw outputs and {len(checks)} reconstruction checks pass.', '',
             'Population differences are elementary charges on the metal center. All ghost-center populations are included in the difference construction.', '',
             '| Parent | Mulliken C | Mulliken P | Mulliken B | Loewdin C | Loewdin P | Loewdin B |', '|---|---:|---:|---:|---:|---:|---:|']
    for gid, parent in parents.items():
        v = parent['metal_population_charge_differences_e']
        lines.append('| ' + gid + ' | ' + ' | '.join(f'{v[s][q]:+.6f}' for s in ('mulliken', 'loewdin') for q in ('C', 'P', 'B')) + ' |')
    lines += ['', 'Positive B means less metal-centered electronic population gain for the full shell than the additive pair construction under that population scheme. It is not a unique charge-transfer observable or an energy attribution.', '',
              'Dipole vectors and all per-center populations are retained in RESULT.json. Nuclear cancellation passes; the neutral difference combinations are independent of a common origin shift. Cross-parent vector comparisons additionally depend on coordinate orientation and are not interpreted here.', '',
              'Original ECP/PBE0 state, basis and SCF limitations remain. No new molecular calls, fitting or model selection.']
    (output / 'REPORT.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps({'status': result['status'], 'rows': len(rows), 'checks': len(checks), 'cost': result['cost']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    main(args.config, args.output)
