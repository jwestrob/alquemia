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
    review = json.loads(readpin(cfg['independent_review']))
    assert review['status'] == 'passed_independent_actual_common_basis_review'
    assert review['collector_result'] == cfg['collection']
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
    expected = set(GIDS)
    for n in (8, 9):
        for metal in ('La', 'Dy'):
            prefix = f'pbe0__aquo{n}__{metal}'
            expected.update(prefix + '__' + role for role in ('ion_ghost_shell', 'shell_ghost_ion'))
            expected.update(prefix + f'__water_{i:02d}__' + role for i in range(1, n + 1) for role in ('metal_water', 'water_only'))
    check('complete_exact80_inventory', len(sources) == 80 and set(sources) == expected)
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
        # Establish the printed origin's units from its equality to the native
        # mass-weighted A.U. coordinates, including zero-mass ghost centers.
        origin_hits = re.findall(r'^Position of electric origin\s*\.\.\.\s*(' + NUMBER + r')\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s*$', text, re.M)
        check(cid + ':one_origin', len(origin_hits) == 1)
        origin = list(map(float, origin_hits[0]))
        coord_block = text.split('CARTESIAN COORDINATES (A.U.)')
        check(cid + ':one_au_coordinate_block', len(coord_block) == 2)
        coords = []
        for line in coord_block[1].splitlines()[1:]:
            match = re.match(r'^\s*(\d+)\s+([A-Za-z]+)\s+(' + NUMBER + r')\*?\s+\d+\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s*$', line)
            if match:
                coords.append((int(match[1]), match[2], float(match[3]), float(match[4]), [float(match[i]) for i in (5, 6, 7)]))
            elif coords:
                break
        check(cid + ':au_centers', [(i, s) for i, s, z, mass, r in coords] == list(enumerate(labels)))
        masses = sum(x[3] for x in coords)
        check(cid + ':positive_mass', masses > 0)
        com = [sum(mass * xyz[j] for i, s, z, mass, xyz in coords) / masses for j in range(3)]
        # Coordinates/origins print six decimals, masses three. Ghost masses
        # are exactly zero; propagate rounding of the positive printed masses.
        mass_error = sum(5e-4 for x in coords if x[3] > 0)
        check(cid + ':mass_rounding_bound', masses > mass_error)
        com_bounds = [1e-6 + sum(5e-4 * abs(xyz[j] - com[j]) for i, s, z, mass, xyz in coords if mass > 0) / (masses - mass_error) for j in range(3)]
        check(cid + ':origin_is_au_COM', all(abs(a - b) <= bound for a, b, bound in zip(com, origin, com_bounds)))
        electron_hits = re.findall(r'^\s*Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)\s*$', text, re.M)
        check(cid + ':one_electron_count', len(electron_hits) == 1)
        ne = int(electron_hits[0])
        ztot = sum(x[2] for x in coords)
        check(cid + ':effective_nuclear_electron_charge', abs(ztot - ne - expected_q) < 1e-8)
        row['native_electric_origin_au'] = origin
        row['origin_vs_rounded_mass_COM_bounds_au'] = com_bounds
        row['native_electron_count'] = ne
        row['native_effective_nuclear_charge'] = ztot
        row['printed_dipoles_at_native_origin'] = {key: row[key] for key in ('electronic_dipole_au', 'nuclear_dipole_au', 'total_dipole_au')}
        row['dipole_roundoff_bounds_au'] = {}
        for key, charge in (('electronic_dipole_au', -ne), ('nuclear_dipole_au', ztot), ('total_dipole_au', expected_q)):
            row[key] = [mu + charge * o for mu, o in zip(row[key], origin)]
            row['dipole_roundoff_bounds_au'][key] = 5e-10 + abs(charge) * 5e-7
        nuclear_from_coords = [sum(z * xyz[j] for i, s, z, mass, xyz in coords) for j in range(3)]
        check(cid + ':nuclear_moment_and_origin_units', max(abs(a - b) for a, b in zip(nuclear_from_coords, row['nuclear_dipole_au'])) <= row['dipole_roundoff_bounds_au']['nuclear_dipole_au'] + ztot * 5e-7 + 1e-10)
        check(cid + ':converted_dipole_components', max(abs(t - e - z) for t, e, z in zip(row['total_dipole_au'], row['electronic_dipole_au'], row['nuclear_dipole_au'])) <= 5e-9)
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
        parent = {'water_count': n, 'combinations': {'C': c_terms, 'P': p_terms, 'B': b_terms}, 'vectors': {}, 'dipole_roundoff_bounds_au': {}}
        for key in ('mulliken', 'loewdin', 'electronic_dipole_au', 'nuclear_dipole_au', 'total_dipole_au'):
            c, p, b = [combine(terms, key) for terms in (c_terms, p_terms, b_terms)]
            check(gid + ':' + key + ':expanded_closure', max(abs(x - y - z) for x, y, z in zip(c, p, b)) <= 1e-10)
            parent['vectors'][key] = {'C': c, 'P': p, 'B': b}
            if key in ('mulliken', 'loewdin'):
                check(gid + ':' + key + ':neutral_differences', max(abs(sum(v)) for v in (c, p, b)) <= 5e-4)
            if key == 'nuclear_dipole_au':
                for quantity, terms, vec in zip(('C', 'P', 'B'), (c_terms, p_terms, b_terms), (c, p, b)):
                    bound = sum(abs(coef) * rows[cid]['dipole_roundoff_bounds_au'][key] for coef, cid in terms)
                    check(gid + ':' + quantity + ':nuclear_cancellation', max(abs(x) for x in vec) <= bound + 1e-10)
            if key.endswith('_dipole_au'):
                parent['dipole_roundoff_bounds_au'][key] = {q: sum(abs(coef) * rows[cid]['dipole_roundoff_bounds_au'][key] for coef, cid in terms) for q, terms in zip(('C', 'P', 'B'), (c_terms, p_terms, b_terms))}
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
    start, cpu = time.monotonic(), time.process_time()
    output_preexisted = args.output.exists()
    try:
        main(args.config, args.output)
    except Exception as exc:
        if not output_preexisted:
            args.output.mkdir(parents=True, exist_ok=True)
            failure = {'status': 'failed_saved_review', 'implementation': pin(__file__), 'config': pin(args.config) if args.config.is_file() else None, 'error': repr(exc), 'new_native_calls': 0,
                       'cost': {'wall_seconds': time.monotonic() - start, 'CPU_seconds': time.process_time() - cpu, 'local_threads': 1}}
            (args.output / 'FAILURE.json').write_text(json.dumps(failure, indent=2) + '\n')
        raise
