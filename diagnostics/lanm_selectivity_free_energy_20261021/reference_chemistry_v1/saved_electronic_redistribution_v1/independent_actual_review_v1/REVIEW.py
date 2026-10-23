#!/usr/bin/env python3
"""Independent saved-only population and full-precision dipole replay."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time

GIDS = ('aquo8_central_La', 'aquo8_central_Dy', 'aquo9_central_La', 'aquo9_central_Dy')
N = r'[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?'


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def load(pin_value):
    assert pin(pin_value['path']) == pin_value, pin_value
    return json.loads(Path(pin_value['path']).read_text())


def section(text, name):
    found = re.findall(r'^\$' + re.escape(name) + r'\s*\n(.*?)^\$End\s*$', text, re.M | re.S)
    assert len(found) == 1, name
    return found[0]


def array(block, name, length):
    tail = block.split('&' + name + ' ', 1)
    assert len(tail) == 2, name
    body = tail[1].split('\n   &', 1)[0]
    entries = re.findall(r'^\s*(\d+)\s+(' + N + r')\s*$', body, re.M)
    assert [int(i) for i, _ in entries] == list(range(length)), name
    values = [float(v) for _, v in entries]
    assert all(math.isfinite(v) for v in values)
    return values


def maxdiff(a, b):
    assert len(a) == len(b)
    return max(abs(x - y) for x, y in zip(a, b))


def main(cfg_path, out):
    started, cpu = time.monotonic(), time.process_time()
    assert all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'))
    cfg = json.loads(cfg_path.read_text())
    assert cfg['implementation'] == pin(__file__)
    actual = load(cfg['actual_result'])
    root_cfg = load(actual['config'])
    assert actual['implementation'] == root_cfg['implementation']
    load(root_cfg['source_review'])
    assert pin(root_cfg['implementation']['path']) == root_cfg['implementation']
    collection = load(root_cfg['collection'])
    independent_collection = load(root_cfg['independent_review'])
    assert independent_collection['status'] == 'passed_independent_actual_common_basis_review'
    assert independent_collection['collector_result'] == root_cfg['collection']
    assert actual['status'] == 'complete_saved_redistribution_diagnostic'
    assert actual['new_native_calls'] == 0 and len(actual['checks']) == 1881
    sources = dict(collection['reused_terms'], **collection['new_fragments'])
    expected = set(GIDS)
    for n in (8, 9):
        for metal in ('La', 'Dy'):
            prefix = f'pbe0__aquo{n}__{metal}'
            expected.update(prefix + '__' + x for x in ('ion_ghost_shell', 'shell_ghost_ion'))
            expected.update(prefix + f'__water_{i:02d}__' + role for i in range(1, n + 1) for role in ('metal_water', 'water_only'))
    assert set(sources) == set(actual['rows']) == expected and len(expected) == 80
    assert not out.exists(), 'Preserve prior review'
    out.mkdir(parents=True)
    checks, rows = [], {}
    maxima = {'raw_population_replay_e': 0.0, 'property_population_printing_e': 0.0,
              'property_dipole_printing_au': 0.0, 'nuclear_COM_reconstruction_au': 0.0,
              'row_total_dipole_vs_printed_au': 0.0, 'parent_total_dipole_vs_printed_au': 0.0,
              'total_dipole_bound_fraction': 0.0}

    def check(name, passed):
        assert passed, name
        checks.append(name)

    for cid in sorted(expected):
        old = actual['rows'][cid]
        check(cid + ':raw_authority', sources[cid]['output'] == old['output'] == pin(old['output']['path']))
        raw_path = Path(old['output']['path'])
        text = raw_path.read_text()
        prop_path = raw_path.with_name('endpoint.property.txt')
        prop = prop_path.read_text()
        check(cid + ':property_version', 'ORCA 6.1.1' in prop)
        labels = old['labels']
        count = len(labels)
        row = {'output': old['output'], 'property_file': pin(prop_path)}
        # Read every indexed population, including all ghost centers.
        for scheme, title in (('mulliken', 'MULLIKEN'), ('loewdin', 'LOEWDIN')):
            blocks = re.split(r'^' + title + r' ATOMIC CHARGES\s*$', text, flags=re.M)
            check(cid + ':' + scheme + ':one_raw_table', len(blocks) == 2)
            values = []
            for line in blocks[1].lstrip('\n- ').splitlines():
                m = re.fullmatch(r'\s*(\d+)\s+([A-Za-z]+)\s*:\s*(' + N + r')\s*', line)
                if not m:
                    break
                check(cid + ':' + scheme + ':center' + m[1], int(m[1]) == len(values) and m[2] == labels[len(values)])
                values.append(float(m[3]))
            check(cid + ':' + scheme + ':complete_raw_vector', len(values) == count and values == old[scheme])
            check(cid + ':' + scheme + ':original_charge_sum', abs(sum(values) - old['expected_total_charge_e']) <= 5e-5)
            props = section(prop, 'SCF_' + ('Mulliken' if scheme == 'mulliken' else 'Loewdin') + '_Population_Analysis')
            full = array(props, 'AtomicCharges', count)
            delta = maxdiff(full, values)
            check(cid + ':' + scheme + ':property_rounds_to_raw', delta <= 0.5e-6)
            maxima['property_population_printing_e'] = max(maxima['property_population_printing_e'], delta)
            row[scheme] = values
            row[scheme + '_property_precision'] = full

        native_geometry = text.split('CARTESIAN COORDINATES (A.U.)')
        check(cid + ':one_native_AU_geometry', len(native_geometry) == 2)
        atoms = []
        for line in native_geometry[1].splitlines():
            m = re.fullmatch(r'\s*(\d+)\s+([A-Za-z]+)\s+(' + N + r')\*?\s+\d+\s+(' + N + r')\s+(' + N + r')\s+(' + N + r')\s+(' + N + r')\s*', line)
            if m:
                atoms.append((int(m[1]), m[2], float(m[3]), float(m[4]), [float(m[i]) for i in (5, 6, 7)]))
            elif atoms:
                break
        check(cid + ':native_centers', [(x[0], x[1]) for x in atoms] == list(enumerate(labels)))
        geometry = section(prop, 'Geometry')
        check(cid + ':explicit_property_Bohr_units', '&Units "Bohr"' in geometry)
        positions = []
        for line in geometry.split('&Units "Bohr"]', 1)[1].splitlines():
            m = re.fullmatch(r'\s*([A-Za-z]+)\s+(' + N + r')\s+(' + N + r')\s+(' + N + r')\s*', line)
            if m:
                check(cid + ':property_center' + str(len(positions)), m[1] == labels[len(positions)])
                positions.append([float(m[i]) for i in (2, 3, 4)])
        check(cid + ':property_geometry_complete', len(positions) == count)
        check(cid + ':property_geometry_rounds_to_raw', max(maxdiff(x[4], r) for x, r in zip(atoms, positions)) <= 0.5e-6)
        mass = sum(x[3] for x in atoms)
        com = [math.fsum(x[3] * r[j] for x, r in zip(atoms, positions)) / mass for j in range(3)]
        za = sum(x[2] for x in atoms)
        ne = old['native_electron_count']
        q = old['expected_total_charge_e']
        check(cid + ':effective_charge', za == old['native_effective_nuclear_charge'] == ne + q)
        check(cid + ':ghost_charge_mass', all((x[2] == 0) == (x[3] == 0) for x in atoms))
        check(cid + ':COM_rounds_to_printed_origin', maxdiff(com, old['native_electric_origin_au']) <= 0.5e-6)
        dipole = section(prop, 'SCF_Dipole_Moment')
        full_native = {}
        row['dipole_roundoff_bounds_au'] = old['dipole_roundoff_bounds_au']
        for key, field, charge in (('total_dipole_au', 'dipoleTotal', q), ('nuclear_dipole_au', 'dipoleNucContrib', za), ('electronic_dipole_au', 'dipoleElecContrib', -ne)):
            full = array(dipole, field, 3)
            delta = maxdiff(full, old['printed_dipoles_at_native_origin'][key])
            check(cid + ':' + key + ':property_rounds_to_raw', delta <= 0.5e-9)
            maxima['property_dipole_printing_au'] = max(maxima['property_dipole_printing_au'], delta)
            full_native[key] = full
            row[key] = [full[j] + charge * com[j] for j in range(3)]
            difference = maxdiff(row[key], old[key])
            bound = old['dipole_roundoff_bounds_au'][key]
            check(cid + ':' + key + ':within_existing_printed_bound', difference <= bound)
            if key == 'total_dipole_au':
                maxima['row_total_dipole_vs_printed_au'] = max(maxima['row_total_dipole_vs_printed_au'], difference)
                maxima['total_dipole_bound_fraction'] = max(maxima['total_dipole_bound_fraction'], difference / bound)
        nuclear_prediction = [math.fsum(x[2] * (r[j] - com[j]) for x, r in zip(atoms, positions)) for j in range(3)]
        nuclear_error = maxdiff(nuclear_prediction, full_native['nuclear_dipole_au'])
        maxima['nuclear_COM_reconstruction_au'] = max(maxima['nuclear_COM_reconstruction_au'], nuclear_error)
        # Geometry prints twelve decimals: two position terms yield ZA*1e-12.
        check(cid + ':nuclear_COM_identity', nuclear_error <= za * 1e-12 + 1e-12)
        row.update({'COM_Bohr': com, 'native_effective_nuclear_charge': za, 'native_electrons': ne, 'charge': q,
                    'nuclear_COM_reconstruction_error_au': nuclear_error, 'property_dipoles_at_native_origin_au': full_native})
        rows[cid] = row

    parents = {}
    for gid in GIDS:
        n = 8 if 'aquo8' in gid else 9
        prefix = 'pbe0__aquo' + str(n) + '__' + gid.rsplit('_', 1)[1]
        ion, shell = prefix + '__ion_ghost_shell', prefix + '__shell_ghost_ion'
        pairs = [prefix + f'__water_{i:02d}__metal_water' for i in range(1, n + 1)]
        waters = [prefix + f'__water_{i:02d}__water_only' for i in range(1, n + 1)]
        terms = {'C': [(1, gid), (-1, ion), (-1, shell)],
                 'P': [(1, x) for x in pairs] + [(-n, ion)] + [(-1, x) for x in waters],
                 'B': [(1, gid), (n-1, ion), (-1, shell)] + [(-1, x) for x in pairs] + [(1, x) for x in waters]}
        check(gid + ':coefficients', {k: [list(t) for t in v] for k, v in terms.items()} == actual['parents'][gid]['combinations'])
        parent = {'vectors': {}, 'metal_charge_e': {}, 'dipole_comparison': {}}
        for key in ('mulliken', 'loewdin', 'total_dipole_au', 'nuclear_dipole_au', 'electronic_dipole_au'):
            vectors = {k: [math.fsum(a * rows[cid][key][i] for a, cid in terms[k]) for i in range(len(rows[gid][key]))] for k in terms}
            parent['vectors'][key] = vectors
            check(gid + ':' + key + ':expanded_algebra', max(abs(c-p-b) for c, p, b in zip(vectors['C'], vectors['P'], vectors['B'])) <= 1e-10)
            if key in ('mulliken', 'loewdin'):
                for k in terms:
                    error = maxdiff(vectors[k], actual['parents'][gid]['vectors'][key][k])
                    check(gid + ':' + key + ':' + k + ':all_center_replay', error <= 1e-10)
                    check(gid + ':' + key + ':' + k + ':original_neutral_sum', abs(sum(vectors[k])) <= 5e-4)
                    maxima['raw_population_replay_e'] = max(maxima['raw_population_replay_e'], error)
                parent['metal_charge_e'][key] = {k: vectors[k][0] for k in terms}
            else:
                parent['dipole_comparison'][key] = {}
                for k in terms:
                    bound = actual['parents'][gid]['dipole_roundoff_bounds_au'][key][k]
                    check(gid + ':' + key + ':' + k + ':bound_coefficients', abs(bound - math.fsum(abs(a)*rows[cid]['dipole_roundoff_bounds_au'][key] for a, cid in terms[k])) <= 1e-15)
                    error = maxdiff(vectors[k], actual['parents'][gid]['vectors'][key][k])
                    check(gid + ':' + key + ':' + k + ':within_existing_bound', error <= bound)
                    if key == 'nuclear_dipole_au':
                        check(gid + ':' + k + ':original_nuclear_cancellation_bound', max(abs(x) for x in vectors[k]) <= bound + 1e-10)
                    parent['dipole_comparison'][key][k] = {'max_difference_from_printed_au': error, 'existing_bound_au': bound}
                    if key == 'total_dipole_au':
                        maxima['parent_total_dipole_vs_printed_au'] = max(maxima['parent_total_dipole_vs_printed_au'], error)
                        maxima['total_dipole_bound_fraction'] = max(maxima['total_dipole_bound_fraction'], error/bound)
        parents[gid] = parent
    paired = {s: {q: math.fsum(sign*parents[g]['metal_charge_e'][s][q] for sign, g in zip((1,-1,-1,1), GIDS)) for q in ('C','P','B')} for s in ('mulliken','loewdin')}
    for scheme in paired:
        check(scheme + ':paired_replay', max(abs(paired[scheme][q]-actual['paired_metal_population_charge_differences_e'][scheme][q]) for q in paired[scheme]) <= 1e-10)
    result = {'status': 'passed_independent_saved_redistribution_review', 'implementation': pin(__file__), 'configuration': pin(cfg_path),
              'actual_result': cfg['actual_result'], 'root_source': actual['implementation'], 'collection': root_cfg['collection'],
              'source_review': root_cfg['source_review'], 'complete_raw_and_property_rows': len(rows), 'rows': rows,
              'parents': parents, 'paired_metal_charge_e': paired, 'maxima': maxima, 'checks': checks,
              'new_native_calls': 0, 'unique_charge_transfer_identified': False, 'error': None,
              'cost': {'wall_seconds': time.monotonic()-started, 'parent_CPU_seconds': time.process_time()-cpu, 'local_threads': 1}}
    (out/'REVIEW.json').write_text(json.dumps(result, indent=2)+'\n')
    lines = ['# Independent actual redistribution review', '', f'Passed: all 80 raw outputs and 80 full-precision property files; {len(checks)} saved checks; zero new molecular calls.', '',
             'Both raw population schemes replay for every center, including ghosts. The full-precision property dipoles, transformed using COM reconstructed from Bohr coordinates and native masses, agree with all printed-data C/P/B vectors within the original coefficient-propagated bounds.', '',
             '| Parent | Mulliken metal B, e | Loewdin metal B, e | Full-precision total dipole Bz, au |', '|---|---:|---:|---:|']
    for gid, parent in parents.items():
        lines.append(f"| {gid} | {parent['metal_charge_e']['mulliken']['B']:+.9f} | {parent['metal_charge_e']['loewdin']['B']:+.9f} | {parent['vectors']['total_dipole_au']['B'][2]:+.12f} |")
    lines += ['', f"Paired metal B: Mulliken {paired['mulliken']['B']:+.9f} e; Loewdin {paired['loewdin']['B']:+.9f} e. The opposite signs are retained.", '',
              f"Largest parent total-dipole difference from the printed-data analysis: {maxima['parent_total_dipole_vs_printed_au']:.9g} au. Largest native nuclear-moment reconstruction residual: {maxima['nuclear_COM_reconstruction_au']:.9g} au.", '',
              'Positive per-parent metal B occurs in both schemes: the complete shell produces less metal-centered population gain than the additive pair construction under those partitions. This is not a uniquely defined charge-transfer observable or an energy/mechanism attribution. Cross-parent dipole components retain the fixed coordinate orientation; they are not rotation-invariant scalars. Frozen-core ECP, finite-basis, native SCF and physical-state limitations remain.', '',
              f"Saved review cost: {result['cost']['wall_seconds']:.3f} s wall, {result['cost']['parent_CPU_seconds']:.3f} CPU-s, one thread. Original results remain unchanged."]
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'status':result['status'],'complete':len(rows),'checks':len(checks),'maxima':maxima,'cost':result['cost'],'review':pin(out/'REVIEW.json')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    existed = args.output.exists()
    start, cpu = time.monotonic(), time.process_time()
    try:
        main(args.config, args.output)
    except Exception as exc:
        if not existed:
            args.output.mkdir(parents=True, exist_ok=True)
            (args.output/'FAILURE.json').write_text(json.dumps({'status':'failed_independent_saved_review','implementation':pin(__file__),'configuration':pin(args.config),'error':repr(exc),'new_native_calls':0,'wall_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu},indent=2)+'\n')
        raise
