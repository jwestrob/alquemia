#!/usr/bin/env python3
"""Saved-data comparison on exact matched La/Dy chemistry fixtures."""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import time

EV_KCAL = 23.060547830619
HA_KCAL = 627.509474


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    if pin(record['path']) != record:
        raise ValueError('Changed artifact: ' + record['path'])


def read_xyz(record):
    verify(record)
    lines = Path(record['path']).read_text().splitlines()
    rows = [r.split() for r in lines[2:] if r.strip()]
    if int(lines[0]) != len(rows):
        raise ValueError('Atom count differs')
    return [r[0] for r in rows], [[float(v) for v in r[1:]] for r in rows]


def rms(values):
    return math.sqrt(sum(v*v for v in values) / len(values))


def flat(v):
    return [x for row in v for x in row]


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--panel', required=True)
    a.add_argument('--mace', required=True)
    a.add_argument('--out', required=True)
    a.add_argument('references', nargs='*')
    args = a.parse_args()
    tick = time.monotonic()
    p = json.loads(Path(args.panel).read_text())
    geometries = {g['geometry_id']: g for g in p['geometries']}
    expected = {c['cell_id']: c for c in p['cells']}
    m = json.loads(Path(args.mace).read_text())
    verify(m['config']); verify(m['implementation'])
    if m['status'] != 'complete' or m['native_calls'] != 28 or len(m['rows']) != 28:
        raise ValueError('Incomplete declared MACE panel')
    data = {'MACE': {}}
    for row in m['rows']:
        gid = row['id']; g = geometries[gid]; task = row['task']
        symbols, xyz = read_xyz(task['xyz'])
        if symbols != g['symbols'] or xyz != g['positions_A'] or task['charge'] != g['charge'] or task['physical_multiplicity'] != g['physical_multiplicity']:
            raise ValueError('MACE physical input differs: ' + gid)
        if row['status'] != 'complete' or gid in data['MACE']:
            raise ValueError('Missing/duplicate MACE row')
        data['MACE'][gid] = {'energy': row['energy_eV'] * EV_KCAL,
                            'force': [[x * EV_KCAL for x in r] for r in row['forces_eV_A']]}
    failed = []; ref_pins = []; seen = set()
    for path in args.references:
        r = json.loads(Path(path).read_text()); cid = r['cell_id']
        if cid in seen or r['panel'] != pin(args.panel):
            raise ValueError('Duplicate/reference panel differs')
        seen.add(cid); c = expected[cid]; ref_pins.append(pin(path))
        if r['status'] != 'complete':
            failed.append({'cell': cid, 'status': r['status']}); continue
        if not all(r['checks'].values()):
            raise ValueError('Reference admission contradictory')
        for name in ('output', 'engrad', 'input', 'geometry'):
            verify(r[name])
        if r['geometry'] != c['geometry']['xyz']:
            raise ValueError('Reference geometry differs')
        data.setdefault(c['recipe'], {})[c['geometry_id']] = {
            'energy': r['energy_hartree'] * HA_KCAL, 'force': r['forces_kcal_mol_A']}
    contrasts = []; deformation = []; force_errors = []
    for method, rows in data.items():
        for donor in ('acetate_bi', 'acetate_mono', 'acetamide'):
            for dist in ('minus', 'central', 'plus'):
                ids = [f'{donor}_{dist}_Dy', f'water_{dist}_Dy', f'{donor}_{dist}_La', f'water_{dist}_La']
                if all(k in rows for k in ids):
                    e = [rows[k]['energy'] for k in ids]
                    contrasts.append({'method': method, 'donor': donor, 'distortion': dist,
                                      'Dy_minus_La_kcal_mol': (e[0]-e[1])-(e[2]-e[3])})
        ids = ['aquo9_central_Dy', 'aquo8_central_Dy', 'aquo9_central_La', 'aquo8_central_La']
        if all(k in rows for k in ids):
            e = [rows[k]['energy'] for k in ids]
            contrasts.append({'method': method, 'donor': 'aquo8_to_9', 'distortion': 'central',
                              'Dy_minus_La_kcal_mol': (e[0]-e[1])-(e[2]-e[3])})
        for gid, row in rows.items():
            g = geometries[gid]
            origin = f"{g['case_id']}_central_{g['metal']}"
            if g['distortion'] != 'central' and origin in rows:
                deformation.append({'method': method, 'geometry': gid,
                                    'delta_energy_kcal_mol': row['energy']-rows[origin]['energy']})
    for method, rows in data.items():
        if method == 'MACE': continue
        for gid in rows:
            g = geometries[gid]
            if g['metal'] != 'La': continue
            dy = gid[:-2] + 'Dy'
            if dy not in rows: continue
            ml, md = flat(data['MACE'][gid]['force']), flat(data['MACE'][dy]['force'])
            rl, rd = flat(rows[gid]['force']), flat(rows[dy]['force'])
            if not len(ml) == len(md) == len(rl) == len(rd):
                raise ValueError('Force dimensions differ')
            err = [(d-l)-(b-a) for l,d,a,b in zip(ml,md,rl,rd)]
            force_errors.append({'reference': method, 'case': g['case_id'], 'distortion': g['distortion'],
                                 'paired_component_rmse_kcal_mol_A': rms(err),
                                 'paired_component_max_abs_kcal_mol_A': max(abs(x) for x in err),
                                 'La_component_rmse_kcal_mol_A': rms([a-b for a,b in zip(ml,rl)]),
                                 'Dy_component_rmse_kcal_mol_A': rms([a-b for a,b in zip(md,rd)])})
    result = {'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'panel': pin(args.panel), 'mace': pin(args.mace), 'references': ref_pins,
              'implementation': pin(__file__), 'available': {k: len(v) for k,v in data.items()},
              'failed_references': failed, 'contrasts': contrasts, 'deformation': deformation,
              'force_errors': force_errors, 'elapsed_seconds': time.monotonic()-tick,
              'local_threads': 1, 'new_molecular_calls': 0, 'binding_free_energy_claimed': False}
    out = Path(args.out); out.mkdir(parents=True, exist_ok=False)
    (out/'RESULT.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    lines = ['# Matched fixed-configuration chemistry comparison', '',
             'These are balanced electronic-energy contrasts, not solution binding free energies.', '',
             '| Method | Donor / coordination change | Radial configuration | Dy − La contrast (kcal/mol) |',
             '|---|---|---|---:|']
    lines += [f"| {r['method']} | {r['donor']} | {r['distortion']} | {r['Dy_minus_La_kcal_mol']:.4f} |" for r in contrasts]
    lines += ['', 'Negative contrasts favor Dy relative to water on the specified fixed configurations.',
              'Missing reference cells remain missing. Individual and paired force errors are in RESULT.json.',
              'Reference SCF, electronic-state, basis and spin–orbit limitations remain unresolved.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'available': result['available'], 'contrasts': contrasts, 'force_errors': force_errors}))


if __name__ == '__main__': main()
