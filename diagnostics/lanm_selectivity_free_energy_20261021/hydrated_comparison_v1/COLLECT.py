#!/usr/bin/env python3
"""Read-only molecular inputs; immutable saved arithmetic for hydrated fixtures."""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import time

EV_KCAL = 23.060547830619
HA_KCAL = 627.509474


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    if pin(p['path']) != p:
        raise ValueError('Artifact changed: ' + p['path'])


def load(p):
    return json.loads(Path(p).read_text())


def xyz(record):
    verify(record)
    lines = Path(record['path']).read_text().splitlines()
    atoms = [s.split() for s in lines[2:] if s.strip()]
    if len(atoms) != int(lines[0]) or any(len(s) != 4 for s in atoms):
        raise ValueError('XYZ shape differs')
    return [s[0] for s in atoms], [[float(x) for x in s[1:]] for s in atoms]


def geometry(record, g):
    s, r = xyz(record)
    if s != g['symbols'] or r != g['positions_A']:
        raise ValueError('Geometry differs: ' + g['geometry_id'])


def add(data, method, gid, energy, forces, g):
    if gid in data.setdefault(method, {}):
        raise ValueError('Duplicate molecular row')
    if len(forces) != len(g['symbols']) or any(len(r) != 3 for r in forces):
        raise ValueError('Force shape differs')
    if not all(math.isfinite(x) for x in [energy] + [v for r in forces for v in r]):
        raise ValueError('Nonfinite energy/force')
    data[method][gid] = {'energy': energy, 'forces': forces}


def rms(v):
    return math.sqrt(sum(x*x for x in v) / len(v))


def derivative(row, g):
    return -sum(row['forces'][i][a] * g['radial_direction'][a]
                for i in g['ligand_indices'] for a in range(3))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    t0 = time.monotonic()
    b = Path(args.root).resolve()
    ref = b / 'reference_chemistry_v1'
    mixed_path = ref / 'hydrated_mixed_v1/prepared_v2/PANEL.json'
    native_path = ref / 'hydrated_mixed_v1/bound_v1/NATIVE_PANEL.json'
    orig_path = ref / 'prepared_v2/PANEL.json'
    mixed, native, orig = map(load, [mixed_path, native_path, orig_path])
    if mixed['contact_failures'] or len(mixed['geometries']) != 18 or len(mixed['reused_controls']) != 2:
        raise ValueError('Mixed preparation differs')
    gs = {g['geometry_id']: g for g in mixed['geometries'] + mixed['reused_controls']}
    if len(gs) != 20:
        raise ValueError('Geometry denominator differs')
    for g in gs.values():
        geometry(g['xyz'], g)
    inputs = [pin(p) for p in [mixed_path, native_path, orig_path]]
    data = {}
    for rel, required in [('mace_chemistry_v1/submission_v1/execution/RESULT.json', 28),
                          ('hydrated_chemistry_run_v1/mace/execution/RESULT.json', 18)]:
        path = b / rel; d = load(path); inputs.append(pin(path))
        verify(d['config']); verify(d['implementation'])
        if d['status'] != 'complete' or len(d['rows']) != required or d['native_calls'] != required:
            raise ValueError('MACE panel incomplete')
        if d['energy_units'] != 'eV' or d['force_units'] != 'eV/Angstrom':
            raise ValueError('MACE units differ')
        for r in d['rows']:
            gid = r['id']
            if gid not in gs: continue
            g = gs[gid]; geometry(r['task']['xyz'], g)
            if r['status'] != 'complete' or r['task']['charge'] != g['charge'] or r['task']['physical_multiplicity'] != g['physical_multiplicity']:
                raise ValueError('MACE state differs')
            add(data, 'MACE', gid, r['energy_eV'] * EV_KCAL,
                [[v*EV_KCAL for v in f] for f in r['forces_eV_A']], g)
    path = b / 'explicit_model_v1/hydrated_mixed_run_v1/RESULT.json'
    d = load(path); inputs.append(pin(path)); verify(d['config']); verify(d['implementation'])
    if d['status'] != 'complete' or d['complete_rows'] != 20 or len(d['rows']) != 20:
        raise ValueError('Polarizable panel incomplete')
    for r in d['rows']:
        g = gs[r['id']]; inp = r['input']
        geometry(inp['source_geometry']['xyz'], g)
        for key in ['xyz', 'key', 'parameters']: verify(inp[key])
        for key in ['analyze_output', 'gradient_output']: verify(r['native'][key])
        if r['status'] != 'complete' or inp['charge'] != g['charge']:
            raise ValueError('Polarizable state differs')
        add(data, 'AMOEBA_transfer', r['id'], r['native']['total_energy_kcal_mol'],
            r['native']['force_kcal_mol_A'], g)
        if 'ligand_radial_derivative_kcal_mol_A' in r:
            actual = derivative(data['AMOEBA_transfer'][r['id']], g)
            if abs(actual-r['ligand_radial_derivative_kcal_mol_A']) > 1e-9:
                raise ValueError('Independent radial derivative disagrees')
    missing, failed, admitted = [], [], []
    expected = [(c, native_path, b/'hydrated_chemistry_run_v1/reference/tasks') for c in native['cells']]
    expected += [(c, orig_path, ref/'ecp_remaining_v1/tasks') for c in orig['cells']
                 if c['recipe'] == 'ecp_pbe0' and c['geometry_id'] in ('aquo8_central_La','aquo8_central_Dy')]
    if len(expected) != 8: raise ValueError('Reference inventory differs')
    for c, panelpath, parent in expected:
        path = parent / c['cell_id'] / 'execution/RESULT.json'
        if not path.exists(): missing.append(str(path)); continue
        r = load(path); inputs.append(pin(path))
        if r['panel'] != pin(panelpath) or r['cell_id'] != c['cell_id']:
            raise ValueError('Reference panel/identity differs')
        if r['status'] != 'complete': failed.append({'path':str(path),'status':r['status']}); continue
        if not r['checks'] or not all(r['checks'].values()):
            raise ValueError('Contradictory reference admission')
        for k in ['output','engrad','input','geometry']: verify(r[k])
        if r['geometry'] != c['geometry']['xyz']: raise ValueError('Reference cell geometry differs')
        g = gs[c['geometry_id']]; geometry(r['geometry'],g)
        if c['charge'] != g['charge'] or c['physical_multiplicity'] != g['physical_multiplicity']:
            raise ValueError('Reference chemical state differs')
        add(data,'ecp_pbe0',g['geometry_id'],r['energy_hartree']*HA_KCAL,r['forces_kcal_mol_A'],g)
        admitted.append(c['cell_id'])
    if len(data['MACE']) != 20 or len(data['AMOEBA_transfer']) != 20:
        raise ValueError('Model inventory differs')
    contrasts, radial, errors = [], [], []
    for method, rows in data.items():
        for donor in ['acetate_mono','acetate_bi','acetamide']:
            for dist in ['minus','central','plus']:
                ids = [f'hydrated_{donor}_{dist}_Dy','aquo8_central_Dy',f'hydrated_{donor}_{dist}_La','aquo8_central_La']
                if all(i in rows for i in ids):
                    e = [rows[i]['energy'] for i in ids]
                    contrasts.append({'method':method,'donor':donor,'distortion':dist,'Dy_minus_La_kcal_mol':(e[0]-e[1])-(e[2]-e[3])})
        for gid,row in rows.items():
            if 'ligand_indices' in gs[gid]:
                radial.append({'method':method,'geometry':gid,'ligand_energy_derivative_kcal_mol_A':derivative(row,gs[gid])})
    refs = data.get('ecp_pbe0',{})
    for donor in ['acetate_mono','acetate_bi','acetamide']:
        la,dy=[f'hydrated_{donor}_central_{m}' for m in ['La','Dy']]
        if not all(i in refs for i in [la,dy]): continue
        for method in ['MACE','AMOEBA_transfer']:
            fl,fd,rl,rd=[[v for f in rows[i]['forces'] for v in f] for rows,i in [(data[method],la),(data[method],dy),(refs,la),(refs,dy)]]
            pair=[(b-a)-(d-c) for a,b,c,d in zip(fl,fd,rl,rd)]
            errors.append({'model':method,'reference':'ecp_pbe0','donor':donor,
                           'La_component_rmse_kcal_mol_A':rms([a-b for a,b in zip(fl,rl)]),
                           'Dy_component_rmse_kcal_mol_A':rms([a-b for a,b in zip(fd,rd)]),
                           'paired_component_rmse_kcal_mol_A':rms(pair),
                           'paired_ligand_derivative_error_kcal_mol_A':
                           (derivative(data[method][dy],gs[dy])-derivative(data[method][la],gs[la]))-
                           (derivative(refs[dy],gs[dy])-derivative(refs[la],gs[la]))})
    result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'implementation':pin(__file__),
            'inputs':inputs,'available':{k:len(v) for k,v in data.items()},'missing_references':missing,
            'failed_references':failed,'native_admitted_cells':admitted,'contrasts':contrasts,'radial_derivatives':radial,
            'force_errors':errors,'new_molecular_calls':0,'local_threads':1,'elapsed_seconds':time.monotonic()-t0,
            'reference_accuracy_qualified':False,'binding_free_energy_claimed':False}
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# Hydrated donor competition: saved comparison','','Fixed constructed configurations; these are not binding free energies.','',
           '| Model | Donor | Radial position | Dy minus La (kcal/mol) |','|---|---|---|---:|']
    lines += [f"| {r['method']} | {r['donor']} | {r['distortion']} | {r['Dy_minus_La_kcal_mol']:.4f} |" for r in contrasts]
    lines += ['',f'Missing references: {len(missing)}. Failed references: {len(failed)}.',
              'Native admission does not establish strict SCF residual, electronic-state or reference accuracy qualification.',
              'Full force errors and rigid-ligand energy derivatives are preserved in RESULT.json.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'available':result['available'],'missing':len(missing),'contrasts':contrasts,'force_errors':errors}))


if __name__ == '__main__': main()
