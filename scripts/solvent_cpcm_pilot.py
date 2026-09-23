"""Fixed-state native GFN2 CPCM transfer pilot; reuse the existing ORCA runner."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re
import shutil
import time

from affordable_common import InvalidArtifact, HA_TO_KCAL, read_json, record, verify, write_new, xyz
from affordable_workflow import dry_run, execute as run_existing
from compact_solvation import completed, diagnostics
from mace_hybrid import EV_TO_KCAL
from nikasha_pool import PROTOCOL as SOURCE_PROTOCOL

PROTOCOL = 'native_OMOL_plus_native_GFN2_CPCM_water_transfer_v1'
PILOT = ('1H4I', '4MAE', 'q9z4j7-pqq-la_model', 'a0a3f2yly8-pqq-la_model')
SETTINGS = {'solvent': 'CPCM(Water)', 'UseXTBMixer': True, 'SmearTemp_K': 300,
            'MaxIter': 500, 'parameter_policy': 'installed_native_GFN2_defaults',
            'cavity_policy': 'installed_default_gaussian_vdw_no_CDS_no_DRACO',
            'epsilon': 80.151, 'refrac': 1.3328, 'surface_density_A_minus2': 5.0,
            'expected_radii_A': {'Ca': 2.772, 'La': 2.4, 'C': 2.04, 'N': 1.86, 'O': 1.824, 'H': 1.32},
            'coordinate_tolerance_A': 1e-12, 'charge_closure_e': 5e-4, 'max_abs_charge_e': 4.0}


def same_atoms(a, b):
    return len(a) == len(b) and all(x[0] == y[0] and max(abs(i-j) for i, j in zip(x[1:], y[1:])) <= 1e-12 for x, y in zip(a, b))


def recipe(charge, mult):
    return ('! Native-GFN2-xTB NoAutostart CPCM(Water)\n%maxcore 2000\n'
            '%method\n WriteXTBParam true\n ReadXTBParam false\nend\n'
            '%scf\n SmearTemp 300\n UseXTBMixer true\n MaxIter 500\nend\n'
            f'* xyzfile {charge} {mult} core.xyz\n')


def sources(collections):
    cases, identity = {}, None
    for path in collections:
        col = read_json(path); parent = read_json(verify(col['manifest']))
        if col['protocol_id'] != SOURCE_PROTOCOL: raise InvalidArtifact('source is not the actual original common pool')
        sig = {k: parent[k] for k in ('model', 'software', 'orca')}
        if identity is not None and identity != sig: raise InvalidArtifact('source method/software differs')
        identity = sig
        for c in col['cases']:
            if c['case_id'] in cases: raise InvalidArtifact('duplicate source case')
            cases[c['case_id']] = (c, record(path))
    return cases, identity


def original_endpoint(case, metal, identity):
    cell = case['matrix'][metal]['origin']
    if cell['status'] != 'complete': raise InvalidArtifact('actual original endpoint unavailable')
    coords = xyz(verify(cell['xyz'])); native = read_json(verify(cell['MACE']))
    request = read_json(verify(native['request']))
    if (native['status'] != 'complete' or native['model'] != identity['model'] or
            not same_atoms(coords, xyz(verify(request['xyz']))) or
            coords[0][0] != metal or request['multiplicity'] != 1 or
            native['state']['charge'] != request['charge'] or
            native['state']['spin_multiplicity'] != request['multiplicity']):
        raise InvalidArtifact('native MACE reuse state differs')
    low = {}
    for medium in ('vacuum', 'alpb'):
        endpoint = cell['low'][medium]; parent = read_json(verify(endpoint['manifest']))
        task = next(t for t in parent['all_tasks'] if t['task_id'] == endpoint['task_id'])
        if (parent['orca'] != identity['orca'] or task['solver'] != 'native' or
                task['charge'] != request['charge'] or task['multiplicity'] != request['multiplicity'] or
                not same_atoms(coords, xyz(verify(task['xyz'])))):
            raise InvalidArtifact('native GFN2 reuse state differs')
        audit = diagnostics(endpoint, task)
        if ('INFO: Using special xTB SCF mixer' not in verify(endpoint['output']).read_text() or
                audit['charge_sanity_status'] != 'pass'):
            raise InvalidArtifact('archived native solver/charge check differs')
        low[medium] = {'endpoint': endpoint, 'task': task, 'audit': audit}
    if read_json(verify(low['vacuum']['audit']['parameter_export'])) != read_json(verify(low['alpb']['audit']['parameter_export'])):
        raise InvalidArtifact('old vacuum and ALPB parameter sets differ')
    return {'xyz': cell['xyz'], 'charge': request['charge'], 'multiplicity': request['multiplicity'],
            'MACE': cell['MACE'], 'MACE_energy_eV': native['energy_eV'], 'GFN2': low,
            'source_preparation': low['vacuum']['task']['source_preparation']}


def prepare(collections, agreement, output):
    start = time.monotonic(); source, identity = sources(collections)
    if set(PILOT) - set(source): raise InvalidArtifact('declared pilot member absent')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    cases, tasks = [], []
    for cid in PILOT:
        case, pin = source[cid]; endpoints = {z: original_endpoint(case, z, identity) for z in ('Ca', 'La')}
        ca, la = (xyz(verify(endpoints[z]['xyz'])) for z in ('Ca', 'La'))
        if not same_atoms([('M', *a[1:]) if i == 0 else a for i, a in enumerate(ca)],
                          [('M', *a[1:]) if i == 0 else a for i, a in enumerate(la)]):
            raise InvalidArtifact('Ca/La nuclear geometry differs')
        if endpoints['La']['charge'] != endpoints['Ca']['charge']+1: raise InvalidArtifact('paired charge policy differs')
        cases.append({'case_id': cid, 'source_collection': pin, 'old_result': case['old_result'], 'endpoints': endpoints})
        for metal, e in endpoints.items():
            tid = f'{cid}__{metal}__cpcm_native'; td = out/'tasks'/tid; td.mkdir(parents=True)
            xp = td/'core.xyz'; shutil.copyfile(verify(e['xyz']), xp)
            ip = td/'endpoint.inp'; ip.write_text(recipe(e['charge'], e['multiplicity']))
            tasks.append({'task_id': tid, 'case_id': cid, 'case': cid, 'metal': metal, 'medium': 'cpcm',
                          'representation': 'context', 'charge': e['charge'], 'multiplicity': e['multiplicity'],
                          'xyz': record(xp), 'input': record(ip), 'output_path': str(td/'endpoint.out')})
    impl = out/'implementation'; impl.mkdir(); pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        target = impl/p.name; shutil.copyfile(p, target); pins[p.name] = record(target)
    manifest = {'protocol_id': PROTOCOL, 'settings': SETTINGS, **identity, 'agreement': record(agreement),
                'source_collections': [record(p) for p in collections], 'cases': cases, 'tasks': tasks,
                'implementation': pins, 'execution_resources': {'mpi_ranks': 8, 'concurrent_tasks': 8},
                'execution_policy': {'task_runner': pins['run_orca_task_manifest.py'],
                                     'runtime_renderer': pins['render_orca_runtime_input.py']},
                'reference': None, 'new_DFT_calls': 0, 'new_MACE_calls': 0, 'new_GFN2_calls': 8,
                'source_preparation_wall_seconds': time.monotonic()-start, 'production_changed': False}
    path = out/'manifest.json'; write_new(path, manifest)
    result = validate(path); write_new(out/'PREFLIGHT.json', result); return result


def validate(manifest):
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL or m['settings'] != SETTINGS or m['new_GFN2_calls'] != 8:
        raise InvalidArtifact('frozen solvent protocol differs')
    if len(m['cases']) != 4 or {c['case_id'] for c in m['cases']} != set(PILOT): raise InvalidArtifact('pilot population changed')
    if len(m['tasks']) != 8 or {(t['case_id'], t['metal']) for t in m['tasks']} != {(c,z) for c in PILOT for z in ('Ca','La')}:
        raise InvalidArtifact('finite task set changed')
    verify(m['agreement']); verify(m['orca'])
    for p in m['source_collections']: verify(p)
    for p in m['implementation'].values(): verify(p)
    cases = {c['case_id']: c for c in m['cases']}
    for t in m['tasks']:
        e = cases[t['case_id']]['endpoints'][t['metal']]
        if (t['charge'], t['multiplicity']) != (e['charge'], e['multiplicity']): raise InvalidArtifact('task state differs')
        if verify(t['xyz']).read_bytes() != verify(e['xyz']).read_bytes(): raise InvalidArtifact('source coordinates changed')
        if verify(t['input']).read_text() != recipe(t['charge'], t['multiplicity']): raise InvalidArtifact('CPCM input recipe differs')
        if not Path(t['output_path']).is_relative_to(Path(manifest).resolve().parent/'tasks'):
            raise InvalidArtifact('task output escapes own workspace')
    dry_run(manifest)
    return {'status': 'validated', 'tasks': 8, 'cases': 4, 'manifest': record(manifest)}


def cpcm_audit(endpoint, task, source):
    audit = diagnostics(endpoint, task); text = verify(endpoint['output']).read_text()
    if audit['charge_sanity_status'] != 'pass': raise InvalidArtifact('CPCM atomic charges outside declared range')
    if 'INFO: Using special xTB SCF mixer' not in text:
        raise InvalidArtifact('effective mixer changed: matched ordinary-SCF vacuum qualification required')
    if read_json(verify(audit['parameter_export'])) != read_json(verify(source['GFN2']['vacuum']['audit']['parameter_export'])):
        raise InvalidArtifact('CPCM/native-vacuum parameters differ')
    if 'CPCM SOLVATION MODEL' not in text or 'GAUSSIAN VDW' not in text:
        raise InvalidArtifact('actual CPCM cavity not confirmed')
    block = text.split('CPCM SOLVATION MODEL', 1)[1]
    scalar = lambda name: float(re.findall(r'^\s*'+re.escape(name)+r'\s+\.\.\.\s+([-+0-9.eE]+)', block, re.M)[0])
    eps, refrac = scalar('Epsilon'), scalar('Refrac')
    if abs(eps-SETTINGS['epsilon']) > 1e-4 or abs(refrac-SETTINGS['refrac']) > 1e-4:
        raise InvalidArtifact('installed CPCM water constants differ from frozen default')
    radii = {z: float(r) for z, r in re.findall(r'Radius for (\w+)\s+used is\s+[-+0-9.eE]+ Bohr \(=\s+([-+0-9.eE]+) Ang\.\)', block)}
    for z in {a[0] for a in xyz(verify(task['xyz']))}:
        if z not in radii: raise InvalidArtifact('CPCM element radius unprinted')
        if z in SETTINGS['expected_radii_A'] and abs(radii[z]-SETTINGS['expected_radii_A'][z]) > 1e-4:
            raise InvalidArtifact('CPCM radius differs from installed default policy')
    hits = re.findall(r'^\s*CPCM Dielectric\s*:\s*([-+0-9.eE]+) Eh', text, re.M)
    if not hits or not math.isfinite(float(hits[-1])) or abs(float(hits[-1])) < 1e-8:
        raise InvalidArtifact('nonzero CPCM reaction-field energy not demonstrated')
    return {**audit, 'effective_native_mixer': True, 'CPCM_dielectric_hartree': float(hits[-1]),
            'CPCM_epsilon': eps, 'CPCM_refrac': refrac, 'CPCM_radii_A': radii,
            'CPCM_setup_text': block.split('Overall time for CPCM initialization',1)[0],
            'warnings': [s for s in text.splitlines() if 'WARNING' in s.upper() or 'IGNOR' in s.upper()]}


def collect(manifest, output):
    validate(manifest); m = read_json(manifest); cases = {c['case_id']: c for c in m['cases']}; endpoints = []
    for t in m['tasks']:
        endpoint = completed(manifest, t['task_id']); source = cases[t['case_id']]['endpoints'][t['metal']]
        row = {'task_id': t['task_id'], 'case_id': t['case_id'], 'metal': t['metal'], 'status': 'unavailable',
               'energy_hartree': None, 'raw_converged_endpoint': endpoint, 'source': source}
        if endpoint:
            try: row.update(status='complete', **endpoint, audit=cpcm_audit(endpoint,t,source))
            except (InvalidArtifact, IndexError, KeyError, ValueError) as exc: row['reason'] = str(exc)
        else:
            op = Path(t['output_path']); row['artifacts'] = [record(p) for p in (op,Path(str(op)+'.execution.json')) if p.exists()]
            row['reason'] = 'failed, nonconverged or unrun'
        endpoints.append(row)
    rows = []
    for cid, c in cases.items():
        es = {z: next(e for e in endpoints if (e['case_id'],e['metal']) == (cid,z)) for z in ('Ca','La')}
        components = {}
        for z, e in es.items():
            s = e['source']; vac = s['GFN2']['vacuum']['endpoint']['energy_hartree']; alpb = s['GFN2']['alpb']['endpoint']['energy_hartree']
            solvent = (e['energy_hartree']-vac)*HA_TO_KCAL if e['status'] == 'complete' else None
            components[z] = {'MACE_vacuum_eV': s['MACE_energy_eV'], 'GFN2_vacuum_hartree': vac,
                             'GFN2_ALPB_hartree': alpb, 'GFN2_CPCM_hartree': e['energy_hartree'],
                             'CPCM_transfer_kcal_mol': solvent,
                             'composite_energy_kcal_mol': s['MACE_energy_eV']*EV_TO_KCAL+solvent if solvent is not None else None}
        good = all(e['status']=='complete' for e in es.values()); old = c['old_result']['R0']['composite_R_model_kcal_mol']
        score = components['Ca']['composite_energy_kcal_mol']-components['La']['composite_energy_kcal_mol'] if good else None
        rows.append({'case_id': cid, 'expected_class': c['old_result']['expected_class'], 'label_scope': c['old_result']['label_scope'],
                     'status': 'available' if good else 'unavailable', 'components': components,
                     'released_ALPB_R': old, 'CPCM_R': score, 'CPCM_minus_ALPB_R': score-old if good else None,
                     'new_reference': None, 'new_decision': None, 'affinity_claim': False})
    result = {'protocol_id': PROTOCOL, 'manifest': record(manifest), 'endpoints': endpoints, 'rows': rows,
              'complete_endpoints': sum(e['status']=='complete' for e in endpoints), 'endpoint_denominator': 8,
              'available_pairs': sum(r['status']=='available' for r in rows), 'case_denominator': 4,
              'new_DFT_MACE_calls': 0, 'production_changed': False, 'reference': None}
    write_new(output,result); return {k:v for k,v in result.items() if k not in ('rows','endpoints')}


def execute(manifest):
    validate(manifest); return run_existing(manifest)


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='operation', required=True)
    a=sub.add_parser('prepare');a.add_argument('--collections',nargs='+',required=True);a.add_argument('--agreement',required=True);a.add_argument('--output',required=True)
    for name in ('validate','execute','collect'):
        a=sub.add_parser(name);a.add_argument('--manifest',required=True)
        if name=='collect':a.add_argument('--output',required=True)
    args=vars(p.parse_args());op=args.pop('operation');print(json.dumps(globals()[op](**args)))


if __name__ == '__main__': main()
