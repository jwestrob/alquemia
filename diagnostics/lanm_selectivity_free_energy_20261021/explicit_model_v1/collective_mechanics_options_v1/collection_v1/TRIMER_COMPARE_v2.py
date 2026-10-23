#!/usr/bin/env python3
"""Read-only complete seven-term energy and electric-moment comparison."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import time

HA = 627.509474
BOHR = 0.529177210544
DEBYE = 4.803206802
CONVERSION = BOHR * DEBYE
NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path']) == record, 'Changed artifact: ' + record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def module(record, name):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def difference(a, b):
    require(len(a) == len(b), 'Vector lengths differ')
    return max(abs(x-y) for x, y in zip(a, b))


def vector(text, label):
    hits = re.findall(r'^\s*' + re.escape(label) + r'\s*:\s*(' + NUMBER +
                      r')\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s*$', text, re.M)
    require(len(hits) == 1, 'Missing/duplicate vector: ' + label)
    values = list(map(float, hits[0]))
    require(all(math.isfinite(x) for x in values), 'Nonfinite moment')
    return values


def quantum_moment(output_pin, expected_charge, helper):
    """Same raw/property and COM bounds as the admitted independent 80-row reader."""
    output = verify(output_pin)
    text = output.read_text()
    prop_path = output.with_name('endpoint.property.txt')
    prop = prop_path.read_text()
    require('ORCA 6.1.1' in prop, 'Unexpected property version')
    blocks = text.split('CARTESIAN COORDINATES (A.U.)')
    require(len(blocks) == 2, 'Missing/duplicate AU coordinates')
    atoms = []
    for line in blocks[1].splitlines():
        match = re.fullmatch(r'\s*(\d+)\s+([A-Za-z]+)\s+(' + NUMBER +
            r')\*?\s+\d+\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s+(' +
            NUMBER + r')\s+(' + NUMBER + r')\s*', line)
        if match:
            atoms.append((int(match[1]), match[2], float(match[3]), float(match[4]),
                          [float(match[j]) for j in (5, 6, 7)]))
        elif atoms:
            break
    require(atoms and [x[0] for x in atoms] == list(range(len(atoms))), 'AU inventory')
    geometry = helper.section(prop, 'Geometry')
    require('&Units "Bohr"' in geometry, 'Property coordinate units differ')
    positions = []
    for line in geometry.split('&Units "Bohr"]', 1)[1].splitlines():
        match = re.fullmatch(r'\s*([A-Za-z]+)(:?)\s+(' + NUMBER + r')\s+(' +
                             NUMBER + r')\s+(' + NUMBER + r')\s*', line)
        if match:
            require(len(positions) < len(atoms), 'Too many property centers')
            atom = atoms[len(positions)]
            require(match[1] == atom[1] and bool(match[2]) == (atom[2] == 0),
                    'Property label/ghost differs')
            positions.append([float(match[j]) for j in (3, 4, 5)])
    require(len(positions) == len(atoms), 'Incomplete property coordinates')
    require(max(difference(x[4], r) for x, r in zip(atoms, positions)) <= 0.5e-6,
            'Property geometry differs from raw')
    origins = re.findall(r'^Position of electric origin\s*\.\.\.\s*(' + NUMBER +
        r')\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s*$', text, re.M)
    electrons = re.findall(r'^\s*Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)\s*$', text, re.M)
    require(len(origins) == len(electrons) == 1, 'Origin/electron denominator differs')
    raw_origin = list(map(float, origins[0]))
    ne = int(electrons[0]); za = math.fsum(x[2] for x in atoms)
    require(za - ne == expected_charge and za > 0, 'Effective charge mismatch')
    require(all((x[2] == 0) == (x[3] == 0) for x in atoms), 'Ghost mass/charge mismatch')
    mass = math.fsum(x[3] for x in atoms)
    mass_error = sum(5e-4 for x in atoms if x[3] > 0)
    require(mass > mass_error, 'Invalid mass denominator')
    com = [math.fsum(x[3]*r[j] for x, r in zip(atoms, positions))/mass for j in range(3)]
    bounds = [1e-6 + math.fsum(5e-4*abs(r[j]-com[j]) for x, r in zip(atoms, positions)
                             if x[3] > 0)/(mass-mass_error) for j in range(3)]
    require(all(abs(a-b) <= t for a, b, t in zip(com, raw_origin, bounds)), 'Native COM differs')
    block = helper.section(prop, 'SCF_Dipole_Moment')
    nuclear = helper.array(block, 'dipoleNucContrib', 3)
    nuclear_zero = [math.fsum(x[2]*r[j] for x, r in zip(atoms, positions)) for j in range(3)]
    precise_origin = [(nuclear_zero[j]-nuclear[j])/za for j in range(3)]
    require(difference(precise_origin, raw_origin) <= 0.5e-6, 'Precise origin differs')
    require(all(abs(a-b) <= t for a, b, t in zip(precise_origin, com, bounds)), 'Precise COM differs')
    row = {'output': output_pin, 'property_file': pin(prop_path), 'charge': expected_charge,
           'native_effective_nuclear_charge': za, 'native_electrons': ne,
           'precise_origin_Bohr': precise_origin, 'printed_origin_Bohr': raw_origin,
           'labels': [x[1] for x in atoms], 'ghost_mask': [x[2] == 0 for x in atoms],
           'positions_Bohr': positions, 'raw_common_origin': {}, 'print_bounds_au': {}}
    for name, field, label, charge in (
        ('total_dipole_au', 'dipoleTotal', 'Total Dipole Moment', expected_charge),
        ('electronic_dipole_au', 'dipoleElecContrib', 'Electronic contribution', -ne),
        ('nuclear_dipole_au', 'dipoleNucContrib', 'Nuclear contribution', za)):
        full = helper.array(block, field, 3); raw = vector(text, label)
        require(difference(full, raw) <= 0.5e-9, 'Property dipole does not round to raw')
        row[name] = [full[j]+charge*precise_origin[j] for j in range(3)]
        row['raw_common_origin'][name] = [raw[j]+charge*raw_origin[j] for j in range(3)]
        bound = 5e-10 + abs(charge)*5e-7
        row['print_bounds_au'][name] = bound
        require(difference(row[name], row['raw_common_origin'][name]) <= bound,
                'Common-origin moment outside original print bound')
    require(difference(row['total_dipole_au'], [a+b for a,b in
            zip(row['electronic_dipole_au'],row['nuclear_dipole_au'])]) <= 5e-9,
            'Dipole component closure failed')
    return row


def combine(rows, terms, field):
    require(set(terms) <= set(rows), 'Missing cycle term')
    return [math.fsum(coef*rows[cid][field][j] for cid,coef in terms.items()) for j in range(3)]


def compare(cfg):
    require(cfg['implementation'] == pin(__file__), 'Collector source changed')
    verify(cfg['declaration']); m = read(cfg['manifest'])
    executor = module(cfg['executor'], 'trimer_executor_for_collection')
    qreview = read(m['sources']['quantum_moments_review'])
    qhelper = module(qreview['implementation'], 'admitted_property_helpers')
    qdefs = {c['cell_id']:c for c in m['reference_cells']}
    require(len(qdefs) == len(m['model_rows']) == 16 and len(m['cycles']) == 8,
            'Trimer denominator differs')
    require(set(cfg['reference_results']) == set(qdefs), 'Incomplete reference inventory')
    energies = {}; qrows = {}; native = {}
    for cid, result_pin in cfg['reference_results'].items():
        result = read(result_pin); rcfg = read(result['config'])
        require(rcfg['implementation'] == cfg['executor'] and result['manifest'] == cfg['manifest'],
                'Reference runtime source/manifest differs')
        executor.admit_reference(rcfg, m, verify(result_pin))
        require(result['cell'] == qdefs[cid], 'Reference cell differs')
        native[cid] = result['parsed']; energies[cid] = result['parsed']['electronic_energy_hartree']
        qrows[cid] = quantum_moment(result['output'], qdefs[cid]['charge'], qhelper)
        require(qrows[cid]['labels'] == qdefs[cid]['symbols'] and
                qrows[cid]['ghost_mask'] == qdefs[cid]['ghost_mask'], 'Reference moment atom inventory differs')
        require(max(difference([v*BOHR for v in r], p) for r,p in
                    zip(qrows[cid]['positions_Bohr'],qdefs[cid]['positions_A'])) <= 1e-6,
                'Moment geometry differs from declared geometry')
    for cid, old in m['reused_reference'].items():
        require(old['record']['output'] == old['moment_row']['output'], 'Reused energy/moment join differs')
        qrows[cid] = quantum_moment(old['record']['output'], old['moment_row']['expected_total_charge_e'], qhelper)
        for field in ('total_dipole_au','nuclear_dipole_au','electronic_dipole_au'):
            require(difference(qrows[cid][field],qreview['rows'][cid][field]) <= 1e-12,
                    'Reused full-precision authority differs')
        energies[cid] = old['record']['electronic_energy_hartree']
    model = read(cfg['model_result']); mcfg = read(model['config'])
    require(mcfg['implementation'] == cfg['executor'] and model['manifest'] == cfg['manifest'],
            'Model runtime source/manifest differs')
    executor.admit_model(mcfg,m,model)
    model_rows = {cid:r['parsed_moments'] for cid,r in model['rows'].items()}
    model_rows.update({cid:r['record']['moments'] for cid,r in m['reused_model'].items()})
    cycles = {}
    for cid, cycle in m['cycles'].items():
        qt = cycle['reference_terms']; mt = cycle['model_terms']
        require(len(qt) == len(mt) == 7, 'Seven-term denominator differs')
        require(math.fsum(c*qrows[k]['charge'] for k,c in qt.items()) == 0, 'Nonneutral reference cycle')
        qe = HA*math.fsum(c*energies[k] for k,c in qt.items())
        me = math.fsum(c*model_rows[k]['energy_kcal_mol'] for k,c in mt.items())
        q = [CONVERSION*x for x in combine(qrows,qt,'total_dipole_au')]
        x = combine(model_rows,mt,'induced_sum_Debye')
        total = combine(model_rows,mt,'total_lab_Debye')
        qbound = CONVERSION*math.fsum(abs(c)*qrows[k]['print_bounds_au']['total_dipole_au'] for k,c in qt.items())
        mbound = cycle['induced_print_bound_Debye_per_component']
        require(difference(total,x) <= cycle['M_print_bound_Debye_per_component']+mbound+2e-9*sum(abs(c) for c in mt.values()),
                'Permanent/origin cancellation fails model total-M check')
        nuclear = combine(qrows,qt,'nuclear_dipole_au')
        nb = math.fsum(abs(c)*qrows[k]['print_bounds_au']['nuclear_dipole_au'] for k,c in qt.items())
        require(max(map(abs,nuclear)) <= nb, 'Reference nuclear cycle failed')
        names = set().union(*(model_rows[k]['components'] for k in mt))
        components = {name:math.fsum(c*model_rows[k]['components'].get(name,{}).get('energy_kcal_mol',0)
                                   for k,c in mt.items()) for name in sorted(names)}
        require(abs(me-math.fsum(components.values())) <= 7e-6, 'Cycle energy component closure failed')
        cycles[cid] = {'selection':cycle['selection'],'reference_terms':qt,'model_terms':mt,
            'reference_electronic_E3_kcal_mol':qe,'model_E3_kcal_mol':me,'model_minus_reference_E3_kcal_mol':me-qe,
            'model_E3_components_kcal_mol':components,'reference_mu3_Debye':q,'model_mu3_Debye':x,
            'model_total_M_mu3_Debye':total,'model_minus_reference_mu3_Debye':[a-b for a,b in zip(x,q)],
            'combined_print_bound_Debye_per_component':qbound+mbound,
            'reference_nuclear_cycle_au':nuclear,'reference_print_bound_Debye_per_component':qbound,
            'model_print_bound_Debye_per_component':mbound}
    return {'status':'complete_saved_trimer_comparison','inputs':cfg,'cycles':cycles,
            'quantum_moment_rows':qrows,'reference_native_metadata':native,'new_native_calls':0,
            'conversion_Debye_per_au':CONVERSION,'hartree_to_kcal_mol':HA,
            'physical_acceptance_threshold':None,'fitting_performed':False,
            'limitations':['Electronic E3 excludes D4.','Eight sampled cycles cannot reconstruct full B.',
                           'Print bounds exclude reference and model uncertainty.',
                           'Quantum E3 is not uniquely polarization or charge transfer.']}


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--config',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();cpu=time.process_time()
    try:
        result=compare(json.loads(Path(args.config).read_text()));filename='RESULT.json'
    except Exception as exc:
        result={'status':'failed_saved_trimer_comparison','error':repr(exc),'implementation':pin(__file__),'new_native_calls':0};filename='FAILURE.json'
    result['cost']={'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'local_threads':1}
    (out/filename).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'output':str(out)}));raise SystemExit(0 if filename=='RESULT.json' else 1)
