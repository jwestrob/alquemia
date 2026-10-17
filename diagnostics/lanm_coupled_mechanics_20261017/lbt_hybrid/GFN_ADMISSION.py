#!/usr/bin/env python3
"""Dimension/charge generalization of the qualified Mex native admission.

No changed convergence tolerances. Returns native analytic gradients; caller
applies the exact source-dependent cap/physical map. Scientific calls use Slurm.
"""
import hashlib
import json
import math
from pathlib import Path
import re

BOHR_A = 0.529177210903
HA_KCAL = 627.509474
Z = {'H':1, 'C':6, 'N':7, 'O':8, 'La':57, 'Nd':60, 'Dy':66}
TOL = {'charge_closure_e':1e-5, 'xtb_json_stdout_hartree':1e-7,
       'engrad_energy_hartree':1e-8, 'engrad_coordinates_A':1e-6,
       'xtb_json_population_e':1e-7}

def require(ok, message):
    if not ok: raise ValueError(message)

def pin(path):
    p=Path(path).resolve()
    return {'path':str(p), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'], 'Changed pin: '+p['path'])
    return Path(p['path'])

def xyz(p):
    rows=verify(p).read_text().splitlines(); n=int(rows[0]); data=[r.split() for r in rows[2:] if r.strip()]
    require(len(data)==n and all(len(r)==4 for r in data), 'XYZ layout differs')
    sy=[r[0] for r in data]; x=[[float(v) for v in r[1:]] for r in data]
    require(all(math.isfinite(v) for r in x for v in r), 'Nonfinite coordinates')
    return sy,x

def admit(directory, task, execution, *, backend, gradients, np):
    nq = task["atom_count"]; nf = task["pointcharge_count"]
    require(execution['returncode'] == 0 and execution['owned_process_group_empty'], 'Process failure')
    data = backend.parse_backend(directory, 'xtb', task)
    require(data['status'] == 'available', 'SCC/population failure: ' + str(data.get('failure_reason')))
    energy, charges = data['energy_hartree'], data['charges_e']
    payload = json.loads((directory/'xtbout.json').read_text())
    require(payload['method'] == 'GFN2-xTB' and payload['xtb version'].startswith('6.7.1'), 'Actual method/version differs')
    require(payload['number of electrons'] == task['valence_electrons']
            and payload['number of unpaired electrons'] == 0, 'Actual effective electronic state differs')
    require(abs(payload['total energy']-energy) <= TOL['xtb_json_stdout_hartree'], 'JSON/stdout energy mismatch')
    require(len(payload['partial charges']) == nq and all(math.isfinite(q) for q in payload['partial charges'])
            and max(abs(a-b) for a,b in zip(payload['partial charges'],charges)) <= TOL['xtb_json_population_e'],
            'JSON/charge-file populations differ')
    require(abs(math.fsum(charges)-task['charge']) <= TOL['charge_closure_e'], 'QM charge closure')
    text = (directory/'stdout.txt').read_text()
    combined = text + '\n' + (directory/'stderr.txt').read_text()
    require('normal termination of xtb' in combined.lower(), 'Normal termination absent')
    require(not re.search(r'failed to converge|NOT CONVERGED|abnormal termination|numerical gradient', combined, re.I),
            'SCC failure or unsupported numerical gradient')
    for name, expected in ((r'net charge',task['charge']),(r'unpaired electrons',0),(r'# electrons',task['valence_electrons']),
                           (r'electronic temp\.',300),(r':  accuracy',.0001)):
        values = re.findall(name+r'\s+([-+\d.Ee]+)', text)
        require(bool(values) and float(values[-1]) == expected, 'Actual setting differs: '+name)
    embedded = task['kind'] == 'embedded'
    for pattern in (r'max\. iterations\s+500',r'restarted\?\s+false',r'GBSA solvation\s+false',
                    r'PC potential\s+'+('true' if embedded else 'false'),
                    r'convergence criteria satisfied after'):
        require(bool(re.search(pattern,text)), 'Actual SCC setting/evidence absent: '+pattern)
    require('Solvation model:               ALPB' not in text, 'Unexpected solvent')
    if embedded:
        values = re.findall(r'# point charges\s+(\d+)',text)
        totals = re.findall(r'sum of PC\s+([-+\d.Ee]+)',text)
        require(bool(values) and int(values[-1]) == nf and bool(totals)
                and abs(float(totals[-1])-task['field_charge_e']) <= TOL['charge_closure_e'], 'Printed field inventory/charge differs')
    else:
        require(not re.search(r'# point charges',text) and not (directory/'pcgrad').exists(), 'Unexpected isolated field')
    engrad = [directory/n for n in ('core.engrad','xtb-orca.engrad') if (directory/n).is_file()]
    require(len(engrad) == 1, 'Expected one native analytic QM gradient output')
    header = engrad[0].read_text()
    require('Eh/bohr' in header and 'coordinates in Bohr' in header, 'Gradient unit header absent')
    g = gradients.read_engrad(engrad[0]); symbols, positions = xyz(task['xyz'])
    require(g['atom_count'] == nq and abs(g['energy_Ha']-energy) <= TOL['engrad_energy_hartree'],
            'Gradient count/energy mismatch')
    require(np.array_equal(g['atomic_numbers'],[Z[s] for s in symbols])
            and np.allclose(g['coordinates_bohr']*BOHR_A,positions,atol=TOL['engrad_coordinates_A'],rtol=0),
            'Gradient atom order/coordinates differ')
    if embedded:
        records = [line for line in (directory/'pcgrad').read_text().splitlines() if line.strip()]
        # Native writes3f12.8: adjacent full-width negative fields need not have
        # whitespace between them. Parse the actual fixed-width format.
        require(len(records) == nf and all(len(line) == 36 for line in records), 'Native pcgrad layout/count mismatch')
        rows = [[line[k:k+12] for k in (0,12,24)] for line in records]
        pc = np.array([[float(v) for v in r] for r in rows],dtype=float)
        require(pc.shape == (nf,3) and np.isfinite(pc).all(), 'Invalid external-site gradient')
    else:
        pc = np.zeros((nf,3))
    cycles = re.findall(r'convergence criteria satisfied after\s+(\d+) iterations',text)
    require(bool(cycles), 'SCC iteration receipt absent')
    return {'energy_hartree':energy,'charges_e':charges,'charge_closure_e':abs(math.fsum(charges)-task['charge']),
            'solver_observation':data,'scc_iterations':int(cycles[-1]),'xtbout':pin(directory/'xtbout.json'),
            'qm_gradient':{'artifact':pin(engrad[0]),'hartree_per_bohr':g['gradient_Ha_per_bohr'].tolist(),
                           'kcal_mol_A':(g['gradient_Ha_per_bohr']*HA_KCAL/BOHR_A).tolist()},
            'pointcharge_gradient':{'artifact':pin(directory/'pcgrad') if embedded else None,
                 'source':'native analytic pcgrad' if embedded else 'exact zero: isolated Hamiltonian has no external field',
                 'hartree_per_bohr':pc.tolist(),'kcal_mol_A':(pc*HA_KCAL/BOHR_A).tolist(),
                 'rounding_bound_per_component_Eh_bohr':5e-9 if embedded else 0.,'field_mapping':task['field']},
            'quantity':'gradient=dE/dR; force=-gradient; native analytic return, no independent finite-difference qualification',
            'gradient_units':'hartree/bohr and kcal/mol/angstrom',
            'translation_gradient_sum_Eh_bohr':(g['gradient_Ha_per_bohr'].sum(axis=0)+pc.sum(axis=0)).tolist()}
