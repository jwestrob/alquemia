#!/usr/bin/env python3
"""Isolated native fixed-charge energy/force test on actual LBT3 artifacts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
import traceback

TERMS = ('Gelec', 'Ghb', 'Gsasa', 'Gshift')
BOHR_A = 0.529177210903


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed input: '+record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, data):
    with Path(path).open('x') as handle:
        json.dump(data, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def xyz(path):
    lines = Path(path).read_text().splitlines()
    require(int(lines[0]) == 256 and len(lines) == 258, '256-atom fixture required')
    atoms = [line.split() for line in lines[2:]]
    return [a[0] for a in atoms], np.array([[float(z) for z in a[1:]] for a in atoms])


def evaluate(binary, task, positions, charges, out, env, calls):
    out.mkdir()
    (out/'input.xyz').write_text('256\nModeled derivative probe; fixed archived charges\n'+''.join(
        '%s %.17g %.17g %.17g\n' % (el, *v) for el, v in zip(task['symbols'], positions)))
    (out/'charges').write_text(''.join('%.17g\n' % q for q in charges))
    command = [str(binary), str(out/'input.xyz'), str(out/'charges')]
    write(out/'START.json', {'command': command, 'source': task['source'], 'positions': pin(out/'input.xyz'),
                           'charges': pin(out/'charges'), 'spin_role': 'No electronic calculation'})
    started = time.monotonic()
    with (out/'stdout.txt').open('x') as stdout, (out/'stderr.txt').open('x') as stderr:
        p = subprocess.run(command, cwd=out, env=env, stdout=stdout, stderr=stderr)
    execution = {'command': command, 'returncode': p.returncode, 'wall_seconds': time.monotonic()-started}
    write(out/'EXECUTION.json', execution)
    calls.append({'directory': str(out), 'execution': execution})
    require(p.returncode == 0, 'Native component call failed: '+str(out))
    text = (out/'stdout.txt').read_text()
    def values(label, n):
        found = re.findall(r'^'+label+r'\s+(.+)$', text, re.M)
        require(len(found) == 1, 'Missing or repeated '+label)
        v = [float(x) for x in found[0].split()]
        require(len(v) == n and np.isfinite(v).all(), 'Bad '+label)
        return v
    parts = dict(zip(TERMS, values('COMPONENTS', 4)))
    bounds = dict(zip(('Gelec', 'Ghb'), values('ROUNDING_BOUNDS', 2)))
    charge = values('CHARGE_SUM', 1)[0]
    gradient = np.loadtxt(out/'gradient.dat')
    require(gradient.shape == (256, 3) and np.isfinite(gradient).all(), 'Invalid native gradient')
    require(abs(charge+2.) <= 1e-5, 'Charge closure failed')
    result = {'components_hartree': parts, 'total_hartree': sum(parts.values()),
              'rounding_bounds_hartree': bounds, 'charge_e': charge,
              'gradient_hartree_per_bohr': gradient.tolist(), 'gradient': pin(out/'gradient.dat'),
              'force_convention': 'force = minus gradient', 'execution': execution}
    write(out/'RESULT.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Actual Slurm allocation required')
    cfg = json.loads(args.config.read_text())
    require(pin(__file__) == cfg['implementation'], 'Changed implementation')
    verify(cfg['plan'])
    global np
    import numpy as np
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    result = {'protocol_id': 'native_frozen_charge_ALPB_force_lbt3_v1', 'status': 'unavailable',
              'config': pin(args.config), 'plan': cfg['plan'], 'slurm_job_id': os.environ['SLURM_JOB_ID'],
              'MACE_calls': 0, 'SCC_calls': 0, 'quantum_calls': 0, 'native_calls': [], 'checks': []}
    write(out/'START.json', result)
    tick = time.monotonic()
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               BLIS_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', LD_LIBRARY_PATH=cfg['library_directory'])
    try:
        for record in cfg['build_inputs']:
            verify(record)
        previous = read(cfg['previous_build'])
        require(previous['status'] == 'built', 'Prior native build unavailable')
        for record in previous['linked_libraries']:
            verify(record)
        source = verify(cfg['fortran'])
        binary = out/'frozen_solvent_force'
        command = cfg['compiler_command']+[str(source), '-o', str(binary)]+cfg['link_arguments']
        with (out/'compile.out').open('x') as stdout, (out/'compile.err').open('x') as stderr:
            p = subprocess.run(command, cwd=out, env=env, stdout=stdout, stderr=stderr)
        result['build'] = {'command': command, 'returncode': p.returncode, 'fortran': cfg['fortran']}
        require(p.returncode == 0, 'Wrapper compilation failed')
        result['build']['binary'] = pin(binary)
        parent = read(cfg['solvent_result'])
        require(parent['status'] == 'complete' and len(parent['rows']) == 8, 'Eight accepted original endpoints required')
        tasks = {}
        for row in parent['rows']:
            t = row['task']
            require(row['status'] == 'accepted' and t['charge'] == -2 and t['atom_count'] == 256,
                    'Invalid real fixture state')
            symbols, positions = xyz(verify(t['xyz']))
            q = np.loadtxt(verify(row['artifacts']['charges']))
            atoms = read(t['physical_source'])
            require(len(atoms) == 256 and symbols[255] == t['metal'] and q.shape == (256,)
                    and np.isfinite(q).all() and abs(q.sum()+2.) <= 1e-5
                    and np.array_equal(q, np.asarray(row['charges_e']))
                    and np.allclose(positions, np.array([a['xyz_A'] for a in atoms]), rtol=0, atol=1e-12),
                    'Actual charge/coordinate identity differs')
            key = (t['source_id'], t['metal'], t['medium'])
            require(key not in tasks, 'Duplicate fixture')
            tasks[key] = {'source': t, 'symbols': symbols, 'x': positions, 'q': q, 'atoms': atoms,
                          'stdout': verify(row['artifacts']['stdout.txt'])}
        require(set(tasks) == {(s,m,v) for s in ('7CCO','7CCN') for m in ('La','Dy')
                               for v in ('vacuum','ALPB(Water)')}, 'Fixture denominator changed')
        for (s,m,v), task in tasks.items():
            if v != 'ALPB(Water)':
                continue
            actual = evaluate(binary,task,task['x'],task['q'],out/(s+'_'+m+'_water_replay'),env,result['native_calls'])
            oldtext = task['stdout'].read_text()
            for term in TERMS:
                matches = re.findall(r'->\s*'+term+r'\s+([-+0-9.Ee]+)\s+Eh', oldtext)
                require(len(matches) == 1, 'Archived native component missing')
                expected = float(matches[0]); bound = actual['rounding_bounds_hartree'].get(term,0.)+1e-9
                residual = abs(actual['components_hartree'][term]-expected)
                result['checks'].append({'kind':'native_energy_replay','source':s,'metal':m,'component':term,
                    'residual_hartree':residual,'bound_hartree':bound,'pass':bound<=1e-6 and residual<=bound})
        require(all(c['pass'] for c in result['checks']), 'Native component replay gate failed')
        bases = {}
        for (s,m,v), task in tasks.items():
            if v == 'vacuum':
                bases[(s,m)] = evaluate(binary,task,task['x'],task['q'],out/(s+'_'+m+'_vacuum_charges'),env,result['native_calls'])
        for m in ('La','Dy'):
            task=tasks[('7CCO',m,'vacuum')]; base=bases[('7CCO',m)]
            donor = next(i for i,a in enumerate(task['atoms']) if a['id']=='A/5//OD1')
            for name,index in [('metal_x',255),('Asn5_OD1_x',donor)]:
                analytic = base['gradient_hartree_per_bohr'][index][0]/BOHR_A
                for step in (1e-4,5e-5):
                    energies=[]
                    for sign in (-1,1):
                        x=task['x'].copy();x[index,0]+=sign*step
                        call=evaluate(binary,task,x,task['q'],out/('%s_%s_%g_%+d'%(m,name,step,sign)),env,result['native_calls'])
                        energies.append(call['total_hartree'])
                    finite=(energies[1]-energies[0])/(2*step); residual=abs(finite-analytic)
                    tolerance=1e-5+1e-4*abs(analytic)
                    result['checks'].append({'kind':'directional_derivative','metal':m,'direction':name,'step_A':step,
                        'analytic_Eh_A':analytic,'finite_difference_Eh_A':finite,'residual_Eh_A':residual,
                        'tolerance_Eh_A':tolerance,'pass':residual<=tolerance})
        task=tasks[('7CCO','La','vacuum')];base=bases[('7CCO','La')]
        for name,shift in [('repeat',(0.,0.,0.)),('joint_translation',(1.25,-.75,.5))]:
            call=evaluate(binary,task,task['x']+np.array(shift),task['q'],out/name,env,result['native_calls'])
            energy_error=abs(call['total_hartree']-base['total_hartree'])
            gradient_error=float(np.max(np.abs(np.array(call['gradient_hartree_per_bohr'])-base['gradient_hartree_per_bohr'])))
            result['checks'].append({'kind':name,'energy_error_Eh':energy_error,'max_gradient_error_Eh_bohr':gradient_error,
                                    'pass':energy_error<=1e-10 and gradient_error<=1e-9})
        require(len(result['native_calls']) == 26, 'Declared component-call count differs')
        result['fixed_vacuum_charge_origins']={s+'_'+m:b for (s,m),b in bases.items()}
        result['status']='passed' if all(c['pass'] for c in result['checks']) else 'failed_qualification'
    except Exception as exc:
        result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick,native_calls_attempted=len(result['native_calls']),
                  allocation={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','SLURM_MEM_PER_CPU')})
    write(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')}),flush=True)
    return 0 if result['status']=='passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
