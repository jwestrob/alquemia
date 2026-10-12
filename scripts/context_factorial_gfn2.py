"""Finite, separately receipted Linux GFN2 context factorial; no MACE or geometry search."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')
for key in THREAD_KEYS:
    os.environ[key] = '1'
Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'S': 16, 'Ca': 20, 'La': 57}
VALENCE = {'H': 1, 'C': 4, 'N': 5, 'O': 6, 'S': 6, 'Ca': 2, 'La': 3}
BOHR_A = 0.529177210903


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(p):
    if isinstance(p, dict):
        assert pin(p['path']) == p, ('pin mismatch', p)
        p = p['path']
    return json.loads(Path(p).read_text())


def write(p, data):
    with Path(p).open('x') as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write('\n')


def coordinates(p):
    if isinstance(p, dict):
        assert pin(p['path']) == p
        p = p['path']
    lines = Path(p).read_text().splitlines()
    rows = [line.split() for line in lines[2:] if line.strip()]
    assert int(lines[0]) == len(rows)
    return [x[0] for x in rows], [[float(v) for v in x[1:4]] for x in rows]


def clean_env(param):
    env = os.environ.copy()
    for key in THREAD_KEYS:
        env[key] = '1'
    env.update(OMP_DYNAMIC='FALSE', OMP_MAX_ACTIVE_LEVELS='1', XTBPATH=str(Path(param).parent))
    for key in ('LD_PRELOAD', 'LD_LIBRARY_PATH', 'PYTHONPATH', 'PYTHONHOME'):
        env.pop(key, None)
    return env


def dynamic_pins(binary, env):
    run = subprocess.run(['ldd', str(binary)], env=env, text=True, capture_output=True)
    assert run.returncode == 0, run.stderr
    paths = re.findall(r'(?:=>\s+)?(/[^\s()]+)', run.stdout)
    return {'ldd': run.stdout, 'libraries': [pin(p) for p in sorted(set(paths)) if Path(p).is_file()]}


def prepare(args):
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False)
    prepared = read(args.preparation)
    assert prepared['status'] == 'prepared_no_energy_calls'
    snapshot = Path(args.snapshot)
    comparison = snapshot / 'la_followup/strict208_recompute/reference_edge_context_comparison.tsv'
    oldrows = list(csv.DictReader(comparison.open(), delimiter='\t'))
    runtimes = {'xtb': pin(args.xtb), 'xtb_parameter_file': pin(args.parameter_file),
                'tblite_python': pin(args.tblite_python), 'tblite_launcher_path': str(Path(args.tblite_python).absolute()),
                'driver': pin(__file__)}
    env = clean_env(args.parameter_file)
    runtimes['xtb_dynamic'] = dynamic_pins(args.xtb, env)
    inspector = '''import importlib.metadata,json,hashlib,pathlib,tblite.interface,tblite.library
def pin(p):
 p=pathlib.Path(p).resolve();return {"path":str(p),"sha256":hashlib.sha256(p.read_bytes()).hexdigest()}
data={"versions":{},"distribution_files":[],"shared_objects":[]}
for d in importlib.metadata.distributions():
 data["versions"][d.metadata["Name"]]=d.version
 for f in d.files or []:
  p=pathlib.Path(d.locate_file(f))
  if p.is_file() and not p.name.endswith('.pyc'):
   data["distribution_files"].append(pin(p))
   if '.so' in p.name:data["shared_objects"].append(str(p.resolve()))
print(json.dumps(data))'''
    found = subprocess.run([args.tblite_python, '-c', inspector], env=env, text=True, capture_output=True, check=True)
    runtimes['tblite_distribution'] = json.loads(found.stdout)
    assert runtimes['tblite_distribution']['versions']['tblite'] == '0.7.0'
    runtimes['tblite_dynamic'] = {p: dynamic_pins(p, env) for p in runtimes['tblite_distribution']['shared_objects']}
    write(out / 'RUNTIME.json', runtimes)
    tasks = []
    for row in prepared['rows']:
        protein, ctx, geom = row['protein'], row['target_context'], row['geometry_id']
        case = next(c for c in prepared['cases'] if c['protein'] == protein)
        pool = read(case['contexts'][ctx]['manifest'])
        pool_case = next(c for c in pool['cases'] if c['case_id'] == (protein if ctx == 'original' else case['source_case']+'__envelope_'+('5eba0c9b39ead1ee' if protein.startswith('q9') else 'aa36ec84c4cbdc62')))
        same_context = row['previously_scored_same_context']
        candidate = 'origin' if geom == 'origin' else geom.removeprefix(ctx+'_')
        for metal in ('Ca', 'La'):
            src = pool_case['matrix'][metal][candidate]['xyz'] if same_context else row['endpoints'][metal]['xyz']
            symbols, positions = coordinates(src)
            assert symbols[0] == metal and all(s in Z for s in symbols)
            charge = -3 if metal == 'Ca' else -2
            electrons = sum(VALENCE[s] for s in symbols) - charge
            assert electrons % 2 == 0
            for medium in ('vacuum', 'water'):
                tid = '__'.join((protein, ctx, geom, metal, medium))
                old = next((r for r in oldrows if r['protein_id'] == protein and r['panel'] == ('strict208' if ctx == 'original' else 'current34_39')), None)
                ref = None
                if old and geom == ctx+'_adaptive_La':
                    assert src['sha256'] == old[metal+'_xyz_sha256'], 'reference input hash mismatch'
                    ref = {'scope': 'archived_corrected_Mac_standalone_selected_cell', 'table': pin(comparison),
                           'energy_hartree': float(old[f'{metal}_{"vacuum" if medium == "vacuum" else "alpb"}_standalone_energy_hartree'])}
                tasks.append({'id': tid, 'protein': protein, 'context': ctx, 'geometry': geom, 'metal': metal,
                              'medium': medium, 'xyz': src, 'charge': charge, 'multiplicity': 1, 'uhf': 0,
                              'valence_electrons': electrons, 'atom_count': len(symbols),
                              'reference': ref, 'scout': protein.startswith('q9') and ctx == 'original' and geom == 'original_adaptive_La'})
    assert len(tasks) == 80 and sum(t['scout'] for t in tasks) == 4
    manifest = {'protocol_id': 'Linux_xtb671_tblite070_boundary_factorial_v1', 'runtime': pin(out/'RUNTIME.json'),
        'preparation': pin(args.preparation), 'implementation': pin(__file__), 'tasks': tasks,
        'settings': {'accuracy': 0.0001, 'maxiter': 500, 'temperature_K': 300, 'temperature_Eh': 300*3.166808578545117e-6,
          'energy_parity_Eh': 1e-6, 'population_parity_e': 1e-5, 'charge_sum_tolerance_e': 1e-5,
          'bohr_A': BOHR_A, 'workers': 2, 'numerical_threads_per_worker': 1, 'solvation': 'ALPB water gsolv',
          'primary': 'xtb6.7.1 native GFN2, not --tblite', 'independent': 'isolated tblite0.7.0',
          'require_independent_parity': True, 'no_restart': True, 'retry_policy': 'none'},
        'reuse': {'corrected_raw_archives_absent': True, 'basis': 'selected snapshot lacks raw per-candidate outputs; literal Mac paths absent',
                 'new_Linux_endpoint_cap_per_backend': 80, 'scouts_counted_and_reused': 4}}
    write(out/'manifest.json', manifest)
    print(json.dumps({'manifest': str(out/'manifest.json'), 'tasks': len(tasks)}))


def tblite_point(args):
    import numpy as np
    from tblite.interface import Calculator
    import importlib.metadata
    assert importlib.metadata.version('tblite') == '0.7.0'
    task = read(args.task); symbols, positions = coordinates('core.xyz')
    calc = Calculator('GFN2-xTB', np.array([Z[s] for s in symbols]), np.array(positions)/BOHR_A,
                      charge=task['charge'], uhf=task['uhf'])
    calc.set('accuracy', .0001); calc.set('max-iter', 500)
    calc.set('temperature', 300*3.166808578545117e-6); calc.set('verbosity', 2)
    if task['medium'] == 'water':
        calc.add('alpb-solvation', 'water')
    result = calc.singlepoint()
    write('scientific.json', {'energy_hartree': float(result.get('energy')), 'charges_e': result.get('charges').tolist(),
          'orbital_occupations_e': result.get('orbital-occupations').tolist(), 'library_returned_without_error': True})


def call(command, directory, env):
    start = time.time()
    with (directory/'stdout.txt').open('x') as out, (directory/'stderr.txt').open('x') as err:
        proc = subprocess.Popen(command, cwd=directory, env=env, stdout=out, stderr=err, start_new_session=True)
        write(directory/'START.json', {'pid': proc.pid, 'process_group': proc.pid, 'command': command,
              'started_epoch': start, 'thread_environment': {k: env[k] for k in THREAD_KEYS}})
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
        try:
            os.killpg(proc.pid, 0); empty = False
        except ProcessLookupError:
            empty = True
        receipt = {'command': command, 'returncode': proc.returncode, 'wall_seconds': time.time()-start,
                   'cpu_user_seconds': usage.ru_utime, 'cpu_system_seconds': usage.ru_stime,
                   'peak_RSS_KiB': usage.ru_maxrss, 'owned_process_group_empty': empty,
                   'stdout': pin(directory/'stdout.txt'), 'stderr': pin(directory/'stderr.txt')}
        write(directory/'EXECUTION.json', receipt)
        if not empty:
            raise RuntimeError('owned child process group survived completion; stop scheduling')
        return receipt


def parse_backend(directory, backend, task):
    receipt = read(directory/'EXECUTION.json'); stdout = (directory/'stdout.txt').read_text()
    data = {'status': 'unavailable', 'energy_hartree': None, 'charges_e': None, 'execution': pin(directory/'EXECUTION.json')}
    try:
        if receipt['returncode'] != 0:
            raise ValueError(f'exit {receipt["returncode"]}')
        if backend == 'xtb':
            energies = re.findall(r'TOTAL ENERGY\s+(-?\d+\.\d+)', stdout)
            assert energies and re.search(r'convergence criteria satisfied|SCC converged', stdout, re.I), 'missing SCC convergence'
            assert 'FAILED TO CONVERGE' not in stdout.upper()
            energy = float(energies[-1]); charges = [float(x) for x in (directory/'charges').read_text().split()]
            tolerance = {'energy': re.findall(r'(?i)energy convergence\s*[:=]?\s*([\d.Ee+-]+)', stdout),
                         'charge': re.findall(r'(?i)charge convergence\s*[:=]?\s*([\d.Ee+-]+)', stdout)}
            residual_lines = [line for line in stdout.splitlines() if re.match(r'^\s*\d+\s+-?\d+\.\d+\s+', line)]
            evidence = {'convergence_declared': True, 'printed_tolerances': tolerance, 'iteration_lines': residual_lines}
        else:
            result = read(directory/'scientific.json')
            energy, charges = result['energy_hartree'], result['charges_e']
            occupations = result['orbital_occupations_e']; occupation_sum = sum(occupations)
            assert abs(occupation_sum-task['valence_electrons']) <= 1e-5, 'occupation electron mismatch'
            assert result['library_returned_without_error']
            residual_lines = [line for line in stdout.splitlines() if re.match(r'^\s*\d+\s+-?\d+\.\d+\s+', line)]
            evidence = {'library_returned_without_error': True, 'orbital_occupation_sum': occupation_sum,
                        'iteration_lines': residual_lines}
        assert len(charges) == task['atom_count'] and all(math.isfinite(x) for x in charges)
        assert math.isfinite(energy)
        charge_sum = sum(charges); assert abs(charge_sum-task['charge']) <= 1e-5, 'charge sum mismatch'
        data.update(status='available', energy_hartree=energy, charges_e=charges, charge_sum_e=charge_sum,
                    charge_sum_residual_e=charge_sum-task['charge'], solver_evidence=evidence)
    except Exception as exc:
        data['failure_reason'] = str(exc)
    return data


def execute_point(task, manifest, base, runtimes):
    cell = base/'points'/task['id']
    if (cell/'RESULT.json').exists():
        old = read(cell/'RESULT.json'); assert old['task'] == task
        return old
    cell.mkdir(parents=True, exist_ok=False)
    write(cell/'TASK.json', task)
    env = clean_env(runtimes['xtb_parameter_file']['path'])
    for backend in ('xtb', 'tblite'):
        directory = cell/backend; directory.mkdir()
        shutil.copyfile(task['xyz']['path'], directory/'core.xyz')
        assert pin(directory/'core.xyz')['sha256'] == task['xyz']['sha256']
        if backend == 'xtb':
            command = [runtimes['xtb']['path'], 'core.xyz', '--gfn', '2', '--chrg', str(task['charge']),
                       '--uhf', '0', '--etemp', '300', '--acc', '0.0001', '--iterations', '500', '--norestart', '--json', '--pop']
            if task['medium'] == 'water': command += ['--alpb', 'water', 'gsolv']
        else:
            command = [runtimes['tblite_launcher_path'], str(Path(__file__).resolve()), 'tblite-point', '--task', str(cell/'TASK.json')]
        call(command, directory, env)
    parsed = {b: parse_backend(cell/b, b, task) for b in ('xtb', 'tblite')}
    result = {'task': task, 'implementation': pin(__file__), 'manifest': pin(base/'manifest.json'),
              'backends': parsed, 'accepted': False, 'energy_hartree': None, 'parity': None}
    if all(r['status'] == 'available' for r in parsed.values()):
        de = parsed['xtb']['energy_hartree']-parsed['tblite']['energy_hartree']
        dq = max(abs(x-y) for x,y in zip(parsed['xtb']['charges_e'],parsed['tblite']['charges_e']))
        result['parity'] = {'energy_difference_hartree': de, 'maximum_population_difference_e': dq}
        result['accepted'] = abs(de)<=1e-6 and dq<=1e-5
        if result['accepted']: result['energy_hartree'] = parsed['xtb']['energy_hartree']
        if task['reference']:
            result['reference_energy_difference_hartree'] = parsed['xtb']['energy_hartree']-task['reference']['energy_hartree']
            if abs(result['reference_energy_difference_hartree'])>1e-6:
                result['accepted']=False;result['energy_hartree']=None
    write(cell/'RESULT.json', result)
    print(json.dumps({'id':task['id'], 'accepted':result['accepted'], 'parity':result['parity']}), flush=True)
    return result


def collect(base, output):
    manifest=read(base/'manifest.json'); rows=[]
    for task in manifest['tasks']:
        path=base/'points'/task['id']/'RESULT.json'
        rows.append(read(path) if path.exists() else {'task':task,'accepted':False,'status':'not_executed','energy_hartree':None})
    write(output, {'protocol_id':manifest['protocol_id'],'manifest':pin(base/'manifest.json'),'rows':rows,
                   'accepted':sum(r['accepted'] for r in rows),'declared':len(rows),'executed':sum('backends' in r for r in rows)})


def execute(args):
    base=Path(args.directory).resolve();manifest=read(base/'manifest.json');runtimes=read(manifest['runtime'])
    assert pin(__file__)==manifest['implementation'], 'implementation changed after freeze'
    for key in ('xtb','xtb_parameter_file','tblite_python'):
        assert pin(runtimes[key]['path'])==runtimes[key]
    if args.stage=='remaining':
        scout=read(base/'SCOUT.json');assert scout['accepted']==scout['declared']==4,'scout failed'
    tasks=[t for t in manifest['tasks'] if t['scout']==(args.stage=='scout')]
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows=list(pool.map(lambda task:execute_point(task,manifest,base,runtimes),tasks))
    if args.stage=='scout':
        write(base/'SCOUT.json',{'declared':4,'accepted':sum(r['accepted'] for r in rows),'rows':rows})
    collect(base,base/('SCOUT_COLLECTION.json' if args.stage=='scout' else 'FINAL_COLLECTION.json'))


def main():
    parser=argparse.ArgumentParser();sub=parser.add_subparsers(dest='op',required=True)
    p=sub.add_parser('prepare')
    for k in ('preparation','snapshot','xtb','parameter-file','tblite-python','output'):p.add_argument('--'+k,required=True)
    p=sub.add_parser('tblite-point');p.add_argument('--task',required=True)
    p=sub.add_parser('execute');p.add_argument('--directory',required=True);p.add_argument('--stage',choices=['scout','remaining'],required=True)
    p=sub.add_parser('collect');p.add_argument('--directory',required=True);p.add_argument('--output',required=True)
    args=parser.parse_args()
    if args.op=='prepare':prepare(args)
    elif args.op=='tblite-point':tblite_point(args)
    elif args.op=='execute':execute(args)
    else:collect(Path(args.directory).resolve(),args.output)


if __name__=='__main__':main()
