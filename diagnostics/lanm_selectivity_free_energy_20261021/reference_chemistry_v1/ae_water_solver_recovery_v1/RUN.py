#!/usr/bin/env python3
"""One Dy-water solver recovery; preparation is saved-only, execution Slurm-only."""
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
WORK = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1'
OLD = WORK/'pilot_v1/tasks/ae_pbe0__water_central_Dy/execution'
CELL = 'ae_pbe0__water_central_Dy'
INSERT = ' MaxIter 300\n trah\n  QuadRegionStart 1e-6\n end\nend\n%basis'
ORIGINAL = ' MaxIter 300\nend\n%basis'
FLOAT = r'([-+0-9.eE]+)'


def check(ok, why):
    if not ok:
        raise ValueError(why)


def pin(path):
    path = Path(path).resolve()
    with path.open('rb') as f:
        sha = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'path': str(path), 'sha256': sha}


def verify(record):
    check(pin(record['path']) == record, 'Changed pin: '+record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, data):
    with Path(path).open('x') as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write('\n')
    return pin(path)


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def default_config():
    return {'implementation': pin(__file__), 'declaration': pin(HERE/'ROOT_DECLARATION.md'),
            'native_adapter': pin(HERE.parent/'REFERENCE.py'),
            'panel': pin(WORK/'prepared_v2/PANEL.json'),
            'original_start': pin(OLD/'START.json'), 'original_input': pin(OLD/'endpoint.inp'),
            'original_output': pin(OLD/'endpoint.out'), 'original_result': pin(OLD/'RESULT.json'),
            'assessment': pin(HERE.parent/'ae_water_recovery_assessment_v1/ASSESSMENT.json'),
            'documentation': pin(HERE.parent/'ae_water_recovery_assessment_v1/PRIMARY_SOURCE.md'),
            'launcher': pin(HERE.parent/'RUN.sbatch'), 'cell_id': CELL,
            'MPI_ranks': 4, 'CPUs_per_task': 1, 'max_native_invocations': 1,
            'QuadRegionStart': 1e-6, 'TolG': 2e-6, 'MaxIter': 300}


def load(config_path):
    cfg = json.loads(Path(config_path).read_text())
    check(cfg == default_config(), 'Configuration must match the declared single cell and exact pins')
    for value in cfg.values():
        if isinstance(value, dict) and 'sha256' in value:
            verify(value)
    spec = importlib.util.spec_from_file_location('ae_water_original_native', cfg['native_adapter']['path'])
    native = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(native)
    panel = native.load_panel(verify(cfg['panel']))
    cell, = [c for c in panel['cells'] if c['cell_id'] == CELL]
    old = read(cfg['original_start'])
    check(cell == old['cell'], 'Original cell changed')
    check((cell['charge'], cell['multiplicity'], cell['physical_multiplicity'], cell['explicit_electrons'],
           cell['ecp_core'], cell['geometry']['atoms']) == (3, 6, 6, 73, 0, 4), 'State/size changed')
    check(old['panel'] == cfg['panel'] and old['job_id'] == '17260', 'Original authority differs')
    check(native.deck(panel, cell, 4, 3000) == verify(cfg['original_input']).read_text(), 'Original deck replay differs')
    check(read(cfg['original_result'])['status'] == 'failed_collection', 'Original failure changed')
    return cfg, native, panel, cell


def recovery_deck(native, panel, cell, ranks, maxcore):
    base = native.deck(panel, cell, ranks, maxcore)
    check(base.count(ORIGINAL) == 1 and 'QuadRegionStart' not in base, 'Unexpected original SCF block')
    changed = base.replace(ORIGINAL, INSERT)
    check(changed.replace(INSERT, ORIGINAL) == base, 'More than the declared solver setting changed')
    return base, changed


def scf_evidence(text):
    lines = text.splitlines()
    def numbers(pattern):
        return [float(x) for x in re.findall(pattern, text)]
    tol = numbers(r'Orbital Gradient\s+TolG\s+\.{2,}\s*'+FLOAT)
    trah_tol = numbers(r'Converg\. threshold\s+\(grad\. norm\)\s+\.{2,}\s*'+FLOAT)
    entry = numbers(r'NR start threshold \(gradient norm\)\s+\.{2,}\s*'+FLOAT)
    final = re.findall(r'Last Orbital Gradient\s+\.{2,}\s*'+FLOAT+r'\s+Tolerance\s*:\s*'+FLOAT, text)
    macro = []
    for line_no, line in enumerate(lines, 1):
        m = re.match(r'^\s*(\d+)\s+(-\d+\.\d+)\s+([0-9.eE+-]+)\s+.*\((TRAH|NR)\s+MAcro\)', line)
        if m:
            macro.append({'line': line_no, 'iteration': int(m[1]), 'orbital_gradient_norm': float(m[3]), 'mode': m[4]})
    used = bool(macro)
    native_success = bool(re.search(r'SCF CONVERGED AFTER', text))
    signals = [line.strip() for line in lines if 'check signals convergence' in line]
    final_gradient = float(final[-1][0]) if final else None
    residual = macro[-1]['orbital_gradient_norm'] if used else final_gradient
    target = (trah_tol[-1] if trah_tol else None) if used else (tol[-1] if tol else None)
    within = residual is not None and target is not None and math.isfinite(residual) and residual <= target
    # TRAH's macro norm is active. For ordinary SCF retain its native stopping path;
    # do not reinterpret stale DIIS/density flags as additional criteria.
    gradient_stop = any('Gradient check' in s for s in signals)
    active = native_success and (within if used or gradient_stop else True)
    criteria = [line.strip() for line in lines if re.search(r'Last (?:Energy|MAX-Density|RMS-Density|DIIS|Orbital)', line)]
    state_markers = [{'line': i, 'text': line.strip()} for i, line in enumerate(lines, 1)
                     if re.search(r'Expectation value of <S\*\*2>|ORBITAL ENERGIES|ATOMIC CHARGES|SPIN POPULATIONS|N\(Alpha\)|N\(Beta\)', line)]
    return {'native_scf_converged': native_success, 'TRAH_used': used, 'macro_rows': len(macro),
            'last_macro': macro[-1] if macro else None, 'declared_TolG': tol[-1] if tol else None,
            'printed_TRAH_convergence_threshold': trah_tol[-1] if trah_tol else None,
            'active_orbital_gradient_threshold': target,
            'printed_ConvCheckMode': re.findall(r'Convergence Check Mode\s+ConvCheckMode\s+\.{2,}\s*([^\n]+)', text),
            'printed_NR_entry_threshold': entry[-1] if entry else None,
            'final_printed_orbital_gradient': final_gradient,
            'final_printed_orbital_gradient_tolerance': float(final[-1][1]) if final else None,
            'active_residual_source': 'last_TRAH_or_NR_macro_norm' if used else 'final_Last_Orbital_Gradient',
            'orbital_residual_within_TolG': within, 'active_scf_convergence': active,
            'native_stopping_signals': signals, 'printed_final_criteria': criteria,
            'micro_residual_used_for_admission': False,
            'solver_change_causality': 'not established' if used else 'TRAH not used; success cannot be attributed to new control',
            'state_analysis_markers_in_raw_output': state_markers,
            'physical_state_qualified': False}


def prepare(config_path, out):
    wall, cpu = time.monotonic(), time.process_time()
    cfg, native, panel, cell = load(config_path)
    out.mkdir(parents=True, exist_ok=False)
    base, changed = recovery_deck(native, panel, cell, 4, 3000)
    (out/'original.template.inp').write_text(base)
    (out/'recovery.template.inp').write_text(changed)
    shutil.copyfile(native.verify(cell['geometry']['xyz']), out/'geometry.xyz')
    manifest = {'status': 'prepared', 'config': pin(config_path), 'implementation': pin(__file__),
                'cell': cell, 'native_adapter': cfg['native_adapter'], 'panel': cfg['panel'],
                'original_template': pin(out/'original.template.inp'), 'recovery_template': pin(out/'recovery.template.inp'),
                'geometry': pin(out/'geometry.xyz'), 'template_MaxCore_MiB': 3000,
                'resource_rule': 'floor(0.75 * actual_granted_MiB / 4)',
                'sole_scientific_delta': 'nested TRAH QuadRegionStart 1e-6',
                'orca': panel['orca'], 'original_native_failure': cfg['original_result'],
                'native_calls': 0, 'physical_state_qualified': False}
    mp = write(out/'MANIFEST.json', manifest)
    rp = write(out/'RESULT.json', {'status': 'preparation_complete', 'manifest': mp,
               'original_deck_exact': True, 'sole_delta_proved': True, 'cells': 1, 'native_calls': 0,
               'local_threads': 1, 'wall_seconds': time.monotonic()-wall, 'CPU_seconds': time.process_time()-cpu})
    release = {'status': None, 'config': pin(config_path), 'manifest': mp, 'preparation': rp,
               'implementation': pin(__file__), 'independent_review': None, 'cell_id': CELL,
               'native_invocations': 1, 'MPI_ranks': 4, 'CPUs_per_task': 1, 'root_release': None}
    write(out/'ROOT_RELEASE.template.json', release)
    print(json.dumps({'manifest': mp, 'preparation': rp}))


def admit(cfg, native, panel, cell, directory, raw):
    text = (directory/'endpoint.out').read_text(errors='replace')
    scf = scf_evidence(text)
    def exact(pattern, expected):
        found = re.findall(pattern, text)
        return bool(found) and found[-1] == expected
    checks = {'original_native_admission': raw.get('status') == 'complete',
              'SCF_active_convergence': scf['active_scf_convergence'],
              'declared_TolG_unchanged': scf['declared_TolG'] == 2e-6,
              'TRAH_convergence_threshold_unchanged': scf['printed_TRAH_convergence_threshold'] == 2e-6,
              'NR_entry_threshold_executed': scf['printed_NR_entry_threshold'] == 1e-6,
              'native_charge': exact(r'Total Charge\s+Charge\s+\.{2,}\s*([0-9]+)', '3'),
              'native_basis_dimension': exact(r'Basis Dimension\s+Dim\s+\.{2,}\s*([0-9]+)', '164'),
              'native_unrestricted': exact(r'Hartree-Fock type\s+HFTyp\s+\.{2,}\s*(\w+)', 'UHF'),
              'native_MaxIter': exact(r'Maximum # iterations\s+MaxIter\s+\.{2,}\s*([0-9]+)', '300'),
              'native_version': 'Program Version 6.1.1' in text,
              'finite_E_and_forces': False}
    if raw.get('energy_hartree') is not None and raw.get('forces_kcal_mol_A'):
        vals = [raw['energy_hartree']] + [x for v in raw['forces_kcal_mol_A'] for x in v]
        checks['finite_E_and_forces'] = len(vals) == 13 and all(math.isfinite(x) for x in vals)
    base, changed = recovery_deck(native, panel, cell, 4, raw['MaxCore_MiB_per_rank'])
    checks['actual_input_sole_solver_delta'] = (directory/'endpoint.inp').read_text() == changed
    return {'status': 'recovery_admitted' if all(checks.values()) else 'failed_recovery_admission',
            'cell_id': CELL, 'checks': checks, 'SCF_evidence': scf, 'native_result': pin(directory/'RESULT.json'),
            'native_output': pin(directory/'endpoint.out'), 'native_input': pin(directory/'endpoint.inp'),
            'original_failure': cfg['original_result'], 'wrapper': pin(__file__), 'native_helper': cfg['native_adapter'],
            'input_relation': 'Original deck plus QuadRegionStart; NOT an unchanged baseline execution',
            'native_invocations': 1, 'physical_state_qualified': False, 'method_accuracy_qualified': False,
            'force_accuracy_qualified': False, 'automatic_retries': 0,
            'energy_hartree': raw.get('energy_hartree') if all(checks.values()) else None,
            'forces_kcal_mol_A': raw.get('forces_kcal_mol_A') if all(checks.values()) else None}


def run(config_path, manifest_path, release_path, out):
    cfg, native, panel, cell = load(config_path)
    manifest = json.loads(Path(manifest_path).read_text()); release = json.loads(Path(release_path).read_text())
    check(release['status'] == 'approved_single_AE_Dy_water_solver_recovery', 'Explicit root release required')
    for key, value in {'config': pin(config_path), 'manifest': pin(manifest_path), 'implementation': pin(__file__),
                       'cell_id': CELL, 'native_invocations': 1, 'MPI_ranks': 4, 'CPUs_per_task': 1}.items():
        check(release[key] == value, 'Release mismatch: '+key)
    check(release.get('root_release'), 'Missing root authority')
    verify(release['independent_review']); prep = read(release['preparation'])
    check(prep['status'] == 'preparation_complete' and prep['manifest'] == pin(manifest_path), 'Preparation failed/differs')
    check(manifest['config'] == pin(config_path) and manifest['cell'] == cell and manifest['orca'] == panel['orca'], 'Manifest mismatch')
    base, template = recovery_deck(native, panel, cell, 4, 3000)
    check(verify(manifest['original_template']).read_text() == base, 'Baseline template differs')
    check(verify(manifest['recovery_template']).read_text() == template, 'Recovery template differs')
    check(verify(manifest['geometry']).read_bytes() == native.verify(cell['geometry']['xyz']).read_bytes(), 'Geometry differs')
    check(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Native execution requires Slurm')
    check(int(os.environ.get('SLURM_NTASKS', '0')) == 4 and int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) == 1, 'Expected four MPI ranks, one thread/rank')
    mem = int(os.environ.get('SLURM_MEM_PER_NODE', '0')) or 4*int(os.environ.get('SLURM_MEM_PER_CPU', '0'))
    check(mem > 0, 'Actual memory grant required'); maxcore = math.floor(.75*mem/4)
    out.mkdir(parents=True, exist_ok=False); directory = out/'native'; directory.mkdir()
    base, changed = recovery_deck(native, panel, cell, 4, maxcore)
    (out/'baseline_at_actual_resources.inp').write_text(base)
    (directory/'endpoint.inp').write_text(changed)
    shutil.copyfile(native.verify(cell['geometry']['xyz']), directory/'geometry.xyz')
    env = os.environ.copy()
    for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
        env[key] = '1'
    env['OMP_DYNAMIC'] = env['MKL_DYNAMIC'] = 'FALSE'
    env['PATH'] = str(Path(panel['orca']['path']).parent)+os.pathsep+env.get('PATH', '')
    cmd = [str(native.verify(panel['orca'])), 'endpoint.inp', '--bind-to none']
    write(directory/'START.json', {'utc': now(), 'host': socket.gethostname(), 'job_id': os.environ['SLURM_JOB_ID'],
          'cell': cell, 'wrapper': pin(__file__), 'native_helper': cfg['native_adapter'], 'panel': cfg['panel'],
          'release': pin(release_path), 'manifest': pin(manifest_path), 'input': pin(directory/'endpoint.inp'),
          'baseline_at_actual_resources': pin(out/'baseline_at_actual_resources.inp'),
          'sole_input_delta': 'nested TRAH QuadRegionStart 1e-6', 'ranks': 4, 'memory_MiB': mem,
          'MaxCore_MiB_per_rank': maxcore, 'command': cmd,
          'thread_environment': {k: env[k] for k in env if k.endswith('NUM_THREADS') or k == 'OMP_THREAD_LIMIT'}})
    wall = time.monotonic(); error = None; returncode = None
    try:
        with (directory/'endpoint.out').open('x') as stdout, (directory/'endpoint.err').open('x') as stderr:
            proc = subprocess.run(cmd, cwd=directory, env=env, stdout=stdout, stderr=stderr, check=False)
        returncode = proc.returncode
        raw = native.collect_one(panel, cell, directory, returncode)
    except Exception as exc:
        error = repr(exc); raw = {'status': 'failed_collection', 'cell_id': CELL, 'exception': error}
    elapsed = time.monotonic()-wall
    raw.update({'returncode': returncode, 'elapsed_seconds': elapsed, 'allocated_cpu_seconds': elapsed*4,
                'job_id': os.environ['SLURM_JOB_ID'], 'finished_utc': now(), 'panel': cfg['panel'],
                'native_collection_helper': cfg['native_adapter'], 'execution_wrapper': pin(__file__),
                'MaxCore_MiB_per_rank': maxcore, 'input_delta': 'QuadRegionStart 1e-6'})
    write(directory/'RESULT.json', raw)
    try:
        result = admit(cfg, native, panel, cell, directory, raw)
    except Exception as exc:
        result = {'status': 'failed_recovery_collection', 'cell_id': CELL, 'exception': repr(exc),
                  'native_result': pin(directory/'RESULT.json'), 'physical_state_qualified': False}
    result.update({'job_id': os.environ['SLURM_JOB_ID'], 'release': pin(release_path), 'manifest': pin(manifest_path),
                   'elapsed_seconds': elapsed, 'allocated_cpu_seconds': elapsed*4, 'finished_utc': now()})
    write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] == 'recovery_admitted' else 1


def main():
    parser = argparse.ArgumentParser(); sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('config'); p.add_argument('--out', type=Path, required=True)
    for name in ('prepare', 'run'):
        p = sub.add_parser(name); p.add_argument('--config', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
        if name == 'run':
            p.add_argument('--manifest', type=Path, required=True); p.add_argument('--release', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'config': write(args.out, default_config())
    elif args.action == 'prepare': prepare(args.config, args.out)
    else: return run(args.config, args.manifest, args.release, args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
