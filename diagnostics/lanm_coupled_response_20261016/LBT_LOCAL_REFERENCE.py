#!/usr/bin/env python3
"""Finite LBT48 reference using the unchanged pinned compact ORCA primitives."""
import argparse
import datetime as dt
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import re
import resource
import socket
import subprocess
import sys
import time
import traceback

PROTOCOL = 'lbt3_local_source_response_reference_v1'
IDS = ('7CCO_La', '7CCN_La', '7CCO_Dy', '7CCN_Dy')
COMPACT_SHA = '18ed7c0f4141b95cba20191bf36867bc839ae8722945b81dfb74fe3eeeda7f25'
METHOD = 'PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    p = Path(record['path'])
    require(pin(p)['sha256'] == record['sha256'], 'Pin changed: ' + str(p))
    return p


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def inputs(cfg):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm is required, including preflight')
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Adapter differs')
    verify(cfg['plan']); verify(cfg['health_monitor'])
    require(pin(sys.executable)['sha256'] == cfg['python']['sha256'], 'Python executable differs')
    compact = read(cfg['compact_manifest'])
    require(cfg['compact_manifest']['sha256'] == COMPACT_SHA, 'Compact parent differs')
    for rec in compact['implementation'].values():
        verify(rec)
    verify(compact['orca'])
    for assets in compact['assets'].values():
        for rec in assets.values():
            verify(rec)
    directory = verify(compact['implementation']['metal_environment_frozen_f.py']).parent
    sys.path.insert(0, str(directory))
    model = importlib.import_module('metal_environment_frozen_f')
    runner = importlib.import_module('run_orca_task_manifest')
    collector = importlib.import_module('ggr_sensitivity')
    residual = importlib.import_module('mace_omol_vacuum')
    runtime = importlib.import_module('render_orca_runtime_input')
    for name, obj in [('metal_environment_frozen_f.py', model), ('run_orca_task_manifest.py', runner),
                      ('ggr_sensitivity.py', collector), ('mace_omol_vacuum.py', residual),
                      ('render_orca_runtime_input.py', runtime)]:
        require(pin(obj.__file__)['sha256'] == compact['implementation'][name]['sha256'], 'Imported helper differs: ' + name)
    require(model.METHOD == METHOD, 'Compact method changed')
    prep = read(cfg['preparation_result']); mace = read(cfg['MACE_result'])
    mcfg = read(mace['config'])
    require(prep['status'] == mace['status'] == 'complete'
            and prep['slurm_job_id'] == '12469' and mace['slurm_job_id'] == '12471'
            and mcfg['preparation'] == cfg['preparation_result'], 'Preparation/MACE parent differs')
    prepared = {t['source_id'] + '_' + t['metal']: t for r in prep['rows'] for t in r['tasks']
                if t['representation'] == 'matched_core'}
    saved = {r['source_id'] + '_' + r['metal']: r for r in mace['rows']
             if r['representation'] == 'matched_core'}
    require(set(prepared) == set(saved) == set(cfg['tasks']) == set(IDS), 'Four-state denominator differs')
    states = {}
    for tid in IDS:
        task = cfg['tasks'][tid]
        require(task == prepared[tid], 'Prepared task changed: ' + tid)
        for key in ('xyz', 'physical_source', 'source_cap_map'):
            verify(task[key])
            require(saved[tid][key] == task[key], 'MACE source/map differs')
        require(saved[tid]['status'] == 'complete' and task['atom_count'] == 48
                and task['charge'] == -1 and task['metal_index'] == 0
                and task['multiplicity'] == (1 if task['metal'] == 'La' else 6), 'Physical state differs')
        atoms = model.xyz(verify(task['xyz']))
        require(len(atoms) == 48 and atoms[0][0] == task['metal'], 'Actual core inventory differs')
        st = model.state(verify(task['xyz']))
        require(st['charge'] == -1 and st['effective_multiplicity'] == 1
                and st['physical_multiplicity'] == task['multiplicity']
                and st['explicit_electrons'] == 200 and st['alpha_electrons'] == st['beta_electrons'] == 100
                and st['all_electron_count'] == task['all_electron_count'], 'Derived valence state differs')
        states[tid] = st
    for source in ('7CCO', '7CCN'):
        la = model.xyz(verify(prepared[source + '_La']['xyz']))
        dy = model.xyz(verify(prepared[source + '_Dy']['xyz']))
        require([a[1:] for a in la] == [a[1:] for a in dy]
                and [a[0] for a in la[1:]] == [a[0] for a in dy[1:]], 'Paired real coordinates differ')
    return compact, model, runner, collector, residual, runtime, states


def allocation(runtime, output):
    require(int(os.environ.get('SLURM_NTASKS', '0')) == 12
            and int(os.environ.get('SLURM_CPUS_PER_TASK', '1')) == 1
            and int(os.environ.get('SLURM_CPUS_ON_NODE', '0')) == 12
            and int(os.environ.get('SLURM_JOB_NUM_NODES', '1')) == 1, 'Exactly12 real MPI slots on one node required')
    text = subprocess.check_output(['scontrol', 'show', 'job', os.environ['SLURM_JOB_ID'], '-o'], text=True)
    (output / 'SCONTROL_JOB.txt').write_text(text)
    match = re.search(r'(?:^|\s)AllocTRES=([^\s]+)', text)
    require(match is not None, 'Actual allocated TRES unavailable')
    tres = dict(field.split('=', 1) for field in match[1].split(','))
    require(int(tres['cpu']) == 12, 'Actual CPU grant differs')
    mem = re.fullmatch(r'([0-9]+(?:\.[0-9]+)?)([KMGT]?)', tres.get('mem', ''), re.I)
    require(mem is not None, 'Actual memory grant unavailable')
    granted = float(mem[1]) * {'K': 1/1024, 'M': 1, 'G': 1024, 'T': 1024**2, '': 1}[mem[2].upper()]
    require(granted > 0 and granted == int(granted), 'Unsupported memory grant')
    # This exports the observed grant to the existing renderer; it requests no memory.
    os.environ['SLURM_MEM_PER_NODE'] = str(int(granted))
    os.environ['METAL_ENV_WORKERS'] = '1'
    policy = runtime.allocation_policy(workers=1, ranks=12)
    require(policy['allocated_memory_MiB'] == int(granted)
            and policy['memory_source'] == 'SLURM_MEM_PER_NODE'
            and policy['workers'] == 1 and policy['ranks_per_worker'] == 12,
            'Renderer did not use the actual shared grant')
    return dict(policy, scheduler_receipt=pin(output / 'SCONTROL_JOB.txt'), allocated_TRES=tres)


def stage(cfg, compact, model, states, output, selected):
    full = dict(protocol_id=PROTOCOL, agreement=cfg['plan'], plan=cfg['plan'], method=METHOD,
                adapter=cfg['implementation'], inputs_config=cfg['_config_pin'],
                compact_manifest=cfg['compact_manifest'], preparation_result=cfg['preparation_result'],
                MACE_result=cfg['MACE_result'], assets=compact['assets'], orca=compact['orca'],
                implementation=compact['implementation'], execution_policy=compact['execution_policy'],
                execution_resources={'mpi_ranks': 12, 'concurrent_tasks': 1}, tasks=[])
    for tid in IDS:
        original = cfg['tasks'][tid]; directory = output / tid; directory.mkdir()
        (directory / 'core.xyz').write_bytes(verify(original['xyz']).read_bytes())
        assets = compact['assets'][original['metal']]
        template = model.input_text(verify(assets['basis']).read_text(), verify(assets['aux']).read_text())
        (directory / 'endpoint.inp').write_text(template)
        require(not re.search(r'%maxcore|%pal|DoEQ|pointcharges|CPCM|\bOpt\b', template, re.I), 'Unexpected template control')
        task = dict(task_id=tid, source_id=original['source_id'], metal=original['metal'],
                    charge=-1, multiplicity=1, physical_multiplicity=original['multiplicity'],
                    electronic_state=states[tid], source_task=original,
                    xyz=pin(directory / 'core.xyz'), input=pin(directory / 'endpoint.inp'),
                    output_path=str(directory / 'endpoint.out'), engrad_path=str(directory / 'endpoint.engrad'),
                    task_type='analytic_gradient')
        require(task['xyz']['sha256'] == original['xyz']['sha256'], 'Copy changed coordinates')
        full['tasks'].append(task)
    write(output / 'FULL_MANIFEST.json', full)
    execution = dict(full, full_manifest=pin(output / 'FULL_MANIFEST.json'),
                     tasks=[t for t in full['tasks'] if t['task_id'] in selected])
    write(output / 'MANIFEST.json', execution)
    return full, execution


def observations(output, residual):
    text = Path(output).read_text(errors='replace')
    values = residual.terminal_scf_residuals(text)
    require(all(math.isfinite(v[k]) for v in values.values() for k in ('value', 'printed_tolerance')),
            'Nonfinite SCF evidence')
    for v in values.values():
        v['within_printed_tolerance'] = abs(v['value']) <= v['printed_tolerance']
    modes = re.findall(r'Convergence Check Mode\s+ConvCheckMode\s*\.{2,}\s*([^\n]+)', text)
    forced = re.findall(r'Convergence forced\s+ConvForced\s*\.{2,}\s*(\d+)', text)
    require(len(modes) == len(forced) == 1, 'Native convergence mode evidence missing/ambiguous')
    warnings = [line for line in text.splitlines() if re.search(r'warn|near.degen|small.*gap|homo.lumo', line, re.I)]
    return {'terminal_residuals': values, 'ConvCheckMode': modes[0].strip(), 'ConvForced': int(forced[0]),
            'all_printed_residuals_pass': all(v['within_printed_tolerance'] for v in values.values()),
            'SOSCF_initialized': 'Initializing SOSCF' in text,
            'native_cycles': re.findall(r'SCF CONVERGED AFTER\s+(\d+)\s+CYCLES', text),
            'warnings_and_gap_lines': warnings,
            'admission_policy': 'unchanged native convergence and compact state/gradient parser; residual flags are observations, not new all-criteria admission'}


def admit(manifest_path, model, collector, residual):
    manifest, rows = collector.executed(manifest_path)
    for task in manifest['tasks']:
        row = rows[task['task_id']]
        row.update(task_id=task['task_id'], source_id=task['source_id'], metal=task['metal'],
                   task=task, physical_multiplicity=task['physical_multiplicity'])
        try:
            if Path(task['output_path']).exists():
                row['SCF_evidence'] = observations(task['output_path'], residual)
            if row['status'] != 'complete':
                continue
            parsed = model.parse(task, {'electronic_state': task['electronic_state']})
            require(math.isfinite(parsed['energy_hartree']), 'Nonfinite native energy')
            row.update(parsed, status='available')
            receipt = read(row['receipt'])
            row['molecular_wall_seconds'] = (dt.datetime.fromisoformat(receipt['finished_at_utc'])
                                           - dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds()
            sidecar = read(receipt['artifacts']['runtime_sidecar'])
            row['runtime_memory_policy'] = sidecar['memory_policy']
        except Exception as exc:
            row.update(status='unavailable', failure=repr(exc), energy_hartree=None)
    return rows


def scout_gate(cfg, model, collector, residual):
    prior = read(cfg['scout_result']); pcfg = read(prior['config'])
    require(prior['stage'] == 'scout' and prior['status'] == 'complete'
            and prior['selected_tasks'] == ['7CCO_La'] and prior['new_calls_attempted'] == 1,
            'Admitted actual La7CCO scout required')
    for key in ('plan', 'compact_manifest', 'preparation_result', 'MACE_result', 'tasks', 'python', 'health_monitor'):
        require(pcfg[key] == cfg[key], 'Scout scientific inputs differ: ' + key)
    verify(prior['implementation']); verify(pcfg['implementation'])
    require(prior['implementation']['sha256'] == pcfg['implementation']['sha256']
            == cfg['implementation']['sha256'], 'Scout implementation differs')
    actual = admit(verify(prior['manifest']), model, collector, residual)['7CCO_La']
    require(actual == prior['rows']['7CCO_La'] and actual['status'] == 'available', 'Scout raw admission does not replay')
    return actual


def run(args):
    output = Path(args.output).resolve(); output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic(); result = dict(protocol_id=PROTOCOL, stage=args.stage, status='execution_failure',
        config=pin(args.config), implementation=pin(__file__), slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        host=socket.gethostname(), selected_tasks=[], new_calls_attempted=0, new_MACE_calls=0,
        new_GFN_calls=0, optimizer_steps=0, rows={tid: {'task_id': tid, 'status': 'not_requested'} for tid in IDS})
    try:
        cfg = json.loads(Path(args.config).read_text()); cfg['_config_pin'] = pin(args.config)
        compact, model, runner, collector, residual, runtime, states = inputs(cfg)
        result.update(plan=cfg['plan'], compact_manifest=cfg['compact_manifest'], preparation_result=cfg['preparation_result'],
                      MACE_result=cfg['MACE_result'], states=states, preflight='passed')
        if args.stage == 'remaining':
            result['rows']['7CCO_La'] = scout_gate(cfg, model, collector, residual)
            result['scout_result'] = cfg['scout_result']
        selected = [] if args.stage == 'preflight' else (['7CCO_La'] if args.stage == 'scout'
                    else ([args.task] if args.task else list(IDS[1:])))
        require(args.task is None or args.stage == 'remaining' and args.task in IDS[1:], 'Unsupported task selection')
        result['selected_tasks'] = selected
        result['allocation'] = allocation(runtime, output)
        full, manifest = stage(cfg, compact, model, states, output, selected)
        result.update(full_manifest=pin(output / 'FULL_MANIFEST.json'), manifest=pin(output / 'MANIFEST.json'))
        write(output / 'START.json', result)
        if selected:
            try:
                runner.run_manifest(output / 'MANIFEST.json', orca_path=verify(compact['orca']),
                    workers=1, nprocs=12, selected_tasks=set(selected),
                    expected_runner_sha256=compact['execution_policy']['task_runner']['sha256'])
            except Exception as exc:
                result['runner_error'] = repr(exc)
            finally:
                result['new_calls_attempted'] = sum(Path(t['output_path']).exists() for t in manifest['tasks'])
                result['rows'].update(admit(output / 'MANIFEST.json', model, collector, residual))
        result['status'] = 'complete' if all(result['rows'][t]['status'] == 'available' for t in selected) else 'unavailable'
        if args.stage == 'preflight':
            result['status'] = 'preflight_passed'
        result['full_four_available'] = all(row['status'] == 'available' for row in result['rows'].values())
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    finally:
        result['elapsed_seconds'] = time.monotonic() - start
        result['application_allocated_CPU_seconds'] = result['elapsed_seconds'] * int(os.environ.get('SLURM_CPUS_ON_NODE', '0'))
        result['peak_parent_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result['parent_RSS_is_not_MPI_peak'] = True
        result['physical_state_qualified'] = result['force_consistency_qualified'] = False
        write(output / 'RESULT.json', result)
    print(json.dumps({k: result[k] for k in ('status', 'selected_tasks', 'new_calls_attempted', 'elapsed_seconds')}))
    return 0 if result['status'] in ('complete', 'preflight_passed') else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('preflight', 'scout', 'remaining'), required=True)
    parser.add_argument('--task', choices=IDS[1:])
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    raise SystemExit(run(parser.parse_args()))
