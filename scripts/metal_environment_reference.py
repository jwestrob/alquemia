"""Finite embedded-electronic response scout; no production or hybrid score.

Reuses the native manifested ORCA runner and strict r2SCAN-3c gradient parser.
Inputs must be an audited real-source preparation; never creates MM parameters.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import numpy as np

from affordable_common import (BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, cache_key,
                              read_json, record, verify, write_new, paired)
from affordable_workflow import dry_run, execute
from mace_omol_vacuum import (METHOD, embedded_input, scientific_input,
                              parse_endpoint)

PROTOCOL = 'nikasha_finite_embedded_electronic_response_v1'
ROOT = Path(__file__).resolve().parents[1]
ORCA = Path('/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca')


def read_pcgrad(path, count):
    lines = Path(path).read_text().splitlines()
    if not lines or int(lines[0]) != count or len(lines) != count + 1:
        raise InvalidArtifact('point-charge gradient inventory differs')
    values = np.array([[float(x.replace('D', 'E')) for x in s.split()] for s in lines[1:]])
    if values.shape != (count, 3) or not np.isfinite(values).all():
        raise InvalidArtifact('invalid native point-charge gradient')
    return values


def prepare(inputs, plan, output):
    config = read_json(inputs)
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    # Snapshot Python sources to protect queued research from concurrent edits.
    impl = out / 'implementation'
    impl.mkdir()
    pins = {}
    for source in sorted(Path(__file__).parent.glob('*.py')):
        target = impl / source.name
        shutil.copyfile(source, target)
        pins[source.name] = record(target)
    # Research-only allocation-aware renderer; production renderer stays intact.
    shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py')
    shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
    pins['_base_render_orca_runtime_input.py']=record(impl/'_base_render_orca_runtime_input.py')
    pins['render_orca_runtime_input.py']=record(impl/'render_orca_runtime_input.py')
    tasks = []
    for metal in ('Ca', 'La'):
        endpoint = config['endpoints'][metal]
        if endpoint['multiplicity'] != 1:
            raise InvalidArtifact('scout only supports physical singlet Ca/La')
        for environment in ('A', 'B', 'isolated'):
            td = out / f'{metal}_{environment}'
            td.mkdir()
            shutil.copyfile(verify(endpoint['xyz']), td / 'core.xyz')
            embedded = environment != 'isolated'
            inp = embedded_input(endpoint['charge']) if embedded else scientific_input(endpoint['charge'])
            (td / 'endpoint.inp').write_text(inp)
            task = dict(task_id=td.name, metal=metal, environment=environment,
                        charge=endpoint['charge'], multiplicity=1,
                        input=record(td/'endpoint.inp'), xyz=record(td/'core.xyz'),
                        output_path=str(td/'endpoint.out'), engrad_path=str(td/'endpoint.engrad'),
                        task_type='analytic_gradient')
            if embedded:
                shutil.copyfile(verify(config['environments'][environment]['pointcharges']), td/'environment.pc')
                task['pointcharges'] = record(td/'environment.pc')
            tasks.append(task)
    m = dict(protocol_id=PROTOCOL, method=METHOD, inputs=record(inputs), agreement=record(plan),
             tasks=tasks, orca=record(ORCA), implementation=pins,
             execution_policy={'task_runner':pins['run_orca_task_manifest.py'],
                               'runtime_renderer':pins['render_orca_runtime_input.py']},
             execution_resources={'mpi_ranks':16, 'concurrent_tasks':4},
             energy_scope='embedded_electronic_component_only',
             full_hybrid_status='unsupported_cross_parameters_and_boundary_reference',
             classification=None, affinity=None, compute_budget=None, wall_time_limit=None)
    for task in tasks:
        task['cache_key'] = cache_key({'task':task, 'protocol':PROTOCOL, 'method':METHOD,
                                      'inputs':m['inputs'], 'orca':m['orca'], 'implementation':pins})
    write_new(out/'manifest.json',m)
    return {'status':'prepared', 'manifest':record(out/'manifest.json'), 'task_count':len(tasks)}


def validate(manifest):
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL or m['method'] != METHOD:
        raise InvalidArtifact('unexpected reference protocol')
    config=read_json(verify(m['inputs']))
    # Recursively verify preparation references without interpreting unrelated
    # historical output paths as compatible computed endpoints.
    def check_pins(value):
        if isinstance(value,dict):
            if 'path' in value and 'sha256' in value: verify(value)
            else:
                for child in value.values(): check_pins(child)
        elif isinstance(value,list):
            for child in value: check_pins(child)
    check_pins(config)
    for pin in m['implementation'].values(): verify(pin)
    if {(t['metal'],t['environment']) for t in m['tasks']} != {
            (metal,env) for metal in ('Ca','La') for env in ('A','B','isolated')} or len(m['tasks']) != 6:
        raise InvalidArtifact('six-cell scout inventory changed')
    for t in m['tasks']:
        expected = embedded_input(t['charge']) if t['environment'] != 'isolated' else scientific_input(t['charge'])
        if verify(t['input']).read_text() != expected or t['multiplicity'] != 1:
            raise InvalidArtifact('reference Hamiltonian or state changed')
        if 'pointcharges' in t: verify(t['pointcharges'])
        bare = {k:v for k,v in t.items() if k != 'cache_key'}
        if t['cache_key'] != cache_key({'task':bare,'protocol':PROTOCOL,'method':METHOD,
                                      'inputs':m['inputs'],'orca':m['orca'],'implementation':m['implementation']}):
            raise InvalidArtifact('reference cache identity changed')
        original=config['endpoints'][t['metal']]
        if t['charge']!=original['charge'] or t['xyz']['sha256']!=original['xyz']['sha256']:
            raise InvalidArtifact('source geometry or charge changed')
        if t['environment']!='isolated' and t['pointcharges']['sha256']!=config['environments'][t['environment']]['pointcharges']['sha256']:
            raise InvalidArtifact('source field changed')
    for env in ('A','B','isolated'):
        ca,la = [next(t for t in m['tasks'] if t['metal']==metal and t['environment']==env) for metal in ('Ca','La')]
        paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
        if env != 'isolated' and ca['pointcharges']['sha256'] != la['pointcharges']['sha256']:
            raise InvalidArtifact('paired environment differs')
    return dry_run(manifest)


def collect(manifest):
    from ggr_sensitivity import executed
    validate(manifest)
    m, rows = executed(manifest)
    for t in m['tasks']:
        row = rows[t['task_id']]
        if row['status'] != 'complete': continue
        try:
            row.update(parse_endpoint(t, verify(row['output']), t['engrad_path'],
                                      permanent_field=t['environment']!='isolated'))
            if t['environment'] != 'isolated':
                path = Path(t['output_path']).with_name('endpoint.runtime.pcgrad')
                count = int(verify(t['pointcharges']).read_text().splitlines()[0])
                g = read_pcgrad(path,count)
                # Preserve raw units and output; no total physical force claim before FD.
                row.update(pointcharge_gradient=record(path), pointcharge_gradient_units='Hartree/bohr',
                           pointcharge_gradient_max_abs=float(abs(g).max()),
                           derivative_qualification='pending_directional_native_checks',
                           physical_boundary_force_status='not_projected')
        except (ValueError,OSError) as exc:
            row.update(status='invalid', reason=str(exc), energy_hartree=None)
    response = {}
    for metal in ('Ca','La'):
        a,b = (rows[f'{metal}_{env}'] for env in ('A','B'))
        response[metal] = (b['energy_hartree']-a['energy_hartree']) if all(
            r['status']=='complete' for r in (a,b)) else None
    delta = response['La']-response['Ca'] if all(v is not None for v in response.values()) else None
    return dict(protocol_id=PROTOCOL, manifest=record(manifest), rows=rows,
                status='complete' if all(r['status']=='complete' for r in rows.values()) else 'incomplete',
                per_metal_response_hartree=response, delta_env_el_hartree=delta,
                delta_env_el_kcal_mol=None if delta is None else delta*HA_TO_KCAL,
                full_hybrid_status=m['full_hybrid_status'], classification=None,
                candidate_comparison_status='unavailable_until_exact_checkpoint_and_interface',
                measured_execution_events=str(Path(manifest).parent/'budget_events.jsonl'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='op',required=True)
    a=sub.add_parser('prepare')
    for k in ('inputs','plan','output'): a.add_argument('--'+k,required=True)
    for op in ('dry-run','execute','collect'):
        a=sub.add_parser(op); a.add_argument('--manifest',required=True); a.add_argument('--output')
    args=p.parse_args()
    if args.op=='prepare': result=prepare(args.inputs,args.plan,args.output)
    elif args.op=='dry-run': result=validate(args.manifest)
    elif args.op=='collect': result=collect(args.manifest)
    else:
        validate(args.manifest)
        # Match 64 allocated task slots; the existing runner partitions four workers.
        if int(os.environ.get('SLURM_NTASKS','0')) < 16:
            raise InvalidArtifact('native MPI requires at least 16 allocated task slots')
        result=execute(args.manifest)
    if args.op!='prepare' and args.output: write_new(args.output,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
