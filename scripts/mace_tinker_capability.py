"""Exact real-framework export and parameter-only Tinker initialization."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new
from mace_amoeba_capability import CASES, topology

PROTOCOL = 'frozen_framework_Tinker26_2_parameter_capability_v1'
STATE_CASES = ('GGR_extended', 'ALPHA_1F6S', 'ALPHA_6IP9')


def prepare(framework_result, parameters, build_receipt, states, plan, output):
    from openmm import app
    parent = read_json(framework_result)
    build = read_json(build_receipt)
    if build['status'] != 'built_not_scientifically_evaluated':
        raise InvalidArtifact('native parameter probe unavailable')
    verify(build['executable'])
    root = Path(output).resolve(); root.mkdir(parents=True, exist_ok=False)
    impl = root/'implementation'; impl.mkdir()
    for name in ('mace_tinker_capability.py','mace_amoeba_capability.py','affordable_common.py'):
        shutil.copyfile(Path(__file__).with_name(name), impl/name)
    tasks = []
    start = time.monotonic(); cpu = time.process_time()
    ffpin = read_json(verify(parent['manifest']))['forcefields'][0]
    ff = app.ForceField(str(verify(ffpin)))
    for row, case, state_case in zip(parent['cases'], CASES, STATE_CASES):
        if row['case_id'] != case or row['status'] != 'framework_parameterized_no_energy':
            raise InvalidArtifact('declared framework does not match')
        prep = read_json(verify(row['source_preparation']))
        mapping = read_json(verify(row['mapping']))
        top, positions, ids, bonds, water_bonds, metals = topology(prep)
        if ids != mapping['system_atom_ids']:
            raise InvalidArtifact('framework atom order changed')
        np.testing.assert_array_equal(positions, mapping['positions_A'])
        bonded = ff._buildBondedToAtomList(top)
        types = []
        for residue in top.residues():
            template, matches = ff._getResidueTemplateMatches(residue, bonded)
            if matches is None:
                raise InvalidArtifact('unmatched real residue')
            types.extend(int(template.atoms[i].type) for i in matches)
        if len(types) != len(ids):
            raise InvalidArtifact('atom type count differs')
        state_path = Path(states)/state_case/'state.json'
        state = read_json(state_path)
        if state['physical_atoms'] != prep['physical_atoms']:
            raise InvalidArtifact('QM support has a different physical system')
        support = set(state['projection_support_ids'])
        if not support <= set(ids)|{'metal'}:
            raise InvalidArtifact('unsupported QM source identity')
        frozen = [i+1 for i, pid in enumerate(ids) if pid in support]
        allowed = [i+1 for i, pid in enumerate(ids) if pid not in support]
        if not frozen or not allowed:
            raise InvalidArtifact('missing frozen/environment source set')
        directory = root/case; directory.mkdir()
        neighbors = {pid: [] for pid in ids}; index = {pid:i+1 for i,pid in enumerate(ids)}
        for a,b in bonds+water_bonds:
            neighbors[a].append(index[b]); neighbors[b].append(index[a])
        lines = [f'{len(ids)}  {case} frozen physical protein/water framework; metal separately unparameterized']
        for i, (atom, pid, coords, kind) in enumerate(zip(top.atoms(), ids, positions, types), 1):
            lines.append(f'{i} {atom.element.symbol} '+ ' '.join(f'{x:.17g}' for x in coords)
                         + f' {kind} '+ ' '.join(map(str, sorted(neighbors[pid]))))
        (directory/'framework.xyz').write_text('\n'.join(lines)+'\n')
        key = ['parameters '+str(Path(parameters).resolve()), 'polarization MUTUAL',
               'polar-eps 0.00001', 'solvate GK']
        # Positive indices only, with short lines; these are actual environmental
        # atoms. Native initialization and the explicit adapter mask are both audited.
        key.extend('polarizable '+' '.join(map(str, allowed[i:i+12])) for i in range(0,len(allowed),12))
        (directory/'framework.key').write_text('\n'.join(key)+'\n')
        (directory/'framework.freeze').write_text(str(len(frozen))+'\n'+'\n'.join(map(str,frozen))+'\n')
        meta = dict(case_id=case, system_atom_ids=ids, physical_atoms=prep['physical_atoms'],
                    tinker_types=types, frozen_indices=frozen, allowed_indices=allowed,
                    source_mapping=row['mapping'], source_preparation=row['source_preparation'],
                    support_state=record(state_path), unparameterized_sources=metals)
        write_new(directory/'mapping.json', meta)
        task = dict(case_id=case, xyz=record(directory/'framework.xyz'), key=record(directory/'framework.key'),
                    mask=record(directory/'framework.freeze'), mapping=record(directory/'mapping.json'))
        task['cache_key'] = cache_key(dict(task=task, parameters=record(parameters), executable=build['executable'], protocol=PROTOCOL))
        tasks.append(task)
    manifest = dict(protocol=PROTOCOL, plan=record(plan), parent=record(framework_result),
        build_receipt=record(build_receipt), executable=build['executable'], parameters=record(parameters),
        forcefield=ffpin, implementation=[record(p) for p in sorted(impl.glob('*.py'))], tasks=tasks,
        new_energy_calls=0, new_force_calls=0, full_model_qualified=False,
        preparation_receipt=dict(wall_seconds=time.monotonic()-start, process_cpu_seconds=time.process_time()-cpu,
            peak_process_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    write_new(root/'manifest.json', manifest)
    return manifest


def validate(path):
    m = read_json(path)
    for pin in [m['plan'],m['parent'],m['build_receipt'],m['executable'],m['parameters'],m['forcefield']]+m['implementation']:
        verify(pin)
    for task in m['tasks']:
        body = {k:v for k,v in task.items() if k!='cache_key'}
        if task['cache_key'] != cache_key(dict(task=body,parameters=m['parameters'],executable=m['executable'],protocol=m['protocol'])):
            raise InvalidArtifact('changed task configuration')
        for k in ('xyz','key','mask','mapping'):
            verify(task[k])
    return m


def inspect_output(text, mapping):
    records = {k: {} for k in ('ATOM','AXIS','POLE','BONDS','FROZEN_BEFORE','FROZEN_AFTER')}
    complete = False
    for line in text.splitlines():
        fields = line.split()
        if fields and fields[0] in records:
            kind, index = fields[0], int(fields[1])
            if index in records[kind]:
                raise InvalidArtifact('duplicate native atom row')
            records[kind][index] = fields[2:]
        if fields and fields[0]=='ALQUEMIA_PARAMETER_ONLY_COMPLETE':
            complete = [int(v) for v in fields[1:]] == [len(mapping['system_atom_ids']),len(mapping['frozen_indices'])]
    if not complete:
        raise InvalidArtifact('native parameter-only completion missing')
    old = read_json(verify(mapping['source_mapping']))
    n = len(mapping['system_atom_ids']); expected = set(range(1,n+1))
    for k in ('ATOM','AXIS','POLE','BONDS'):
        if set(records[k]) != expected:raise InvalidArtifact('missing physical parameter rows')
    differences = dict(max_coordinate_error_A=0.,max_charge_error_e=0.,max_polarizability_error_A3=0.,
                       max_damping_error_A_half=0.,max_local_multipole_error_MD_units=0.,axis_mismatches=0)
    by_id = {pid:i+1 for i,pid in enumerate(mapping['system_atom_ids'])}
    bondset = {tuple(sorted((by_id[a],by_id[b]))) for a,b in old['protein_bonds']+old['retained_water_bonds']}
    native_bonds = set()
    for i, params in enumerate(old['parameters'], 1):
        row = records['ATOM'][i]; mp = params['multipole_md_units']
        if int(row[0]) != mapping['tinker_types'][i-1]:raise InvalidArtifact('native atom type differs')
        diffs = np.abs(np.array(list(map(float,row[2:5])))-old['positions_A'][i-1])
        differences['max_coordinate_error_A'] = max(differences['max_coordinate_error_A'],float(max(diffs)))
        pole = list(map(float,records['POLE'][i])); converted = [pole[0]]+[v*.1 for v in pole[1:4]]+[v*.01 for v in pole[4:13]]
        target = [mp[0]]+mp[1]+mp[2]
        differences['max_local_multipole_error_MD_units'] = max(differences['max_local_multipole_error_MD_units'],max(abs(a-b) for a,b in zip(converted,target)))
        differences['max_charge_error_e'] = max(differences['max_charge_error_e'],abs(pole[0]-mp[0]))
        differences['max_polarizability_error_A3'] = max(differences['max_polarizability_error_A3'],abs(float(row[5])-mp[9]*1000))
        differences['max_damping_error_A_half'] = max(differences['max_damping_error_A_half'],abs(float(row[6])-mp[8]*np.sqrt(10)))
        axis = records['AXIS'][i]
        expected_axes = [v+1 if v>=0 else 0 for v in mp[4:7]]
        if [abs(int(v)) for v in axis[1:]] != expected_axes:differences['axis_mismatches'] += 1
        b = list(map(int,records['BONDS'][i]))
        if b[0]!=len(b)-1:raise InvalidArtifact('native bond row count')
        native_bonds.update(tuple(sorted((i,j))) for j in b[1:])
    if native_bonds != bondset:raise InvalidArtifact('native connectivity differs')
    frozen = set(mapping['frozen_indices'])
    if any(set(records[k])!=frozen for k in ('FROZEN_BEFORE','FROZEN_AFTER')):
        raise InvalidArtifact('source-freeze rows differ')
    reenabled = 0
    for i in frozen:
        before, after = records['FROZEN_BEFORE'][i], records['FROZEN_AFTER'][i]
        if after[0]!='F' or before[1:]!=after[1:]:
            raise InvalidArtifact('freeze adapter changed damping/parameters or failed')
        reenabled += before[0]=='T'
    # No tolerance is a prediction/energy gate: these compare serialized input
    # parameters in explicitly stated units. Differences remain fully visible.
    return dict(status='native_parameters_inspected', differences=differences,
        coordinates_exact=differences['max_coordinate_error_A']==0., connectivity_exact=True,
        native_keyword_reenabled_sources=reenabled, requested_frozen_sources=len(frozen),
        adapter_mask_pass=True, energy=None, full_model_qualified=False, native_records=records)


def execute(manifest_path):
    m = validate(manifest_path); root = Path(manifest_path).resolve().parent
    results = []
    for task in m['tasks']:
        directory = root/task['case_id']; target = directory/'native_result.json'
        if target.exists():
            result = read_json(target)
            if result['cache_key']!=task['cache_key']:raise InvalidArtifact('incompatible partial-job result')
            results.append(result); continue
        before=time.monotonic(); cpu=resource.getrusage(resource.RUSAGE_CHILDREN)
        env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
        cmd=[str(verify(m['executable'])),str(verify(task['xyz'])),str(verify(task['mask']))]
        with (directory/'native.log').open('x') as out:
            p=subprocess.run(cmd,cwd=directory,env=env,stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT)
        result=dict(status='failed',case_id=task['case_id'],cache_key=task['cache_key'],command=cmd,returncode=p.returncode)
        try:
            if p.returncode:raise InvalidArtifact('native parameter probe failed')
            result.update(inspect_output((directory/'native.log').read_text(),read_json(verify(task['mapping']))))
        except Exception as error:result['failure_reason']=str(error)
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        result['receipt']=dict(wall_seconds=time.monotonic()-before,child_CPU_seconds=after.ru_utime+after.ru_stime-cpu.ru_utime-cpu.ru_stime,peak_child_RSS_KiB=after.ru_maxrss,new_energy_calls=0,new_force_calls=0)
        result['native_log']=record(directory/'native.log')
        write_new(target,result); results.append(result)
    write_new(root/'result.json',dict(manifest=record(manifest_path),cases=results,full_model_qualified=False,numerical_score=None))
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    prep=sub.add_parser('prepare')
    for name in ('framework-result','parameters','build-receipt','states','plan','output'):prep.add_argument('--'+name,required=True)
    for name in ('dry-run','execute'):
        s=sub.add_parser(name);s.add_argument('--manifest',required=True)
    a=p.parse_args()
    if a.command=='prepare':
        m=prepare(a.framework_result,a.parameters,a.build_receipt,a.states,a.plan,a.output);print(json.dumps({'prepared':len(m['tasks']),'energy_calls':0}))
    elif a.command=='dry-run':
        m=validate(a.manifest);print(json.dumps({'validated':len(m['tasks']),'energy_calls':0}))
    else:
        results=execute(a.manifest);print(json.dumps([{k:r.get(k) for k in ('case_id','status','differences','native_keyword_reenabled_sources','failure_reason')} for r in results]))


if __name__=='__main__':main()
