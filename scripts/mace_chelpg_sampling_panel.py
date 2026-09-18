"""Declared same-density sampling extension; no new SCF or biological score."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
import mace_chelpg_sampling as replay

PROTOCOL = replay.PROTOCOL + '_five_geometry_panel_v1'
CASES = ['GGR_extended', 'GGR_2FW0', 'GGR_2FVY', 'ALPHA_1F6S', 'ALPHA_6IP9']
SETTINGS = {'original': {'grid_A': .3, 'extent_A': 2.8},
            'finer': {'grid_A': .2, 'extent_A': 2.8},
            'finest': {'grid_A': .15, 'extent_A': 2.8},
            'extent': {'grid_A': .3, 'extent_A': 3.5}}


def parse_sampling(text, elements, charge, sampling):
    """Explicit variant parser; the released default parser stays strict."""
    for pattern, expected in [(r'Grid spacing\s+\.{3}\s+([\d.]+)', sampling['grid_A']),
                              (r'Point Cut-Off\s+\.{3}\s+([\d.]+)', sampling['extent_A'])]:
        values = re.findall(pattern, text)
        if len(values) != 1 or float(values[0]) != expected:
            raise InvalidArtifact('declared native sampling control differs')
    if not re.search(r'Van-der-Waals Radii\s+\.{3}\s+COSMO', text) or not re.search(r'Dipole moment constraint\s+\.{3}\s+FALSE', text):
        raise InvalidArtifact('CHELPG radii/dipole convention differs')
    if 'CHELPG charges calculated...' not in text or text.count('CHELPG Charges') != 1:
        raise InvalidArtifact('CHELPG calculation incomplete')
    values = re.findall(r'^\s*(\d+)\s+(\w+)\s+:\s+([-+\d.eE]+)\s*$', text, re.M)
    if [(int(i), e) for i, e, q in values] != list(enumerate(elements)):
        raise InvalidArtifact('CHELPG atom order/count differs')
    q = np.array([float(v) for i, e, v in values]); total = re.findall(r'Total charge:\s+([-+\d.eE]+)', text)
    if not np.isfinite(q).all() or abs(math.fsum(q)-charge) > 5e-5 or len(total) != 1 or abs(float(total[0])-charge) > 1e-6:
        raise InvalidArtifact('CHELPG formal/ECP charge mismatch')
    return q


def source_inventory(c):
    old = read_json(verify(c['old_audit'])); obs = read_json(verify(c['observations_manifest']))
    new = read_json(verify(c['density_report'])); qr = read_json(verify(obs['quantum'])); qm = read_json(verify(qr['manifest']))
    oq = read_json(verify(old['quantum'])); oc = read_json(verify(old['charges']))
    if not old['complete'] or not old['gates_pass'] or new['manifest'] != c['observations_manifest'] or new['status'] != 'complete' or not new['projection_gate_pass']:
        raise InvalidArtifact('qualified source observations required')
    result = {}
    for case in CASES:
        for metal in ('Ca', 'La'):
            tid = case+'_'+metal
            if tid in old['rows']:
                src = old['rows'][tid]['source_task']; qt = old['rows'][tid]['quantum_task']; qrow = oq['rows'][tid]; drow = oc['rows'][tid]
                group = read_json(verify(src['probe_groups'])); exterior = group['exterior_indices']
            else:
                src = next(t for t in obs['tasks'] if t['task_id'] == tid)
                qt = next(t for t in qm['tasks'] if t['task_id'] == tid); qrow = qr['rows'][tid]; drow = new['rows'][tid]
                group = read_json(verify(src['probes'])); exterior = list(range(group['exterior_start'], group['exterior_start']+group['exterior_count']))
            for pin in [qt['input'], qt['xyz'], qt['pointcharges'], src['files']['gbw'], qrow['receipt'], qrow['output'], drow['potential'], src['points'], src['projection']]: verify(pin)
            if qrow['energy_hartree'] is None or qt['multiplicity'] != 1: raise InvalidArtifact('source state unsupported')
            result[tid] = dict(quantum_task=qt, observation_task=src, density_row=drow, quantum_row=qrow, exterior_indices=exterior)
    return result, qm, obs


def qualified(c):
    pop = read_json(verify(c['qualification_population'])); identity = read_json(verify(c['qualification_identity']))
    im = replay.validate_identity(verify(identity['manifest']))
    if im['population'] != c['qualification_population'] or not identity['complete'] or not identity['qualification_pass']:
        raise InvalidArtifact('actual default density identity not qualified')
    mp = verify(pop['manifest']); replay.validate(mp)
    from density_embedding import parse_potential
    for t in read_json(mp)['tasks']:
        if replay.native_replay(t, mp) != pop['rows'][t['task_id']]: raise InvalidArtifact('qualification population changed')
        it = next(x for x in im['tasks'] if x['task_id'] == t['task_id']); row = identity['rows'][t['task_id']]
        verify(row['receipt']['log']); verify(row['receipt']['resource_usage'])
        values = parse_potential(verify(row['potential']), np.loadtxt(verify(it['points']), skiprows=1))
        error = float(np.max(abs(values-np.array(read_json(verify(it['expected']))['potential_au']))))
        if row['receipt']['returncode'] or error > 1e-8 or error != row['maximum_identity_error_au']:
            raise InvalidArtifact('qualification potential identity changed')
    return pop, identity


def prepare(config, output):
    c = read_json(config); verify(c['plan']); pop, identity = qualified(c)
    if c['protocol'] != PROTOCOL or c['cases'] != CASES or c['settings'] != SETTINGS:
        raise InvalidArtifact('declared sampling panel differs')
    sources, qm, obs = source_inventory(c)
    root = Path(output).resolve(); root.mkdir(parents=True, exist_ok=False); impl = root/'implementation'; impl.mkdir()
    for n, pin in qm['implementation'].items(): shutil.copyfile(verify(pin), impl/n)
    for n, pin in obs['implementation'].items():
        if not (impl/n).exists(): shutil.copyfile(verify(pin), impl/n)
    for name in ('mace_chelpg_sampling.py', Path(__file__).name): shutil.copyfile(Path(__file__).with_name(name), impl/name)
    pins = {p.name: record(p) for p in impl.glob('*.py')}
    m = dict(protocol_id=PROTOCOL, method=replay.METHOD, config=record(config), agreement=c['plan'], implementation=pins,
             orca=qm['orca'], utility=obs['utilities']['orca_vpot'], tasks=[], reused=[], phase='sampling',
             execution_policy={k: pins[n] for k, n in [('task_runner', 'run_orca_task_manifest.py'), ('runtime_renderer', 'render_orca_runtime_input.py')]},
             planned_SCF_optimizations=0, requested_property_starts=38, total_property_starts_including_qualification=40,
             new_score=None, baseline_changed=False, compute_budget=None, wall_time_limit=None)
    for tid, s in sources.items():
        qt=s['quantum_task']; src=s['observation_task']; qr=s['quantum_row']
        for name, setting in SETTINGS.items():
            task_id=tid+'_'+name
            if task_id in pop['rows']:
                original=next(t for t in read_json(verify(pop['manifest']))['tasks'] if t['task_id']==task_id)
                if original['source_gbw']!=src['files']['gbw'] or original['source_quantum_receipt']!=qr['receipt'] or original['sampling']!=setting:
                    raise InvalidArtifact('default reuse source differs')
                m['reused'].append(dict(task_id=task_id, manifest=pop['manifest'], population=c['qualification_population'], identity=c['qualification_identity']))
                continue
            d=root/task_id; d.mkdir()
            for pin, filename in [(qt['xyz'], 'core.xyz'), (qt['pointcharges'], 'environment.pc'), (src['files']['gbw'], 'seed.gbw')]: shutil.copyfile(verify(pin), d/filename)
            (d/'endpoint.inp').write_text(replay.scientific_input(qt['charge'], setting['grid_A'], setting['extent_A']))
            t=dict(task_id=task_id, case_id=qt['case_id'], metal=qt['metal'], charge=qt['charge'], multiplicity=1,
                   task_type='fixed_orbital_population', setting=name, sampling=setting, input=record(d/'endpoint.inp'), xyz=record(d/'core.xyz'),
                   pointcharges=record(d/'environment.pc'), seed_gbw=record(d/'seed.gbw'), source_gbw=src['files']['gbw'],
                   source_quantum_receipt=qr['receipt'], source_quantum_output=qr['output'], source_energy_hartree=qr['energy_hartree'],
                   source_density_row=s['density_row'], source_observation_task=src, exterior_indices=s['exterior_indices'], output_path=str(d/'endpoint.out'))
            t['cache_key']=replay.key(t,m); m['tasks'].append(t)
    write_new(root/'manifest.json',m); return validate(root/'manifest.json')


def validate(path):
    from affordable_workflow import dry_run
    m=read_json(path); c=read_json(verify(m['config'])); qualified(c); sources, qm, obs=source_inventory(c)
    if m['protocol_id']!=PROTOCOL or m['method']!=replay.METHOD or c['settings']!=SETTINGS or c['cases']!=CASES or len(m['tasks'])!=38 or len(m['reused'])!=2 or m['planned_SCF_optimizations']:
        raise InvalidArtifact('declared sampling inventory differs')
    for pin in [m['agreement'],m['orca'],m['utility'],*m['implementation'].values()]: verify(pin)
    expected={case+'_'+metal+'_'+name for case in CASES for metal in ('Ca','La') for name in SETTINGS}
    if {t['task_id'] for t in m['tasks']+m['reused']}!=expected: raise InvalidArtifact('sampling cases missing or duplicated')
    if m['orca']!=qm['orca'] or m['utility']!=obs['utilities']['orca_vpot']: raise InvalidArtifact('native executables differ')
    for t in m['tasks']:
        s=sources[t['case_id']+'_'+t['metal']]; qt=s['quantum_task']; qr=s['quantum_row']; setting=SETTINGS[t['setting']]
        if t['cache_key']!=replay.key(t,m) or t['sampling']!=setting or verify(t['input']).read_text()!=replay.scientific_input(t['charge'],setting['grid_A'],setting['extent_A']): raise InvalidArtifact('changed sampling or method')
        if t['charge']!=qt['charge'] or t['multiplicity']!=qt['multiplicity'] or t['source_energy_hartree']!=qr['energy_hartree'] or t['source_quantum_receipt']!=qr['receipt'] or t['source_density_row']!=s['density_row'] or t['source_observation_task']!=s['observation_task'] or t['exterior_indices']!=s['exterior_indices']: raise InvalidArtifact('source state differs')
        for name in ('xyz','pointcharges'):
            if verify(t[name]).read_bytes()!=verify(qt[name]).read_bytes(): raise InvalidArtifact('geometry or generating field changed')
        if verify(t['seed_gbw']).read_bytes()!=verify(s['observation_task']['files']['gbw']).read_bytes(): raise InvalidArtifact('imported orbitals changed')
    return dry_run(path)


def population_row(t,path):
    row=replay.native_replay(t,path,lambda text,elements,charge:parse_sampling(text,elements,charge,t['sampling']))
    error=row.pop('default_charge_max_error_e'); default=t['setting']=='original'
    row.update(charge_change_from_original_max_e=error, default_charge_identity_pass=(error<=2e-6 if default else None),
               preliminary_identity_pass=abs(row['energy_identity_error_hartree'])<=1e-7 and (not default or error<=2e-6),
               qualification_status='requires_actual_saved_density_potential_comparison')
    return row


def collect(path,output):
    """Collect retained native receipts without starting a job or property call."""
    validate(path);m=read_json(path);rows={}
    for t in m['tasks']:
        try:rows[t['task_id']]=population_row(t,path)
        except (ValueError,OSError,KeyError) as exc:rows[t['task_id']]=dict(status='invalid_or_unavailable',reason=str(exc),preliminary_identity_pass=False)
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,root/Path(__file__).name)
    r=dict(protocol=PROTOCOL,manifest=record(path),collection_implementation=record(root/Path(__file__).name),rows=rows,
           new_native_calls=0,complete=all(v['status']=='native_NoIter_population_completed' for v in rows.values()),
           preliminary_identity_pass=all(v['preliminary_identity_pass'] for v in rows.values()),new_score=None)
    write_new(root/'result.json',r);return r


def execute(path,output):
    """Reuse the native task runner; retain its SCF failure separately from NoIter."""
    from run_orca_task_manifest import run_manifest
    validate(path); m=read_json(path); root=Path(path).resolve().parent
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<64: raise InvalidArtifact('64CPU declared allocation required')
    with (root/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB); pending=[]; rows={}
        for t in m['tasks']:
            out=Path(t['output_path']); receipt=Path(str(out)+'.execution.json')
            if out.exists() or receipt.exists():
                if not(out.exists() and receipt.exists()): raise InvalidArtifact('partial attempt retained; needs explicit fresh retry')
                rows[t['task_id']]=population_row(t,path)
            else: pending.append(t['task_id'])
        start=time.monotonic(); generic=None; error=None
        if pending:
            try: generic=run_manifest(Path(path),orca_path=verify(m['orca']),workers=4,nprocs=16,selected_tasks=set(pending),expected_runner_sha256=m['execution_policy']['task_runner']['sha256'])
            except ValueError as exc: error=str(exc)
        for t in m['tasks']:
            try: rows[t['task_id']]=population_row(t,path)
            except (ValueError,OSError,KeyError) as exc: rows[t['task_id']]=dict(status='invalid_or_unavailable',reason=str(exc),preliminary_identity_pass=False)
        r=dict(protocol=PROTOCOL,manifest=record(path),rows=rows,generic_runner_result=generic,generic_runner_error=error,
               admitted_property_starts=len(pending),SCF_optimizations=0,slurm_job_id=os.environ['SLURM_JOB_ID'],wall_seconds=time.monotonic()-start,
               complete=all(v['status']=='native_NoIter_population_completed' for v in rows.values()),
               preliminary_identity_pass=all(v['preliminary_identity_pass'] for v in rows.values()),new_score=None)
        write_new(output,r)
    return r


def prepare_identity(population,output):
    from density_embedding import parse_potential
    pop=read_json(population); mp=verify(pop['manifest']); validate(mp); parent=read_json(mp)
    if not pop['complete'] or not pop['preliminary_identity_pass']: raise InvalidArtifact('population identity gate failed')
    root=Path(output).resolve(); root.mkdir(parents=True,exist_ok=False); impl=root/'implementation'; impl.mkdir()
    for name,pin in parent['implementation'].items(): shutil.copyfile(verify(pin),impl/name)
    shutil.copyfile(__file__,impl/Path(__file__).name)
    m=dict(protocol=PROTOCOL+'_density_identity',population=record(population),plan=parent['agreement'],utility=parent['utility'],
           implementation={p.name:record(p) for p in impl.glob('*.py')},tasks=[],potential_tolerance_au=1e-8,workers=8,requested_queries=38,requested_DFT_calls=0)
    for t in parent['tasks']:
        row=population_row(t,mp)
        if row!=pop['rows'][t['task_id']]: raise InvalidArtifact('population record changed')
        src=t['source_observation_task']; coords=np.loadtxt(verify(src['points']),skiprows=1); values=parse_potential(verify(t['source_density_row']['potential']),coords)
        ix=np.array(t['exterior_indices']); points=coords[ix]; expected=values[ix]
        d=root/'tasks'/t['task_id']; d.mkdir(parents=True); pp=d/'points_bohr.xyz'
        pp.write_text(str(len(points))+'\n'+''.join(' '.join(format(float(v),'.14f') for v in point)+'\n' for point in points))
        write_new(d/'expected.json',dict(potential_au=expected.tolist(),source_potential=t['source_density_row']['potential'],source_points=src['points'],exterior_indices=t['exterior_indices']))
        task=dict(task_id=t['task_id'],files=row['files'],points=record(pp),expected=record(d/'expected.json'),source_replay_receipt=row['execution_receipt'])
        task['cache_key']=cache_key(dict(task=task,population=m['population'],utility=m['utility'],implementation=m['implementation'],tolerance=m['potential_tolerance_au']));m['tasks'].append(task)
    write_new(root/'manifest.json',m);return validate_identity(root/'manifest.json')


def validate_identity(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL+'_density_identity' or len(m['tasks'])!=38 or m['requested_DFT_calls'] or m['potential_tolerance_au']!=1e-8: raise InvalidArtifact('identity query scope changed')
    for pin in [m['population'],m['plan'],m['utility'],*m['implementation'].values()]:verify(pin)
    for t in m['tasks']:
        original={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key(dict(task=original,population=m['population'],utility=m['utility'],implementation=m['implementation'],tolerance=m['potential_tolerance_au'])):raise InvalidArtifact('query cache changed')
        for pin in [t['points'],t['expected'],t['source_replay_receipt'],*t['files'].values()]:verify(pin)
        e=read_json(verify(t['expected']));verify(e['source_potential']);verify(e['source_points'])
    return m


def execute_identity(path,output):
    from affordable_solver import run_command
    from density_embedding import parse_potential
    m=validate_identity(path);root=Path(path).resolve().parent
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<m['workers']:raise InvalidArtifact('declared utility allocation required')
    def work(t):
        d=root/'tasks'/t['task_id']; attempt=d/'attempt_0001'
        if attempt.exists():
            saved=attempt/'result.json'
            if not saved.exists():raise InvalidArtifact('partial query retained; explicit recovery needed')
            r=read_json(saved)
            if r['cache_key']!=t['cache_key'] or r['status']!='computed':raise InvalidArtifact('failed query retained; explicit retry needed')
            for pin in (r['potential'],r['receipt']['log'],r['receipt']['resource_usage']):verify(pin)
            return r
        attempt.mkdir();r=dict(task_id=t['task_id'],cache_key=t['cache_key'],status='failed',receipt=None)
        try:
            for pin in t['files'].values():shutil.copyfile(verify(pin),attempt/Path(pin['path']).name)
            gbw=attempt/Path(t['files']['gbw']['path']).name; potential=attempt/'potential.out'
            r['receipt']=run_command([str(verify(m['utility'])),str(gbw),'endpoint.runtime.scfp',str(verify(t['points'])),str(potential)],attempt,attempt/'vpot.log',attempt/'resources.txt')
            if r['receipt']['returncode']:raise InvalidArtifact('density query failed')
            values=parse_potential(potential,np.loadtxt(verify(t['points']),skiprows=1));expected=np.array(read_json(verify(t['expected']))['potential_au']);err=float(np.max(abs(values-expected)))
            r.update(status='computed',potential=record(potential),maximum_identity_error_au=err,pass_=err<=1e-8)
        except Exception as exc:r.update(failure_reason=str(exc),pass_=False)
        write_new(attempt/'result.json',r);return r
    with (root/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with ThreadPoolExecutor(max_workers=m['workers']) as pool: rows=list(pool.map(work,m['tasks']))
    r=dict(protocol=m['protocol'],manifest=record(path),rows={r['task_id']:r for r in rows},complete=all(r['status']=='computed' for r in rows),identity_pass=all(r['pass_'] for r in rows),new_score=None)
    write_new(output,r);return r


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    q=sub.add_parser('prepare-identity');q.add_argument('--population',required=True);q.add_argument('--output',required=True)
    for op in ('dry-run','execute','collect','dry-run-identity','execute-identity'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op.startswith('execute') or op=='collect':q.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.config,a.output)
    elif a.command=='prepare-identity':r=prepare_identity(a.population,a.output);r={'status':'prepared','tasks':len(r['tasks'])}
    elif a.command=='dry-run':r=validate(a.manifest)
    elif a.command=='dry-run-identity':r=validate_identity(a.manifest);r={'status':'pass','tasks':len(r['tasks'])}
    elif a.command=='execute-identity':r=execute_identity(a.manifest,a.output)
    elif a.command=='collect':r=collect(a.manifest,a.output)
    else:r=execute(a.manifest,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','generic_runner_result')}))

if __name__=='__main__':main()
