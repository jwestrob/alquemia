"""Admit declared physical-coordinate candidates to the existing common-pool runner.

The specification pins common pilot inputs and a branch agreement. Each case has
case_id, status (prepared/unavailable), and candidates [{id, full_q, optional
coordinate, optional native_reuse}]. native_reuse maps a metal to its real saved
MACE receipt. Chemical states and source maps always come from the pinned inputs.
"""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import re
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from accommodation_torsion_profiles import geometry_check
from mace_hybrid import check_atoms, write_xyz
from mace_site_kinematics import Kinematics
import nikasha_pool as pool

PROTOCOL = 'nikasha_declared_finite_geometry_native_OMOL_GFN2_ALPB_v1'
BRANCHES = ('structure_informed_starts', 'solvent_guided', 'local_basin_breadth')


def sources(inputs, cid):
    row = next(r for r in inputs['cases'] if r['case_id'] == cid)
    collection = pool.pinned(row['pool_collection'])
    base = next(c for c in collection['cases'] if c['case_id'] == cid)
    bm = pool.pinned(row['pool_manifest'])
    if collection['manifest'] != row['pool_manifest'] or base['pool']['status'] != 'available':
        raise InvalidArtifact('complete pinned base pool required')
    proposals = pool.pinned(row['proposal_collection']); pm = pool.pinned(proposals['manifest'])
    tasks = {z: next(t for t in pm['tasks'] if (t['case_id'], t['metal']) == (cid, z)) for z in ('Ca', 'La')}
    for z, task in tasks.items():
        if xyz(verify(task['xyz'])) != xyz(verify(base['matrix'][z]['origin']['xyz'])):
            raise InvalidArtifact('source proposal origin differs from base pool')
        if task['multiplicity'] != 1: raise InvalidArtifact('finite pilot supports original singlets only')
        verify(task['mapping']); verify(task['source_preparation'])
    if tasks['La']['charge'] - tasks['Ca']['charge'] != 1:
        raise InvalidArtifact('paired charge convention differs')
    if not pool.same_geometry(xyz(verify(tasks['Ca']['xyz'])), xyz(verify(tasks['La']['xyz']))):
        raise InvalidArtifact('paired origins differ')
    if tasks['Ca']['source_preparation'] != tasks['La']['source_preparation']:
        raise InvalidArtifact('paired preparation differs')
    return row, base, bm, tasks


def physical_candidate(candidate, tasks, settings):
    if not re.fullmatch(r'[A-Za-z0-9_+-]+', candidate['id']):
        raise InvalidArtifact('unsafe candidate identifier')
    q = np.asarray(candidate['full_q'], dtype=float)
    rows = {}; checks = {}
    for z, task in tasks.items():
        kin = Kinematics(pool.pinned(task['mapping'])['context'])
        if q.shape != (len(kin.modes),) or not np.isfinite(q).all():
            raise InvalidArtifact('candidate physical coordinate shape differs')
        if any(q[i] != 0 for i in set(range(len(q))) - set(task['active_indices'])):
            raise InvalidArtifact('candidate moves an undeclared physical mode')
        if any(kin.modes[i]['unit'] != 'radian' for i in task['active_indices']):
            raise InvalidArtifact('pilot only admits existing physical angular modes')
        if np.max(np.abs(q)) > settings['maximum_angle_radian'] + 1e-12:
            raise InvalidArtifact('candidate outside declared angle domain')
        atoms = xyz(verify(task['xyz'])); symbols = [a[0] for a in atoms]
        coords = kin.evaluate(q)[1]
        checks[z] = geometry_check(kin, q, symbols)
        if not checks[z]['pass'] or checks[z]['maximum_heavy_displacement_A'] > settings['maximum_heavy_displacement_A'] + 1e-7:
            raise InvalidArtifact('candidate physical geometry/domain rejected')
        rows[z] = [(a, *map(float, p)) for a, p in zip(symbols, coords)]
        check_atoms(rows[z], task['charge'])
    if not pool.same_geometry(rows['Ca'], rows['La']):
        raise InvalidArtifact('candidate paired coordinates differ')
    if candidate.get('coordinate'):
        actual = xyz(verify(candidate['coordinate']))
        if not pool.same_geometry(rows['Ca'], actual):
            raise InvalidArtifact('declared candidate coordinate does not reproduce physical map')
        # Preserve the actually evaluated coordinates after checking their map.
        # Regenerating floating-point copies would defeat exact native reuse.
        rows = {z: [(z,*actual[0][1:]),*actual[1:]] for z in ('Ca','La')}
    return rows, checks


def load_spec(path):
    spec = read_json(path); inputs = pool.pinned(spec['inputs'])
    verify(spec['agreement']); verify(inputs['agreement'])
    if spec['branch'] not in BRANCHES: raise InvalidArtifact('unsupported pilot branch')
    expected = {'maximum_angle_radian': .8, 'maximum_heavy_displacement_A': .8}
    if spec['coordinate_limits'] != expected: raise InvalidArtifact('current approved physical domain changed')
    ids = [c['case_id'] for c in spec['cases']]
    declared = [r['case_id'] for r in inputs['cases'] if spec['branch'] != 'local_basin_breadth' or r['basin_four']]
    if ids != declared: raise InvalidArtifact('declared common pilot denominator changed')
    maximum = spec['maximum_candidates_per_case']
    if not isinstance(maximum, int) or maximum < 1: raise InvalidArtifact('finite candidate count required')
    if spec['branch'] != 'local_basin_breadth' and maximum > 4:
        raise InvalidArtifact('search pilot exceeds four new candidates per source')
    for c in spec['cases']:
        if c['status'] not in ('prepared', 'unavailable'): raise InvalidArtifact('unknown candidate status')
        candidates = c.get('candidates', [])
        if len(candidates) > maximum or len({p['id'] for p in candidates}) != len(candidates):
            raise InvalidArtifact('candidate count/identities differ')
        if c['status'] == 'prepared' and not candidates: raise InvalidArtifact('prepared source has no proposal')
    return spec, inputs


def prepare(specification, output, shards=1):
    began = time.monotonic(); spec, inputs = load_spec(specification)
    if shards not in (1, 4): raise InvalidArtifact('unsupported existing shard layout')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    cases = []; tasks = []; reference = None
    for proposed in spec['cases']:
        cid = proposed['case_id']; source, base, bm, endpoints = sources(inputs, cid)
        if reference is None: reference = bm
        if any(bm[k] != reference[k] for k in ('model','software','orca','gpu_python','cpu_python')):
            raise InvalidArtifact('pilot scientific executable/model mismatch')
        c = copy.deepcopy(base); c.pop('pool')
        c.update(base_collection=source['pool_collection'], prior_pool=base['pool'],
                 proposal_metadata=proposed, status=proposed['status'], reason=proposed.get('reason'))
        for cells in c['matrix'].values():
            for cell in cells.values(): cell['reused'] = True
        if proposed['status'] != 'prepared': cases.append(c); continue
        known = {p['id']: xyz(verify(p['xyz'])) for p in c['candidates']}
        for candidate in proposed['candidates']:
            name = candidate['id']
            if name in c['aliases'] or name in known: raise InvalidArtifact('new candidate shadows archived identity')
            paired, checks = physical_candidate(candidate, endpoints, spec['coordinate_limits'])
            match = next((k for k, a in known.items() if pool.same_geometry(a, paired['Ca'])), None)
            if match is not None:
                c['aliases'][name] = {'representative': match, 'full_q': candidate['full_q'], 'physical_checks': checks}
                continue
            match = name; known[name] = paired['Ca']
            cd = out/'candidates'/cid/name; cd.mkdir(parents=True)
            xp = cd/'geometry.xyz'; write_xyz(xp, paired['Ca'])
            c['candidates'].append({'id': name, 'xyz': record(xp), 'full_q': candidate['full_q']})
            c['aliases'][name] = {'representative': match, 'source_xyz': record(xp), 'full_q': candidate['full_q'], 'physical_checks': checks}
            for z, endpoint in endpoints.items():
                tid = cid+'__'+z+'__at_'+name; td = out/'cells'/tid; td.mkdir(parents=True)
                cp = td/'context.xyz'; write_xyz(cp, paired[z])
                task = {'task_id':tid,'case_id':cid,'metal':z,'candidate':name,'xyz':record(cp),
                        'charge':endpoint['charge'],'multiplicity':endpoint['multiplicity'],
                        'source_mapping':endpoint['mapping'],'source_preparation':endpoint['source_preparation']}
                reuse = candidate.get('native_reuse', {}).get(z)
                if reuse: task['native_reuse'] = reuse; pool.native_reuse(task, bm)
                tasks.append(task)
                c['matrix'][z][name] = {'status':'pending','task_id':tid,'xyz':task['xyz'],'reused':False}
        cases.append(c)
    impl = out/'implementation'; impl.mkdir(); pins = {}
    for path in Path(__file__).parent.glob('*.py'):
        dst = impl/path.name; shutil.copyfile(path,dst); pins[path.name] = record(dst)
    m = {k:reference[k] for k in ('model','software','orca','cpu_python','gpu_python','resources')}
    m.update(protocol_id=PROTOCOL, branch=spec['branch'], specification=record(specification),
             agreement=spec['agreement'], inputs=spec['inputs'], cases=cases, tasks=tasks,
             declared_case_ids=[c['case_id'] for c in cases], implementation=pins, shard_count=shards,
             new_MACE_cells=sum(not t.get('native_reuse') for t in tasks), maximum_new_GFN2_calls=2*len(tasks),
             new_DFT_calls=0, new_optimizations=0, production_changed=False, reference=None,
             GFN2_maxiter=reference['GFN2_maxiter'], numerical_policy_id=reference['numerical_policy_id'],
             numerical_qualification=reference['numerical_qualification'],
             source_preparation_wall_seconds=time.monotonic()-began)
    mp = out/'manifest.json'; write_new(mp,m); pool.low_prepare(mp)
    result = validate(mp); write_new(out/'PREFLIGHT.json',result); return result


def validate(manifest):
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL: raise InvalidArtifact('finite pool protocol differs')
    spec, inputs = load_spec(verify(m['specification']))
    if m['inputs'] != spec['inputs'] or m['agreement'] != spec['agreement'] or m['branch'] != spec['branch']:
        raise InvalidArtifact('finite pilot specification differs')
    for pin in m['implementation'].values(): verify(pin)
    verify(m['software']); verify(m['orca'])
    ids = [c['case_id'] for c in spec['cases']]
    if m['declared_case_ids'] != ids or [c['case_id'] for c in m['cases']] != ids:
        raise InvalidArtifact('finite pilot denominator differs')
    lookup = {t['task_id']:t for t in m['tasks']}
    if len(lookup) != len(m['tasks']): raise InvalidArtifact('duplicate finite cell')
    used = set()
    for proposed, c in zip(spec['cases'], m['cases']):
        source, base, bm, endpoints = sources(inputs, c['case_id'])
        if any(m[k] != bm[k] for k in ('model','software','orca','gpu_python','cpu_python','GFN2_maxiter','numerical_policy_id','numerical_qualification')):
            raise InvalidArtifact('finite pool changed source scientific recipe')
        if c['prior_pool'] != base['pool'] or c['proposal_metadata'] != proposed or c['status'] != proposed['status']:
            raise InvalidArtifact('prior result or declared proposal state differs')
        for z in ('Ca','La'):
            for name, old in base['matrix'][z].items():
                if c['matrix'][z][name] != {**old,'reused':True}: raise InvalidArtifact('archived cell changed')
        if proposed['status'] != 'prepared': continue
        for candidate in proposed['candidates']:
            paired, checks = physical_candidate(candidate,endpoints,spec['coordinate_limits'])
            alias = c['aliases'][candidate['id']]; name = alias['representative']
            if alias['full_q'] != candidate['full_q']: raise InvalidArtifact('alias coordinate differs')
            for z in ('Ca','La'):
                cell = c['matrix'][z][name]
                if not pool.same_geometry(paired[z],xyz(verify(cell['xyz']))):
                    raise InvalidArtifact('finite cell no longer matches physical candidate')
                if cell['reused']: continue
                t = lookup[cell['task_id']]; used.add(t['task_id'])
                if (t['case_id'],t['metal'],t['candidate'],t['charge'],t['multiplicity'],t['source_mapping'],t['source_preparation']) != (c['case_id'],z,name,endpoints[z]['charge'],1,endpoints[z]['mapping'],endpoints[z]['source_preparation']):
                    raise InvalidArtifact('finite cell source or state differs')
                atoms = xyz(verify(t['xyz']))
                if atoms[0][0] != z or not pool.same_geometry(atoms,paired[z]):
                    raise InvalidArtifact('finite endpoint metal or coordinates differ')
                check_atoms(atoms,t['charge'])
                if t.get('native_reuse'): pool.native_reuse(t,m)
    if used != set(lookup) or m['maximum_new_GFN2_calls'] != 2*len(lookup):
        raise InvalidArtifact('undeclared or missing finite calculation')
    if m['new_MACE_cells'] != sum(not t.get('native_reuse') for t in lookup.values()):
        raise InvalidArtifact('finite native-call count differs')
    return {'status':'validated','manifest':record(manifest),'denominator':len(ids),
            'prepared':sum(c['status']=='prepared' for c in m['cases']),
            'new_MACE_cells':m['new_MACE_cells'],'new_GFN2_calls':2*len(lookup),'new_DFT_calls':0}


def main():
    p = argparse.ArgumentParser(description=__doc__); s = p.add_subparsers(dest='operation',required=True)
    q=s.add_parser('prepare'); q.add_argument('--specification',required=True);q.add_argument('--output',required=True);q.add_argument('--shards',type=int,default=1)
    for name in ('validate','execute-mace','collect'):
        q=s.add_parser(name);q.add_argument('--manifest',required=True)
        if name=='collect':q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('operation').replace('-','_')
    fn={'prepare':prepare,'validate':validate,'execute_mace':pool.execute_mace,'collect':pool.collect}[op]
    print(json.dumps(fn(**a)))


if __name__ == '__main__': main()
