#!/usr/bin/env python3
"""Slurm-only isolated native component build and finite saved-charge replay."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import socket
import subprocess
import tarfile
import time

TERMS = ('Gelec', 'Ghb', 'Gsasa', 'Gshift')


def read(p):
    return json.loads(Path(p).read_text())


def pin(p):
    p = Path(p).resolve()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return {'path': str(p), 'sha256': h.hexdigest()}


def verify(p):
    if pin(p['path']) != p:
        raise ValueError('changed pin: '+p['path'])
    return Path(p['path'])


def dump(p, obj):
    with Path(p).open('x') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')


def command(args, cwd, stem, env):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.monotonic()
    with (cwd/(stem+'.out')).open('x') as out, (cwd/(stem+'.err')).open('x') as err:
        result = subprocess.run(args, cwd=cwd, env=env, stdout=out, stderr=err)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    row = {'command': args, 'cwd': str(cwd), 'returncode': result.returncode,
           'wall_seconds': time.monotonic()-start,
           'child_CPU_seconds': after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
           'cumulative_child_peak_RSS_KiB': after.ru_maxrss,
           'stdout': pin(cwd/(stem+'.out')), 'stderr': pin(cwd/(stem+'.err'))}
    dump(cwd/(stem+'_EXECUTION.json'), row)
    if result.returncode:
        raise RuntimeError('failed command: '+str(args)+'; see '+str(cwd/(stem+'_EXECUTION.json')))
    return row


def preflight(m, mp):
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('actual Slurm allocation required')
    if not 1 <= int(os.environ.get('SLURM_NTASKS', '0')) <= 12:
        raise ValueError('requires one to twelve allocated task CPUs')
    for p in m['implementation'].values():
        verify(p)
    if pin(__file__) != m['implementation']['REPLAY.py']:
        raise ValueError('unfrozen runner')
    for p in m['provenance'].values():
        verify(p)
    for p in m['build']['tools'].values():
        verify(p)
    for p in m['build']['libraries'].values():
        verify(p)
    for source in m['build']['sources'].values():
        verify(source['archive'])
    if m['protocol'] != 'xtb671_native_fixed_charge_alpb_replay_v1' or len(m['cases']) != 4:
        raise ValueError('changed finite scope')
    seen = set()
    for c in m['cases']:
        if c['case_id'] in seen:
            raise ValueError('duplicate case')
        seen.add(c['case_id'])
        verify(c['state_signature'])
        xyz = verify(c['xyz'])
        lines = xyz.read_text().splitlines()
        if int(lines[0]) != 192 or len(lines[2:]) != 192:
            raise ValueError('changed atom inventory')
        atoms = [line.split() for line in lines[2:]]
        if any(len(a) != 4 or not all(math.isfinite(float(v)) for v in a[1:]) for a in atoms):
            raise ValueError('invalid exact saved geometry')
        if atoms[0][0] != c['metal'] or any(a[0] not in ('H','C','N','O','S',c['metal']) for a in atoms):
            raise ValueError('unexpected atom identity/order')
        for medium in ('water','vacuum'):
            e = c['endpoints'][medium]
            point = read(verify(e['point']))
            receipt = read(verify(e['execution']))
            stdout = verify(e['stdout']).read_text()
            local_xyz = verify(e['executed_xyz'])
            qs = verify(e['charges']).read_text().split()
            if local_xyz.read_bytes() != xyz.read_bytes() or receipt['returncode'] != 0:
                raise ValueError('executed geometry or completion differs')
            if len(qs) != 192 or any(not re.fullmatch(r'[-+]?\d+\.\d{8}', q) for q in qs):
                raise ValueError('saved eight-decimal charge format differs')
            q = [float(v) for v in qs]
            if q != point['backends']['xtb']['charges_e'] or point['backends']['xtb']['status'] != 'available':
                raise ValueError('charge array does not match accepted endpoint')
            if abs(math.fsum(q)-c['charge']) > 192*m['tolerances']['charge_rounding_e']+1e-12:
                raise ValueError('rounded charges do not close to prescribed state')
            t = point['task']
            if any(t[k] != c[k] for k in ('metal','charge','multiplicity','valence_electrons','atom_count')) or t['xyz'] != c['xyz'] or t['medium'] != medium:
                raise ValueError('source state or medium differs')
            if medium == 'water':
                for expected in m['required_setup_lines']:
                    if expected not in stdout:
                        raise ValueError('archived setup differs: '+expected)
    expected = {f'sample{s}_{metal}' for s in (1,3) for metal in ('Ca','La')}
    if seen != expected:
        raise ValueError('finite source IDs differ')
    return {'status':'passed', 'manifest':pin(mp), 'slurm_job_id':os.environ['SLURM_JOB_ID'],
            'hostname':socket.gethostname(), 'affinity':sorted(os.sched_getaffinity(0))}


def build(m, out, env):
    source = out/'sources'
    source.mkdir()
    for name, record in m['build']['sources'].items():
        with tarfile.open(verify(record['archive'])) as archive:
            for member in archive.getmembers():
                destination = (source/member.name).resolve()
                if not destination.is_relative_to(source.resolve()) or member.issym() or member.islnk():
                    raise ValueError('unsupported source archive member')
            archive.extractall(source)
    wrapper = out/'wrapper'
    wrapper.mkdir()
    for name in ('CMakeLists.txt','REPLAY.f90'):
        shutil.copyfile(verify(m['implementation'][name]), wrapper/name)
    tools = {k:str(verify(v)) for k,v in m['build']['tools'].items()}
    libs = {k:str(verify(v)) for k,v in m['build']['libraries'].items()}
    directory = out/'build'
    args = [tools['cmake'],'-S',str(wrapper),'-B',str(directory),'-G','Unix Makefiles',
            '-DCMAKE_BUILD_TYPE=Release','-DCMAKE_Fortran_COMPILER='+tools['gfortran'],
            '-DCMAKE_C_COMPILER='+tools['gcc'],'-DCMAKE_MAKE_PROGRAM='+tools['make'],
            '-DSOURCE_ROOT='+str(source),'-DBLAS_LIBRARIES='+libs['blas'],
            '-DLAPACK_LIBRARIES='+libs['blas'],'-DBLAS_RUNTIME_DIR='+str(Path(libs['blas']).parent),
            '-DBLA_VENDOR=OpenBLAS','-DBLAS_openblas_LIBRARY='+libs['blas'],
            '-DCMAKE_LIBRARY_PATH='+str(Path(libs['blas']).parent)]
    configure = command(args,out,'configure',env)
    compile_ = command([tools['cmake'],'--build',str(directory),'--target','replay','--parallel',
                        str(int(os.environ['SLURM_NTASKS']))],out,'compile',env)
    binary = pin(directory/'replay')
    dep = command(['/usr/bin/ldd',binary['path']],out,'libraries',env)
    linked = []
    for line in verify(dep['stdout']).read_text().splitlines():
        found = re.search(r'(?:=>\s+)?(/\S+)\s+\(',line)
        if found:
            linked.append(pin(found.group(1)))
    record = {'status':'built', 'binary':binary, 'linked_libraries':linked, 'implementation':m['implementation'],
              'source_archives':m['build']['sources'], 'configure':configure, 'compile':compile_}
    dump(out/'BUILD.json', record)
    return record


def evaluate(m, c, medium, binary, out, env):
    task = out/'tasks'/(c['case_id']+'__'+medium+'_charges')
    task.mkdir(parents=True)
    for name, record in (('core.xyz',c['xyz']),('charges',c['endpoints'][medium]['charges'])):
        shutil.copyfile(verify(record),task/name)
    local_env = dict(env, XTBPATH=str(task), XTBHOME=str(task))
    receipt = command([str(verify(binary)),str(task/'core.xyz'),str(task/'charges')],task,'native',local_env)
    text = verify(receipt['stdout']).read_text()
    values = re.findall(r'^COMPONENTS\s+(.+)$',text,re.M)
    bounds = re.findall(r'^ROUNDING_BOUNDS\s+(.+)$',text,re.M)
    if len(values) != 1 or len(bounds) != 1:
        raise ValueError('missing unique component output')
    energy = dict(zip(TERMS,map(float,values[0].split())))
    rounding = dict(zip(('Gelec','Ghb'),map(float,bounds[0].split())))
    if len(energy) != 4 or len(rounding) != 2 or not all(math.isfinite(v) for v in [*energy.values(),*rounding.values()]):
        raise ValueError('invalid native result')
    if 'internal GFN2-xTB/ALPB' not in text:
        raise ValueError('native parameter factory not proven')
    row = {'case_id':c['case_id'],'charge_source':medium,'components_Eh':energy,
           'rounding_bounds_Eh':rounding,'G_Eh':math.fsum(energy.values()),
           'execution':pin(task/'native_EXECUTION.json'), 'xyz':pin(task/'core.xyz'),
           'charges':pin(task/'charges'),'cavity':pin(task/'cavity.tsv'),
           'born_matrix_with_HB':pin(task/'born_matrix.dat')}
    dump(task/'RESULT.json',row)
    row['result'] = pin(task/'RESULT.json')
    return row


def qualify(m, rows, mp, build_record, out):
    checks = []
    for c in m['cases']:
        r = rows[c['case_id']]
        for term in TERMS:
            bound = r['rounding_bounds_Eh'].get(term,0.0)
            tolerance = bound+m['tolerances']['native_component_Eh']
            residual = r['components_Eh'][term]-c['archived_components_Eh'][term]
            informative = tolerance*m['HA_TO_KCAL'] <= m['tolerances']['maximum_informative_bound_kcal_mol']
            checks.append({'case_id':c['case_id'],'component':term,'residual_Eh':residual,
                           'rounding_bound_Eh':bound,'tolerance_Eh':tolerance,'informative':informative,
                           'passed':informative and abs(residual)<=tolerance})
    passed = all(r['passed'] for r in checks)
    q = {'status':'passed' if passed else 'failed','manifest':pin(mp),'build':build_record,
         'rows':rows,'checks':checks}
    dump(out/'QUALIFICATION.json',q)
    if not passed:
        raise ValueError('diagonal replay failed or rounding uncertainty uninformative; no vacuum-charge evaluations')
    return q


def main():
    p = argparse.ArgumentParser()
    p.add_argument('action',choices=('build-test-execute','build','test','execute'))
    p.add_argument('manifest')
    p.add_argument('--sha',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--build-receipt')
    p.add_argument('--qualification')
    a = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Slurm is mandatory; nothing created')
    mp = Path(a.manifest).resolve()
    if pin(mp)['sha256'] != a.sha:
        raise SystemExit('manifest pin differs')
    m = read(mp)
    out = Path(a.output).resolve()
    out.mkdir(parents=True,exist_ok=False)
    start = time.monotonic()
    result = {'status':'failed','terminal':True,'action':a.action,'manifest':pin(mp),
              'molecular_SCC_calls':0,'DFT_calls':0,'component_results':[]}
    env = {k:v for k,v in os.environ.items() if k not in ('LD_PRELOAD','PYTHONPATH','CPATH','LIBRARY_PATH','CMAKE_PREFIX_PATH')}
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',
               LD_LIBRARY_PATH=m['build']['library_directory'])
    try:
        dump(out/'PREFLIGHT.json',preflight(m,mp))
        if a.action in ('build','build-test-execute'):
            b = build(m,out,env)
        else:
            if not a.build_receipt:
                raise ValueError('existing exact build receipt required')
            b = read(a.build_receipt)
            if b['source_archives'] != m['build']['sources'] or b['implementation'] != m['implementation'] or b['status'] != 'built':
                raise ValueError('different source build')
        verify(b['binary'])
        for linked in b['linked_libraries']:
            verify(linked)
        result['build'] = b
        if a.action in ('test','build-test-execute'):
            water = {}
            for c in m['cases']:
                water[c['case_id']] = evaluate(m,c,'water',b['binary'],out,env)
                result['component_results'].append(water[c['case_id']])
            q = qualify(m,water,mp,b,out)
        elif a.action == 'execute':
            if not a.qualification:
                raise ValueError('actual diagonal qualification required')
            q = read(a.qualification)
            if q['status'] != 'passed' or q['manifest'] != pin(mp) or q['build']['binary'] != b['binary'] or not all(x['passed'] for x in q['checks']):
                raise ValueError('unqualified build/manifest')
            for r in q['rows'].values():
                verify(r['result'])
                verify(r['cavity'])
                verify(r['born_matrix_with_HB'])
        if a.action in ('execute','build-test-execute'):
            result['qualification'] = q
            result['cases'] = []
            for c in m['cases']:
                vac = evaluate(m,c,'vacuum',b['binary'],out,env)
                result['component_results'].append(vac)
                water = q['rows'][c['case_id']]
                for key in ('cavity','born_matrix_with_HB'):
                    if vac[key]['sha256'] != water[key]['sha256']:
                        raise ValueError('fixed cavity differs between charge distributions')
                # Preserve archived intrinsic remainder; no deletion or refit.
                delta_g = water['G_Eh']-vac['G_Eh']
                bracket = delta_g+c['intrinsic_remainder_Eh']
                closure = vac['G_Eh']+bracket-c['U_Eh']
                budget = sum(x['tolerance_Eh'] for x in q['checks'] if x['case_id']==c['case_id'])
                if abs(closure)>budget:
                    raise ValueError('U component accounting exceeds propagated diagonal budget')
                result['cases'].append({'case_id':c['case_id'],'source_case_id':c['source_case_id'],
                    'metal':c['metal'],'water_charge_replay':water,'vacuum_charge_replay':vac,
                    'U_archived_Eh':c['U_Eh'],'I_archived_Eh':c['intrinsic_remainder_Eh'],
                    'G_qvac_Eh':vac['G_Eh'],'G_qwater_minus_G_qvac_Eh':delta_g,
                    'redistribution_plus_I_Eh':bracket,'U_closure_residual_Eh':closure,
                    'U_closure_tolerance_Eh':budget})
            by_id = {r['case_id']:r for r in result['cases']}
            signs = {'sample1_Ca':1,'sample1_La':-1,'sample3_Ca':-1,'sample3_La':1}
            contrasts = {key:math.fsum(sign*by_id[name][key] for name,sign in signs.items())
                         for key in ('G_qvac_Eh','redistribution_plus_I_Eh','U_archived_Eh','U_closure_residual_Eh')}
            contrasts['sum_of_two_terms_Eh'] = contrasts['G_qvac_Eh']+contrasts['redistribution_plus_I_Eh']
            result['source_metal_double_contrast'] = {
                'definition':'sample1_Ca−sample1_La−sample3_Ca+sample3_La',
                'coefficients':signs,'hartree':contrasts,'HA_TO_KCAL':m['HA_TO_KCAL'],
                'kcal_mol':{k.removesuffix('_Eh'):v*m['HA_TO_KCAL'] for k,v in contrasts.items()}}
        result['status'] = 'complete'
    except Exception as exc:
        result['failure'] = repr(exc)
    finally:
        result['elapsed_seconds'] = time.monotonic()-start
        result['completed_component_evaluations'] = len(result['component_results'])
        result['slurm_job_id'] = os.environ['SLURM_JOB_ID']
        result['hostname'] = socket.gethostname()
        dump(out/'FINAL_COLLECTION.json',result)
    if result['status'] != 'complete':
        raise SystemExit(result['failure'])


if __name__ == '__main__':
    main()
