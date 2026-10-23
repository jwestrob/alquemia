#!/usr/bin/env python3
"""Nine consumed fixtures; metadata preparation or one root-bound numerical run."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shutil
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
W = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
D = HERE.parent
THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
PROFILES = {'absent': None, 'zero': (0., 0.), 'direct': (.05, 0.), 'mutual': (0., -.05), 'combined': (.05, -.05)}
IONS = ('ion_La', 'ion_Dy')
PURE = 'aquo9_central_La__near__two_water'
DIMERS = tuple('aquo9_central_' + m + '__dimer_01' for m in ('La', 'Dy'))
TRIMERS = tuple('aquo9_central_' + m + '__' + k + '__metal_two_water' for m in ('La', 'Dy') for k in ('near', 'far'))
NEAR = tuple(g for g in TRIMERS if '__near__' in g)
CELLS = IONS + (PURE,) + DIMERS + TRIMERS
MAX_CALLS, MAX_FD = 56, 384
SOURCES = {
    'model78': (W/'saved_redistribution_counterpart_v1/recovery_preparation_v1/MANIFEST.json', '141048e19aed8337ead156d2aa9025fc546f8c2252cb64c23cb44f30a2239e92'),
    'trimer': (W/'collective_mechanics_options_v1/preparation_v1/MANIFEST.json', 'a8fdf0c1ed724657bd642862d6af52632854e48ad8006f68cb0d87c2aef6f34c'),
    'native': (D/'RUN_PANEL.py', '4189602145853bc39d8dfaf30915b2a6b1d4632984e0c50fc0bf321644fa306c'),
    'gradient': (D/'paired_water_transfer_v1/RUN.py', '27591f28fbbfcfc7f7ac8f607711fbdacc3b783b643bf0a849e9149a7b5150e0'),
    'induction': (D/'collective_response_design_v1/RUN.py', '9eb552e6ec65da01b3e62ac803740405ee7cc05ab45d1371294fdea5b3ca422d'),
    'original_testgrad': (W/'tinker_build/testgrad', 'd82b2872136785cadf0c0763aafef60ba1bdd89283f685850485a07a875f6dc2'),
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path']) == p, 'Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def save(path, obj):
    p = Path(path)
    tmp = p.with_name(p.name + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    os.replace(tmp, p)


def module(p, name):
    s = importlib.util.spec_from_file_location(name, verify(p))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def source_pins():
    result = {}
    for name, (path, digest) in SOURCES.items():
        result[name] = pin(path)
        require(result[name]['sha256'] == digest, 'Frozen source differs: ' + name)
    return result


def xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    rows = [line.split() for line in lines[1:] if line.strip()]
    require(len(rows) == n and all(int(r[0]) == i+1 for i, r in enumerate(rows)), 'Typed XYZ inventory')
    coords = [[float(x) for x in r[2:5]] for r in rows]
    require(all(math.isfinite(x) for r in coords for x in r), 'Finite XYZ')
    return lines[0], rows, coords


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def unit(v):
    length = norm(v)
    require(length > 1e-12, 'Undefined geometry direction')
    return [x/length for x in v]


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def path_geometry(coords, kind, q):
    """Unit-norm tangent at q=0; exact rigid rotation/translation path."""
    out = copy.deepcopy(coords)
    tangent = [[0., 0., 0.] for _ in coords]
    if kind == 'radial':
        u = unit([coords[1][j]-coords[0][j] for j in range(3)])
        tangent[0] = [-v/2 for v in u]
        for a in (1, 2, 3):
            tangent[a] = [v/2 for v in u]
        out = [[v+q*t for v, t in zip(r, tr)] for r, tr in zip(coords, tangent)]
    else:
        require(kind == 'rotation', 'Unknown directional path')
        axis = unit([1., 2., 3.])
        center = coords[1]
        relative = [[coords[a][j]-center[j] for j in range(3)] for a in (1, 2, 3)]
        velocity = [cross(axis, r) for r in relative]
        length = norm([v for r in velocity for v in r])
        require(length > 1e-12, 'Undefined rotational normalization')
        theta = q/length
        for a, r, v in zip((1, 2, 3), relative, velocity):
            tangent[a] = [x/length for x in v]
            axial = sum(x*y for x, y in zip(axis, r))
            out[a] = [center[j]+r[j]*math.cos(theta)+v[j]*math.sin(theta)+axis[j]*axial*(1-math.cos(theta)) for j in range(3)]
    require(abs(norm([v for r in tangent for v in r])-1.) < 1e-12, 'Path tangent is not normalized')
    return out, tangent


def inventory():
    jobs = []
    for profile in ('original', 'absent', 'zero'):
        jobs += [dict(id=profile+'__'+g, cell=g, profile=profile, program='testgrad', fd=False) for g in CELLS]
    jobs += [dict(id='combined__'+g, cell=g, profile='combined', program='testgrad', fd=g in DIMERS+TRIMERS) for g in CELLS]
    for profile in ('direct', 'mutual'):
        jobs += [dict(id=profile+'__'+g, cell=g, profile=profile, program='testgrad', fd=True) for g in NEAR]
    for g in NEAR:
        for kind in ('radial', 'rotation'):
            for h in (1e-4, 5e-5):
                for sign in (-1, 1):
                    jobs.append(dict(id=f'directional__{g}__{kind}__{h:g}__{sign:+d}', cell=g, profile='combined', program='analyze', fd=False, kind=kind, h=h, q=sign*h))
    require(len(jobs) == MAX_CALLS and len({r['id'] for r in jobs}) == MAX_CALLS, 'Program denominator')
    return jobs


def prepare(out):
    sources = source_pins()
    a, b = read(sources['model78']), read(sources['trimer'])
    rows = {}
    for g in CELLS:
        r = (a['rows'] if g in IONS+DIMERS else b['model_rows'])[g]
        directory = out/'fixtures'/g
        directory.mkdir(parents=True)
        for k, filename in (('xyz', 'system.xyz'), ('key', 'system.key'), ('parameters', 'amoeba09.prm')):
            shutil.copyfile(verify(r['input'][k]), directory/filename)
        files = {k: pin(directory/f) for k, f in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
        header, atoms, coords = xyz(files['xyz']['path'])
        require([int(t[5]) for t in atoms] == r['types'], 'Consumed types differ')
        require(len(atoms) == r['atoms'] == (1 if g in IONS else 6 if g == PURE else 4 if g in DIMERS else 7), 'Fixture atom count')
        require(all(files[k]['sha256'] == r['input'][k]['sha256'] for k in files), 'Fixture copy differs')
        key = Path(files['key']['path']).read_text()
        require('polar-pair-response' not in key.lower() and 'polar-print' in key.lower(), 'Unexpected response or missing residual readback')
        rows[g] = {'input': files, 'source_input': r['input'], 'atoms': r['atoms'], 'types': r['types'], 'frame_signatures': r['frame_signatures'], 'positions_A': coords}
    jobs = inventory()
    require(sum(6*rows[j['cell']]['atoms'] for j in jobs if j['fd']) == MAX_FD, 'Internal Cartesian count')
    p = {'status':'prepared_no_native_calls', 'implementation':pin(__file__), 'plan':pin(HERE/'QUALIFICATION_PLAN.md'), 'sources':sources, 'rows':rows, 'programs':jobs, 'maximum_native_programs':MAX_CALLS, 'maximum_internal_FD_energies':MAX_FD, 'native_calls':0}
    save(out/'MANIFEST.json', p)
    template = {'implementation':pin(__file__), 'manifest':pin(out/'MANIFEST.json'), 'prototype_testgrad':None, 'prototype_analyze':None, 'prototype_source_manifest':None, 'source_audit':None, 'pair_scope_fallback_verified':False, 'dual_field_paths_verified':False, 'unsupported_routes_rejected':False, 'root_authorized':False, 'maximum_native_programs':MAX_CALLS, 'maximum_internal_FD_energies':MAX_FD, 'numerical_threads':1, 'fit':False, 'physical_qualification':False}
    save(out/'RUN_CONFIG_TEMPLATE.json', template)
    return p


def identity(a, b):
    require(len(a['forces']) == len(b['forces']), 'Force replay shape')
    de = abs(a['energy']-b['energy'])
    df = max(abs(x-y) for ax, bx in zip(a['forces'], b['forces']) for x, y in zip(ax, bx))
    require(math.isfinite(de) and math.isfinite(df), 'Nonfinite identity')
    return {'energy_error_kcal_mol':de, 'force_max_error_kcal_mol_A':df, 'printed_exact':de == 0 and df == 0, 'passed':de <= 1e-4 and df <= 1e-3}


def run(config_path, out):
    start, cpu = time.monotonic(), time.process_time()
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    config_pin = pin(config_path)
    c = read(config_pin)
    b = read(c['manifest'])
    require(c['implementation'] == pin(__file__) == b['implementation'] and b['plan'] == pin(HERE/'QUALIFICATION_PLAN.md'), 'Source/plan binding')
    require(c['root_authorized'] is True and all(c[x] is True for x in ('pair_scope_fallback_verified','dual_field_paths_verified','unsupported_routes_rejected')), 'Root source audit/execution binding required')
    require(c['fit'] is False and c['physical_qualification'] is False and c['maximum_native_programs'] == MAX_CALLS and c['maximum_internal_FD_energies'] == MAX_FD and c['numerical_threads'] == 1, 'Finite numerical-only scope')
    for k in ('prototype_testgrad','prototype_analyze','prototype_source_manifest','source_audit'):
        verify(c[k])
    require(c['prototype_testgrad']['sha256'] != b['sources']['original_testgrad']['sha256'], 'Prototype not isolated')
    require(b['status'] == 'prepared_no_native_calls' and b['native_calls'] == 0 and b['sources'] == source_pins() and set(b['rows']) == set(CELLS) and b['programs'] == inventory(), 'Prepared inventory differs')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1')) == 1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0')) == 1 and all(os.environ.get(k) == '1' for k in THREADS), 'One CPU Slurm task and one numerical thread required')
    for row in b['rows'].values():
        for k, p in row['input'].items():
            verify(p)
            require(p['sha256'] == row['source_input'][k]['sha256'], 'Original fixture binding')
        for p in row['source_input'].values():
            verify(p)
    native = module(b['sources']['native'], 'pair_qualification_native')
    helper = module(b['sources']['gradient'], 'pair_qualification_gradient')
    induct = module(b['sources']['induction'], 'pair_qualification_induction')
    result = {'status':'running', 'implementation':pin(__file__), 'config':config_pin, 'manifest':c['manifest'], 'job_id':os.environ['SLURM_JOB_ID'], 'rows':{j['id']:{'status':'not_attempted'} for j in b['programs']}, 'fixtures':list(CELLS), 'native_program_invocations':0, 'internal_FD_energies_requested':0, 'internal_FD_energies_completed':0, 'directional_checks':{}, 'error':None, 'numerical_checks_passed':False, 'physical_accuracy_qualified':False}
    def persist():
        save(out/'RESULT.json', result)
    try:
        persist()
        for job in b['programs']:
            g, label = job['cell'], job['id']
            base, target = b['rows'][g], result['rows'][label]
            cell = out/'cells'/label
            cell.mkdir(parents=True)
            for k, name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
                shutil.copyfile(verify(base['input'][k]), cell/name)
            profile = job['profile']
            if profile not in ('original','absent'):
                dd, dm = PROFILES[profile]
                with (cell/'system.key').open('a') as f:
                    for metal in (500, 501):
                        f.write(f'\nPOLAR-PAIR-RESPONSE {metal} 36 {dd:.17g} {dm:.17g}\n')
            if 'q' in job:
                header, atoms, coords = xyz(cell/'system.xyz')
                displaced, tangent = path_geometry(coords, job['kind'], job['q'])
                lines = [header]
                for atom, position in zip(atoms, displaced):
                    lines.append(' '.join(atom[:2]+[f'{v:.17g}' for v in position]+atom[5:]))
                (cell/'system.xyz').write_text('\n'.join(lines)+'\n')
                target['direction'] = {'kind':job['kind'], 'q_A':job['q'], 'tangent':tangent}
            target.update(status='started', job=job, input={k:pin(cell/f) for k,f in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}, receipts=[])
            fd = 6*base['atoms'] if job['fd'] else 0
            require(result['native_program_invocations'] < MAX_CALLS and result['internal_FD_energies_requested']+fd <= MAX_FD, 'Native finite ceiling')
            result['native_program_invocations'] += 1
            result['internal_FD_energies_requested'] += fd
            persist()
            if job['program'] == 'testgrad':
                exe = b['sources']['original_testgrad'] if profile == 'original' else c['prototype_testgrad']
                argv = [str(verify(exe)), 'system.xyz', 'Y', 'Y', '0.00001'] if job['fd'] else [str(verify(exe)), 'system.xyz', 'Y', 'N']
            else:
                argv = [str(verify(c['prototype_analyze'])), 'system.xyz', 'E']
            output = cell/(job['program']+'.out')
            native.invoke(argv, cell, output, target['receipts'])
            target['induction'] = induct.induction(output)
            if job['program'] == 'testgrad':
                value = helper.gradient(output, base['atoms'], job['fd'], native)
                require(math.isfinite(value['energy']), 'Finite energy')
                target['value'] = value
                result['internal_FD_energies_completed'] += fd
                require(not job['fd'] or value['fd_passed'], 'Cartesian numerical derivative failed')
                if profile in ('absent','zero') or (profile == 'combined' and g in IONS+(PURE,)):
                    target['identity'] = identity(value, result['rows']['original__'+g]['value'])
                    require(target['identity']['passed'], 'Disabled/zero/unrelated-system E/F identity failed')
            else:
                helper.program_admission(output, 'analyze', native)
                import re
                energy = float(re.findall(r'Total Potential Energy\s*:\s*('+native.NUMBER+')', output.read_text())[0])
                require(math.isfinite(energy), 'Finite directional energy')
                target['value'] = {'energy':energy, 'output':pin(output)}
            target['status'] = 'complete'
            persist()
            if job['program'] == 'analyze' and job['q'] > 0:
                kind, h = job['kind'], job['h']
                _, tangent = path_geometry(b['rows'][g]['positions_A'], kind, 0.)
                center = result['rows']['combined__'+g]['value']
                analytic = sum(v*t for ar,tr in zip(center['gradient'],tangent) for v,t in zip(ar,tr))
                projection_bound = 5e-9*sum(abs(t) for tr in tangent for t in tr)
                minus = result['rows'][f'directional__{g}__{kind}__{h:g}__-1']['value']['energy']
                plus = result['rows'][f'directional__{g}__{kind}__{h:g}__+1']['value']['energy']
                slope = (plus-minus)/(2*h)
                rounding = 1e-8/(2*h)+projection_bound
                check = {'half_step_A':h,'finite_difference_slope':slope,'analytic_slope':analytic,'absolute_error':abs(slope-analytic),'printing_bound':rounding,'passed':abs(slope-analytic) <= .001+rounding}
                row = result['directional_checks'].setdefault(g+'__'+kind, {'steps':[], 'passed':None})
                row['steps'].append(check)
                persist()
                require(check['passed'], 'Directional derivative check failed')
                if len(row['steps']) == 2:
                    checks = row['steps']
                    difference = abs(checks[0]['finite_difference_slope']-checks[1]['finite_difference_slope'])
                    tolerance = .001+sum(x['printing_bound'] for x in checks)
                    row.update(two_step_difference=difference, two_step_tolerance=tolerance, passed=difference <= tolerance)
                    persist()
                    require(row['passed'], 'Directional two-step check failed')
        require(result['native_program_invocations'] == MAX_CALLS and result['internal_FD_energies_completed'] == MAX_FD and all(v['status'] == 'complete' for v in result['rows'].values()) and len(result['directional_checks']) == 4 and all(x['passed'] is True for x in result['directional_checks'].values()), 'Full numerical denominator required')
        result.update(status='metal_pair_numerical_qualification_complete', numerical_checks_passed=True)
    except Exception as exc:
        result.update(status='numerical_qualification_failed', error=repr(exc))
        for target in result['rows'].values():
            if target['status'] == 'started':
                target.update(status='failed', error=repr(exc))
    finally:
        now = resource.getrusage(resource.RUSAGE_CHILDREN)
        result['cost'] = {'wall_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':now.ru_utime-child.ru_utime,'child_system_CPU_seconds':now.ru_stime-child.ru_stime,'numerical_threads':1}
        persist()
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('prepare','run'))
    parser.add_argument('--out', required=True)
    parser.add_argument('--config')
    args = parser.parse_args()
    require(args.action != 'run' or args.config, 'Run config required')
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    r = prepare(out) if args.action == 'prepare' else run(args.config, out)
    print(json.dumps({'status':r['status'],'out':str(out)}))
    if r['status'] == 'numerical_qualification_failed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
