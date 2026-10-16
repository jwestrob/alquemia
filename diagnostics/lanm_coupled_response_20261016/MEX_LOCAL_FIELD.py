#!/usr/bin/env python3
"""Four fixed local Mex GFN2 gradient calls; separate Slurm-only prepare/run.

CLI: --stage prepare|run --config CONFIG --output FRESH_DIRECTORY
Config pins: implementation, plan, inputs, solvent_result (12510),
force_mapping_helper, mapping_dependencies{name:pin}, api_evidence{name:pin};
workers=4, repository.
run additionally requires preparation_result. Runtime/Python/parsers/settings
inherit unchanged from12510/12488. No model import/evaluation on login.
"""
from concurrent.futures import ThreadPoolExecutor
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import socket
import sys
import time
import traceback

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
PROTOCOL = 'mex209_fixed_local_gfn2_field_scout_v1'
BOHR_A = 0.529177210903
HA_KCAL = 627.509474
Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'La': 57, 'Nd': 60, 'Dy': 66}
TOL = {'charge_closure_e': 1e-5, 'xtb_json_stdout_hartree': 1e-7,
       'engrad_energy_hartree': 1e-8, 'engrad_coordinates_A': 1e-6,
       'xtb_json_population_e': 1e-7}
CONTROL = '$embedding\n interface=orca\n input=external.pc\n gradient=pcgrad\n$end\n'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return {'path': str(path), 'sha256': h.hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def module(name, path):
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def xyz(p):
    lines = verify(p).read_text().splitlines()
    rows = [s.split() for s in lines[2:] if s.strip()]
    require(int(lines[0]) == len(rows) == 209 and all(len(r) == 4 for r in rows), 'Local XYZ layout')
    symbols = [r[0] for r in rows]
    positions = [[float(v) for v in r[1:]] for r in rows]
    require(all(math.isfinite(v) for r in positions for v in r), 'Nonfinite XYZ')
    return symbols, positions


def inputs(config):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required for every stage')
    cfg = json.loads(config.read_text())
    require(Path(cfg['repository']).resolve() == ROOT and cfg['workers'] == 4, 'Scope/workers differ')
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation changed')
    verify(cfg['plan']); verify(cfg['inputs']); verify(cfg['force_mapping_helper'])
    require(verify(cfg['inputs']).resolve() == ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Mex8FNS/INPUTS.json',
            'Declared repaired Mex source differs')
    require(set(cfg['mapping_dependencies']) == {'metal_environment_reference.py', 'metal_environment_electronic_state.py'},
            'Mapping dependency pins absent')
    for name, p in cfg['mapping_dependencies'].items():
        require(verify(p).resolve() == ROOT/'scripts'/name, 'Mapping dependency source differs')
    expected_api = {'header', 'library', 'embedding_source', 'control_source', 'calculator_source', 'api_calculator_source', 'api_results_source'}
    require(set(cfg['api_evidence']) == expected_api, 'Incomplete native interface evidence')
    for p in cfg['api_evidence'].values():
        verify(p)
    solvent = read(cfg['solvent_result']); scfg = read(solvent['config'])
    require(solvent['status'] == 'complete' and solvent['slurm_job_id'] == '12510'
            and solvent['accepted'] == solvent['attempted'] == 8, 'Parent standalone result differs')
    require(scfg['tolerances'] == TOL, 'Parent tolerances differ')
    require(pin(sys.executable)['sha256'] == scfg['python']['sha256'], 'Python differs')
    qualification = read(scfg['qualification']); qcfg = read(qualification['config'])
    qmanifest = read(qualification['manifest'])
    require(qualification['status'] == 'passed' and qualification['slurm_job_id'] == '12488'
            and qualification['accepted'] == qualification['attempted'] == 8
            and all(p['passed'] for p in qualification['parity'] + qualification['solvent_component_parity']),
            'Qualified standalone lineage differs')
    runtime = scfg['runtime']
    for p in [runtime['executable'], runtime['parameters'], *runtime['runtime_libraries']]:
        verify(p)
    accepted = [r['task'] for r in qualification['rows'] if r['task']['backend'] == 'standalone' and r['status'] == 'accepted']
    require(len(accepted) == 4 and all({k:t[k] for k in runtime} == runtime for t in accepted), 'Qualified runtime differs')
    require(qcfg['staging']['standalone_settings'] == {'accuracy': 0.0001, 'maxiter': 500,
                'electronic_temperature_K': 300, 'alpb_reference': 'gsolv'}, 'Qualified settings differ')
    profile = read(qmanifest['profile'])
    sys.path.insert(0, str(ROOT / 'scripts'))
    import numpy as np
    backend = module('mex_field_backend', verify(scfg['backend_helper']))
    gradients = module('mex_field_gradients', verify(scfg['gradient_helper']))
    common = module('mex_field_common', verify(scfg['common_helper']))
    mapping = module('mex_field_mapping', verify(cfg['force_mapping_helper']))
    require(common.BOHR_TO_A == BOHR_A and common.HA_TO_KCAL == HA_KCAL, 'Unit convention changed')
    env = os.environ.copy()
    native = [r['task']['environment'] for r in qualification['rows'] if r['task']['backend'] == 'native']
    require(len(native) == 4 and all(e == native[0] for e in native), 'Inherited environment differs')
    for setting in (native[0], runtime['environment']):
        for k, v in setting.items():
            if v is None:
                env.pop(k, None)
            else:
                env[k] = v
    for k in ('OMP_PLACES', 'OMP_PROC_BIND'):
        env.pop(k, None)
    return cfg, scfg, qmanifest, profile, runtime, np, backend, gradients, mapping, env


def prepare(config, out, world):
    cfg, scfg, qmanifest, profile, runtime, np, backend, gradients, mapper, env = world
    source = read(cfg['inputs'])
    require(source['charge_accounting'] == {'MM_e': 2.0, 'QM_e': 0, 'all_site_ions_e': 12.0,
                'protein_water_e': -10, 'total_e': 2.0}, 'Physical source charge accounting differs')
    c = source['configurations']['A']; field = read(source['environment_atoms'])
    boundary = read(source['boundary_mapping'])
    rows = verify(c['pointcharges']).read_text().splitlines()
    require(int(rows[0]) == len(rows)-1 == len(field) == 1840, 'Actual field count differs')
    require(len({a['id'] for a in field}) == 1840, 'Duplicate field source IDs')
    require(all(a['kind'] in {'protein_source', 'retained_crystal_water', 'spectator_formal_monopole'} for a in field),
            'Unresolved virtual field site')
    require({a['element'] for a in field} == {'H', 'C', 'N', 'O', 'Nd'}, 'Field element inventory differs')
    require(abs(math.fsum(a['charge_e'] for a in field)-2) <= TOL['charge_closure_e'], 'Field charge closure')
    gam = {}
    for block in re.split(r'(?m)^\$Z=\s*', verify(runtime['parameters']).read_text())[1:]:
        number = int(block.split()[0])
        value = re.search(r'(?m)^\s*GAM=\s*([-+\d.Ee]+)\s*$', block)
        if value:
            gam[number] = float(value.group(1))
    expected_gam = {'H': .405771, 'C': .538015, 'N': .461493, 'O': .451896, 'Nd': .658531}
    require(all(gam[Z[s]] == v for s,v in expected_gam.items()), 'Author hardness differs')
    external = ['1840']; metadata = []
    for i, (line, a) in enumerate(zip(rows[1:], field)):
        tokens = line.split()
        require(len(tokens) == 4, 'Source point-charge format differs')
        values = [float(v) for v in tokens]
        require(abs(values[0]-a['charge_e']) <= 1e-12
                and np.allclose(values[1:], a['xyz_A'], atol=1e-8, rtol=0), 'Field identity/order/coordinates differ')
        hardness = gam[Z[a['element']]]
        external.append(' '.join(tokens) + ' ' + format(hardness, '.17g'))
        metadata.append({'pc_index': i, 'source_id': a['id'], 'element': a['element'], 'charge_e': values[0],
                         'xyz_A': values[1:], 'GAM_atomic_units': hardness, 'source_kind': a['kind']})
    external_text = '\n'.join(external) + '\n'
    (out/'external.pc').write_text(external_text)
    write(out/'FIELD.json', {'source': c['pointcharges'], 'source_mapping': source['environment_atoms'],
          'boundary_mapping': source['boundary_mapping'], 'parameters': runtime['parameters'], 'rows': metadata,
          'convention': 'ORCA-format Angstrom coordinates; explicit element-specific author GAM; no virtual sites',
          'pc_gradient_layout': 'headerless1840x3,input-order,Eh/bohr,8decimal', 'pc_gradient_rounding_Eh_bohr': 5e-9,
          'external_reference': 'nonperiodic hardness-damped shell potential tending to zero at infinity; no PC-PC self-energy'})
    maps = {}; paired = {}; tasks = []
    for metal in ('La', 'Dy'):
        state = c['endpoints'][metal]
        require(state['charge'] == 0 and state['oxidation_state'] == 3
                and state['multiplicity'] == (1 if metal == 'La' else 6), 'Local physical state differs')
        symbols, positions = xyz(state['xyz']); paired[metal] = (symbols, positions)
        require(set(symbols) <= {'H','C','N','O',metal} and symbols.count(metal) == 1, 'Local element inventory')
        electrons = sum(sum(profile['element'][s]['refocc']) for s in symbols)
        require(int(electrons) == electrons and int(electrons) % 2 == 0, 'Effective f-in-core electron parity')
        mapping = mapper.prepare_mapping(verify(cfg['inputs']), 'A', metal)
        require(len(mapping['core']) == 209 and len(mapping['field']) == 1840
                and len(mapping['physical_atoms']) == 2049 and len(mapping['caps']) == 4, 'Full source/cap map inventory')
        missing = set(a['id'] for a in mapping['physical_atoms']) - set(a['id'] for a in field)
        real = set(mapper.target_source_id(a) if a['kind']=='metal' else a['source_id']
                   for a in mapping['core'] if a['kind'] != 'cap')
        require(missing-real == {f'A/{r}//CA' for r in (66,67,83,96)}, 'Zero-field physical boundary inventory differs')
        write(out/f'MAPPING_{metal}.json', mapping); maps[metal] = pin(out/f'MAPPING_{metal}.json')
        for kind in ('isolated', 'embedded'):
            directory = out/f'{metal}_{kind}'; directory.mkdir()
            (directory/'core.xyz').write_bytes(verify(state['xyz']).read_bytes())
            command = [runtime['executable']['path'], 'core.xyz', '--gfn','2','--chrg','0','--uhf','0',
                       '--etemp','300','--acc','0.0001','--iterations','500','--norestart','--json','--pop','--grad']
            staged = {'core.xyz': pin(directory/'core.xyz')}
            if kind == 'embedded':
                (directory/'external.pc').write_text(external_text)
                (directory/'control.inp').write_text(CONTROL)
                command += ['--input','control.inp']
                staged.update({name:pin(directory/name) for name in ('external.pc','control.inp')})
            task = {'task_id': f'Mex8FNS_{metal}_{kind}', 'metal': metal, 'kind':kind,
                    'atom_count':209,'pointcharge_count':1840 if kind=='embedded' else 0,
                    'charge':0,'multiplicity':1,'uhf':0,'physical_multiplicity':state['multiplicity'],
                    'physical_electronic_state':state['electronic_state'],'oxidation_state':3,
                    'electronic_state':'original published spin-independent f-in-core GFN2 effective singlet',
                    'valence_electrons':int(electrons), 'source_xyz':state['xyz'], 'xyz':staged['core.xyz'],
                    'mapping':maps[metal],'source_pointcharges':c['pointcharges'], 'field':pin(out/'FIELD.json'),
                    'staged_files':staged,'command':command}
            write(directory/'TASK.json', task); tasks.append(dict(task, task_record=pin(directory/'TASK.json')))
    require(paired['La'][1] == paired['Dy'][1]
            and [s if s!='La' else 'Ln' for s in paired['La'][0]] == [s if s!='Dy' else 'Ln' for s in paired['Dy'][0]],
            'La/Dy do not share exact coordinates/identity')
    return {'status':'prepared','tasks':tasks,'mappings':maps,'field':pin(out/'FIELD.json'),
            'validation':'all4states,209QM,1840PC,2049physical,4caps,exactpairedXYZ,authorGAM,charge/mapping passed',
            'molecular_calls':0,'runtime':runtime,'solvent_result':cfg['solvent_result'],
            'profile':qmanifest['profile'],'source_inputs':cfg['inputs'],'api_evidence':cfg['api_evidence'],
            'source_boundary':source['boundary_mapping']}


def admit(directory, task, execution, world):
    cfg, scfg, qmanifest, profile, runtime, np, backend, gradients, mapper, env = world
    require(execution['returncode'] == 0 and execution['owned_process_group_empty'], 'Process failure')
    data = backend.parse_backend(directory, 'xtb', task)
    require(data['status'] == 'available', 'SCC/population failure: ' + str(data.get('failure_reason')))
    energy, charges = data['energy_hartree'], data['charges_e']
    payload = json.loads((directory/'xtbout.json').read_text())
    require(payload['method'] == 'GFN2-xTB' and payload['xtb version'].startswith('6.7.1'), 'Actual method/version differs')
    require(payload['number of electrons'] == task['valence_electrons']
            and payload['number of unpaired electrons'] == 0, 'Actual effective electronic state differs')
    require(abs(payload['total energy']-energy) <= TOL['xtb_json_stdout_hartree'], 'JSON/stdout energy mismatch')
    require(len(payload['partial charges']) == 209 and all(math.isfinite(q) for q in payload['partial charges'])
            and max(abs(a-b) for a,b in zip(payload['partial charges'],charges)) <= TOL['xtb_json_population_e'],
            'JSON/charge-file populations differ')
    require(abs(math.fsum(charges)) <= TOL['charge_closure_e'], 'QM charge closure')
    text = (directory/'stdout.txt').read_text()
    combined = text + '\n' + (directory/'stderr.txt').read_text()
    require('normal termination of xtb' in combined.lower(), 'Normal termination absent')
    require(not re.search(r'failed to converge|NOT CONVERGED|abnormal termination|numerical gradient', combined, re.I),
            'SCC failure or unsupported numerical gradient')
    for name, expected in ((r'net charge',0),(r'unpaired electrons',0),(r'# electrons',task['valence_electrons']),
                           (r'electronic temp\.',300),(r':  accuracy',.0001)):
        values = re.findall(name+r'\s+([-+\d.Ee]+)', text)
        require(bool(values) and float(values[-1]) == expected, 'Actual setting differs: '+name)
    embedded = task['kind'] == 'embedded'
    for pattern in (r'max\. iterations\s+500',r'restarted\?\s+false',r'GBSA solvation\s+false',
                    r'PC potential\s+'+('true' if embedded else 'false'),
                    r'convergence criteria satisfied after'):
        require(bool(re.search(pattern,text)), 'Actual SCC setting/evidence absent: '+pattern)
    require('Solvation model:               ALPB' not in text, 'Unexpected solvent')
    if embedded:
        values = re.findall(r'# point charges\s+(\d+)',text)
        totals = re.findall(r'sum of PC\s+([-+\d.Ee]+)',text)
        require(bool(values) and int(values[-1]) == 1840 and bool(totals)
                and abs(float(totals[-1])-2) <= TOL['charge_closure_e'], 'Printed field inventory/charge differs')
    else:
        require(not re.search(r'# point charges',text) and not (directory/'pcgrad').exists(), 'Unexpected isolated field')
    engrad = [directory/n for n in ('core.engrad','xtb-orca.engrad') if (directory/n).is_file()]
    require(len(engrad) == 1, 'Expected one native analytic QM gradient output')
    header = engrad[0].read_text()
    require('Eh/bohr' in header and 'coordinates in Bohr' in header, 'Gradient unit header absent')
    g = gradients.read_engrad(engrad[0]); symbols, positions = xyz(task['xyz'])
    require(g['atom_count'] == 209 and abs(g['energy_Ha']-energy) <= TOL['engrad_energy_hartree'],
            'Gradient count/energy mismatch')
    require(np.array_equal(g['atomic_numbers'],[Z[s] for s in symbols])
            and np.allclose(g['coordinates_bohr']*BOHR_A,positions,atol=TOL['engrad_coordinates_A'],rtol=0),
            'Gradient atom order/coordinates differ')
    if embedded:
        records = [line for line in (directory/'pcgrad').read_text().splitlines() if line.strip()]
        # Native writes3f12.8: adjacent full-width negative fields need not have
        # whitespace between them. Parse the actual fixed-width format.
        require(len(records) == 1840 and all(len(line) == 36 for line in records), 'Native pcgrad layout/count mismatch')
        rows = [[line[k:k+12] for k in (0,12,24)] for line in records]
        pc = np.array([[float(v) for v in r] for r in rows],dtype=float)
        require(pc.shape == (1840,3) and np.isfinite(pc).all(), 'Invalid external-site gradient')
    else:
        pc = np.zeros((1840,3))
    mapping = read(task['mapping'])
    physical = mapper.scatter(mapping,g['gradient_Ha_per_bohr'],pc)
    cycles = re.findall(r'convergence criteria satisfied after\s+(\d+) iterations',text)
    require(bool(cycles), 'SCC iteration receipt absent')
    return {'energy_hartree':energy,'charges_e':charges,'charge_closure_e':abs(math.fsum(charges)),
            'solver_observation':data,'scc_iterations':int(cycles[-1]),'xtbout':pin(directory/'xtbout.json'),
            'qm_gradient':{'artifact':pin(engrad[0]),'hartree_per_bohr':g['gradient_Ha_per_bohr'].tolist(),
                           'kcal_mol_A':(g['gradient_Ha_per_bohr']*HA_KCAL/BOHR_A).tolist()},
            'pointcharge_gradient':{'artifact':pin(directory/'pcgrad') if embedded else None,
                 'source':'native analytic pcgrad' if embedded else 'exact zero: isolated Hamiltonian has no external field',
                 'hartree_per_bohr':pc.tolist(),'kcal_mol_A':(pc*HA_KCAL/BOHR_A).tolist(),
                 'rounding_bound_per_component_Eh_bohr':5e-9 if embedded else 0.,'field_mapping':task['field']},
            'mapped_physical_gradient':physical,
            'quantity':'gradient=dE/dR; force=-gradient; native analytic return, no independent finite-difference qualification',
            'gradient_units':'hartree/bohr and kcal/mol/angstrom',
            'translation_gradient_sum_Eh_bohr':(g['gradient_Ha_per_bohr'].sum(axis=0)+pc.sum(axis=0)).tolist()}


def execute(config, out, world):
    cfg, scfg, qmanifest, profile, runtime, np, backend, gradients, mapper, env = world
    prepared = read(cfg['preparation_result'])
    require(prepared['stage'] == 'prepare' and prepared['status'] == 'prepared' and prepared['molecular_calls'] == 0,
            'Zero-call preparation unavailable')
    prep_cfg = read(prepared['config'])
    require({k:v for k,v in cfg.items() if k != 'preparation_result'} == prep_cfg,
            'Scientific config changed after preparation')
    require(prepared['implementation']['sha256'] == pin(__file__)['sha256'] and prepared['plan'] == cfg['plan']
            and prepared['runtime'] == runtime and prepared['source_inputs'] == cfg['inputs']
            and prepared['api_evidence'] == cfg['api_evidence'], 'Preparation lineage differs')
    tasks = prepared['tasks']
    require(len(tasks) == 4 and {(t['metal'],t['kind']) for t in tasks}
            == {(m,k) for m in ('La','Dy') for k in ('isolated','embedded')}, 'Four-state denominator differs')
    require(int(os.environ['SLURM_CPUS_PER_TASK']) >= cfg['workers'], 'Four one-thread workers require four allocated CPUs')
    # Validate every prepared task before any backend starts.
    for task in tasks:
        expected = dict(task); expected.pop('task_record')
        require(read(task['task_record']) == expected, 'Prepared TASK content differs')
        verify(task['mapping']); verify(task['field']); verify(task['source_xyz']); verify(task['source_pointcharges'])
        for p in task['staged_files'].values():
            verify(p)
        require(task['xyz']['sha256'] == task['source_xyz']['sha256'] and task['charge'] == 0
                and task['multiplicity'] == 1 and task['uhf'] == 0, 'Prepared local state differs')

    def cell(task):
        directory = out/task['task_id']; directory.mkdir()
        actual = dict(task)
        staged = {}
        for name,p in task['staged_files'].items():
            (directory/name).write_bytes(verify(p).read_bytes())
            staged[name] = pin(directory/name)
            require(staged[name]['sha256'] == p['sha256'], 'Execution copy differs')
        actual.update(xyz=staged['core.xyz'],staged_files=staged,preparation_task=task['task_record'])
        write(directory/'TASK.json',actual)
        row = {'task':actual,'status':'unavailable','attempted':False,'energy_hartree':None}
        try:
            for p in [runtime['executable'],runtime['parameters'],*runtime['runtime_libraries']]:
                verify(p)
            row['attempted'] = True
            execution = backend.call(actual['command'],directory,env)
            row['execution'] = pin(directory/'EXECUTION.json')
            for p in staged.values():
                verify(p)
            row.update(admit(directory,actual,execution,world),status='accepted')
        except Exception as exc:
            row.update(failure=repr(exc),traceback=traceback.format_exc())
        row['artifacts'] = {p.name:pin(p) for p in directory.iterdir() if p.is_file()}
        write(directory/'RESULT.json',row)
        row['result'] = pin(directory/'RESULT.json')
        print(json.dumps({'task':actual['task_id'],'status':row['status'],'failure':row.get('failure')}),flush=True)
        return row

    with ThreadPoolExecutor(max_workers=cfg['workers']) as pool:
        rows = list(pool.map(cell,tasks))
    by = {(r['task']['metal'],r['task']['kind']):r for r in rows}
    differences = []
    for metal in ('La','Dy'):
        isolated,embedded = (by[(metal,k)] for k in ('isolated','embedded'))
        item = {'metal':metal,'status':'unavailable','definition':'embedded minus isolated local GFN2 only',
                'isolated':isolated['result'],'embedded':embedded['result']}
        if isolated['status'] == embedded['status'] == 'accepted':
            de = math.fsum([embedded['energy_hartree'],-isolated['energy_hartree']])
            dg = np.array(embedded['mapped_physical_gradient']['gradient_kcal_mol_A']) - np.array(isolated['mapped_physical_gradient']['gradient_kcal_mol_A'])
            item.update(status='available',energy_hartree=de,energy_kcal_mol=de*HA_KCAL,
                        gradient_kcal_mol_A=dg.tolist(),force_kcal_mol_A=(-dg).tolist(),
                        physical_mapping=prepared['mappings'][metal],
                        qm_population_change_e=[a-b for a,b in zip(embedded['charges_e'],isolated['charges_e'])])
        differences.append(item)
    observed = [read(r['execution']) for r in rows if r.get('execution')]
    return {'status':'complete' if all(r['status']=='accepted' for r in rows) else 'unavailable',
            'preparation_result':cfg['preparation_result'],'rows':rows,'field_differences':differences,
            'attempted':sum(r['attempted'] for r in rows),'accepted':sum(r['status']=='accepted' for r in rows),
            'molecular_calls':sum(r['attempted'] for r in rows),'runtime':runtime,
            'observed_child_CPU_seconds':math.fsum(r['cpu_user_seconds']+r['cpu_system_seconds'] for r in observed),
            'peak_child_RSS_KiB':max((r['peak_RSS_KiB'] for r in observed),default=0),
            'MACE_calls':0,'wholeprotein_SCC_calls':0,'solvent_calls':0,'optimizer_steps':0,
            'scope':'local209 field correction only; no total hybrid/solvent/affinity/classifier; fixed actual source'}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage',choices=('prepare','run'),required=True)
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    config = args.config.resolve(); out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    start = time.monotonic(); cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
    result = {'protocol_id':PROTOCOL,'stage':args.stage,'status':'running','config':pin(config),
              'implementation':pin(__file__),'slurm_job_id':os.environ['SLURM_JOB_ID'],
              'host':socket.gethostname(),'allocated_cpus':cpus,'declared_cells':4}
    write(out/'START.json',result)
    try:
        world = inputs(config); result['plan'] = world[0]['plan']
        result.update(prepare(config,out,world) if args.stage == 'prepare' else execute(config,out,world))
    except Exception as exc:
        result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic()-start
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds']*cpus
    write(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'result':str(out/'RESULT.json')}),flush=True)
    return result['status'] in ('prepared','complete')


if __name__ == '__main__':
    if not main():
        raise SystemExit(1)
