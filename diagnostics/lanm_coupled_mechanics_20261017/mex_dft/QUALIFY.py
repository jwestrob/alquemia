#!/usr/bin/env python3
"""Eight fixed Mex DFT derivative endpoints; existing calculators, Slurm only."""
import argparse
import copy
import datetime as dt
from decimal import Decimal
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import sys
import time
import traceback

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
PROTOCOL = 'mex_embedded_PBE0_total_derivative_eight_SP_v1'
SETTINGS = {'h_A': .001, 'absolute_kcal_mol_A': .005, 'relative': .0001,
            'maximum_displaced_SP': 8, 'new_central_gradients': 0,
            'directions': ['MM_O97', 'boundary_CA83']}
METALS = ('La', 'Dy')
COMPONENTS = ('PBE0_embedded', 'MM_ledger', 'fixed_FF_charge_ALPB')
IDS = tuple(f'{d}_{m}_{s}' for d in SETTINGS['directions'] for m in METALS for s in ('plus', 'minus'))
SCOUT = IDS[:2]
HA = 627.509474


def require(ok, message):
    if not ok: raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: '+p['path'])
    return Path(p['path'])


def read(p): return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n')
    return pin(path)


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj); return obj


def canonical(value): return json.loads(json.dumps(value, allow_nan=False))


def inputs(cfg):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    require(cfg['protocol'] == PROTOCOL and cfg['settings'] == SETTINGS, 'Protocol differs')
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Adapter differs')
    require(pin(sys.executable)['sha256'] == cfg['python']['sha256'], 'Interpreter differs')
    for p in (cfg['plan'], cfg['specification'], cfg['origin_result'], *cfg['helpers'].values()): verify(p)
    # Load the proven allocation-aware runner before other helpers can import it.
    sys.path.insert(0, str(ROOT/'scripts'))
    sys.path.insert(0, str(verify(cfg['native_implementation']['run_orca_task_manifest.py']).parent))
    native = {}
    for name in ('render_orca_runtime_input', 'run_orca_task_manifest', 'ggr_sensitivity', 'mace_omol_vacuum', 'affordable_common'):
        obj = importlib.import_module(name)
        require(pin(obj.__file__)['sha256'] == cfg['native_implementation'][name+'.py']['sha256'], 'Imported native helper differs: '+name)
        native[name] = obj
    for p in cfg['native_implementation'].values(): verify(p)
    original = read(cfg['origin_result'])
    require(original['status'] == 'complete' and original['slurm_job_id'] == '12849'
            and original['molecular_calls'] == 0 and original['accepted'] == 2, 'Compatible saved origins required')
    helper = module('mex_dft_origin_replay', original['implementation'])
    acfg, prior, ref, world = helper.inputs(verify(original['config']))
    np = world[-1]; origins = {}
    for old in original['rows']:
        metal = old['metal']; replay = helper.assemble(acfg, ref, world, metal)
        require(old['status'] == 'complete' and old['result'] and read(old['result']) == {k:v for k,v in old.items() if k != 'result'}, 'Origin row receipt differs')
        require(old['native_energies'] == replay['native_energies'] and old['physical_atoms'] == replay['physical_atoms']
                and old['electronic_state'] == replay['electronic_state'], 'Origin energy/state/source differs')
        for k in COMPONENTS:
            ref.close(np, replay['component_gradients_kcal_mol_A'][k], old['component_gradients_kcal_mol_A'][k], 'Origin '+k)
        ref.close(np, replay['gradient_kcal_mol_A'], old['gradient_kcal_mol_A'], 'Origin total')
        origins[metal] = old  # Admitted archived derivatives remain authoritative.
    require(set(origins) == set(METALS) and origins['La']['source_ids'] == origins['Dy']['source_ids']
            and np.array_equal([a['xyz_A'] for a in origins['La']['physical_atoms']], [a['xyz_A'] for a in origins['Dy']['physical_atoms']]), 'Origin pair differs')
    manifest = world[7]
    require(manifest['orca'] == cfg['orca'] and manifest['assets'] == cfg['basis_assets'], 'Native executable/basis differs')
    require(cfg['ALPB_binary'] == original['solvent_build']['binary']
            and cfg['helpers']['ALPB'] == world[1]['implementation'], 'Shared solvent implementation differs')
    verify(cfg['orca'])
    for assets in cfg['basis_assets'].values():
        for p in assets.values(): verify(p)
    return original, ref, world, origins, native


def directions(world, origins):
    np = world[-1]; ids = origins['La']['source_ids']; ix = {s:i for i,s in enumerate(ids)}
    x = np.asarray([a['xyz_A'] for a in origins['La']['physical_atoms']])
    require([ix[s] for s in ('A/97//O','A/83//CA','A/83//C','A/83//N')] == [1019,801,803,799], 'Fixture IDs differ')
    bonds = {frozenset(b) for b in world[9]}
    for b in [('A/97//C','A/97//O'),('A/83//C','A/83//CA'),('A/83//CA','A/83//N')]:
        require(frozenset(b) in bonds, 'Fixture source bond absent')
    mp = read(origins['La']['physical_mapping']); cap = next(a for a in mp['core'] if a['qm_index'] == 205)
    require(mp['field'][832]['id'] == 'A/97//O' and mp['field'][832]['charge_e'] == -.5679
            and all(a['id'] != 'A/83//CA' for a in mp['field'])
            and cap['retained_source_id'] == 'A/83//C' and cap['omitted_source_id'] == 'A/83//CA'
            and cap['length_A'] == 1.09, 'Fixture cap/field relation differs')
    def unit(v):
        n = float(np.linalg.norm(v)); require(n > 1e-8, 'Degenerate direction'); return v/n
    u = unit(x[801]-x[803]); w = x[799]-x[801]
    return x, {'MM_O97': {'index':1019, 'id':'A/97//O', 'direction':unit(x[1019]-x[ix['A/97//C']]).tolist()},
               'boundary_CA83': {'index':801, 'id':'A/83//CA', 'direction':unit(w-u*np.dot(u,w)).tolist()}}


def update_geometry(world, origin, y, out, fixture):
    """Same physical/core/field/cap update as MEX_RESPONSE, without unused GFN tasks."""
    np = world[-1]; mapper = world[12]; source = world[5]; export = read(source['export'])
    mp = read(origin['physical_mapping']); ids = origin['source_ids']; ix = {s:i for i,s in enumerate(ids)}
    x = np.asarray([a['xyz_A'] for a in origin['physical_atoms']]); changed = set(np.flatnonzero(np.any(y != x, axis=1)))
    require(changed == {fixture['index']}, 'Single physical displacement changed unexpected atoms')
    def moved(atoms):
        a = copy.deepcopy(atoms)
        for row in a: row['xyz_A'] = y[ix[row['id']]].tolist()
        return a
    physical = write(out/'physical.json', moved(origin['physical_atoms']))
    atom_pin = write(out/'atoms.json', moved(read(export['atoms'])))
    ep = write(out/'EXPORT.json', dict(export, atoms=atom_pin, origin_export=source['export']))
    field = copy.deepcopy(mp['field']); lines = verify(mp['pointcharges']).read_text().splitlines(keepends=True)
    field_changes = []
    for i,a in enumerate(field):
        if ix[a['id']] not in changed: continue
        field_changes.append(i); a['xyz_A'] = y[ix[a['id']]].tolist()
        if 'source_xyz_A' in a: a['source_xyz_A'] = a['xyz_A']
        if 'source' in a: a['source']['xyz_A'] = a['xyz_A']
        tok = lines[i+1].split(); require(len(tok) == 4, 'Source field format differs')
        lines[i+1] = tok[0]+' '+' '.join(format(v,'.17g') for v in a['xyz_A'])+'\n'
    env = write(out/'environment_atoms.json', field); (out/'environment.pc').write_text(''.join(lines))
    boundary = copy.deepcopy(read(source['boundary_mapping']))
    boundary.update(original_boundary_mapping=source['boundary_mapping'], dipole_diagnostics_scope='historical source only; fixed shift membership')
    for sid in boundary['retained_physical_CA_coordinates']: boundary['retained_physical_CA_coordinates'][sid] = y[ix[sid]].tolist()
    bp = write(out/'boundary_mapping.json', boundary)
    core = copy.deepcopy(mp['core']); xyz_lines = verify(mp['xyz']).read_text().splitlines(keepends=True); core_changes = []
    for a in core:
        qi = a['qm_index']
        if a['kind'] == 'cap':
            ri,oi = ix[a['retained_source_id']],ix[a['omitted_source_id']]
            if not ({ri,oi} & changed): continue
            r,o = y[ri],y[oi]; delta = o-r; dist = float(np.linalg.norm(delta)); require(dist > 0, 'Coincident cap anchors')
            u = delta/dist; jo = a['length_A']/dist*(np.eye(3)-np.outer(u,u))
            a.update(xyz_A=(r+a['length_A']*u).tolist(), retained_xyz_A=r.tolist(), omitted_xyz_A=o.tolist(),
                     jacobian_retained=(np.eye(3)-jo).tolist(), jacobian_omitted=jo.tolist())
        else:
            sid = mapper.target_source_id(a) if a['kind'] == 'metal' else a['source_id']
            if ix[sid] not in changed: continue
            a['xyz_A'] = y[ix[sid]].tolist()
            if a['kind'] == 'metal': a['source_site']['xyz_A'] = a['xyz_A']
        core_changes.append(qi); tok = xyz_lines[qi+2].split()
        xyz_lines[qi+2] = tok[0]+' '+' '.join(format(v,'.17g') for v in a['xyz_A'])+'\n'
    (out/'core.xyz').write_text(''.join(xyz_lines)); cp = write(out/'core_mapping.json', core)
    require((field_changes,core_changes) == ([832],[]) if fixture['id'] == 'A/97//O' else (field_changes,core_changes) == ([],[205]), 'Fixture electronic updates differ')
    endpoint = copy.deepcopy(source['configurations']['A']['endpoints'][origin['metal']]); endpoint['xyz'] = pin(out/'core.xyz')
    inp = copy.deepcopy(source); inp.update(origin_inputs=origin['source_inputs'], export=ep, environment_atoms=env, boundary_mapping=bp,
        physical_source=atom_pin, physical_sources={'A':atom_pin}, core_mapping={'A':cp},
        configurations={'A':{'endpoints':{origin['metal']:endpoint},'core_mapping':cp,'pointcharges':pin(out/'environment.pc')}})
    ip = write(out/'INPUTS.json',inp); mapping = mapper.prepare_mapping(verify(ip),'A',origin['metal'])
    require([a['id'] for a in mapping['physical_atoms']] == ids and np.array_equal([a['xyz_A'] for a in mapping['physical_atoms']],y), 'Current physical mapping differs')
    require([a['charge_e'] for a in mapping['field']] == [a['charge_e'] for a in mp['field']], 'External charge vector changed')
    return {'physical':physical,'source_inputs':ip,'mapping':write(out/'MAPPING.json',mapping),
            'xyz':pin(out/'core.xyz'),'pointcharges':pin(out/'environment.pc'),'core_mapping':cp,'boundary_mapping':bp}


def prepare(cfg, out, result):
    original,ref,world,origins,native = inputs(cfg); np = world[-1]; x,fixtures = directions(world,origins)
    manifest = world[7]; tasks=[]; rows=[]
    origin_pin = write(out/'ORIGINS.json',origins); fixture_pin = write(out/'DIRECTIONS.json',fixtures)
    for name,f in fixtures.items():
        for metal in METALS:
            oldtask = next(t for t in manifest['tasks'] if t['task_id'] == metal+'_A')
            text = verify(oldtask['input']).read_text(); require(text.startswith('! '+manifest['method']+'\n'), 'Original method differs')
            template = text.replace(' VeryTightSCF EnGrad\n',' VeryTightSCF\n',1)
            require(template != text and 'EnGrad' not in template and 'DoEQ false' in template
                    and not re.search(r'%maxcore|%pal|\bOpt\b|CPCM|NumGrad',template,re.I), 'Unexpected new template option')
            for sign in ('plus','minus'):
                tid=f'{name}_{metal}_{sign}'; d=out/tid;d.mkdir();h=SETTINGS['h_A']*(1 if sign=='plus' else -1)
                y=x.copy();y[f['index']]+=h*np.asarray(f['direction'])
                require(abs(float(np.linalg.norm(y-x,axis=1).max())-abs(h))<1e-12,'Displacement length differs')
                maps=update_geometry(world,origins[metal],y,d,f);(d/'endpoint.inp').write_text(template)
                task=dict(task_id=tid,metal=metal,charge=0,multiplicity=1,physical_multiplicity=oldtask['physical_multiplicity'],
                    electronic_state=oldtask['electronic_state'],task_type='single_point',input=pin(d/'endpoint.inp'),
                    output_path=str(d/'endpoint.out'),origin_task=oldtask,fixture=name,signed_h_A=h,**maps)
                task['cache_key']=hashlib.sha256(json.dumps({'protocol':PROTOCOL,'state':task['electronic_state'],
                    'input':task['input']['sha256'],'xyz':task['xyz']['sha256'],'field':task['pointcharges']['sha256'],
                    'physical':task['physical']['sha256'],'orca':cfg['orca'],'basis_ECP':cfg['basis_assets'],
                    'native_helpers':cfg['native_implementation'],'MM_template':cfg['MM_template'],
                    'origin_MM_ledger':read(origins[metal]['classical_reused'])['ledger'],
                    'ALPB_binary':cfg['ALPB_binary'],'ALPB_library_directory':cfg['ALPB_library_directory'],
                    'ALPB_source_charge_vector':origins[metal]['solvent_charges'],'helpers':cfg['helpers']},sort_keys=True).encode()).hexdigest()
                rows.append(write(d/'TASK.json',task));tasks.append(task)
    require(tuple(t['task_id'] for t in tasks)==IDS,'Eight-state denominator differs')
    for name in fixtures:
        for sign in ('plus','minus'):
            pair=[t for t in tasks if t['fixture']==name and t['task_id'].endswith('_'+sign)]
            require(np.array_equal([a['xyz_A'] for a in read(pair[0]['physical'])],[a['xyz_A'] for a in read(pair[1]['physical'])]),'Paired physical geometry differs')
    header={'protocol_id':PROTOCOL,'agreement':cfg['specification'],'plan':cfg['plan'],'source_origin_result':cfg['origin_result'],
        'method':manifest['method'].removesuffix(' EnGrad'),'orca':cfg['orca'],'assets':cfg['basis_assets'],
        'implementation':cfg['native_implementation'],'execution_policy':{k:cfg['native_implementation'][v] for k,v in
            [('task_runner','run_orca_task_manifest.py'),('runtime_renderer','render_orca_runtime_input.py')]},
        'execution_resources':{'mpi_ranks':12,'concurrent_tasks':1},'tasks':tasks}
    result.update(status='prepared',origins=origin_pin,directions=fixture_pin,tasks=rows,
        full_manifest=write(out/'FULL_MANIFEST.json',header),origin_replay='passed',new_molecular_calls=0)


def prepared(cfg):
    p = read(cfg['preparation_result']); pcfg = read(p['config'])
    require(p['status'] == 'prepared' and p['protocol'] == PROTOCOL and p['new_molecular_calls'] == 0
            and p['origin_replay'] == 'passed', 'Actual zero-call preparation required')
    for k in ('protocol','settings','implementation','plan','specification','origin_result','python','helpers',
              'native_implementation','orca','basis_assets','MM_template','ALPB_binary','ALPB_library_directory'):
        require(pcfg[k] == cfg[k], 'Prepared recipe differs: '+k)
    m = read(p['full_manifest']); tasks = {read(q)['task_id']:read(q) for q in p['tasks']}
    require(tuple(tasks) == IDS and list(tasks.values()) == m['tasks'], 'Prepared eight-task inventory differs')
    for t in tasks.values():
        for k in ('physical','source_inputs','mapping','xyz','pointcharges','core_mapping','boundary_mapping','input'): verify(t[k])
    return p,m,tasks


def admit_sp(manifest_pin, cfg, native):
    manifest, rows = native['ggr_sensitivity'].executed(verify(manifest_pin))
    require(len(manifest['tasks']) == 1, 'Execution manifest must contain exactly one actual task')
    task = manifest['tasks'][0]; row = rows[task['task_id']]; row.update(task=task)
    if row['status'] != 'complete': return row
    text = verify(row['output']).read_text(); state = task['electronic_state']
    for pat,wanted in [(r'Program Version\s+(6\.1\.1)\b','6.1.1'),
        (r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)','0'), (r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)','1'),
        (r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)','832'), (r'Hartree-Fock type\s+HFTyp\s+\.{2,}\s+(\w+)','RHF')]:
        require(re.findall(pat,text) == [wanted], 'Actual native state/version differs: '+pat)
    ecp = re.findall(r'Type\s+(\w+)\s+ECP(?:\s+(\S+))?\s+\(replacing\s+(\d+)\s+core electrons',text)
    require(ecp == [(task['metal'],'',str(state['ecp_core_electrons']))], 'Executed ECP differs')
    counts = re.findall(r'Reading point charge file\s+\.{2,}\s+ok\s+\((\d+) point charges\)',text)
    require(counts and all(v == '1840' for v in counts) and 'environment.pc' in text
            and 'DFTD4' in text and 'INITIAL GUESS: MODEL POTENTIAL' in text, 'Executed field/dispersion/guess differs')
    require(not re.search(r'gCP correction\s+[-+0-9.]|^\s*(?:CPCM|SMD|COSMO) SOLVATION MODEL',text,re.M), 'Undeclared energy component')
    obs = module('mex_dft_observation_helper',cfg['helpers']['reference_runtime']).observations(row['output']['path'],native['mace_omol_vacuum'])
    require(obs['ConvCheckMode'] == 'Total+1el-Energy' and obs['ConvForced'] == 0, 'Native convergence semantics differ')
    values = re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+0-9.Ee]+)',text)
    require(len(values) == 1 and float(values[0]) == row['energy_hartree'], 'Native energy text differs')
    receipt = read(row['receipt']); sidecar = read(receipt['artifacts']['runtime_sidecar'])
    require(receipt['orca_executable'] == cfg['orca'] and receipt['returncode'] == 0
            and receipt['normal_termination'] and receipt['scf_converged'], 'Native receipt failed')
    policy=sidecar['memory_policy']
    require(policy['workers']==1 and policy['ranks_per_worker']==12 and policy['memory_source']=='SLURM_MEM_PER_NODE',
            'Native allocation-derived memory policy differs')
    row.update(status='available',executed_ecp=canonical(ecp),SCF=obs,
        printed_energy=values[0],printed_quantum_hartree=float(Decimal(1).scaleb(Decimal(values[0]).as_tuple().exponent)),
        runtime_memory_policy=sidecar['memory_policy'],
        molecular_wall_seconds=(dt.datetime.fromisoformat(receipt['finished_at_utc'])-dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds(),
        physical_state_qualified=False)
    return row


def singlepoint(cfg, out, result, tid):
    require(tid in IDS, 'Unknown finite task')
    original,ref,world,origins,native = inputs(cfg); np = world[-1]; prep,full,tasks = prepared(cfg)
    if tid not in SCOUT:
        scout = read(cfg['scout_collection'])
        require(scout['status'] == 'complete' and scout['stage'] == 'collect' and scout['scope'] == 'scout_pair'
                and scout['scout_pair_passed'] and scout['preparation_result'] == cfg['preparation_result']
                and scout['selected_tasks'] == list(SCOUT) and scout['settings'] == SETTINGS,
                'Actual admitted La external-charge pair required before remaining six')
    task = copy.deepcopy(tasks[tid]); origin = origins[task['metal']]
    for k,name in [('xyz','core.xyz'),('pointcharges','environment.pc'),('input','endpoint.inp')]:
        (out/name).write_bytes(verify(task[k]).read_bytes()); task[k] = pin(out/name)
    task['prepared_task'] = next(p for p in prep['tasks'] if read(p)['task_id'] == tid)
    task['output_path'] = str(out/'endpoint.out'); task_pin=write(out/'TASK.json',task)
    manifest = dict(full,tasks=[task],full_prepared_manifest=prep['full_manifest'],preparation_result=cfg['preparation_result'])
    manifest_pin = write(out/'MANIFEST.json',manifest); result.update(task_id=tid,task=task_pin,manifest=manifest_pin,
        preparation_result=cfg['preparation_result'],selected_tasks=[tid],electronic_state=task['electronic_state'])
    runtime_helper = module('mex_dft_runtime',cfg['helpers']['reference_runtime'])
    result['allocation'] = runtime_helper.allocation(native['render_orca_runtime_input'],out)
    # Use the existing one-cell MM evaluator, including its original exclusions.
    mmcfg = read(cfg['MM_template']); mmcfg.update(trials=[task_pin],plan=cfg['specification'])
    require(mmcfg['ledgers'] == {m:read(origins[m]['classical_reused'])['ledger'] for m in METALS}, 'MM ledger differs')
    mmout=out/'MM'; mmout.mkdir(); mmresult={}; mmhelper=module('mex_dft_mm',cfg['helpers']['MM'])
    result['MM_config']=write(out/'MM_CONFIG.json',mmcfg)
    try: mmhelper.mm_stage(mmcfg,mmout,mmresult)
    finally:
        result['MM']=write(mmout/'RESULT.json',mmresult)
        result['new_MM_configurations_attempted']=mmresult.get('configurations_attempted',0)
        result['new_MM_group_queries_attempted']=mmresult.get('group_queries_attempted',0)
    require(mmresult.get('status') == 'complete' and len(mmresult['rows']) == 1, 'MM unavailable')
    mmrow=mmresult['rows'][0]; require(mmrow['task_id']==tid and mmrow['metal']==task['metal'], 'MM state differs')
    # Same native global solvent executable and exact unshifted physical charges.
    phys=read(task['physical']); x=np.asarray([a['xyz_A'] for a in phys]); q=np.loadtxt(verify(origin['solvent_charges']))
    require(x.shape==(2049,3) and q.shape==(2049,) and np.isfinite(x).all() and np.isfinite(q).all(), 'Current solvent inventory differs')
    env=dict(os.environ,LD_LIBRARY_PATH=cfg['ALPB_library_directory'],OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
    static=module('mex_dft_solvent',cfg['helpers']['ALPB']); verify(cfg['ALPB_binary'])
    result['new_ALPB_calls_attempted']+=1
    sol=static.solvent_cell(verify(cfg['ALPB_binary']),{'physical':phys,'x':x,'q':q,'mapping_pin':task['mapping']},out/'ALPB',env,np)
    result['ALPB']=pin(out/'ALPB/RESULT.json'); result['ALPB_stdout']=pin(out/'ALPB/stdout.txt')
    for k in ('xyz','pointcharges','input'): verify(task[k])
    write(out/'START.json',result)
    try:
        result['new_DFT_calls_requested']=1
        native['run_orca_task_manifest'].run_manifest(verify(manifest_pin),orca_path=verify(cfg['orca']),workers=1,nprocs=12,
            selected_tasks={tid},expected_runner_sha256=cfg['native_implementation']['run_orca_task_manifest.py']['sha256'])
    except Exception as exc: result['runner_failure']=repr(exc)
    finally:
        result['new_DFT_calls_attempted']=int((out/'endpoint.out').exists())
        for k in ('xyz','pointcharges','input'): verify(task[k])
        row=admit_sp(manifest_pin,cfg,native); result['DFT']=row
    if row['status']=='available':
        result.update(status='complete',native_energies={'PBE0_embedded_hartree':row['energy_hartree'],
            'MM_kcal_mol':mmrow['total_kcal_mol'],'MM_groups_kcal_mol':mmrow['components_kcal_mol'],'ALPB_hartree':sol['total_hartree']})
    else: result['status']='unavailable'


def compare_pair(cfg, origins, fixtures, rows, name, metal):
    import numpy as np
    r={'fixture':name,'metal':metal,'status':'unavailable'}
    p=rows[f'{name}_{metal}_plus'];m=rows[f'{name}_{metal}_minus']
    if p['status']!='complete' or m['status']!='complete': return r
    f=fixtures[name]; v=np.asarray(f['direction']); origin=origins[metal]; h=SETTINGS['h_A']; den=2*h
    analytic={k:float(np.dot(np.asarray(origin['component_gradients_kcal_mol_A'][k])[f['index']],v)) for k in COMPONENTS}
    analytic['total']=float(np.dot(np.asarray(origin['gradient_kcal_mol_A'])[f['index']],v))
    ep,em=p['native_energies'],m['native_energies']
    fd={'PBE0_embedded':float((Decimal(p['DFT']['printed_energy'])-Decimal(m['DFT']['printed_energy']))*Decimal(str(HA))/Decimal(str(den))),
        'MM_ledger':(ep['MM_kcal_mol']-em['MM_kcal_mol'])/den,
        'fixed_FF_charge_ALPB':(ep['ALPB_hartree']-em['ALPB_hartree'])*HA/den}
    fd['total']=math.fsum(fd.values()); checks={}
    for k in (*COMPONENTS,'total'):
        tol=SETTINGS['absolute_kcal_mol_A']+SETTINGS['relative']*max(abs(fd[k]),abs(analytic[k]))
        checks[k]={'analytic_kcal_mol_A':analytic[k],'FD_kcal_mol_A':fd[k],
            'residual_kcal_mol_A':fd[k]-analytic[k],'tolerance_kcal_mol_A':tol,'passed':abs(fd[k]-analytic[k])<=tol}
    classical=read(origin['classical_reused']); ci=classical['source_ids'].index(f['id']); groupchecks={}
    for k,g in classical['component_gradients_kcal_mol_A'].items():
        a=float(np.dot(np.asarray(g)[ci],v)); d=(ep['MM_groups_kcal_mol'][k]-em['MM_groups_kcal_mol'][k])/den
        groupchecks[k]={'analytic_kcal_mol_A':a,'FD_kcal_mol_A':d,'residual_kcal_mol_A':d-a}
    def energy_residual(row):
        found=[v for k,v in row['DFT']['SCF']['terminal_residuals'].items() if k=='Energy change']
        require(len(found)==1,'Printed SCF energy residual key differs');return found[0]
    rp,rm=energy_residual(p),energy_residual(m)
    precision={'printed_energy_half_quantum_sum_kcal_mol_A':.5*(p['DFT']['printed_quantum_hartree']+m['DFT']['printed_quantum_hartree'])*HA/den,
        'float64_operand_ulp_sum_kcal_mol_A':(math.ulp(ep['PBE0_embedded_hartree'])+math.ulp(em['PBE0_embedded_hartree']))*HA/den,
        'SCF_energy_tolerance_sum_sensitivity_proxy_kcal_mol_A':(rp['printed_tolerance']+rm['printed_tolerance'])*HA/den,
        'SCF_last_energy_change_sum_sensitivity_proxy_kcal_mol_A':(abs(rp['value'])+abs(rm['value']))*HA/den,
        'SCF_proxies_are_error_bounds':False,'truncation_or_step_halving_qualified':False,
        'DFT_subtraction':'Decimal subtraction of exact printed energy strings; native MM/ALPB differences before summation'}
    r.update(status='available',passed=all(c['passed'] for c in checks.values()),checks=checks,MM_groups=groupchecks,
        numerical_precision=precision,origin_SCF=origin['SCF'],plus_SCF=p['DFT']['SCF'],minus_SCF=m['DFT']['SCF'])
    return r


def collect(cfg,out,result):
    original,ref,world,origins,native=inputs(cfg); prep,manifest,tasks=prepared(cfg)
    scope=cfg['collection_scope'];require(scope in ('scout_pair','all_eight'),'Unsupported collection scope')
    selected=list(SCOUT if scope=='scout_pair' else IDS); endpoint_pins=cfg['endpoint_results']
    require(set(endpoint_pins)==set(selected),'Declared collection denominator differs')
    rows={tid:{'task_id':tid,'status':'not_requested' if tid not in selected else 'unavailable'} for tid in IDS}
    result.update(scope=scope,selected_tasks=selected,preparation_result=cfg['preparation_result'],rows=rows,checks=[])
    for tid in selected:
        row=rows[tid];row['result']=endpoint_pins[tid]
        try:
            parent=read(endpoint_pins[tid]);pcfg=read(parent['config'])
            require(parent['status']=='complete' and parent['stage']=='singlepoint' and parent['task_id']==tid
                and parent['preparation_result']==cfg['preparation_result'] and parent['implementation']['sha256']==cfg['implementation']['sha256'], 'Endpoint execution unavailable/different')
            for k in ('protocol','settings','implementation','plan','specification','origin_result','python','helpers','native_implementation','orca','basis_assets','MM_template','ALPB_binary','ALPB_library_directory'):
                require(pcfg[k]==cfg[k],'Endpoint recipe differs: '+k)
            task=read(parent['task']); expected=tasks[tid]
            require(task['prepared_task']==next(p for p in prep['tasks'] if read(p)['task_id']==tid), 'Prepared endpoint identity differs')
            for k in ('xyz','input','pointcharges'):
                verify(task[k]);require(task[k]['sha256']==expected[k]['sha256'],'Executed input bytes differ: '+k)
            for k in ('physical','mapping','source_inputs','core_mapping','boundary_mapping','electronic_state','cache_key','fixture','signed_h_A'):
                require(task[k]==expected[k],'Current map/state differs: '+k)
            dft=admit_sp(parent['manifest'],cfg,native);require(canonical(dft)==parent['DFT'] and dft['status']=='available'
                and dft['task']==task,'Raw DFT replay differs')
            mm=read(parent['MM']);mc=read(parent['MM_config']);sol=read(parent['ALPB'])
            require(mm['status']=='complete' and len(mm['rows'])==1 and mc['trials']==[parent['task']], 'MM cell binding differs')
            mr=mm['rows'][0];require(mr['task_id']==tid and mr['metal']==task['metal'] and mr['status']=='complete','MM state differs')
            require(mc['ledgers']=={m:read(origins[m]['classical_reused'])['ledger'] for m in METALS},'MM ledger changed')
            for k in ('gradient','xyz','charges','execution'):verify(sol[k])
            ex=read(sol['execution']);require(ex['returncode']==0 and ex['command']==[cfg['ALPB_binary']['path'],sol['xyz']['path'],sol['charges']['path']], 'Native solvent execution differs')
            stdout=verify(parent['ALPB_stdout']).read_text()
            require(verify(parent['ALPB_stdout']).parent==verify(sol['execution']).parent,'Solvent stdout location differs')
            parts_lines=re.findall(r'^COMPONENTS\s+(.+)$',stdout,re.M)
            require(len(parts_lines)==1 and dict(zip(('Gelec','Ghb','Gsasa','Gshift'),map(float,parts_lines[0].split())))==sol['components_hartree'], 'Raw solvent energy differs')
            require(sol['charges']['sha256']==origins[task['metal']]['solvent_charges']['sha256'],'Physical FF charge vector differs')
            # Full current solvent coordinate identity, not only the core.
            raw=verify(sol['xyz']).read_text().splitlines(); parts=[s.split() for s in raw[2:]]; physical=read(task['physical']);np=world[-1]
            require(int(raw[0])==len(parts)==2049 and [s[0] for s in parts]==[a['element'] for a in physical]
                and np.array_equal(np.asarray([s[1:] for s in parts],float),[a['xyz_A'] for a in physical]),'Solvent physical coordinates differ')
            native_values={'PBE0_embedded_hartree':dft['energy_hartree'],'MM_kcal_mol':mr['total_kcal_mol'],
                'MM_groups_kcal_mol':mr['components_kcal_mol'],'ALPB_hartree':sol['total_hartree']}
            require(native_values==parent['native_energies'] and math.fsum(mr['components_kcal_mol'].values())==mr['total_kcal_mol']
                and math.fsum(sol['components_hartree'].values())==sol['total_hartree'],'Native component closure differs')
            row.update(status='complete',DFT=dft,native_energies=native_values,MM=parent['MM'],ALPB=parent['ALPB'])
        except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
    fixtures=read(prep['directions'])
    for name in SETTINGS['directions']:
        for metal in METALS:
            if f'{name}_{metal}_plus' in selected: result['checks'].append(compare_pair(cfg,origins,fixtures,rows,name,metal))
    result['common_and_differential']=[]
    for name in SETTINGS['directions']:
        pair={r['metal']:r for r in result['checks'] if r['fixture']==name and r['status']=='available'}
        if set(pair)!=set(METALS): continue
        item={'fixture':name,'values':{}}
        for k in (*COMPONENTS,'total'):
            a,b=pair['La']['checks'][k],pair['Dy']['checks'][k]
            item['values'][k]={label:{v:(a[v]+b[v])/2 if label=='common' else b[v]-a[v]
                for v in ('analytic_kcal_mol_A','FD_kcal_mol_A','residual_kcal_mol_A')} for label in ('common','Dy_minus_La')}
        result['common_and_differential'].append(item)
    complete=all(rows[t]['status']=='complete' for t in selected)
    passed=complete and all(r.get('passed',False) for r in result['checks'])
    result.update(status='complete' if passed else ('failed_derivative_gate' if complete else 'unavailable'),
        scout_pair_passed=passed if scope=='scout_pair' else None,
        all_eight_directional_derivatives_passed=passed if scope=='all_eight' else None,
        physical_state_qualified=False,global_force_qualified=False,SCF_error_bound_available=False,
        new_molecular_calls=0,interpretation='Directional numerical consistency only; unresolved density/state and single-step truncation uncertainty remain.')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stage',choices=('prepare','singlepoint','collect'),required=True)
    ap.add_argument('--task',choices=IDS);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required');out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);tick=time.monotonic()
    cfg=json.loads(a.config.read_text());result={'protocol':PROTOCOL,'stage':a.stage,'status':'execution_failure',
        'implementation':pin(__file__),'config':pin(a.config),'settings':SETTINGS,'slurm_job_id':os.environ['SLURM_JOB_ID'],
        'declared_displaced_SP':list(IDS),'new_DFT_calls_requested':0,'new_DFT_calls_attempted':0,
        'new_MM_configurations_attempted':0,'new_MM_group_queries_attempted':0,'new_ALPB_calls_attempted':0,
        'new_central_gradients':0,'new_MACE_calls':0,'new_GFN_calls':0,'optimization_steps':0}
    try:
        if a.stage=='prepare':prepare(cfg,out,result)
        elif a.stage=='singlepoint':singlepoint(cfg,out,result,a.task)
        else:collect(cfg,out,result)
    except Exception as exc:result.update(failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-tick
    result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ.get('SLURM_CPUS_ON_NODE','1'))
    result['parent_peak_RSS_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['parent_RSS_is_not_MPI_peak']=True;write(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')}));return 0 if result['status'] in ('prepared','complete') else 1


if __name__=='__main__':raise SystemExit(main())
