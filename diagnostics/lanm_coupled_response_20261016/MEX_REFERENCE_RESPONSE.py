#!/usr/bin/env python3
"""Slurm saved-data comparison of six fixed Mex electronic-response directions.

CLI --config CONFIG --output FRESH. No molecular engine is called.
Config: repository, implementation, plan, reference_collection, static_result,
dependencies {each DEPENDENCIES filename: {path,sha256}}. Root binds static_result
only after actual admission; this adapter does not generate a replacement.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
PROTOCOL = 'mex_saved_pbe0_six_direction_response_v1'
SOURCE_SHA = '4b3b98983af0cb5fcdbf7b8f52e35063308201ce9df777a30484acf5146b7930'
METHOD = 'PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad'
DEPENDENCIES = ('metal_environment_force_assembly.py', 'affordable_common.py',
    'affordable_response.py', 'metal_environment_reference.py',
    'metal_environment_electronic_state.py', 'affordable_workflow.py',
    'mace_omol_vacuum.py', 'mace_file_checks.py', 'mace_site_kinematics.py')
COMPONENTS = ('MACE_vacuum', 'GFN2_field_difference', 'MM_ledger', 'fixed_FF_charge_ALPB')
METALS = ('La', 'Dy')
HA, BOHR, EV = 627.509474, .529177210903, 23.06054783061903


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Missing/changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def close(np, actual, expected, label):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
            and np.allclose(a, b, atol=1e-10, rtol=1e-12), label + ' differs')
    return float(np.max(np.abs(a-b)))


def inputs(config):
    cfg = json.loads(config.read_text())
    require(Path(cfg['repository']).resolve() == ROOT, 'Repository differs')
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation differs')
    verify(cfg['plan'])
    require(set(cfg['dependencies']) == set(DEPENDENCIES), 'Dependency inventory differs')
    for name, p in cfg['dependencies'].items():
        require(verify(p).resolve() == ROOT/'scripts'/name, 'Dependency location differs: ' + name)
    sys.path.insert(0, str(ROOT/'scripts'))
    import numpy as np
    import metal_environment_force_assembly as mapper
    from mace_site_kinematics import Kinematics
    require(mapper.HA_TO_KCAL == HA and mapper.BOHR_TO_A == BOHR, 'Unit constants differ')

    static = read(cfg['static_result']); scfg = read(static['config'])
    require(static['status'] == 'complete' and static['protocol_id'] ==
            'mex_static_local_MACE_field_difference_MM_fixedFF_ALPB_v1', 'Admitted static model required')
    require(static['implementation']['sha256'] == scfg['implementation']['sha256']
            and static['plan']['sha256'] == scfg['plan']['sha256'],
            'Static source/plan binding differs')
    for p in (static['implementation'], scfg['implementation'], static['plan'], scfg['plan']):
        verify(p)
    require(static['new_MACE_calls_attempted'] == static['native_model_calls'] == 2
            and static['new_ALPB_calls_attempted'] == 2 and static['new_SCC_calls'] == static['new_MM_calls'] == 0,
            'Static call inventory differs')
    field = read(scfg['field_result']); fcfg = read(field['config'])
    prep = read(field['preparation_result']); source_pin = prep['source_inputs']; source = read(source_pin)
    require(field['status'] == 'complete' and field['slurm_job_id'] == '12570'
            and field['attempted'] == field['accepted'] == 4, 'Local field origin unavailable')
    require(prep['status'] == 'prepared' and prep['slurm_job_id'] == '12569'
            and source_pin['sha256'] == SOURCE_SHA and fcfg['inputs'] == source_pin,
            'Repaired source binding differs')
    require(fcfg['force_mapping_helper']['sha256'] == cfg['dependencies']['metal_environment_force_assembly.py']['sha256'],
            'Static force-mapping helper differs')
    for obj in (field, prep):
        verify(obj['implementation']); verify(obj['plan'])
    mf = read(scfg['mace_force_result']); mcfg = read(mf['config'])
    require(static['model'] == mcfg['model'] and static['versions'] == mcfg['expected_versions'],
            'Static MACE model/version differs')
    verify(mcfg['model']['checkpoint']); verify(scfg['native_helper'])

    ref = read(cfg['reference_collection']); manifest = read(ref['manifest'])
    require(ref['source_id'] == 'Mex_8FNS' and set(ref['rows']) == {'La_A', 'Dy_A'}
            and manifest['inputs'] == source_pin and manifest['method'] == METHOD,
            'Archived reference source/method differs')
    require(ref['force_consistency_qualified'] is False, 'Archived reference qualification changed')
    verify(ref['collector']); verify(manifest['orca']); verify(manifest['agreement'])
    for assets in manifest['assets'].values():
        for p in assets.values():
            verify(p)
    export = read(source['export']); atoms = read(export['atoms']); bonds = read(export['bonds'])
    verify(export['system'])
    require(export['atom_count'] == len(atoms) == 2045, 'Native source inventory differs')
    require(len(source['spectator_occupancy']) == 3 and all(a['element'] == 'Nd' and a['charge_e'] == 3
            for a in source['spectator_occupancy']), 'Spectator occupancy differs')
    rows = {r['metal']: r for r in static['rows']}
    require(len(static['rows']) == len(rows) == 2 and set(rows) == set(METALS), 'Static pair incomplete')
    fields = {(r['task']['metal'], r['task']['kind']): r for r in field['rows']}
    require(len(fields) == 4, 'Duplicate local field cells')
    return cfg, static, scfg, prep, source_pin, source, ref, manifest, atoms, bonds, rows, fields, mapper, Kinematics, np


def gradients(world, metal):
    cfg, static, scfg, prep, source_pin, source, ref, manifest, atoms, bonds, rows, fields, mapper, _, np = world
    row = rows[metal]
    require(row['status'] == 'complete' and read(row['result']) == {k:v for k,v in row.items() if k != 'result'},
            'Static row receipt differs')
    assembled = read(row['assembled']); mp = read(assembled['physical_mapping'])
    require(assembled['physical_mapping'] == prep['mappings'][metal] and assembled['metal'] == metal
            and mp['inputs'] == source_pin and mp['metal'] == metal, 'Physical mapping binding differs')
    physical = mp['physical_atoms']; ids = [a['id'] for a in physical]
    require(len(ids) == len(set(ids)) == 2049 and len(mp['core']) == 209 and len(mp['field']) == 1840
            and len(mp['caps']) == 4 and assembled['source_ids'] == ids
            and assembled['physical_atoms'] == physical, 'Physical map/inventory differs')
    endpoint = source['configurations']['A']['endpoints'][metal]
    require(assembled['source_xyz'] == endpoint['xyz'] and mp['xyz'] == endpoint['xyz'], 'Static local source differs')
    verify(endpoint['xyz']); verify(mp['pointcharges'])
    task = next(t for t in manifest['tasks'] if t['task_id'] == metal+'_A')
    r = ref['rows'][metal+'_A']; receipt = read(r['receipt'])
    require(r['status'] == 'complete' and receipt['returncode'] == 0 and receipt['normal_termination']
            and receipt['scf_converged'] and receipt['expected_engrad_complete']
            and receipt['manifest'] == ref['manifest'] and receipt['task_id'] == task['task_id']
            and receipt['orca_version'] == '6.1.1', 'Archived reference execution unavailable')
    for key in ('input', 'xyz', 'pointcharges', 'core_mapping', 'boundary_mapping'):
        verify(task[key])
    template = verify(task['input']).read_text()
    require(template.startswith('! '+METHOD+'\n') and 'DoEQ false' in template
            and '* xyzfile 0 1 core.xyz' in template, 'Archived executed recipe differs')
    for key in ('engrad', 'pointcharge_gradient', 'output'):
        verify(r[key])
    for p in receipt['artifacts'].values():
        verify(p)
    require(task['xyz']['sha256'] == endpoint['xyz']['sha256']
            and task['pointcharges']['sha256'] == mp['pointcharges']['sha256']
            and task['core_mapping'] == source['configurations']['A']['core_mapping']
            and task['boundary_mapping'] == source['boundary_mapping'], 'Reference/static geometry or field differs')
    state = r['electronic_state']; physical_mult = 1 if metal == 'La' else 6; core = 46 if metal == 'La' else 55
    require(r['quantity'] == 'gradient_not_force' and r['pointcharge_gradient_units'] == 'Hartree/bohr',
            'Reference units/sign differ')
    require(state == task['electronic_state'] and state['metal'] == metal and state['charge'] == 0
            and state['effective_multiplicity'] == task['multiplicity'] == 1
            and state['physical_multiplicity'] == task['physical_multiplicity'] == physical_mult
            and state['ecp_core_electrons'] == core and state['explicit_electrons'] == 832
            and state['representation'] == 'spin_free_frozen_f_valence_model'
            and state['spin_orbit_included'] is False and state['physical_state_qualified'] is False
            and r['physical_state_qualified'] is False and r['force_consistency_qualified'] is False,
            'Reference electronic-state distinction differs')
    require(r['executed_ecp'] == [[metal, '', str(core)]], 'Executed reference ECP differs')
    reference = mapper.assemble(verify(source_pin), 'A', metal, verify(cfg['reference_collection']), metal+'_A')
    require(reference['physical_atoms'] == physical, 'Reference physical ordering/coordinates differ')
    rg = np.asarray(reference['electronic']['gradient_kcal_mol_A'])
    components = {k: np.asarray(v) for k,v in assembled['component_gradients_kcal_mol_A'].items()}
    require(set(components) == set(COMPONENTS), 'Static component inventory differs')
    residuals = {}
    for k, g in components.items():
        require(g.shape == (2049,3) and np.isfinite(g).all(), 'Invalid static component: '+k)
    residuals['static_total'] = close(np, sum(components.values()), assembled['gradient_kcal_mol_A'], 'Static sum')
    residuals['static_force_sign'] = close(np, -np.asarray(assembled['gradient_kcal_mol_A']), assembled['force_kcal_mol_A'], 'Force sign')
    require(assembled['MACE'] == row['MACE'] and assembled['ALPB'] == row['ALPB'], 'Static raw component binding differs')
    ml = read(assembled['MACE']); check = ml['input_state_check']
    require(ml['xyz'] == endpoint['xyz'] and ml['charge'] == 0 and ml['physical_multiplicity'] == physical_mult
            and check['atoms'] == 209 and check['charge'] == 0 and check['spin_multiplicity'] == physical_mult
            and check['selected_metal_index'] == 0 and check['head'] == 'omol' and check['graph_count'] == 1,
            'Native MACE physical state differs')
    require(ml['force_convention'] == 'negative vacuum energy derivative', 'MACE force convention differs')
    f = np.asarray(ml['forces_eV_A']); require(f.shape == (209,3) and np.isfinite(f).all(), 'Invalid MACE forces')
    mapped = mapper.scatter(mp, -f*EV*BOHR/HA, np.zeros((1840,3)))
    residuals['MACE_cap_mapping'] = close(np, mapped['gradient_kcal_mol_A'], components['MACE_vacuum'], 'MACE physical replay')
    fg = {}
    for kind in ('embedded', 'isolated'):
        fieldrow = fields[(metal, kind)]; saved = read(assembled['GFN2_reused'][kind])
        require(assembled['GFN2_reused'][kind] == fieldrow['result']
                and saved == {k:v for k,v in fieldrow.items() if k != 'result'}
                and saved['status'] == 'accepted', 'GFN2 reused row differs')
        t = saved['task']
        require(t['charge'] == 0 and t['multiplicity'] == 1 and t['uhf'] == 0
                and t['physical_multiplicity'] == physical_mult
                and t['metal'] == metal and t['xyz']['sha256'] == endpoint['xyz']['sha256']
                and t['mapping'] == assembled['physical_mapping'], 'GFN2 source/effective state differs')
        for p in saved['artifacts'].values():
            verify(p)
        fg[kind] = np.asarray(saved['mapped_physical_gradient']['gradient_kcal_mol_A'])
        require(fg[kind].shape == (2049,3) and np.isfinite(fg[kind]).all(), 'Invalid field physical gradient')
    residuals['field_difference'] = close(np, fg['embedded']-fg['isolated'], components['GFN2_field_difference'], 'Field difference')
    require(assembled['classical_reused'] == scfg['classical'][metal], 'Classical binding differs')
    classical = read(assembled['classical_reused']); ci = {sid:i for i,sid in enumerate(classical['source_ids'])}
    require(classical['inputs'] == source_pin and classical['metal'] == metal and set(ci) == set(ids), 'Classical state/map differs')
    residuals['MM'] = close(np, np.asarray(classical['gradient_kcal_mol_A'])[[ci[sid] for sid in ids]], components['MM_ledger'], 'MM replay')
    solvent = read(assembled['ALPB']); verify(solvent['gradient']); verify(solvent['charges']); verify(solvent['xyz'])
    residuals['ALPB'] = close(np, np.asarray(solvent['gradient_hartree_per_bohr'])*HA/BOHR, components['fixed_FF_charge_ALPB'], 'ALPB replay')
    approx = components['MACE_vacuum'] + components['GFN2_field_difference']
    shared = components['MM_ledger'] + components['fixed_FF_charge_ALPB']
    arrays = dict(components, PBE0_embedded=rg, approximation_electronic=approx,
                  approximation_minus_PBE0=approx-rg, approximation_total=approx+shared,
                  PBE0_plus_shared_MM_ALPB=rg+shared)
    evidence = {'reference_task': task, 'reference_row': {k:r[k] for k in ('receipt','engrad','pointcharge_gradient','output','electronic_state')},
                'assembled': row['assembled'], 'replay_residuals_kcal_mol_A': residuals,
                'reference_quantity': 'gradient=dE/dR; force=-gradient', 'reference_physical_state_qualified': False,
                'reference_force_consistency_qualified': False, 'static_total_derivative_FD_qualified': assembled['total_derivative_FD_qualified']}
    return mp, arrays, evidence


def directions(mp, atoms, bonds, Kinematics, np):
    ids = [a['id'] for a in mp['physical_atoms']]; ix = {sid:i for i,sid in enumerate(ids)}
    x = np.asarray([a['xyz_A'] for a in mp['physical_atoms']]); metadata = {a['id']:a for a in atoms}
    graph = {sid:set() for sid in ids}
    for a,b in bonds:
        require(a in ix and b in ix, 'Native graph atom absent')
        graph[a].add(b); graph[b].add(a)
    targets = [a for a in mp['core'] if a['kind'] == 'metal']; require(len(targets) == 1, 'One local target required')
    t = targets[0]; s = t['source_site']; target = f"{s['chain']}/{s['resid']}/{s.get('insertion_code','')}/{s['name']}"
    require(t['qm_index'] == 0 and not graph[target], 'Target index/graph differs')
    modes = []
    for k, name in enumerate(('x','y','z')):
        axis = [0.,0.,0.]; axis[k] = 1.
        modes.append({'name':'target_'+name, 'unit':'angstrom', 'axis':axis,
                      'moving_indices':[ix[target]], 'moving_source_ids':[target]})
    for resid, resname, names in ((84,'ASP',('CB','HB2','HB3','CG','OD1','OD2')),
            (86,'ASP',('CB','HB2','HB3','CG','OD1','OD2')),
            (95,'GLU',('CB','HB2','HB3','CG','HG2','HG3','CD','OE1','OE2'))):
        a,b = f'A/{resid}//CA', f'A/{resid}//CB'
        require(b in graph[a] and metadata[a]['resname'] == metadata[b]['resname'] == resname, 'Donor axis differs')
        seen, pending = {b}, [b]
        while pending:
            current = pending.pop()
            for other in graph[current]-seen:
                if {current,other} != {a,b}:
                    seen.add(other); pending.append(other)
        require(a not in seen and seen == {f'A/{resid}//{n}' for n in names}, 'Donor graph component differs')
        moving = sorted(seen, key=ix.get)
        modes.append({'name':f'A{resid}_CA_CB', 'unit':'radian', 'axis_indices':[ix[a],ix[b]],
                      'axis_source_ids':[a,b], 'moving_indices':[ix[sid] for sid in moving], 'moving_source_ids':moving})
    links = []
    for a in mp['core']:
        if a['kind'] == 'metal':
            continue
        if a['kind'] == 'source':
            links.append([a['qm_index'],'source',ix[a['source_id']]])
        else:
            require(a['kind'] == 'cap', 'Unknown local mapping kind')
            links.append([a['qm_index'],'cap',ix[a['retained_source_id']],ix[a['omitted_source_id']],a['length_A'],a['xyz_A']])
    data = {'positions_A':x.tolist(), 'core_positions_A':[a['xyz_A'] for a in mp['core']],
            'modes':modes, 'core_links':links, 'metal_index':ix[target],
            'heavy_indices':[i for i,a in enumerate(mp['physical_atoms']) if a['element'] != 'H'],
            'physical_ids':ids, 'bonds':[[ix[a],ix[b]] for a,b in bonds]}
    physical, core, j, cj = Kinematics(data).evaluate(np.zeros(6))
    close(np, physical, x, 'Origin kinematic replay'); close(np, core, data['core_positions_A'], 'Origin core replay')
    scales = np.linalg.norm(j.reshape(6,-1), axis=1); require(np.all(scales > 0), 'Zero tangent')
    unit = j/scales[:,None,None]; gram = unit.reshape(6,-1)@unit.reshape(6,-1).T
    eig = np.linalg.eigvalsh(gram); require(np.all(eig > eig[-1]*1e-10), 'Six-direction span is rank deficient')
    description = {'physical_source_ids':ids, 'target_source_id':target, 'modes':modes,
        'column_norms_A_per_native_unit':scales.tolist(), 'normalized_Gram':gram.tolist(),
        'Gram_eigenvalues':eig.tolist(), 'rank':6, 'relative_rank_cutoff':1e-10,
        'Gram_condition_number':float(eig[-1]/eig[0]),
        'normalization':'unit Euclidean physical displacement over all 2049 atoms including H',
        'native_Jacobian_shape':[6,2049,3], 'finite_geometries_generated':0}
    return j, cj, unit, description


def run(config, out, result):
    world = inputs(config); cfg, static, _, _, _, _, _, manifest, atoms, bonds, _, _, _, Kinematics, np = world
    gradients_by_metal, mappings, evidence = {}, {}, {}
    for metal in METALS:
        mappings[metal], gradients_by_metal[metal], evidence[metal] = gradients(world, metal)
    a,b = mappings['La'],mappings['Dy']
    require([v['id'] for v in a['physical_atoms']] == [v['id'] for v in b['physical_atoms']]
            and np.array_equal([v['xyz_A'] for v in a['physical_atoms']], [v['xyz_A'] for v in b['physical_atoms']]),
            'Paired physical coordinates/order differ')
    j, cj, unit, description = directions(a, atoms, bonds, Kinematics, np)
    path = out/'JACOBIANS.npz'; np.savez_compressed(path, physical_native=j, core_native=cj, physical_unit=unit)
    description['jacobians'] = pin(path); write(out/'DIRECTIONS.json', description)
    all_arrays = dict(gradients_by_metal)
    all_arrays['common'] = {k:(gradients_by_metal['La'][k]+gradients_by_metal['Dy'][k])/2 for k in gradients_by_metal['La']}
    all_arrays['Dy_minus_La'] = {k:gradients_by_metal['Dy'][k]-gradients_by_metal['La'][k] for k in gradients_by_metal['La']}
    projections = {}
    for state, arrays in all_arrays.items():
        projections[state] = {}
        for component, g in arrays.items():
            load = np.einsum('mij,ij->m',unit,g)
            projections[state][component] = {'gradient_kcal_mol_A':load.tolist(), 'force_kcal_mol_A':(-load).tolist()}
    result.update(status='complete', static_result=cfg['static_result'], reference_collection=cfg['reference_collection'],
                  reference_manifest=world[6]['manifest'], reference_method=manifest['method'],
                  reference_basis_assets=manifest['assets'], static_model=static['model'], evidence=evidence,
                  directions=pin(out/'DIRECTIONS.json'), projections=projections,
                  mode_order=[m['name'] for m in description['modes']], dependencies=cfg['dependencies'],
                  interpretation='Local reference-response diagnostic; no physical accuracy pass threshold, affinity or classifier')
    lines = ['# Mex saved electronic-response comparison', '',
        'Signed energy-gradient projections in kcal/mol/Å; force loads have the opposite sign.', '',
        '| Direction | PBE0 La | Approx La | PBE0 Dy | Approx Dy | PBE0 Dy−La | Approx Dy−La |',
        '|---|---:|---:|---:|---:|---:|---:|']
    def value(state, method, i):
        return projections[state][method]['gradient_kcal_mol_A'][i]
    for i,name in enumerate(result['mode_order']):
        values = [value(s,m,i) for s in ('La','Dy','Dy_minus_La') for m in ('PBE0_embedded','approximation_electronic')]
        lines.append('| '+name+' | '+' | '.join(f'{v:+.6f}' for v in values)+' |')
    lines += ['', '| Direction | PBE0 common | Approx common | Approx−PBE0 common | Approx−PBE0 Dy−La |',
              '|---|---:|---:|---:|---:|']
    for i,name in enumerate(result['mode_order']):
        values = [value('common',m,i) for m in ('PBE0_embedded','approximation_electronic','approximation_minus_PBE0')]
        values.append(value('Dy_minus_La','approximation_minus_PBE0',i))
        lines.append('| '+name+' | '+' | '.join(f'{v:+.6f}' for v in values)+' |')
    lines += ['', 'Approximation = MACE vacuum + GFN2 embedded − isolated. Full own/common/differential',
        'and shared MM/ALPB loads are retained in RESULT.json. Directions were fixed before projections.', '',
        'PBE0 uses La ECP46/Dy ECP55 and effective restricted singlets; MACE uses physical',
        'multiplicities 1/6. GFN2 is a separate effective-singlet model. Archived physical-state',
        'and force-consistency qualification remain false; no spin–orbit or stable-state claim.',
        'Agreement is conditional local response evidence. Disagreement cannot uniquely separate',
        'MACE error from the GFN2 field approximation. No affinity, finite-energy prediction or',
        'gradient qualification follows. No molecular calls, new geometries or optimization.', '']
    (out/'REPORT.md').write_text('\n'.join(lines)); result['report'] = pin(out/'REPORT.md')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); require(os.environ.get('SLURM_JOB_ID','').isdigit(), 'Saved scientific analysis requires Slurm')
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    cfg = json.loads(a.config.read_text())
    result = {'protocol_id':PROTOCOL, 'status':'running', 'config':pin(a.config), 'implementation':pin(__file__),
        'plan':cfg['plan'], 'slurm_job_id':os.environ['SLURM_JOB_ID'], 'molecular_calls':0,
        'new_geometries':0, 'optimization_steps':0, 'affinity':None, 'classification':None}
    write(out/'START.json', result)
    try:
        run(a.config.resolve(), out, result)
    except Exception as exc:
        result.update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
    result['wall_seconds'] = time.monotonic()-start
    result['allocated_cpus'] = int(os.environ.get('SLURM_CPUS_PER_TASK','1'))
    result['application_allocated_CPU_seconds'] = result['wall_seconds']*result['allocated_cpus']
    write(out/'RESULT.json', result)
    print(json.dumps({'status':result['status'],'result':str(out/'RESULT.json')}),flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    sys.exit(main())
