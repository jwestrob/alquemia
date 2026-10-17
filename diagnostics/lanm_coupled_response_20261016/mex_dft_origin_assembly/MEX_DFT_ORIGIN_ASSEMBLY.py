#!/usr/bin/env python3
"""Saved-only two-state embedded-PBE0 + MM + fixed-FF ALPB assembly; Slurm only."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import sys
import time
import traceback

PROTOCOL = 'mex_saved_embedded_PBE0_MM_fixedFF_ALPB_origins_v1'
METALS = ('La', 'Dy')
GROUPS = ('retained_bonded', 'MM_Coulomb', 'MM_LJ', 'QM_MM_LJ')
COMPONENTS = ('PBE0_embedded', 'MM_ledger', 'fixed_FF_charge_ALPB')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Missing/changed pin: '+p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, data):
    with Path(path).open('x') as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write('\n')


def inputs(config):
    cfg = json.loads(config.read_text())
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation differs')
    verify(cfg['implementation']); verify(cfg['plan'])
    previous = read(cfg['reference_response_result'])
    require(previous['status'] == 'complete' and previous['slurm_job_id'] == '12686',
            'Admitted saved reference response required')
    oldcfg = read(previous['config'])
    require(previous['implementation']['sha256'] == oldcfg['implementation']['sha256']
            and previous['static_result'] == oldcfg['static_result']
            and previous['reference_collection'] == oldcfg['reference_collection'], 'Parent lineage differs')
    verify(oldcfg['implementation'])
    spec = importlib.util.spec_from_file_location('saved_mex_reference', verify(previous['implementation']))
    helper = importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
    world = helper.inputs(verify(previous['config']))
    require(helper.METALS == METALS, 'Parent metal denominator differs')
    require(set(cfg['solvent_stdout']) == set(METALS), 'Solvent evidence denominator differs')
    verify(world[1]['build']['binary'])
    for p in cfg['solvent_stdout'].values():
        verify(p)
    return cfg, previous, helper, world


def assemble(cfg, helper, world, metal):
    from affordable_common import energy
    from affordable_response import read_engrad
    from mace_omol_vacuum import terminal_scf_residuals
    np = world[-1]
    mp, arrays, evidence = helper.gradients(world, metal)
    old = read(world[10][metal]['assembled'])
    rr = world[6]['rows'][metal+'_A']
    output = verify(rr['output']).read_text()
    dft = energy(verify(rr['output']))
    raw = read_engrad(verify(rr['engrad']))
    require(dft == rr['energy_hartree'] and abs(raw['energy_Ha']-dft) <= 1e-8,
            'Native DFT energy replay differs')
    scf = terminal_scf_residuals(output)
    for v in scf.values():
        require(math.isfinite(v['value']) and math.isfinite(v['printed_tolerance']), 'Nonfinite SCF evidence')
        v['within_printed_tolerance'] = abs(v['value']) <= v['printed_tolerance']
    scf_lines = [s.strip() for s in output.splitlines() if 'ConvCheckMode' in s or 'ConvForced' in s]
    require(any('ConvCheckMode' in s for s in scf_lines) and any('ConvForced' in s for s in scf_lines),
            'Missing native convergence semantics')

    classical = read(old['classical_reused']); ledger = read(classical['ledger'])
    for p in ledger['artifacts'].values():
        verify(p)
    verify(ledger['native_system'])
    require(classical['QM_MM_Coulomb_included'] is False and classical['C4_induction_included'] is False
            and set(classical['components_kcal_mol']) == set(GROUPS)
            and set(classical['component_gradients_kcal_mol_A']) == set(GROUPS), 'MM exclusions/groups differ')
    mm = classical['total_kcal_mol']
    mm_residual = abs(math.fsum(classical['components_kcal_mol'].values())-mm)
    require(math.isfinite(mm) and mm_residual <= 1e-8, 'MM component energy closure differs')
    ids = [a['id'] for a in mp['physical_atoms']]
    ix = {sid:i for i,sid in enumerate(classical['source_ids'])}
    order = [ix[sid] for sid in ids]
    require(len(ix) == len(ids) == 2049 and set(ix) == set(ids), 'MM physical IDs differ')
    require(np.allclose(np.asarray(classical['coordinates_A'])[order],
                       [a['xyz_A'] for a in mp['physical_atoms']], atol=1e-8, rtol=0), 'MM coordinates differ')
    mm_gradient_residual = helper.close(np,
        sum(np.asarray(classical['component_gradients_kcal_mol_A'][k])[order] for k in GROUPS),
        arrays['MM_ledger'], 'MM component gradient closure')

    solvent = read(old['ALPB']); execution = read(solvent['execution'])
    require(execution['returncode'] == 0, 'Archived solvent execution failed')
    stdout_path = verify(cfg['solvent_stdout'][metal])
    require(stdout_path.parent == verify(solvent['execution']).parent, 'Solvent stdout parent differs')
    text = stdout_path.read_text()
    def values(label, n):
        hits = re.findall(r'^'+label+r'\s+(.+)$', text, re.M)
        require(len(hits) == 1, 'Missing/duplicate solvent '+label)
        value = [float(v) for v in hits[0].split()]
        require(len(value) == n and all(math.isfinite(v) for v in value), 'Invalid solvent '+label)
        return value
    parts = dict(zip(('Gelec','Ghb','Gsasa','Gshift'), values('COMPONENTS',4)))
    require(parts == solvent['components_hartree'] and values('CHARGE_SUM',1)[0] == solvent['charge_e']
            and abs(solvent['charge_e']-2) <= 1e-5 and values('ROUNDING_BOUNDS',2) == [0.,0.]
            and values('CAVITY_GLOBAL',4) == solvent['cavity_global'], 'Solvent raw replay differs')
    alpb = solvent['total_hartree']; alpb_residual = abs(math.fsum(parts.values())-alpb)
    require(math.isfinite(alpb) and alpb_residual == 0., 'ALPB component energy closure differs')
    g = np.loadtxt(verify(solvent['gradient']))
    require(g.shape == (2049,3) and np.isfinite(g).all()
            and np.array_equal(g, solvent['gradient_hartree_per_bohr']), 'Raw ALPB gradient replay differs')
    q = np.loadtxt(verify(solvent['charges']))
    require(q.shape == (2049,) and np.isfinite(q).all() and abs(math.fsum(q)-2) <= 1e-5,
            'Physical FF charge vector differs')
    xyz = verify(solvent['xyz']).read_text().splitlines(); coords = [s.split() for s in xyz[2:]]
    require(int(xyz[0]) == len(coords) == 2049
            and [s[0] for s in coords] == [a['element'] for a in mp['physical_atoms']]
            and np.array_equal(np.asarray([s[1:] for s in coords], float),
                               [a['xyz_A'] for a in mp['physical_atoms']]), 'ALPB physical coordinates differ')
    require(execution['command'] == [world[1]['build']['binary']['path'], solvent['xyz']['path'],
                                     solvent['charges']['path']], 'Solvent execution input differs')
    energies = {'PBE0_embedded': dft*helper.HA, 'MM_ledger':mm, 'fixed_FF_charge_ALPB':alpb*helper.HA}
    require(energies['MM_ledger'] == old['component_energies_kcal_mol']['MM_ledger']
            and energies['fixed_FF_charge_ALPB'] == old['component_energies_kcal_mol']['fixed_FF_charge_ALPB'],
            'Unchanged shared component energy differs')
    gradients = {k:np.asarray(arrays[k], float) for k in COMPONENTS}
    for k,g in gradients.items():
        require(g.shape == (2049,3) and np.isfinite(g).all(), 'Invalid gradient '+k)
    total = sum(gradients.values())
    total_residual = helper.close(np, total, arrays['PBE0_plus_shared_MM_ALPB'], 'Reference total replay')
    return {'metal':metal, 'status':'complete', 'physical_atoms':mp['physical_atoms'], 'source_ids':ids,
        'source_inputs':world[4], 'target_source_id':classical['target_source_id'],
        'spectator_occupancy':world[5]['spectator_occupancy'], 'physical_holo_charge_e':2,
        'inventory':{'physical_atoms':2049, 'electronic_atoms_including_caps':209,
                     'real_QM_atoms':205, 'caps':4, 'MM_atoms':1844, 'embedding_sites':1840},
        'physical_mapping':old['physical_mapping'], 'electronic_state':rr['electronic_state'],
        'native_energies':{'PBE0_embedded_hartree':dft, 'MM_kcal_mol':mm, 'ALPB_hartree':alpb},
        'component_energies_kcal_mol':energies, 'total_energy_kcal_mol':math.fsum(energies.values()),
        'MM_components_kcal_mol':classical['components_kcal_mol'], 'ALPB_components_hartree':parts,
        'component_gradients_kcal_mol_A':{k:g.tolist() for k,g in gradients.items()},
        'gradient_kcal_mol_A':total.tolist(), 'force_kcal_mol_A':(-total).tolist(),
        'translation_gradient_sum_kcal_mol_A':total.sum(axis=0).tolist(),
        'component_translation_gradient_sums_kcal_mol_A':{k:g.sum(axis=0).tolist() for k,g in gradients.items()},
        'quantity':'gradient=dE/dR; force=-gradient', 'reference_evidence':evidence,
        'classical_reused':old['classical_reused'], 'ALPB_reused':old['ALPB'],
        'solvent_stdout':cfg['solvent_stdout'][metal], 'solvent_charges':solvent['charges'],
        'SCF':{'native_converged':True, 'printed_residuals':scf, 'convergence_semantics_lines':scf_lines,
               'all_printed_residuals_within_tolerance':all(v['within_printed_tolerance'] for v in scf.values()),
               'warning_lines':[s.strip() for s in output.splitlines() if re.search(r'warning|small.*gap',s,re.I)],
               'all_residuals_flag_is_admission_gate':False},
        'replay_residuals':{'DFT_engrad_energy_hartree':abs(raw['energy_Ha']-dft),
            'MM_energy_sum_kcal_mol':mm_residual, 'ALPB_energy_sum_hartree':alpb_residual,
            'MM_gradient_sum_kcal_mol_A':mm_gradient_residual, 'total_gradient_kcal_mol_A':total_residual},
        'QM_MM_classical_Coulomb_included':False, 'C4_induction_included':False,
        'solvent_charge_response_included':False, 'total_derivative_FD_qualified':False,
        'physical_state_qualified':False, 'affinity':None}


def run(config, out, result):
    cfg, previous, helper, world = inputs(config); np = world[-1]
    result.update(reference_response_result=cfg['reference_response_result'],
        static_result=previous['static_result'], reference_collection=previous['reference_collection'],
        reference_manifest=previous['reference_manifest'], reference_method=previous['reference_method'],
        source_inputs=world[4], reference_basis_assets=previous['reference_basis_assets'],
        solvent_build=world[1]['build'], dependencies=world[0]['dependencies'],
        runtime={'python':pin(sys.executable), 'python_version':sys.version, 'numpy_version':np.__version__})
    for i,metal in enumerate(METALS):
        try:
            row = assemble(cfg, helper, world, metal)
            write(out/(metal+'_ASSEMBLED.json'), row)
            row['result'] = pin(out/(metal+'_ASSEMBLED.json')); result['rows'][i] = row
        except Exception as exc:
            result['rows'][i].update(failure=repr(exc), traceback=traceback.format_exc())
    complete = all(r['status'] == 'complete' for r in result['rows'])
    if complete:
        a,b = result['rows']
        require(a['source_ids'] == b['source_ids'] and np.array_equal(
            [v['xyz_A'] for v in a['physical_atoms']], [v['xyz_A'] for v in b['physical_atoms']])
            and a['solvent_charges']['sha256'] == b['solvent_charges']['sha256'], 'Paired source/charge identity differs')
    result.update(status='complete' if complete else 'unavailable', accepted=sum(r['status']=='complete' for r in result['rows']))
    lines = ['# Saved Mex DFT origin assembly', '',
        'Ledger: embedded PBE0 + existing MM + fixed-FF-charge global ALPB.', '',
        '| Metal | Status | All printed SCF residuals within tolerance |', '|---|---|---|']
    for row in result['rows']:
        lines.append('| '+row['metal']+' | '+row['status']+' | '+str(row.get('SCF',{}).get('all_printed_residuals_within_tolerance'))+' |')
    lines += ['', 'Full component energies, physical gradients, forces and raw receipt pins are in each',
        'ASSEMBLED.json. The two absolute metal energies are not an affinity comparison.',
        'Native SCF convergence is retained; printed residual failures and warnings remain visible.',
        'Physical electronic states and the replacement total derivative remain unqualified.',
        'Zero molecular calls, new geometries and optimization steps; no default change.', '']
    (out/'REPORT.md').write_text('\n'.join(lines)); result['report'] = pin(out/'REPORT.md')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); require(os.environ.get('SLURM_JOB_ID','').isdigit(), 'Saved scientific assembly requires Slurm')
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    cfg = json.loads(a.config.read_text())
    result = {'protocol_id':PROTOCOL, 'status':'running', 'implementation':pin(__file__),
        'config':pin(a.config), 'plan':cfg['plan'], 'slurm_job_id':os.environ['SLURM_JOB_ID'],
        'declared_states':list(METALS), 'rows':[{'metal':m,'status':'unavailable'} for m in METALS],
        'molecular_calls':0, 'new_geometries':0, 'optimization_steps':0, 'affinity':None, 'classification':None}
    write(out/'START.json', result)
    try:
        run(a.config.resolve(), out, result)
    except Exception as exc:
        result.update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
    result['wall_seconds'] = time.monotonic()-start
    result['allocated_cpus'] = int(os.environ.get('SLURM_CPUS_PER_TASK','1'))
    result['application_allocated_CPU_seconds'] = result['wall_seconds']*result['allocated_cpus']
    write(out/'RESULT.json', result)
    print(json.dumps({'status':result['status'],'result':str(out/'RESULT.json')}), flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    sys.exit(main())
