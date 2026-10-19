#!/usr/bin/env python3
"""Saved-only conditional next-direction audit; no search or molecular engine."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback

PROTOCOL = 'lbt_7ccn_La_saved_conditional_direction_v1'
FIXED = {
    'collector': 'fff95c2f235947acb24e91991915b320731992f48c2b3d6c9ca3912a6e927f38',
    'force_report': 'e6b1cb9311b7e0cfe2fce2010a5839077c3814b5a1202ca4d03e2cbecd56c619',
    'checkpoint': '4e3a612380f93e9ba7f374636ed883a659c8b0ade6b35e2f7099f5d7643e9ecf',
    'point': 'cf59ef7add95948eb46e3b61f86014e473e0f04035d7a06c15885163e4042022',
    'optimizer': 'd93821da485f8d078e42bf7031018931632858ea720af617f26da2dab77617ee',
    'direction_helper': '07f29ca60adefa185ad89ec68a51fe3221149230c80d016c00d6994ceec98226',
    'search_implementation': 'f2419e88d3c43d12f9b476bb2b5175308657f723c70a3e52eb08398f064e8f4b',
    'reference_implementation': 'c0cca6dd4f44619fc2e3201f024b7913f6b7e9d8e8c4fca694fe17bd942d0ecf',
}
EV, HA, BOHR = 23.06054783061903, 627.509474, 0.529177210903


def require(ok, message):
    if not ok: raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as handle: sha = hashlib.file_digest(handle, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: '+p['path'])
    return Path(p['path'])


def read(p): return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')
    return pin(path)


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def saved_inputs(cfg, result):
    require(cfg['protocol'] == PROTOCOL and pin(__file__) == cfg['implementation'], 'Audit source differs')
    for key, sha in FIXED.items(): require(cfg[key]['sha256'] == sha, 'Frozen selection differs: '+key)
    for value in cfg.values():
        if isinstance(value, dict) and set(value) == {'path', 'sha256'}: verify(value)
    require(pin(sys.executable)['sha256'] == cfg['python']['sha256'] and
            str(Path(sys.prefix).resolve()) == cfg['python_prefix'], 'Saved-replay Python environment differs')
    import numpy as np
    import scipy
    result['runtime'] = dict(python=sys.version, executable=pin(sys.executable), prefix=sys.prefix,
        numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform(),
        threads={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
    sys.path.insert(0, str(verify(cfg['direction_helper']).parent))
    core = module('saved_direction_optimizer', cfg['optimizer'])
    helper = module('saved_direction_lbfgs', cfg['direction_helper'])
    require(pin(sys.modules['affordable_common'].__file__) == cfg['common_helper'], 'Helper dependency differs')
    collector, report = read(cfg['collector']), read(cfg['force_report'])
    search, run, cp, point = (read(cfg[k]) for k in ('search','search_result','checkpoint','point'))
    ccfg, rcfg = read(collector['config']), read(run['config'])
    require(collector['status'] == 'complete' and collector['stage'] == 'collect' and collector['full_four_available']
            and collector['slurm_job_id'] == '15173' and collector['implementation'] == cfg['reference_implementation']
            and ccfg['implementation'] == collector['implementation'] and ccfg['plan'] == cfg['plan']
            and collector['force_reports']['after_La'] == cfg['force_report'], 'Complete exact reference collector required')
    require(run['status'] == 'complete' and run['stage'] == 'run' and run['slurm_job_id'] == '12968'
            and run['implementation'] == cfg['search_implementation'] and rcfg['implementation'] == run['implementation']
            and run['optimizer'] == rcfg['optimizer'] == cfg['optimizer']
            and run['geometry_helper'] == cfg['direction_helper'] and run['rows'][0]['search'] == cfg['search'],
            'Exact completed Cartesian search lineage differs')
    require(search['status'] == 'complete' and search['source_id'] == '7CCN' and search['metal'] == 'La'
            and search['termination'] == 'accepted_step_limit' and not search['converged']
            and search['physical_state_admitted'] and search['final_geometry']['pass']
            and search['final_geometry']['state_compatible'] and search['checkpoint'] == cfg['checkpoint']
            and search['terminal_point'] == cp['point'] == cfg['point'] and search['start_point'] == cp['initial_point']
            and search['accepted_steps'] == cp['accepted_steps'] == 80
            and search['distinct_evaluations'] == cp['distinct_evaluations'] == 80
            and search['settings'] == cp['settings'] == core.SETTINGS
            and search['limits'] == cp['limits'] == {'accepted_steps':80,'evaluations':200}
            and cp['failure'] is None, 'Exact finite checkpoint differs')
    row = collector['rows']['after_La']; task = row['task']; sol = collector['solvent_admission']['after_La']
    require(all(collector['rows'][k]['status'] == 'available' for k in ('before_La','after_La','before_Dy','after_Dy')),
            'Reference denominator unavailable')
    endpoint = read(ccfg['endpoint_results']['after_La'])
    require(endpoint['status'] == 'complete' and endpoint['rows']['after_La'] == row, 'Admitted reference row differs')
    for k in ('receipt','output','engrad'): verify(row[k])
    for k in ('xyz','input'): verify(task[k])
    receipt = read(row['receipt'])
    require(receipt['returncode'] == 0 and receipt['normal_termination'] and receipt['scf_converged']
            and receipt['expected_engrad_complete'] and receipt['task_id'] == 'after_La'
            and receipt['artifacts']['engrad'] == row['engrad'] and receipt['artifacts']['output'] == row['output'],
            'Native analytic receipt differs')
    require(point['status'] == 'complete' and point['source_id'] == '7CCN' and point['metal'] == 'La'
            and point['compute_force'] and point['physical_source'] == report['physical'] == task['physical_parent']
            and point['positions'] == report['positions'] == task['positions']
            and point['ALPB'] == sol['receipt'] == task['ALPB'] and sol['status'] == 'available',
            'Same-geometry/common-solvent join differs')
    atoms = read(point['physical_source']); ids = [a['id'] for a in atoms]
    x = np.load(verify(point['positions']), allow_pickle=False)
    g = np.load(verify(point['gradient_kcal_mol_A']), allow_pickle=False)
    require(x.shape == g.shape == (256,3) and np.isfinite(x).all() and np.isfinite(g).all()
            and np.array_equal(x,[a['xyz_A'] for a in atoms]) and len(set(ids)) == 256
            and ids == point['recipe']['physical_ids'] == task['physical_ids'] and ids[-1] == 'metal'
            and [a['physical_id'] for a in report['force_table']] == ids, 'Full physical ordering differs')
    state = row['electronic_state']; native_task = read(point['task']); verify(point['xyz'])
    require(point['recipe']['physical_charge'] == native_task['charge'] == state['charge'] == -2
            and point['recipe']['physical_multiplicity'] == native_task['multiplicity'] == state['physical_multiplicity'] == 1
            and state['effective_multiplicity'] == 1 and state['ecp_core_electrons'] == 46
            and state['explicit_electrons'] == 1030 and task['target_index'] == point['recipe']['target_index'] == 255
            and native_task['physical_source'] == point['physical_source'] and native_task['xyz'] == point['xyz'],
            'Physical/native state differs')
    ml, alpb = read(point['MACE']), read(point['ALPB'])
    for value in alpb['artifacts'].values(): verify(value)
    raw_s = np.asarray([[float(v) for v in line.split()] for line in verify(alpb['artifacts']['gradient.dat']).read_text().splitlines()])
    require(raw_s.shape == (256,3) and np.isfinite(raw_s).all()
            and np.array_equal(raw_s,sol['gradient_hartree_per_bohr'])
            and np.array_equal(raw_s*HA/BOHR,sol['gradient_kcal_mol_A'])
            and np.array_equal(-np.asarray(ml['forces_eV_A'])*EV+raw_s*(HA/BOHR),g), 'Original component arithmetic differs')
    gr = np.asarray(row['gradient_kcal_mol_per_A'])+np.asarray(sol['gradient_kcal_mol_A'])
    require(gr.shape == (256,3) and np.isfinite(gr).all() and report['quantity'] == 'force=-gradient'
            and np.array_equal(-g,[a['force_kcal_mol_A']['whole_MACE_plus_fixed_ALPB'] for a in report['force_table']])
            and np.array_equal(-gr,[a['force_kcal_mol_A']['whole_PBE0_plus_fixed_ALPB'] for a in report['force_table']]),
            'Saved full-force representations differ')
    require(len(cp['history']) == 8 and core.SETTINGS['minimum_trust_A'] <= cp['trust_A'] <= core.SETTINGS['maximum_trust_A'],
            'Checkpoint history/trust differs')
    history = []
    for h in cp['history']:
        s,y = np.asarray(h['s']),np.asarray(h['y']); rho = h['rho']
        require(s.shape == y.shape == (768,) and np.isfinite(s).all() and np.isfinite(y).all()
                and np.isfinite(rho) and rho > 0, 'Invalid original history pair')
        history.append((s,y,rho))
    result.update(admitted_cells=1, state=state, SCF_evidence=row['SCF_evidence'], caveats=collector['caveats'],
        archived_force_metrics=search['final_full_force'], archived_force_errors=report['MACE_minus_PBE0_force_errors'],
        model_recipe=point['recipe'], saved_component_receipts={'MACE':point['MACE'],'ALPB':point['ALPB'],'PBE0':row['receipt']},
        other_cases={'before_La':'not_selected','before_Dy':'not_selected; model force unavailable',
                     'after_Dy':'not_selected; model force unavailable'})
    return np,core,helper,atoms,x,g,gr,cp,history


def analyze(cfg, result):
    np,core,helper,atoms,x,g,gr,cp,history = saved_inputs(cfg,result)
    pg,projection = core.remove_rigid(x,g)
    raw_d = helper.lbfgs_direction(pg.ravel(),history).reshape(x.shape)
    d,dprojection = core.remove_rigid(x,raw_d)
    before_fallback = d.copy(); test = float(np.sum(pg*d)); fallback = test >= -1e-14
    if fallback:
        history = []; d = -pg
    maximum = float(np.linalg.norm(d,axis=1).max())
    result['direction'] = dict(meaning='Conditional direction if a separately declared continuation were permitted; original80-step stop remains.',
        actually_proposed_or_executed=False, gradient_projection=projection, direction_projection=dprojection,
        projected_gradient=pg.tolist(), raw_LBFGS_direction=raw_d.tolist(), projected_LBFGS_direction=before_fallback.tolist(),
        fallback_test_sum_pg_d=test, fallback_threshold=-1e-14, fallback_taken=fallback,
        fallback_history_length=len(history), selected_before_trust=d.tolist(), maximum_atom_before_trust=maximum,
        original_zero_direction_threshold=1e-14, trust_A=cp['trust_A'])
    require(maximum >= 1e-14, 'Original formula would stop at zero_rigid_free_direction')
    scale = min(1.,cp['trust_A']/maximum); d *= scale
    slope = float(np.sum(g*d)); require(slope < 0,'Original non-descent direction gate failed')
    length = float(np.linalg.norm(d)); u = d/length
    derivatives = {name:{'dE_dalpha_kcal_mol':float(np.sum(grad*d)),
                         'dE_ds_kcal_mol_per_A':float(np.sum(grad*u))}
                   for name,grad in (('MACE_plus_fixed_ALPB',g),('PBE0_plus_fixed_ALPB',gr))}
    result['direction'].update(trust_scale=scale, proposed_displacement_A=d.tolist(),
        full_L2_length_A=length, normalized_full_direction=u.tolist(),
        original_model_slope_kcal_mol=slope, derivatives=derivatives,
        parameterization='X(alpha)=X+alpha*d, alpha dimensionless; s=alpha*||d||_2 in Angstrom. No coordinates evaluated.',
        nominal_local_descent={k:v['dE_dalpha_kcal_mol']<0 for k,v in derivatives.items()},
        signs_are_descriptive_not_accuracy_or_continuation_gates=True)
    error = gr-g  # F_model-F_reference = -g-(-gr).
    norms = np.linalg.norm(error,axis=1)
    def metrics(v):
        n=np.linalg.norm(v,axis=1)
        return {'max_atom_kcal_mol_A':float(n.max()),'RMS_atom_vector_kcal_mol_A':float(np.sqrt(np.mean(n*n))),
                'L2_kcal_mol_A':float(np.linalg.norm(v))}
    atom_rows=[]; groups={}
    for i,a in enumerate(atoms):
        rid=a['id'].rsplit('/',1)[0] if '/' in a['id'] else 'target'; groups.setdefault(rid,[]).append(i)
        atom_rows.append(dict(index=i,physical_id=a['id'],residue_id=rid,resname=a.get('resname'),
            model_force_kcal_mol_A=(-g[i]).tolist(),reference_force_kcal_mol_A=(-gr[i]).tolist(),
            model_minus_reference_force_kcal_mol_A=error[i].tolist(),force_error_norm_kcal_mol_A=float(norms[i]),
            model_dE_ds_contribution=float(np.sum(g[i]*u[i])),reference_dE_ds_contribution=float(np.sum(gr[i]*u[i]))))
    residues=[dict(residue_id=rid,indices=ix,atom_count=len(ix),force_error=metrics(error[ix]),
                  model_force=metrics(g[ix]),reference_force=metrics(gr[ix]),
                  model_dE_ds_contribution=float(np.sum(g[ix]*u[ix])),
                  reference_dE_ds_contribution=float(np.sum(gr[ix]*u[ix]))) for rid,ix in groups.items()]
    result.update(force_metrics={'model':metrics(g),'reference':metrics(gr),'model_minus_reference_error':metrics(error)},
        atom_rows=atom_rows,residue_rows=residues,
        atom_error_ranking=[a['index'] for a in sorted(atom_rows,key=lambda a:(-a['force_error_norm_kcal_mol_A'],a['index']))],
        residue_error_ranking=[a['residue_id'] for a in sorted(residues,key=lambda a:(-a['force_error']['L2_kcal_mol_A'],a['residue_id']))],
        ranking_is_display_only=True, residue_summaries_are_force_metrics_not_energy_decomposition=True, status='complete')


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--config',required=True); parser.add_argument('--output',required=True)
    args=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required before scientific imports')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK','0')) >= 1,'Allocated CPU required')
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=False); tick=time.monotonic()
    result=dict(protocol=PROTOCOL,status='unavailable',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],declared_cells=1,admitted_cells=0,declared_primary_directions=1,
        molecular_calls=0,model_loads=0,geometry_trials=0,optimizer_steps=0,continuation_authorized=False,
        physical_state_qualified=False,force_accuracy_qualified=False)
    try:
        cfg=read(result['config']); result['inputs']=cfg; analyze(cfg,result)
        verify(result['config'])
    except Exception as exc: result.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-tick
    write(out/'RESULT.json',result)
    return 0 if result['status']=='complete' else 1


if __name__=='__main__': raise SystemExit(main())
