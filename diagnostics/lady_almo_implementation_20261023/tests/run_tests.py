#!/usr/bin/env python3
"""Scheduled independent ALMO checks; never launches a subprocess or retries.

Example: python run_tests.py --almo-source ../almo.py --mode core --output NEW.json
Modes: core=PBE0 correctness; benchmark=PBE independent-implementation values.
All molecular and linear-algebra work occurs only after the Slurm environment
check. Source preparation/importing this module does not import numerical tools.
"""
import argparse
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import traceback


HERE = Path(__file__).resolve().parent


def pin(path):
    p = Path(path).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def require(value, message):
    if not value:
        raise AssertionError(message)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--almo-source', required=True, type=Path)
    ap.add_argument('--mode', choices=['core', 'benchmark'], required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    require(os.environ.get('SLURM_JOB_ID'), 'Numerical tests require a real Slurm job')
    require(os.environ.get('SLURM_JOB_QOS') == 'standard_expand', 'Require standard_expand')
    require(not args.output.exists(), 'Refuse to overwrite an existing result')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Serial numerical work within this job. Independent jobs may run in parallel.
    for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS']:
        os.environ[key] = '1'
    result = {'status': 'failed', 'mode': args.mode, 'checks': [],
              'implementation': pin(__file__), 'almo_source': pin(args.almo_source),
              'fixture': pin(HERE/'FIXTURE.json'),
              'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'runtime': {'host': os.uname().nodename, 'pid': os.getpid(),
                          'job_id': os.environ['SLURM_JOB_ID'],
                          'qos': os.environ['SLURM_JOB_QOS'], 'numerical_threads': 1},
              'error': None}
    start, cpu = time.monotonic(), time.process_time()
    stream = None
    try:
        import numpy as np
        import pyscf
        from pyscf import dft, gto, lib
        lib.num_threads(1)
        result['versions'] = {'numpy': np.__version__, 'pyscf': pyscf.__version__}
        spec = importlib.util.spec_from_file_location('tested_almo', args.almo_source)
        almo = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = almo
        spec.loader.exec_module(almo)
        fixture = json.loads((HERE/'FIXTURE.json').read_text())
        for key in ['cp2k_basis', 'cp2k_potential']:
            require(pin(fixture[key]['path']) == fixture[key], 'Literal fixture pin changed')
        atoms = [(r[0], tuple(map(float, r[1:]))) for r in fixture['geometry_A_tokens']]
        basis = {s: [[sh[0], *[list(map(float, row)) for row in sh[1:]]]
                     for sh in shells] for s, shells in fixture['basis_pyscf_tokens'].items()}
        log = args.output.with_suffix('.pyscf.log')
        require(not log.exists(), 'Refuse to overwrite a backend log')
        stream = log.open('x')
        xc = 'PBE0' if args.mode == 'core' else 'PBE'
        atom_grid = (99, 590) if args.mode == 'core' else (150, 974)
        result['method'] = {'xc': xc, 'basis': 'literal CP2K 6-31G*', 'cart': False,
                            'atom_grid': list(atom_grid), 'prune': None,
                            'density_fitting': False, 'ecp': None, 'dispersion': None,
                            'fragment_charges': [0, 0], 'fragment_spins': [0, 0]}

        def check(name, passed, **evidence):
            result['checks'].append({'name': name, 'passed': bool(passed), **evidence})
            require(passed, name)

        def make_mf(atom_rows):
            mol = gto.M(atom=atom_rows, basis=basis, unit='Angstrom', charge=0,
                        spin=0, cart=False, verbose=4, max_memory=3000)
            mol.stdout = stream
            mf = dft.RKS(mol)
            mf.stdout = stream
            mf.xc = xc
            mf.grids.atom_grid = atom_grid
            mf.grids.prune = None
            mf.small_rho_cutoff = 0.0
            mf.conv_tol = 1e-12
            mf.conv_tol_grad = 1e-9
            mf.max_cycle = 150
            mf.chkfile = None
            mf.grids.build()
            return mf

        mf = make_mf(atoms)
        S = mf.get_ovlp()
        check('AO_and_electron_count', mf.mol.nao_nr() == 36 and mf.mol.nelectron == 20,
              nao=mf.mol.nao_nr(), electrons=mf.mol.nelectron)
        slices = mf.mol.aoslice_by_atom()
        aos = [[i for atom in frag for i in range(slices[atom, 2], slices[atom, 3])]
               for frag in fixture['fragment_atoms_zero_based']]
        nocc = [5, 5]
        mask = np.zeros((36, 10), dtype=bool)
        C0 = np.zeros((36, 10))
        fragment_energies = []
        for f, atom_ids in enumerate(fixture['fragment_atoms_zero_based']):
            fm = make_mf([atoms[i] for i in atom_ids])
            fm.kernel()
            check('isolated_fragment_%d_converged' % (f+1), fm.converged)
            occ = fm.mo_coeff[:, fm.mo_occ > 0]
            check('fragment_%d_integer_occupation' % (f+1),
                  occ.shape == (18, 5) and np.all(fm.mo_occ[fm.mo_occ > 0] == 2))
            cols = list(range(5*f, 5*f+5))
            mask[np.ix_(aos[f], cols)] = True
            C0[np.ix_(aos[f], cols)] = occ
            fragment_energies.append(float(fm.e_tot))
        L = np.linalg.cholesky(S)

        def independent_density(C):
            # QR of AO-metric-whitened occupied vectors constructs the same
            # occupied subspace without reusing the implementation's inverse.
            q, _ = np.linalg.qr(L.T @ C, mode='reduced')
            orth = np.linalg.solve(L.T, q)
            return 2.0 * orth @ orth.T

        def backend_energy(C):
            return float(mf.energy_tot(dm=independent_density(C)))

        def physical_checks(tag, C, density):
            dref = independent_density(C)
            electron_error = abs(float(np.trace(density @ S))-20)
            idem = float(np.max(np.abs(density @ S @ density - 2*density)))
            support = float(np.max(np.abs(C[~mask])))
            derror = float(np.max(np.abs(density-dref)))
            check(tag+'_electron_idempotency_support', electron_error <= 2e-9 and
                  idem <= 2e-9 and support <= 1e-13 and derror <= 2e-9,
                  electron_error=electron_error, idempotency_max=idem,
                  forbidden_support_max=support, independent_density_max=derror)

        problem = almo.ALMOProblem(mf, aos, nocc)
        frozen = problem.evaluate(C0)
        physical_checks('frozen', C0, frozen.density)
        Efrozen = backend_energy(C0)
        check('frozen_backend_energy', abs(frozen.energy-Efrozen) <= 2e-10,
              difference_hartree=float(frozen.energy-Efrozen))

        if args.mode == 'core':
            rng = np.random.default_rng(20261023)
            Ctest = C0 + 0.025 * rng.standard_normal(C0.shape) * mask
            ev = problem.evaluate(Ctest)
            physical_checks('perturbed', Ctest, ev.density)
            G = np.asarray(ev.coefficient_gradient)
            fdrows = []
            for k in range(4):
                direction = rng.standard_normal(C0.shape) * mask
                if k < 2:
                    direction[:, 5*(1-k):5*(1-k)+5] = 0
                direction /= np.linalg.norm(direction)
                analytic = float(np.sum(G*direction))
                for h in (1e-4, 5e-5):
                    slope = (backend_energy(Ctest+h*direction)-backend_energy(Ctest-h*direction))/(2*h)
                    allowance = 2e-7 + 1e-6*abs(analytic)
                    fdrows.append({'direction': k, 'h': h, 'analytic': analytic,
                                   'finite_difference': slope, 'error': abs(slope-analytic),
                                   'allowance': allowance})
            result['orbital_gradient_finite_differences'] = fdrows
            check('orbital_gradient_two_step_FD', all(r['error'] <= r['allowance'] for r in fdrows))
            check('FD_control_has_signal', max(abs(r['finite_difference']) for r in fdrows) > 1e-4)
            check('FD_rejects_factor_two_corruption', any(
                abs(r['finite_difference']-2*r['analytic']) > r['allowance'] for r in fdrows))
            Cmix = Ctest.copy()
            for f in range(2):
                b = np.diag(np.linspace(0.7, 1.3, 5))
                b += np.diag(np.full(4, 0.17), 1)
                Cmix[:, 5*f:5*f+5] = Cmix[:, 5*f:5*f+5] @ b
            mixed = problem.evaluate(Cmix)
            check('nonsingular_within_block_occupied_gauge',
                  abs(mixed.energy-ev.energy) <= 2e-10 and
                  np.max(np.abs(mixed.density-ev.density)) <= 2e-9,
                  energy_difference=float(mixed.energy-ev.energy),
                  density_max=float(np.max(np.abs(mixed.density-ev.density))))
            for label, bad in [('rank_deficient', C0.copy()), ('outside_fragment', C0.copy())]:
                if label == 'rank_deficient':
                    bad[:, 0] = 0
                else:
                    bad[aos[1][0], 0] = 0.01
                rejected = False
                try:
                    problem.evaluate(bad)
                except (ValueError, np.linalg.LinAlgError):
                    rejected = True
                check(label+'_rejected', rejected)

        optimum = almo.solve_almo(mf, aos, nocc, initial_coefficients=C0,
                                 max_iterations=300, max_evaluations=1500,
                                 gradient_tolerance=1e-8, energy_tolerance=1e-12)
        result['optimizer'] = {'status': optimum.status, 'converged': bool(optimum.converged),
                               'gradient_norm': float(optimum.gradient_norm),
                               'history': optimum.history}
        check('ALMO_converged', optimum.converged)
        physical_checks('optimized', optimum.coefficients, optimum.density)
        Ealmo = backend_energy(optimum.coefficients)
        check('optimized_backend_energy', abs(Ealmo-optimum.energy) <= 2e-10)
        mf.kernel(dm0=optimum.density)
        check('ordinary_full_SCF_converged', mf.converged)
        Efull = float(mf.e_tot)
        check('variational_order', Efull <= Ealmo+2e-8 and Ealmo <= Efrozen+2e-8,
              Efull=Efull, Ealmo=Ealmo, Efrozen=Efrozen, allowance_hartree=2e-8)
        if args.mode == 'core':
            one = almo.solve_almo(mf, [list(range(36))], [10], initial_coefficients=C0,
                                 max_iterations=300, max_evaluations=1500,
                                 gradient_tolerance=1e-8, energy_tolerance=1e-12)
            check('one_fragment_reduces_to_ordinary_SCF', one.converged and
                  abs(one.energy-Efull) <= 2e-7,
                  difference_hartree=float(one.energy-Efull), tolerance_hartree=2e-7)
            released = almo.release_constraints(mf, optimum, conv_tol=1e-12,
                                                conv_tol_grad=1e-9, max_cycle=150)
            check('release_constraints_ordinary_SCF', released['converged'] and
                  abs(released['energy']-Efull) <= 2e-8,
                  difference_hartree=float(released['energy']-Efull))
        R = sum(fragment_energies)
        result['energies_hartree'] = {'fragment_1': fragment_energies[0],
            'fragment_2': fragment_energies[1], 'R': R, 'F': Efrozen, 'P': Ealmo, 'U': Efull}
        result['components_hartree'] = {'FRZ': Efrozen-R, 'POL': Ealmo-Efrozen,
                                       'CT': Efull-Ealmo, 'interaction': Efull-R}
        result['status'] = 'passed_core_tests' if args.mode == 'core' else 'passed_pyscf_benchmark'
        result['independent_native_comparison'] = 'pending separate CP2K outputs'
    except Exception:
        result['error'] = traceback.format_exc()
    finally:
        if stream is not None:
            stream.close()
        result['cost'] = {'wall_seconds': time.monotonic()-start,
                          'process_cpu_seconds': time.process_time()-cpu,
                          'scheduler_accounting': 'root collects actual parent allocation separately'}
        with args.output.open('x') as f:
            json.dump(result, f, indent=2, allow_nan=False, default=lambda x: x.item())
            f.write('\n')
    return 0 if result['status'].startswith('passed_') else 1


if __name__ == '__main__':
    sys.exit(main())
