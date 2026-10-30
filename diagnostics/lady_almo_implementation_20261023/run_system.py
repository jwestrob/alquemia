"""Slurm-only molecular runner for the exact archived La/Dy quartet.

Baseline mode is independent of the custom ALMO implementation. Every job has
fresh artifacts; failures retain native logs, history and an explicit result.
"""
from __future__ import annotations

import argparse
import ast
import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time
import traceback

CALL_COUNTS = {"scf_kernel": 0}

def write_json(path, value):
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def metadata_versions():
    from pyscf import gto, dft, scf
    import pyscf
    import archive
    modules = [gto, gto.mole, gto.ecp, dft, dft.rks, scf.hf]
    pins = [archive.pin(m.__file__) for m in modules]
    pins += [archive.pin(Path(pyscf.__file__).parent / "gto/basis/parse_nwchem_ecp.py")]
    try:
        from threadpoolctl import threadpool_info
        threadpools = threadpool_info()
    except ImportError:
        threadpools = {"status": "threadpoolctl_not_installed",
                       "observed_pyscf_openmp_threads": pyscf.lib.num_threads(),
                       "effective_blas_thread_count_measured": False}
    return {"python": sys.version, "executable": sys.executable,
            "platform": platform.platform(),
            "packages": {name: importlib.metadata.version(name)
                         for name in ("pyscf", "numpy", "scipy", "h5py")},
            "pyscf_module_pins": pins,
            "libxc_version": str(dft.libxc.__version__),
            "threadpools": threadpools,
            "thread_environment": {name: os.environ.get(name) for name in
                                   ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                    "NUMEXPR_NUM_THREADS", "OMP_MAX_ACTIVE_LEVELS")}}


def configure_mean_field(mol, checkpoint):
    from pyscf import dft
    mf = dft.RKS(mol)
    mf.xc = "PBE0"
    mf.grids.atom_grid = (150, 974)
    mf.grids.prune = None
    mf.small_rho_cutoff = 0
    mf.conv_tol = 1e-10
    mf.conv_tol_grad = 1e-7
    mf.max_cycle = 300
    mf.max_memory = mol.max_memory
    mf.direct_scf = True
    mf.direct_scf_tol = 1e-13
    mf.init_guess = "minao"
    mf.chkfile = str(checkpoint)
    return mf


def physics_signature(info):
    """Exact actual subset inputs, allowing identical fragments across contexts."""
    subset = info["runtime_subset"]
    atoms = [info["atoms"][i] for i in subset["atom_indices"]]
    elements = sorted({a[0] for a in atoms})
    payload = {"atoms_angstrom": atoms, "charge": subset["charge"],
               "spin": subset["spin"], "nelectron": subset["nelectron"],
               "basis": {e: info["basis"][e] for e in elements},
               "ecp": {e: info["ecp"][e] for e in elements if e in info["ecp"]},
               "cart": False}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def common_energy_settings(settings):
    """Admit only the two declared starting guesses, with every other key exact."""
    expected = numerical_settings()
    if set(settings) != set(expected) or settings["initial_guess"] not in ("minao", "1e"):
        raise ValueError("Unknown/incomplete numerical settings or initial guess")
    common = {key: value for key, value in settings.items() if key != "initial_guess"}
    if common != {key: value for key, value in expected.items() if key != "initial_guess"}:
        raise ValueError("Electronic model, quadrature or SCF criteria changed")
    return common


def validate_initialization(run, info):
    """A narrow Dy55 initial-guess exception; never a solver/model bypass."""
    import archive
    common_energy_settings(run["settings"])
    guess = run["settings"]["initial_guess"]
    if guess == "minao":
        if "initialization" in run and run["initialization"]["selected_guess"] != "minao":
            raise ValueError("Declared and recorded starting guesses differ")
        return {"selected_guess": "minao", "recovery": False}
    subset = info["runtime_subset"]
    selected = [info["atoms"][i][0] for i in subset["atom_indices"]]
    if "Dy" not in selected or info["ecp"].get("Dy", [None])[0] != 55:
        raise ValueError("Native 1e recovery is restricted to Dy55-containing states")
    if run["mode"] not in ("baseline", "fragment") or (run["mode"] == "fragment" and subset["fragment"] != 0):
        raise ValueError("Native 1e recovery is restricted to Dy baselines and the Dy metal fragment")
    record = run["initialization"]
    if record["selected_guess"] != "1e" or record["recovery_kind"] != "Dy55_minao_rectangular_projection_to_native_1e":
        raise ValueError("Unrecognized initialization recovery declaration")
    failed_path = archive.verify_pin(record["previous_failure"])
    previous = json.loads(failed_path.read_text())
    failed_info = json.loads(archive.verify_pin(record["previous_input"]).read_text())
    previous_start = json.loads(archive.verify_pin(record["previous_start"]).read_text())
    archive.verify_pin(record["pyscf_hf_source"])
    if previous["status"] != "failed" or previous["mode"] != run["mode"] or previous["system_id"] != run["system_id"]:
        raise ValueError("Initialization recovery does not bind the matching failed task")
    if previous["settings"] != numerical_settings("minao") or previous_start["settings"] != previous["settings"]:
        raise ValueError("Failed task used a different scientific protocol")
    if previous.get("error") != "LinAlgError('Last 2 dimensions of the array must be square')":
        raise ValueError("Recovery is not the observed minao rectangular-determinant failure")
    if "minao_basis" not in previous.get("traceback", "") or "s12[occ4ecp>0][:,occ>0]" not in previous["traceback"]:
        raise ValueError("The bound failure did not occur in the reviewed minao branch")
    if physics_signature(failed_info) != physics_signature(info):
        raise ValueError("Initialization recovery changes the failed state's physical inputs")
    if previous["runtime"]["packages"] != run["runtime"]["packages"]:
        raise ValueError("Initialization recovery changes the electronic runtime")
    old_hf = [p for p in previous["runtime"]["pyscf_module_pins"] if p["path"].endswith("/scf/hf.py")]
    if len(old_hf) != 1 or old_hf[0]["sha256"] != record["pyscf_hf_source"]["sha256"]:
        raise ValueError("Initialization recovery changed the reviewed PySCF source")
    return record


def prepare_initialization(args, runtime):
    """Bind explicit CLI recovery to the original failed receipt before SCF."""
    import archive
    if args.initial_guess == "minao":
        if args.initialization_recovery_from is not None:
            raise ValueError("A recovery receipt requires explicit --initial-guess 1e")
        return {"selected_guess": "minao", "recovery": False,
                "used_for_scf": args.mode in ("baseline", "fragment")}
    if args.initialization_recovery_from is None:
        raise ValueError("Native 1e requires --initialization-recovery-from FAILED_RESULT_JSON")
    if args.mode not in ("baseline", "fragment") or (args.mode == "fragment" and args.fragment != 0):
        raise ValueError("Native 1e is only a Dy baseline/metal-fragment initialization recovery")
    info = archive.load_system(args.system)
    subset = (info["fragments"][args.fragment] if args.mode == "fragment" else
              {"atom_indices": list(range(len(info["atoms"]))),
               "charge": info["charge"], "spin": info["spin"], "nelectron": info["nelectron"]})
    info["runtime_subset"] = {**subset, "fragment": args.fragment if args.mode == "fragment" else None}
    failed_path = args.initialization_recovery_from.resolve()
    from pyscf.scf import hf
    record = {"selected_guess": "1e", "recovery": True, "used_for_scf": True,
              "recovery_kind": "Dy55_minao_rectangular_projection_to_native_1e",
              "previous_failure": archive.pin(failed_path),
              "previous_input": archive.pin(failed_path.parent / "INPUT.json"),
              "previous_start": archive.pin(failed_path.parent / "START.json"),
              "pyscf_hf_source": archive.pin(hf.__file__),
              "native_path": "SCF.init_guess_by_1e: get_hcore, get_ovlp, generalized eigensystem, integer occupations",
              "scope": "Only the initial density changes; Hamiltonian, bases, ECPs, geometry, electron count, grid, SCF criteria and all admission gates stay fixed.",
              "state_limit": "A different initial density may find a different SCF stationary solution; convergence and archived-energy gates remain mandatory."}
    candidate = {"mode": args.mode, "system_id": args.system,
                 "settings": numerical_settings("1e"), "initialization": record, "runtime": runtime}
    validate_initialization(candidate, info)
    return record


def restore_state(directory, label, mf, info):
    """Reuse only an admitted, hashed state with identical inputs/settings."""
    import archive
    import numpy as np
    directory = Path(directory).resolve()
    run = json.loads((directory / "RESULT.json").read_text())
    report = json.loads((directory / f"{label}_RESULT.json").read_text())
    previous_info = json.loads((directory / "INPUT.json").read_text())
    start = json.loads((directory / "START.json").read_text())
    if run["status"] != "passed" or report["status"] != "passed" or not all(report["checks"].values()):
        raise ValueError(f"Cannot reuse an unadmitted {label} state: {directory}")
    if physics_signature(previous_info) != physics_signature(info):
        raise ValueError(f"Saved {label} input signature differs")
    if start["settings"] != run["settings"]:
        raise ValueError("Saved startup/final settings differ")
    initialization = validate_initialization(run, previous_info)
    if common_energy_settings(start["settings"]) != common_energy_settings(numerical_settings()):
        raise ValueError(f"Saved {label} electronic/grid/SCF settings differ")
    current_versions = metadata_versions()
    if start["runtime"]["packages"] != current_versions["packages"] or start["runtime"]["libxc_version"] != current_versions["libxc_version"]:
        raise ValueError("Saved-state backend package or LibXC version differs")
    for record in start["runtime"]["pyscf_module_pins"]:
        archive.verify_pin(record)
    runner_pin = next(p for p in start["source_pins"] if Path(p["path"]).name == "run_system.py")
    saved_source = archive.verify_pin(runner_pin).read_text()
    def configuration_ast(text):
        node = next(n for n in ast.parse(text).body
                    if isinstance(n, ast.FunctionDef) and n.name == "configure_mean_field")
        return ast.dump(node, include_attributes=False)
    if configuration_ast(saved_source) != configuration_ast(Path(__file__).read_text()):
        raise ValueError("Saved-state mean-field/grid constructor differs")
    state_path = directory / f"{label}_state.npz"
    if archive.sha256(state_path) != report["state_sha256"]:
        raise ValueError("Saved state hash mismatch")
    with np.load(state_path, allow_pickle=False) as state:
        mf.mo_coeff = state["mo_coeff"].copy()
        mf.mo_occ = state["mo_occ"].copy()
        mf.mo_energy = state["mo_energy"].copy()
        mf.e_tot = float(state["e_tot"])
        saved_overlap = state["overlap"].copy()
        saved_density = state["density"].copy()
    if mf.mo_coeff.shape != (mf.mol.nao_nr(), mf.mol.nao_nr()):
        raise ValueError("Saved coefficient dimensions differ")
    if not np.all(np.isfinite(mf.mo_coeff)) or not np.all((mf.mo_occ == 0) | (mf.mo_occ == 2)):
        raise ValueError("Invalid saved coefficients or occupations")
    if int(mf.mo_occ.sum()) != mf.mol.nelectron or mf.e_tot != report["energy_hartree"]:
        raise ValueError("Saved electron/energy identity mismatch")
    if not np.allclose(mf.get_ovlp(), saved_overlap, rtol=1e-11, atol=1e-12):
        raise ValueError("Saved AO overlap/order differs")
    if not np.allclose(mf.make_rdm1(), saved_density, rtol=1e-10, atol=1e-11):
        raise ValueError("Saved density does not reproduce from orbitals")
    mf.converged = True
    return {"result": archive.pin(directory / f"{label}_RESULT.json"),
            "run_result": archive.pin(directory / "RESULT.json"),
            "input": archive.pin(directory / "INPUT.json"),
            "state": archive.pin(state_path), "energy_hartree": mf.e_tot,
            "physics_signature": physics_signature(info), "new_scf_calls": 0,
            "grid_points": report["grid_points"],
            "initialization": initialization,
            "same_mean_field_constructor_and_backend": True}


def numerical_settings(initial_guess="minao"):
    return {"xc": "PBE0", "cart": False, "density_fitting": False,
            "grid_radial_angular": [150, 974], "grid_prune": None,
            "small_rho_cutoff": 0, "conv_tol": 1e-10,
            "conv_tol_grad": 1e-7, "direct_scf_tol": 1e-13,
            "max_cycle": 300, "initial_guess": initial_guess,
            "d4": "reuse archived printed scalar once"}


def scf_state(mf, output, label, *, dm0=None):
    """One declared SCF call, with native output and incremental cycle ledger."""
    import numpy as np
    start = time.monotonic()
    history_path = output / f"{label}_scf_history.jsonl"

    def callback(env):
        record = {"cycle": int(env["cycle"]),
                  "elapsed_seconds": time.monotonic() - start}
        for name in ("e_tot", "norm_gorb", "norm_ddm"):
            if name in env:
                record[name] = float(env[name])
        with history_path.open("a") as stream:
            stream.write(json.dumps(record, allow_nan=False) + "\n")
    mf.callback = callback
    CALL_COUNTS["scf_kernel"] += 1
    energy = mf.kernel(dm0=dm0)
    dm = mf.make_rdm1()
    overlap = mf.get_ovlp()
    fock = mf.get_fock(dm=dm)
    orbital_gradient = mf.get_grad(mf.mo_coeff, mf.mo_occ, fock)
    occ_coeff = mf.mo_coeff[:, mf.mo_occ > 1.0]
    identity_error = np.max(np.abs(occ_coeff.T @ overlap @ occ_coeff -
                                   np.eye(occ_coeff.shape[1])))
    electron_trace = np.einsum("ij,ji", dm, overlap)
    idempotency_error = np.linalg.norm(dm @ overlap @ dm - 2 * dm)
    gradient_norm = np.linalg.norm(orbital_gradient)
    state_path = output / f"{label}_state.npz"
    np.savez_compressed(state_path, mo_coeff=mf.mo_coeff, mo_occ=mf.mo_occ,
                        mo_energy=mf.mo_energy, density=dm,
                        e_tot=np.array(float(energy)), overlap=overlap)
    checks = {"native_scf_converged": bool(mf.converged),
              "finite_energy_and_density": bool(np.isfinite(energy) and np.isfinite(dm).all()),
              "electron_trace": bool(abs(electron_trace - mf.mol.nelectron) < 1e-7),
              "occupied_orthonormality": bool(identity_error < 1e-8),
              "density_idempotency": bool(idempotency_error < 1e-6),
              "orbital_gradient": bool(gradient_norm <= mf.conv_tol_grad)}
    report = {"status": "passed" if all(checks.values()) else "failed_scf_admission",
              "energy_hartree": float(energy), "checks": checks,
              "initial_guess": mf.init_guess,
              "orbital_gradient_norm": float(gradient_norm),
              "electron_trace": float(electron_trace),
              "occupied_orthonormality_max_error": float(identity_error),
              "density_idempotency_frobenius_error": float(idempotency_error),
              "nelectron": mf.mol.nelectron, "nao": mf.mol.nao_nr(),
              "grid_points": len(mf.grids.weights),
              "elapsed_seconds": time.monotonic() - start,
              "state_file": str(state_path),
              "state_sha256": hashlib.sha256(state_path.read_bytes()).hexdigest()}
    write_json(output / f"{label}_RESULT.json", report)
    if not all(checks.values()):
        raise RuntimeError(f"{label} SCF admission failed: {checks}")
    return report


def baseline(args, output):
    import archive
    mol, info = archive.build_molecule(args.system, output / "baseline_pyscf.log",
                                      args.max_memory_mb)
    write_json(output / "INPUT.json", info)
    mf = configure_mean_field(mol, output / "baseline.chk")
    mf.init_guess = args.initial_guess
    report = scf_state(mf, output, "baseline")
    reference = info["archive_energies"]
    discrepancy = report["energy_hartree"] - reference["electronic_hartree"]
    composite = report["energy_hartree"] + reference["d4_printed_hartree"]
    report["archive_comparison"] = {
        "electronic_difference_hartree": discrepancy,
        "electronic_difference_kcal_mol": discrepancy * archive.HARTREE_TO_KCAL,
        "pbe0_plus_archived_printed_d4_hartree": composite,
        "full_difference_hartree": composite - reference["full_pbe0_d4_hartree"],
        "full_energy_tolerance_hartree": 1e-4,
        "electronic_agreement": abs(discrepancy) <= 1e-4,
        "full_energy_agreement": abs(composite - reference["full_pbe0_d4_hartree"]) <= 1e-4,
        "d4_recomputed": False,
        "numerical_scheme_note": "Direct PySCF J/K and fixed unpruned 150x974 atom grid differ from archived ORCA RIJCOSX/DefGrid3; agreement is measured, not assumed."}
    report["checks"]["archived_full_energy_agreement"] = report["archive_comparison"]["full_energy_agreement"]
    report["checks"]["archived_electronic_energy_agreement"] = report["archive_comparison"]["electronic_agreement"]
    report["status"] = "passed" if all(report["checks"].values()) else "failed_archive_agreement"
    write_json(output / "baseline_RESULT.json", report)
    return report


def fragment(args, output):
    import archive
    if args.fragment not in (0, 1):
        raise ValueError("Fragment mode requires --fragment 0 or 1")
    mol, info = archive.build_molecule(args.system, output / "fragment_pyscf.log",
                                      args.max_memory_mb, fragment=args.fragment)
    info["physics_signature"] = physics_signature(info)
    write_json(output / "INPUT.json", info)
    mf = configure_mean_field(mol, output / "fragment.chk")
    mf.init_guess = args.initial_guess
    report = scf_state(mf, output, "fragment")
    report["fragment"] = info["runtime_subset"]
    report["physics_signature"] = info["physics_signature"]
    report["fragment_dispersion"] = "not evaluated; R is the isolated electronic reference"
    write_json(output / "fragment_RESULT.json", report)
    return report


def constrained(args, output):
    import archive
    import almo
    import numpy as np
    if args.baseline_dir is None or args.fragment_dirs is None:
        raise ValueError("ALMO mode requires an admitted baseline and two fragment directories")
    mol, info = archive.build_molecule(args.system, output / "almo_pyscf.log",
                                      args.max_memory_mb)
    write_json(output / "INPUT.json", info)
    mf = configure_mean_field(mol, output / "unused_full_state.chk")
    reuse = {"baseline": restore_state(args.baseline_dir, "baseline", mf, info), "fragments": []}
    fragment_mfs = []
    for i, source in enumerate(args.fragment_dirs):
        frag_mol, frag_info = archive.build_molecule(
            args.system, output / f"fragment_{i}_reuse_pyscf.log",
            args.max_memory_mb, fragment=i)
        frag_mf = configure_mean_field(frag_mol, output / f"unused_fragment_{i}.chk")
        reuse["fragments"].append(restore_state(source, "fragment", frag_mf, frag_info))
        fragment_mfs.append(frag_mf)
    write_json(output / "REUSED_STATES.json", reuse)
    atom_slices = mol.aoslice_by_atom()
    fragment_aos = [[j for i in f["atom_indices"]
                     for j in range(int(atom_slices[i, 2]), int(atom_slices[i, 3]))]
                    for f in info["fragments"]]
    fragment_nocc = [f["nelectron"] // 2 for f in info["fragments"]]
    problem = almo.ALMOProblem(mf, fragment_aos, fragment_nocc)
    if len(mf.grids.weights) != reuse["baseline"]["grid_points"]:
        raise ValueError("ALMO full-system grid count differs from saved baseline")
    grid_record = {"points": len(mf.grids.weights),
                   "coordinates_sha256": hashlib.sha256(mf.grids.coords.tobytes()).hexdigest(),
                   "weights_sha256": hashlib.sha256(mf.grids.weights.tobytes()).hexdigest(),
                   "baseline_grid_identity_evidence": "Identical inputs, pinned backend, exact constructor AST, no pruning, same native build call and point count; baseline did not archive point hashes."}
    write_json(output / "FIXED_GRID.json", grid_record)
    initial = problem.initial_coefficients(fragment_mfs=fragment_mfs)
    metric = initial.T @ problem.overlap @ initial
    frozen_evaluation = problem.evaluate(initial)
    initial_density = frozen_evaluation.density
    np.savez_compressed(output / "frozen_state.npz", coefficients=initial,
                        density=initial_density, occupied_overlap=metric,
                        energy_hartree=np.array(frozen_evaluation.energy))

    def progress(row):
        with (output / "almo_history.jsonl").open("a") as stream:
            stream.write(json.dumps(row, allow_nan=False) + "\n")

    solution = almo.solve_almo(mf, fragment_aos, fragment_nocc,
                               initial_coefficients=initial,
                               max_iterations=300, max_evaluations=1500,
                               gradient_tolerance=1e-7, energy_tolerance=1e-10,
                               callback=progress)
    np.savez_compressed(output / "polarized_state.npz", coefficients=solution.coefficients,
                        density=solution.density, coefficient_gradient=solution.coefficient_gradient)
    summary = solution.summary()
    write_json(output / "ALMO_SOLVER_RESULT.json", summary)
    R = sum(f["energy_hartree"] for f in reuse["fragments"])
    F, P, U = solution.initial_energy, solution.energy, reuse["baseline"]["energy_hartree"]
    frz, pol, ct = F - R, P - F, U - P
    closure = (frz + pol + ct) - (U - R)
    checks = {"solver_converged": bool(solution.converged),
              "frozen_saved_state_matches_solver": abs(frozen_evaluation.energy - F) <= 1e-10,
              "density_and_support_invariants": bool(solution.invariants["passed"]),
              "polarization_variational_order": P <= F + 1e-8,
              "full_relaxation_variational_order": U <= P + 1e-8,
              "energy_decomposition_closure": abs(closure) <= 1e-10}
    checks["fixed_grid_preserved"] = (
        hashlib.sha256(mf.grids.coords.tobytes()).hexdigest() == grid_record["coordinates_sha256"]
        and hashlib.sha256(mf.grids.weights.tobytes()).hexdigest() == grid_record["weights_sha256"])
    d4 = info["archive_energies"]["d4_printed_hartree"]
    report = {"status": "passed" if all(checks.values()) else "failed_almo_admission",
              "checks": checks, "solver": summary,
              "energies_hartree": {"R": R, "F": F, "P": P, "U": U,
                                    "FRZ": frz, "POL": pol, "CT": ct,
                                    "closure_residual": closure,
                                    "archived_printed_D4": d4,
                                    "full_PBE0_D4": U + d4,
                                    "CT_suppressed_PBE0_D4": P + d4},
              "reused_states": reuse, "fragment_aos": fragment_aos,
              "fixed_grid": grid_record,
              "fragment_nocc": fragment_nocc,
              "new_scf_calls": 0, "almo_backend_evaluations": solution.evaluations,
              "frozen_state_backend_evaluations": 1,
              "total_constrained_backend_evaluations": solution.evaluations + 1,
              "saved_frozen_energy_minus_solver_initial_hartree": frozen_evaluation.energy - F,
              "frozen_normalization_note": "Saved frozen orbitals have their own evaluated energy; solver re-normalization preserves the occupied subspace, with the residual checked explicitly.",
              "frozen_state": archive.pin(output / "frozen_state.npz"),
              "polarized_state": archive.pin(output / "polarized_state.npz"),
              "interpretation_limits": [
                  "AO-block variational stationary state; global constrained minimum not certified",
                  "Full U is the admitted saved ordinary SCF, not a repeated release calculation",
                  "No electrostatics/Pauli split; FRZ is a combined frozen term",
                  "No isolated-fragment D4 is invented; R and decomposition components are electronic",
                  "Effective singlet model is not a qualification of physical Dy multiplet or method accuracy"]}
    write_json(output / "almo_RESULT.json", report)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--system", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--max-memory-mb", type=int, default=28000)
    parser.add_argument("--mode", choices=("baseline", "fragment", "almo"), default="baseline")
    parser.add_argument("--fragment", type=int, choices=(0, 1))
    parser.add_argument("--baseline-dir", type=Path)
    parser.add_argument("--fragment-dirs", type=Path, nargs=2)
    parser.add_argument("--initial-guess", choices=("minao", "1e"), default="minao")
    parser.add_argument("--initialization-recovery-from", type=Path)
    args = parser.parse_args()
    if args.repo:
        os.environ["LADY_REPO_ROOT"] = str(args.repo.resolve())
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Numerical work requires a root-submitted Slurm allocation")
    allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
    if not 1 <= args.threads <= 24 or args.threads > allocated:
        raise RuntimeError("Requested threads exceed the actual scheduler grant")
    granted_memory = int(os.environ.get("SLURM_MEM_PER_NODE", "0"))
    if args.max_memory_mb <= 0 or (granted_memory and args.max_memory_mb >= granted_memory):
        raise RuntimeError("Mean-field memory must fit inside the actual scheduler grant")
    if os.environ.get("SLURM_JOB_QOS") not in (None, "standard_expand"):
        raise RuntimeError("This diagnostic requires standard_expand")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    start = time.monotonic()
    result = {"status": "started", "system_id": args.system, "mode": args.mode,
              "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "job_id": os.environ["SLURM_JOB_ID"], "allocated_cpus": allocated,
              "requested_threads": args.threads, "max_memory_mb": args.max_memory_mb,
              "new_scf_calls": 0}
    exit_code = 1
    try:
        from pyscf import lib
        import archive
        if args.system not in archive.SYSTEMS:
            raise ValueError("Outside exact declared quartet")
        lib.num_threads(args.threads)
        result["runtime"] = metadata_versions()
        result["source_pins"] = [archive.pin(__file__), archive.pin(archive.__file__),
                                  archive.pin(archive.MANIFEST)]
        result["settings"] = numerical_settings(args.initial_guess)
        result["initialization"] = prepare_initialization(args, result["runtime"])
        if args.mode == "almo":
            import almo
            result["source_pins"].append(archive.pin(almo.__file__))
            result["almo_settings"] = {"max_iterations": 300, "max_evaluations": 1500,
                                       "gradient_tolerance": 1e-7, "energy_tolerance": 1e-10,
                                       "variational_order_tolerance_hartree": 1e-8,
                                       "decomposition_closure_tolerance_hartree": 1e-10}
        write_json(output / "START.json", result)
        runner = {"baseline": baseline, "fragment": fragment, "almo": constrained}[args.mode]
        result[args.mode] = runner(args, output)
        result["status"] = result[args.mode]["status"]
        exit_code = 0 if result["status"] == "passed" else 2
    except Exception as error:
        result["status"] = "failed"
        result["error"] = repr(error)
        result["traceback"] = traceback.format_exc()
        print(result["traceback"], file=sys.stderr)
    finally:
        result["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        result["elapsed_seconds"] = time.monotonic() - start
        result["new_scf_calls"] = CALL_COUNTS["scf_kernel"]
        result["allocated_cpu_seconds"] = allocated * result["elapsed_seconds"]
        result["process_cpu_seconds"] = resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime
        result["peak_process_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result["physical_state_qualified"] = False
        result["method_accuracy_qualified"] = False
        write_json(output / "RESULT.json", result)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
