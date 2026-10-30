#!/usr/bin/env python3
"""Three declared XC evaluations on the exact saved La-aquo density; zero SCF.

This diagnoses a cross-code energy mismatch. It neither admits the failed
baseline nor changes a production/default Hamiltonian. The documented ORCA PBE
parameter substitution is not a claim of complete ORCA XC equivalence.
"""
from __future__ import annotations

import argparse
import ast
import ctypes
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback


ORCA_DEFAULTS = {101: {"_kappa": 0.804, "_mu": 0.219520}, 130: {"_beta": 0.066725}}
LIBXC_DEFAULTS = {
    101: {"_kappa": 0.804, "_mu": 0.2195149727645171},
    130: {"_beta": 0.06672455060314922, "_gamma": 0.031090690869654895034, "_B": 1.0},
}
EXPLICIT_PBE0 = "0.25*HF + 0.75*GGA_X_PBE, GGA_C_PBE"
# RKS.do_nlc calls parse_dft, which lowercases names before LibXC lookup.
CUSTOM_CODE = "lady_orca_documented_pbe_parameters_v1"
SOURCE_EVIDENCE = {
    "orca_defaults": "https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html#customization-of-basic-functionals",
    "orca_table": "Table 3.3: XMuePBE=0.21952, XKappa=0.804, CBetaPBE=0.066725",
    "libxc_exchange": "https://gitlab.com/libxc/libxc/-/raw/7.0.0/src/gga_x_pbe.c",
    "libxc_correlation": "https://gitlab.com/libxc/libxc/-/raw/7.0.0/src/gga_c_pbe.c",
    "libxc_mu": "https://gitlab.com/libxc/libxc/-/raw/7.0.0/src/util.h",
    "libxc_pbe0": "https://gitlab.com/libxc/libxc/-/raw/7.0.0/src/hyb_gga_xc_pbeh.c",
    "limitation": "ORCA PWLDA internal constants are not proven identical to LibXC modified PW92; this is only the documented mu/kappa/beta substitution.",
}

# Exact immutable artifacts from the observed serialization-only failure. This
# recovery is deliberately restricted to that completed first evaluation.
FIRST_EVALUATION_RECOVERY_PINS = {
    "CONFIG.json": "8db44f62971cf97a90ab25b962ac776670920151f92eba8378a2efb3e27a046e",
    "SUBMISSION.json": "f8983d54c52e544dfb98a61aa768d285cb22e8890ed510d5409e928e1230cd51",
    "slurm-22416.out": "b22fa36f89ec812ebd6f038b7d7e7fa8db0a9cdfd8c6ba042ab9b1f2d664bc7f",
    "source/diagnose_xc.py": "84453aaff4ec1780e17e464abf7499a258c981215a58e9ae788f193716957620",
    "result/START.json": "9f96ea9a70ad397730ab7215e3d50ddcb8e980bdfab23f7bcbf536423ac16f15",
    "result/INPUT.json": "30e704002bd814f54c8db669db9042245a6069102ba76c7df65de85c2f5cb883",
    "result/FIXED_INPUT_ADMISSION.json": "2516d4a22d96e5644a6af1c33e91b23b640b03f67d52fa13677f78dd41a04ed1",
    "result/builtin_PBE0.json": "7168bc98ffc1ed3b63e5a4b7d6de92b84d0ea7bf2e71289c36ec1b8266e2902b",
    "result/FIXED_DENSITY_ARRAYS.npz": "eea423a7447c35f23dcb770bf853a99a7e2f43388b2fae90eaa3359caa1b6fdf",
}
FIRST_EVALUATION_TERMINAL = {
    "job_id": "22416", "state": "FAILED", "exit_code": "1:0",
    "start": "2026-10-23T22:45:45", "end": "2026-10-23T22:47:14",
    "elapsed_seconds": "89", "allocated_cpus": "8", "nodes": "node-64-768g-14",
}
SECOND_EVALUATION_RECOVERY_PINS = {
    "CONFIG.json": "fc5c1ac1dde39a0e9ff71f8c87d5ba6bf2905360facc16b809a3348bea09c0f2",
    "SUBMISSION.json": "36b2fc62ebabe4505306174bcb76cdb15855276963cf94e02ee8a11591795c42",
    "slurm-22419.out": "95ffec6ec5a6711560c5620efff4b6b6a9b5d79779b2874b78018e2076ef2c4c",
    "source/diagnose_xc.py": "b3cf5524d412ac6fef8fb59e43f57face5f57441f59cd20ccc8cf817f8c95d24",
    "result/START.json": "00bdc3de24cc6fea51c2a43fb9b0dae69d0cf1e2efcc4e3e3823a474ca4691ef",
    "result/INPUT.json": "30e704002bd814f54c8db669db9042245a6069102ba76c7df65de85c2f5cb883",
    "result/FIXED_INPUT_ADMISSION.json": "6ba84b60f1f94b7c7f60093ec297059dd8a3be971c249b6c220fbf1b3e80a46d",
    "result/RECOVERY_ADMISSION.json": "19b4cbfcce7d65a50ef2c4569d76ceb83be9b5214f08ba4597150b3f9b7f79c8",
    "result/RESULT.json": "a16bb87fae6f8d0e77aae8e03feed72152f1f7eacb025701faa9f54ead2b11d3",
    "result/explicit_native_components.json": "a7e1d248e249ca9b68d29edfca6007836c8fe19659dd29738414cf3cb51ad440",
    "result/builtin_PBE0.json": "7168bc98ffc1ed3b63e5a4b7d6de92b84d0ea7bf2e71289c36ec1b8266e2902b",
    "result/FIXED_DENSITY_ARRAYS.npz": "d605fe3c98c5c7c13a80dfba281dae6ae570abf32690f0131e6030f3693aedcc",
}
SECOND_EVALUATION_TERMINAL = {
    "job_id": "22419", "state": "FAILED", "exit_code": "1:0",
    "start": "2026-10-23T22:53:03", "end": "2026-10-23T22:55:17",
    "elapsed_seconds": "134", "allocated_cpus": "8", "nodes": "node-64-768g-4",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def source_ast(text, name):
    node = next(n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name == name)
    return ast.dump(node, include_attributes=False)


def functional_dispatch_metadata(libxc, code):
    """Metadata-only real dispatcher smoke; no Mole, integrals, density or SCF."""
    from types import SimpleNamespace
    from pyscf.dft import rks
    from pyscf.scf.dispersion import parse_dft
    parsed = parse_dft(code)
    require(parsed[0] == code == code.lower(), "Custom registration must survive parse_dft normalization")

    def inspect(name):
        return {"xc_type": libxc.xc_type(name), "is_nlc": bool(libxc.is_nlc(name)),
                "is_hybrid_xc": bool(libxc.is_hybrid_xc(name)),
                "hybrid_coeff": float(libxc.hybrid_coeff(name)),
                "rsh_coeff": [float(value) for value in libxc.rsh_coeff(name)]}

    raw, normalized = inspect(code), inspect(parsed[0])
    require(raw == normalized and raw["xc_type"] == "GGA" and not raw["is_nlc"] and
            raw["is_hybrid_xc"] and raw["hybrid_coeff"] == 0.25 and
            raw["rsh_coeff"] == [float(value) for value in libxc.rsh_coeff(EXPLICIT_PBE0)],
            "Custom functional dispatcher metadata differs from the intended hybrid")
    lightweight = SimpleNamespace(xc=code, nlc="", _numint=SimpleNamespace(libxc=libxc))
    do_nlc = bool(rks.KohnShamDFT.do_nlc(lightweight))
    require(not do_nlc, "Real RKS.do_nlc dispatcher must resolve the custom name without NLC")
    return {"custom_code": code, "parse_dft": list(parsed), "raw_lookup": raw,
            "parsed_lookup": normalized, "RKS_do_nlc": do_nlc,
            "molecular_integrals": 0, "density_evaluations": 0, "scf_calls": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume-from-first-evaluation", type=Path,
                        help="Exact original 22416 result directory; reuse its pinned builtin evaluation and run only the remaining two")
    parser.add_argument("--resume-from-second-evaluation", type=Path,
                        help="Exact original 22419 result directory; also reuse its explicit-components evaluation and run only the last variant")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--max-memory-mb", type=int, default=28000)
    args = parser.parse_args()
    require(not args.resume_from_second_evaluation or args.resume_from_first_evaluation,
            "Second-evaluation reuse requires the independently verified first-evaluation chain")
    require(os.environ.get("SLURM_JOB_ID"), "Diagnostic requires a real Slurm allocation")
    require(os.environ.get("SLURM_JOB_QOS") in (None, "standard_expand"), "Require standard_expand")
    allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
    require(1 <= args.threads <= 8 and args.threads <= allocated, "Requested threads exceed declared/actual allocation")
    granted_memory = int(os.environ.get("SLURM_MEM_PER_NODE", "0"))
    require(0 < args.max_memory_mb <= 28000 and (not granted_memory or args.max_memory_mb < granted_memory),
            "Memory exceeds declared/actual allocation")
    if args.repo:
        os.environ["LADY_REPO_ROOT"] = str(args.repo.resolve())
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[name] = "1"
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    result = {
        "status": "failed", "schema": "saved_density_pbe_parameter_diagnostic_v1",
        "started_utc": stamp(), "job_id": os.environ["SLURM_JOB_ID"],
        "allocated_cpus": allocated, "requested_threads": args.threads,
        "source_evidence": SOURCE_EVIDENCE, "new_scf_calls": 0,
        "full_potential_evaluations_attempted": 0, "variants": [], "checks": {},
        "full_potential_evaluations_reused": 0,
        "resume_from_first_evaluation": str(args.resume_from_first_evaluation.resolve()) if args.resume_from_first_evaluation else None,
        "resume_from_second_evaluation": str(args.resume_from_second_evaluation.resolve()) if args.resume_from_second_evaluation else None,
        "previous_attempt_accounting": [],
        "custom_name_compatibility": {
            "previous_name": "LADY_ORCA_DOCUMENTED_PBE_PARAMETERS_V1",
            "registered_name": CUSTOM_CODE,
            "reason": "RKS.do_nlc uses parse_dft, whose normalized name is lowercase; LibXC custom-name lookup is case-sensitive",
            "functional_definition_changed": False,
        },
        "baseline_admitted": False, "production_or_default_changed": False,
        "density_relaxed_for_substituted_parameters": False,
        "full_orca_xc_equivalence_established": False,
    }
    write_json = None
    try:
        import numpy as np
        import pyscf
        from pyscf import lib
        from pyscf.dft import libxc
        import archive
        import run_system
        write_json = run_system.write_json
        lib.num_threads(args.threads)
        result["runtime"] = run_system.metadata_versions()
        result["source_pins"] = [archive.pin(__file__), archive.pin(archive.__file__),
                                 archive.pin(run_system.__file__), archive.pin(archive.MANIFEST),
                                 archive.pin(libxc.__file__), archive.pin(libxc._itrf._name)]
        libxc_library = Path(pyscf.__file__).parent / "lib/deps/lib/libxc.so"
        result["source_pins"].append(archive.pin(libxc_library))
        result["source_pins"].append(archive.pin(Path(pyscf.__file__).parent / "scf/dispersion.py"))
        result["parameters"] = {"libxc_defaults": LIBXC_DEFAULTS,
                                "documented_orca_substitution": ORCA_DEFAULTS,
                                "explicit_pbe0": EXPLICIT_PBE0,
                                "hybrid_fraction": 0.25,
                                "unchanged_correlation_parameters": ["_gamma", "_B"],
                                "parameter_array_order": {101: ["_kappa", "_mu"],
                                                           130: ["_beta", "_gamma", "_B"]}}
        result["thresholds"] = {"original_density_energy_replay_hartree": 1e-8,
                                "explicit_vs_builtin_energy_hartree": 1e-9,
                                "explicit_vs_builtin_fock_max_abs": 1e-9,
                                "original_density_orbital_gradient_norm": 1e-7}
        write_json(output / "START.json", result)

        saved = args.baseline_dir.resolve()
        require(saved != output and saved not in output.parents and output not in saved.parents,
                "Diagnostic and baseline artifact directories must be disjoint")
        run = json.loads((saved / "RESULT.json").read_text())
        report = json.loads((saved / "baseline_RESULT.json").read_text())
        previous_info = json.loads((saved / "INPUT.json").read_text())
        start = json.loads((saved / "START.json").read_text())
        require(run["system_id"] == "aquo8_central_La" and run["mode"] == "baseline",
                "Only the exact saved La-aquo baseline is declared")
        require(run["status"] in ("passed", "failed_archive_agreement") and
                report["status"] == run["status"], "Baseline is not a terminal native-SCF result")
        require(run["baseline"] == report, "Baseline parent and per-state reports differ")
        native_checks = {"native_scf_converged", "finite_energy_and_density", "electron_trace",
                         "occupied_orthonormality", "density_idempotency", "orbital_gradient"}
        parity_checks = {"archived_full_energy_agreement", "archived_electronic_energy_agreement"}
        require(set(report["checks"]) == native_checks | parity_checks, "Unexpected baseline admission schema")
        require(all(report["checks"][key] is True for key in native_checks), "Saved native-SCF gate failed")
        result["baseline_admission_exception"] = {
            "purpose": "diagnose archive disagreement; baseline remains unadmitted",
            "only_permitted_false_checks": sorted(parity_checks),
            "actual_checks": report["checks"], "native_scf_admitted": True,
            "original_status": report["status"],
        }
        for record in start["source_pins"] + start["runtime"]["pyscf_module_pins"]:
            archive.verify_pin(record)
        source_record = next(p for p in start["source_pins"] if Path(p["path"]).name == "run_system.py")
        require(source_ast(Path(source_record["path"]).read_text(), "configure_mean_field") ==
                source_ast(Path(run_system.__file__).read_text(), "configure_mean_field"),
                "Mean-field/grid constructor differs from original baseline")
        require(start["settings"] == run_system.numerical_settings(), "Baseline numerical settings differ")
        for key in ("packages", "libxc_version"):
            require(start["runtime"][key] == result["runtime"][key], "Backend version differs: " + key)
        result["input_pins"] = [archive.pin(saved / filename) for filename in
                                ("RESULT.json", "baseline_RESULT.json", "INPUT.json", "START.json", "baseline_state.npz")]
        require(archive.sha256(saved / "baseline_state.npz") == report["state_sha256"], "Saved-state hash differs")
        mol, info = archive.build_molecule("aquo8_central_La", output / "pyscf.log", args.max_memory_mb)
        require(run_system.physics_signature(info) == run_system.physics_signature(previous_info),
                "Exact saved geometry/basis/ECP/charge/state differs")
        write_json(output / "INPUT.json", info)
        mf = run_system.configure_mean_field(mol, output / "unused_no_scf.chk")
        with np.load(saved / "baseline_state.npz", allow_pickle=False) as state:
            coefficients, occupations = state["mo_coeff"].copy(), state["mo_occ"].copy()
            orbital_energies, density = state["mo_energy"].copy(), state["density"].copy()
            saved_overlap, saved_energy = state["overlap"].copy(), float(state["e_tot"])
        require(coefficients.shape == (mol.nao_nr(), mol.nao_nr()), "Saved coefficient shape mismatch")
        require(np.isfinite(coefficients).all() and np.isfinite(density).all(), "Nonfinite saved state")
        require(np.all((occupations == 0) | (occupations == 2)) and int(occupations.sum()) == mol.nelectron,
                "Saved occupations differ")
        require(saved_energy == report["energy_hartree"], "Saved/report energy mismatch")
        overlap, hcore = mf.get_ovlp(), mf.get_hcore()
        require(np.allclose(overlap, saved_overlap, rtol=1e-11, atol=1e-12), "Saved overlap differs")
        mf.mo_coeff, mf.mo_occ, mf.mo_energy, mf.e_tot = coefficients, occupations, orbital_energies, saved_energy
        require(np.allclose(mf.make_rdm1(), density, rtol=1e-10, atol=1e-11), "Saved density/orbitals disagree")
        occupied = coefficients[:, occupations == 2]
        density_record = {
            "sha256": hashlib.sha256(density.tobytes()).hexdigest(),
            "dtype": str(density.dtype), "shape": list(density.shape),
            "electron_trace": float(np.einsum("ij,ji->", density, overlap)),
            "idempotency_frobenius": float(np.linalg.norm(density @ overlap @ density - 2 * density)),
            "occupied_orthonormality_max_abs": float(np.max(np.abs(occupied.T @ overlap @ occupied - np.eye(len(occupied.T))))),
        }
        require(abs(density_record["electron_trace"] - mol.nelectron) < 1e-7 and
                density_record["idempotency_frobenius"] < 1e-6 and
                density_record["occupied_orthonormality_max_abs"] < 1e-8,
                "Saved density invariants failed")
        result["density"] = density_record
        mf.grids.build(with_non0tab=True)
        require(len(mf.grids.weights) == report["grid_points"], "Reconstructed original grid count differs")
        grid_record = {
            "points": len(mf.grids.weights),
            "coordinates_sha256": hashlib.sha256(mf.grids.coords.tobytes()).hexdigest(),
            "weights_sha256": hashlib.sha256(mf.grids.weights.tobytes()).hexdigest(),
            "identity_evidence": "Identical input, constructor AST, pinned backend, unpruned grid, build call and original point count; original baseline did not save grid hashes.",
        }
        result["fixed_grid"] = grid_record
        write_json(output / "FIXED_INPUT_ADMISSION.json", result)

        # Runtime parameter metadata only: no evaluation of densities here.
        get_name = libxc._itrf.xc_func_info_get_ext_params_name
        get_name.argtypes, get_name.restype = (ctypes.c_void_p, ctypes.c_int), ctypes.c_char_p
        get_value = libxc._itrf.xc_func_get_ext_params_value
        get_value.argtypes, get_value.restype = (ctypes.c_void_p, ctypes.c_int), ctypes.c_double

        def parameter_metadata(code):
            cache = libxc._get_xc(code, spin=0)
            entries = {}
            for xid, pointer in cache.obj_by_id().items():
                func_info = libxc._itrf.xc_func_get_info(pointer)
                count = libxc._itrf.xc_func_info_get_n_ext_params(func_info)
                entries[int(xid)] = {get_name(func_info, i).decode(): float(get_value(pointer, i)) for i in range(count)}
            return {"components": entries,
                    "component_factors": {int(xid): float(factor) for xid, factor in zip(cache.fn_ids, cache.facs)},
                    "hybrid_coeff": float(libxc.hybrid_coeff(code, spin=0)),
                    "xc_type": libxc.xc_type(code)}

        explicit_metadata = parameter_metadata(EXPLICIT_PBE0)
        require(explicit_metadata["components"] == LIBXC_DEFAULTS, "Installed LibXC defaults differ from declared source")
        require(explicit_metadata["hybrid_coeff"] == 0.25 and
                explicit_metadata["component_factors"] == {101: 0.75, 130: 1.0}, "Explicit PBE0 factors differ")
        libxc.register_custom_functional_(CUSTOM_CODE, EXPLICIT_PBE0,
                                         ext_params=ORCA_DEFAULTS, restricted=True, unrestricted=False)
        custom_metadata = parameter_metadata(CUSTOM_CODE)
        expected_custom = {xid: dict(values) for xid, values in LIBXC_DEFAULTS.items()}
        for xid, overrides in ORCA_DEFAULTS.items():
            expected_custom[xid].update(overrides)
        require(custom_metadata["components"] == expected_custom and
                custom_metadata["hybrid_coeff"] == 0.25 and
                custom_metadata["component_factors"] == {101: 0.75, 130: 1.0}, "Custom functional parameter readback differs")
        result["runtime_parameter_readback"] = {"builtin": parameter_metadata("PBE0"),
                                                 "explicit": explicit_metadata,
                                                 "substitution": custom_metadata}
        result["functional_dispatch_metadata"] = functional_dispatch_metadata(libxc, CUSTOM_CODE)
        # Catch metadata serialization defects before any expensive evaluation.
        json.dumps(result, sort_keys=True, allow_nan=False)
        arrays = {"density": density, "overlap": overlap, "occupied_coefficients": occupied}
        baseline_fock = baseline_gradient = None
        reference = info["archive_energies"]
        variants = [("builtin_PBE0", "PBE0"), ("explicit_native_components", EXPLICIT_PBE0),
                    ("documented_ORCA_parameter_substitution", CUSTOM_CODE)]
        if args.resume_from_first_evaluation:
            recovered = args.resume_from_first_evaluation.resolve()
            old_task = recovered.parent
            require(recovered.name == "result" and recovered != output and
                    recovered not in output.parents and output not in recovered.parents,
                    "Recovery and output paths must be disjoint original result directories")
            recovery_pins = []
            for relative, digest in FIRST_EVALUATION_RECOVERY_PINS.items():
                path = old_task / relative
                require(archive.sha256(path) == digest, "Original first-evaluation artifact changed: " + relative)
                recovery_pins.append(archive.pin(path))
            old_submission = json.loads((old_task / "SUBMISSION.json").read_text())
            old_config = json.loads((old_task / "CONFIG.json").read_text())
            terminal = json.loads((old_task / "completion_event.json").read_text())
            require(old_submission["job_id"] == "22416" and
                    old_submission["config_sha256"] == FIRST_EVALUATION_RECOVERY_PINS["CONFIG.json"] and
                    old_config["source_sha256"]["diagnose_xc.py"] == FIRST_EVALUATION_RECOVERY_PINS["source/diagnose_xc.py"],
                    "Original submission/source identity differs")
            require(terminal["records"]["22416"] == FIRST_EVALUATION_TERMINAL,
                    "Original failed terminal evidence differs")
            log = (old_task / "slurm-22416.out").read_text()
            require("TypeError: keys must be str, int, float, bool or None, not int32" in log and
                    'write_json(output / "PROGRESS.json", result)' in log,
                    "Original failure is not the admitted first-evaluation serialization failure")
            old_admission = json.loads((recovered / "FIXED_INPUT_ADMISSION.json").read_text())
            old_start = json.loads((recovered / "START.json").read_text())
            for record in old_admission["source_pins"] + old_admission["input_pins"]:
                archive.verify_pin(record)
            require(old_admission["input_pins"] == result["input_pins"],
                    "Recovery baseline inputs differ")
            require(old_admission["density"]["sha256"] == density_record["sha256"] and
                    old_admission["density"]["dtype"] == density_record["dtype"] and
                    old_admission["density"]["shape"] == density_record["shape"],
                    "Recovery density hash/type/shape differs")
            require(old_admission["fixed_grid"] == grid_record, "Recovery grid identity differs")
            require(old_admission["parameters"] == json.loads(json.dumps(result["parameters"])) and
                    old_admission["thresholds"] == result["thresholds"] and
                    old_admission["source_evidence"] == SOURCE_EVIDENCE,
                    "Recovery functional parameters or fixed gates differ")
            for key in ("packages", "libxc_version"):
                require(old_start["runtime"][key] == result["runtime"][key], "Recovery backend differs: " + key)
            require(old_start["job_id"] == "22416" and old_admission["job_id"] == "22416" and
                    old_admission["new_scf_calls"] == 0 and old_admission["variants"] == [] and
                    old_admission["full_potential_evaluations_attempted"] == 0,
                    "Recovery initial-admission provenance differs")
            row = json.loads((recovered / "builtin_PBE0.json").read_text())
            require(row["label"] == "builtin_PBE0" and row["xc_code"] == "PBE0" and
                    row["new_scf_calls"] == 0 and row["archived_D4_reused_once"] is True,
                    "Recovery first variant differs")
            with np.load(recovered / "FIXED_DENSITY_ARRAYS.npz", allow_pickle=False) as old_arrays:
                require(set(old_arrays.files) == {"density", "overlap", "occupied_coefficients",
                                                  "builtin_PBE0_fock", "builtin_PBE0_orbital_gradient"},
                        "Recovery archive is not exactly one completed variant")
                for key in ("density", "overlap", "occupied_coefficients"):
                    require(np.array_equal(old_arrays[key], arrays[key]), "Recovery array differs: " + key)
                baseline_fock = old_arrays["builtin_PBE0_fock"].copy()
                baseline_gradient = old_arrays["builtin_PBE0_orbital_gradient"].copy()
            require(np.isfinite(baseline_fock).all() and np.isfinite(baseline_gradient).all(),
                    "Recovery Fock/gradient is nonfinite")
            replay_gradient = mf.get_grad(coefficients, occupations, baseline_fock)
            require(np.max(np.abs(replay_gradient - baseline_gradient)) <= 1e-12 and
                    abs(float(np.linalg.norm(baseline_gradient)) - row["orbital_gradient_norm"]) <= 1e-12 and
                    abs(float(np.max(np.abs(baseline_gradient))) - row["orbital_gradient_max_abs"]) <= 1e-12,
                    "Recovery gradient/row does not reconstruct from the saved Fock")
            require(float(np.max(np.abs(baseline_fock - baseline_fock.T))) == row["fock_symmetry_max_abs"],
                    "Recovery Fock symmetry/row differs")
            algebraic_energy = (float(np.einsum("ij,ji->", density, hcore)) + float(mol.energy_nuc()) +
                                row["coulomb_energy_hartree"] + row["xc_energy_including_exact_exchange_hartree"])
            require(abs(algebraic_energy - row["energy_hartree"]) <= 1e-10 and
                    abs(row["semilocal_xc_energy_hartree"] + row["exact_exchange_energy_hartree"] -
                        row["xc_energy_including_exact_exchange_hartree"]) <= 1e-12,
                    "Recovery energy components do not close")
            require(row["difference_from_saved_energy_hartree"] == row["energy_hartree"] - saved_energy and
                    row["difference_from_ORCA_electronic_hartree"] == row["energy_hartree"] - reference["electronic_hartree"] and
                    row["plus_archived_printed_D4_hartree"] == row["energy_hartree"] + reference["d4_printed_hartree"] and
                    row["difference_from_ORCA_full_hartree"] == row["plus_archived_printed_D4_hartree"] - reference["full_pbe0_d4_hartree"],
                    "Recovery baseline/archive/D4 comparisons differ")
            result["checks"]["original_density_energy_replay"] = abs(row["energy_hartree"] - saved_energy) <= 1e-8
            result["checks"]["original_density_orbital_gradient"] = row["orbital_gradient_norm"] <= 1e-7
            require(all(result["checks"].values()), "Recovered original-density replay fails unchanged gates")
            result["full_potential_evaluations_reused"] = 1
            result["previous_attempt_accounting"].append({
                "job_id": "22416", "full_potential_evaluations_attempted": 1,
                "completed_new_evaluations": 1, "incomplete_potential_evaluations": 0,
                "failure": "JSON serialization after completed builtin evaluation"})
            result["variants"].append(row)
            arrays["builtin_PBE0_fock"] = baseline_fock
            arrays["builtin_PBE0_orbital_gradient"] = baseline_gradient
            result["recovery"] = {
                "status": "admitted_exact_first_evaluation", "original_failure_preserved": True,
                "original_job_id": "22416", "artifact_pins": recovery_pins,
                "terminal_record": FIRST_EVALUATION_TERMINAL,
                "terminal_receipt_observed": archive.pin(old_task / "completion_event.json"),
                "new_full_potential_evaluations_for_replay": 0,
                "new_scf_calls": 0, "reused_full_potential_evaluations": 1,
                "remaining_declared_variants": [pair[0] for pair in variants[1:]],
            }
            np.savez_compressed(output / "FIXED_DENSITY_ARRAYS.npz", **arrays)
            write_json(output / "builtin_PBE0.json", row)
            write_json(output / "RECOVERY_ADMISSION.json", result)
            variants = variants[1:]
        if args.resume_from_second_evaluation:
            recovered_second = args.resume_from_second_evaluation.resolve()
            second_task = recovered_second.parent
            require(recovered_second.name == "result" and recovered_second != output and
                    recovered_second not in output.parents and output not in recovered_second.parents,
                    "Second recovery and output paths must be disjoint")
            second_pins = []
            for relative, digest in SECOND_EVALUATION_RECOVERY_PINS.items():
                path = second_task / relative
                require(archive.sha256(path) == digest, "Second-evaluation artifact changed: " + relative)
                second_pins.append(archive.pin(path))
            second_submission = json.loads((second_task / "SUBMISSION.json").read_text())
            second_config = json.loads((second_task / "CONFIG.json").read_text())
            second_terminal = json.loads((second_task / "completion_event.json").read_text())
            require(second_submission["job_id"] == "22419" and
                    second_submission["config_sha256"] == SECOND_EVALUATION_RECOVERY_PINS["CONFIG.json"] and
                    second_config["source_sha256"]["diagnose_xc.py"] == SECOND_EVALUATION_RECOVERY_PINS["source/diagnose_xc.py"] and
                    second_terminal["records"]["22419"] == SECOND_EVALUATION_TERMINAL,
                    "Second recovery submission/source/terminal identity differs")
            previous = json.loads((recovered_second / "RESULT.json").read_text())
            second_row = json.loads((recovered_second / "explicit_native_components.json").read_text())
            require(previous["status"] == "failed" and previous["job_id"] == "22419" and
                    "ks.do_nlc()" in previous["error"] and
                    "LADY_ORCA_DOCUMENTED_PBE_PARAMETERS_V1' not found" in previous["error"],
                    "Second recovery is not the declared custom-name dispatch failure")
            require(previous["full_potential_evaluations_attempted"] == 2 and
                    previous["full_potential_evaluations_reused"] == 1 and
                    previous["new_scf_calls"] == 0 and len(previous["variants"]) == 2 and
                    previous["variants"][0] == result["variants"][0] and
                    previous["variants"][1] == second_row,
                    "Second recovery completed/failed evaluation accounting differs")
            expected_old_checks = {"original_density_energy_replay", "original_density_orbital_gradient",
                                   "explicit_native_energy_equivalence", "explicit_native_fock_equivalence"}
            require(set(previous["checks"]) == expected_old_checks and all(previous["checks"].values()),
                    "Second recovery did not pass the unchanged original and explicit-component gates")
            for record in previous["source_pins"] + previous["input_pins"]:
                archive.verify_pin(record)
            require(previous["input_pins"] == result["input_pins"] and
                    previous["fixed_grid"] == grid_record and
                    previous["density"]["sha256"] == density_record["sha256"] and
                    previous["parameters"] == json.loads(json.dumps(result["parameters"])) and
                    previous["thresholds"] == result["thresholds"] and
                    previous["runtime_parameter_readback"] == json.loads(json.dumps(result["runtime_parameter_readback"])) and
                    previous["source_evidence"] == SOURCE_EVIDENCE,
                    "Second recovery fixed inputs/parameters/gates differ")
            require(previous["recovery"]["original_job_id"] == "22416" and
                    previous["recovery"]["artifact_pins"] == result["recovery"]["artifact_pins"] and
                    previous["recovery"]["terminal_record"] == FIRST_EVALUATION_TERMINAL,
                    "Second recovery first-evaluation provenance chain differs")
            for key in ("packages", "libxc_version"):
                require(previous["runtime"][key] == result["runtime"][key], "Second recovery backend differs: " + key)
            with np.load(recovered_second / "FIXED_DENSITY_ARRAYS.npz", allow_pickle=False) as second_arrays:
                require(set(second_arrays.files) == set(arrays) | {"explicit_native_components_fock", "explicit_native_components_orbital_gradient"},
                        "Second recovery archive contains an unexpected evaluation panel")
                for key, value in arrays.items():
                    require(np.array_equal(second_arrays[key], value), "Second recovery inherited array differs: " + key)
                second_fock = second_arrays["explicit_native_components_fock"].copy()
                second_gradient = second_arrays["explicit_native_components_orbital_gradient"].copy()
            require(np.isfinite(second_fock).all() and np.isfinite(second_gradient).all(), "Second recovery arrays are nonfinite")
            require(second_row["label"] == "explicit_native_components" and second_row["xc_code"] == EXPLICIT_PBE0,
                    "Second recovery functional definition differs")
            second_gradient_replay = mf.get_grad(coefficients, occupations, second_fock)
            require(np.max(np.abs(second_gradient_replay - second_gradient)) <= 1e-12 and
                    abs(float(np.linalg.norm(second_gradient)) - second_row["orbital_gradient_norm"]) <= 1e-12 and
                    abs(float(np.linalg.norm(second_gradient - baseline_gradient)) - second_row["orbital_gradient_change_norm"]) <= 1e-12,
                    "Second recovery gradient does not reconstruct from saved Fock")
            require(float(np.max(np.abs(second_fock - baseline_fock))) == second_row["fock_change_max_abs"] and
                    abs(second_row["energy_hartree"] - result["variants"][0]["energy_hartree"]) <= 1e-9 and
                    second_row["fock_change_max_abs"] <= 1e-9,
                    "Second recovery explicit/native equivalence fails unchanged gates")
            second_algebraic_energy = (float(np.einsum("ij,ji->", density, hcore)) + float(mol.energy_nuc()) +
                                       second_row["coulomb_energy_hartree"] + second_row["xc_energy_including_exact_exchange_hartree"])
            require(abs(second_algebraic_energy - second_row["energy_hartree"]) <= 1e-10 and
                    second_row["difference_from_saved_energy_hartree"] == second_row["energy_hartree"] - saved_energy and
                    second_row["difference_from_ORCA_electronic_hartree"] == second_row["energy_hartree"] - reference["electronic_hartree"] and
                    second_row["plus_archived_printed_D4_hartree"] == second_row["energy_hartree"] + reference["d4_printed_hartree"] and
                    second_row["difference_from_ORCA_full_hartree"] == second_row["plus_archived_printed_D4_hartree"] - reference["full_pbe0_d4_hartree"],
                    "Second recovery energy/archive/D4 bookkeeping differs")
            result["checks"]["explicit_native_energy_equivalence"] = True
            result["checks"]["explicit_native_fock_equivalence"] = True
            result["full_potential_evaluations_reused"] = 2
            result["previous_attempt_accounting"].append({
                "job_id": "22419", "full_potential_evaluations_attempted": 2,
                "completed_new_evaluations": 1, "incomplete_potential_evaluations": 1,
                "failure": "Custom-name lookup in RKS.do_nlc during last get_veff; no completed custom energy/Fock"})
            result["variants"].append(second_row)
            arrays["explicit_native_components_fock"] = second_fock
            arrays["explicit_native_components_orbital_gradient"] = second_gradient
            result["second_recovery"] = {
                "status": "admitted_exact_second_evaluation", "original_failure_preserved": True,
                "original_job_id": "22419", "artifact_pins": second_pins,
                "terminal_record": SECOND_EVALUATION_TERMINAL,
                "terminal_receipt_observed": archive.pin(second_task / "completion_event.json"),
                "first_evaluation_chain_retained": True, "new_scf_calls": 0,
                "new_full_potential_evaluations_for_replay": 0,
                "remaining_declared_variants": ["documented_ORCA_parameter_substitution"],
            }
            np.savez_compressed(output / "FIXED_DENSITY_ARRAYS.npz", **arrays)
            write_json(output / "explicit_native_components.json", second_row)
            write_json(output / "SECOND_RECOVERY_ADMISSION.json", result)
            variants = variants[1:]
        for label, code in variants:
            mf.xc = code
            result["full_potential_evaluations_attempted"] += 1
            potential = mf.get_veff(mol, density, hermi=1)
            energy = float(mf.energy_tot(dm=density, h1e=hcore, vhf=potential))
            fock = hcore + np.asarray(potential)
            gradient = mf.get_grad(coefficients, occupations, fock)
            require(np.isfinite(energy) and np.isfinite(fock).all() and np.isfinite(gradient).all(),
                    "Nonfinite fixed-density evaluation: " + label)
            exact_exchange = float(-0.25 * np.einsum("ij,ji->", density, potential.vk))
            xc_energy = float(potential.exc)
            row = {
                "label": label, "xc_code": code, "energy_hartree": energy,
                "xc_energy_including_exact_exchange_hartree": xc_energy,
                "semilocal_xc_energy_hartree": xc_energy - exact_exchange,
                "exact_exchange_energy_hartree": exact_exchange,
                "coulomb_energy_hartree": float(potential.ecoul),
                "orbital_gradient_norm": float(np.linalg.norm(gradient)),
                "orbital_gradient_max_abs": float(np.max(np.abs(gradient))),
                "fock_symmetry_max_abs": float(np.max(np.abs(fock - fock.T))),
                "difference_from_saved_energy_hartree": energy - saved_energy,
                "difference_from_ORCA_electronic_hartree": energy - reference["electronic_hartree"],
                "plus_archived_printed_D4_hartree": energy + reference["d4_printed_hartree"],
                "difference_from_ORCA_full_hartree": energy + reference["d4_printed_hartree"] - reference["full_pbe0_d4_hartree"],
                "archived_D4_reused_once": True, "new_scf_calls": 0,
            }
            if baseline_fock is None:
                baseline_fock, baseline_gradient = fock.copy(), gradient.copy()
                result["checks"]["original_density_energy_replay"] = abs(energy - saved_energy) <= 1e-8
                result["checks"]["original_density_orbital_gradient"] = row["orbital_gradient_norm"] <= 1e-7
            else:
                row["energy_change_from_builtin_hartree"] = energy - result["variants"][0]["energy_hartree"]
                row["xc_change_from_builtin_hartree"] = xc_energy - result["variants"][0]["xc_energy_including_exact_exchange_hartree"]
                row["fock_change_max_abs"] = float(np.max(np.abs(fock - baseline_fock)))
                row["fock_change_frobenius"] = float(np.linalg.norm(fock - baseline_fock))
                row["orbital_gradient_change_norm"] = float(np.linalg.norm(gradient - baseline_gradient))
                if label == "explicit_native_components":
                    result["checks"]["explicit_native_energy_equivalence"] = abs(row["energy_change_from_builtin_hartree"]) <= 1e-9
                    result["checks"]["explicit_native_fock_equivalence"] = row["fock_change_max_abs"] <= 1e-9
            arrays[label + "_fock"] = fock.copy()
            arrays[label + "_orbital_gradient"] = gradient.copy()
            result["variants"].append(row)
            np.savez_compressed(output / "FIXED_DENSITY_ARRAYS.npz", **arrays)
            write_json(output / (label + ".json"), row)
            write_json(output / "PROGRESS.json", result)
            if label == "builtin_PBE0":
                require(all(result["checks"].values()), "Original-density replay gate failed; later variants not evaluated")
        result["checks"]["density_unchanged"] = hashlib.sha256(density.tobytes()).hexdigest() == density_record["sha256"]
        result["checks"]["grid_unchanged"] = (
            hashlib.sha256(mf.grids.coords.tobytes()).hexdigest() == grid_record["coordinates_sha256"] and
            hashlib.sha256(mf.grids.weights.tobytes()).hexdigest() == grid_record["weights_sha256"])
        result["checks"]["exactly_three_fixed_density_evaluations"] = (
            result["full_potential_evaluations_attempted"] + result["full_potential_evaluations_reused"] == 3 and
            len(result["variants"]) == 3)
        expected_new = 1 if args.resume_from_second_evaluation else (2 if args.resume_from_first_evaluation else 3)
        result["checks"]["declared_new_evaluation_count"] = result["full_potential_evaluations_attempted"] == expected_new
        result["array_archive"] = archive.pin(output / "FIXED_DENSITY_ARRAYS.npz")
        result["status"] = "passed_fixed_density_diagnostic" if all(result["checks"].values()) else "failed_diagnostic_gate"
        result["interpretation"] = [
            "A fixed-density energy response is diagnostic; no density was reconverged for the substituted constants.",
            "The original archive-parity failure is preserved; this test does not admit a corrected baseline or choose a functional.",
            "Even closer agreement would not establish full ORCA XC equivalence; internal PWLDA constants, quadrature and RIJCOSX remain distinct issues.",
        ]
        libxc.unregister_custom_functional_(CUSTOM_CODE)
    except Exception:
        result["error"] = traceback.format_exc()
        print(result["error"], file=sys.stderr)
    finally:
        result["finished_utc"] = stamp()
        result["elapsed_seconds"] = time.monotonic() - started
        result["allocated_cpu_seconds"] = allocated * result["elapsed_seconds"]
        usage = resource.getrusage(resource.RUSAGE_SELF)
        result["process_cpu_seconds"] = usage.ru_utime + usage.ru_stime
        result["peak_process_rss_kib"] = usage.ru_maxrss
        if write_json is None:
            (output / "RESULT.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        else:
            write_json(output / "RESULT.json", result)
    return 0 if result["status"] == "passed_fixed_density_diagnostic" else 1


if __name__ == "__main__":
    raise SystemExit(main())
