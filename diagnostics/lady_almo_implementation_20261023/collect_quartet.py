#!/usr/bin/env python3
"""Admit and collect four saved ALMO results; never execute molecular work.

Run in a root-submitted standard_expand Slurm job. A failed or absent
prerequisite produces a failed RESULTS.json without scientific contrasts.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
import datetime
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback

getcontext().prec = 40
SYSTEMS = ("aquo8_central_La", "aquo8_central_Dy",
           "hydrated_acetate_mono_central_La", "hydrated_acetate_mono_central_Dy")
SIGNS = dict(zip(SYSTEMS, (1, -1, -1, 1)))
STAGES = ("R", "F", "P", "U", "FRZ", "POL", "CT", "archived_printed_D4",
          "full_PBE0_D4", "CT_suppressed_PBE0_D4")
ROUNDING_GATE = Decimal("1e-10")
ORDER_GATE = Decimal("1e-8")
ARCHIVE_GATE = Decimal("1e-4")


def need(ok, message):
    if not ok:
        raise ValueError(message)


def dec(value):
    value = Decimal(str(value))
    need(value.is_finite(), "Nonfinite saved scalar")
    return value


def pin(path):
    path = Path(path).resolve()
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return {"path": str(path), "sha256": h.hexdigest()}


def verified(record):
    need(pin(record["path"])["sha256"] == record["sha256"],
         f"Pinned artifact changed: {record['path']}")
    return Path(record["path"])


def read(path):
    return json.loads(Path(path).read_text())


def all_checks(checks, label):
    need(isinstance(checks, dict) and checks, label + ": checks absent")
    need(all(value is True for value in checks.values()), label + ": failed/nonboolean check")


def list_checks(checks, label):
    need(isinstance(checks, list) and checks, label + ": checks absent")
    names = [row["name"] for row in checks]
    need(len(names) == len(set(names)), label + ": duplicate checks")
    need(all(row["passed"] is True for row in checks), label + ": failed check")
    return set(names)


def common_energy_settings(settings):
    expected = {"xc": "PBE0", "cart": False, "density_fitting": False,
                "grid_radial_angular": [150, 974], "grid_prune": None,
                "small_rho_cutoff": 0, "conv_tol": 1e-10, "conv_tol_grad": 1e-7,
                "direct_scf_tol": 1e-13, "max_cycle": 300,
                "d4": "reuse archived printed scalar once"}
    need(set(settings) == set(expected) | {"initial_guess"}, "Unknown/incomplete scientific settings")
    need(settings["initial_guess"] in ("minao", "1e"), "Unrecognized starting guess")
    common = {key: value for key, value in settings.items() if key != "initial_guess"}
    need(common == expected, "Hamiltonian, quadrature or SCF criteria changed")
    return common


def subset_signature(info):
    subset = info["runtime_subset"]
    atoms = [info["atoms"][i] for i in subset["atom_indices"]]
    elements = sorted({a[0] for a in atoms})
    payload = {"atoms_angstrom": atoms, "charge": subset["charge"],
               "spin": subset["spin"], "nelectron": subset["nelectron"],
               "basis": {e: info["basis"][e] for e in elements},
               "ecp": {e: info["ecp"][e] for e in elements if e in info["ecp"]},
               "cart": False}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def admit_initialization(run, info, stage=None):
    """Permit exactly the reviewed Dy55 native-1e initialization recovery."""
    common_energy_settings(run["settings"])
    guess = run["settings"]["initial_guess"]
    if stage is not None and "initial_guess" in stage:
        need(stage["initial_guess"] == guess, "Saved SCF actual/declared initialization differs")
    if guess == "minao":
        if "initialization" in run:
            need(run["initialization"]["selected_guess"] == "minao", "Starting guess provenance differs")
        return
    subset = info["runtime_subset"]
    need("Dy" in [info["atoms"][i][0] for i in subset["atom_indices"]]
         and info["ecp"].get("Dy", [None])[0] == 55, "Native 1e outside Dy55 recovery scope")
    need(run["mode"] in ("baseline", "fragment")
         and (run["mode"] != "fragment" or subset["fragment"] == 0), "Wrong native-1e recovery task")
    need(stage is not None and stage.get("initial_guess") == "1e", "Native-1e actual SCF record absent")
    record = run["initialization"]
    need(record["selected_guess"] == "1e"
         and record["recovery_kind"] == "Dy55_minao_rectangular_projection_to_native_1e", "Unknown initialization exception")
    previous = read(verified(record["previous_failure"]))
    previous_info = read(verified(record["previous_input"]))
    previous_start = read(verified(record["previous_start"]))
    verified(record["pyscf_hf_source"])
    need(previous["status"] == "failed" and previous["mode"] == run["mode"]
         and previous["system_id"] == run["system_id"], "Wrong failed initialization predecessor")
    need(previous["settings"]["initial_guess"] == "minao"
         and previous_start["settings"] == previous["settings"], "Failed predecessor was not original MINAO")
    need(common_energy_settings(previous["settings"]) == common_energy_settings(run["settings"]),
         "Recovery changed noninitialization settings")
    need(previous.get("error") == "LinAlgError('Last 2 dimensions of the array must be square')"
         and "minao_basis" in previous.get("traceback", "")
         and "s12[occ4ecp>0][:,occ>0]" in previous["traceback"], "Wrong initialization failure branch")
    need(subset_signature(previous_info) == subset_signature(info), "Recovery changed physical inputs")
    need(previous["runtime"]["packages"] == run["runtime"]["packages"], "Recovery changed electronic runtime")
    old_hf = [p for p in previous["runtime"]["pyscf_module_pins"] if p["path"].endswith("/scf/hf.py")]
    need(len(old_hf) == 1 and old_hf[0]["sha256"] == record["pyscf_hf_source"]["sha256"],
         "Reviewed PySCF initial-guess source changed")


def source_pin(run, name):
    matches = [r for r in run["source_pins"] if Path(r["path"]).name == name]
    need(len(matches) == 1, f"Expected one saved source pin for {name}")
    verified(matches[0])
    return matches[0]


def snapshot_inputs(info, expected):
    need(info["system_id"] == expected["system_id"], "Wrong input system ID")
    for key in ("atoms", "charge", "spin", "nelectron", "basis", "ecp", "fragments",
                "physical_multiplicity_metadata", "archive_energies"):
        need(info[key] == expected[key], f"Archive quartet input differs: {key}")
    for record in info["provenance"].values():
        verified(record)


def admit_reused_state(reuse, mode, system_id=None):
    for key in ("result", "run_result", "input", "state"):
        verified(reuse[key])
    run = read(reuse["run_result"]["path"])
    stage = read(reuse["result"]["path"])
    info = read(reuse["input"]["path"])
    need(run["status"] == "passed" and run["mode"] == mode, "Unadmitted reused " + mode)
    need(stage["status"] == "passed", "Failed reused stage")
    all_checks(stage["checks"], "Reused " + mode)
    admit_initialization(run, info, stage)
    need(run[mode] == stage, "Reused top-level and stage results disagree")
    need(reuse["state"]["sha256"] == stage["state_sha256"], "Reused state digest differs")
    need(dec(reuse["energy_hartree"]) == dec(stage["energy_hartree"]), "Reused energy differs")
    need(reuse["new_scf_calls"] == 0, "Reused state unexpectedly reran SCF")
    if system_id is not None:
        need(run["system_id"] == system_id and info["system_id"] == system_id,
             "Baseline belongs to another system")
    return run, stage, info


def admit_prerequisites(core_path, comparison_path):
    core = read(core_path)
    need(core["status"] == "passed_core_tests" and core["mode"] == "core", "Core tests did not pass")
    names = list_checks(core["checks"], "Core tests")
    required = {"AO_and_electron_count", "frozen_backend_energy",
                "orbital_gradient_two_step_FD", "FD_control_has_signal",
                "FD_rejects_factor_two_corruption", "nonsingular_within_block_occupied_gauge",
                "rank_deficient_rejected", "outside_fragment_rejected", "ALMO_converged",
                "optimized_backend_energy", "ordinary_full_SCF_converged", "variational_order",
                "one_fragment_reduces_to_ordinary_SCF", "release_constraints_ordinary_SCF"}
    need(required <= names, "Core validation coverage missing")
    need(core["method"]["xc"] == "PBE0", "Wrong core-test functional")
    for key in ("implementation", "almo_source", "fixture"):
        verified(core[key])
    comparison = read(comparison_path)
    need(comparison["status"] == "passed_independent_cp2k_pbe_comparison",
         "Independent CP2K comparison did not pass")
    all_checks(comparison["gates"], "Independent comparison")
    expected_gates = {component + suffix for component in ("FRZ", "POL", "CT", "interaction")
                      for suffix in ("_agreement", "_resolution_change")}
    need(set(comparison["gates"]) == expected_gates, "Incomplete independent comparison gates")
    for key in ("implementation", "manifest", "pyscf", "coarse_input", "fine_input"):
        verified(comparison[key])
    benchmark = read(comparison["pyscf"]["path"])
    need(benchmark["status"] == "passed_pyscf_benchmark" and benchmark["mode"] == "benchmark",
         "Independent PySCF benchmark did not pass")
    need(benchmark["method"]["xc"] == "PBE", "Wrong independent-test functional")
    list_checks(benchmark["checks"], "Independent PySCF benchmark")
    for key in ("implementation", "almo_source", "fixture"):
        verified(benchmark[key])
    need(core["almo_source"]["sha256"] == benchmark["almo_source"]["sha256"],
         "Core and independent benchmark used different ALMO sources")
    need(core["fixture"]["sha256"] == benchmark["fixture"]["sha256"], "Validation fixture changed")
    receipts = {}
    for name, cutoff in (("coarse", 400), ("fine", 600)):
        native = comparison[name + "_native"]
        output = verified(native["output"])
        receipt_path = output.parent / "RUN_RESULT.json"
        receipt = read(receipt_path)
        need(receipt["status"] == "native_complete_needs_comparison" and receipt["returncode"] == 0
             and not receipt.get("timed_out", False), "Native comparison runner failed")
        need(receipt["cutoff_Ry"] == cutoff and receipt["qos"] == "standard_expand"
             and receipt["native_calls"] == 1, "Wrong native benchmark execution")
        need(receipt["input"]["sha256"] == comparison[name + "_input"]["sha256"],
             "Native actual-input linkage differs")
        need(native["output"] in receipt["artifacts"], "Native output not bound by actual runner")
        for record in receipt["artifacts"]:
            verified(record)
        verified(receipt["input"])
        verified(receipt["runtime"])
        verified(receipt["manifest"])
        receipts[name] = pin(receipt_path)
    return {"core": pin(core_path), "comparison": pin(comparison_path),
            "pyscf_benchmark": comparison["pyscf"], "native_receipts": receipts,
            "almo_source_sha256": core["almo_source"]["sha256"],
            "core_versions": core["versions"], "benchmark_versions": benchmark["versions"]}


def admit_endpoint(directory, system_id, expected, tested_hash):
    directory = Path(directory).resolve()
    run, stage = read(directory / "RESULT.json"), read(directory / "almo_RESULT.json")
    start, info = read(directory / "START.json"), read(directory / "INPUT.json")
    need(run["status"] == "passed" and run["mode"] == "almo" and run["system_id"] == system_id,
         "Unadmitted/wrong quartet endpoint: " + system_id)
    need(start["mode"] == "almo" and start["system_id"] == system_id, "Startup identity differs")
    need(stage["status"] == "passed" and run["almo"] == stage, "ALMO stage/top-level disagreement")
    all_checks(stage["checks"], "ALMO " + system_id)
    need(stage["checks"].get("frozen_saved_state_matches_solver") is True,
         "Saved frozen state was not independently matched to solver F")
    need(abs(dec(stage["saved_frozen_energy_minus_solver_initial_hartree"])) <= ROUNDING_GATE,
         "Saved frozen state differs from solver F")
    need(stage["solver"]["converged"] is True and stage["solver"]["invariants"]["passed"] is True,
         "Solver not converged/invariant-admitted")
    need(dec(stage["solver"]["gradient_norm"]) <= Decimal("1e-7"), "ALMO gradient gate failed")
    need(run["new_scf_calls"] == 0 and stage["new_scf_calls"] == 0, "ALMO endpoint repeated admitted SCF")
    snapshot_inputs(info, expected)
    admit_initialization(run, info)
    need(source_pin(run, "almo.py")["sha256"] == tested_hash, "Endpoint uses untested ALMO source")
    for record in run["source_pins"]:
        verified(record)
    need(start["source_pins"] == run["source_pins"], "Startup/final source pins differ")
    need(start["settings"] == run["settings"] and start["almo_settings"] == run["almo_settings"],
         "Startup/final numerical settings differ")
    need(dec(run["almo_settings"]["gradient_tolerance"]) == Decimal("1e-7"), "Changed ALMO gradient tolerance")
    need(dec(run["almo_settings"]["energy_tolerance"]) == Decimal("1e-10"), "Changed ALMO energy tolerance")
    for label in ("frozen_state", "polarized_state"):
        verified(stage[label])
    need(stage["fixed_grid"] == read(directory / "FIXED_GRID.json"), "Fixed-grid evidence differs")
    need(stage["reused_states"] == read(directory / "REUSED_STATES.json"), "Reuse evidence differs")
    reuse = stage["reused_states"]
    baseline_run, baseline_stage, baseline_info = admit_reused_state(reuse["baseline"], "baseline", system_id)
    snapshot_inputs(baseline_info, expected)
    need(common_energy_settings(baseline_run["settings"]) == common_energy_settings(run["settings"]),
         "Baseline/ALMO electronic model or convergence settings differ")
    need(baseline_run["runtime"]["packages"] == run["runtime"]["packages"], "Baseline/ALMO backend differs")
    comparison = baseline_stage["archive_comparison"]
    need(dec(comparison["full_energy_tolerance_hartree"]) == ARCHIVE_GATE, "Archive gate changed")
    need(comparison["electronic_agreement"] is True and comparison["full_energy_agreement"] is True,
         "Ordinary baseline did not reproduce archive")
    fragment_rows = []
    need(len(reuse["fragments"]) == 2, "Wrong fragment count")
    for index, saved in enumerate(reuse["fragments"]):
        frun, fstage, finfo = admit_reused_state(saved, "fragment")
        need(common_energy_settings(frun["settings"]) == common_energy_settings(run["settings"]),
             "Fragment electronic model or convergence settings differ")
        need(frun["runtime"]["packages"] == run["runtime"]["packages"], "Fragment backend differs")
        expected_fragment = expected["fragments"][index]
        need(finfo["runtime_subset"]["charge"] == expected_fragment["charge"]
             and finfo["runtime_subset"]["spin"] == expected_fragment["spin"]
             and finfo["runtime_subset"]["nelectron"] == expected_fragment["nelectron"],
             "Wrong fragment state")
        need(saved["physics_signature"] == fstage["physics_signature"] == finfo["physics_signature"],
             "Fragment identity signature differs")
        fragment_rows.append(saved)
    e = {key: dec(stage["energies_hartree"][key]) for key in STAGES}
    need(e["U"] == dec(reuse["baseline"]["energy_hartree"]), "U is not the admitted baseline")
    residuals = {"R_minus_fragments": e["R"] - sum(dec(f["energy_hartree"]) for f in fragment_rows),
                 "FRZ_minus_F_minus_R": e["FRZ"] - (e["F"] - e["R"]),
                 "POL_minus_P_minus_F": e["POL"] - (e["P"] - e["F"]),
                 "CT_minus_U_minus_P": e["CT"] - (e["U"] - e["P"]),
                 "component_closure": e["FRZ"] + e["POL"] + e["CT"] - (e["U"] - e["R"]),
                 "full_composite": e["full_PBE0_D4"] - (e["U"] + e["archived_printed_D4"]),
                 "suppressed_composite": e["CT_suppressed_PBE0_D4"] - (e["P"] + e["archived_printed_D4"])}
    need(all(abs(value) <= ROUNDING_GATE for value in residuals.values()), "Endpoint energy identities failed")
    need(e["U"] <= e["P"] + ORDER_GATE and e["P"] <= e["F"] + ORDER_GATE, "Endpoint variational order failed")
    need(e["archived_printed_D4"] == dec(expected["archive_energies"]["d4_printed_hartree"]),
         "Parent D4 scalar differs")
    need(abs(e["U"] - dec(expected["archive_energies"]["electronic_hartree"])) <= ARCHIVE_GATE,
         "Saved U fails literal archive electronic comparison")
    need(abs(e["full_PBE0_D4"] - dec(expected["archive_energies"]["full_pbe0_d4_hartree"])) <= ARCHIVE_GATE,
         "Saved full composite fails literal archive comparison")
    return {"system_id": system_id, "energy": e, "checks": stage["checks"],
            "identity_residuals_hartree": residuals,
            "pins": {filename: pin(directory / filename) for filename in
                     ("RESULT.json", "almo_RESULT.json", "START.json", "INPUT.json", "FIXED_GRID.json", "REUSED_STATES.json")},
            "source_almo": source_pin(run, "almo.py"),
            "settings": run["settings"], "almo_settings": run["almo_settings"],
            "backend_packages": run["runtime"]["packages"],
            "reused_states": reuse,
            "solver_gradient_norm": stage["solver"]["gradient_norm"],
            "archive_comparison": comparison}


def scientific_collection(args, result):
    repo = args.repo.resolve()
    manifest_path = repo / "diagnostics/lady_almo_implementation_20261023/INPUT_PROVENANCE.json"
    manifest = read(manifest_path)
    for record in manifest["input_pins"]:
        verified(record)
    ledger_path = repo / "diagnostics/lady_electronic_decomposition_20261023/cp2k_inputs_v1/ENERGY_LEDGER.json"
    ledger = read(ledger_path)
    need(ledger["status"] == "complete_saved_archived_energy_ledger", "Archive ledger incomplete")
    need(pin(ledger_path) in manifest["input_pins"], "Energy ledger not pinned by input declaration")
    conversion = dec(ledger["units"]["hartree_to_kcal_mol_literal_from_original_collector"])
    prerequisites = admit_prerequisites(args.core_tests, args.independent_comparison)
    rows = {}
    directories = (args.aquo_la, args.aquo_dy, args.mixed_la, args.mixed_dy)
    for name, directory in zip(SYSTEMS, directories):
        rows[name] = admit_endpoint(directory, name, manifest["systems"][name],
                                    prerequisites["almo_source_sha256"])
    first = rows[SYSTEMS[0]]
    for row in rows.values():
        for key in ("settings", "almo_settings", "backend_packages"):
            need(row[key] == first[key], "Quartet does not share " + key)
        for package in ("pyscf", "numpy"):
            need(row["backend_packages"][package] == prerequisites["core_versions"][package]
                 == prerequisites["benchmark_versions"][package], "Quartet and tests use different " + package + " versions")
    # These are exact reused reference states, not an assumed cancellation.
    pairs = ((SYSTEMS[0], SYSTEMS[2], 0), (SYSTEMS[1], SYSTEMS[3], 0),
             (SYSTEMS[0], SYSTEMS[1], 1), (SYSTEMS[2], SYSTEMS[3], 1))
    shared = []
    for left, right, index in pairs:
        a, b = (rows[n]["reused_states"]["fragments"][index] for n in (left, right))
        need(a["physics_signature"] == b["physics_signature"], "Matched fragment identity differs")
        need(a["state"]["sha256"] == b["state"]["sha256"], "Matched fragment was not reused identically")
        need(dec(a["energy_hartree"]) == dec(b["energy_hartree"]), "Shared fragment energy differs")
        shared.append({"left": left, "right": right, "fragment": index,
                       "state_sha256": a["state"]["sha256"], "physics_signature": a["physics_signature"]})
    contrasts = {key: sum(SIGNS[name] * rows[name]["energy"][key] for name in SYSTEMS)
                 for key in STAGES}
    closure = {"reference_R_cancellation": contrasts["R"],
               "FRZ_POL_CT_minus_U_minus_R": contrasts["FRZ"] + contrasts["POL"] + contrasts["CT"] - (contrasts["U"] - contrasts["R"]),
               "full_minus_suppressed_minus_CT": contrasts["full_PBE0_D4"] - contrasts["CT_suppressed_PBE0_D4"] - contrasts["CT"],
               "full_minus_U_minus_D4": contrasts["full_PBE0_D4"] - contrasts["U"] - contrasts["archived_printed_D4"],
               "suppressed_minus_P_minus_D4": contrasts["CT_suppressed_PBE0_D4"] - contrasts["P"] - contrasts["archived_printed_D4"]}
    need(all(abs(value) <= 4 * ROUNDING_GATE for value in closure.values()), "Quartet closure or R cancellation failed")
    old = ledger["double_difference"]
    model = dec(old["original_amoeba_transfer_value_kcal_mol"])
    original = {"electronic_kcal_mol": dec(old["electronic_kcal_mol_from_tokens"]),
                "printed_D4_kcal_mol": dec(old["direct_printed_d4_kcal_mol"]),
                "effective_D4_token_difference_kcal_mol": dec(old["effective_d4_kcal_mol_from_full_minus_electronic"]),
                "full_token_kcal_mol": dec(old["full_kcal_mol_from_tokens"]),
                "published_full_kcal_mol": dec(old["published_reference_value_kcal_mol"]),
                "AMOEBA_transfer_kcal_mol": model}
    comparisons = {"new_full_minus_archived_full_kcal_mol": contrasts["full_PBE0_D4"] * conversion - original["full_token_kcal_mol"],
                   "new_electronic_minus_archived_electronic_kcal_mol": contrasts["U"] * conversion - original["electronic_kcal_mol"],
                   "AMOEBA_minus_new_full_kcal_mol": model - contrasts["full_PBE0_D4"] * conversion,
                   "AMOEBA_minus_CT_suppressed_kcal_mol": model - contrasts["CT_suppressed_PBE0_D4"] * conversion,
                   "AMOEBA_minus_published_archive_kcal_mol": model - original["published_full_kcal_mol"],
                   "full_minus_CT_suppressed_kcal_mol": contrasts["CT"] * conversion}
    result.update(status="passed_saved_quartet_collection", prerequisites=prerequisites,
                  input_provenance=pin(manifest_path), archived_energy_ledger=pin(ledger_path),
                  hartree_to_kcal_mol=conversion, sign_convention="(Dy-La)_mixed - (Dy-La)_aquo8; negative favors Dy in this fixed-geometry substitution contrast",
                  coefficients=SIGNS, endpoints=rows,
                  double_differences={key: {"hartree": value, "kcal_mol": value * conversion}
                                      for key, value in contrasts.items()},
                  shared_reference_states=shared, closure_residuals_hartree=closure,
                  archived_reference=original, comparisons=comparisons,
                  numerical_limits={"per_endpoint_archive_gate_hartree": ARCHIVE_GATE,
                                    "arithmetic_endpoint_tolerance_hartree": ROUNDING_GATE,
                                    "arithmetic_quartet_tolerance_hartree": 4 * ROUNDING_GATE,
                                    "worst_case_four_archive_gates_kcal_mol": 4 * ARCHIVE_GATE * conversion,
                                    "estimated_CT_uncertainty": None,
                                    "materiality_threshold": None},
                  interpretation_limits=[
                      "Fixed archived geometries and frozen-f effective singlets; no equilibrium or protein binding free energy",
                      "AO-block constrained stationary solutions, not certified global minima or experimentally unique energy components",
                      "Core tests and independent PBE water comparison validate sampled implementation behavior, not physical Dy multiplets",
                      "Archive compatibility tolerances and closure residuals are not confidence intervals for CT",
                      "The same archived parent D4 scalar enters full and constrained totals once; no isolated-fragment D4 is invented",
                      "CT contribution/sign is reported directly without choosing a post hoc materiality cutoff or assigning force-field causality"])


def report_markdown(result):
    if result["status"] != "passed_saved_quartet_collection":
        return "# Quartet collection failed\n\nNo scientific contrast was admitted.\n\n" + result.get("error", "Unknown admission failure") + "\n"
    rows = result["double_differences"]
    lines = ["# La/Dy ALMO quartet", "", "All supplied implementation, independent comparison, baseline and quartet admission gates passed.", "",
             "Sign: `(Dy−La)mixed − (Dy−La)aquo8`; negative values favor Dy within this fixed-geometry substitution comparison.", "",
             "| Quantity | kcal/mol |", "|---|---:|"]
    for key in STAGES:
        lines.append(f"| {key} | {rows[key]['kcal_mol']:+.9f} |")
    lines += ["", "## Comparison with existing results", "",
              "| Quantity | kcal/mol |", "|---|---:|"]
    for key, value in result["comparisons"].items():
        lines.append(f"| {key} | {value:+.9f} |")
    lines += ["", f"Archived full PBE0-D4 contrast: {result['archived_reference']['published_full_kcal_mol']:+.9f} kcal/mol. Transferred AMOEBA contrast: {result['archived_reference']['AMOEBA_transfer_kcal_mol']:+.9f} kcal/mol.", "",
              "`CT = U−P` is the change when the metal–shell constraint is released. Polarization remains allowed inside each fragment. The whole shell is one fragment.", "",
              "## Numerical scope", "",
              f"Reference cancellation residual: {result['closure_residuals_hartree']['reference_R_cancellation']:.3E} Eh. Largest reported quartet identity residual: {max(abs(v) for v in result['closure_residuals_hartree'].values()):.3E} Eh.", "",
              "The four per-endpoint archive gates permit a worst-case sum of 0.251003790 kcal/mol; this is a compatibility bound, not a CT error bar. No CT uncertainty or materiality threshold has been inferred.", "",
              "Results concern the exact frozen geometries, effective singlet Hamiltonians and AO-block partition. They do not establish physical Dy-state accuracy, solution/protein selectivity, a global constrained minimum or a unique force-field error. See RESULTS.json for every artifact pin, raw endpoint scalar and admission residual.", ""]
    return "\n".join(lines)


def json_default(value):
    if isinstance(value, Decimal):
        return str(value)
    raise TypeError(type(value).__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(os.environ.get("LADY_REPO_ROOT", Path(__file__).resolve().parent.parent.parent)))
    for flag in ("aquo-la", "aquo-dy", "mixed-la", "mixed-dy", "core-tests", "independent-comparison", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    need(os.environ.get("SLURM_JOB_ID"), "Saved numerical collection requires Slurm")
    need(os.environ.get("SLURM_JOB_QOS") == "standard_expand", "Require standard_expand")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    result = {"status": "failed_saved_quartet_collection", "collector": pin(__file__),
              "job_id": os.environ["SLURM_JOB_ID"], "new_molecular_calls": 0,
              "new_fragment_calls": 0, "new_solver_calls": 0,
              "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    result["requested_inputs"] = {key: str(value.resolve()) for key, value in vars(args).items()}
    try:
        scientific_collection(args, result)
    except Exception as error:
        result["error"] = str(error)
        result["traceback"] = traceback.format_exc()
    result["elapsed_seconds"] = time.monotonic() - started
    result["allocated_cpu_seconds"] = int(os.environ.get("SLURM_CPUS_PER_TASK", "1")) * result["elapsed_seconds"]
    usage = resource.getrusage(resource.RUSAGE_SELF)
    result["process_cpu_seconds"] = usage.ru_utime + usage.ru_stime
    result["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    (output / "REPORT.md").write_text(report_markdown(result))
    with (output / "RESULTS.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False, default=json_default)
        stream.write("\n")
    return 0 if result["status"] == "passed_saved_quartet_collection" else 1


if __name__ == "__main__":
    raise SystemExit(main())
