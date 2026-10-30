#!/usr/bin/env python3
"""Collect exactly four terminal saved ALMO receipts; never run a parser or CP2K.

  python COLLECT.py --results-root /absolute/cp2k_results_v1 \
      --config /absolute/cp2k_inputs_v1/CONFIG.json --output /fresh/collection

The results root may instead be the configured workspace tasks directory.
Writes RESULTS.json, table.csv and REPORT.md exclusively in a fresh directory.
Exit 0 means all four original output checks and lineage checks passed. Missing,
failed or changed records produce a failure report without a quartet contrast.
All arithmetic uses saved Decimal tokens. No package/native code is imported.
"""

import argparse
import csv
import datetime
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import sys
import time


LEDGER_SHA = "24c699cc97c3c55857c56a5ee07dda6213eeaf12e59af69395f23bbd05b732c5"
PARSER_SHA = "743ddda796bb469ad21eb5329e0b7df5f02ebd29e34f319b11419f4f3cdf6202"
SIGNS = {"aquo8_central_La": 1, "aquo8_central_Dy": -1,
         "hydrated_acetate_mono_central_La": -1,
         "hydrated_acetate_mono_central_Dy": 1}
HA = Decimal("627.509474")
FOUR_CELL_BOUND = Decimal("4e-4")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(value):
    result = Decimal(str(value).replace("D", "E").replace("d", "e"))
    require(result.is_finite(), "Nonfinite saved scalar")
    return result


class Saved:
    def __init__(self):
        self.pins = {}

    def pin(self, path, expected=None):
        path = Path(path).resolve()
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024*1024), b""):
                digest.update(block)
        result = {"path": str(path), "sha256": digest.hexdigest()}
        require(expected is None or result["sha256"] == expected, f"Hash mismatch: {path}")
        old = self.pins.get(str(path))
        require(old is None or old == result, f"Input changed during collection: {path}")
        self.pins[str(path)] = result
        return result

    def read(self, path, expected=None):
        p = self.pin(path, expected)
        raw = Path(p["path"]).read_bytes()
        require(hashlib.sha256(raw).hexdigest() == p["sha256"], f"Input changed while reading: {path}")
        return json.loads(raw), p


def raw_records(value):
    if isinstance(value, dict):
        if "line" in value and "raw" in value:
            yield value
        for child in value.values():
            yield from raw_records(child)
    elif isinstance(value, list):
        for child in value:
            yield from raw_records(child)


def required_gates():
    names = {"native_exit_code_zero", "native_no_failure_or_early_stop", "native_stage_order",
             "native_full_scf", "two_domains", "block_native_settings", "full_native_settings",
             "no_outer_iteration_extension", "parent_AO_count_from_domains", "fragment_reference_sum",
             "full_trace_matches_final_energy", "component_arithmetic_closure",
             "variational_block_not_above_frozen", "variational_full_not_above_block",
             "archived_electronic_agreement", "archived_total_agreement"}
    for phase in ("block_diis", "full_pcg"):
        names.update({phase+"_iteration_sequence", phase+"_printed_residual"})
    for i in (1, 2):
        names.update({f"domain_{i}_state", f"fragment_{i}_state", f"fragment_{i}_membership",
                      f"fragment_{i}_native_convergence"})
        names.update(f"fragment_{i}_{label}_if_emitted" for label in
                     ("electrons", "occupied orbitals", "orbital functions"))
    return names


def cell(saved, root, config, config_pin, ledger, gid, parser_pin):
    spec = config["systems"][gid]
    task = Path(spec["task"]).resolve()
    result, result_pin = saved.read(root/gid/"RESULT.json")
    run, run_pin = saved.read(root/gid/"RUN_RESULT.json")
    # Diagnostic copies must be byte-identical to their original terminal receipts.
    saved.pin(task/"RESULT.json", result_pin["sha256"])
    saved.pin(task/"RUN_RESULT.json", run_pin["sha256"])
    require(run["status"] == "native_finished_review_required" and run.get("error") is None,
            f"{gid}: native runner not successfully terminal")
    require(run["geometry_id"] == result["geometry_id"] == gid, f"{gid}: geometry label mismatch")
    require(run["config_sha256"] == config_pin["sha256"], f"{gid}: different configuration")
    require(run["molecular_parent_calls"] == 1, f"{gid}: native call count differs")
    require(Path(run["parser_result"]).resolve() == task/"RESULT.json", f"{gid}: parser output path differs")
    require(result["implementation"] == parser_pin, f"{gid}: parser source identity differs")
    require(result["ledger"]["sha256"] == LEDGER_SHA, f"{gid}: energy ledger differs")
    require(result["passed"] is True and result["status"] == "passed_cp2k_almo_output_checks"
            and result["error"] is None and result["exit_code"] == 0, f"{gid}: output admission failed")
    gates = result["gates"]
    require(required_gates() <= gates.keys(), f"{gid}: missing mandatory parser gates")
    require(all(g.get("passed") is True for g in gates.values()), f"{gid}: nonpassing parser gate")
    require(result["archive"] == ledger["systems"][gid]["provenance"], f"{gid}: archive lineage differs")
    require(spec["geometry_state"] == ledger["systems"][gid]["geometry_state"], f"{gid}: electronic state differs")
    require(int(run["allocation"]["SLURM_CPUS_PER_TASK"]) == config["cpus"] == 24,
            f"{gid}: allocation differs")
    commands = run["commands"]
    require(len(commands) == 3 and all(c["returncode"] == 0 and not c.get("timed_out", False)
                                      for c in commands), f"{gid}: incomplete or failed invocation chain")
    require(commands[1]["cap_seconds"] == config["native_cap_seconds"] and
            0 < number(config["native_cap_seconds"]) <= 7200,
            f"{gid}: native process cap differs")
    require(number(commands[1]["elapsed_seconds"]) <= number(config["native_cap_seconds"]),
            f"{gid}: native process exceeded its configured cap")
    require(commands[1]["argv"][-5:] == ["/work/build/bin/cp2k.psmp", "-i", "/task/input.inp", "-o", "/task/endpoint.out"],
            f"{gid}: unexpected native command endpoint")
    require(str(config["parser"]) in commands[2]["argv"], f"{gid}: parser command differs")
    if gid != config["pilot"]:
        admission_pin = run["pilot_admission"]
        require(Path(admission_pin["path"]).resolve() == Path(config["pilot_admission"]).resolve(),
                f"{gid}: pilot admission path differs")
        admission, _ = saved.read(admission_pin["path"], admission_pin["sha256"])
        require(admission["admitted"] is True and admission["config_sha256"] == config_pin["sha256"],
                f"{gid}: pilot was not admitted for this configuration")
    artifacts = {Path(a["path"]).name: a for a in run["artifacts"]}
    require({"input.inp", "endpoint.out", "parse.log", "native_launcher.log", "syntax.log"} <= artifacts.keys(),
            f"{gid}: incomplete executed-artifact inventory")
    for name, record in artifacts.items():
        require(Path(record["path"]).resolve() == task/"execution"/name, f"{gid}: artifact escaped its task")
        saved.pin(record["path"], record["sha256"])
        if root/gid != task:
            saved.pin(root/gid/name, record["sha256"])
    require(artifacts["input.inp"]["sha256"] == spec["deck_sha256"], f"{gid}: executed deck differs")
    output_pin = artifacts["endpoint.out"]
    require(Path(result["output"]["path"]).resolve() == task/"execution/endpoint.out"
            and result["output"]["sha256"] == output_pin["sha256"], f"{gid}: parser raw output identity differs")
    lines = Path(output_pin["path"]).read_text().splitlines()
    for record in raw_records(result["evidence"]):
        line = record["line"]
        require(isinstance(line, int) and 1 <= line <= len(lines) and lines[line-1] == record["raw"],
                f"{gid}: saved native line evidence differs")
    ev = result["evidence"]
    tokens = {"R": ev["R"][0]["tokens"][0], "F": ev["F"][0]["tokens"][0],
              "P": ev["block_diis"][-1]["tokens"][2], "U": ev["U"][0]["tokens"][0]}
    vals = {key: number(value) for key, value in tokens.items()}
    require({k: number(v) for k, v in result["energies_hartree"].items()} == vals,
            f"{gid}: reported energies differ from saved native tokens")
    vals.update(FRZ=vals["F"]-vals["R"], POL=vals["P"]-vals["F"], CT=vals["U"]-vals["P"])
    require(all(number(result["components_hartree"][k]) == vals[k] for k in ("FRZ", "POL", "CT")),
            f"{gid}: component bookkeeping differs")
    vals["metal_reference"] = number(ev["fragment_1_energy"][0]["tokens"][0])
    vals["shell_reference"] = number(ev["fragment_2_energy"][0]["tokens"][0])
    a = ledger["systems"][gid]["energy_hartree"]
    vals["D4"] = number(a["d4_printed"]["token"])
    require(number(result["dispersion_bookkeeping"]["archived_D4_hartree"]) == vals["D4"], f"{gid}: D4 differs")
    vals["full_total"] = vals["U"]+vals["D4"]
    vals["suppressed_total"] = vals["P"]+vals["D4"]
    vals["interaction_electronic"] = vals["U"]-vals["R"]
    vals["archive_electronic"] = number(a["electronic_after_final_cosx_integration"]["token"])
    vals["archive_total"] = number(a["full_pbe0_d4"]["token"])
    vals["full_minus_archive"] = vals["full_total"]-vals["archive_total"]
    elapsed = number(run["elapsed_seconds"])
    require(abs(number(run["allocated_cpu_seconds_worker_elapsed"])-elapsed*24) <= Decimal("1e-8"),
            f"{gid}: allocated-capacity arithmetic inconsistent")
    return {"result": result_pin, "run_result": run_pin, "job_id": run["job_id"],
            "gates": gates, "energy_tokens": tokens, "energies_hartree": {k: str(v) for k, v in vals.items()},
            "fragment_reference_sum_residual_hartree": str(vals["metal_reference"]+vals["shell_reference"]-vals["R"]),
            "convergence_evidence": {key: ev[key] for key in
                ("fragment_1_converged", "fragment_2_converged", "fragment_1_last_scf_row", "fragment_2_last_scf_row")},
            "block_final": ev["block_diis"][-1], "full_final": ev["full_pcg"][-1],
            "native_elapsed_seconds": commands[1]["elapsed_seconds"], "worker_elapsed_seconds": run["elapsed_seconds"],
            "allocated_cpu_seconds_worker_elapsed": run["allocated_cpu_seconds_worker_elapsed"]}


def collect(args):
    started = time.monotonic()
    cpu = time.process_time()
    saved = Saved()
    result = {"status": "incomplete_or_failed_quartet", "passed": False, "systems": {}, "errors": [],
              "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "implementation": saved.pin(__file__), "double_differences": None,
              "new_molecular_calls": 0, "parser_reexecutions": 0}
    try:
        # The explicit configuration is bound to every actual terminal RUN_RESULT
        # by its SHA256 below. No directory/version or initial-attempt hash is
        # hardcoded; this permits a separately declared resource-only recovery.
        config, cp = saved.read(args.config)
        result["configuration"] = cp
        pins = {str(Path(p["path"]).resolve()): p for p in config["pins"]}
        parser_path = Path(config["parser"]).resolve()
        parser_pin = pins[str(parser_path)]
        require(parser_pin["sha256"] == PARSER_SHA, "Unexpected frozen parser")
        ledger_path = parser_path.with_name("ENERGY_LEDGER.json")
        ledger_pin = pins[str(ledger_path)]
        require(ledger_pin["sha256"] == LEDGER_SHA, "Configuration energy ledger pin differs")
        ledger, lp = saved.read(ledger_path, LEDGER_SHA)
        result["energy_ledger"] = lp
        require(set(config["systems"]) == set(SIGNS) and ledger["double_difference"]["coefficients"] == SIGNS,
                "Quartet or sign convention differs")
        require(config["mpi_ranks"] == 12 and config["omp_threads"] == 2, "MPI/OpenMP configuration differs")
        # These are exact configured source/assets/binary/image pins; streaming
        # hashing does not load native libraries or execute molecular software.
        for p in config["pins"]:
            saved.pin(p["path"], p["sha256"])
        root = args.results_root.resolve()
        taskroots = {Path(s["task"]).resolve().parent for s in config["systems"].values()}
        require(root == Path(config["diagnostic_results"]).resolve() or root in taskroots,
                "Results root is outside immutable configuration")
        result["results_root"] = str(root)
        for gid in SIGNS:
            try:
                result["systems"][gid] = {"status": "passed_saved_receipt_joins",
                                          **cell(saved, root, config, cp, ledger, gid, parser_pin)}
            except (OSError, ValueError, KeyError, TypeError, IndexError, ArithmeticError) as error:
                result["systems"][gid] = {"status": "missing_or_failed", "error": f"{type(error).__name__}: {error}"}
                result["errors"].append({"geometry_id": gid, "error": result["systems"][gid]["error"]})
        if result["errors"]:
            raise ValueError("At least one exact quartet slot is missing, failed or inconsistent; no quartet inference")
        require(len({r["job_id"] for r in result["systems"].values()}) == 4, "Repeated parent job identity")
        rows = {g: {k: number(v) for k, v in r["energies_hartree"].items()} for g, r in result["systems"].items()}
        def dd(key):
            return sum((SIGNS[g]*rows[g][key] for g in SIGNS), Decimal(0))
        contrasts = {key: dd(key) for key in next(iter(rows.values()))}
        result["double_differences"] = {"definition": "(Dy mixed - La mixed) - (Dy aquo8 - La aquo8)",
                                        "signs": SIGNS, "hartree": {k: str(v) for k, v in contrasts.items()},
                                        "kcal_mol": {k: str(v*HA) for k, v in contrasts.items()}}
        aLa, aDy, mLa, mDy = (rows[g] for g in SIGNS)
        residuals = {
            "La_ion_mixed_minus_aquo": mLa["metal_reference"]-aLa["metal_reference"],
            "Dy_ion_mixed_minus_aquo": mDy["metal_reference"]-aDy["metal_reference"],
            "aquo_shell_Dy_minus_La": aDy["shell_reference"]-aLa["shell_reference"],
            "mixed_shell_Dy_minus_La": mDy["shell_reference"]-mLa["shell_reference"],
            "DD_R_reported": contrasts["R"],
            "DD_sum_actual_fragments": contrasts["metal_reference"]+contrasts["shell_reference"],
            "DD_R_minus_actual_fragment_sum": contrasts["R"]-contrasts["metal_reference"]-contrasts["shell_reference"],
            "full_component_closure": contrasts["full_total"]-sum(contrasts[k] for k in ("FRZ", "POL", "CT", "R", "D4")),
            "suppressed_component_closure": contrasts["suppressed_total"]-sum(contrasts[k] for k in ("FRZ", "POL", "R", "D4")),
            "full_minus_suppressed_minus_CT": contrasts["full_total"]-contrasts["suppressed_total"]-contrasts["CT"]}
        result["reference_and_closure_residuals_hartree"] = {k: str(v) for k, v in residuals.items()}
        original = ledger["double_difference"]
        model = number(original["original_amoeba_transfer_value_kcal_mol"])
        oldref = number(original["published_reference_value_kcal_mol"])
        result["comparison_kcal_mol"] = {
            "archived_published_full_DD": str(oldref), "archived_decimal_full_DD": str(contrasts["archive_total"]*HA),
            "native_full_DD": str(contrasts["full_total"]*HA), "native_CT_suppressed_DD": str(contrasts["suppressed_total"]*HA),
            "native_minus_archived_decimal_full_DD": str(contrasts["full_minus_archive"]*HA),
            "native_minus_archived_published_full_DD": str(contrasts["full_total"]*HA-oldref),
            "original_AMOEBA_transfer_DD": str(model), "model_minus_native_full_DD": str(model-contrasts["full_total"]*HA),
            "model_minus_CT_suppressed_DD": str(model-contrasts["suppressed_total"]*HA),
            "original_model_minus_archived_published_DD": str(model-oldref),
            "CT_DD": str(contrasts["CT"]*HA)}
        result["numerical_resolution"] = {
            "per_cell_archive_gate_hartree": "1e-4", "four_cell_worst_case_hartree": str(FOUR_CELL_BOUND),
            "four_cell_worst_case_kcal_mol": str(FOUR_CELL_BOUND*HA),
            "observed_full_DD_disagreement_within_four_cell_bound": abs(contrasts["full_minus_archive"]) <= FOUR_CELL_BOUND,
            "observed_absolute_CT_DD_above_four_cell_comparison_scale": abs(contrasts["CT"]) > FOUR_CELL_BOUND,
            "interpretation": "4e-4 Eh is the triangle-inequality bound implied by four chosen archive-agreement gates. It is NOT a universal CT/component uncertainty, SCF error bar, or physical accuracy bound. Repeated-reference residuals are observed and retained, never forced to zero. No new reference-cancellation threshold is imposed."}
        result["archive_cautions"] = {g: {"scf": ledger["systems"][g]["scf"],
                                          "final_cosx_reintegration": ledger["systems"][g]["energy_hartree"]["final_cosx_reintegration_change"]}
                                      for g in SIGNS}
        result["limitations"] = [
            "Only the combined frozen FRZ term is available; electrostatic/Pauli separation is unavailable.",
            "True absolute D4-corrected isolated-fragment interactions are unavailable without isolated-shell D4 references; no missing constants are invented. Matched reference D4 cancels algebraically in this quartet DD.",
            "The raw electronic interaction DD is DD[U]-DD[R]; full absolute-energy DD includes the observed DD[R] when reconstructed from FRZ/POL/CT. These distinct quantities are not silently equated.",
            "Density/SCF flags and COSX integration shifts in the archive are retained; decimal printing is not their uncertainty. Native ALMO residuals are separate from telescoping closure.",
            "Own-fragment-basis AO-block ALMO at fixed geometries/effective singlets is a conditional decomposition, not a unique mechanism, a proof of missing force-field CT, an explicit-f validation, or a physical binding free energy.",
            "Root must inspect actual input/output/convergence and decide whether any repeated-reference discrepancy contaminates interpretation. Passing this saved collector does not itself release further calculations."]
        result["cost"] = {"native_parent_calls_in_these_four_successful_receipts": 4, "observed_internal_fragment_SCFs": 8,
                          "observed_block_optimizations": 4, "observed_full_optimizations": 4,
                          "native_wall_seconds_sum": str(sum(number(r["native_elapsed_seconds"]) for r in result["systems"].values())),
                          "worker_elapsed_allocated_cpu_seconds_sum": str(sum(number(r["allocated_cpu_seconds_worker_elapsed"]) for r in result["systems"].values())),
                          "job_ids": [r["job_id"] for r in result["systems"].values()],
                          "configured_native_cap_seconds": config["native_cap_seconds"],
                          "scope": "These four successful receipts only. Prior failed startup/OOM attempts, builds, scheduler parent accounting and peak RSS must be retained in the root cumulative cost ledger; they are not zero or erased by this collector. Runner wall and 24*worker wall are elapsed time/allocated capacity, not measured CPU utilization. This collector does not independently admit a resource-only recovery or its cumulative budget."}
        # Receipt mutation cannot turn a previously read running job into a
        # silently mixed terminal collection. Native output hashes were checked.
        for row in result["systems"].values():
            saved.pin(row["result"]["path"], row["result"]["sha256"])
            saved.pin(row["run_result"]["path"], row["run_result"]["sha256"])
        result["status"] = "complete_saved_quartet_comparison"
        result["passed"] = True
    except (OSError, ValueError, KeyError, TypeError, IndexError, ArithmeticError) as error:
        result["errors"].append({"error": f"{type(error).__name__}: {error}"})
        result["double_differences"] = None
        for key in ("comparison_kcal_mol", "reference_and_closure_residuals_hartree", "numerical_resolution"):
            result.pop(key, None)
    result["verified_pins"] = list(saved.pins.values())
    result["postprocessing_cost"] = {"wall_seconds": time.monotonic()-started, "process_cpu_seconds": time.process_time()-cpu}
    return result


def write_outputs(destination, result):
    destination.mkdir(parents=True, exist_ok=False)
    (destination/"RESULTS.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    keys = ["R", "F", "P", "U", "D4", "FRZ", "POL", "CT", "full_total", "suppressed_total",
            "metal_reference", "shell_reference", "full_minus_archive"]
    with (destination/"table.csv").open("x", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["geometry_id", "status"]+[k+"_hartree" for k in keys])
        for gid, row in result["systems"].items():
            writer.writerow([gid, row["status"]]+[row.get("energies_hartree", {}).get(k, "") for k in keys])
        if result["passed"]:
            writer.writerow(["DD", "complete"]+[result["double_differences"]["hartree"][k] for k in keys])
    text = ["# Saved ALMO quartet comparison", "", f"Status: **{result['status']}**.", ""]
    if not result["passed"]:
        text += ["At least one required record or join failed. No quartet contrast is admitted.", ""]
        text += [f"- {e.get('geometry_id', 'collection')}: {e['error']}" for e in result["errors"]]
    else:
        cmp = result["comparison_kcal_mol"]
        text += ["All values below are fixed-geometry kcal/mol; DD means (Dy−La)mixed − (Dy−La)aquo.", "",
                 "| Quantity | DD (kcal/mol) |", "|---|---:|"]
        for label, key in [("Native full PBE0 + archived D4", "native_full_DD"),
                           ("CT-suppressed P + archived D4", "native_CT_suppressed_DD"),
                           ("CT contribution U−P", "CT_DD"),
                           ("Archived published PBE0-D4", "archived_published_full_DD"),
                           ("Native − archived full", "native_minus_archived_decimal_full_DD"),
                           ("Original AMOEBA transfer", "original_AMOEBA_transfer_DD"),
                           ("AMOEBA − native full", "model_minus_native_full_DD")]:
            text.append(f"| {label} | {cmp[key]} |")
        text += ["", "| Electronic/dispersion contribution | DD (kcal/mol) |", "|---|---:|"]
        text += [f"| {key} | {result['double_differences']['kcal_mol'][key]} |" for key in ("FRZ", "POL", "CT", "R", "D4")]
        text += ["", "Full DD = FRZ + POL + CT + R + D4. The observed reference contribution R is retained.", "",
                 "Repeated-reference residuals (Eh):", ""]
        text += [f"- {key}: {value}" for key, value in result["reference_and_closure_residuals_hartree"].items()]
        text += ["", "The four-cell archive-agreement worst-case scale is 0.0004 Eh = 0.2510037896 kcal/mol. "
                 "It is not a universal CT uncertainty. Actual SCF residuals, raw printing, arithmetic closure and "
                 "archive differences remain separate in RESULTS.json.", "",
                 "Only combined FRZ is available. Absolute isolated-shell D4 reference energies are missing; their "
                 "matched quartet DD cancels, while absolute corrected interactions remain unavailable.", "",
                 "Interpretation: this compares release of the declared metal–whole-shell orbital constraint. "
                 "Even a substantial CT term does not identify a force-field defect or establish a physical binding free energy. "
                 "Root review of raw convergence and repeated references remains necessary; no further run is released.", "",
                 "Cost for these four successful native parents: wall sum "+result["cost"]["native_wall_seconds_sum"]+" s; worker elapsed allocated capacity "+
                 result["cost"]["worker_elapsed_allocated_cpu_seconds_sum"]+" CPU-s. Prior failed/OOM attempts and scheduler/peak-RSS accounting remain separate, nonzero consumed work; see the root cumulative ledger."]
    (destination/"REPORT.md").write_text("\n".join(text)+"\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-root", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    if args.output.exists():
        ap.error("Output directory exists; preserve prior collection and choose a fresh directory")
    with localcontext() as context:
        context.prec = 50
        result = collect(args)
    write_outputs(args.output, result)
    print(json.dumps({"status": result["status"], "output": str(args.output.resolve())}))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
