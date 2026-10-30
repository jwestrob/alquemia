#!/usr/bin/env python3
"""Saved-only, exact-artifact correction of one ORCA backend parser predicate.

No executable is launched. Original results/source remain immutable. The
existing parser is reused; only its direct auxiliary-metric veto is corrected
after independent native SCF-timing and Hamiltonian checks pass.
"""
import argparse
from decimal import Decimal
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import time
import traceback

SOURCE_SHA = "eb1cf5247f98e982282497dfcddccac717be4080c45a039147898619411646f9"
MANIFEST_SHA = "2bc29798a5ace80a91909ee53171b26ea68b600d83b6bf1dba7ac472d4afc8ba"
PINS = {
    "direct": {
        "RESULT.json": "300dc77735f95f8ccb7b1cf973414936528f4821f14ca5ade88105104fb28f87",
        "START.json": "3545522de607de038567f729a628c68dd35a1152645ad42d42cd480ed6f8f58c",
        "frozen_direct.out": "ead114d6e53fea6aad3b77da1acc5ef5acb90136a3e0b2fe8c31d8f2868c1566",
        "frozen_direct.err": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        "frozen_direct.property.txt": "1376ac19a6f5c0d634f5f9d1b02901b6c3f4662e4a93091356480135a5307da0",
        "frozen_direct.bibtex": "cd0b51f29484fc04741d306638c62fc5925f21810b86692297179c47127d666f",
    },
    "rijcosx": {
        "RESULT.json": "fbaba70a2e3a029f0b2d5d670e00b18964d36f67be668107c48d9bba4a8297aa",
        "START.json": "32116b5478f4c164176e412c9417278c85001e18a661573202a3ea33dabeb0a2",
        "frozen_rijcosx.out": "0df28a8bb73c05486227d5f47d8a2d59cec556f5a2581a7ef7b2d8cdace2f61a",
        "frozen_rijcosx.err": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        "frozen_rijcosx.property.txt": "ffc81b3de81e5fbb59631aa3f106818f0a068eb00811d48636f7eb5e0325d763",
        "frozen_rijcosx.bibtex": "2d4d762bab2ba5a8843920bf1d8318fcb5d7110fc9781d442ced4ef3d370b7a3",
    },
}


def backend_evidence(m, variant, text, property_text, bibtex):
    """Require independent native timing channels, exact flags and log context."""
    blocks = re.findall(r"(?ms)^\$SCF_Timings\s*\n(.*?)^\$End\s*$", property_text)
    m.require(len(blocks) == 1, "Expected one native SCF timing block")
    fields = re.findall(r'(?m)^\s*&([A-Za-z0-9_]+)\s+\[&Type "Double"\]\s+(' + m.NUM + ')', blocks[0])
    timings = dict(fields)
    m.require(len(fields) == len(timings), "Duplicate native timing channel")
    hamiltonians = re.findall(r"(?s)SCF SETTINGS\s*\n-+\s*\nHamiltonian:\s*\n(.*?)\nGeneral Settings:", text)
    m.require(len(hamiltonians) == 1, "Expected one complete native Hamiltonian block")
    hamiltonian = hamiltonians[0]
    ri_on = "RI-approximation to the Coulomb term is turned on" in hamiltonian
    cosx_on = bool(re.search(r"RIJ-COSX.*\.{2,}\s+on", hamiltonian))
    citations = ("10.1002/jcc.10318", "10.1016/j.chemphys.2008.10.036", "10.1063/5.0058766")
    metric = "Calculating RI/J V-Matrix + Cholesky decomp"
    checks = {
        "metric_prep_occurs_only_in_startup": text.count(metric) == 1 and
            text.index("SHARK INTEGRAL PACKAGE") < text.index(metric) < text.index("SCF SETTINGS"),
        "auxiliary_basis_remains_supplied": bool(re.search(r"# of basis functions in Aux-J\s+\.{2,}\s+919", text)),
        "explicit_noiter_contract": all(text.count(p) == 1 for p in m.NOITER_PHRASES),
    }
    if variant == "direct":
        checks.update(
            conventional_J_timing_present="J" in timings and m.decimal(timings["J"]) > 0,
            no_RI_or_COSX_timing=not (set(timings) & {"SPLITRIJ", "RIJ", "RIJK", "COSX"}),
            no_active_RI_or_COSX_hamiltonian=not ri_on and not cosx_on and
                not re.search(r"RIJ-ON|RI-JK|RIJK", hamiltonian),
            no_COSX_grid="COSX GRID GENERATION" not in text,
            no_RI_or_COSX_method_citations=all(doi not in bibtex for doi in citations),
        )
    else:
        checks.update(
            RI_and_COSX_timing_present=all(k in timings and m.decimal(timings[k]) > 0 for k in ("SPLITRIJ", "COSX")),
            no_conventional_J_timing="J" not in timings,
            active_RI_and_COSX_hamiltonian=ri_on and cosx_on,
            COSX_grid_present="COSX GRID GENERATION" in text,
            RI_and_COSX_method_citations=all(doi in bibtex for doi in citations),
        )
    m.require(all(checks.values()), "Independent backend evidence failed: " + repr(checks))
    return {"checks": checks, "SCF_timing_tokens_seconds": timings,
            "native_hamiltonian_block": hamiltonian, "native_timing_block": blocks[0],
            "metric_prep_line": next({"line": i, "text": line} for i, line in enumerate(text.splitlines(), 1) if metric in line),
            "interpretation": "The startup auxiliary metric is distinct from the density-dependent Fock backend; actual SCF timing channels corroborate the requested backend."}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--rijcosx", type=Path, required=True)
    p.add_argument("--direct", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    manifest_path = args.manifest.resolve()
    source = manifest_path.parent / "source/orca_frozen_diagnosis.py"
    # Importing the pinned stdlib adapter defines functions; it launches nothing.
    import hashlib
    if hashlib.sha256(source.read_bytes()).hexdigest() != SOURCE_SHA:
        raise ValueError("Original parser source changed")
    spec = importlib.util.spec_from_file_location("original_frozen_parser", source)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    allocation = m.slurm_allocation()
    m.require(allocation["allocated_cpus"] == 1, "Saved-only recovery requires one CPU")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    result = {"schema": "orca_frozen_backend_saved_recovery_v1", "status": "failed",
              "allocation": allocation, "created_utc": m.stamp(), "new_native_calls": 0,
              "source": m.pin(__file__), "original_parser": m.pin(source),
              "baseline_admitted": False, "scf_optimization_performed": False,
              "interpretation": "Fixed-density one-Fock pair, not converged SCF or ALMO admission.",
              "documentation": {
                  "metric_and_RI_off": "https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/RI.html#ri-j",
                  "backend_controls": "https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/RI.html#keywords"}}
    try:
        m.require(m.sha(manifest_path) == MANIFEST_SHA, "Only the exact MPI-recovery manifest is admitted")
        manifest = m.validate(manifest_path)
        result["manifest"] = m.pin(manifest_path)
        rows, raw, evidence = {}, {}, {}
        for variant, job in (("rijcosx", "22436"), ("direct", "22437")):
            directory = getattr(args, variant).resolve()
            for name, digest in PINS[variant].items():
                m.verify({"path": str(directory / name), "sha256": digest})
            row = json.loads((directory / "RESULT.json").read_text())
            expected_status = "fixed_density_completed" if variant == "rijcosx" else "failed_native_admission"
            m.require(row["status"] == expected_status and row["allocation"]["job_id"] == job and
                      row["variant"] == variant and row["system_id"] == m.SYSTEM, "Unexpected original result")
            m.require(row["manifest"] == m.pin(manifest_path) and row["implementation"]["sha256"] == SOURCE_SHA and
                      row["native_calls_attempted"] == 1 and row["returncode"] == 0 and not row["timed_out"],
                      "Original native execution differs")
            expected_false = {"requested_jk_backend"} if variant == "direct" else set()
            m.require({k for k, v in row["checks"].items() if v is not True} == expected_false,
                      "Recovery is limited to the single auxiliary-metric parser veto")
            for record in list(row["input_pins"].values()) + list(row["output_pins"].values()):
                m.verify(record)
            m.require(row["input_pins"]["archive_seed.gbw"]["sha256"] == manifest["copied_gbw"]["sha256"] and
                      row["input_pins"]["geometry.xyz"]["sha256"] == manifest["copied_xyz"]["sha256"], "Imported state differs")
            v = manifest["variants"][variant]
            input_text = m.verify(v["input"]).read_text()
            m.require(row["input_pins"][v["basename"] + ".inp"]["sha256"] == v["input"]["sha256"], "Native input differs")
            raw[variant] = (directory / (v["basename"] + ".out")).read_text()
            parsed = m.parse_native(raw[variant], variant, input_text, directory)
            m.require(parsed["energies"] == row["energies"] and parsed["grid_evidence"] == row["grid_evidence"], "Saved native parse differs")
            m.require(parsed["checks"] == {k: val for k, val in row["checks"].items() if k != "process_success"},
                      "Original parser checks differ from receipt")
            property_text = (directory / (v["basename"] + ".property.txt")).read_text()
            bibtex = (directory / (v["basename"] + ".bibtex")).read_text()
            evidence[variant] = backend_evidence(m, variant, raw[variant], property_text, bibtex)
            # The sole corrected predicate. No native input, numerical gate,
            # source snapshot or old result is changed. Independent evidence
            # above is required even for these exact pinned native artifacts.
            if variant == "direct":
                parsed["checks"]["requested_jk_backend"] = True
            m.require(all(parsed["checks"].values()), "A separate native contract failed")
            rows[variant] = {"original_result": m.pin(directory / "RESULT.json"),
                             "original_status": row["status"], "original_checks": row["checks"],
                             "recovered_checks": parsed["checks"], "energies": parsed["energies"],
                             "grid_evidence": parsed["grid_evidence"], "native_backend_evidence": evidence[variant],
                             "new_native_calls": 0, "original_elapsed_seconds": row["elapsed_seconds"],
                             "original_allocated_cpu_seconds": row["allocated_cpu_seconds"],
                             "artifact_pins": {name: m.pin(directory / name) for name in PINS[variant]}}
        dft = [[g["settings"] for g in rows[v]["grid_evidence"]["printed_grid_sections"]
                if g["section"] == "DFT GRID GENERATION"] for v in m.VARIANTS]
        m.require(dft[0] == dft[1] and len(dft[0]) == 1, "Common DFT grid differs")
        same_density_observables = {}
        for label, pattern in (("one_electron_energy", r"^One Electron Energy:\s*(" + m.NUM + ")"),
                               ("xc_energy", r"^E\(XC\)\s*:\s*(" + m.NUM + ")"),
                               ("integrated_electrons", r"^N\(Total\)\s*:\s*(" + m.NUM + ")")):
            observed = {v: m.one_token(pattern, raw[v]) for v in m.VARIANTS}
            m.require(m.decimal(observed["direct"]) == m.decimal(observed["rijcosx"]), "Matched density/grid observable differs")
            same_density_observables[label] = observed
        contrasts = {}
        for key in ("electronic", "full_pbe0_d4", "d4_printed"):
            values = {v: m.decimal(rows[v]["energies"][key]["token"]) for v in m.VARIANTS}
            archive = m.decimal(manifest["archive_energies"][key]["token"])
            delta = values["direct"] - values["rijcosx"]
            contrasts[key] = {"direct_minus_rijcosx_hartree": str(delta),
                              "direct_minus_rijcosx_kcal_mol": str(delta * Decimal("627.509474")),
                              "rijcosx_minus_archive_hartree": str(values["rijcosx"] - archive),
                              "direct_minus_archive_hartree": str(values["direct"] - archive)}
        gate = Decimal(manifest["settings"]["archive_comparison_tolerance_hartree"])
        result.update(status="fixed_density_pair_completed_after_parser_recovery", rows=rows, contrasts=contrasts,
                      same_printed_density_observables=same_density_observables,
                      original_archive_energies=manifest["archive_energies"],
                      archive_grid_evidence=manifest["archive_grid_evidence"],
                      existing_archive_gate_hartree=str(gate),
                      replay_within_existing_archive_gate={k: abs(m.decimal(contrasts[k]["rijcosx_minus_archive_hartree"])) <= gate
                          for k in ("electronic", "full_pbe0_d4")},
                      exact_archived_final_cosx_grid_reproduction_established=False,
                      recovered_predicate="Only requested_jk_backend for exact job22437; startup V-metric preparation is not Fock-backend activation.",
                      limitation="This pair measures the native backend difference at imported orbitals under NOITER conventions. The RIJCOSX control differs from the archived post-final-COSX energy; report that residual separately. No converged direct SCF, physical-state qualification, baseline admission or complete cross-code equivalence follows.")
    except Exception as exc:
        result["exception"] = repr(exc)
        result["traceback"] = traceback.format_exc()
    result["elapsed_seconds"] = time.monotonic() - started
    result["allocated_cpu_seconds"] = result["elapsed_seconds"] * allocation["allocated_cpus"]
    m.write(output / "RESULTS.json", result)
    report = ["# Saved ORCA backend admission recovery", "", result["interpretation"], "", "Status: " + result["status"] + ".", ""]
    if result["status"] == "fixed_density_pair_completed_after_parser_recovery":
        report += ["The old failed parser result is preserved. The correction uses native SCF timings: direct has J and no SPLITRIJ/COSX; the control has SPLITRIJ/COSX and no J. Startup auxiliary-metric preparation is not an active-Fock marker.", "",
                   "| Energy | RIJCOSX (Eh) | Direct (Eh) | Direct − RIJCOSX (kcal/mol) |", "|---|---:|---:|---:|"]
        for key in ("electronic", "full_pbe0_d4", "d4_printed"):
            report.append(f"| {key} | {rows['rijcosx']['energies'][key]['token']} | {rows['direct']['energies'][key]['token']} | {contrasts[key]['direct_minus_rijcosx_kcal_mol']} |")
        report += ["", "Control − archived electronic energy: " + contrasts["electronic"]["rijcosx_minus_archive_hartree"] + " Eh.",
                   "Control within unchanged 1e-4 Eh electronic/full archive gates: " + json.dumps(result["replay_within_existing_archive_gate"]) + ".",
                   "", result["limitation"], "", "Zero new native calls. All original pins, failed checks, actual grids, timing evidence and costs are retained in RESULTS.json."]
    else:
        report += ["No admitted paired contrast. " + result.get("exception", "")]
    (output / "REPORT.md").write_text("\n".join(report) + "\n")
    print(json.dumps({"status": result["status"], "new_native_calls": 0}))
    return 0 if result["status"] == "fixed_density_pair_completed_after_parser_recovery" else 1


if __name__ == "__main__":
    sys.exit(main())
