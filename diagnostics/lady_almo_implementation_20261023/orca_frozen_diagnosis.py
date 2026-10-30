#!/usr/bin/env python3
"""Exactly two ORCA one-Fock energies on the archived La-aquo orbitals.

Preparation is metadata-only. Native execution and saved collection require
Slurm. These are fixed-density diagnostics, never converged-SCF references.
"""
from __future__ import annotations

import argparse
import datetime as dt
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import socket
import subprocess
import sys
import time
import traceback

SCHEMA = "orca_frozen_laaquo_ri_diagnostic_v1"
SYSTEM = "aquo8_central_La"
VARIANTS = ("rijcosx", "direct")
OLD = "diagnostics/lady_electronic_decomposition_20261023"
PROVENANCE_SHA256 = "677c759c1045e2d7d5ab2dd05977e4a516ca04ff8d520c40111669eabb452d66"
LEDGER_SHA256 = "24c699cc97c3c55857c56a5ee07dda6213eeaf12e59af69395f23bbd05b732c5"
SETTINGS = {"ranks": 8, "threads_per_rank": 1, "maxcore_mb_per_rank": 3000,
            "max_allocation_cpus": 8, "max_allocation_memory_mib": 32768,
            "max_job_wall_seconds": 3600, "native_timeout_seconds": 900,
            "method": "PBE0 D4", "grid": "DefGrid3", "convergence_keyword": "VeryTightSCF",
            "max_scf_iterations": 0, "maximum_native_calls_per_variant": 1,
            "archive_comparison_tolerance_hartree": "0.0001"}
NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
NOITER_PHRASES = (
    "INITIAL GUESS: MOREAD",
    "Guess MOs are being read from file: archive_seed.gbw",
    "Input Geometry matches current geometry (good)",
    "Input basis set matches current basis set (good)",
    "Occupation numbers from the input are taken because of NoIter",
    "MaxIter==0 has been chosen - skipping the SCF",
    "CalcGuessEnergy requested: will build one Fock matrix to compute the energy",
)
DOCS = {
    "ri_controls": "https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/RI.html#keywords",
    "fixed_density": "https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html#calculating-only-the-energy-of-an-input-density",
    "same_basis_restart": "https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html#restarting-scf-calculations",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": sha(path)}


def verify(record):
    path = Path(record["path"])
    require(path.is_file() and sha(path) == record["sha256"], "Missing/changed pinned file: " + str(path))
    return path


def write(path, value):
    with Path(path).open("x") as f:
        json.dump(value, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")


def decimal(token):
    value = Decimal(token.replace("D", "E").replace("d", "e"))
    require(value.is_finite(), "Nonfinite printed energy")
    return value


def energy_record(token):
    return {"token": token, "hartree": float(decimal(token))}


def tokens(pattern, text):
    return re.findall(pattern, text, re.M)


def one_token(pattern, text):
    values = tokens(pattern, text)
    require(len(values) == 1, "Expected exactly one printed field: " + pattern)
    return values[0]


def raw_grid_evidence(text):
    records, current = [], None
    for number, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if stripped == "DFT GRID GENERATION" or re.fullmatch(r"GRIDX\s+\d+", stripped):
            current = {"section": stripped, "line": number, "settings": []}
            records.append(current)
        if current is not None and re.search(
                r"IntAcc|RadialGrid|AngularGrid|GridPruning|WeightScheme|BFCut|WCut|"
                r"Partially contracted basis|Rotationally invariant|Total number of grid points|"
                r"Total number of batches|UseSFitting", line):
            current["settings"].append(line.strip())
        if "Grids setup in" in line:
            current = None
    phase = [{"line": i, "text": line} for i, line in enumerate(text.splitlines(), 1)
             if re.search(r"final.*integrat|integrat.*final|switch.*grid|grid.*switch", line, re.I)]
    return {"printed_grid_sections": records, "printed_phase_markers": phase,
            "active_cosx_energy_grid_inferred": False,
            "note": "Generated GRIDX 1/2/3 sections do not by themselves identify the grid used for this one-Fock energy."}


def render(original, variant):
    require(variant in VARIANTS, "Undeclared variant")
    old_first = "! PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad"
    require(original.splitlines()[0] == old_first, "Archived keyword line differs")
    first = old_first.replace("RIJCOSX", "RIJCOSX" if variant == "rijcosx" else "NOCOSX NORI")
    first = first.replace("EnGrad", "SP MOREAD CALCGUESSENERGY NOITER")
    result = original.replace(old_first, first, 1)
    require(result.count("%pal nprocs 4 end") == 1, "Archived MPI declaration differs")
    result = result.replace("%pal nprocs 4 end", "%pal nprocs 8 end", 1)
    old_scf = "%scf\n Guess PModel\n MaxIter 300\nend\n"
    require(result.count(old_scf) == 1, "Archived SCF block differs")
    new_scf = '%scf\n Guess MORead\n MaxIter 0\nend\n%moinp "archive_seed.gbw"\n'
    if variant == "direct":
        new_scf += "%method\n RI off\nend\n"
    result = result.replace(old_scf, new_scf, 1)
    require(result.split("%basis\n", 1)[1] == original.split("%basis\n", 1)[1],
            "Basis/ECP/auxiliary basis/geometry declaration changed")
    require(result.endswith("* xyzfile 3 1 geometry.xyz\n"), "Charge/state/geometry declaration differs")
    require("GuessMode" not in result and "%maxcore 3000\n" in result, "Unexpected projection or memory settings")
    return result


def prepare(repo, output):
    provenance_path = repo / OLD / "inputs/PROVENANCE.json"
    verify({"path": str(provenance_path), "sha256": PROVENANCE_SHA256})
    provenance = json.loads(provenance_path.read_text())
    saved = provenance["systems"][SYSTEM]
    selected = {name: saved[name] for name in (
        "xyz", "original_xyz", "original_executed_xyz", "original_native_input",
        "archived_native_input_copy", "original_full_result", "original_output",
        "original_engrad", "original_panel", "original_adapter")}
    selected["provenance"] = pin(provenance_path)
    selected["energy_ledger"] = {"path": str(repo / OLD / "cp2k_inputs_v1/ENERGY_LEDGER.json"),
                                  "sha256": LEDGER_SHA256}
    for kind, entries in provenance["files"]["metal_assets"]["La"].items():
        for role, record in entries.items():
            selected[f"La_{kind}_{role}"] = record
    for record in selected.values():
        verify(record)
    original = verify(selected["original_native_input"]).read_text()
    basis = "\n".join(verify(selected[f"La_{kind}_copy"]).read_text() for kind in ("basis", "aux"))
    require(original.split("%basis\n", 1)[1] == basis + "\nend\n* xyzfile 3 1 geometry.xyz\n",
            "Executed basis/ECP/auxiliary text differs from pinned assets")
    require(len({selected[k]["sha256"] for k in ("xyz", "original_xyz", "original_executed_xyz")}) == 1,
            "Geometry copies differ")
    require(selected["original_native_input"]["sha256"] == selected["archived_native_input_copy"]["sha256"],
            "Native input archive differs")
    original_result = json.loads(verify(selected["original_full_result"]).read_text())
    require(original_result["status"] == "complete" and all(original_result["checks"].values()),
            "Original reference was not admitted")
    archived = verify(selected["original_output"]).read_text()
    require("Program Version 6.1.1" in archived and "ORCA TERMINATED NORMALLY" in archived,
            "Original ORCA version/termination differs")
    require(tokens(r"Number of eigenvalues below threshold\s+\.{2,}\s+(\d+)", archived) == ["0"],
            "NOITER is unsupported when overlap components were removed")
    gbw = verify(selected["original_native_input"]).with_suffix(".gbw")
    require(gbw.is_file() and gbw.stat().st_size > 100000, "Saved archived orbitals are missing/truncated")
    selected["original_gbw"] = pin(gbw)
    panel = json.loads(verify(selected["original_panel"]).read_text())
    native = verify(panel["orca"])
    runtime = {"orca": panel["orca"]}
    for name in ("orca_startup", "orca_startup_mpi", "orca_guess", "orca_guess_mpi",
                 "orca_leanscf", "orca_leanscf_mpi", "orca_prop", "orca_prop_mpi",
                 "orca_pop", "orca_pop_mpi"):
        runtime[name] = pin(native.parent / name)
    for path in sorted((native.parent / "lib").glob("liborca_tools*")):
        runtime[path.name] = pin(path)
    ledger = json.loads(verify(selected["energy_ledger"]).read_text())["systems"][SYSTEM]["energy_hartree"]
    archive_energies = {
        "electronic": energy_record(ledger["electronic_after_final_cosx_integration"]["token"]),
        "full_pbe0_d4": energy_record(ledger["full_pbe0_d4"]["token"]),
        "d4_printed": energy_record(ledger["d4_printed"]["token"]),
    }
    require(decimal(one_token(r"^Total Energy\s*:\s*(" + NUM + r")\s+Eh", archived)) ==
            decimal(archive_energies["electronic"]["token"]), "Archived electronic energy/ledger differs")
    require(decimal(one_token(r"^FINAL SINGLE POINT ENERGY\s+(" + NUM + r")", archived)) ==
            decimal(archive_energies["full_pbe0_d4"]["token"]), "Archived full energy/ledger differs")
    output.mkdir(parents=True, exist_ok=False)
    (output / "source").mkdir()
    shutil.copyfile(__file__, output / "source/orca_frozen_diagnosis.py")
    shutil.copyfile(gbw, output / "archive_seed.gbw")
    shutil.copyfile(verify(selected["original_executed_xyz"]), output / "geometry.xyz")
    variants = {}
    for name in VARIANTS:
        path = output / f"frozen_{name}.inp"
        path.write_text(render(original, name))
        variants[name] = {"input": pin(path), "basename": path.stem}
    manifest = {"schema": SCHEMA, "created_utc": stamp(), "system_id": SYSTEM,
                "settings": SETTINGS, "source_inputs": selected, "runtime_files": runtime,
                "implementation": pin(output / "source/orca_frozen_diagnosis.py"),
                "copied_gbw": pin(output / "archive_seed.gbw"), "copied_xyz": pin(output / "geometry.xyz"),
                "variants": variants, "archive_energies": archive_energies,
                "archive_grid_evidence": raw_grid_evidence(archived),
                "archive_elapsed_seconds": original_result["elapsed_seconds"],
                "archive_allocated_cpu_seconds": original_result["allocated_cpu_seconds"],
                "documentation": DOCS, "new_native_calls": 0,
                "interpretation": "Two fixed-density diagnostics; no converged-SCF, ALMO-admission, or production claim.",
                "grid_phase_policy": "Preserve DefGrid3; record actual grids. Do not force a COSX phase or assert archived final-energy replay from the input alone."}
    require(manifest["copied_gbw"]["sha256"] == selected["original_gbw"]["sha256"], "Copied GBW differs")
    write(output / "MANIFEST.json", manifest)
    validate(output / "MANIFEST.json")
    print(json.dumps({"status": "prepared_zero_native_calls", "manifest": pin(output / "MANIFEST.json")}))


def validate(manifest_path):
    m = json.loads(manifest_path.read_text())
    require(m["schema"] == SCHEMA and m["system_id"] == SYSTEM and m["settings"] == SETTINGS,
            "Diagnostic scope/settings differ")
    require(set(m["variants"]) == set(VARIANTS), "Variant inventory differs")
    for record in list(m["source_inputs"].values()) + list(m["runtime_files"].values()):
        verify(record)
    verify(m["implementation"])
    require(sha(__file__) == m["implementation"]["sha256"], "Runner differs from prepared source")
    for name in ("copied_gbw", "copied_xyz"):
        verify(m[name])
    require(m["copied_gbw"]["sha256"] == m["source_inputs"]["original_gbw"]["sha256"], "Seed GBW differs")
    require(m["copied_xyz"]["sha256"] == m["source_inputs"]["original_executed_xyz"]["sha256"], "Copied XYZ differs")
    original = verify(m["source_inputs"]["original_native_input"]).read_text()
    for name in VARIANTS:
        v = m["variants"][name]
        require(v["basename"] == "frozen_" + name and verify(v["input"]).read_text() == render(original, name),
                "Prepared input differs from exact transformation")
    return m


def slurm_allocation(native=False):
    require(os.environ.get("SLURM_JOB_ID", "").isdigit(), "Execution/collection requires Slurm")
    require(os.environ.get("SLURM_JOB_QOS") in (None, "standard_expand"), "Require standard_expand QoS")
    ntasks = int(os.environ.get("SLURM_NTASKS", "1"))
    cpus = ntasks * int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    memory = int(os.environ.get("SLURM_MEM_PER_NODE", "0"))
    if not memory:
        memory = int(os.environ.get("SLURM_MEM_PER_CPU", "0")) * cpus
    if native:
        require(cpus == 8, "This fixed input requires exactly eight allocated CPUs")
        require(24000 < memory <= 32768, "Require >24000 and <=32768 MiB granted memory")
    return {"job_id": os.environ["SLURM_JOB_ID"], "allocated_cpus": cpus,
            "allocated_memory_mib": memory, "host": socket.gethostname()}


def parse_native(text, variant, input_text, output):
    echoed = tokens(r"^\|\s*\d+> ?(.*)$", text)
    expected = input_text.splitlines()
    checks = {"normal_termination": text.count("ORCA TERMINATED NORMALLY") == 1,
              "orca_version": tokens(r"Program Version\s+(\S+)", text) == ["6.1.1"],
              "native_fixed_density_contract": all(text.count(p) == 1 for p in NOITER_PHRASES),
              "no_scf_iterations_or_gradient": not re.search(
                  r"SCF NOT CONVERGED|SCF CONVERGED AFTER|^\s*SCF ITERATIONS\s*$|ORCA SCF GRADIENT CALCULATION",
                  text, re.M) and not list(output.glob("*.engrad")),
              "no_basis_or_geometry_projection": not re.search(
                  r"Input basis set (?:is compatible with but different|needs projection)|"
                  r"Input geometry compatible with but different|Occupation numbers will be reassigned", text),
              "unchanged_input_echo": echoed[:len(expected)] == expected and
                  [line.strip() for line in echoed[len(expected):]] == ["", "****END OF INPUT****"]}
    for name, pattern, expected in (
        ("charge", r"Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)", "3"),
        ("multiplicity", r"Multiplicity\s+Mult\s+\.{2,}\s+(\d+)", "1"),
        ("electrons", r"Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)", "88"),
        ("basis_dimension", r"Basis Dimension\s+Dim\s+\.{2,}\s+(\d+)", "391"),
        ("shell_count", r"Number of shells\s+\.{2,}\s+(\d+)", "167"),
        ("zero_iterations", r"Maximum # iterations\s+MaxIter\s+\.{2,}\s+(\d+)", "0"),
        ("no_removed_overlap_modes", r"Number of eigenvalues below threshold\s+\.{2,}\s+(\d+)", "0"),
    ):
        checks[name] = tokens(pattern, text) == [expected]
    checks["ecp46"] = tokens(r"Type\s+La\s+ECP\s*\(replacing\s+(\d+)\s+core electrons,\s*lmax=(\d+)\)", text) == [("46", "4")]
    checks["pbe0_exchange"] = tokens(r"Exchange Functional\s+Exchange\s+\.{2,}\s+(\S+)", text) == ["PBE"]
    checks["pbe0_correlation"] = tokens(r"Correlation Functional\s+Correlation\s+\.{2,}\s+(\S+)", text) == ["PBE"]
    checks["pbe0_hybrid_fraction"] = tokens(r"Fraction HF Exchange\s+ScalHFX\s+\.{2,}\s+(\S+)", text) == ["0.250000"]
    ri_on = "RI-approximation to the Coulomb term is turned on" in text
    cosx_on = bool(re.search(r"RIJ-COSX.*\.{2,}\s+on", text))
    checks["requested_jk_backend"] = ((ri_on and cosx_on) if variant == "rijcosx" else
        (not ri_on and not cosx_on and "Calculating RI/J V-Matrix" not in text and "COSX GRID GENERATION" not in text))
    grids = raw_grid_evidence(text)
    checks["dft_grid_reported"] = any(g["section"] == "DFT GRID GENERATION" and
        any("Total number of grid points" in s for s in g["settings"]) for g in grids["printed_grid_sections"])
    energies = {}
    for name, pattern in (
        ("electronic", r"^Total Energy\s*:\s*(" + NUM + r")\s+Eh"),
        ("full_pbe0_d4", r"^FINAL SINGLE POINT ENERGY\s+(" + NUM + r")"),
        ("d4_printed", r"^\s*Dispersion correction\s+(" + NUM + r")"),
    ):
        values = tokens(pattern, text)
        checks["one_finite_" + name] = len(values) == 1
        if len(values) == 1:
            energies[name] = energy_record(values[0])
    if len(energies) == 3:
        numbers = {k: decimal(v["token"]) for k, v in energies.items()}
        rounding = sum(Decimal("0.5") * Decimal(10) ** v.as_tuple().exponent for v in numbers.values())
        residual = numbers["full_pbe0_d4"] - numbers["electronic"] - numbers["d4_printed"]
        checks["d4_closure_within_printed_precision"] = abs(residual) <= rounding
    else:
        checks["d4_closure_within_printed_precision"] = False
    checks["native_d4_present"] = "DFTD4" in text
    return {"checks": checks, "energies": energies, "grid_evidence": grids,
            "orbital_import_evidence": [{"line": i, "text": line} for i, line in enumerate(text.splitlines(), 1)
                if re.search(r"INITIAL GUESS:|Guess MOs|Input [Gg]eometry|Input basis|Occupation numbers|MOs were|MaxIter==0|CalcGuessEnergy requested", line)],
            "native_component_lines": [{"line": i, "text": line} for i, line in enumerate(text.splitlines(), 1)
                if re.search(r"^\s*(Nuclear Repulsion|One Electron Energy|Two Electron Energy|E\(X\)|E\(C\)|E\(XC\)|N\(Total\)|RI-approximation|RIJ-COSX|TOTAL RUN TIME)", line)],
            "scf_converged": False, "scf_optimization_performed": False,
            "density_bitwise_identity_claimed": False,
            "density_note": "Same-basis native MORead can renormalize/reorthogonalize; no orbital optimization or projection to a changed basis is admitted."}


def run(manifest_path, variant, output):
    allocation = slurm_allocation(native=True)
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    result = {"schema": SCHEMA, "status": "failed", "variant": variant, "system_id": SYSTEM,
              "interpretation": "Fixed-density one-Fock diagnostic, not converged SCF.",
              "started_utc": stamp(), "allocation": allocation, "native_calls_attempted": 0,
              "manifest": pin(manifest_path), "implementation": pin(__file__),
              "baseline_admitted": False, "production_changed": False}
    try:
        m = validate(manifest_path)
        v = m["variants"][variant]
        for record, name in ((m["copied_gbw"], "archive_seed.gbw"), (m["copied_xyz"], "geometry.xyz"),
                             (v["input"], v["basename"] + ".inp")):
            shutil.copyfile(verify(record), output / name)
            require(sha(output / name) == record["sha256"], "Staged input changed: " + name)
        native = verify(m["runtime_files"]["orca"])
        env = os.environ.copy()
        for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
            env[name] = "1"
        env.update(OMP_DYNAMIC="FALSE", MKL_DYNAMIC="FALSE", OMP_MAX_ACTIVE_LEVELS="1")
        env["PATH"] = str(native.parent) + os.pathsep + env.get("PATH", "")
        command = [str(native), v["basename"] + ".inp", "--bind-to none"]
        result["command"] = command
        result["environment"] = {k: env.get(k) for k in ("PATH", "LD_LIBRARY_PATH", "OMP_NUM_THREADS",
            "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_DYNAMIC", "OMP_MAX_ACTIVE_LEVELS")}
        result["input_pins"] = {name: pin(output / name) for name in
            ("archive_seed.gbw", "geometry.xyz", v["basename"] + ".inp")}
        result["runtime_files"] = m["runtime_files"]
        write(output / "START.json", result)
        # One durable reservation per manifest/variant, also if native startup fails.
        write(manifest_path.parent / ("ATTEMPT_" + variant + ".json"),
              {"utc": stamp(), "manifest": pin(manifest_path), "output": str(output),
               "variant": variant, "job_id": allocation["job_id"], "implementation": pin(__file__)})
        result["native_calls_attempted"] = 1
        result["timed_out"] = False
        native_started = time.monotonic()
        with (output / (v["basename"] + ".out")).open("x") as stdout, (output / (v["basename"] + ".err")).open("x") as stderr:
            proc = subprocess.Popen(command, cwd=output, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
            try:
                result["returncode"] = proc.wait(timeout=SETTINGS["native_timeout_seconds"])
            except subprocess.TimeoutExpired:
                result["timed_out"] = True
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                result["returncode"] = proc.returncode
        result["native_elapsed_seconds"] = time.monotonic() - native_started
        for record in result["input_pins"].values():
            verify(record)
        validate(manifest_path)  # Recheck originals, copied seed and native files after execution.
        raw = (output / (v["basename"] + ".out")).read_text(errors="replace")
        parsed = parse_native(raw, variant, verify(v["input"]).read_text(), output)
        result.update(parsed)
        result["checks"]["process_success"] = result["returncode"] == 0 and not result["timed_out"]
        result["output_pins"] = {p.name: pin(p) for p in sorted(output.iterdir()) if p.is_file() and
            p.name not in ("START.json", "RESULT.json")}
        result["status"] = "fixed_density_completed" if all(result["checks"].values()) else "failed_native_admission"
    except Exception as exc:
        result["exception"] = repr(exc)
        result["traceback"] = traceback.format_exc()
    result["finished_utc"] = stamp()
    result["elapsed_seconds"] = time.monotonic() - started
    result["allocated_cpu_seconds"] = result["elapsed_seconds"] * allocation["allocated_cpus"]
    write(output / "RESULT.json", result)
    print(json.dumps({"status": result["status"], "variant": variant, "native_calls_attempted": result["native_calls_attempted"]}))
    return 0 if result["status"] == "fixed_density_completed" else 1


def collect(manifest_path, paths, output):
    allocation = slurm_allocation()
    output.mkdir(parents=True, exist_ok=False)
    result = {"schema": SCHEMA, "status": "failed", "created_utc": stamp(), "system_id": SYSTEM,
              "allocation": allocation, "new_native_calls": 0, "manifest": pin(manifest_path),
              "interpretation": "Paired fixed-density diagnostics, not converged SCF or ALMO admission.",
              "baseline_admitted": False, "production_changed": False}
    try:
        m = validate(manifest_path)
        rows = {}
        for variant in VARIANTS:
            path = paths[variant]
            row = json.loads((path / "RESULT.json").read_text())
            require(row["schema"] == SCHEMA and row["variant"] == variant and row["system_id"] == SYSTEM,
                    "Unexpected saved variant")
            require(row["status"] == "fixed_density_completed" and all(row["checks"].values()), "Unadmitted native diagnostic")
            require(row["manifest"] == pin(manifest_path) and row["native_calls_attempted"] == 1 and
                    row["implementation"]["sha256"] == m["implementation"]["sha256"], "Saved run/source mismatch")
            for record in list(row["input_pins"].values()) + list(row["output_pins"].values()):
                verify(record)
            require(row["input_pins"]["archive_seed.gbw"]["sha256"] == m["copied_gbw"]["sha256"], "Saved seed differs")
            v = m["variants"][variant]
            parsed = parse_native((path / (v["basename"] + ".out")).read_text(errors="replace"),
                                  variant, verify(v["input"]).read_text(), path)
            require(all(parsed["checks"].values()) and parsed["energies"] == row["energies"] and
                    parsed["grid_evidence"] == row["grid_evidence"], "Native saved replay differs")
            rows[variant] = row
        dft = [[g["settings"] for g in rows[v]["grid_evidence"]["printed_grid_sections"]
                if g["section"] == "DFT GRID GENERATION"] for v in VARIANTS]
        require(dft[0] == dft[1] and len(dft[0]) == 1, "Printed common DFT grid settings differ")
        contrasts = {}
        for key in ("electronic", "full_pbe0_d4", "d4_printed"):
            values = {v: decimal(rows[v]["energies"][key]["token"]) for v in VARIANTS}
            archive_value = decimal(m["archive_energies"][key]["token"])
            delta = values["direct"] - values["rijcosx"]
            contrasts[key] = {"direct_minus_rijcosx_hartree": str(delta),
                              "direct_minus_rijcosx_kcal_mol": str(delta * Decimal("627.509474")),
                              "rijcosx_minus_archive_hartree": str(values["rijcosx"] - archive_value),
                              "direct_minus_archive_hartree": str(values["direct"] - archive_value)}
        gate = Decimal(SETTINGS["archive_comparison_tolerance_hartree"])
        result.update(status="fixed_density_pair_completed", rows=rows, contrasts=contrasts,
                      archive_energies=m["archive_energies"], archive_grid_evidence=m["archive_grid_evidence"],
                      result_pins={v: pin(paths[v] / "RESULT.json") for v in VARIANTS},
                      replay_within_existing_archive_gate={k: abs(decimal(contrasts[k]["rijcosx_minus_archive_hartree"])) <= gate
                          for k in ("electronic", "full_pbe0_d4")},
                      existing_archive_gate_hartree=str(gate),
                      exact_archived_final_cosx_grid_reproduction_established=False,
                      limitation="The pair measures the operator change under native NOITER conventions at imported orbitals. Generated grid lists alone do not establish active COSX phase. Archive-gate agreement is reported separately; it does not qualify a converged direct SCF, explain other backend differences, or admit the failed PySCF baseline.")
    except Exception as exc:
        result["exception"] = repr(exc)
        result["traceback"] = traceback.format_exc()
    write(output / "RESULTS.json", result)
    report = ["# ORCA frozen-density diagnosis", "", result["interpretation"], "", "Status: " + result["status"] + ".", ""]
    if result["status"] == "fixed_density_pair_completed":
        report += ["| Energy | RIJCOSX (Eh) | Direct J/K (Eh) | Direct − RIJCOSX (kcal/mol) |", "|---|---:|---:|---:|"]
        for key in ("electronic", "full_pbe0_d4", "d4_printed"):
            report.append(f"| {key} | {rows['rijcosx']['energies'][key]['token']} | {rows['direct']['energies'][key]['token']} | {contrasts[key]['direct_minus_rijcosx_kcal_mol']} |")
        report += ["", "RIJCOSX − archived electronic energy: " + contrasts["electronic"]["rijcosx_minus_archive_hartree"] + " Eh.",
                   "Replay within existing 1e-4 Eh electronic/full archive gates: " + json.dumps(result["replay_within_existing_archive_gate"]) + ".",
                   "", result["limitation"], "", "Actual grid sections, phase markers, orbital-import messages, raw artifact pins and costs are retained in RESULTS.json."]
    else:
        report += ["No admitted paired contrast is available. " + result.get("exception", "")]
    (output / "REPORT.md").write_text("\n".join(report) + "\n")
    return 0 if result["status"] == "fixed_density_pair_completed" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--repo", type=Path, default=Path(os.environ.get("LADY_REPO_ROOT", Path(__file__).resolve().parents[2])))
    p.add_argument("--output", type=Path, required=True)
    for action in ("validate", "run", "collect"):
        p = sub.add_parser(action)
        p.add_argument("--manifest", type=Path, required=True)
        if action != "validate":
            p.add_argument("--output", type=Path, required=True)
        if action == "run":
            p.add_argument("--variant", choices=VARIANTS, required=True)
        if action == "collect":
            p.add_argument("--rijcosx", type=Path, required=True)
            p.add_argument("--direct", type=Path, required=True)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.repo.resolve(), args.output.resolve())
    elif args.action == "validate":
        validate(args.manifest.resolve())
        print("Validated exact source/input pins; zero native calls")
    elif args.action == "run":
        return run(args.manifest.resolve(), args.variant, args.output.resolve())
    else:
        return collect(args.manifest.resolve(), {v: getattr(args, v).resolve() for v in VARIANTS}, args.output.resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main())
