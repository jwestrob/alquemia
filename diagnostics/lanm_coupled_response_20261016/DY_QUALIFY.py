#!/usr/bin/env python3
"""Stage only the eight 7CCO La/Dy native/standalone vacuum/ALPB controls.

No engine, scheduler, package discovery, or molecular calculation is invoked.
Usage: python DY_QUALIFY.py --config ROOT_CONFIG.json --output FRESH_DIRECTORY

Required root config keys (no numerical defaults):
  protocol_id; plan: {path, sha256}; source_xyz: {La: pin, Dy: pin};
  native_settings: {maxcore_mb, maxiter, tole, electronic_temperature_K};
  standalone_settings: {accuracy, maxiter, electronic_temperature_K,
                        alpb_reference: "gsolv"};
  tolerances: nonempty root-declared acceptance specification.
The fixed scope is 48 atoms, charge -1, effective multiplicity 1 / UHF 0;
physical Dy(III) multiplicity 6 is retained separately. Tolerances are copied
and pinned, never selected or applied by this preparation adapter.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import copy
import json
import math
import os
from pathlib import Path
import sys
import tomllib

REPO = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs")
sys.path.insert(0, str(REPO / "scripts"))
import native_gfn2_parameters as guard

BASE_SHA = "a4f83ca72d58d6e80f634e61bb8e99b7af0f26c38cc4a9295e85215ecf72d05c"
BASE = REPO / "params/native_gfn2_parameters_v1/corrected_la.xtb.json"
RUNTIME = REPO / "workspaces/orca_parameter_audit_20261012/runtime_qualification_v1/RUNTIME_BINDING.json"
OBSERVED = REPO / "diagnostics/orca_parameter_audit_20261012/parameter_guard/OBSERVED_RUNTIME_SNAPSHOT.json"
XTB_RUNTIME = REPO / "workspaces/orca_parameter_audit_20261012/context_gfn2_v1/RUNTIME.json"
CORE = REPO / "workspaces/lanm_coupled_response_20261016/lbt3_preparation_v1/execution/7CCO"
FIELDS = ["kcn", "shpoly", "lgam"]


def write_json(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


@contextmanager
def pinned_environment(environment):
    saved = {key: os.environ.get(key) for key in environment}
    try:
        for key, value in environment.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        yield
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def positive(value, name, integer=False):
    guard.require(type(value) in ((int,) if integer else (int, float))
                  and math.isfinite(value) and value > 0, f"invalid {name}")


def stage(config_path, output):
    cfg = guard.load_json(config_path)
    guard.require(isinstance(cfg.get("protocol_id"), str) and cfg["protocol_id"], "protocol_id required")
    guard.check_pin(cfg["plan"])
    guard.require(isinstance(cfg.get("tolerances"), dict) and cfg["tolerances"], "root tolerances required")
    native, standalone = cfg["native_settings"], cfg["standalone_settings"]
    for name in ("maxcore_mb", "maxiter"):
        positive(native[name], "native " + name, integer=True)
    positive(native["tole"], "native tole")
    positive(native["electronic_temperature_K"], "native electronic temperature")
    positive(standalone["accuracy"], "standalone accuracy")
    positive(standalone["maxiter"], "standalone maxiter", integer=True)
    positive(standalone["electronic_temperature_K"], "standalone electronic temperature")
    guard.require(native["electronic_temperature_K"] == standalone["electronic_temperature_K"],
                  "electronic temperatures must match")
    guard.require(standalone["alpb_reference"] == "gsolv", "only declared gsolv reference supported")
    guard.require(set(cfg["source_xyz"]) == {"La", "Dy"}, "exact La/Dy sources required")
    source = {}
    for metal in ("La", "Dy"):
        source[metal] = guard.check_pin(cfg["source_xyz"][metal])
        guard.require(source[metal].resolve() == (CORE / f"matched_core_{metal}.xyz").resolve(),
                      "source outside fixed 7CCO matched-core pair")
        lines = source[metal].read_text().splitlines()
        guard.require(lines[0].strip() == "48" and len(lines) == 50, "expected exactly 48 atoms")
        symbols = [line.split()[0] for line in lines[2:]]
        guard.require(symbols[0] == metal and symbols.count(metal) == 1
                      and set(symbols[1:]) <= {"H", "C", "N", "O"}, "unexpected core inventory")
        guard.xyz_elements(source[metal])
    la = [line.split() for line in source["La"].read_text().splitlines()[2:]]
    dy = [line.split() for line in source["Dy"].read_text().splitlines()[2:]]
    guard.require(la[0][1:] == dy[0][1:] and la[1:] == dy[1:], "metal pair must share exact coordinates/order")
    guard.require(guard.sha256(BASE) == BASE_SHA, "corrected-La full profile changed")
    guard.require(guard.sha256(guard.AUTHOR) == guard.AUTHOR_SHA, "author TOML changed")
    base = guard.load_json(BASE)
    candidate = copy.deepcopy(base)
    author = tomllib.loads(guard.AUTHOR.read_text())["element"]["Dy"]
    shells = candidate["element"]["Dy"]["shells"]
    guard.require(shells == ["5d", "6s", "6p"], "unexpected declared Dy shells")
    for field in FIELDS:
        candidate["element"]["Dy"][field] = [author[field][author["shells"].index(shell)] for shell in shells]
    changes = guard._diff(base, candidate)
    guard.require(len(changes) == 9 and all(change["path"][:2] == ["element", "Dy"]
                  and change["path"][2] in FIELDS for change in changes), "only nine Dy coefficients may change")
    xtb = guard.load_json(XTB_RUNTIME)
    for record in [xtb["xtb"], xtb["xtb_parameter_file"], *xtb["xtb_dynamic"]["libraries"]]:
        guard.check_pin(record)
    output.mkdir(parents=True, exist_ok=False)
    frozen_config = output / "ROOT_CONFIG.json"
    with frozen_config.open("xb") as handle:
        handle.write(config_path.read_bytes())
    profile, policy_path = output / "corrected_la_dy.xtb.json", output / "la_dy_shell_three_arrays.json"
    write_json(profile, candidate)
    policy = {"schema_version": 1, "correction_id": cfg["protocol_id"] + "__La_Dy_author_shell_arrays",
              "original_sha256": guard.ORIGINAL_SHA, "author_sha256": guard.AUTHOR_SHA,
              "candidate_sha256": guard.sha256(profile), "corrections": {"La": FIELDS, "Dy": FIELDS}}
    write_json(policy_path, policy)
    audit = guard.audit_parameters(profile, policy)
    guard.require(all(audit["elements"][e]["parameter_consistent"] for e in ("H", "C", "N", "O", "La", "Dy")),
                  "a used element is inconsistent with author shell labels")
    write_json(output / "PARAMETER_AUDIT.json", audit)
    # Match scripts/context_factorial_gfn2.py::clean_env without importing its
    # module-level environment mutations or any numerical backend.
    xtb_environment = {key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
    xtb_environment.update(OMP_DYNAMIC="FALSE", OMP_MAX_ACTIVE_LEVELS="1",
                           XTBPATH=str(Path(xtb["xtb_parameter_file"]["path"]).parent),
                           LD_PRELOAD=None, LD_LIBRARY_PATH=None, PYTHONPATH=None, PYTHONHOME=None)
    common_sources = [frozen_config, config_path, Path(cfg["plan"]["path"]), Path(__file__),
                      BASE, guard.AUTHOR, RUNTIME, OBSERVED, XTB_RUNTIME,
                      REPO / "scripts/context_factorial_gfn2.py"]
    tasks = []
    for metal in ("La", "Dy"):
        for medium in ("vacuum", "ALPB(Water)"):
            for backend in ("native", "standalone"):
                task_id = "__".join((metal, "vacuum" if medium == "vacuum" else "water", backend))
                directory = output / task_id
                directory.mkdir()
                xyz = directory / "core.xyz"
                with xyz.open("xb") as handle:
                    handle.write(source[metal].read_bytes())
                task = {"id": task_id, "backend": backend, "metal": metal, "medium": medium,
                        "directory": str(directory), "xyz": guard.pin(xyz), "source_xyz": cfg["source_xyz"][metal],
                        "atom_count": 48, "charge": -1, "multiplicity": 1, "uhf": 0,
                        "physical_multiplicity": 6 if metal == "Dy" else 1,
                        "state_convention": guard.STATE_CONVENTION, "runtime_qualified": False,
                        "tolerances": cfg["tolerances"], "status": "staged_not_executed"}
                if backend == "native":
                    inp = directory / "endpoint.inp"
                    solvent = " ALPB(water)" if medium != "vacuum" else ""
                    with inp.open("x") as handle:
                        handle.write(f"! Native-GFN2-xTB NoAutostart{solvent}\n%maxcore {native['maxcore_mb']}\n"
                                     + guard.import_directives()
                                     + f"%scf\n MaxIter {native['maxiter']}\n TolE {native['tole']}\n"
                                     + f" SmearTemp {native['electronic_temperature_K']}\n UseXTBMixer true\nend\n"
                                     + "%pal\n nprocs 1\nend\n* xyzfile -1 1 core.xyz\n")
                    binding = guard.build_binding(input_path=inp, xyz=xyz, parameters=profile,
                        runtime_binding=RUNTIME, observed_runtime=OBSERVED,
                        source_inputs=[*common_sources, source[metal]], protocol_id=cfg["protocol_id"],
                        charge=-1, multiplicity=1, medium=medium)
                    with pinned_environment(binding["environment"]):
                        task["guard"] = guard.stage(binding, profile, policy, mode="diagnostic", variant="corrected")
                    task.update(command=[binding["executable"]["path"], inp.name, "--bind-to none"],
                                environment=binding["environment"], input=guard.pin(inp))
                else:
                    command = [xtb["xtb"]["path"], xyz.name, "--gfn", "2", "--chrg", "-1", "--uhf", "0",
                               "--etemp", str(standalone["electronic_temperature_K"]), "--acc", str(standalone["accuracy"]),
                               "--iterations", str(standalone["maxiter"]), "--norestart", "--json", "--pop"]
                    if medium != "vacuum":
                        command += ["--alpb", "water", standalone["alpb_reference"]]
                    task.update(command=command, environment=xtb_environment,
                                executable=xtb["xtb"], parameters=xtb["xtb_parameter_file"],
                                runtime_libraries=xtb["xtb_dynamic"]["libraries"])
                write_json(directory / "TASK.json", task)
                tasks.append({**task, "task_record": guard.pin(directory / "TASK.json")})
    manifest = {"schema": "lbt3_dy_gfn2_qualification_stage.v1", "protocol_id": cfg["protocol_id"],
                "status": "staged_not_executed", "numerical_calls": 0, "runtime_qualified": False,
                "implementation": guard.pin(__file__), "config": guard.pin(frozen_config),
                "plan": cfg["plan"], "profile": guard.pin(profile), "policy": guard.pin(policy_path),
                "base_profile": guard.pin(BASE), "changes_from_corrected_La": changes,
                "source_pins": [guard.pin(path) for path in common_sources], "tasks": tasks,
                "scope": "eight fixed conditional runtime controls; no affinity, solvent score, or production admission"}
    guard.require(len(tasks) == 8, "exact eight-cell matrix required")
    write_json(output / "MANIFEST.json", manifest)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    stage(args.config.resolve(), args.output.resolve())
