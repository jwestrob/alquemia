"""Pinned TABI-PB adapter with a physical cavity independent of QM charge sites.

The native solver calls NanoShaper through PATH.  A preserved bridge replaces
only its generated XYZR with the pinned source-atom cavity, runs NanoShaper,
and retains the mesh before TABI deletes it.  No solver equations are changed.
"""
from __future__ import annotations

import argparse
import copy
import csv
import fcntl
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

import numpy as np

from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new
from affordable_environment import KJ_PER_KCAL, physical_boundary_key, pqr_lines

# Capture the bytes associated with this import, before any long-running
# subprocess. A later live-file edit must not relabel the loaded executor.
LOADED_IMPLEMENTATION = record(__file__)

PROTOCOL = "source_cavity_frozen_charge_tabi_v1"
TABI_SOURCE_COMMIT = "fe1c237b057418fed48535db125394607040d9de"
BINARY_HASHES = {
    "tabipb": "2e3b2ef3016982438401896ad0817fda5095b3042dc23a601e4e9e96406e0a9b",
    "NanoShaper": "6d41b6f0ddd400aa5713d08b3a178884bfff4c7ae826530ab73549aa93584baf",
}
PHYSICS = {"solute_dielectric": 1.0, "solvent_dielectric": 78.54,
           "salt_molar": 0.0, "temperature_K": 298.15,
           "probe_radius_A": 1.4, "radii_policy": "Bondi_CHNOS_common_1p8A_metal_v1"}
LEVELS = {
    "primary": {"grid_scale_inverse_A": 2.0, "tree_degree": 5, "tree_theta": 0.5},
    "refined": {"grid_scale_inverse_A": 3.0, "tree_degree": 5, "tree_theta": 0.5},
    "tree": {"grid_scale_inverse_A": 2.0, "tree_degree": 7, "tree_theta": 0.3},
}
NUMERICAL_CONSTANTS = {"tree_max_per_leaf": 500, "precondition": False,
                       "gmres_relative_tolerance": 1e-4, "gmres_max_iterations": 100,
                       "gmres_restart": 10, "mesh": "ses", "threads": 1}
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def install_software(source_bin, output):
    """Copy exact archived executables; never chmod or modify their sources."""
    output = Path(output).resolve()
    if output.exists():
        raise InvalidArtifact(f"refusing existing software directory: {output}")
    originals = {}
    for name, expected in BINARY_HASHES.items():
        item = record(Path(source_bin) / name)
        if item["sha256"] != expected:
            raise InvalidArtifact(f"unreviewed executable: {name}")
        originals[name] = item
    output.mkdir(parents=True)
    executables = {}
    for name, rec in originals.items():
        target = output / name
        shutil.copyfile(verify(rec), target)
        target.chmod(0o755)
        executables[name] = record(target)
    result = {"status": "copied_not_scientifically_executed", "executables": executables,
              "archived_sources": originals, "tabi_source_commit": TABI_SOURCE_COMMIT}
    write_new(output / "software_manifest.json", result)
    return result


def transform_definition(name):
    rotation = np.eye(3)
    translation = np.zeros(3)
    if name == "translated":
        translation = np.array([0.173, 0.271, 0.389])
    elif name == "rotated":
        axis = np.array([1., 2., 3.]) / math.sqrt(14.)
        x, y, z = axis
        skew = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
        angle = math.radians(37.)
        rotation = np.eye(3) * math.cos(angle) + (1 - math.cos(angle)) * np.outer(axis, axis) + math.sin(angle) * skew
    elif name != "identity":
        raise InvalidArtifact(f"unsupported rigid transform: {name}")
    return {"name": name, "rotation_matrix": rotation.tolist(), "translation_A": translation.tolist(),
            "convention": "column coordinate r_transformed = rotation @ r + translation; origin fixed"}


def transformed(atoms, transform):
    result = copy.deepcopy(atoms)
    matrix, shift = np.array(transform["rotation_matrix"]), np.array(transform["translation_A"])
    for atom in result:
        # The original source convention is ten-decimal Angstrom coordinates.
        point = np.array([round(float(v), 10) for v in atom["xyz_A"]])
        atom["xyz_A"] = [round(float(v), 10) for v in matrix @ point + shift]
    return result


def physical_xyzr(atoms):
    """Charge and cap inventories cannot affect the source-only meshing bytes."""
    return "".join(" ".join(f"{float(v):.10f}" for v in a["xyz_A"]) + f" {a['radius_A']:.6f}\n"
                   for a in sorted(atoms, key=lambda a: a["id"]))


def validate_state(state):
    required = ("source", "assembly", "microstate", "explicit_waters", "boundary_mapping",
                "charge_quality", "physical_atoms", "core_atoms", "environment_atoms",
                "core_total_charge_e", "expected_environment_charge_e", "settings")
    for key in required:
        if key not in state:
            raise InvalidArtifact(f"missing environmental state field: {key}")
    verify(state["source"])
    verify(state["boundary_mapping"])
    if state["charge_quality"].get("status") != "passed":
        raise InvalidArtifact("QM charge quality unavailable or failed")
    verify(state["charge_quality"]["receipt"])
    for key, value in PHYSICS.items():
        if state["settings"].get(key) != value:
            raise InvalidArtifact(f"physical setting differs from approved protocol: {key}")
    collections = {}
    for field in ("physical_atoms", "core_atoms", "environment_atoms"):
        atoms = state[field]
        if not atoms and field != "environment_atoms":
            raise InvalidArtifact(f"empty {field}")
        ids = set()
        for atom in atoms:
            if atom["id"] in ids or len(atom["xyz_A"]) != 3:
                raise InvalidArtifact(f"duplicate/invalid atom in {field}")
            ids.add(atom["id"])
            values = [*atom["xyz_A"], atom["radius_A"], atom.get("charge_e", float("nan"))]
            if not all(math.isfinite(v) for v in values) or atom["radius_A"] <= 0:
                raise InvalidArtifact(f"nonfinite/invalid atom in {field}")
            if field == "physical_atoms" and atom.get("kind") == "cap" and state.get("cavity_role") != "isolated_reduction":
                raise InvalidArtifact("synthetic cap in physical cavity")
        collections[field] = {a["id"]: a for a in atoms}
    physical, core, environment = (collections[k] for k in ("physical_atoms", "core_atoms", "environment_atoms"))
    if core.keys() & environment.keys():
        raise InvalidArtifact("force-field charge on a QM source atom")
    for key, expected, tolerance in (("core_atoms", "core_total_charge_e", 1e-4),
                                      ("environment_atoms", "expected_environment_charge_e", 1e-6)):
        if not math.isfinite(state[expected]) or abs(sum(a["charge_e"] for a in state[key]) - state[expected]) > tolerance:
            raise InvalidArtifact(f"charge closure failed: {key}")
    pcoords = np.array([a["xyz_A"] for a in physical.values()])
    pradii = np.array([a["radius_A"] for a in physical.values()])
    for atom in [*core.values(), *environment.values()]:
        source = physical.get(atom["id"])
        if source is None:
            if atom["id"] in environment or atom.get("kind") != "cap":
                raise InvalidArtifact("charge site absent from physical source mapping")
            # Cap point, not cap radius, participates in PB. Require the stronger
            # historical contained-sphere condition as the approved boundary audit.
            if not np.any(np.linalg.norm(pcoords - atom["xyz_A"], axis=1) + atom["radius_A"] <= pradii + 1e-6):
                raise InvalidArtifact("cap not contained in source cavity")
        elif (not np.allclose(source["xyz_A"], atom["xyz_A"], atol=1e-9, rtol=0)
              or source["radius_A"] != atom["radius_A"]):
            raise InvalidArtifact("charge site changes physical source geometry/radius")
    if environment:
        from scipy.spatial import cKDTree
        distances, _ = cKDTree([a["xyz_A"] for a in environment.values()]).query([a["xyz_A"] for a in core.values()])
        if np.any(distances < 1.0):
            raise InvalidArtifact("environment charge too close to QM atom/cap")


def prepare(state_path, output, software_manifest, level="primary", transform="identity", component="total", numerics=None):
    """Prepare only. Components all share the same physical dielectric surface."""
    state = read_json(state_path)
    validate_state(state)
    if level not in LEVELS or component not in ("total", "core", "environment", "zero"):
        raise InvalidArtifact("unsupported numerical level/component")
    if numerics is None:
        raise InvalidArtifact("explicit frozen NUMERICS.md path required")
    software = read_json(software_manifest)
    for name, expected in BINARY_HASHES.items():
        rec = software["executables"][name]
        verify(rec)
        if rec["sha256"] != expected or not os.access(rec["path"], os.X_OK):
            raise InvalidArtifact(f"unreviewed or nonexecutable solver: {name}")
    definition = transform_definition(transform)
    physical = transformed(state["physical_atoms"], definition)
    core = transformed(state["core_atoms"], definition)
    environment = transformed(state["environment_atoms"], definition)
    charges = core + environment
    if component in ("core", "environment", "zero"):
        active = {a["id"] for a in (core if component == "core" else environment if component == "environment" else [])}
        charges = [dict(a, charge_e=a["charge_e"] if a["id"] in active else 0.0) for a in charges]
    output = Path(output).resolve()
    if output.exists():
        raise InvalidArtifact(f"refusing existing preparation: {output}")
    output.mkdir(parents=True)
    (output / "physical.xyzr").write_text(physical_xyzr(physical))
    (output / "charges.pqr").write_text(pqr_lines(charges))
    setting = {**PHYSICS, **NUMERICAL_CONSTANTS, **LEVELS[level]}
    lines = ["pqr charges.pqr", f"pdie {setting['solute_dielectric']}", f"sdie {setting['solvent_dielectric']}",
             f"bulk {setting['salt_molar']}", f"temp {setting['temperature_K']}", "mesh ses",
             f"sdens {setting['grid_scale_inverse_A']}", f"srad {setting['probe_radius_A']}",
             f"tree_degree {setting['tree_degree']}", f"tree_theta {setting['tree_theta']}",
             f"tree_max_per_leaf {setting['tree_max_per_leaf']}", "precondition false", "nonpolar false",
             "outdata csv", "outdata csv_headers", "outdata timers"]
    (output / "tabi.in").write_text("\n".join(lines) + "\n")
    # Preserve the exact adapter and dependencies used by the later bridge.
    implementation = {}
    for name in ("affordable_tabi.py", "affordable_common.py", "affordable_environment.py"):
        source = Path(__file__).resolve().parent / name
        destination = output / "implementation" / name
        destination.parent.mkdir(exist_ok=True)
        shutil.copyfile(source, destination)
        implementation[name] = record(destination)
    inputs = {key: record(output / name) for key, name in (("input", "tabi.in"), ("charges", "charges.pqr"), ("physical_xyzr", "physical.xyzr"))}
    scientific = {"protocol_id": PROTOCOL, "state": record(state_path), "software": record(software_manifest),
                  "settings": setting, "transform": definition, "component": component,
                  "implementation": implementation, "numerics": record(numerics),
                  "input_hashes": {k: v["sha256"] for k, v in inputs.items()}}
    result = {**scientific, **inputs, "cache_key": cache_key(scientific), "status": "prepared_not_executed",
              "physical_boundary_hash": physical_boundary_key(physical),
              "original_physical_boundary_hash": physical_boundary_key(state["physical_atoms"]),
              "physical_atom_count": len(physical), "charge_site_count": len(charges),
              "total_charge_e": sum(a["charge_e"] for a in charges), "level": level}
    write_new(output / "tabi_manifest.json", result)
    return result


def parse_output(output_path, csv_path=None, headers_path=None):
    """Read native reaction energy; retain Coulomb only as an audit component."""
    text = Path(output_path).read_text()
    if (text.count("*** OUTPUT FOR TABI-PB RUN ***") != 1 or "GMRES error code" in text
            or re.search(r"(?i)(?<![A-Za-z])(?:[-+]?nan|[-+]?inf)(?![A-Za-z])", text)):
        raise InvalidArtifact("incomplete, ambiguous, or nonfinite TABI output")
    found = re.findall(r"GMRES completed\.\s+(\d+) iterations,\s*(" + FLOAT + r") residual\.", text)
    if len(found) != 1:
        raise InvalidArtifact("missing/ambiguous GMRES convergence")
    iterations, residual = int(found[0][0]), float(found[0][1])
    if not 0 <= iterations <= 100 or not 0 <= residual <= 1e-4:
        raise InvalidArtifact("TABI GMRES nonconvergence")
    solvation = re.findall(r"^\s*Solvation energy\s*=\s*(" + FLOAT + r")\s+kJ/mol\s*$", text, re.M)
    free = re.findall(r"^\s*Free energy\s*=\s*(" + FLOAT + r")\s+kJ/mol\s*$", text, re.M)
    if len(solvation) != 1 or len(free) != 1:
        raise InvalidArtifact("missing/ambiguous native TABI energies")
    rf, total = float(solvation[0]), float(free[0])
    coulomb = total - rf
    precision = "native_stdout_six_decimals_kJ_mol"
    csv_values = None
    if csv_path is not None or headers_path is not None:
        if csv_path is None or headers_path is None:
            raise InvalidArtifact("both native CSV and headers required")
        with Path(headers_path).open() as stream:
            header_rows = list(csv.reader(stream))
        if len(header_rows) != 1 or not header_rows[0]:
            raise InvalidArtifact("empty/truncated or ambiguous native TABI CSV header")
        headers = [x.strip() for x in header_rows[0]]
        with Path(csv_path).open() as stream:
            rows = list(csv.reader(stream))
        if len(rows) != 1 or len(rows[0]) != len(headers) or len(set(headers)) != len(headers):
            raise InvalidArtifact("invalid native TABI CSV")
        csv_values = {k: float(v.strip()) for k, v in zip(headers, rows[0]) if k}
        if not all(math.isfinite(v) for v in csv_values.values()):
            raise InvalidArtifact("nonfinite native TABI CSV")
        for key in ("solvation_energy", "coulombic_energy", "free_energy", "residual", "num_iterations"):
            if key not in csv_values:
                raise InvalidArtifact(f"missing CSV component: {key}")
        if (abs(csv_values["solvation_energy"] - rf) > 1e-6 or abs(csv_values["free_energy"] - total) > 1e-6
                or abs(csv_values["residual"] - residual) > 1e-8 or csv_values["num_iterations"] != iterations):
            raise InvalidArtifact("stdout/CSV energy or convergence mismatch")
        rf, coulomb, total = (csv_values[k] for k in ("solvation_energy", "coulombic_energy", "free_energy"))
        residual = csv_values["residual"]
        precision = "native_CSV_twelve_digits_scientific_kJ_mol"
    else:
        # Actual bundled APBS fixtures additionally expose full-precision totals.
        direct = re.findall(r"Global net COULOMBIC energy\s*=\s*(" + FLOAT + r")", text)
        reaction = re.findall(r"Global net ELEC energy\s*=\s*(" + FLOAT + r")", text)
        if len(direct) == len(reaction) == 1:
            coulomb = float(direct[0])
            refined_rf = float(reaction[0])
            if abs(refined_rf - rf) > 1e-6:
                raise InvalidArtifact("APBS/native solvation energy mismatch")
            rf = refined_rf
            precision = "archived_APBS_wrapper_full_precision_kJ_mol"
    if abs(rf + coulomb - total) > max(2e-6, abs(total) * 2e-12):
        raise InvalidArtifact("TABI component accounting failed")
    return {"reaction_field_kJ_mol": rf, "reaction_field_kcal_mol": rf / KJ_PER_KCAL,
            "coulomb_kJ_mol": coulomb,
            "audit_only_coulomb_kJ_mol": coulomb, "audit_only_free_energy_kJ_mol": total,
            "coulomb_added_to_descriptor": False, "gmres_iterations": iterations,
            "gmres_residual": residual, "energy_precision": precision, "native_csv_values": csv_values}


def _verified_manifest(path):
    manifest = read_json(path)
    if manifest["protocol_id"] != PROTOCOL:
        raise InvalidArtifact("wrong TABI protocol")
    scientific_keys = ("protocol_id", "state", "software", "settings", "transform", "component", "implementation", "numerics", "input_hashes")
    if cache_key({key: manifest[key] for key in scientific_keys}) != manifest["cache_key"]:
        raise InvalidArtifact("TABI manifest/cache key mismatch")
    if (manifest.get("level") not in LEVELS
            or manifest["settings"] != {**PHYSICS, **NUMERICAL_CONSTANTS, **LEVELS[manifest["level"]]}):
        raise InvalidArtifact("TABI settings differ from frozen numerical schedule")
    for key in ("state", "software", "input", "charges", "physical_xyzr", "numerics"):
        verify(manifest[key])
    for rec in manifest["implementation"].values():
        verify(rec)
    software = read_json(verify(manifest["software"]))
    for name, rec in software["executables"].items():
        verify(rec)
        if BINARY_HASHES.get(name) != rec["sha256"]:
            raise InvalidArtifact("unreviewed executable hash")
    return manifest, software


def mesh_bridge(manifest_path, attempt):
    """Called only by TABI's subprocess. No charges enter NanoShaper input."""
    manifest, software = _verified_manifest(manifest_path)
    attempt = Path(attempt).resolve()
    if Path.cwd().resolve() != attempt:
        raise InvalidArtifact("mesh bridge working directory mismatch")
    original = attempt / "surfaceConfiguration.prm"
    config = original.read_text()
    expected = {"Grid_scale": str(manifest["settings"]["grid_scale_inverse_A"]),
                "XYZR_FileName": "molecule.xyzr", "Surface": "ses", "Probe_Radius": "1.4"}
    entries = dict(re.findall(r"^\s*(\w+)\s*=\s*(.*?)\s*$", config, re.M))
    for key, value in expected.items():
        if key in ("Grid_scale", "Probe_Radius"):
            matches = float(entries.get(key, "nan")) == float(value)
        else:
            matches = entries.get(key) == value
        if not matches:
            raise InvalidArtifact(f"native mesher setting mismatch: {key}")
    saved = attempt / "mesh"
    saved.mkdir(exist_ok=False)
    (saved / "native_surfaceConfiguration.prm").write_text(config)
    # These control CPU use only; the physical/mathematical mesh settings stay native.
    effective = config + "Number_thread = 1\n"
    original.write_text(effective)
    (saved / "effective_surfaceConfiguration.prm").write_text(effective)
    shutil.copyfile(verify(manifest["physical_xyzr"]), attempt / "molecule.xyzr")
    shutil.copyfile(attempt / "molecule.xyzr", saved / "physical.xyzr")
    start = time.monotonic()
    with (saved / "nanoshaper.out").open("x") as log:
        completed = subprocess.run([str(verify(software["executables"]["NanoShaper"]))], cwd=attempt,
                                   stdout=log, stderr=subprocess.STDOUT, check=False)
    result = {"command": [software["executables"]["NanoShaper"]["path"]], "returncode": completed.returncode,
              "wall_seconds": time.monotonic() - start, "source_xyzr": record(saved / "physical.xyzr"),
              "native_configuration": record(saved / "native_surfaceConfiguration.prm"),
              "effective_configuration": record(saved / "effective_surfaceConfiguration.prm"),
              "output": record(saved / "nanoshaper.out"), "executable": software["executables"]["NanoShaper"],
              "mesh": {}}
    if completed.returncode == 0:
        for name in ("triangulatedSurf.vert", "triangulatedSurf.face"):
            source = attempt / name
            if source.is_file():
                shutil.copyfile(source, saved / name)
                result["mesh"][name] = record(saved / name)
    write_new(saved / "execution.json", result)
    if completed.returncode != 0 or len(result["mesh"]) != 2:
        raise InvalidArtifact("NanoShaper failed or did not produce its mesh")
    return result


def validate_native_controls(values, settings, charged_atom_count):
    """Check controls echoed by the real binary, not just its input keywords."""
    expected = {
        "mesh_density": settings["grid_scale_inverse_A"],
        "mesh_probe_radius": settings["probe_radius_A"],
        "tree_degree": settings["tree_degree"],
        "tree_theta": settings["tree_theta"],
        "tree_max_per_leaf": settings["tree_max_per_leaf"],
        "precondition": float(settings["precondition"]),
        "num_atoms": charged_atom_count,
    }
    actual = values["native_csv_values"]
    for key, wanted in expected.items():
        got = actual.get(key)
        if got is None or not math.isfinite(got) or abs(got - wanted) > 1e-12:
            raise InvalidArtifact(f"native solver control differs from manifest: {key}")
    return {"status": "passed", "native_echoed_controls": expected}


def collect(manifest_path, receipt_path):
    manifest, _ = _verified_manifest(manifest_path)
    receipt = read_json(receipt_path)
    verify(receipt["manifest"])
    if receipt["manifest"] != record(manifest_path) or receipt["cache_key"] != manifest["cache_key"]:
        raise InvalidArtifact("execution receipt belongs to another TABI task")
    if receipt["returncode"] != 0:
        raise InvalidArtifact("TABI execution failed; reaction energy unavailable")
    for rec in receipt["artifacts"].values():
        verify(rec)
    for key in ("input", "charges"):
        if receipt["artifacts"][key]["sha256"] != manifest[key]["sha256"]:
            raise InvalidArtifact(f"executed {key} differs from pinned preparation")
    bridge_path = verify(receipt["artifacts"]["mesh_receipt"])
    bridge = read_json(bridge_path)
    if bridge["returncode"] != 0:
        raise InvalidArtifact("mesher execution failed")
    for key in ("source_xyzr", "native_configuration", "effective_configuration", "output", "executable"):
        verify(bridge[key])
    if bridge["source_xyzr"]["sha256"] != manifest["physical_xyzr"]["sha256"]:
        raise InvalidArtifact("executed mesh did not use physical source cavity")
    for rec in bridge["mesh"].values():
        verify(rec)
    mesh_geometry = {}
    for name, rec in bridge["mesh"].items():
        # Preserve full original files separately; omit only two descriptive
        # MSMS header lines for the common geometry fingerprint.
        mesh_geometry[name] = verify(rec).read_text().splitlines()[2:]
    vertices = mesh_geometry.get("triangulatedSurf.vert", [])
    faces = mesh_geometry.get("triangulatedSurf.face", [])
    if not vertices or not faces:
        raise InvalidArtifact("incomplete retained MSMS mesh")
    vertex_count, face_count = int(vertices[0].split()[0]), int(faces[0].split()[0])
    if vertex_count != len(vertices) - 1 or face_count != len(faces) - 1 or min(vertex_count, face_count) <= 0:
        raise InvalidArtifact("retained mesh vertex/face count mismatch")
    for line in vertices[1:]:
        values = line.split()
        if len(values) < 6 or not all(math.isfinite(float(value)) for value in values[:6]):
            raise InvalidArtifact("invalid/nonfinite mesh coordinates or normals")
    for line in faces[1:]:
        values = line.split()
        if len(values) < 3 or not all(1 <= int(value) <= vertex_count for value in values[:3]):
            raise InvalidArtifact("mesh face index outside vertex inventory")
    values = parse_output(verify(receipt["artifacts"]["output"]), verify(receipt["artifacts"]["csv"]),
                          verify(receipt["artifacts"]["headers"]))
    charged_atom_count = sum(line.startswith(("ATOM ", "HETATM "))
                             for line in verify(manifest["charges"]).read_text().splitlines())
    controls = validate_native_controls(values, manifest["settings"], charged_atom_count)
    if values["native_csv_values"].get("num_particles") != vertex_count:
        raise InvalidArtifact("retained mesh differs from native solver particle count")
    return {"protocol_id": PROTOCOL, "status": "computed_unvalidated_reaction_field", **values,
            "source_manifest": record(manifest_path), "execution_receipt": record(receipt_path),
            "mesh": bridge["mesh"], "physical_boundary_hash": manifest["physical_boundary_hash"],
            "mesh_geometry_sha256": cache_key(mesh_geometry),
            "mesh_vertex_count": vertex_count, "mesh_face_count": face_count,
            "executed_controls": controls,
            "physical_xyzr_sha256": manifest["physical_xyzr"]["sha256"],
            "collector": record(__file__), "decision": "uncalibrated_protocol", "correction": None}


def execute(manifest_path, reuse=True):
    """One explicitly prepared task, immutable attempts, no automatic retries."""
    if not os.environ.get("SLURM_JOB_ID"):
        raise InvalidArtifact("TABI execution requires an existing SLURM allocation")
    manifest_path = Path(manifest_path).resolve()
    manifest, software = _verified_manifest(manifest_path)
    root = manifest_path.parent
    with (root / "execution.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if reuse:
            for path in sorted(root.glob("attempt_*/execution.json"), reverse=True):
                try:
                    result = collect(manifest_path, path)
                except (InvalidArtifact, KeyError):
                    continue
                return {"status": "cached_verified", "result": result, "execution_receipt": record(path)}
        count = len(list(root.glob("attempt_*"))) + 1
        attempt = root / f"attempt_{count:04d}"
        attempt.mkdir()
        for key, name in (("input", "tabi.in"), ("charges", "charges.pqr")):
            shutil.copyfile(verify(manifest[key]), attempt / name)
        bridge = attempt / "bridge_bin"
        bridge.mkdir()
        adapter = verify(manifest["implementation"]["affordable_tabi.py"])
        command = [sys.executable, str(adapter), "mesh-bridge", "--manifest", str(manifest_path), "--attempt", str(attempt)]
        (bridge / "NanoShaper").write_text("#!/bin/sh\nexec " + shlex.join(command) + "\n")
        (bridge / "NanoShaper").chmod(0o755)
        env = {**os.environ, "PATH": str(bridge) + os.pathsep + os.environ.get("PATH", ""),
               "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        invocation = ["/usr/bin/time", "-v", "-o", str(attempt / "resources.txt"),
                      str(verify(software["executables"]["tabipb"])), "tabi.in"]
        started = time.monotonic()
        with (attempt / "tabi.out").open("x") as log:
            process = subprocess.run(invocation, cwd=attempt, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
        artifacts = {}
        for key, name in (("input", "tabi.in"), ("charges", "charges.pqr"),
                          ("output", "tabi.out"), ("resources", "resources.txt"), ("csv", "output.csv"),
                          ("headers", "headers.csv"), ("mesh_receipt", "mesh/execution.json"), ("bridge", "bridge_bin/NanoShaper")):
            if (attempt / name).is_file():
                artifacts[key] = record(attempt / name)
        receipt = {"manifest": record(manifest_path), "cache_key": manifest["cache_key"], "command": invocation,
                   "returncode": process.returncode, "wall_seconds": time.monotonic() - started,
                   "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "allocated_cpus": os.environ.get("SLURM_CPUS_ON_NODE"),
                   "thread_environment": {k: env[k] for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
                   "mesh_bridge_command": command, "artifacts": artifacts, "executor": dict(LOADED_IMPLEMENTATION),
                   "executor_identity_semantics": "source_identity_captured_at_module_load",
                   "python_executable": record(sys.executable), "status": "attempt_finished"}
        receipt_path = attempt / "execution.json"
        write_new(receipt_path, receipt)
        try:
            result = collect(manifest_path, receipt_path)
        except (InvalidArtifact, KeyError) as exc:
            result = {"status": "failed", "reason": str(exc), "reaction_field_kcal_mol": None}
        write_new(attempt / "result.json", result)
        return {"status": result["status"], "result": result, "execution_receipt": record(receipt_path)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    cmd = sub.add_parser("software"); cmd.add_argument("--source-bin", required=True); cmd.add_argument("--output", required=True)
    cmd = sub.add_parser("prepare")
    for key in ("state", "output", "software-manifest", "numerics"):
        cmd.add_argument("--" + key, required=True)
    cmd.add_argument("--level", choices=LEVELS, default="primary")
    cmd.add_argument("--transform", choices=("identity", "translated", "rotated"), default="identity")
    cmd.add_argument("--component", choices=("total", "core", "environment", "zero"), default="total")
    cmd = sub.add_parser("execute"); cmd.add_argument("--manifest", required=True); cmd.add_argument("--no-reuse", action="store_true")
    cmd = sub.add_parser("collect"); cmd.add_argument("--manifest", required=True); cmd.add_argument("--receipt", required=True); cmd.add_argument("--output", required=True)
    cmd = sub.add_parser("mesh-bridge"); cmd.add_argument("--manifest", required=True); cmd.add_argument("--attempt", required=True)
    args = parser.parse_args()
    if args.operation == "software":
        result = install_software(args.source_bin, args.output)
    elif args.operation == "prepare":
        result = prepare(args.state, args.output, args.software_manifest, args.level, args.transform, args.component, args.numerics)
    elif args.operation == "execute":
        result = execute(args.manifest, not args.no_reuse)
    elif args.operation == "collect":
        result = collect(args.manifest, args.receipt); write_new(args.output, result)
    else:
        result = mesh_bridge(args.manifest, args.attempt)
    print(result.get("status", "completed"))


if __name__ == "__main__":
    main()
