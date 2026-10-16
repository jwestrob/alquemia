"""Saved Hans source-gradient profiles; numerical execution requires Slurm.

No molecular engine, geometry proposal, optimization, or finite energy is run.
Cartesian source differences are diagnostic vectors, not admissible motions.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import sys
import time
import traceback

ROOT = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs")
SOURCES = ("Hans8DQ2", "Hans8FNR")
METALS = ("La", "Dy")


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def write_new(path, data):
    # Serialize before creating the immutable output, preserving failed attempts.
    text = json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n"
    with Path(path).open("x") as handle:
        handle.write(text)


def profile(provenance_path, output):
    # Deliberately after the Slurm guard in main(). No scientific imports on login.
    import numpy as np
    sys.path.insert(0, str(ROOT / "scripts"))
    from affordable_common import read_json, verify, InvalidArtifact
    from metal_environment_force_assembly import assemble, target_source_id

    dependencies = [pin(module.__file__) for module in list(sys.modules.values())
                    if getattr(module, "__file__", None)
                    and Path(module.__file__).resolve().parent == ROOT / "scripts"]

    consumed = {}

    def read_pinned(record):
        path = verify(record)
        consumed[str(path.resolve())] = pin(path)
        return read_json(path)

    def require(condition, message):
        if not condition:
            raise InvalidArtifact(message)

    def index_unique(rows):
        result = {row["id"]: row for row in rows}
        require(len(result) == len(rows), "Duplicate physical/source atom ID")
        return result

    def stats(array):
        array = np.asarray(array, dtype=float)
        require(array.ndim == 2 and array.shape[1] == 3 and np.isfinite(array).all(),
                "Nonfinite or malformed gradient")
        return {"atom_count": len(array), "gradient_norm_kcal_mol_A": float(np.linalg.norm(array)),
                "force_norm_kcal_mol_A": float(np.linalg.norm(array)),
                "gradient_cartesian_RMS_kcal_mol_A": float(np.sqrt(np.mean(array**2))),
                "net_gradient_kcal_mol_A": array.sum(axis=0).tolist(),
                "net_force_kcal_mol_A": (-array.sum(axis=0)).tolist()}

    provenance = read_json(provenance_path)
    final = read_pinned(provenance["result"])
    require(final["complete_cells"] == final["declared_cells"] == 6,
            "Expected the complete six-cell final comparison")
    rows = {(row["source"], row["metal"]): row for row in final["rows"]}
    require(len(rows) == len(final["rows"]) == 6, "Duplicate/incomplete final cells")
    require(set(rows) == {(s, m) for s in (*SOURCES, "Mex8FNS") for m in METALS},
            "Unexpected final source/metal inventory")
    assemblies, source_data, checks = {}, {}, {}
    for source in SOURCES:
        inputs_record = provenance["source_states"][source]["inputs"]
        inputs = read_pinned(inputs_record)
        exported = read_pinned(inputs["export"])
        metadata = index_unique(read_pinned(exported["atoms"]))
        physical = index_unique(read_pinned(inputs["physical_sources"]["A"]))
        bonds = read_pinned(exported["bonds"])
        core = read_pinned(inputs["configurations"]["A"]["core_mapping"])
        target = [target_source_id(a) for a in core if a["kind"] == "metal"]
        require(len(target) == 1, "Expected exactly one physical target")
        spectators = [row["id"] for row in inputs["spectator_occupancy"]]
        protein = {sid for sid, a in metadata.items() if a["kind"] == "protein_source"}
        waters = {sid for sid, a in metadata.items() if a["kind"] == "retained_crystal_water"}
        require(protein | waters == set(metadata) and not protein & waters,
                "Unclassified native protein/water export atoms")
        require(set(physical) == set(metadata), "Repaired physical/export inventories differ")
        grouped_ids = {"protein": sorted(protein), "waters": sorted(waters),
                       "target": target, "spectators": spectators}
        source_data[source] = {"inputs": inputs, "metadata": metadata,
                               "protein": protein, "groups": grouped_ids, "bonds": bonds}
        for metal in METALS:
            row = rows[(source, metal)]
            require(row["status"] == row["classical_status"] == "complete",
                    "Unavailable electronic/classical cell")
            require(row["inputs"] == inputs_record and row["task_id"] == metal + "_A",
                    "Final cell does not bind the exact repaired origin")
            collection = read_pinned(row["collection"])
            manifest = read_pinned(row["manifest"])
            classical = read_pinned(row["classical"])
            require(classical["metal"] == metal and classical["inputs"] == inputs_record,
                    "Classical metal/input binding differs")
            current = collection["rows"][row["task_id"]]
            require(current["status"] == "complete" and current["energy_hartree"] == row["energy_hartree"],
                    "Final comparison/collection energy differs")
            receipt = read_pinned(current["receipt"])
            require(receipt["manifest"] == collection["manifest"] == row["manifest"],
                    "Collection manifest differs")
            task = next(t for t in manifest["tasks"] if t["task_id"] == row["task_id"])
            require(current["electronic_state"] == row["electronic_state"] == task["electronic_state"],
                    "Accepted electronic-state records differ")
            require(task["metal"] == metal and task["charge"] ==
                    inputs["configurations"]["A"]["endpoints"][metal]["charge"],
                    "Task metal or charge differs from prepared endpoint")
            source_xyz = inputs["configurations"]["A"]["endpoints"][metal]["xyz"]
            for coordinate_record in (task["xyz"], source_xyz):
                coordinate_path = verify(coordinate_record)
                consumed[str(coordinate_path.resolve())] = pin(coordinate_path)
            require(task["xyz"]["sha256"] == source_xyz["sha256"],
                    "Task geometry bytes differ from repaired origin")
            for key in ("engrad", "pointcharge_gradient"):
                path = verify(current[key])
                consumed[str(path.resolve())] = pin(path)
            assembled = assemble(inputs_record["path"], "A", metal,
                                 row["collection"]["path"], row["task_id"], row["classical"]["path"])
            atoms = index_unique(assembled["physical_atoms"])
            ids = [a["id"] for a in assembled["physical_atoms"]]
            require(set(ids) == set(metadata) | set(target) | set(spectators),
                    "Physical categories do not close")
            for sid, atom in physical.items():
                require(np.allclose(atom["xyz_A"], atoms[sid]["xyz_A"], atol=1e-8, rtol=0),
                        "Assembly differs from current repaired physical coordinates: " + sid)
            total = np.asarray(assembled["total_gradient_kcal_mol_A"])
            require(total.shape == (len(ids), 3) and np.isfinite(total).all(), "Invalid full gradient")
            require(np.array_equal(-total, assembled["total_force_kcal_mol_A"]), "Force sign differs")
            checks[source + "/" + metal] = {
                "source_input_and_metal_binding": True, "current_repaired_coordinates": True,
                "electronic_state_bound_to_final_collection_and_task": row["electronic_state"],
                "native_gradient_coordinate_energy_collection_checks": "passed by existing assemble()",
                "exact_classical_ID_coordinate_units_exclusions": "passed by existing assemble()",
                "physical_atom_count": len(ids), "force_is_negative_gradient": True,
                "mapping_translation_conservation_residual_kcal_mol_A":
                    assembled["electronic"]["mapping_translation_conservation_residual_kcal_mol_A"],
                "total_translation_gradient_sum_kcal_mol_A": total.sum(axis=0).tolist(),
                "full_hybrid_derivative_qualified": False}
            write_new(output / (source + "_" + metal + "_ASSEMBLED.json"), assembled)
            assemblies[(source, metal)] = assembled
        la, dy = (assemblies[(source, m)] for m in METALS)
        require([a["id"] for a in la["physical_atoms"]] == [a["id"] for a in dy["physical_atoms"]],
                "Within-source target substitutions have different physical ordering")
        require(np.allclose([a["xyz_A"] for a in la["physical_atoms"]],
                            [a["xyz_A"] for a in dy["physical_atoms"]], atol=1e-8, rtol=0),
                "Within-source target substitutions have different coordinates")

    common = sorted(source_data[SOURCES[0]]["protein"] & source_data[SOURCES[1]]["protein"])
    require(len(common) == 1665 and all(set(common) == source_data[s]["protein"] for s in SOURCES),
            "Expected exact common 1665-atom Hans protein inventory")
    for sid in common:
        a, b = (source_data[s]["metadata"][sid] for s in SOURCES)
        require(all(a[key] == b[key] for key in ("element", "resname", "charge_e")),
                "Protein identity or protonation/charge differs: " + sid)
    heavy = [i for i, sid in enumerate(common)
             if source_data[SOURCES[0]]["metadata"][sid]["element"] != "H"]
    protein_coordinates = {}
    for source in SOURCES:
        atoms = index_unique(assemblies[(source, "La")]["physical_atoms"])
        protein_coordinates[source] = np.asarray([atoms[sid]["xyz_A"] for sid in common])
    reference, mobile = (protein_coordinates[s] for s in SOURCES)
    ref_center, mob_center = reference[heavy].mean(axis=0), mobile[heavy].mean(axis=0)
    left, singular, right_t = np.linalg.svd((mobile[heavy] - mob_center).T @ (reference[heavy] - ref_center))
    require(singular[-1] > 1e-10, "Degenerate heavy-atom alignment")
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(left @ right_t)
    rotation = left @ correction @ right_t
    require(np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-12, rtol=0)
            and abs(np.linalg.det(rotation) - 1) < 1e-12, "Improper/nonorthogonal alignment")
    aligned = (mobile - mob_center) @ rotation + ref_center
    displacement = aligned - reference
    norm = float(np.linalg.norm(displacement))
    require(norm > 0 and np.isfinite(norm), "Missing source-difference direction")
    direction = displacement / norm
    result = {
        "protocol": "lanm_saved_Hans_source_gradient_profile_v1", "new_molecular_calls": 0,
        "final_comparison": provenance["result"], "assembly_checks": checks,
        "units": {"coordinates_displacement": "angstrom", "gradients_forces": "kcal/mol/angstrom",
                  "normalized_direction_projection": "kcal/mol/angstrom"},
        "convention": "gradient=dE/dR; force=-gradient; differential=Dy-La; common=(La+Dy)/2",
        "scope": "Source-difference tangent diagnostics only; no finite motion or energy prediction",
        "alignment": {"reference": SOURCES[0], "mobile": SOURCES[1],
                      "rule": "equal weight all exactly mapped common protein heavy atoms",
                      "anchor_source_ids": [common[i] for i in heavy],
                      "row_vector_rotation": rotation.tolist(), "reference_center_A": ref_center.tolist(),
                      "mobile_center_A": mob_center.tolist(),
                      "aligned_heavy_RMSD_A": float(np.sqrt(np.mean(np.sum(displacement[heavy]**2, axis=1)))),
                      "all_protein_displacement_norm_A": norm,
                      "water_target_spectator_atoms_used_for_alignment": False},
        "matched_protein": [{"id": sid, "element": source_data[SOURCES[0]]["metadata"][sid]["element"],
                             "resname": source_data[SOURCES[0]]["metadata"][sid]["resname"],
                             "DQ2_xyz_A": reference[i].tolist(), "aligned_FNR_xyz_A": aligned[i].tolist(),
                             "FNR_minus_DQ2_displacement_A": displacement[i].tolist()}
                            for i, sid in enumerate(common)],
        "source_profiles": {}, "finite_matrix": None, "relaxation_energy": None,
        "classifier": None, "physical_motion_admissible": False, "full_hybrid_qualified": False,
        "limitations": ["Cartesian interpolation is not a covalent-preserving path",
                        "Projection is local load, not finite work or a relaxed-state energy",
                        "Waters and target/spectators are not part of the protein-difference direction",
                        "Waters are grouped only within each source; no cross-source water matching",
                        "Frozen-f electronic state, fixed occupancy and absent bulk solvent remain",
                        "No Hessian, state populations, force correction or new derivative qualification"]}
    table = ["source\tresidue\tresname\tatom_count\tdisplacement_norm_A\tLa_local_gradient_projection\tDy_local_gradient_projection\tcommon_local_gradient_projection\tdifferential_local_gradient_projection\tdifferential_local_force_projection\tLa_additive_global_contribution\tDy_additive_global_contribution\tcommon_additive_global_contribution\tdifferential_additive_global_contribution"]
    residues = sorted({sid.rsplit("/", 1)[0] for sid in common},
                      key=lambda sid: (sid.split("/")[0], int(sid.split("/")[1]), sid))
    for source in SOURCES:
        la, dy = (assemblies[(source, m)] for m in METALS)
        ids = [a["id"] for a in la["physical_atoms"]]
        ix = {sid: i for i, sid in enumerate(ids)}
        ordered = [ix[sid] for sid in common]
        transform = np.eye(3) if source == SOURCES[0] else rotation
        # Each source's direction points toward the other source in the same aligned frame.
        own_direction = direction if source == SOURCES[0] else -direction
        components = {}
        for component, field in (("electronic", "gradient_kcal_mol_A"),
                                 ("classical", "classical_gradient_kcal_mol_A"),
                                 ("total", "total_gradient_kcal_mol_A")):
            values = {m: np.asarray(assemblies[(source, m)]["electronic"][field]
                                   if component == "electronic" else assemblies[(source, m)][field])
                      for m in METALS}
            values["Dy_minus_La"] = values["Dy"] - values["La"]
            values["common"] = (values["Dy"] + values["La"]) / 2
            group_stats = {group: {m: stats(g[[ix[sid] for sid in group_ids]]) for m, g in values.items()}
                           for group, group_ids in source_data[source]["groups"].items()}
            projections = {m: float(np.sum((g[ordered] @ transform) * own_direction))
                           for m, g in values.items()}

            def local_load(indices):
                local_norm = float(np.linalg.norm(displacement[indices]))
                loads = {}
                for label, gradient in values.items():
                    part = gradient[[ordered[i] for i in indices]] @ transform
                    additive = float(np.sum(part * own_direction[indices]))
                    local = additive * norm / local_norm if local_norm > 0 else None
                    loads[label] = dict(stats(part),
                        local_gradient_projection_kcal_mol_A=local,
                        local_force_projection_kcal_mol_A=-local if local is not None else None,
                        additive_global_gradient_contribution_kcal_mol_A=additive,
                        additive_global_force_contribution_kcal_mol_A=-additive)
                return {"atom_count": len(indices), "source_displacement_norm_A": local_norm,
                        "local_direction_available": local_norm > 0, "loads": loads}

            segments = {
                "whole_protein": list(range(len(common))),
                "EF3_A83_94": [i for i, sid in enumerate(common)
                               if sid.split("/")[0] == "A" and 83 <= int(sid.split("/")[1]) <= 94],
                "remaining_protein": [i for i, sid in enumerate(common)
                                      if not (sid.split("/")[0] == "A" and 83 <= int(sid.split("/")[1]) <= 94)],
                "protein_heavy": heavy,
                "protein_hydrogens": [i for i, sid in enumerate(common)
                    if source_data[source]["metadata"][sid]["element"] == "H"]}
            per_residue = []
            for residue in residues:
                indices = [i for i, sid in enumerate(common) if sid.rsplit("/", 1)[0] == residue]
                per_residue.append({"residue": residue,
                    "resname": source_data[source]["metadata"][common[indices[0]]]["resname"],
                    **local_load(indices)})
            for metal in values:
                require(abs(sum(r["loads"][metal]["additive_global_gradient_contribution_kcal_mol_A"] for r in per_residue)
                            - projections[metal]) < 1e-8, "Per-residue projection sum mismatch")
            components[component] = {"groups": group_stats,
                "gradient_projection_kcal_mol_A": projections,
                "force_projection_kcal_mol_A": {m: -v for m, v in projections.items()},
                "protein_segments": {name: local_load(indices) for name, indices in segments.items()},
                "per_residue": per_residue}
        for row in components["total"]["per_residue"]:
            loads = row["loads"]
            table.append("\t".join(map(str, [source, row["residue"], row["resname"], row["atom_count"],
                row["source_displacement_norm_A"], *[loads[m]["local_gradient_projection_kcal_mol_A"]
                for m in (*METALS, "common", "Dy_minus_La")],
                loads["Dy_minus_La"]["local_force_projection_kcal_mol_A"],
                *[loads[m]["additive_global_gradient_contribution_kcal_mol_A"]
                  for m in (*METALS, "common", "Dy_minus_La")]])))
        # Native graph inventory only: no torsion angles, candidate coordinates or force-based selection.
        graph = {sid: set() for sid in common}
        for a, b in source_data[source]["bonds"]:
            if a in graph and b in graph:
                graph[a].add(b)
                graph[b].add(a)
        axes = []
        for name, a, b in (("Asp85_CA_CB", "A/85//CA", "A/85//CB"),
                           ("Glu91_CB_CG", "A/91//CB", "A/91//CG"),
                           ("N83_CA83", "A/83//N", "A/83//CA")):
            require(a in graph and b in graph[a], "Declared archived axis is not a native protein bond")
            seen, pending = {b}, [b]
            while pending:
                current = pending.pop()
                for neighbor in graph[current]:
                    if {current, neighbor} == {a, b} or neighbor in seen:
                        continue
                    seen.add(neighbor)
                    pending.append(neighbor)
            axes.append({"name": name, "axis_source_ids": [a, b], "native_bond_present": True,
                         "is_graph_bridge": a not in seen,
                         "distal_component_source_ids": sorted(seen) if a not in seen else None,
                         "motion_generated_by_this_profile": False,
                         "finite_response_generated_by_this_profile": False})
        result["source_profiles"][source] = {"direction": "toward other Hans source after rigid alignment",
            "normalization": "Local loads use each subset's Euclidean displacement norm; additive contributions use the whole-protein norm. Only additive contributions sum to total.",
            "group_vector_frame": "native source frame; norms are rotation invariant",
            "protein_projection_frame": "aligned DQ2 frame; FNR gradients rotated with coordinates",
            "components": components, "archived_axis_graph_inventory": axes}
    for record in [*consumed.values(), *dependencies]:
        verify(record)
    result["consumed_pins"] = list(consumed.values())
    result["runtime_dependencies"] = dependencies
    result["runtime_versions"] = {"python": sys.version, "numpy": np.__version__}
    write_new(output / "PROFILE.json", result)
    with (output / "RESIDUES.tsv").open("x") as handle:
        handle.write("\n".join(table) + "\n")
    return {"result": pin(output / "PROFILE.json"), "residue_table": pin(output / "RESIDUES.tsv"),
            "assembled_cells": 4, "common_protein_atoms": len(common), "new_molecular_calls": 0}


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Saved-force numerical profiling requires Slurm")
    if int(os.environ.get("SLURM_CPUS_PER_TASK", "1")) != 1 or int(os.environ.get("SLURM_NTASKS", "1")) != 1:
        raise RuntimeError("Expected the declared one-task, one-CPU allocation")
    for path in (args.plan, args.provenance):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.output.mkdir(parents=True, exist_ok=False)
    receipt = {"status": "running", "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "host": socket.gethostname(), "plan": pin(args.plan), "provenance": pin(args.provenance),
        "implementation": pin(__file__), "started_epoch": time.time(), "allocated_CPUs": 1,
        "new_molecular_calls": 0}
    write_new(args.output / "START.json", receipt)
    status = 1
    try:
        receipt.update(profile(args.provenance, args.output))
        for name in ("plan", "provenance", "implementation"):
            if pin(receipt[name]["path"]) != receipt[name]:
                raise RuntimeError("Changed input/implementation during execution: " + name)
        receipt["status"] = "complete"
        status = 0
    except Exception as error:
        receipt.update(status="failed", result_accepted=False, error=str(error), traceback=traceback.format_exc())
    finally:
        receipt["finished_epoch"] = time.time()
        receipt["elapsed_seconds"] = receipt["finished_epoch"] - receipt["started_epoch"]
        receipt["measured_worker_CPU_allocation_seconds"] = receipt["elapsed_seconds"]
        receipt["terminal_scheduler_accounting"] = "collect separately after job termination"
        write_new(args.output / "RECEIPT.json", receipt)
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    arguments.plan = arguments.plan.resolve()
    arguments.provenance = arguments.provenance.resolve()
    arguments.output = arguments.output.resolve()
    raise SystemExit(main(arguments))
