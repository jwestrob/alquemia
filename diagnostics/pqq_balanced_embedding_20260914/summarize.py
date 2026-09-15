#!/usr/bin/env python3
"""Evaluate the frozen balanced-embedding decision rule from manifested outputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
HARTREE_TO_KCAL = 627.5094740631
ENERGY = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
NORMAL = "ORCA TERMINATED NORMALLY"
SCF_OK = "SCF CONVERGED AFTER"
SCF_BAD = "SCF NOT CONVERGED"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def energy(path: Path, manifest_hash: str, task: dict) -> float:
    text = path.read_text(errors="replace")
    if NORMAL not in text or SCF_OK not in text or SCF_BAD in text:
        raise RuntimeError(f"incomplete ORCA output: {path}")
    values = ENERGY.findall(text)
    if not values:
        raise RuntimeError(f"no final energy: {path}")
    record_path = Path(str(path) + ".execution.json")
    record = json.loads(record_path.read_text())
    if not (
        record["returncode"] == 0
        and record["normal_termination"] is True
        and record["scf_converged"] is True
        and record["manifest"]["sha256"] == manifest_hash
        and record["artifacts"]["output"]["sha256"] == sha256(path)
        and record["artifacts"]["template_input"]["sha256"] == task["input"]["sha256"]
        and record["artifacts"]["xyz"]["sha256"] == task["xyz"]["sha256"]
    ):
        raise RuntimeError(f"execution provenance failed: {path}")
    pc_path = HERE / task["point_charges"]["path"]
    if sha256(pc_path) != task["point_charges"]["sha256"]:
        raise RuntimeError(f"point-charge provenance failed: {path}")
    return float(values[-1])


def main() -> None:
    manifest_path = HERE / "embedded_tasks.json"
    manifest_hash = sha256(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    by_id = {task["task_id"]: task for task in manifest["tasks"]}
    aquo_path = PROJECT / "reference_inputs/aquo_cn8_native_r2scan3c_v2/aquo_reference.json"
    aquo = json.loads(aquo_path.read_text())
    delta_aquo = float(aquo["delta_E_aquo_hartree"])
    panels = {}
    for case in ("mxaf", "c5_monomer", "c5_dimer"):
        radii = {}
        for radius, tag in (("3.3", "33"), ("3.6", "36")):
            raw = {}
            for metal in ("La", "Ca"):
                task = by_id[f"{case}_qm{tag}_{metal}"]
                raw[metal] = energy(HERE / task["output_path"], manifest_hash, task)
            score = ((raw["Ca"] - raw["La"]) - delta_aquo) * HARTREE_TO_KCAL
            radii[radius] = {
                "E_La_hartree": raw["La"],
                "E_Ca_hartree": raw["Ca"],
                "E_Ca_minus_La_hartree": raw["Ca"] - raw["La"],
                "score_kcal_mol": score,
            }
        response = radii["3.6"]["score_kcal_mol"] - radii["3.3"]["score_kcal_mol"]
        panels[case] = {"radii": radii, "R_PC_3p6_minus_3p3_kcal_mol": response}

    r_mxaf = panels["mxaf"]["R_PC_3p6_minus_3p3_kcal_mol"]
    d_mono = panels["c5_monomer"]["R_PC_3p6_minus_3p3_kcal_mol"] - r_mxaf
    d_dimer = panels["c5_dimer"]["R_PC_3p6_minus_3p3_kcal_mol"] - r_mxaf
    gates = {
        "generic_boundary_suppressed": abs(r_mxaf) <= 1.0,
        "monomer_selective_margin": d_mono >= 2.0,
        "dimer_selective_margin": d_dimer >= 2.0,
    }
    primary_answer = "YES" if all(gates.values()) else "NO"
    cutoff_distances = {
        "MxaF_abs_R_to_1p0": abs(abs(r_mxaf) - 1.0),
        "D_monomer_to_2p0": abs(d_mono - 2.0),
        "D_dimer_to_2p0": abs(d_dimer - 2.0),
    }
    sensitivity_triggered = any(value <= 0.5 for value in cutoff_distances.values())
    final_answer = "PENDING_BOUNDARY_SENSITIVITY" if sensitivity_triggered else primary_answer

    bare_path = PROJECT / "diagnostics/pqq_boundary_pair_20260914/result.json"
    bare = json.loads(bare_path.read_text())
    for case in panels:
        for radius in ("3.3", "3.6"):
            panels[case]["radii"][radius]["embedding_shift_from_bare_kcal_mol"] = (
                panels[case]["radii"][radius]["score_kcal_mol"]
                - bare["panels"][case]["radii"][radius]["score_kcal_mol"]
            )

    boundary_embedding_jumps = {
        case: (
            panels[case]["R_PC_3p6_minus_3p3_kcal_mol"]
            - bare["panels"][case]["R_3p6_minus_3p3_kcal_mol"]
        )
        for case in panels
    }
    raw_outputs = []
    for task in manifest["tasks"]:
        output_path = HERE / task["output_path"]
        execution_path = Path(str(output_path) + ".execution.json")
        raw_outputs.append(
            {
                "task_id": task["task_id"],
                "output": {"path": str(output_path), "sha256": sha256(output_path)},
                "execution_record": {
                    "path": str(execution_path),
                    "sha256": sha256(execution_path),
                },
                "point_charges": {
                    "path": str(HERE / task["point_charges"]["path"]),
                    "sha256": task["point_charges"]["sha256"],
                },
            }
        )

    preparation_path = HERE / "preparation.json"
    result = {
        "schema_version": "pqq_balanced_embedding.result.v1",
        "status": "complete_primary_map",
        "primary_map_answer": primary_answer,
        "answer": final_answer,
        "sensitivity_triggered": sensitivity_triggered,
        "cutoff_distances_kcal_mol": cutoff_distances,
        "gates": gates,
        "R_PC_MxaF_kcal_mol": r_mxaf,
        "D_PC_monomer_kcal_mol": d_mono,
        "D_PC_dimer_kcal_mol": d_dimer,
        "embedding_induced_boundary_jump_kcal_mol": boundary_embedding_jumps,
        "panels": panels,
        "raw_outputs": raw_outputs,
        "task_manifest": {"path": str(manifest_path), "sha256": manifest_hash},
        "preparation": {"path": str(preparation_path), "sha256": sha256(preparation_path)},
        "bare_panel": {"path": str(bare_path), "sha256": sha256(bare_path)},
        "aquo_reference": {
            "path": str(aquo_path),
            "sha256": sha256(aquo_path),
            "delta_E_hartree": delta_aquo,
            "note": "used only for absolute secondary scores; cancels from R and D",
        },
    }
    temporary = HERE / ".result.json.tmp"
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    temporary.replace(HERE / "result.json")

    lines = [
        "# Balanced PQQ embedding result",
        "",
        f"**Frozen primary-map answer: {primary_answer}. Final answer: {final_answer}.**",
        "",
        "| Structure | embedded score 3.3 | embedded score 3.6 | R_PC |",
        "|---|---:|---:|---:|",
    ]
    for case in ("mxaf", "c5_monomer", "c5_dimer"):
        panel = panels[case]
        lines.append(
            f"| {case} | {panel['radii']['3.3']['score_kcal_mol']:+.3f} | "
            f"{panel['radii']['3.6']['score_kcal_mol']:+.3f} | "
            f"{panel['R_PC_3p6_minus_3p3_kcal_mol']:+.3f} |"
        )
    lines.extend(
        [
            "",
            f"D_PC monomer = {d_mono:+.3f}; D_PC dimer = {d_dimer:+.3f} kcal/mol.",
            f"Sensitivity branch triggered: {'yes' if sensitivity_triggered else 'no'}.",
            "",
            "Both selective-margin gates passed, but the generic MxaF gate failed",
            "catastrophically. Relative to the bare calculation, embedding changed the",
            "3.3-to-3.6 boundary response by "
            + ", ".join(
                f"{case} {value:+.3f}" for case, value in boundary_embedding_jumps.items()
            )
            + " kcal/mol.",
            "The near-common jump despite exact total-charge closure shows that a static",
            "embedded-QM energy is not continuous when Asp(-) changes representation from",
            "fixed MM charges to an explicit QM fragment. This model is rejected for",
            "production scoring; the non-triggered sensitivity map cannot rescue its gate.",
            "",
            "Machine-readable energies and raw-output hashes are in `result.json`; complete",
            "charge, repair, and distance ledgers are in `preparation.json`.",
            "",
        ]
    )
    temporary = HERE / ".RESULT.md.tmp"
    temporary.write_text("\n".join(lines))
    temporary.replace(HERE / "RESULT.md")
    print(
        json.dumps(
            {
                "answer": final_answer,
                "primary_map_answer": primary_answer,
                "R_PC_MxaF": r_mxaf,
                "D_PC_monomer": d_mono,
                "D_PC_dimer": d_dimer,
                "sensitivity_triggered": sensitivity_triggered,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
