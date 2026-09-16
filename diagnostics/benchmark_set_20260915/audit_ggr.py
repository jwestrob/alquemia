"""Read-only bookkeeping audit of the four completed GGR endpoint outputs.

No calculation, altered score, classification or causal energy decomposition.
Run from any directory with explicit --result and --output paths.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from affordable_common import HA_TO_KCAL, record, verify, write_new


def require(condition, message):
    if not condition:
        raise ValueError(message)


def xyz(artifact):
    lines = verify(artifact).read_text().splitlines()
    atoms = [(p[0], *map(float, p[1:])) for p in map(str.split, lines[2:]) if p]
    require(len(atoms) == int(lines[0]), "XYZ atom count")
    require(all(len(a) == 4 for a in atoms), "XYZ shape")
    return atoms


def components(endpoint):
    text = verify(endpoint["output"]).read_text()
    receipt = json.loads(verify(endpoint["receipt"]).read_text())
    require(endpoint["status"] == "complete", "endpoint incomplete")
    require(receipt["normal_termination"] and receipt["scf_converged"], "receipt failure")
    require("ORCA TERMINATED NORMALLY" in text and "SCF CONVERGED" in text,
            "output not converged/terminated")
    patterns = {
        "scf_total": r"^Total Energy\s*:\s*([-+\d.]+)\s+Eh",
        "cpcm_printed": r"^CPCM Dielectric\s*:\s*([-+\d.]+)\s+Eh",
        "dispersion": r"^Dispersion correction\s+([-+\d.]+)\s*$",
        "gcp": r"^gCP correction\s+([-+\d.]+)\s*$",
        "final": r"^FINAL SINGLE POINT ENERGY\s+([-+\d.]+)\s*$",
    }
    values = {}
    for name, pattern in patterns.items():
        matches = re.findall(pattern, text, re.M)
        require(len(matches) == 1, f"ambiguous/missing {name}")
        values[name] = float(matches[0])
    require(values["final"] == endpoint["energy_hartree"], "recorded energy mismatch")
    residual = values["final"] - sum(values[k] for k in ("scf_total", "dispersion", "gcp"))
    # Dispersion/gCP are printed to nine decimal places in these archived files.
    require(abs(residual) < 2e-9, "printed component closure exceeds precision")
    values["scf_minus_cpcm_bookkeeping_only"] = values["scf_total"] - values["cpcm_printed"]
    return {"hartree": values, "closure_residual_hartree": residual,
            "output": endpoint["output"], "receipt": endpoint["receipt"]}


def audit(result_path):
    result = json.loads(result_path.read_text())
    rows = {x["lane"]: x for x in result["scores"] if x["case"] == "ggr_1glg_GGR"}
    require(set(rows) == {"baseline", "repaired"}, "expected both GGR representations")
    manifests = {lane: json.loads(verify(row["source_manifest"]).read_text())
                 for lane, row in rows.items()}
    old, new = manifests["baseline"], manifests["repaired"]
    for key in ("source_structure", "protonation_manifest", "coordination"):
        require(old[key] == new[key], f"changed {key}")
    require(old["charge_ledger"]["fragment_formal_charge_sum"] ==
            sum(f["formal_charge"] for f in new["charge_ledger"]) == -3, "ligand charge")
    for metal in ("Ca", "La"):
        a, b = old["charge_ledger"]["clusters"][metal], new["outputs"][metal]
        require(a["total_charge"] == b["charge"] and a["multiplicity"] == b["multiplicity"] == 1,
                "charge/multiplicity changed")
    verify(new["source_structure"])
    verify(new["protonation_manifest"])
    differences = {}
    for metal in ("Ca", "La"):
        original = xyz(old["outputs"][metal]["xyz"])
        repaired = xyz(new["outputs"][metal]["xyz"])
        a, b = Counter(original), Counter(repaired)
        # v2 writes caps to six decimals; v3 preserves ten. Compare at the
        # original recorded precision without claiming bit-identical caps.
        rounded = lambda atom: (atom[0], *(round(v, 6) for v in atom[1:]))
        ar, br = Counter(map(rounded, original)), Counter(map(rounded, repaired))
        matching = [(u, v) for u in original for v in repaired if rounded(u) == rounded(v)]
        differences[metal] = {
            "original_atoms": len(original), "repaired_atoms": len(repaired),
            "exactly_unchanged_atoms": sum((a & b).values()),
            "unchanged_at_original_six_decimal_precision": sum((ar & br).values()),
            "maximum_matched_coordinate_difference_A": max(math.dist(u[1:], v[1:]) for u, v in matching),
            "removed_atoms_element_xyz_A_at_original_precision": list((ar - br).elements()),
            "added_atoms_element_xyz_A_at_original_precision": list((br - ar).elements()),
        }
        require(len(original) == 50 and len(repaired) == 52, "unexpected inventory")
        require(sum((ar & br).values()) == 49, "unexpected coordinate changes")
        require(all(t[0] == "H" for t in (a - b)), "original heavy atom changed")
    for manifest in manifests.values():
        ca = xyz(manifest["outputs"]["Ca"]["xyz"])
        la = xyz(manifest["outputs"]["La"]["xyz"])
        require([a[1:] for a in ca] == [a[1:] for a in la], "unpaired coordinates")
        require([(a[0], b[0]) for a, b in zip(ca, la) if a[0] != b[0]] == [("Ca", "La")],
                "unexpected paired atom substitution")
    endpoints = {lane: {metal: components(ep) for metal, ep in row["endpoints"].items()}
                 for lane, row in rows.items()}
    contrasts = {lane: {key: (pair["Ca"]["hartree"][key] - pair["La"]["hartree"][key])
                            for key in pair["Ca"]["hartree"]}
                 for lane, pair in endpoints.items()}
    shifts = {key: (contrasts["repaired"][key] - contrasts["baseline"][key]) * HA_TO_KCAL
              for key in contrasts["baseline"]}
    require(abs(shifts["final"] - (rows["repaired"]["score"]["S_kcal_mol"] -
                                    rows["baseline"]["score"]["S_kcal_mol"])) < 1e-7,
            "score difference algebra")
    return {
        "schema_version": "alquemia.ggr_archived_audit.v1", "status": "verified",
        "scientific_endpoint_evaluations_launched": 0, "implementation": record(__file__),
        "source_result": record(result_path), "hartree_to_kcal_mol": HA_TO_KCAL,
        "manifests": {lane: row["source_manifest"] for lane, row in rows.items()},
        "protocols": {lane: row["protocol_id"] for lane, row in rows.items()},
        "coordinate_comparison": differences, "endpoints": endpoints,
        "reported_S_kcal_mol": {lane: row["score"]["S_kcal_mol"] for lane, row in rows.items()},
        "repair_effect_on_Ca_minus_La_kcal_mol": shifts,
        "previously_reported_all_six_repair_shifts": result["repair_differences"],
        "interpretation": "SCF minus printed CPCM is bookkeeping on a polarized density, not a vacuum energy. Components do not isolate causal mechanisms. Repaired v3 has no calibrated zero/bands; do not remove components or refit a threshold.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    write_new(args.output, audit(args.result.resolve()))
    print(args.output)
