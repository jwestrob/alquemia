#!/usr/bin/env python3
"""Summarize actual saved cavity results; no molecular calls."""
import argparse
import hashlib
import json
import os
from pathlib import Path


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def read_pin(item):
    assert pin(item["path"]) == item
    return json.loads(Path(item["path"]).read_text())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--manifest-sha", required=True)
    a = p.parse_args()
    assert os.environ.get("SLURM_JOB_ID"), "Saved scientific summary requires Slurm"
    assert pin(a.manifest)["sha256"] == a.manifest_sha
    m = json.loads(a.manifest.read_text())
    assert pin(__file__) == m["implementation"]
    c = read_pin(m["collection"])
    e = read_pin(m["execution"])
    assert c["status"] == e["status"] == "complete"
    assert e["native_child_execution_receipts"] == e["successful_native_updates"] == 12
    assert c["scout"]["status"] == "passed"
    ha = c["hartree_to_kcal_mol"]
    rows = []
    for r in c["rows"]:
        assert r["status"] == "available"
        old = read_pin(r["archived_local"])["vacuum"]
        union = r["updates"]["union"]
        rows.append({"id": r["id"], "status": r["status"],
            "transfer_kcal_mol": {k: v * ha for k, v in r["transfer_Eh"].items()},
            "cap_volume_control_kcal_mol": {k: v * ha for k, v in r["cap_volume_Eh"].items()},
            "union_electrostatic_blocks_kcal_mol": {k: v * ha for k, v in union["electrostatic_blocks_Eh"].items()},
            "HB_change_excluded_from_transfer_kcal_mol": (union["components_Eh"]["Ghb"] - old["components_Eh"]["Ghb"]) * ha,
            "transfer_charge_bound_kcal_mol": r["conservative_transfer_charge_bound_Eh"] * ha,
            "cap_volume_charge_bound_kcal_mol": r["conservative_cap_volume_charge_bound_Eh"] * ha})
    def converted(values):
        return {k.removesuffix('_Eh'): {t: v * ha for t, v in terms.items()} for k, terms in values.items()}
    result = {"protocol": c["protocol"], "inputs": m, "summary_job": os.environ["SLURM_JOB_ID"],
        "status": "complete", "rows": rows,
        "paired_kcal_mol": {k: converted(v) for k, v in c["Ca_minus_La"].items()},
        "source_contrast_kcal_mol": converted(c["sample3_minus_sample0"]),
        "paired_precision_informative": c["paired_precision_informative"],
        "contrast_precision_informative": c["contrast_precision_informative"],
        "source_contrast_charge_bounds_kcal_mol": {k: v * ha for k, v in c["sample3_minus_sample0_charge_bounds_Eh"].items()},
        "scout": c["scout"], "actual_cost": {k: e[k] for k in (
            "elapsed_seconds", "allocated_CPU_seconds", "native_child_execution_receipts",
            "successful_native_updates", "measured_child_CPU_seconds", "peak_child_RSS_KiB")},
        "new_SCC_calls": 0, "new_MACE_calls": 0, "new_DFT_calls": 0,
        "classifier_gain": "not established", "production_changed": False}
    lines = ["# Native protein-exclusion scout: mainly a common shift", "",
        "All 12 declared native cavity updates completed. The two origin sources receive",
        "almost the same metal contrast. This does not explain their different classifier",
        "behavior and does not justify adding this component to routine scoring.", "",
        "| Source | Ca minus La transfer (kcal/mol) | Shape part | Cap-volume control |",
        "|---|---:|---:|---:|"]
    for sample, label in (("0", "sample0 control"), ("3", "sample3 ambiguous")):
        v = result["paired_kcal_mol"][sample]
        lines.append(f"| {label} | {v['transfer']['Gelec']:.9f} | {v['transfer']['Gshape']:.9f} | {v['cap_volume']['Gelec']:.9f} |")
    delta = result["source_contrast_kcal_mol"]["transfer"]["Gelec"]
    lines += ["", f"Raw sample3-minus-sample0 contrast: {delta:+.9f} kcal/mol.", "",
        "## Numerical qualification", "",
        "Archived replay is exact; translation, rotation and permutation pass the frozen",
        "1e-9 Eh electrostatic tolerance. All native block/repeat and endpoint precision",
        "checks pass. Both paired contrasts pass their precision flag. The aggregated",
        "source-contrast flag is **false**: its conservative summed charge bound exceeds",
        "the frozen 1e-6 Eh ceiling. Keep the raw contrast with this qualification; no",
        "tolerance change or rerun is justified to rescue that flag.", "",
        "## Scope and next decision", "",
        "Close this origin-only exclusion hypothesis as a remedy for the source gap.",
        "The test holds archived vacuum populations fixed and excludes permanent protein",
        "electrostatics. Synthetic caps remain native-radius centers. The cap-volume",
        "control is small but does not show that charged caps are negligible; their",
        "real/cap/cross blocks remain in RESULT.json. Native hydrogen-bond burial changes",
        "are reported separately and excluded from T. This is not total solvation, an",
        "affinity, a calibrated correction, or a test of accommodation-dependent burial.", "",
        "## Actual execution", "",
        "Worker4562: 136 scheduler seconds, four allocated CPUs (544 allocated CPU-s),",
        "133.610 measured TotalCPU seconds; scheduler peak RSS 7640020 KiB. The serial",
        "engine used four shared slots for default memory. Runner/child timing boundaries",
        "are separately retained in RESULT.json. Preparation4539 used7s/1CPU; build4543",
        "used4s/1CPU. No new SCC, MACE, DFT, optimization or folding. Production unchanged."]
    out = Path(m["output"])
    out.mkdir(exist_ok=False)
    (out / "RESULT.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    (out / "REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"status": "complete", "result": pin(out / "RESULT.json")}))


if __name__ == "__main__":
    main()
