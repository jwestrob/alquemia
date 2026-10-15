#!/usr/bin/env python3
"""Compare saved common pools after the declared Gelec-only transfer diagnostic."""
import argparse
import hashlib
import json
import os
from pathlib import Path


def pin(p):
    p = Path(p).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def read(p):
    assert pin(p["path"]) == p
    return json.loads(Path(p["path"]).read_text())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--manifest-sha", required=True)
    a = p.parse_args()
    assert os.environ.get("SLURM_JOB_ID"), "Saved numerical comparison requires Slurm"
    assert pin(a.manifest)["sha256"] == a.manifest_sha
    m = json.loads(a.manifest.read_text())
    assert pin(__file__) == m["implementation"]
    new, old, baseline = [read(m[k]) for k in ("new_collection", "origin_collection", "baseline")]
    assert new["status"] == old["status"] == "complete"
    assert new["origin_collection"] == m["origin_collection"]
    new_by = {r["id"]: r for r in new["rows"]}
    origin_by = {r["id"]: r for r in old["rows"]}
    order = ("origin", "adaptive_Ca", "adaptive_La")
    ha = new["hartree_to_kcal_mol"]
    assert ha == old["hartree_to_kcal_mol"]
    rows = []
    for base in baseline["rows"]:
        sample = base["sample"]
        assert sample in (0, 3)
        metals = {}
        for metal in ("Ca", "La"):
            origin = origin_by[f"sample{sample}_{metal}"]
            records = {"origin": origin, **{
                q: new_by[f"sample{sample}_{q}_{metal}"] for q in order[1:]}}
            assert all(r["status"] == "available" for r in records.values())
            # Baseline values are exact saved works, each relative to its own origin.
            works = {q: base["all_origin_relative_works_kcal_mol"][metal][q]["total"] for q in order}
            historical = base["pools"]["historical3"]["mathematical"][metal]
            assert all(abs(works[q] - historical["all_works_kcal_mol"][q]) <= 1e-8 for q in order)
            before = min(order, key=works.get)
            assert before == historical["candidate"]
            assert abs(works[before] - historical["work_kcal_mol"]) <= 1e-8
            transfer = {q: records[q]["transfer_Eh"]["Gelec"] * ha for q in order}
            diagnostic = {q: works[q] + transfer[q] for q in order}
            relative = {q: diagnostic[q] - diagnostic["origin"] for q in order}
            after = min(order, key=diagnostic.get)
            change = diagnostic[after] - works[before]
            work_change = relative[after] - works[before]
            assert abs(change - transfer["origin"] - work_change) <= 1e-8
            metals[metal] = {"old_works_kcal_mol": works, "transfers_kcal_mol": transfer,
                "diagnostic_works_kcal_mol": relative,
                "old_mathematical_selected": before, "diagnostic_selected": after,
                "old_operational_selected": base["pools"]["historical3"]["operational"][metal]["candidate"],
                "pool_energy_change_kcal_mol": change,
                "accommodation_work_change_kcal_mol": work_change,
                "input_local_pins": {q: records[q]["archived_local"] for q in order}}
        rchange = metals["Ca"]["pool_energy_change_kcal_mol"] - metals["La"]["pool_energy_change_kcal_mol"]
        accommodation = metals["Ca"]["accommodation_work_change_kcal_mol"] - metals["La"]["accommodation_work_change_kcal_mol"]
        rows.append({"sample": sample, "source_case_id": base["source_case_id"], "metals": metals,
            "paired_pool_change_kcal_mol": rchange,
            "paired_accommodation_change_kcal_mol": accommodation,
            "origin_paired_change_kcal_mol": rchange - accommodation,
            "classification": None})
    by = {r["sample"]: r for r in rows}
    assert set(by) == {0, 3}
    result = {"protocol": "native_capped_union_accommodation_common3_diagnostic_v1",
        "status": "complete", "inputs": m, "job_id": os.environ["SLURM_JOB_ID"],
        "candidate_order": order, "rows": rows,
        "source_contrast_change_kcal_mol": by[3]["paired_pool_change_kcal_mol"] - by[0]["paired_pool_change_kcal_mol"],
        "source_accommodation_change_kcal_mol": by[3]["paired_accommodation_change_kcal_mol"] - by[0]["paired_accommodation_change_kcal_mol"],
        "response_precision_flags": {q: {
            sample: v["precision_informative"] if v else None
            for sample, v in source.items()}
            for q, source in new["Ca_minus_La_accommodation_changes"].items()},
        "source_response_precision_flags": {q: v["precision_informative"] if v else None for q, v in new["sample3_minus_sample0"].items()},
        "uncertainty_scope": "Original conservative charge-rounding bounds/flags stay in source collections; no full model uncertainty or confidence is inferred.",
        "classification": None, "production_changed": False,
        "interpretation": "Mathematical row minima of Gelec-only augmentation, not complete protein solvation or calibrated classification. No thresholds fitted."}
    lines = ["# Cavity-dependent accommodation: saved common-pool comparison", "",
        "All values below describe the declared electrostatic-only diagnostic. They",
        "are not new class assignments or a complete protein-solvent correction.", "",
        "| Source | Paired pool change | Accommodation part | Ca selection | La selection |",
        "|---|---:|---:|---|---|"]
    for row in rows:
        selected = [row["metals"][z] for z in ("Ca", "La")]
        lines.append(f"| sample{row['sample']} | {row['paired_pool_change_kcal_mol']:+.6f} | {row['paired_accommodation_change_kcal_mol']:+.6f} | "
            + " | ".join(f"{v['old_mathematical_selected']} → {v['diagnostic_selected']}" for v in selected) + " |")
    lines += ["", "Energy changes are kcal/mol; positive Ca-minus-La changes are La-directed.",
        f"Source3-minus-source0 paired change: {result['source_contrast_change_kcal_mol']:+.6f} kcal/mol.",
        f"Accommodation contribution to that change: {result['source_accommodation_change_kcal_mol']:+.6f} kcal/mol.", "",
        "Conservative aggregate precision flags are retained in RESULT.json and source",
        "collections. Passing molecular execution/per-endpoint gates does not make those",
        "aggregate flags pass. Native HB/nonpolar burial and permanent environment charges",
        "remain excluded. Synthetic-cap controls do not establish a source-only cavity.", "",
        "No new electronic endpoints, geometries, thresholds or production changes."]
    out = Path(m["output"])
    out.mkdir(exist_ok=False)
    (out / "RESULT.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    (out / "REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"status": "complete", "result": pin(out / "RESULT.json")}))


if __name__ == "__main__":
    main()
