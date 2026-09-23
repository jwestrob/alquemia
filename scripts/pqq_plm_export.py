#!/usr/bin/env python3
"""Join disjoint completed PLM batches, preserving original biological evidence."""
import argparse
import csv
import json
from pathlib import Path

from nikasha_plm_export import parse_tsv, require
from plm_candidate_overlay import FIELDS, build_overlay, pin


def export(proteins, results, preparations, output):
    require(results or preparations, "At least one actual result or preparation ledger required")
    original = parse_tsv(Path(proteins).read_bytes())
    require(not set(FIELDS) & set(original["fields"]), "Use the original protein table without an old overlay")
    rows = [{**r, **{k: "" for k in FIELDS}, "candidate_status": "candidate_not_scored"}
            for r in original["rows"]]
    by_id = {r["target_id"]: r for r in rows}
    require(len(by_id) == len(rows), "Duplicate original protein IDs")
    receipts, seen, references = [], set(), set()
    for result in results:
        _, overlay, receipt = build_overlay(proteins, result)
        scored = {r["target_id"] for r in overlay if r["candidate_status"] != "candidate_not_scored"}
        require(not seen & scored, "Overlapping result proteins; select one declared run per protein explicitly")
        seen.update(scored)
        references.add(receipt["reference"]["sha256"])
        require(len(references) == 1, "Mixed candidate references cannot form one cohort export")
        for r in overlay:
            if r["target_id"] in scored:
                by_id[r["target_id"]].update({k: r[k] for k in FIELDS})
        receipts.append(receipt)
    unavailable = set()
    for preparation in preparations:
        data = json.loads(Path(preparation).read_text())
        require(data["schema"] == "Nikasha_PLM_three_source_operations_v1", "Unexpected preparation ledger")
        for g in data["groups"]:
            pid = g["protein_id"]
            require(pid in by_id, "Preparation protein absent from the biological join: " + pid)
            if g["status"] != "prepared":
                require(pid not in seen and pid not in unavailable, "Conflicting failed/scored or duplicate preparation")
                unavailable.add(pid)
                by_id[pid].update(candidate_status="candidate_unavailable_preparation",
                    candidate_declared_sources="3", candidate_limitation=json.dumps(g, separators=(",", ":")))
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    table = out / "proteins.tsv"
    with table.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=original["fields"] + FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    receipt = {"schema": "Nikasha_PLM_disjoint_batch_export_v1", "implementation": pin(__file__),
               "original_proteins": pin(proteins), "results": receipts,
               "preparations": [pin(p) for p in preparations], "output": pin(table),
               "proteins": len(rows), "result_groups": len(seen), "preparation_unavailable": len(unavailable),
               "not_scored": len(rows) - len(seen) - len(unavailable),
               "original_columns_unchanged": original["fields"], "new_expression_normalization": False,
               "new_molecular_calls": 0}
    (out / "RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return {k: receipt[k] for k in ("output", "proteins", "result_groups", "preparation_unavailable", "not_scored")}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--proteins", required=True)
    p.add_argument("--results", nargs="*", default=[])
    p.add_argument("--preparations", nargs="*", default=[])
    p.add_argument("--output", required=True)
    print(json.dumps(export(**vars(p.parse_args())), indent=2))


if __name__ == "__main__":
    main()
