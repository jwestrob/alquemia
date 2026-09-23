#!/usr/bin/env python3
"""Append experimental three-source results to the existing PLM protein export.

Copies all original strings. This is a join, not a rescore, recalibration or
transcript normalization. Biological preferences remain unknown.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from nikasha_plm_export import digest, parse_tsv, read_pin, require, unique


SCHEMA = "nikasha_PLM_candidate_overlay_v1"
FIELDS = [
    "candidate_status", "candidate_protocol_id", "candidate_result_sha256",
    "candidate_reference_sha256", "candidate_evidence_status",
    "candidate_declared_sources", "candidate_available_sources",
    "candidate_source_ids_json", "candidate_origin_R_model_kcal_mol_json",
    "candidate_origin_median_R_model_kcal_mol", "candidate_origin_range_kcal_mol",
    "candidate_R_model_kcal_mol_json", "candidate_median_R_model_kcal_mol",
    "candidate_range_kcal_mol", "candidate_decision",
    "candidate_mathematical_median_R_model_kcal_mol",
    "candidate_mathematical_decision", "candidate_source_details_json",
    "candidate_limitation",
]


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": digest(path.read_bytes())}


def compact(value):
    return json.dumps(value, separators=(",", ":"), allow_nan=False)


def scalar(value):
    return "" if value is None else str(value)


def build_overlay(proteins, result):
    table_pin, result_pin = pin(proteins), pin(result)
    table = parse_tsv(read_pin(table_pin))
    data = json.loads(read_pin(result_pin))
    require(not set(FIELDS) & set(table["fields"]), "Candidate fields already exist")
    original = unique(table["rows"], "target_id", "PLM protein")
    require(data["candidate_status"] == "experimental_not_promoted", "Expected a separately named experimental candidate")
    require(data["biological_accuracy_evaluated"] is False, "Unknown PLM labels are required")
    require(data["production_changed"] is False and data["original_DFT_and_production_results_changed"] is False,
            "Expected original results to remain unchanged")
    plan = json.loads(read_pin(data["plan"]))
    read_pin(data["reference"])
    require(plan["reference"] == data["reference"], "Plan/reference mismatch")
    request = json.loads(read_pin(plan["request"]))
    declared = unique(request["groups"], "protein_id", "declared group")
    groups = unique(data["groups"], "protein_id", "candidate group")
    sources = unique(data["rows"], "case_id", "candidate source")
    require(set(groups) == set(declared) <= set(original), "Candidate/request/PLM identity mismatch")
    require(len(groups) == data["group_denominator"], "Group denominator mismatch")
    require(len(sources) == data["source_denominator"], "Source denominator mismatch")
    seen, identities, additions = set(), {}, {}
    for protein, group in groups.items():
        members = group["members"]
        require(members == declared[protein]["members"] and len(members) == len(set(members)) == 3,
                f"Not the declared exact three-source group: {protein}")
        require(not seen.intersection(members), "A source occurs in multiple groups")
        seen.update(members)
        require(group["biological_accuracy_evaluated"] is False and group["classification_is_prediction"] is True,
                f"Unresolved prediction evidence changed: {protein}")
        details, sequence_pins = [], []
        for case in members:
            require(case in sources, f"Missing result row: {case}")
            row = sources[case]
            require(row["protein_id"] == protein, f"Source protein mismatch: {case}")
            evidence = group["source_evidence"][case]
            require(evidence["expected_class"] is None, f"Unexpected biological label: {case}")
            sequence_pin = evidence["source_provenance"]["native_AF3_input"]
            af3 = json.loads(read_pin(sequence_pin))
            sequences = [s["protein"]["sequence"] for s in af3["sequences"] if "protein" in s]
            require(len(sequences) == 1, f"Ambiguous protein sequence: {case}")
            sequence = sequences[0]
            require(digest(sequence.encode()) == original[protein]["sequence_sha256"] and
                    len(sequence) == int(original[protein]["full_length_aa"]),
                    f"Exact full-sequence identity mismatch: {case}")
            sequence_pins.append(sequence_pin)
            selected = {}
            if row["status"] == "available":
                for metal, choice in row["pool"]["rows"].items():
                    candidate = choice["operational_candidate"]
                    selected[metal] = {
                        "operational_candidate": candidate,
                        "mathematical_candidate": choice["mathematical_candidate"],
                        "work_from_origin": choice["work_from_origin_kcal_mol"][candidate],
                    }
            details.append({"case_id": case, "status": row["status"], "reason": row["reason"],
                            "variants": row["variants"], "selected": selected})
        available = sum(sources[c]["status"] == "available" for c in members)
        require(available == group["available"] and group["denominator"] == 3,
                f"Group coverage mismatch: {protein}")
        operational, mathematical = (group["variants"][v] for v in ("operational", "mathematical"))
        if available != 3:
            require(operational.get("R_model_kcal_mol") is None and mathematical.get("R_model_kcal_mol") is None,
                    f"Incomplete triple has a numerical aggregate: {protein}")
        origin = group["union_origin"]
        additions[protein] = dict(zip(FIELDS, [
            "experimental_available" if available == 3 else "experimental_unavailable",
            data["protocol_id"], result_pin["sha256"], data["reference"]["sha256"],
            "unknown_PLM_metal_labels_protocol_prediction_only", "3", str(available), compact(members),
            compact(origin.get("member_R_model_kcal_mol")), scalar(origin.get("median_R_model_kcal_mol")),
            scalar(origin.get("range_kcal_mol")), compact(operational.get("member_R_model_kcal_mol")),
            scalar(operational.get("R_model_kcal_mol")), scalar(operational.get("range_kcal_mol")),
            scalar(operational.get("decision")), scalar(mathematical.get("R_model_kcal_mol")),
            scalar(mathematical.get("decision")), compact(details), data["calibration_limitation"],
        ]))
        identities[protein] = {"sequence_sha256": original[protein]["sequence_sha256"],
                               "AF3_inputs": sequence_pins, "members": members}
    require(seen == set(sources), "Unassigned candidate source")
    rows = []
    for row in table["rows"]:
        extra = {field: "" for field in FIELDS}
        extra["candidate_status"] = "candidate_not_scored"
        extra.update(additions.get(row["target_id"], {}))
        rows.append(dict(row, **extra))
    receipt = {
        "schema": SCHEMA, "proteins_input": table_pin, "candidate_input": result_pin,
        "plan": data["plan"], "reference": data["reference"], "request": plan["request"],
        "implementation": pin(__file__), "protein_denominator": len(rows),
        "candidate_group_denominator": len(groups), "candidate_source_denominator": len(sources),
        "unscored_proteins": len(rows) - len(groups), "identities": identities,
        "original_columns_copied_verbatim": table["fields"], "added_columns": FIELDS,
        "new_molecular_calls": 0, "new_expression_normalization": False,
        "biological_accuracy_evaluated": False,
        "scale_warning": "Original DFT S and candidate raw R have different gauges; do not subtract them.",
    }
    return table["fields"] + FIELDS, rows, receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", required=True, type=Path)
    parser.add_argument("--result", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path, help="New output directory; existing directories are rejected")
    args = parser.parse_args()
    fields, rows, receipt = build_overlay(args.proteins, args.result)
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output / "proteins.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    receipt["proteins_output"] = pin(output)
    (args.output / "RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"output": str(output.resolve()), "proteins": len(rows),
                      "candidate_groups": receipt["candidate_group_denominator"]}))


if __name__ == "__main__":
    main()
