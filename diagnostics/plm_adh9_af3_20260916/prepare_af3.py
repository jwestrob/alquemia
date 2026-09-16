#!/usr/bin/env python3
"""Adapt explicitly selected AF3 files to the unchanged PLM fixed-core wrapper.

Only filename adaptation is performed: raw CIF and confidence JSON bytes are
preserved. This module never selects predictions, changes confidence values,
runs ORCA, or submits a job. Failed selected models remain unsupported.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

import gemmi

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / "diagnostics/plm_adh9_fixed_core_20260916"
PROTOCOL = "pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3"
ALLOWED_TARGETS = {"PQQSEQ_48f861015fad150af40a", "PQQSEQ_13d74836d4b7a3e02140"}


def record(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(item):
    observed = record(item["path"])
    if observed["sha256"] != item["sha256"]:
        raise ValueError(f"File hash mismatch: {item['path']}")
    return Path(observed["path"])


def read(path):
    return json.loads(Path(path).read_text())


def write(path, payload):
    with Path(path).open("x") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def load_wrapper(path):
    spec = importlib.util.spec_from_file_location("unchanged_plm_fixed_core_wrapper", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def verify_pins(path):
    pins = read(path)
    if pins["schema_version"] != "plm.adh9.af3_bridge_pins.v1" or pins["protocol_id"] != PROTOCOL:
        raise ValueError("Wrong AF3 bridge pin schema/protocol")
    if verify(pins["bridge"]) != Path(__file__).resolve():
        raise ValueError("Executing bridge differs from pinned bridge")
    for name in ("original_wrapper", "original_wrapper_pins", "prior_candidate_manifest",
                 "calibration_implementation_pins", "calibration_result", "preparation_python"):
        verify(pins[name])
    wrapper = load_wrapper(verify(pins["original_wrapper"]))
    # This checks every original helper, package version, and frozen calibration
    # record, without modifying or preparing any historical structure.
    previous_pins, original, _ = wrapper.verify_pins(verify(pins["original_wrapper_pins"]))
    if previous_pins["wrapper"] != pins["original_wrapper"]:
        raise ValueError("Previous wrapper authority differs")
    if previous_pins["candidate_manifest"] != pins["prior_candidate_manifest"]:
        raise ValueError("Previous reviewed roles differ")
    if set(pins["allowed_target_ids"]) != ALLOWED_TARGETS:
        raise ValueError("Approved target set changed")
    return pins, previous_pins, original


def clean_roles(roles):
    return {role: {key: value[key] for key in ("chain", "resname", "resnum", "icode")}
            for role, value in roles.items()}


def protein_residues(path, chain):
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ValueError("Expected one raw structural model")
    matches = [c for c in structure[0] if c.name == chain]
    if len(matches) != 1:
        raise ValueError("Protein chain must resolve exactly once")
    return [(r.seqid.num, r.seqid.icode.strip(), r.name) for r in matches[0]]


def validate_confidence(summary, chain_order):
    if chain_order != ["A", "B", "C"]:
        raise ValueError("This campaign requires the explicit native A/B/C chain order")
    values = summary.get("chain_pair_iptm")
    if (not isinstance(values, list) or len(values) != 3
            or any(not isinstance(row, list) or len(row) != 3 for row in values)):
        raise ValueError("AF3 native chain_pair_iptm must be a 3-by-3 matrix")
    # AF3 legitimately emits null for undefined entries such as an ion's
    # diagonal. Preserve those nulls; only the tested protein–La entry must
    # contain a finite score. Never impute or manufacture a confidence value.
    if any(x is not None and (not isinstance(x, (int, float)) or isinstance(x, bool)
                             or not math.isfinite(x) or not 0 <= x <= 1)
           for row in values for x in row):
        raise ValueError("AF3 native chain_pair_iptm contains an invalid value")
    if values[0][1] is None:
        raise ValueError("AF3 protein–La confidence is undefined")
    return float(values[0][1])


def adapt_case(case, prior, output, selection_record):
    source = verify(case["source_cif"])
    confidence = verify(case["summary_confidences"])
    roles = clean_roles(case["roles"])
    if roles != prior["roles"]:
        raise ValueError("AF3 role selectors differ from the frozen same-sequence comparison")
    if case["protein_chain"] != "A" or case["metal"]["chain"] != "B" or case["pqq"]["chain"] != "C":
        raise ValueError("Expected protein A / La B / PQQ C from the approved AF3 input")
    raw_model = gemmi.read_structure(str(source))
    order = case.get("confidence_chain_order", ["A", "B", "C"])
    if len(raw_model) != 1 or [c.name for c in raw_model[0]] != order:
        raise ValueError("Raw CIF order does not match the declared native confidence chain order")
    if protein_residues(source, "A") != protein_residues(verify(prior["source_cif"]), "A"):
        raise ValueError("AF3 protein sequence/numbering differs from the original scored monomer")
    native_summary = read(confidence)
    iptm = validate_confidence(native_summary, order)
    destination = output / "adapter_inputs"
    destination.mkdir(exist_ok=True)
    cif_copy = destination / (case["case_id"] + ".cif")
    summary_copy = destination / (case["case_id"] + "_summary.json")
    for src, dst in ((source, cif_copy), (confidence, summary_copy)):
        if dst.exists():
            raise ValueError("Refusing to overwrite adapter input")
        shutil.copyfile(src, dst)
        if record(src)["sha256"] != record(dst)["sha256"]:
            raise ValueError("Byte-preserving copy failed")
    provenance = {"selected_models": selection_record, "raw_AF3_model": record(source),
                  "raw_AF3_summary_confidences": record(confidence),
                  "original_Protenix_model": prior["source_cif"]}
    for name, item in case.get("source_provenance", {}).items():
        provenance["selector_" + name] = record(verify(item))
    target = {"case_id": case["case_id"], "target_id": case["target_id"], "rank": case["rank"],
        "source_cif": record(cif_copy), "summary_confidence": record(summary_copy),
        "source_provenance": provenance, "protein_chain": "A", "metal": case["metal"],
        "pqq": case["pqq"], "roles": roles}
    receipt = {"case_id": case["case_id"], "raw_AF3_source_cif": record(source),
        "raw_AF3_summary_confidences": record(confidence), "adapter_source_cif": record(cif_copy),
        "adapter_summary_confidence": record(summary_copy), "model_bytes_identical": True,
        "summary_bytes_identical": True, "confidence_values_changed": False,
        "confidence_chain_order": order, "native_protein_La_chain_pair_iptm": iptm,
        "same_protein_sequence_and_numbering": True, "roles_identical_to_previous_review": True}
    return target, receipt


def prepare(selected_models, agreement, implementation_pins, output):
    pins, previous_pins, original = verify_pins(implementation_pins)
    selection_record, agreement_record = record(selected_models), record(agreement)
    selected = read(selected_models)
    cases = selected["cases"]
    if len(cases) != 2 or {c["target_id"] for c in cases} != ALLOWED_TARGETS:
        raise ValueError("Exactly the two approved targets, each once, are required")
    if len({c["case_id"] for c in cases}) != 2 or any(not re.fullmatch(r"[A-Za-z0-9_-]+", c["case_id"]) for c in cases):
        raise ValueError("Duplicate or unsafe case IDs")
    prior = {c["target_id"]: c for c in read(verify(pins["prior_candidate_manifest"]))["targets"]}
    output = Path(output).resolve()
    if not output.is_relative_to((ROOT / "workspaces/plm_adh9_af3_20260916").resolve()):
        raise ValueError("AF3 preparation must remain in its separate candidate workspace")
    output.mkdir(parents=True, exist_ok=False)
    adapted, bridges, failed = [], [], {}
    for case in cases:
        if case["status"] == "unsupported":
            if not case.get("reason"):
                raise ValueError("Unsupported selected target lacks a reason")
            failed[case["case_id"]] = case["reason"]
            continue
        if case["status"] != "selected":
            raise ValueError("Selection status must be selected or unsupported")
        try:
            target, bridge = adapt_case(case, prior[case["target_id"]], output, selection_record)
            adapted.append(target)
            bridges.append(bridge)
        except Exception as exc:
            failed[case["case_id"]] = f"AF3 adapter rejected input: {type(exc).__name__}: {exc}"
    candidate_path = output / "candidate_manifest.json"
    write(candidate_path, {"schema_version": "plm.adh9.candidate_manifest.v1", "protocol_id": PROTOCOL,
        "review_authority": selection_record, "approval": agreement_record, "targets": adapted})
    candidate_pins = dict(previous_pins)
    candidate_pins["candidate_manifest"] = record(candidate_path)
    candidate_pins["authorities"] = {"AF3_selection": selection_record, "approval": agreement_record,
        "bridge_pins": record(implementation_pins), "raw_provenance": record(pins["prior_candidate_manifest"]["path"])}
    candidate_pin_path = output / "candidate_implementation_pins.json"
    write(candidate_pin_path, candidate_pins)
    write(output / "adapter_receipt.json", {"schema_version": "plm.adh9.af3_adapter.v1",
        "protocol_id": PROTOCOL, "implementation_pins": record(implementation_pins),
        "selected_models": selection_record, "agreement": agreement_record,
        "bridges": bridges, "adapter_failures": failed,
        "unchanged_wrapper": pins["original_wrapper"], "model_selection_performed": False,
        "confidence_values_changed": False, "orca_executed": False})
    command = [str(verify(pins["preparation_python"])), str(verify(pins["original_wrapper"])),
        "--manifest", str(candidate_path), "--implementation-pins", str(candidate_pin_path),
        "--output", str(output / "pairs")]
    # The unchanged preparation wrapper has no quantum or scheduler entry point.
    process = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output / "preparation.log").write_text(process.stdout)
    if process.returncode:
        raise RuntimeError(f"Unchanged preparation wrapper failed (exit {process.returncode}); retained log")
    inner = read(output / "pairs/prepared_pairs.json")
    ready_by_id = {c["case_id"]: c for c in inner["cases"]}
    final = []
    for case in cases:
        cid = case["case_id"]
        if cid in failed:
            result = {k: case[k] for k in ("case_id", "target_id", "rank")}
            result.update(status="unsupported", reason=failed[cid], carve_manifest=None,
                          source_cif=case.get("source_cif"))
        else:
            result = dict(ready_by_id[cid])
        result["predictor"] = "AlphaFold3"
        result["raw_AF3_source_cif"] = case.get("source_cif")
        result["raw_AF3_summary_confidences"] = case.get("summary_confidences")
        result["selected_sample"] = case.get("sample")
        final.append(result)
    report = {"schema_version": "plm.adh9.prepared_pairs.v1", "protocol_id": PROTOCOL,
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "cases": final,
        "selected_models": selection_record, "agreement": agreement_record,
        "adapter_receipt": record(output / "adapter_receipt.json"),
        "parent_preparation": record(output / "pairs/prepared_pairs.json"),
        "implementation_pins": record(implementation_pins),
        "prepared_pair_count": sum(c["status"] == "ready_for_orca" for c in final),
        "unsupported_count": sum(c["status"] == "unsupported" for c in final), "orca_executed": False}
    write(output / "prepared_pairs.json", report)
    return output / "prepared_pairs.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-models", required=True, type=Path)
    parser.add_argument("--agreement", required=True, type=Path)
    parser.add_argument("--implementation-pins", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(prepare(args.selected_models, args.agreement, args.implementation_pins, args.output))


if __name__ == "__main__":
    main()
