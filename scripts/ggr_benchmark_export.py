"""Append approved GGR A/B development records to an immutable benchmark release.

All inherited fields keep their original values. New scores live separately in
``development_scores``; ``development_export`` describes this derived snapshot.
The exporter launches no calculations and assigns no experimental labels.
"""
from pathlib import Path
from collections import Counter
from copy import deepcopy
import argparse
import json
import shutil

from affordable_common import (InvalidArtifact, cache_key, contrast, energy,
                               read_json, record, verify, write_new)
from ggr_workflow import ALLOWED_PROTOCOLS, checked_rows


ROOT = Path(__file__).resolve().parents[1]
GROUP_ALIASES = {
    "GGR_1GLG": ("ggr", "GGR_MglB"),
    "AEQUORIN_1SL8_VECTOR": ("aequorin", "aequorin_P07164"),
    "ALACTA_BOVINE_STRONG_SITE": ("alacta_bovine", "bovine_alpha_lactalbumin"),
}


def validate_parent(parent):
    rows = parent["rows"]
    if parent["record_count"] != 68 or len(rows) != 68:
        raise InvalidArtifact("the approved parent must retain all 68 records")
    targets = {row["target_id"]: row for row in rows}
    if len(targets) != len(rows):
        raise InvalidArtifact("duplicate parent target_id")
    if "development_export" in parent or any(
            "development_scores" in row or "development_response_records" in row for row in rows):
        raise InvalidArtifact("use the original release, not an earlier development export")
    for target, (_, canonical) in GROUP_ALIASES.items():
        if target not in targets or targets[target]["biological_group"] != canonical:
            raise InvalidArtifact("unexpected existing biological group")
    return targets


def checked_collection(path):
    """Check recorded results against their pinned tasks; never recollect jobs."""
    result = checked_rows(path)
    manifest = read_json(verify(result["manifest"]))
    if manifest.get("schema_version") != "alquemia.ggr_tasks.v1" or manifest["stage"] not in ("A", "B"):
        raise InvalidArtifact("only approved stage A/B task collections are supported")
    for key in ("approved_plan", "agreement", "configuration", "release"):
        verify(manifest[key])
    verify(manifest["reference"]["manifest"])
    cases = {row["case"]: row for row in manifest["cases"]}
    tasks = {task["task_id"]: task for task in manifest["tasks"]}
    if len(cases) != len(manifest["cases"]) or len(tasks) != len(manifest["tasks"]):
        raise InvalidArtifact("duplicate task or case")
    if [row["case"] for row in result["rows"]] != [row["case"] for row in manifest["cases"]]:
        raise InvalidArtifact("collection must retain every manifest case in order")
    completed = 0
    for row in result["rows"]:
        case = cases[row["case"]]
        if any(row.get(key) != value for key, value in case.items()):
            raise InvalidArtifact("collection case differs from its task manifest")
        if row["protocol_id"] not in ALLOWED_PROTOCOLS or row["evaluation_role"] != "method_development_already_consumed":
            raise InvalidArtifact("unapproved protocol or evaluation role")
        if set(row["endpoints"]) != {"La", "Ca"} or len(row["tasks"]) != 2:
            raise InvalidArtifact("exactly two paired endpoints required")
        for metal, task_id in zip(("La", "Ca"), row["tasks"]):
            task = tasks[task_id]
            if task["metal"] != metal or task["case"] != row["case"]:
                raise InvalidArtifact("task/endpoint identity mismatch")
            for key in ("input", "xyz", "source_input", "source_xyz", "source_manifest"):
                verify(task[key])
            endpoint = row["endpoints"][metal]
            if endpoint["status"] == "complete":
                output = verify(endpoint["output"])
                receipt = read_json(verify(endpoint["receipt"]))
                if output != Path(task["output_path"]) or energy(output) != endpoint["energy_hartree"]:
                    raise InvalidArtifact("endpoint output/energy mismatch")
                if receipt["manifest"] != result["manifest"] or receipt["task_id"] != task_id:
                    raise InvalidArtifact("receipt task/manifest mismatch")
                if not receipt["normal_termination"] or not receipt["scf_converged"] or receipt["returncode"] != 0:
                    raise InvalidArtifact("unsuccessful receipt presented as complete")
                completed += 1
            elif endpoint.get("energy_hartree") is not None:
                raise InvalidArtifact("unavailable endpoint must not supply an energy")
        pair_complete = all(ep["status"] == "complete" for ep in row["endpoints"].values())
        if pair_complete:
            expected = contrast(row["endpoints"]["Ca"]["energy_hartree"],
                                row["endpoints"]["La"]["energy_hartree"],
                                manifest["reference"]["delta_E_aquo_hartree"])
            if row["status"] != "complete" or row["score"] != expected or row["decision"] != "uncalibrated_protocol":
                raise InvalidArtifact("score algebra, units, or uncalibrated decision mismatch")
        elif row["score"] is not None or row["status"] == "complete" or row["decision"] != "unavailable":
            raise InvalidArtifact("incomplete pair must retain unavailable score/decision")
    if result["completed_endpoints"] != completed or result["total_endpoints"] != len(tasks):
        raise InvalidArtifact("collection endpoint denominator mismatch")
    return result, manifest, tasks


def development_record(row, collection_path, collection, manifest, tasks, target):
    alias, canonical = GROUP_ALIASES[row["target_id"]]
    if row["biological_group"] != alias or target["biological_group"] != canonical:
        raise InvalidArtifact("execution biological group does not match the existing target")
    if row["observation_group"] != target["target_id"]:
        raise InvalidArtifact("structural replicate must retain the existing observation group")
    preparation_path = verify(row["source_manifest"])
    preparation = read_json(preparation_path)
    if preparation["protocol_id"] != row["protocol_id"]:
        raise InvalidArtifact("preparation protocol mismatch")
    for key in ("source_structure", "topology_definition"):
        verify(preparation[key])
    if preparation.get("protonation_manifest") is not None:
        verify(preparation["protonation_manifest"])
    graph = preparation["atom_graph"]
    if not graph["source_to_qm"] or not graph["retained_bonds"] or not graph["cut_bonds_and_caps"]:
        raise InvalidArtifact("preparation lacks source/cap mapping")
    result = deepcopy(row)
    result.update(
        development_id=cache_key({"manifest": collection["manifest"], "case": row["case"]}),
        biological_group=canonical,
        execution_biological_group=alias,
        independent_biological_observation=False,
        site_specific_affinity_label_added=False,
        calibrated_decision=None,
        experimental_label_policy="existing target evidence retained; no new site-specific label",
        preparation_record=preparation,
        protonation_manifest_status=("recorded_and_verified" if preparation.get("protonation_manifest") is not None
                                     else "not_recorded_in_inherited_preparation"),
        endpoint_tasks={metal: deepcopy(tasks[task_id]) for metal, task_id in zip(("La", "Ca"), row["tasks"])},
        provenance={
            "collection": record(collection_path), "task_manifest": collection["manifest"],
            "preparation_manifest": record(preparation_path),
            "reference": deepcopy(manifest["reference"]),
            "agreement": manifest["agreement"], "approved_plan": manifest["approved_plan"],
            "configuration": manifest["configuration"],
            "implementation_snapshot": manifest["implementation_snapshot"],
            "ggr_implementation_snapshot": manifest["ggr_implementation_snapshot"],
            "execution_policy": deepcopy(manifest["execution_policy"]), "orca": manifest["orca"],
        },
    )
    return result


def response_record(path):
    """Preserve Stage C checks and all actual endpoints without creating scores."""
    path = Path(path).resolve()
    result = read_json(path)
    if result.get("schema_version") != "alquemia.ggr_sensitivity_comparison.v1":
        raise InvalidArtifact("unsupported sensitivity collection schema")
    manifest = read_json(verify(result["manifest"]))
    if (manifest.get("schema_version") != "alquemia.ggr_sensitivity.v1" or
            manifest.get("protocol_id") != "ggr_local_sensitivity_tightscf_dev_v1" or
            result["phase"] != manifest["phase"] or result["phase"] not in ("nominal", "half")):
        raise InvalidArtifact("sensitivity protocol/phase mismatch")
    for key in ("agreement", "plan"):
        verify(manifest[key])
    for item in (result, *result["comparisons"]):
        if (item.get("response_status") != "response_model_not_validated" or
                item.get("relaxation_correction_kcal_mol") is not None or
                item.get("entropy_correction_kcal_mol") is not None):
            raise InvalidArtifact("sensitivity diagnostics cannot supply mechanical corrections")
    if result["absolute_classification"] != "uncalibrated_protocol":
        raise InvalidArtifact("sensitivity diagnostic cannot create a calibrated classification")
    normal = read_json(verify(manifest["normal_manifest"]))
    if set(manifest["representations"]) != {"extended", "connected"}:
        raise InvalidArtifact("sensitivity requires the two approved GGR representations")
    for representation in manifest["representations"].values():
        preparation_pin = representation["preparation"]
        verify(preparation_pin)
        matches = [case for case in normal["cases"] if case["source_manifest"] == preparation_pin]
        if (len(matches) != 1 or matches[0]["target_id"] != "GGR_1GLG" or
                matches[0]["source_structure_id"] != "1GLG"):
            raise InvalidArtifact("sensitivity preparation is not the approved GGR 1GLG case")
    tasks = {task["task_id"]: task for task in manifest["tasks"]}
    if set(tasks) != set(result["energies"]):
        raise InvalidArtifact("sensitivity collection must retain every task endpoint")
    for task in tasks.values():
        for key in ("input", "xyz"):
            verify(task[key])
    for key, endpoint in result["energies"].items():
        for artifact in ("receipt", "output"):
            if artifact in endpoint:
                verify(endpoint[artifact])
        if endpoint["status"] == "complete":
            receipt = read_json(verify(endpoint["receipt"]))
            if (receipt["manifest"] != result["manifest"] or receipt["task_id"] != key or
                    energy(verify(endpoint["output"])) != endpoint["energy_hartree"]):
                raise InvalidArtifact("sensitivity endpoint/receipt mismatch")
    for gradient in result["gradients"].values():
        for artifact in gradient["artifacts"].values():
            verify(artifact)
    if result["phase"] == "half":
        for key in ("nominal_manifest", "trigger_comparison"):
            verify(manifest[key])
    return {
        "collection": record(path), "manifest": result["manifest"],
        "target_id": "GGR_1GLG", "biological_group": "GGR_MglB",
        "evaluation_role": "method_development_already_consumed",
        "independent_biological_observation": False,
        "collection_record": result, "task_manifest_record": manifest,
        "phase": result["phase"],
        "completed_endpoints": sum(endpoint["status"] == "complete" for endpoint in result["energies"].values()),
        "total_endpoints": len(tasks),
        "physical_check_statuses": dict(Counter(item["status"] for item in result["comparisons"])),
        "response_status": "response_model_not_validated",
        "relaxation_correction_kcal_mol": None, "entropy_correction_kcal_mol": None,
        "calibrated_decision": None,
        "interpretation": "Physical gradient/displacement diagnostics of the isolated CPCM model; no added benchmark affinity labels or displacement scores.",
    }


def export_ledger(parent_release, collections, output, sensitivities=()):
    """Write a new directory; every original row field remains exactly intact."""
    parent_release, output = Path(parent_release).resolve(), Path(output).resolve()
    if output.exists() or not output.is_relative_to(ROOT / "workspaces"):
        raise InvalidArtifact("a new output directory under workspaces is required")
    if not collections:
        raise InvalidArtifact("at least one explicit collection path required")
    parent = read_json(parent_release)
    validate_parent(parent)
    result = deepcopy(parent)
    targets = {row["target_id"]: row for row in result["rows"]}
    seen_stages, sources, summaries = set(), [], []
    for path in collections:
        path = Path(path).resolve()
        collection, manifest, tasks = checked_collection(path)
        stage = manifest["stage"]
        if stage in seen_stages:
            raise InvalidArtifact("duplicate stage collection; select one explicit snapshot per stage")
        seen_stages.add(stage)
        for row in collection["rows"]:
            if row["target_id"] not in GROUP_ALIASES:
                raise InvalidArtifact("development case does not map to an approved existing target")
            target = targets[row["target_id"]]
            target.setdefault("development_scores", []).append(
                development_record(row, path, collection, manifest, tasks, target))
        sources.append({
            "stage": stage, "collection": record(path), "task_manifest": collection["manifest"],
            "collection_metadata": {key: deepcopy(value) for key, value in collection.items() if key != "rows"},
            "unsupported_preparations": deepcopy(manifest.get("unsupported_preparations", [])),
        })
    response_records = [response_record(path) for path in sensitivities]
    if len({item["phase"] for item in response_records}) != len(response_records):
        raise InvalidArtifact("select one explicit sensitivity snapshot per phase")
    if response_records:
        plans = {(read_json(verify(source["task_manifest"]))["approved_plan"]["sha256"],
                  read_json(verify(source["task_manifest"]))["agreement"]["sha256"]) for source in sources}
        for response in response_records:
            manifest = response["task_manifest_record"]
            if (manifest["plan"]["sha256"], manifest["agreement"]["sha256"]) not in plans:
                raise InvalidArtifact("sensitivity and A/B collections have different approved scopes")
        targets["GGR_1GLG"]["development_response_records"] = response_records
    for row, old in zip(result["rows"], parent["rows"]):
        if {key: row[key] for key in old} != old:
            raise InvalidArtifact("an inherited benchmark row was changed")
        for score in row.get("development_scores", []):
            summaries.append({key: score[key] for key in (
                "case", "target_id", "biological_group", "stage", "source_structure_id",
                "lane", "protocol_id", "status", "score", "decision")})
    # Validate everything before writing the new derived release.
    output.mkdir(parents=True, exist_ok=False)
    implementation = {}
    for name in ("ggr_benchmark_export.py", "ggr_workflow.py", "affordable_common.py"):
        dest = output / "implementation" / name
        dest.parent.mkdir(exist_ok=True)
        shutil.copyfile(ROOT / "scripts" / name, dest)
        implementation[name] = record(dest)
    result["development_export"] = {
        "schema_version": "alquemia.ggr_benchmark_development_export.v1",
        "parent_release": record(parent_release), "implementation": implementation,
        "status": ("complete" if seen_stages == {"A", "B"} and
                   all(source["collection_metadata"]["status"] == "complete" for source in sources) else "incomplete"),
        "stages_included": sorted(seen_stages), "stages_not_included": sorted({"A", "B"} - seen_stages),
        "collections": sources, "record_count_preserved": len(result["rows"]),
        "development_pair_count": len(summaries),
        "development_pairs_by_biological_group": dict(Counter(row["biological_group"] for row in summaries)),
        "completed_endpoints": sum(source["collection_metadata"]["completed_endpoints"] for source in sources),
        "total_endpoints": sum(source["collection_metadata"]["total_endpoints"] for source in sources),
        "sensitivity_collections": [item["collection"] for item in response_records],
        "sensitivity_phases_included": [item["phase"] for item in response_records],
        "sensitivity_completed_endpoints": sum(item["completed_endpoints"] for item in response_records),
        "sensitivity_total_endpoints": sum(item["total_endpoints"] for item in response_records),
        "summary": summaries, "default_changed": False, "experimental_labels_changed": False,
        "new_independent_biological_observations": 0, "new_calibrated_decisions": 0,
        "response_status": "response_model_not_validated",
        "relaxation_correction_kcal_mol": None,
        "interpretation": (
            "All inherited top-level fields describe the pinned parent release; this block describes the appended development snapshot. "
            "Original row fields, biological groups, scores, labels and thresholds are unchanged. "
            "Aequorin EF3 is a boundary control without a site-specific affinity label; its original ordered site vector remains intact. "
            "GGR structures and both alpha-lactalbumin structures remain grouped by their existing targets. "
            "S is the existing aquo reporting gauge, with no calibrated zero or inherited PQQ bands. "
            "Stage C gradients and displacements are retained separately on GGR as sensitivity diagnostics, not additional benchmark scores. "
            "A/B status and sensitivity endpoint/check statuses are reported separately."
        ),
    }
    write_new(output / "benchmark_manifest.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-release", type=Path, required=True)
    parser.add_argument("--collection", type=Path, action="append", required=True)
    parser.add_argument("--sensitivity", type=Path, action="append", default=[],
                        help="optional nominal or half-step Stage C collection; retained separately")
    parser.add_argument("--output", type=Path, required=True, help="new directory under workspaces")
    args = parser.parse_args()
    result = export_ledger(args.parent_release, args.collection, args.output, args.sensitivity)["development_export"]
    print(json.dumps({key: result[key] for key in (
        "status", "record_count_preserved", "development_pair_count", "completed_endpoints", "total_endpoints")}))


if __name__ == "__main__":
    main()
