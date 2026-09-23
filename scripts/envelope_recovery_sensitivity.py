"""Join the qualified two-start diagnostic separately from the frozen full100 result."""
import argparse
import copy
import json
from collections import Counter

from affordable_common import InvalidArtifact, HA_TO_KCAL, read_json, record, verify, write_new
from motion_envelope_transfer_compare import aggregate, call
from nikasha_pool import choose_rows


PROTOCOL = "Nikasha_envelope100_single_cell_seeded_sensitivity_v1"


def require(condition, message):
    if not condition:
        raise InvalidArtifact(message)


def compare(primary, recovery):
    old, new = read_json(primary), read_json(recovery)
    require(old["denominators"]["triples"] == 100 and len(old["rows"]) == 104, "Expected original104/100 ledger")
    require(new["status"] == "qualified_recovery_sensitivity" and new["agreement_pass"] and
            new["reported_seed"] == "origin" and new["primary_result_changed"] is False,
            "Qualified fixed-origin recovery required")
    require(new["denominator"] == len(new["rows"]) == 2, "Exactly two declared starts required")
    starts = {r["seed_kind"]: r for r in new["rows"]}
    require(set(starts) == {"origin", "adaptive_La"} and
            all(r["status"] == "confirmed_restart" for r in starts.values()), "Both real native restarts required")
    delta = abs(starts["origin"]["energy_hartree"] - starts["adaptive_La"]["energy_hartree"]) * HA_TO_KCAL
    require(delta <= 0.1 and abs(delta - new["seed_difference_kcal_mol"]) < 1e-12,
            "Two-start agreement failed")
    require(new["recovery_energy_hartree"] == starts["origin"]["energy_hartree"], "Predeclared seed changed")
    manifest = read_json(verify(new["manifest"]))
    inventory = read_json(verify(manifest["inventory"]))
    target = inventory["target"]
    for start in starts.values():
        for key in ("output", "receipt", "manifest"):
            verify(start["actual"][key])
    for task in manifest["tasks"]:
        require(verify(task["xyz"]).read_bytes() == verify(target["task"]["xyz"]).read_bytes(), "Recovery geometry differs")
        require((task["charge"], task["multiplicity"], task["metal"], task["medium"]) ==
                (target["task"]["charge"], 1, "La", "vacuum"), "Recovery state differs")
    candidates = [r for r in old["rows"] if r["status"] == "unavailable" and r["matrix"] and
                  r["matrix"]["La"]["adaptive_Ca"]["xyz"]["sha256"] == target["task"]["xyz"]["sha256"]]
    require(len(candidates) == 1, "Exactly one matching failed primary source required")
    failed = candidates[0]
    matrix = copy.deepcopy(failed["matrix"])
    cell = matrix["La"]["adaptive_Ca"]
    require(cell["status"] == "unavailable" and cell["low"]["vacuum"]["reason"] == "SCF_not_converged" and
            cell["low"]["alpb"]["status"] == "complete", "Expected the actual failed vacuum cell")
    mace = read_json(verify(cell["MACE"]))
    require(mace["status"] == "complete", "Archived MACE component unavailable")
    verify(cell["low"]["alpb"]["output"])
    cell.update(status="complete", components={
        "MACE_eV": mace["energy_eV"],
        "GFN2_ALPB_hartree": cell["low"]["alpb"]["energy_hartree"],
        "GFN2_vacuum_hartree": new["recovery_energy_hartree"],
    }, numerical_recovery=record(recovery))
    cell["low"]["vacuum"] = {**starts["origin"]["actual"], "status": "complete",
                                "initialization": "predeclared_origin_seed_after_two_start_agreement"}
    pool = choose_rows(matrix, ["origin", "adaptive_Ca", "adaptive_La"])
    require(pool["status"] == "available", "Recovered common pool remains unavailable")
    ref = read_json(verify(old["reference"]))
    bands = old["bands"]["envelope_operational"]
    require(bands == ref["variants"]["operational"]["bands"], "Frozen bands differ")
    selected = call(pool["operational"]["composite_R_model_kcal_mol"], bands, failed["expected_class"])
    by_pair, rows = {}, []
    for row in old["rows"]:
        recovered = row["case_id"] == failed["case_id"] and row["selection_id"] == failed["selection_id"]
        value = selected if recovered else row["methods"]["envelope_operational"]
        entry = {"case_id": row["case_id"], "selection_id": row["selection_id"],
                 "primary": row["methods"]["envelope_operational"], "sensitivity": value,
                 "seeded_recovery_used": recovered}
        rows.append(entry)
        by_pair[(row["selection_id"], row["case_id"])] = value
    triples = []
    for triple in old["triples"]:
        members = [by_pair.get((triple["selection_id"], cid), call(None, bands, triple["expected_class"]))
                   for cid in triple["members"]]
        sensitivity = aggregate(members, bands, triple["expected_class"], triple["status"] != "prepared")
        triples.append({"triple_id": triple["triple_id"], "protein_id": triple["protein_id"],
                        "members": triple["members"], "primary": triple["methods"]["envelope_operational"],
                        "sensitivity": sensitivity})
    return {"protocol_id": PROTOCOL, "primary": record(primary), "recovery": record(recovery),
            "reference": old["reference"], "bands": bands, "rows": rows, "triples": triples,
            "changed_source": {"case_id": failed["case_id"], "selection_id": failed["selection_id"],
                               "matrix": matrix, "pool": pool},
            "counts": {name: {variant: dict(Counter(r[variant]["outcome"] for r in entries))
                              for variant in ("primary", "sensitivity")}
                       for name, entries in (("pairs", rows), ("triples", triples))},
            "new_molecular_calls": 0, "new_reference_fitted": False,
            "primary_result_changed": False, "production_changed": False,
            "interpretation": "Separate post-failure numerical sensitivity; original primary unavailability retained."
                              " No general automatic restart policy is implemented or qualified here.",
            "implementation": record(__file__)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("primary", "recovery", "output"):
        parser.add_argument("--" + key, required=True)
    args = vars(parser.parse_args())
    output = args.pop("output")
    result = compare(**args)
    write_new(output, result)
    print(json.dumps({"output": record(output), "counts": result["counts"]}, indent=2))
