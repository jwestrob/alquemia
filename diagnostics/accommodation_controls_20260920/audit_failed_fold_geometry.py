#!/usr/bin/env python3
"""Read existing failed preparations only; do not prepare, repair or score."""
from __future__ import annotations
import argparse
from collections import Counter
import datetime as dt
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree


def pin(path):
    p = Path(path).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def atoms(path):
    return [{"chain": c.name, "resname": r.name, "resnum": r.seqid.num,
             "icode": r.seqid.icode.strip(), "atom": a.name,
             "element": a.element.name, "xyz_A": list(a.pos)}
            for c in gemmi.read_structure(str(path))[0] for r in c for a in r]


def distance(a, b):
    return float(np.linalg.norm(np.asarray(a["xyz_A"]) - b["xyz_A"]))


def closest(center, candidates, count=1):
    return sorted([{"atom": a, "distance_A": distance(center, a)}
                   for a in candidates], key=lambda x: x["distance_A"])[:count]


def overlaps(rows, cutoff=.45):
    xyz = np.asarray([a["xyz_A"] for a in rows])
    pairs = cKDTree(xyz).query_pairs(cutoff)
    return sorted([{"a": rows[i], "b": rows[j], "distance_A": distance(rows[i], rows[j])}
                   for i, j in pairs], key=lambda x: x["distance_A"])


def c5_case(case, prepared):
    source = case["source_structure"]["path"]
    if pin(source)["sha256"] != case["source_structure"]["sha256"]:
        raise ValueError("frozen source hash changed")
    rows = atoms(source)
    metal = next(a for a in rows if a["element"] in ("La", "Ca"))
    pqq = [a for a in rows if a["chain"] == "C"]
    protein = [a for a in rows if a["chain"] == "A"]
    result = {"case_id": case["case_id"], "conditioning": case["source_conditioning_metal"],
              "source": pin(source), "metal": metal,
              "PQQ_atom_count": len(pqq), "PQQ_nearest": closest(metal, pqq),
              "PQQ_named_donor_distances_A": {a["atom"]: distance(metal, a)
                  for a in pqq if a["atom"] in ("O1", "N2", "O5", "O6")},
              "nearest_protein_ON": closest(metal, [a for a in protein
                  if a["element"] in ("O", "N")], 10), "roles": {}}
    for role, selector in case["roles"].items():
        selected = [a for a in protein if a["resnum"] == selector["resnum"]
                    and a["icode"] == selector["icode"] and a["element"] in ("O", "N")]
        result["roles"][role] = closest(metal, selected, 2)
    for filename in ("normalized.pdb", "protonated.pdb"):
        p = prepared/"cases"/case["case_id"]/filename
        if p.exists():
            m = next(a for a in atoms(p) if a["element"] in ("La", "Ca"))
            result[filename] = {"file": pin(p), "metal": m,
                                "metal_displacement_from_raw_A": distance(metal, m)}
    failure = prepared/"failures"/(case["case_id"]+".json")
    if failure.exists():
        result["failure"] = {"file": pin(failure), "data": json.loads(failure.read_text())}
    return result


def mmol_case(case, prepared):
    folder = prepared/"cases"/case["case_id"]
    cp = folder/"core_preparation.json"
    core = json.loads(cp.read_text())
    raw = atoms(case["source_structure"]["path"])
    protonated = atoms(folder/"protonated.pdb")
    qm = [{"qm_index": 0, "fragment": "metal", "element": "La",
           "origin": "source_heavy_atom", "xyz_A": case["raw_source_metal"]["xyz_A"]}]
    for frag in core["qm_fragments"]:
        for a in frag["atom_records"]:
            qm.append(dict(a, fragment=frag["id"], qm_index=len(qm)))
    raw_heavy = {(a["chain"], a["resnum"], a["icode"], a["atom"]): a
                 for a in raw if a["element"] not in ("H", "D")}
    prep_heavy = {(a["chain"], a["resnum"], a["icode"], a["atom"]): a
                  for a in protonated if a["element"] not in ("H", "D")}
    shared = set(raw_heavy) & set(prep_heavy)
    pm = folder/"protonated_protonation_manifest.json"
    pairs = overlaps(qm)
    return {"case_id": case["case_id"], "source": pin(case["source_structure"]["path"]),
            "core_manifest": pin(cp), "protonation_manifest": pin(pm),
            "protonation": json.loads(pm.read_text()),
            "source_H_count": sum(a["element"] == "H" for a in raw),
            "core_atom_count": len(qm), "core_overlap_pairs_below_0p45A": pairs,
            "core_overlap_pair_origin_counts": dict(Counter(" / ".join(sorted((p["a"]["origin"], p["b"]["origin"]))) for p in pairs)),
            "protonated_source": pin(folder/"protonated.pdb"),
            "protonated_overlap_pairs_below_0p45A": overlaps(protonated),
            "raw_heavy_overlap_pairs_below_0p45A": overlaps(list(raw_heavy.values())),
            "shared_heavy_atom_count": len(shared),
            "unmatched_raw_heavy_keys": sorted(set(raw_heavy)-set(prep_heavy)),
            "unmatched_prepared_heavy_keys": sorted(set(prep_heavy)-set(raw_heavy)),
            "max_shared_heavy_displacement_A": max(distance(raw_heavy[k],prep_heavy[k]) for k in shared)}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--sources", type=Path, required=True)
    p.add_argument("--prepared", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args=p.parse_args()
    cases=json.loads(args.sources.read_text())["cases"]
    c5=[c5_case(c,args.prepared) for c in cases if c["root_case_id"]=="c5axv8-pqq-la_model"]
    mmol=next(c for c in cases if c["case_id"]=="mmol_1770-pqq-la_model__conditioned_La__seed-1_sample-4")
    result={"scope":"read_only_geometry_of_existing_failures_no_preparation_or_energy_calls",
            "created_at":dt.datetime.now(dt.timezone.utc).isoformat(),
            "frozen_population":pin(args.sources), "population_denominator":len(cases),
            "PQQ_detection_cutoff_A":4.0, "expanded_overlap_cutoff_A":.45,
            "C5AXV8":c5, "MMOL1770_La_sample4":mmol_case(mmol,args.prepared)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x") as f: json.dump(result,f,indent=2);f.write("\n")
    for row in c5:
        print(row["case_id"], "nearest PQQ",row["PQQ_nearest"][0]["distance_A"],
              "nearest protein",row["nearest_protein_ON"][0])
    m=result["MMOL1770_La_sample4"]
    print("MMOL overlaps",len(m["core_overlap_pairs_below_0p45A"]),m["core_overlap_pair_origin_counts"])
    print("MMOL whole protein overlaps",len(m["protonated_overlap_pairs_below_0p45A"]),
          "raw heavy overlaps",len(m["raw_heavy_overlap_pairs_below_0p45A"]),
          "heavy max displacement",m["max_shared_heavy_displacement_A"])


if __name__=="__main__": main()
