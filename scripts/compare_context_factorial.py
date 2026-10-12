"""Collect the real two-source context x saved-geometry factorial; no calibration."""
import argparse
import csv
import json
from pathlib import Path
from context_factorial_mace import EV_TO_KCAL, pin, read, verify, write, xyz

HA_TO_KCAL = 627.509474  # retain the released protocol's conversion exactly
ORDER = ["origin", "original_adaptive_Ca", "original_adaptive_La",
         "expanded_adaptive_Ca", "expanded_adaptive_La"]
POOLS = {"origin": ["origin"], "original_three": ORDER[:3],
         "expanded_three": [ORDER[0], *ORDER[3:]], "common_five": ORDER}


def collect(mace_path, solvent_path):
    mace, solvent = read(mace_path), read(solvent_path)
    if mace["status"] != "complete" or solvent["protocol_id"] != "Linux_xtb671_tblite070_component_factorial_v2":
        raise ValueError("unexpected/incomplete calculation protocol")
    mrows = {r["task_id"]: r for r in mace["rows"]}
    endpoints = []
    for low in solvent["rows"]:
        high = mrows[low["id"]]
        verify(high["receipt"])
        verify(high["xyz"])
        a, b = xyz(verify(high["xyz"])), xyz(verify(low["xyz"]))
        if len(a) != len(b) or [r[0] for r in a] != [r[0] for r in b]:
            raise ValueError("mismatched composite atom membership/order")
        coordinate_error = max(abs(x-y) for u, v in zip(a, b) for x, y in zip(u[1:], v[1:]))
        # Archive files and replay files can differ in comments/decimal spelling.
        # Use the same predeclared 1e-12-A serialization gate as MACE archive reuse.
        if coordinate_error > 1e-12 or high["charge"] != low["charge"] or high["multiplicity"] != low["uhf"] + 1:
            raise ValueError("mismatched composite coordinates or state")
        for receipt in low["inputs"].values():
            if receipt:
                verify(receipt)
        r = {k: high[k] for k in ("task_id", "protein", "target_context", "geometry_id", "metal", "charge", "multiplicity", "xyz")}
        accepted = low["accepted"]
        r.update(accepted=accepted, failure_reason=low.get("failure_reason"),
                 solvent_xyz=low["xyz"], composite_coordinate_max_error_A=coordinate_error,
                 mace_eV=high["energy_eV"], mace_reused=high["reused"],
                 gfn2_vacuum_hartree=low["primary_vacuum_energy_hartree"],
                 gfn2_water_hartree=low["primary_water_energy_hartree"],
                 gfn2_solvation_hartree=low["primary_solvation_hartree"],
                 E_model_kcal_mol=(high["energy_eV"] * EV_TO_KCAL + low["primary_solvation_hartree"] * HA_TO_KCAL) if accepted else None,
                 solver_component_error_hartree=low.get("component_error_hartree"),
                 all_atom_population_error_e=low.get("max_population_error_e"))
        endpoints.append(r)
    if len(endpoints) != 40 or len(mrows) != 40:
        raise ValueError("declared forty-endpoint denominator differs")
    lookup = {(r["protein"], r["target_context"], r["geometry_id"], r["metal"]): r for r in endpoints}
    proteins = sorted(set(r["protein"] for r in endpoints))
    contrasts, pools = [], []
    for protein in proteins:
        for context in ("original", "expanded"):
            for geometry in ORDER:
                ca, la = [lookup[protein, context, geometry, metal] for metal in ("Ca", "La")]
                origins = {metal: lookup[protein, context, "origin", metal] for metal in ("Ca", "La")}
                accepted = ca["accepted"] and la["accepted"]
                r = {"protein": protein, "context": context, "geometry": geometry, "accepted": accepted,
                     "R_MACE_kcal_mol": (ca["mace_eV"] - la["mace_eV"]) * EV_TO_KCAL,
                     "R_solvent_kcal_mol": ((ca["gfn2_solvation_hartree"]-la["gfn2_solvation_hartree"]) * HA_TO_KCAL) if accepted else None,
                     "R_model_kcal_mol": ca["E_model_kcal_mol"] - la["E_model_kcal_mol"] if accepted else None}
                for metal, endpoint in (("Ca", ca), ("La", la)):
                    r["work_" + metal + "_kcal_mol"] = (endpoint["E_model_kcal_mol"]-origins[metal]["E_model_kcal_mol"]
                        if endpoint["accepted"] and origins[metal]["accepted"] else None)
                contrasts.append(r)
            for name, geometries in POOLS.items():
                available = all(lookup[protein, context, g, metal]["accepted"] for g in geometries for metal in ("Ca", "La"))
                r = {"protein": protein, "context": context, "pool": name, "geometry_order": geometries,
                     "available": available, "mathematical": None, "operational": None}
                if available:
                    for policy in ("mathematical", "operational"):
                        selected = {}
                        for metal in ("Ca", "La"):
                            best = min(geometries, key=lambda g: lookup[protein, context, g, metal]["E_model_kcal_mol"])
                            e0 = lookup[protein, context, "origin", metal]["E_model_kcal_mol"]
                            work = lookup[protein, context, best, metal]["E_model_kcal_mol"]-e0
                            if policy == "operational" and work >= -.1:
                                best, work = "origin", 0.0
                            selected[metal] = {"geometry": best, "E_model_kcal_mol": lookup[protein, context, best, metal]["E_model_kcal_mol"],
                                               "work_kcal_mol": work}
                        selected["R_model_kcal_mol"] = selected["Ca"]["E_model_kcal_mol"]-selected["La"]["E_model_kcal_mol"]
                        r[policy] = selected
                pools.append(r)
    cdict = {(r["protein"], r["context"], r["geometry"]): r for r in contrasts}
    effects = []
    for protein in proteins:
        for geometry in ORDER:
            r = {"protein": protein, "geometry": geometry}
            for component in ("MACE", "solvent", "model"):
                field = "R_" + component + "_kcal_mol"
                org, exp = [cdict[protein, context, geometry][field] for context in ("original", "expanded")]
                org0, exp0 = [cdict[protein, context, "origin"][field] for context in ("original", "expanded")]
                r[component + "_context_effect_kcal_mol"] = None if org is None or exp is None else exp-org
                r[component + "_motion_work_original_kcal_mol"] = None if org is None or org0 is None else org-org0
                r[component + "_motion_work_expanded_kcal_mol"] = None if exp is None or exp0 is None else exp-exp0
                r[component + "_interaction_kcal_mol"] = None if any(x is None for x in (org, exp, org0, exp0)) else (exp-exp0)-(org-org0)
            effects.append(r)
    anchors = []
    pdict = {(r["protein"], r["context"], r["pool"]): r for r in pools}
    for context in ("original", "expanded"):
        for pool in POOLS:
            ca = pdict["q9z4j7-pqq-la_model", context, pool]
            la = pdict["c5axv8-pqq-la_model", context, pool]
            anchors.append({"context": context, "pool": pool,
                "two_anchor_La_minus_Ca_kcal_mol": la["operational"]["R_model_kcal_mol"]-ca["operational"]["R_model_kcal_mol"] if ca["available"] and la["available"] else None,
                "is_full_canonical_gap": False})
    return {"protocol_id": "corrected_context_geometry_factorial_v2", "MACE": pin(mace_path), "solvent": pin(solvent_path),
            "implementation": pin(__file__), "declared_endpoints": 40, "accepted_endpoints": sum(r["accepted"] for r in endpoints),
            "units": "eV and hartree for original components; model kcal/mol after one conversion",
            "EV_TO_KCAL": EV_TO_KCAL, "HA_TO_KCAL": HA_TO_KCAL,
            "expression": "E=MACE_vac + GFN2_water-GFN2_vac; R=E_Ca-E_La",
            "bands": None, "classification_status": "no_calibration_for_two_development_cases", "population_claim": False,
            "endpoints": endpoints, "contrasts": contrasts, "pools": pools, "effects": effects, "anchor_ordering": anchors}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mace", required=True)
    p.add_argument("--solvent", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    result = collect(args.mace, args.solvent)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    write(out / "RESULT.json", result)
    for name in ("contrasts", "effects", "anchor_ordering"):
        with (out / (name + ".tsv")).open("x") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(result[name][0]), delimiter="\t")
            writer.writeheader(); writer.writerows(result[name])
    print(json.dumps({k: result[k] for k in ("declared_endpoints", "accepted_endpoints", "anchor_ordering")}, indent=2))


if __name__ == "__main__":
    main()
