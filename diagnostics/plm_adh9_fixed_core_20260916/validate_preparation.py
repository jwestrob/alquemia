#!/usr/bin/env python3
"""Read-only validation of the two prepared pairs and four explicit exclusions."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("plm_candidate_wrapper", HERE / "prepare_candidates.py")
wrapper = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = wrapper
spec.loader.exec_module(wrapper)
OUT = wrapper.ROOT / "workspaces/plm_adh9_fixed_core_20260916"


def check(condition, reason):
    if not condition:
        raise AssertionError(reason)


def run():
    pins, original, cp = wrapper.verify_pins(HERE / "implementation_pins.json")
    inventory = wrapper.read(OUT / "prepared_pairs.json")
    check(len(inventory["cases"]) == 6, "Six selected models must remain visible")
    check(inventory["prepared_pair_count"] == 2 and inventory["unsupported_count"] == 4,
          "Unexpected preparation totals")
    checks = []
    for case in inventory["cases"]:
        wrapper.verify(case["source_cif"])
        if case["status"] != "ready_for_orca":
            check(case["carve_manifest"] is None and case["reason"], "Unsupported case must be unscored with reason")
            checks.append({"case_id": case["case_id"], "status": "unsupported", "reason": case["reason"]})
            continue
        mp = wrapper.verify(case["carve_manifest"])
        m = wrapper.read(mp)
        check(m["protocol_id"] == wrapper.PROTOCOL, "Wrong scientific protocol")
        for name in ("source_cif", "source_structure", "normalized_source_structure", "normalization_manifest", "protonation_manifest"):
            wrapper.verify(m[name])
        h = m["heavy_coordinate_check"]
        hp = Path(h["path"])
        check(hashlib.sha256(hp.read_bytes()).hexdigest() == h["sha256"], "Heavy-check hash differs")
        hj = wrapper.read(hp)
        check(hj["no_heavy_atoms_added"] and hj["all_source_heavy_atoms_retained"], "Heavy changes")
        for item in hj["checks"].values():
            check(item["passes"] and item["maximum_displacement_A"] <= 0.001, "Heavy-coordinate failure")
        pm = wrapper.read(m["protonation_manifest"]["path"])
        check(pm["repaired_missing_atom_count"] == pm["repaired_missing_terminal_atom_count"] == 0,
              "No heavy-atom repair may occur")
        check(m["charge_ledger"]["La_total"] == -2 and m["charge_ledger"]["Ca_total"] == -3,
              "Acidic-D+2 charge ledger mismatch")
        check(m["confidence"]["value"] >= .9 and m["coordination"]["coordination_number"] >= 7,
              "Production admission gate failure")
        check(m["fixed_core"]["catalytic_Asp_partner_geometry"]["maximum_allowed_A"] == 3.5,
              "Partner gate altered")
        xyzs = {}
        for metal in ("La", "Ca"):
            for suffix in ("xyz", "input"):
                rr = m["outputs"][metal + "_" + suffix]
                path = mp.parent / rr["path"]
                check(hashlib.sha256(path.read_bytes()).hexdigest() == rr["sha256"], "Endpoint hash differs")
            inp = (mp.parent / m["outputs"][metal + "_input"]["path"]).read_text()
            check("! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3" in inp, "Electronic method differs")
            check("%pal" not in inp.lower() and "newgto" not in inp.lower() and "new_ecp" not in inp.lower(),
                  "Unexpected resource/basis override")
            xyzs[metal] = (mp.parent / m["outputs"][metal + "_xyz"]["path"]).read_text().splitlines()
            check(int(xyzs[metal][0]) == len(xyzs[metal]) - 2, "XYZ count differs")
        check(xyzs["La"][3:] == xyzs["Ca"][3:], "Nonmetal coordinates differ")
        check(xyzs["La"][2].split()[1:] == xyzs["Ca"][2].split()[1:], "Metal coordinates differ")
        # Verify every final nonmetal atom against the retained provenance ledger,
        # independently of the atom-count and byte-equality assertions in preparation.
        expected_atoms = [a for f in m["qm_fragments"] for a in f["atom_records"]]
        check(len(expected_atoms) == len(xyzs["La"]) - 3, "Fragment/XYZ count differs")
        for line, atom in zip(xyzs["La"][3:], expected_atoms):
            parts = line.split()
            check(parts[0].upper() == atom["element"].upper(), "Fragment element differs")
            check([float(x) for x in parts[1:]] == atom["xyz_A"], "Fragment coordinate differs")
        pqq = m["qm_fragments"][0]
        check(pqq["atom_count"] == 27 and pqq["formal_charge"] == -3, "PQQ incomplete")
        # Hash rejection must fail closed before any preparation occurs.
        broken = dict(case["source_cif"], sha256="0" * 64)
        try:
            wrapper.verify(broken)
        except ValueError:
            pass
        else:
            raise AssertionError("Changed source hash was accepted")
        checks.append({"case_id": case["case_id"], "status": "PASS", "carve_manifest": case["carve_manifest"],
            "qm_atom_count": int(xyzs["La"][0]), "PQQ_atom_count": 27,
            "CN": m["coordination"]["coordination_number"], "confidence": m["confidence"]["value"],
            "raw_Asp_partner_A": m["fixed_core"]["raw_catalytic_Asp_partner_geometry"]["heavy_atom_distance_A"],
            "prepared_Asp_partner_A": m["fixed_core"]["catalytic_Asp_partner_geometry"]["heavy_atom_distance_A"]})
    check(not list(OUT.rglob("*.out")), "ORCA must not have run during preparation")
    receipt = {"status": "PASS", "checks": checks, "prepared_pair_count": 2, "unsupported_count": 4,
        "candidate_manifest": pins["candidate_manifest"], "wrapper": pins["wrapper"],
        "validation_script": wrapper.record(__file__), "orca_executed": False,
        "scientific_helpers_match_original_calibration_pins": True}
    wrapper.write(HERE / "preparation_validation.json", receipt)
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    run()
