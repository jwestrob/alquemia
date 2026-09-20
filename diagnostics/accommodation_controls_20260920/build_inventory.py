#!/usr/bin/env python3
"""Export source mappings for consumed accommodation controls; no energy calls.

Run with the pinned lanm_qmmm Python from any directory. Inputs are existing
immutable preparations. Output contains geometry inventories, not new QM inputs.
"""
from pathlib import Path
import hashlib
import json
import re
from decimal import Decimal
import numpy as np
import gemmi

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "workspaces/accommodation_controls_20260920"


def pin(path):
    p = Path(path).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def write(name, data):
    p = OUT / name
    p.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    return pin(p)


def atomkey(s):
    return (s["chain"], s["resnum"], s.get("insertion_code", ""),
            s["resname"], s["atom"])


def mapkey(row):
    if row["kind"] == "source":
        return ("source", atomkey(row["source"]))
    return (row["kind"], atomkey(row["retained"]), atomkey(row["omitted"]))


def source(path):
    b = gemmi.cif.read_file(str(path)).sole_block()
    s = gemmi.make_structure_from_block(b)
    # Do not use model['A']: mmCIF can split author chain A into polymer,
    # ion and water chain objects with the same author identifier.
    atoms, alternates = {}, []
    protein_residues = {}
    for chain in s[0]:
        if chain.name != "A":
            continue
        for res in chain:
            if res.het_flag == "A":
                protein_residues[res.seqid.num] = res.name
            for a in res:
                if a.element.is_hydrogen:
                    continue
                key = (chain.name, res.seqid.num, res.seqid.icode.strip(), res.name, a.name)
                row = {"source_key": list(key), "element": a.element.name,
                       "xyz_A": list(a.pos), "occupancy": a.occ,
                       "altloc": a.altloc.strip("\x00 "), "protein": res.het_flag == "A"}
                if key in atoms:
                    alternates.append({"key": list(key), "rows": [atoms[key], row]})
                    # Deterministic inventory selection; record every ambiguity.
                    row = sorted([atoms[key], row], key=lambda x: (-x["occupancy"], x["altloc"]))[0]
                atoms[key] = row
    return {"pin": pin(path), "sequence": "".join(gemmi.cif.as_string(
                b.find_values("_entity_poly.pdbx_seq_one_letter_code_can")[0]).split()),
            "atoms": atoms, "observed_residues": protein_residues,
            "alternates": alternates,
            "experiment": list(b.find_values("_exptl.method")),
            "resolution_A": s.resolution,
            "crystallization_details": [gemmi.cif.as_string(x) for x in b.find_values("_exptl_crystal_grow.details")],
            "crystallization_pH": list(b.find_values("_exptl_crystal_grow.pH"))}


def load_xyz(path):
    lines = Path(path).read_text().splitlines()
    rows = [x.split() for x in lines[2:] if x.strip()]
    assert len(rows) == int(lines[0])
    return [x[0] for x in rows], np.array([[float(t) for t in x[1:4]] for x in rows])


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    sources, repairs, cases = {}, {}, {}
    for code in ("1F6S", "6IP9"):
        sp = ROOT / f"diagnostics/laca_benchmark_expansion_20260915/sources/pdb/{code}.cif"
        mp = ROOT / f"workspaces/benchmark_set_20260915/prepared/alacta_{code.lower()}_v2/strong_site/amide_v3/repair_manifest.json"
        sources[code] = source(sp)
        d = json.loads(mp.read_text())
        repairs[code] = d
        arms = {}
        for metal in ("Ca", "La"):
            record = d["outputs"][metal]
            xp = record["xyz"]["path"]
            assert pin(xp) == record["xyz"]
            elements, coords = load_xyz(xp)
            assert elements[0] == metal
            for a in d["atom_graph"]["source_to_qm"]:
                assert np.max(np.abs(coords[a["qm_index"]] - a["xyz_A"])) < 1e-7
            arms[metal] = {"xyz": pin(xp), "n_atoms": len(elements),
                           "charge": record["charge"], "multiplicity": record["multiplicity"]}
        ca = load_xyz(arms["Ca"]["xyz"]["path"])
        la = load_xyz(arms["La"]["xyz"]["path"])
        assert ca[0][1:] == la[0][1:] and np.array_equal(ca[1], la[1])
        heavy = []
        for a in d["atom_graph"]["source_to_qm"]:
            if a["kind"] != "source" or a["source"]["element"] == "H":
                continue
            original = sources[code]["atoms"][atomkey(a["source"])]
            assert np.array_equal(a["xyz_A"], original["xyz_A"])
            heavy.append({"qm_index_zero_based": a["qm_index"], **original})
        cases[code] = {"source": pin(sp), "repair_manifest": pin(mp),
                       "protocol_id": d["protocol_id"], "native_metal": d["selected_site"],
                       "endpoints": arms, "water_inventory": d["explicit_water_inventory"],
                       "charge_ledger": d["charge_ledger"], "source_heavy_atoms": heavy,
                       "observed_chain_A_residues": sources[code]["observed_residues"],
                       "alternates": sources[code]["alternates"],
                       "crystallization_details": sources[code]["crystallization_details"],
                       "crystallization_pH": sources[code]["crystallization_pH"],
                       "coordination_annotation": d["coordination"],
                       "exposure": "previously_scored_consumed_development", "new_scientific_evaluations": 0}
    assert sources["1F6S"]["sequence"] == sources["6IP9"]["sequence"]
    common_ca = sorted(k for k, v in sources["1F6S"]["atoms"].items()
                       if v["protein"] and k[-1] == "CA" and k in sources["6IP9"]["atoms"])
    x = np.array([sources["6IP9"]["atoms"][k]["xyz_A"] for k in common_ca])
    y = np.array([sources["1F6S"]["atoms"][k]["xyz_A"] for k in common_ca])
    u, _, vt = np.linalg.svd((x-x.mean(0)).T @ (y-y.mean(0)))
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    translation = y.mean(0) - x.mean(0) @ rotation
    aligned = x @ rotation + translation
    assert np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-12)
    assert np.linalg.det(rotation) > 0
    alignment = {"moving": "6IP9", "reference": "1F6S", "selection": "all common source chain-A C-alpha atoms; no core or water optimization",
                 "row_vector_expression": "x_1F6S_frame = x_6IP9 @ rotation + translation",
                 "rotation": rotation.tolist(), "translation_A": translation.tolist(),
                 "n_anchors": len(common_ca), "anchor_keys": [list(k) for k in common_ca],
                 "rmsd_A": float(np.sqrt(np.mean(np.sum((aligned-y)**2, axis=1))))}
    maps = {}
    for code, d in repairs.items():
        maps[code] = {mapkey(a): a for a in d["atom_graph"]["source_to_qm"]
                      if not (a["kind"] == "source" and a["source"]["resname"] == "HOH")}
    assert maps["1F6S"].keys() == maps["6IP9"].keys()
    common = []
    for key in maps["1F6S"]:
        a, b = maps["1F6S"][key], maps["6IP9"][key]
        common.append({"kind": a["kind"], "identity": key,
                       "1F6S_qm_index_zero_based": a["qm_index"], "6IP9_qm_index_zero_based": b["qm_index"],
                       "1F6S_xyz_A": a["xyz_A"], "6IP9_xyz_A": b["xyz_A"],
                       "6IP9_aligned_xyz_A": (np.array(b["xyz_A"]) @ rotation+translation).tolist(),
                       "source_heavy": a["kind"] == "source" and a["source"]["element"] != "H"})
    assert sum(a["source_heavy"] for a in common) == 18
    assert len(common) == 33
    water_positions = {}
    for code in cases:
        water_positions[code] = [a for a in cases[code]["source_heavy_atoms"] if a["source_key"][3] == "HOH"]
    water_distance_matrix = [[float(np.linalg.norm(np.array(a["xyz_A"]) - (np.array(b["xyz_A"]) @ rotation+translation)))
                             for b in water_positions["6IP9"]] for a in water_positions["1F6S"]]
    mapping = {"schema_version": "alquemia.accommodation_atom_mapping.v1", "sequence": sources["1F6S"]["sequence"],
               "sequence_length": len(sources["1F6S"]["sequence"]), "alignment": alignment,
               "common_nonmetal_core_atoms": common, "common_protein_heavy_atom_count": 18,
               "common_nonmetal_core_atom_count_including_prepared_H_and_caps": 33,
               "metal_qm_index_zero_based": 0,
               "metal_sites": {k: d["selected_site"] for k, d in repairs.items()},
               "water_positions": water_positions, "water_distance_matrix_A_rows_1F6S_cols_6IP9": water_distance_matrix,
               "water_identity_claim": "none; matrix is only positional proximity after protein alignment",
               "symmetry_warning": "ASP OD1/OD2 are equivalent labels for unprotonated carboxylates; names are traceability, not unique physical displacement identities",
               "hydrogen_warning": "hydrogens and caps were prepared, not resolved experimental atoms; historical water-H coordinates require current versioned water preparation before new scoring",
               "new_scientific_evaluations": 0}
    write("alpha_sources.json", cases)
    write("alpha_atom_mapping.json", mapping)
    # Reuse the already archived, explicitly water-deleted geometry. This is
    # an intervention, not a claim that both native crystals contain 2 waters.
    matched = {}
    for state in ("1F6S__full", "6IP9__minus_A_322"):
        folder = ROOT / "workspaces/hydration_square_20260918/repaired_v1" / state
        preparation = json.loads((folder / "preparation.json").read_text())
        code = preparation["case"]
        assert preparation["parent_preparation"] == cases[code]["repair_manifest"]
        arm_records = {}
        for metal in ("Ca", "La"):
            xp = folder / metal / "core.xyz"
            el, xyz = load_xyz(xp)
            old_el, old_xyz = load_xyz(cases[code]["endpoints"][metal]["xyz"]["path"])
            assert len(el) == 40
            for new_i, old_i in enumerate(preparation["retained_parent_indices"]):
                assert el[new_i] == old_el[old_i]
                if el[new_i] != "H":
                    assert np.array_equal(xyz[new_i], old_xyz[old_i])
            arm_records[metal] = {"xyz": pin(xp), "n_atoms": len(el),
                                  "charge": cases[code]["endpoints"][metal]["charge"], "multiplicity": 1}
        assert np.array_equal(load_xyz(arm_records["Ca"]["xyz"]["path"])[1], load_xyz(arm_records["La"]["xyz"]["path"])[1])
        matched[state] = {"preparation": pin(folder / "preparation.json"),
                          "endpoints": arm_records, "water_geometry_policy": preparation["water_geometry_policy"],
                          "removed_water": preparation["removed_water"],
                          "source_mapping": preparation["source_mapping"],
                          "exposure": "all_single_water_deletions_already_scored_development"}
    write("alpha_equal_two_water_pair.json", {
        "status": "existing_prepared_intervention_pair_not_new_native_pair", "states": matched,
        "composition": "same two peptide-amide units plus three ASP sidechain fragments, metal, two neutral waters; same Ca -1 / La 0 charge",
        "water_oxygen_correspondence": [{"1F6S_source": "A:HOH211/O", "6IP9_source": "A:HOH326/O", "1F6S_prepared_index0": 34, "6IP9_prepared_index0": 37, "distance_after_alignment_A": water_distance_matrix[0][2]},
                                        {"1F6S_source": "A:HOH212/O", "6IP9_source": "A:HOH310/O", "1F6S_prepared_index0": 37, "6IP9_prepared_index0": 34, "distance_after_alignment_A": water_distance_matrix[1][0]}],
        "water_correspondence_basis": "nearest distinct oxygen positions under all-common-C-alpha alignment, without energy or score input",
        "caveats": ["6IP9 water322 was explicitly deleted in an earlier experiment; do not call the result native La-bound geometry",
                    "retained heavy coordinates unchanged; hydrogen geometry follows archived hydration-square normalization, not latest contextual water policy",
                    "hydrogen labels of distinct water molecules do not specify physical identity or a unique displacement path",
                    "equal composition controls water number but does not isolate metal causality from crystal packing or donor/water orientation"],
        "new_energy_calls": 0})
    # Parent's later explicit scope extension: reuse the four compatible
    # archived energies for this exact intervention matrix; never rescore.
    collection_path = ROOT / "workspaces/hydration_square_20260918/repaired_v1/collection_1201801.json"
    collection = json.loads(collection_path.read_text())
    energies, endpoint_records = {}, {}
    for state in matched:
        row = next(r for r in collection["rows"] if r["state_id"] == state)
        assert row["status"] == "complete" and not row["failures"]
        energies[state], endpoint_records[state] = {}, {}
        for metal in ("Ca", "La"):
            ep = row["endpoints"][metal]
            assert pin(ep["output"]["path"]) == ep["output"]
            assert pin(ep["receipt"]["path"]) == ep["receipt"]
            receipt = json.loads(Path(ep["receipt"]["path"]).read_text())
            assert receipt["normal_termination"] and receipt["scf_converged"] and receipt["returncode"] == 0
            assert receipt["orca_version"] == "6.1.1"
            for artifact in receipt["artifacts"].values():
                assert pin(artifact["path"]) == artifact
            assert receipt["artifacts"]["xyz"] == matched[state]["endpoints"][metal]["xyz"]
            template = Path(receipt["artifacts"]["template_input"]["path"]).read_text()
            assert template.splitlines()[0] == "! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3"
            text = Path(ep["output"]["path"]).read_text()
            assert "ORCA TERMINATED NORMALLY" in text
            final = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
            assert len(final) == 1
            energy = Decimal(final[0])
            assert abs(float(energy)-ep["energy_hartree"]) < 1e-12
            energies[state][metal] = energy
            endpoint_records[state][metal] = {**ep, "energy_hartree_decimal": str(energy),
                                              "xyz": receipt["artifacts"]["xyz"], "input": receipt["artifacts"]["template_input"]}
    ca_geometry, la_geometry = "1F6S__full", "6IP9__minus_A_322"
    conversion = Decimal("627.509474")
    contrasts = {s: v["Ca"]-v["La"] for s,v in energies.items()}
    shifts = {m: energies[la_geometry][m]-energies[ca_geometry][m] for m in ("Ca", "La")}
    delta_r = (contrasts[la_geometry]-contrasts[ca_geometry])*conversion
    assert delta_r == (shifts["Ca"]-shifts["La"])*conversion
    energy_matrix = {"status": "compatible_archived_matrix", "collection": pin(collection_path),
                     "protocol_id": collection["protocol_id"], "method": "native ORCA6.1.1 r2SCAN-3c/CPCM(Water)/DefGrid3; archived NormalSCF default; singlets; Ca -1 / La 0",
                     "water_geometry_policy": "reference_internal_geometry; identical in all four endpoints",
                     "hartree_to_kcal_mol": str(conversion), "endpoints": endpoint_records,
                     "R_ECa_minus_ELa_hartree_decimal": {k: str(v) for k,v in contrasts.items()},
                     "geometry_shift_ELaboundsource_minus_ECaboundsource_kcal_mol": {k: str(v*conversion) for k,v in shifts.items()},
                     "delta_R_kcal_mol_decimal": str(delta_r),
                     "algebra": "[ECa(6IP9 minus322)-ELa(6IP9 minus322)]-[ECa(1F6S full)-ELa(1F6S full)]",
                     "interpretation": "negative shifts relative contrast toward Ca at the explicitly two-water La-source geometry; not an affinity classifier or a causal metal-only effect",
                     "absolute_affinity_or_population": None, "new_endpoint_evaluations": 0,
                     "historical_selected_endpoint_wall_seconds_sum": sum(ep["wall_seconds"] for row in endpoint_records.values() for ep in row.values()),
                     "baseline_changed": False}
    write("alpha_archived_dft_matrix.json", energy_matrix)
    (OUT / "alpha_common_heavy.tsv").write_text("source_atom\t1F6S_qm_index0\t6IP9_qm_index0\t1F6S_xyz_A\t6IP9_xyz_A\n" + "".join(
        f"{r['identity'][1]}\t{r['1F6S_qm_index_zero_based']}\t{r['6IP9_qm_index_zero_based']}\t{r['1F6S_xyz_A']}\t{r['6IP9_xyz_A']}\n"
        for r in common if r["source_heavy"]))
    print(json.dumps({"sequence_length": mapping["sequence_length"], "observed_residue_counts": {k: len(v["observed_chain_A_residues"]) for k,v in cases.items()},
                      "common_heavy_atoms": 18, "common_nonmetal_core_atoms": len(common), "alignment_rmsd_A": alignment["rmsd_A"],
                      "water_distance_matrix_A": water_distance_matrix, "checks": "passed", "new_energy_calls": 0}, indent=2))


if __name__ == "__main__":
    main()
