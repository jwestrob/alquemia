"""Declared symmetric donor development objective; saved scalar/array arithmetic only.

The caller admits input/source provenance and passes the frozen rigid-body helper.
No file access, engines, optimizer, parameter renderer, or candidate selection.
"""
import math

FAMILIES = ("acetate_mono", "acetate_bi", "acetamide")
METALS = ("La", "Dy")
WATER = ("water_central", "aquo8_central")
DONOR = tuple("hydrated_" + f + "_" + s for f in FAMILIES
              for s in ("minus", "central", "plus")) + tuple(
                  "matched_dry_" + f + "_central" for f in FAMILIES)
CELLS = tuple(q + "_" + m for q in WATER + DONOR for m in METALS)
WATER_CELLS = tuple(q + "_" + m for q in WATER for m in METALS)
WEIGHTS = {"energy": .25, "La": .25, "Dy": .25, "paired": .25}
ENERGY_WEIGHTS = {"balanced_exchange": .50, "water_addition": .30,
                  "individual_radial_work": .20}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite_vector(vector, length):
    return (isinstance(vector, (list, tuple)) and len(vector) == length
            and all(isinstance(x, (int, float)) and not isinstance(x, bool)
                    and math.isfinite(x) for x in vector))


def check_values(values, atoms, label):
    energy = values.get("energy")
    require(isinstance(energy, (int, float)) and not isinstance(energy, bool)
            and math.isfinite(energy), label + ": nonfinite/absent energy")
    forces = values.get("forces")
    require(isinstance(forces, (list, tuple)) and len(forces) == atoms
            and all(finite_vector(v, 3) for v in forces),
            label + ": wrong/nonfinite force shape")


def check_inputs(data, model, body, water_baseline):
    rows = data["rows"]
    require(set(rows) == set(model) == set(CELLS),
            "Exactly the declared 28 data/model cells are required")
    require(set(water_baseline) == set(WATER_CELLS),
            "Explicit frozen baseline for all four water centers is required")
    for q in WATER + DONOR:
        pair = [rows[q + "_" + m] for m in METALS]
        for metal, row in zip(METALS, pair):
            cell = q + "_" + metal
            g = row["geometry"]
            n = g["atoms"]
            expected_charge = 2 if "acetate_" in q else 3
            expected_n = (4 if q == "water_central" else 25 if q == "aquo8_central"
                          else 10 if q.startswith("matched_dry_acetamide")
                          else 8 if q.startswith("matched_dry_")
                          else 31 if "acetamide" in q
                          else 26 if "acetate_bi" in q else 29)
            require(row["status"] == "complete" and g["geometry_id"] == cell
                    and g["metal"] == metal and g["metal_index"] == 0
                    and g["charge"] == expected_charge and n == expected_n
                    and g["physical_multiplicity"] == (1 if metal == "La" else 6),
                    cell + ": incomplete or wrong chemical state")
            require(len(g["positions_A"]) == n
                    and all(finite_vector(v, 3) for v in g["positions_A"])
                    and len(g["symbols"]) == n and g["symbols"][0] == metal,
                    cell + ": invalid coordinates/symbols")
            bonds = row["input"]["bonds_zero_based"]
            require(len(bonds) == n and all(
                isinstance(j, int) and not isinstance(j, bool) and 0 <= j < n
                and j != i and i in bonds[j]
                for i, neighbors in enumerate(bonds) for j in neighbors),
                cell + ": invalid/asymmetric covalent graph")
            require(row["bodies"] == body.bodies(row), cell + ": body mapping differs")
            require(sum(b["category"] == "ligand" for b in row["bodies"])
                    == (0 if q in WATER else 1), cell + ": ligand body missing/extra")
            check_values(row["reference"], n, cell + " reference")
            check_values(model[cell], n, cell + " model")
            if cell in WATER_CELLS:
                check_values(water_baseline[cell], n, cell + " frozen water")
                require(model[cell]["energy"] == water_baseline[cell]["energy"]
                        and model[cell]["forces"] == water_baseline[cell]["forces"],
                        cell + ": donor model changed frozen water E/forces")
        a, b = pair
        require(a["geometry"]["positions_A"] == b["geometry"]["positions_A"]
                and a["geometry"]["symbols"][1:] == b["geometry"]["symbols"][1:]
                and a["input"]["bonds_zero_based"] == b["input"]["bonds_zero_based"]
                and a["bodies"] == b["bodies"], q + ": paired coordinates/bodies differ")


def evaluate(data, model, body, water_baseline):
    """Score exact28 E[kcal/mol]/F[kcal/mol/A] rows against saved references.

    water_baseline must be supplied explicitly from the future frozen candidate;
    the historical training baseline is never selected automatically. This checks
    its four center records only. Full eight-water qualification is external.
    """
    check_inputs(data, model, body, water_baseline)
    rows = data["rows"]
    reference = {k: r["reference"] for k, r in rows.items()}

    def D(values, q):
        return values[q + "_Dy"]["energy"] - values[q + "_La"]["energy"]

    def exchange_error(q):
        return ((D(model, q) - D(model, "aquo8_central"))
                - (D(reference, q) - D(reference, "aquo8_central")))

    exchange = {q: exchange_error(q) for q in DONOR}
    addition = {f: ((D(model, "hydrated_" + f + "_central")
                    - D(model, "matched_dry_" + f + "_central"))
                   - (D(reference, "hydrated_" + f + "_central")
                      - D(reference, "matched_dry_" + f + "_central")))
                for f in FAMILIES}
    work = {}
    for metal in METALS:
        for family in FAMILIES:
            prefix = "hydrated_" + family + "_"
            for side in ("minus", "plus"):
                endpoint = prefix + side + "_" + metal
                center = prefix + "central_" + metal
                work[metal + "_" + family + "_" + side] = (
                    (model[endpoint]["energy"] - model[center]["energy"])
                    - (reference[endpoint]["energy"] - reference[center]["energy"]))
    paired_work = {f + "_" + s: work["Dy_" + f + "_" + s]
                   - work["La_" + f + "_" + s]
                   for f in FAMILIES for s in ("minus", "plus")}
    energies = {"balanced_exchange": exchange, "water_addition": addition,
                "individual_radial_work": work, "paired_radial_work": paired_work}
    vectors = {"energy": [v * math.sqrt(weight / len(energies[name]))
                           for name, weight in ENERGY_WEIGHTS.items()
                           for v in energies[name].values()]}
    per, controls, water_per, atomic = {}, {}, {}, {}
    for kind in ("La", "Dy", "paired"):
        vec, raw = [], []
        per[kind], water_per[kind] = {}, {}
        for q in DONOR + WATER:
            row = rows[q + "_La"]
            if kind == "paired":
                mf = body.subtract(model[q + "_Dy"]["forces"], model[q + "_La"]["forces"])
                rf = body.subtract(reference[q + "_Dy"]["forces"], reference[q + "_La"]["forces"])
            else:
                mf, rf = model[q + "_" + kind]["forces"], reference[q + "_" + kind]["forces"]
            weighted, components = body.rigid_residual(row, mf, rf)
            rms = 2 * math.sqrt(sum(v * v for v in weighted))
            maximum = max(abs(v) for v in components)
            record = {"rigid_RMS_kcal_mol_A": rms,
                      "maximum_body_component_error_kcal_mol_A": maximum,
                      "passed": rms <= 2 and maximum <= 5,
                      "raw_body_components_kcal_mol_A": components,
                      "model_body_response": body.response(row, mf),
                      "reference_body_response": body.response(row, rf)}
            if q in DONOR:
                per[kind][q] = record
                vec.extend(v / math.sqrt(len(DONOR)) for v in weighted)
                raw.extend(components)
            else:
                water_per[kind][q] = record
            residual = body.subtract(mf, rf)
            atomic[q + "_" + kind] = {
                "atomic_force_residual_kcal_mol_A": residual,
                "atom_vector_RMS_kcal_mol_A": math.sqrt(
                    sum(x * x for v in residual for x in v) / len(residual)),
                "used_in_loss": False}
        vectors[kind] = vec
        controls[kind] = {
            "category_balanced_rigid_RMS_kcal_mol_A": 2 * math.sqrt(sum(v * v for v in vec)),
            "max_body_component_error_kcal_mol_A": max(abs(v) for v in raw),
            "passed": all(v["passed"] for v in per[kind].values()),
            "every_fixture_required": True}
    energy_pass = all(abs(v) <= 1 for group in energies.values() for v in group.values())
    water_exchange = exchange_error("water_central")
    water_pass = (abs(water_exchange) <= 1 and all(
        r["passed"] for group in water_per.values() for r in group.values()))
    donor_pass = energy_pass and all(r["passed"] for r in controls.values())
    weighted = [math.sqrt(WEIGHTS[name]) * v for name in WEIGHTS for v in vectors[name]]
    return {
        "objective": sum(v * v for v in weighted), "weighted_residuals": weighted,
        "blocks": {name: sum(v * v for v in vectors[name]) for name in WEIGHTS},
        "block_weights": dict(WEIGHTS), "energy_subblock_weights": dict(ENERGY_WEIGHTS),
        "energy_errors_kcal_mol": energies, "energy_checks_passed": energy_pass,
        "controls": controls, "per_fixture": per, "atomic_diagnostics": atomic,
        "denominators": {"input_cells": 28, "donor_pairs_per_force_block": 12,
                         "balanced_exchange": 12, "water_addition": 3,
                         "individual_radial_work": 12, "paired_radial_work_controls_only": 6},
        "frozen_water_centers": {"unchanged": True, "per_fixture": water_per,
            "balanced_exchange_error_kcal_mol": water_exchange, "passed": water_pass},
        "donor_development_checks_passed": donor_pass,
        "development_checks_passed": donor_pass and water_pass,
        "development_scope": "Declared28 only; full eight-water qualification remains external",
        "full_water_qualification_assumed": False,
        "reference_eligibility": "Provisional development targets; central amide translation numerically supported",
        "physical_accuracy_qualified": False, "production_released": False}
