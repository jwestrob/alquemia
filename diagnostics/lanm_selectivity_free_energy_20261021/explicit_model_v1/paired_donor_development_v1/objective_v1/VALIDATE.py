"""One-thread, saved-only arithmetic/failure checks. No molecular imports/calls."""
import copy
import hashlib
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
D = ROOT / "diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1"
W = ROOT / "workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1"
TRAINING = W / "matched_dry_v1/training_v1/RESULT.json"
TRAINING_SHA = "99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb"
BODY = D / "matched_dry_v1/OBJECTIVE.py"
BODY_SHA = "bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819"
HERE = Path(__file__).resolve().parent


def pin(path):
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(output):
    started = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    checks = []

    def check(name, condition, details=None):
        checks.append({"name": name, "passed": bool(condition), "details": details})
        if not condition:
            raise AssertionError(name)

    def near(name, actual, expected, tolerance=1e-8):
        check(name, math.isfinite(actual) and abs(actual - expected) <= tolerance,
              {"actual": actual, "expected": expected, "tolerance": tolerance})

    def rejects(name, fn):
        try:
            fn()
        except (ValueError, KeyError):
            check(name, True)
        else:
            check(name, False)

    check("training_pin", pin(TRAINING)["sha256"] == TRAINING_SHA)
    check("body_helper_pin", pin(BODY)["sha256"] == BODY_SHA)
    obj, body = load("donor_objective", HERE / "OBJECTIVE.py"), load("rigid_body", BODY)
    data = json.loads(TRAINING.read_text())
    ref = {k: {"energy": r["reference"]["energy"], "forces": r["reference"]["forces"]}
           for k, r in data["rows"].items()}
    baseline = {k: {"energy": r["baseline"]["energy"], "forces": r["baseline"]["forces"]}
                for k, r in data["rows"].items()}

    def water(values):
        return {k: copy.deepcopy(values[k]) for k in obj.WATER_CELLS}

    def score(values, rows=data):
        return obj.evaluate(rows, values, body, water(values))

    saved = score(baseline)
    check("exact28_and_fixed_denominators", saved["denominators"] == {
        "input_cells": 28, "donor_pairs_per_force_block": 12, "balanced_exchange": 12,
        "water_addition": 3, "individual_radial_work": 12,
        "paired_radial_work_controls_only": 6})
    near("weighted_vector_closure", saved["objective"], sum(x*x for x in saved["weighted_residuals"]))
    near("four_block_closure", saved["objective"], .25*sum(saved["blocks"].values()))

    # Independent signed arithmetic from per-cell model-reference errors.
    errors = {k: baseline[k]["energy"] - ref[k]["energy"] for k in ref}
    paired_errors = {q: errors[q+"_Dy"] - errors[q+"_La"] for q in obj.WATER + obj.DONOR}
    ee = saved["energy_errors_kcal_mol"]
    for q in obj.DONOR:
        near("saved_exchange_"+q, ee["balanced_exchange"][q],
             paired_errors[q]-paired_errors["aquo8_central"])
    for f in obj.FAMILIES:
        near("saved_addition_"+f, ee["water_addition"][f],
             paired_errors["hydrated_"+f+"_central"]-paired_errors["matched_dry_"+f+"_central"])
        for m in obj.METALS:
            for side in ("minus", "plus"):
                key = m+"_"+f+"_"+side
                near("saved_radial_"+key, ee["individual_radial_work"][key],
                     errors["hydrated_"+f+"_"+side+"_"+m]-errors["hydrated_"+f+"_central_"+m])
    expected_energy = (.5*sum(x*x for x in ee["balanced_exchange"].values())/12
                       +.3*sum(x*x for x in ee["water_addition"].values())/3
                       +.2*sum(x*x for x in ee["individual_radial_work"].values())/12)
    near("saved_energy_weights", saved["blocks"]["energy"], expected_energy)

    zero = score(ref)
    near("reference_self_zero", zero["objective"], 0.)
    check("reference_self_controls", zero["development_checks_passed"])
    shifted = copy.deepcopy(baseline)
    for k in shifted:
        shifted[k]["energy"] += 17. if k.endswith("_La") else -11.
    same = score(shifted)
    near("independent_metal_energy_zero_cancellation", same["objective"], saved["objective"])
    check("independent_metal_offsets_preserve_controls", same["controls"] == saved["controls"])
    shifted = copy.deepcopy(ref)
    for q in obj.DONOR:
        offset = 5. if q.startswith("hydrated_") else 2.
        for m in obj.METALS:
            shifted[q+"_"+m]["energy"] += offset
    near("common_wet_dry_composition_offsets_cancel", score(shifted)["objective"], 0.)

    for metal, sign in (("Dy", 1.), ("La", -1.)):
        trial = copy.deepcopy(ref)
        trial["hydrated_acetamide_central_"+metal]["energy"] += 1.
        s = score(trial); e = s["energy_errors_kcal_mol"]
        near(metal+"_central_exchange_sign", e["balanced_exchange"]["hydrated_acetamide_central"], sign)
        near(metal+"_central_wet_dry_sign", e["water_addition"]["acetamide"], sign)
        for side in ("minus", "plus"):
            near(metal+"_individual_work_"+side, e["individual_radial_work"][metal+"_acetamide_"+side], -1.)
            near(metal+"_paired_work_"+side, e["paired_radial_work"]["acetamide_"+side], -sign)
        near(metal+"_central_energy_weight", s["objective"], .04375)
    trial = copy.deepcopy(ref)
    trial["aquo8_central_Dy"]["energy"] += 1.
    s = score(trial)
    check("aquo_anchor_all12_signs", all(abs(x+1.)<1e-8 for x in s["energy_errors_kcal_mol"]["balanced_exchange"].values()))
    near("aquo_anchor_weight", s["objective"], .125)

    # One fixture may fail although the 12-fixture average passes. Shared force
    # errors must cancel only from the paired block, never individual blocks.
    trial = copy.deepcopy(ref)
    for metal in obj.METALS:
        f = trial["hydrated_acetamide_minus_"+metal]["forces"]
        f[0] = [x+4. for x in f[0]]
    s = score(trial)
    for metal in obj.METALS:
        near(metal+"_known_ion_category_weight", s["per_fixture"][metal]["hydrated_acetamide_minus"]["rigid_RMS_kcal_mol_A"], math.sqrt(16./3.))
        near(metal+"_fixture_mean_weight", s["blocks"][metal], 1./9.)
        check(metal+"_one_fixture_not_hidden", not s["controls"][metal]["passed"]
              and s["controls"][metal]["category_balanced_rigid_RMS_kcal_mol_A"] < 2.)
    near("shared_force_error_cancels_paired_only", s["blocks"]["paired"], 0.)
    near("symmetric_force_block_weights", s["objective"], 1./18.)
    check("individual_failure_blocks_development", not s["development_checks_passed"])

    missing = copy.deepcopy(ref); missing.pop("hydrated_acetamide_plus_Dy")
    rejects("missing_model_cell_rejected", lambda: score(missing))
    extra = copy.deepcopy(ref); extra["extra"] = copy.deepcopy(ref["water_central_La"])
    rejects("extra_model_cell_rejected", lambda: score(extra))
    bad_data = copy.deepcopy(data); bad_data["rows"].pop("matched_dry_acetate_bi_central_La")
    rejects("missing_data_cell_rejected", lambda: score(ref, bad_data))
    bad = copy.deepcopy(ref); bad["hydrated_acetamide_plus_Dy"]["forces"][0][0] = float("nan")
    rejects("nan_force_rejected", lambda: score(bad))
    bad = copy.deepcopy(ref); bad["hydrated_acetamide_plus_Dy"]["forces"].pop()
    rejects("short_force_array_rejected", lambda: score(bad))
    bad_data = copy.deepcopy(data)
    bad_data["rows"]["hydrated_acetamide_plus_Dy"]["geometry"]["positions_A"][1][0] += .001
    rejects("mismatched_paired_coordinates_rejected", lambda: score(ref, bad_data))
    bad_data = copy.deepcopy(data)
    bad_data["rows"]["hydrated_acetamide_plus_Dy"]["status"] = "failed"
    rejects("incomplete_reference_row_rejected", lambda: score(ref, bad_data))
    rejects("missing_frozen_water_baseline_rejected", lambda: obj.evaluate(data, ref, body, {}))
    altered = copy.deepcopy(ref); altered["water_central_La"]["energy"] += .001
    rejects("donor_cannot_change_water_anchor", lambda: obj.evaluate(data, altered, body, water(ref)))

    result = {"status": "saved_objective_validation_passed", "checks": checks,
              "check_count": len(checks), "source": pin(HERE/"OBJECTIVE.py"),
              "validation_source": pin(Path(__file__)), "body_helper": pin(BODY),
              "training": pin(TRAINING), "original_design": pin(HERE.parent/"DESIGN.md"),
              "baseline_role": "Historical original28 only; no future water candidate bound",
              "saved_original_baseline_score": saved,
              "future_water_candidate_bound": False, "new_native_program_invocations": 0,
              "new_molecular_calls": 0, "local_threads": 1,
              "reserved_validation_accessed": False, "fit_released": False,
              "elapsed_seconds": time.monotonic()-started}
    (output/"RESULT.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"status": result["status"], "checks": len(checks),
                      "baseline_objective": saved["objective"], "result": pin(output/"RESULT.json"),
                      "elapsed_seconds": result["elapsed_seconds"], "new_native_program_invocations": 0}))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
