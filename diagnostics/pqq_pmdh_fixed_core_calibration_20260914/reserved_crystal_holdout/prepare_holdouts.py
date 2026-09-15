#!/usr/bin/env python3
"""Prepare the preregistered PQQ-MDH crystal holdouts after a passing gate.

This command is intentionally dormant until it is given the immutable JSON
result from a successful 25-member v3 calibration.  It prepares La/Ca vertical
pairs but never starts ORCA.  The atomic primary release is 1H4I plus 4MAE;
6OC6 is included only with an explicit secondary-target flag.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import datetime as dt
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Mapping, Sequence

import gemmi


HERE = Path(__file__).resolve().parent
CALIBRATION_DIR = HERE.parent
PROJECT = CALIBRATION_DIR.parents[1]
SCRIPTS = PROJECT / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(CALIBRATION_DIR))

import carve_with_pqq as base  # noqa: E402
import fixed_core_carver as fixed  # noqa: E402
import protonate_standard_only as standard_only  # noqa: E402
from coordination_policy import (  # noqa: E402
    ALTLOC_POLICY_ID,
    COORDINATION_CUTOFF_A,
    MAX_DIRECT_N,
    WATER_NAMES,
    donor_type,
    selected_residue_atoms,
)
from pqq_microstates import DEFAULT_PQQ_MICROSTATE, prepare_pqq  # noqa: E402


PROTOCOL_ID = "pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3"
HOLDOUT_PREPARATION_SCHEMA = "alchemical_bvs.pqq_fixed_core_holdout_preparation.v1"
HOLDOUT_NORMALIZATION_SCHEMA = "alchemical_bvs.pqq_fixed_core_holdout_normalization.v1"
HOLDOUT_CARVE_SCHEMA = "alchemical_bvs.carve_manifest.v1"
HOLDOUT_PINS_SCHEMA = "alchemical_bvs.pqq_fixed_core_holdout_pins.v1"
HEAVY_CHECK_SCHEMA = "alchemical_bvs.pqq_fixed_core_holdout_heavy_coordinate_check.v1"
CALIBRATION_RESULT_SCHEMA = "alchemical_bvs.pqq_fixed_core_calibration_result.v1"
CORE_POLICY_ID = "pqq_pmdh_E_N_D_Dplus2acidic_Dnetwork_cation_v1"
COORDINATION_POLICY_ID = "typed_cn6_3p1A_maxN2_crystal_transfer_v1"
DRY_POLICY_ID = "dry_fixed_core_v3"
METHOD_ID = "orca_r2scan3c_cpcm_water_native_basis_v2"
EXECUTION_POLICY_ID = "atomic_primary_pair_optional_secondary_full_node_manifested_orca_v1"
EXPECTED_ORCA_DIRECTIVE = "! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3"
COORDINATE_TOLERANCE_A = 0.001
AUDIT_DISTANCE_TOLERANCE_A = 0.011
STANDARD_AA = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
)
EXPECTED_HOLDOUTS = {
    "1H4I": ("primary", "Ca/MxaF"),
    "4MAE": ("primary", "Ln/XoxF"),
    "6OC6": ("secondary", "Ln/XoxF"),
}
ROLE_FIELDS = {
    "anchor_glutamate": "glu_selector",
    "anchor_asparagine": "asn_selector",
    "catalytic_aspartate": "catalytic_asp_selector",
    "extra_acidic_ligand_homolog": "plus2_selector",
    "catalytic_asp_cationic_partner": "cationic_partner_selector",
}
ROLE_DISTANCE_FIELDS = {
    "anchor_glutamate": "glu_metal_A",
    "anchor_asparagine": "asn_metal_A",
    "catalytic_aspartate": "catalytic_asp_metal_A",
    "extra_acidic_ligand_homolog": "plus2_metal_A",
}


class HoldoutPreparationError(RuntimeError):
    """A frozen holdout selector, chemistry rule, or provenance check failed."""


@dataclass(frozen=True)
class ResidueSelector:
    chain: str
    resname: str
    resnum: int
    icode: str = ""

    def label(self) -> str:
        return f"{self.chain}:{self.resname}{self.resnum}{self.icode}"


@dataclass(frozen=True)
class AtomSelector:
    residue: ResidueSelector
    atom: str

    def label(self) -> str:
        return f"{self.residue.label()}/{self.atom}"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise HoldoutPreparationError(f"cannot read JSON object {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise HoldoutPreparationError(f"JSON root is not an object: {path}")
    return value


def write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def write_text_new(path: Path, text: str) -> None:
    with path.open("x") as handle:
        handle.write(text)


def file_record(path: Path) -> dict[str, str]:
    path = path.resolve()
    if not path.is_file():
        raise HoldoutPreparationError(f"required file is absent: {path}")
    return {"path": str(path), "sha256": sha256_file(path)}


def verify_file_record(record: Any, label: str) -> Path:
    if not isinstance(record, dict):
        raise HoldoutPreparationError(f"{label} is not a file record")
    raw_path = record.get("path")
    expected = record.get("sha256")
    if not isinstance(raw_path, str) or not isinstance(expected, str):
        raise HoldoutPreparationError(f"{label} lacks path/SHA-256")
    path = Path(raw_path).resolve()
    if not path.is_file():
        raise HoldoutPreparationError(f"{label} is missing: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise HoldoutPreparationError(
            f"{label} SHA-256 mismatch: expected {expected}, observed {observed}"
        )
    return path


def parse_residue_selector(value: str) -> ResidueSelector:
    match = re.fullmatch(r"([^:]+):(.+?)(-?[0-9]+)([A-Za-z]?)", value)
    if match is None:
        raise HoldoutPreparationError(f"invalid residue selector {value!r}")
    return ResidueSelector(
        chain=match.group(1),
        resname=match.group(2).upper(),
        resnum=int(match.group(3)),
        icode=match.group(4).upper(),
    )


def parse_atom_selector(value: str) -> AtomSelector:
    if value.count("/") != 1:
        raise HoldoutPreparationError(f"invalid atom selector {value!r}")
    residue, atom = value.rsplit("/", 1)
    if not re.fullmatch(r"[A-Za-z0-9_]+", atom):
        raise HoldoutPreparationError(f"invalid atom name in selector {value!r}")
    return AtomSelector(parse_residue_selector(residue), atom.upper())


def parse_bool(value: str, label: str) -> bool:
    if value not in {"true", "false"}:
        raise HoldoutPreparationError(f"{label} must be true or false, not {value!r}")
    return value == "true"


def parse_float(value: str, label: str, *, allow_na: bool = False) -> float | None:
    if value == "NA" and allow_na:
        return None
    try:
        result = float(value)
    except ValueError as exc:
        raise HoldoutPreparationError(f"invalid {label}: {value!r}") from exc
    if not math.isfinite(result):
        raise HoldoutPreparationError(f"nonfinite {label}: {value!r}")
    return result


def insertion_code(residue: gemmi.Residue) -> str:
    value = str(residue.seqid.icode).strip()
    return "" if value in {"", "\x00", ".", "?"} else value.upper()


def atom_altloc(atom: gemmi.Atom) -> str:
    value = str(atom.altloc).strip()
    return "" if value in {"", "\x00", ".", "?"} else value.upper()


def residue_selector(chain: gemmi.Chain, residue: gemmi.Residue) -> ResidueSelector:
    return ResidueSelector(
        chain.name,
        residue.name.upper().strip(),
        residue.seqid.num,
        insertion_code(residue),
    )


def verify_runtime(calibration_pins: Mapping[str, Any]) -> None:
    runtime = calibration_pins.get("preparation_runtime")
    if not isinstance(runtime, dict):
        raise HoldoutPreparationError("calibration pins lack preparation runtime")
    python = runtime.get("python")
    if not isinstance(python, dict):
        raise HoldoutPreparationError("calibration pins lack Python runtime")
    expected_python = Path(str(python.get("path", ""))).resolve()
    if Path(sys.executable).resolve() != expected_python:
        raise HoldoutPreparationError(
            f"run with pinned interpreter {expected_python}, not {Path(sys.executable).resolve()}"
        )
    version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    if version != python.get("version"):
        raise HoldoutPreparationError("Python version differs from frozen runtime")
    observed_packages = {
        "gemmi": getattr(gemmi, "__version__", None),
        "openmm": importlib.metadata.version("openmm"),
        "pdbfixer": importlib.metadata.version("pdbfixer"),
    }
    if observed_packages != runtime.get("packages"):
        raise HoldoutPreparationError(
            f"package versions {observed_packages} differ from frozen runtime"
        )


def verify_holdout_pins(pins_path: Path) -> tuple[dict[str, Any], dict[str, Any], Path]:
    pins_path = pins_path.resolve()
    pins = read_object(pins_path)
    if pins.get("schema_version") != HOLDOUT_PINS_SCHEMA:
        raise HoldoutPreparationError("unsupported holdout-pins schema")
    if pins.get("protocol_id") != PROTOCOL_ID:
        raise HoldoutPreparationError("holdout pins carry another protocol ID")
    calibration_pins_path = verify_file_record(
        pins.get("calibration_implementation_pins"), "calibration implementation pins"
    )
    calibration_pins = fixed.verify_pins(calibration_pins_path)
    verify_file_record(pins.get("calibration_result"), "locked calibration result")
    spec_path = verify_file_record(pins.get("holdout_spec"), "holdout selector authority")
    verify_file_record(pins.get("holdout_selector_audit"), "holdout selector audit")
    if pins.get("holdout_spec") != calibration_pins.get("reserved_holdout_spec"):
        raise HoldoutPreparationError("holdout spec differs from calibration's frozen record")
    if pins.get("holdout_selector_audit") != calibration_pins.get("reserved_holdout_audit"):
        raise HoldoutPreparationError("holdout audit differs from calibration's frozen record")

    implementation = pins.get("holdout_implementation")
    if not isinstance(implementation, dict):
        raise HoldoutPreparationError("holdout pins lack implementation records")
    required_implementation = {
        "prepare_holdouts",
        "run_holdouts",
        "score_holdouts",
        "run_sbatch",
        "design_audit",
    }
    if set(implementation) != required_implementation:
        raise HoldoutPreparationError("holdout implementation file set changed")
    verified_implementation = {
        name: verify_file_record(record, f"holdout implementation {name}")
        for name, record in implementation.items()
    }
    own_path = verified_implementation["prepare_holdouts"]
    if own_path != Path(__file__).resolve():
        raise HoldoutPreparationError("executing holdout preparer is not its pinned file")

    inherited = pins.get("inherited_helpers")
    expected_inherited = {
        "fixed_core_carver": calibration_pins["experiment_helpers"]["fixed_core_carver"],
        "protonate_standard_only": calibration_pins["experiment_helpers"]["protonate_standard_only"],
        "carve_with_pqq": calibration_pins["canonical_helpers"]["carve_with_pqq"],
        "pqq_microstates": calibration_pins["canonical_helpers"]["pqq_microstates"],
        "coordination_policy": calibration_pins["canonical_helpers"]["coordination_policy"],
        "run_orca_task_manifest": calibration_pins["canonical_helpers"]["run_orca_task_manifest"],
        "render_orca_runtime_input": calibration_pins["canonical_helpers"]["render_orca_runtime_input"],
        "calibration_score": calibration_pins["experiment_helpers"]["score"],
    }
    if inherited != expected_inherited:
        raise HoldoutPreparationError("inherited helper records differ from calibration pins")
    for name, record in expected_inherited.items():
        verify_file_record(record, f"inherited helper {name}")
    if Path(fixed.__file__).resolve() != Path(expected_inherited["fixed_core_carver"]["path"]).resolve():
        raise HoldoutPreparationError("imported fixed-core carver differs from frozen helper")
    if Path(standard_only.__file__).resolve() != Path(
        expected_inherited["protonate_standard_only"]["path"]
    ).resolve():
        raise HoldoutPreparationError("imported protonation wrapper differs from frozen helper")
    if pins.get("protonation_subprotocol_id") != calibration_pins.get(
        "protonation_subprotocol", {}
    ).get("id"):
        raise HoldoutPreparationError("holdout protonation subprotocol is not the calibration one")
    if pins.get("selection_policy") != {
        "model": "1",
        "chain": "A",
        "alternate_conformer_policy": ALTLOC_POLICY_ID,
        "metal_preprotonation_policy": "residue_name_to_LA_preserve_native_atom_name_and_element",
        "water_policy": "exclude_all_source_and_synthetic_waters",
        "noncore_heterogen_policy": "exclude_all_including_4MAE_A:15P603",
        "primary_pair_atomic": True,
        "primary_pair_release_runnable": True,
        "secondary_6OC6_requires_explicit_flag": True,
    }:
        raise HoldoutPreparationError("holdout selection policy pins changed")
    if pins.get("frozen_transfer_rule") != {
        "band_scale": "S_aquo_gauge_kcal_mol",
        "band_source": "locked_calibration_result.released_supported_bands",
        "1H4I": "S<=U_max_Ca_S_kcal_mol",
        "4MAE": "S>=L_min_La_S_kcal_mol",
        "open_gap": "indeterminate_and_primary_failure",
        "wrong_band_or_invalid": "primary_failure",
        "6OC6": "secondary_only_never_affects_primary_verdict",
        "threshold_refit_allowed": False,
    }:
        raise HoldoutPreparationError("frozen holdout transfer rule pins changed")
    verify_runtime(calibration_pins)
    return pins, calibration_pins, spec_path


def load_holdout_rows(spec_path: Path) -> list[dict[str, str]]:
    with spec_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_id = {row.get("pdb_id", ""): row for row in rows}
    if len(rows) != 3 or set(by_id) != set(EXPECTED_HOLDOUTS):
        raise HoldoutPreparationError("holdout spec is not the frozen three-structure set")
    for pdb_id, (role, biological_class) in EXPECTED_HOLDOUTS.items():
        row = by_id[pdb_id]
        if row.get("holdout_role") != role or row.get("biological_class") != biological_class:
            raise HoldoutPreparationError(f"frozen disposition changed for {pdb_id}")
        if row.get("model") != "1" or row.get("selected_chain") != "A":
            raise HoldoutPreparationError(f"{pdb_id} is not frozen to model 1, chain A")
        if row.get("pqq_schema") != "pdb_ccd_pqq_v1":
            raise HoldoutPreparationError(f"{pdb_id} does not use frozen CCD PQQ")
        if not row.get("dry_policy", "").startswith(f"{DRY_POLICY_ID}:"):
            raise HoldoutPreparationError(f"{pdb_id} lacks the frozen dry policy")
    return [by_id["1H4I"], by_id["4MAE"], by_id["6OC6"]]


def verify_calibration_gate(
    result_path: Path,
    calibration_pins_path: Path,
    calibration_pins: Mapping[str, Any],
) -> dict[str, Any]:
    result_path = result_path.resolve()
    result = read_object(result_path)
    if result.get("schema_version") != CALIBRATION_RESULT_SCHEMA:
        raise HoldoutPreparationError("gate record is not a v3 calibration result")
    if result.get("protocol_id") != PROTOCOL_ID:
        raise HoldoutPreparationError("gate record carries another protocol")
    expected_pins = {
        "path": str(calibration_pins_path.resolve()),
        "sha256": sha256_file(calibration_pins_path),
    }
    if result.get("implementation_pins") != expected_pins:
        raise HoldoutPreparationError("gate result is not bound to the frozen calibration pins")
    preregistration = calibration_pins.get("preregistration")
    if not isinstance(preregistration, dict) or result.get("preregistration") != {
        "path": preregistration.get("path"),
        "sha256": preregistration.get("sha256"),
    }:
        raise HoldoutPreparationError("gate result is not bound to the frozen preregistration")
    preparation_path = verify_file_record(result.get("preparation"), "calibration preparation")
    preparation = read_object(preparation_path)
    if (
        preparation.get("schema_version") != "alchemical_bvs.pqq_fixed_core_preparation.v1"
        or preparation.get("protocol_id") != PROTOCOL_ID
        or preparation.get("status") != "ready_for_orca"
        or preparation.get("target_count") != 25
        or preparation.get("task_count") != 50
        or preparation.get("implementation_pins") != expected_pins
    ):
        raise HoldoutPreparationError("gate result points to an invalid calibration preparation")
    calibration = result.get("calibration")
    if not isinstance(calibration, dict):
        raise HoldoutPreparationError("gate result lacks calibration data")
    gates = calibration.get("gates")
    required_gates = {
        "all_25_valid": True,
        "strict_separation": True,
        "gap_at_least_5_kcal_mol": True,
        "strict_LOO_midpoint_25_of_25": True,
    }
    scores = result.get("scores")
    if not isinstance(scores, list) or len(scores) != 25:
        raise HoldoutPreparationError("gate result does not contain 25 calibration scores")
    core_map_path = verify_file_record(calibration_pins.get("core_map"), "calibration core map")
    with core_map_path.open(newline="") as handle:
        core_rows = list(csv.DictReader(handle, delimiter="\t"))
    frozen_classes = {row["id"]: row["class"] for row in core_rows}
    if len(frozen_classes) != 25:
        raise HoldoutPreparationError("calibration core-map identity set changed")
    score_by_id: dict[str, tuple[str, float]] = {}
    for score in scores:
        if not isinstance(score, dict):
            raise HoldoutPreparationError("malformed score in calibration gate record")
        panel_id = score.get("panel_id")
        biological_class = score.get("class")
        value = score.get("R_kcal_mol")
        if (
            not isinstance(panel_id, str)
            or panel_id in score_by_id
            or biological_class not in {"La", "Ca"}
            or isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
        ):
            raise HoldoutPreparationError("invalid or duplicate calibration score")
        if frozen_classes.get(panel_id) != biological_class:
            raise HoldoutPreparationError(f"calibration class/identity changed for {panel_id}")
        score_by_id[panel_id] = (biological_class, float(value))
    if set(score_by_id) != set(frozen_classes):
        raise HoldoutPreparationError("calibration result does not cover the frozen panel")
    ca_values = [value for biological_class, value in score_by_id.values() if biological_class == "Ca"]
    la_values = [value for biological_class, value in score_by_id.values() if biological_class == "La"]
    if len(ca_values) != 14 or len(la_values) != 11:
        raise HoldoutPreparationError("calibration score class counts changed")
    recomputed_upper = max(ca_values)
    recomputed_lower = min(la_values)
    recomputed_gap = recomputed_lower - recomputed_upper
    recomputed_threshold = (recomputed_upper + recomputed_lower) / 2.0
    auc_wins = sum(
        1.0 if la > ca else 0.5 if la == ca else 0.0
        for la in la_values
        for ca in ca_values
    )
    recomputed_auc = auc_wins / (len(la_values) * len(ca_values))
    loo_correct = 0
    loo_training_strict = True
    for held_id, (held_class, held_value) in score_by_id.items():
        training = {
            panel_id: item for panel_id, item in score_by_id.items() if panel_id != held_id
        }
        training_ca = [value for biological_class, value in training.values() if biological_class == "Ca"]
        training_la = [value for biological_class, value in training.values() if biological_class == "La"]
        training_upper = max(training_ca)
        training_lower = min(training_la)
        strict = training_lower > training_upper
        loo_training_strict = loo_training_strict and strict
        if strict:
            midpoint = (training_upper + training_lower) / 2.0
            predicted = "La" if held_value > midpoint else "Ca" if held_value < midpoint else "failure"
            loo_correct += predicted == held_class
    recomputed_gates = {
        "all_25_valid": len(score_by_id) == 25,
        "strict_separation": recomputed_lower > recomputed_upper,
        "gap_at_least_5_kcal_mol": recomputed_gap >= 5.0,
        "strict_LOO_midpoint_25_of_25": loo_training_strict and loo_correct == 25,
    }
    if (
        calibration.get("verdict") != "CALIBRATABLE"
        or calibration.get("valid_pair_count") != 25
        or calibration.get("strict_separation") is not True
        or calibration.get("gap_floor_pass") is not True
        or calibration.get("loo_correct_count") != 25
        or gates != required_gates
        or calibration.get("threshold_R_kcal_mol") is None
        or calibration.get("released_supported_bands") is None
        or recomputed_gates != required_gates
    ):
        raise HoldoutPreparationError("all preregistered calibration gates have not passed")
    try:
        upper = float(calibration["U_max_Ca_kcal_mol"])
        lower = float(calibration["L_min_La_kcal_mol"])
        gap = float(calibration["gap_kcal_mol"])
        threshold = float(calibration["threshold_R_kcal_mol"])
        reported_auc = float(calibration["auroc"])
    except (KeyError, TypeError, ValueError) as exc:
        raise HoldoutPreparationError("released calibration boundaries are malformed") from exc
    if not all(math.isfinite(value) for value in (upper, lower, gap, threshold, reported_auc)):
        raise HoldoutPreparationError("released calibration boundaries are nonfinite")
    if not lower > upper or not gap >= 5.0:
        raise HoldoutPreparationError("released calibration boundaries are internally invalid")
    if (
        abs(upper - recomputed_upper) > 1.0e-9
        or abs(lower - recomputed_lower) > 1.0e-9
        or abs(gap - recomputed_gap) > 1.0e-9
        or abs(threshold - recomputed_threshold) > 1.0e-9
        or abs(reported_auc - recomputed_auc) > 1.0e-12
        or abs(gap - (lower - upper)) > 1.0e-9
        or abs(threshold - (upper + lower) / 2.0) > 1.0e-9
    ):
        raise HoldoutPreparationError("released calibration boundary arithmetic is inconsistent")
    bands = calibration["released_supported_bands"]
    if (
        not isinstance(bands, dict)
        or bands.get("Ca", {}).get("R_kcal_mol_max") != upper
        or bands.get("La", {}).get("R_kcal_mol_min") != lower
        or bands.get("indeterminate", {}).get("R_kcal_mol_open_interval") != [upper, lower]
    ):
        raise HoldoutPreparationError("released R bands do not match recomputed U/L")
    try:
        aquo_delta = float(calibration_pins["aquo_reference"]["delta_E_aquo_hartree"])
        upper_s = float(calibration["U_max_Ca_S_kcal_mol"])
        lower_s = float(calibration["L_min_La_S_kcal_mol"])
        threshold_s = float(calibration["threshold_S_kcal_mol"])
    except (KeyError, TypeError, ValueError) as exc:
        raise HoldoutPreparationError("released S bands or aquo gauge are malformed") from exc
    aquo_a = aquo_delta * 627.509474
    if (
        not all(math.isfinite(value) for value in (aquo_delta, aquo_a, upper_s, lower_s, threshold_s))
        or abs(upper_s - (upper - aquo_a)) > 1.0e-9
        or abs(lower_s - (lower - aquo_a)) > 1.0e-9
        or abs(threshold_s - (threshold - aquo_a)) > 1.0e-9
        or bands.get("Ca", {}).get("S_kcal_mol_max") != upper_s
        or bands.get("La", {}).get("S_kcal_mol_min") != lower_s
        or bands.get("indeterminate", {}).get("S_kcal_mol_open_interval")
        != [upper_s, lower_s]
    ):
        raise HoldoutPreparationError("released S bands disagree with frozen R/aquo values")
    interpretation = result.get("interpretation")
    if not isinstance(interpretation, dict) or interpretation.get(
        "reserved_crystal_pair_consumed"
    ) is not False:
        raise HoldoutPreparationError("gate record does not preserve the reserved holdouts")
    return {
        "result": file_record(result_path),
        "verdict": "CALIBRATABLE",
        "valid_pair_count": 25,
        "U_max_Ca_kcal_mol": upper,
        "L_min_La_kcal_mol": lower,
        "gap_kcal_mol": gap,
        "threshold_R_kcal_mol": threshold,
        "U_max_Ca_S_kcal_mol": upper_s,
        "L_min_La_S_kcal_mol": lower_s,
        "threshold_S_kcal_mol": threshold_s,
        "aquo_A_kcal_mol": aquo_a,
        "gates": required_gates,
        "individual_calibration_energies_consumed_by_holdout_preparation": False,
    }


def exact_model(structure: gemmi.Structure, model_name: str) -> gemmi.Model:
    matches = [model for model in structure if str(model.num) == model_name]
    if len(matches) != 1:
        observed = [str(model.num) for model in structure]
        raise HoldoutPreparationError(
            f"model selector {model_name!r} matched {len(matches)} models; observed {observed}"
        )
    return matches[0]


def exact_chain(model: gemmi.Model, chain_name: str) -> gemmi.Chain:
    matches = [chain for chain in model if chain.name == chain_name]
    if len(matches) != 1:
        raise HoldoutPreparationError(
            f"chain selector {chain_name!r} matched {len(matches)} chains"
        )
    return matches[0]


def resolve_residue(model: gemmi.Model, selector: ResidueSelector) -> gemmi.Residue:
    matches = [
        residue
        for chain in model
        for residue in chain
        if residue_selector(chain, residue) == selector
    ]
    if len(matches) != 1:
        raise HoldoutPreparationError(
            f"residue selector {selector.label()} matched {len(matches)} residues"
        )
    return matches[0]


def resolve_atom(model: gemmi.Model, selector: AtomSelector) -> gemmi.Atom:
    residue = resolve_residue(model, selector.residue)
    matches = [
        atom
        for atom in selected_residue_atoms(residue)
        if atom.name.upper().strip() == selector.atom
    ]
    if len(matches) != 1:
        raise HoldoutPreparationError(
            f"atom selector {selector.label()} matched {len(matches)} selected atoms"
        )
    return matches[0]


def residue_has_altloc(residue: gemmi.Residue) -> bool:
    return any(atom_altloc(atom) for atom in residue)


def chosen_altloc(residue: gemmi.Residue) -> str:
    labels = sorted({atom_altloc(atom) for atom in residue} - {""})
    scores = {
        label: sum(
            float(atom.occ)
            if math.isfinite(float(atom.occ)) and float(atom.occ) >= 0.0
            else 0.0
            for atom in residue
            if atom_altloc(atom) == label
        )
        for label in labels
    }
    return min(labels, key=lambda label: (-scores[label], label != "A", label)) if labels else ""


def altloc_decision(residue: gemmi.Residue) -> dict[str, Any]:
    labels = sorted({atom_altloc(atom) for atom in residue} - {""})
    choice = chosen_altloc(residue)
    selected = selected_residue_atoms(residue)
    selected_signatures = {
        (
            atom.name.upper().strip(),
            atom.element.name.upper(),
            atom_altloc(atom),
            round(atom.pos.x, 6),
            round(atom.pos.y, 6),
            round(atom.pos.z, 6),
        )
        for atom in selected
    }
    discarded = []
    for atom in residue:
        signature = (
            atom.name.upper().strip(),
            atom.element.name.upper(),
            atom_altloc(atom),
            round(atom.pos.x, 6),
            round(atom.pos.y, 6),
            round(atom.pos.z, 6),
        )
        if signature in selected_signatures:
            continue
        discarded.append(
            {
                "atom": atom.name.upper().strip(),
                "element": atom.element.name.title(),
                "altloc": atom_altloc(atom),
                "occupancy": float(atom.occ),
                "xyz_A": [round(atom.pos.x, 6), round(atom.pos.y, 6), round(atom.pos.z, 6)],
            }
        )
    return {
        "policy_id": ALTLOC_POLICY_ID,
        "available_nonblank_labels": labels,
        "chosen_nonblank_label": choice,
        "output_altloc_blank": True,
        "discarded_atoms": discarded,
    }


def assert_close(observed: float, expected: float, label: str, tolerance: float) -> None:
    if abs(observed - expected) > tolerance:
        raise HoldoutPreparationError(
            f"{label} mismatch: observed {observed:.6f}, expected {expected:.6f}"
        )


def expected_donor_ledger(row: Mapping[str, str]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    raw = row.get("typed_donor_ledger_3p1_A", "")
    for item in raw.split(";") if raw else []:
        if item.count("=") != 1:
            raise HoldoutPreparationError(f"malformed frozen donor item {item!r}")
        selector, distance = item.split("=", 1)
        records.append(
            {
                "selector": parse_atom_selector(selector),
                "distance_A": float(distance),
            }
        )
    return records


def typed_source_contacts(
    model: gemmi.Model,
    metal_selector: AtomSelector,
    metal: gemmi.Atom,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for chain in model:
        for residue in chain:
            key = residue_selector(chain, residue)
            prepared = None
            for atom in selected_residue_atoms(residue):
                atom_name = atom.name.upper().strip()
                element = atom.element.name.upper()
                if key == metal_selector.residue and atom_name == metal_selector.atom:
                    continue
                if element in {"H", "D"}:
                    continue
                distance = metal.pos.dist(atom.pos)
                if distance > COORDINATION_CUTOFF_A:
                    continue
                kind = None
                if key.resname == "PQQ":
                    if prepared is None:
                        prepared = prepare_pqq(residue, DEFAULT_PQQ_MICROSTATE)
                    if atom_name in prepared.direct_donors:
                        kind = "cofactor_O" if element == "O" else "cofactor_N"
                else:
                    kind = donor_type(key.resname, atom_name, element)
                if kind is not None:
                    records.append(
                        {
                            "selector": AtomSelector(key, atom_name),
                            "element": "H" if element == "D" else element.title(),
                            "donor_type": kind,
                            "distance_A": distance,
                        }
                    )
    records.sort(
        key=lambda record: (
            record["distance_A"],
            record["selector"].residue.chain,
            record["selector"].residue.resnum,
            record["selector"].atom,
        )
    )
    return records


def validate_donor_ledger(
    row: Mapping[str, str], observed: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    expected = expected_donor_ledger(row)
    expected_cn = int(row["fixed_typed_CN_3p1"])
    if len(expected) != expected_cn or len(observed) != expected_cn:
        raise HoldoutPreparationError(
            f"{row['pdb_id']} typed-CN mismatch: expected {expected_cn}, observed {len(observed)}"
        )
    expected_by_selector = {item["selector"]: item for item in expected}
    observed_by_selector = {item["selector"]: item for item in observed}
    if len(expected_by_selector) != len(expected) or len(observed_by_selector) != len(observed):
        raise HoldoutPreparationError(f"{row['pdb_id']} donor ledger contains duplicates")
    if set(expected_by_selector) != set(observed_by_selector):
        missing = sorted(item.label() for item in set(expected_by_selector) - set(observed_by_selector))
        added = sorted(item.label() for item in set(observed_by_selector) - set(expected_by_selector))
        raise HoldoutPreparationError(
            f"{row['pdb_id']} donor identities changed: missing={missing}, added={added}"
        )
    result: list[dict[str, Any]] = []
    # Preserve the deliberately frozen TSV order in emitted provenance.  It is
    # residue-grouped rather than globally sorted by distance.
    for want in expected:
        got = observed_by_selector[want["selector"]]
        assert_close(
            float(got["distance_A"]),
            float(want["distance_A"]),
            f"{row['pdb_id']} donor {want['selector'].label()}",
            AUDIT_DISTANCE_TOLERANCE_A,
        )
        result.append(
            {
                "selector": got["selector"].label(),
                "element": got["element"],
                "donor_type": got["donor_type"],
                "distance_A": round(float(got["distance_A"]), 3),
            }
        )
    if sum(record["element"].upper() == "N" for record in result) > MAX_DIRECT_N:
        raise HoldoutPreparationError(f"{row['pdb_id']} exceeds the two-N donor cap")
    return result


def validate_role_selectors(
    model: gemmi.Model,
    row: Mapping[str, str],
    metal: gemmi.Atom,
) -> tuple[dict[str, ResidueSelector], dict[str, gemmi.Residue], bool, dict[str, Any]]:
    selectors = {
        role: parse_residue_selector(row[field]) for role, field in ROLE_FIELDS.items()
    }
    residues: dict[str, gemmi.Residue] = {}
    for role in fixed.ROLE_ORDER:
        selector = selectors[role]
        if selector.chain != "A":
            raise HoldoutPreparationError(f"{row['pdb_id']} {role} is not on chain A")
        allowed = fixed.ROLE_ALLOWED_RESNAMES[role]
        if selector.resname not in allowed:
            raise HoldoutPreparationError(
                f"{row['pdb_id']} {role} has unsupported identity {selector.resname}"
            )
        residues[role] = resolve_residue(model, selector)
        if residue_has_altloc(residues[role]):
            raise HoldoutPreparationError(
                f"{row['pdb_id']} core selector {selector.label()} has an alternate conformer"
            )
    catalytic = selectors["catalytic_aspartate"]
    plus2 = selectors["extra_acidic_ligand_homolog"]
    if plus2.chain != catalytic.chain or plus2.resnum != catalytic.resnum + 2:
        raise HoldoutPreparationError(f"{row['pdb_id']} D+2 homolog is not catalytic-Asp+2")
    include_plus2 = parse_bool(row["include_plus2"], f"{row['pdb_id']} include_plus2")
    if include_plus2 != (plus2.resname in {"ASP", "GLU"}):
        raise HoldoutPreparationError(
            f"{row['pdb_id']} D+2 inclusion flag disagrees with residue chemistry"
        )
    for role, field in ROLE_DISTANCE_FIELDS.items():
        expected = parse_float(row[field], field, allow_na=True)
        key = fixed.residue_key(
            {
                "chain": selectors[role].chain,
                "resname": selectors[role].resname,
                "resnum": selectors[role].resnum,
                "icode": selectors[role].icode,
            },
            role,
        )
        distances = fixed.donor_distances(residues[role], key, metal.pos)
        if expected is None:
            if distances:
                raise HoldoutPreparationError(
                    f"{row['pdb_id']} {role} unexpectedly has a metal-facing donor"
                )
        elif not distances:
            raise HoldoutPreparationError(f"{row['pdb_id']} {role} has no expected donor")
        else:
            assert_close(
                distances[0]["distance_A"],
                expected,
                f"{row['pdb_id']} {role} distance",
                AUDIT_DISTANCE_TOLERANCE_A,
            )

    partner_geometry = fixed.asp_partner_geometry(
        residues["catalytic_aspartate"],
        residues["catalytic_asp_cationic_partner"],
    )
    left_raw, right_raw = row["cat_asp_partner_contact"].split("--", 1)
    left = parse_atom_selector(left_raw)
    right = parse_atom_selector(right_raw)
    if (
        left.residue != selectors["catalytic_aspartate"]
        or right.residue != selectors["catalytic_asp_cationic_partner"]
        or partner_geometry["asp_atom"] != left.atom
        or partner_geometry["partner_atom"] != right.atom
    ):
        raise HoldoutPreparationError(f"{row['pdb_id']} frozen Asp/partner contact changed")
    assert_close(
        partner_geometry["heavy_atom_distance_A"],
        float(row["cat_asp_partner_A"]),
        f"{row['pdb_id']} Asp/partner distance",
        AUDIT_DISTANCE_TOLERANCE_A,
    )
    return selectors, residues, include_plus2, partner_geometry


def altloc_inventory(model: gemmi.Model) -> list[str]:
    return sorted(
        residue_selector(chain, residue).label()
        for chain in model
        for residue in chain
        if residue_has_altloc(residue)
    )


def water_inventory(
    model: gemmi.Model, selected_chain: str, metal: gemmi.Atom
) -> dict[str, Any]:
    waters: list[tuple[float, str]] = []
    source_count = 0
    chain_count = 0
    for chain in model:
        for residue in chain:
            if residue.name.upper().strip() not in WATER_NAMES:
                continue
            source_count += 1
            if chain.name == selected_chain:
                chain_count += 1
            for atom in selected_residue_atoms(residue):
                if atom.element.name.upper() != "O":
                    continue
                selector = AtomSelector(
                    residue_selector(chain, residue), atom.name.upper().strip()
                )
                waters.append((metal.pos.dist(atom.pos), selector.label()))
    waters.sort()
    return {
        "source_water_count": source_count,
        "selected_chain_water_count": chain_count,
        "waters_within_3p6_A": sum(distance <= 3.6 for distance, _ in waters),
        "nearest_water_selector": waters[0][1] if waters else None,
        "nearest_water_metal_A": waters[0][0] if waters else None,
    }


def validate_water_inventory(row: Mapping[str, str], observed: Mapping[str, Any]) -> dict[str, Any]:
    for field in ("source_water_count", "selected_chain_water_count", "waters_within_3p6_A"):
        if observed[field] != int(row[field]):
            raise HoldoutPreparationError(
                f"{row['pdb_id']} {field} changed: {observed[field]} != {row[field]}"
            )
    expected_nearest = None if row["nearest_water_selector"] == "NA" else row["nearest_water_selector"]
    if observed["nearest_water_selector"] != expected_nearest:
        raise HoldoutPreparationError(f"{row['pdb_id']} nearest-water selector changed")
    expected_distance = parse_float(
        row["nearest_water_metal_A"], "nearest_water_metal_A", allow_na=True
    )
    if expected_distance is None:
        if observed["nearest_water_metal_A"] is not None:
            raise HoldoutPreparationError(f"{row['pdb_id']} unexpectedly has a nearest water")
    else:
        assert_close(
            float(observed["nearest_water_metal_A"]),
            expected_distance,
            f"{row['pdb_id']} nearest-water distance",
            AUDIT_DISTANCE_TOLERANCE_A,
        )
    return {
        **observed,
        "nearest_water_metal_A": (
            round(float(observed["nearest_water_metal_A"]), 3)
            if observed["nearest_water_metal_A"] is not None
            else None
        ),
        "all_source_waters_excluded_from_normalized_and_QM_structures": True,
        "synthetic_waters_added": 0,
    }


def noncore_heterogen_contacts(
    model: gemmi.Model,
    metal_selector: AtomSelector,
    pqq_selector: ResidueSelector,
    metal: gemmi.Atom,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for chain in model:
        for residue in chain:
            selector = residue_selector(chain, residue)
            if (
                selector.resname in STANDARD_AA
                or selector.resname in WATER_NAMES
                or selector == pqq_selector
                or selector == metal_selector.residue
            ):
                continue
            for atom in selected_residue_atoms(residue):
                element = atom.element.name.upper()
                if element not in {"O", "N", "S"}:
                    continue
                distance = metal.pos.dist(atom.pos)
                if distance <= COORDINATION_CUTOFF_A:
                    records.append(
                        {
                            "selector": AtomSelector(selector, atom.name.upper().strip()),
                            "element": element.title(),
                            "distance_A": distance,
                        }
                    )
    records.sort(key=lambda item: (item["distance_A"], item["selector"].label()))
    return records


def validate_noncore_contact(
    row: Mapping[str, str], observed: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    raw = row["noncore_direct_ligand"]
    expected = [] if raw == "none" else [parse_atom_selector(raw)]
    if [item["selector"] for item in observed] != expected:
        raise HoldoutPreparationError(
            f"{row['pdb_id']} noncore direct-ligand inventory changed"
        )
    expected_distance = parse_float(
        row["noncore_direct_ligand_metal_A"],
        "noncore_direct_ligand_metal_A",
        allow_na=True,
    )
    if expected:
        assert expected_distance is not None
        assert_close(
            float(observed[0]["distance_A"]),
            expected_distance,
            f"{row['pdb_id']} noncore ligand distance",
            AUDIT_DISTANCE_TOLERANCE_A,
        )
    elif expected_distance is not None:
        raise HoldoutPreparationError(f"{row['pdb_id']} has distance but no noncore ligand")
    return [
        {
            "selector": item["selector"].label(),
            "element": item["element"],
            "distance_A": round(float(item["distance_A"]), 3),
            "disposition": "excluded_without_replacement",
        }
        for item in observed
    ]


def selected_atom_identity(
    chain: str,
    residue: gemmi.Residue,
    atom: gemmi.Atom,
    *,
    normalized_metal_residue: bool,
) -> tuple[str, int, str, str, str, str]:
    resname = residue.name.upper().strip()
    atom_name = atom.name.upper().strip()
    element = atom.element.name.upper()
    if normalized_metal_residue:
        # The standard-only wrapper keys its allowlist by residue name.  Keep
        # the native Ca/Ce/La atom name and element through protonation so this
        # preprocessing step changes no physical atom identity; only the
        # nonstandard residue label is canonicalized to LA.
        resname = "LA"
    return (
        chain,
        residue.seqid.num,
        insertion_code(residue),
        resname,
        atom_name,
        "H" if element == "D" else element,
    )


def coordinate_digest(
    atoms: Mapping[tuple[str, int, str, str, str, str], tuple[float, float, float]]
) -> str:
    serial = [
        [*identity, *[round(float(value), 6) for value in atoms[identity]]]
        for identity in sorted(atoms)
    ]
    return hashlib.sha256(
        json.dumps(serial, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


def compare_heavy_maps(
    reference: Mapping[tuple[str, int, str, str, str, str], tuple[float, float, float]],
    query: Mapping[tuple[str, int, str, str, str, str], tuple[float, float, float]],
    label: str,
) -> dict[str, Any]:
    if set(reference) != set(query):
        missing = sorted(set(reference) - set(query))
        added = sorted(set(query) - set(reference))
        raise HoldoutPreparationError(
            f"{label} changed heavy identities: missing={missing[:5]}, added={added[:5]}"
        )
    displacements = {key: math.dist(reference[key], query[key]) for key in reference}
    if not displacements:
        raise HoldoutPreparationError(f"{label} has no retained heavy atoms")
    worst, maximum = max(displacements.items(), key=lambda item: item[1])
    if maximum > COORDINATE_TOLERANCE_A:
        raise HoldoutPreparationError(
            f"{label} moved {worst} by {maximum:.8f} A beyond {COORDINATE_TOLERANCE_A:.8f} A"
        )
    return {
        "label": label,
        "atom_count": len(reference),
        "identity_sets_equal": True,
        "reference_coordinate_digest": coordinate_digest(reference),
        "query_coordinate_digest": coordinate_digest(query),
        "maximum_displacement_A": maximum,
        "maximum_displacement_atom": list(worst),
        "allowed_serialization_tolerance_A": COORDINATE_TOLERANCE_A,
        "passes": True,
    }


def heavy_map_from_file(path: Path) -> dict[tuple[str, int, str, str, str, str], tuple[float, float, float]]:
    structure = gemmi.read_structure(str(path))
    model = exact_model(structure, "1")
    chain = exact_chain(model, "A")
    result: dict[tuple[str, int, str, str, str, str], tuple[float, float, float]] = {}
    for residue in chain:
        if residue_has_altloc(residue):
            raise HoldoutPreparationError(f"normalized/protonated residue has altloc: {residue.name}")
        for atom in residue:
            element = atom.element.name.upper()
            if element in {"H", "D"}:
                continue
            key = selected_atom_identity(
                chain.name, residue, atom, normalized_metal_residue=False
            )
            if key in result:
                raise HoldoutPreparationError(f"duplicate heavy identity in {path}: {key}")
            result[key] = (atom.pos.x, atom.pos.y, atom.pos.z)
    return result


def clone_atom(atom: gemmi.Atom) -> gemmi.Atom:
    clone = gemmi.Atom()
    clone.name = atom.name
    clone.element = gemmi.Element(atom.element.name)
    clone.pos = gemmi.Position(atom.pos.x, atom.pos.y, atom.pos.z)
    clone.occ = float(atom.occ)
    clone.b_iso = float(atom.b_iso)
    clone.altloc = "\x00"
    clone.serial = int(atom.serial)
    return clone


def clone_residue_selected(
    residue: gemmi.Residue,
    *,
    normalize_to_la: bool,
) -> gemmi.Residue:
    clone = gemmi.Residue()
    clone.name = "LA" if normalize_to_la else residue.name.upper().strip()
    clone.seqid = gemmi.SeqId(residue.seqid.num, insertion_code(residue) or " ")
    clone.het_flag = "H" if normalize_to_la or clone.name == "PQQ" else residue.het_flag
    atoms = selected_residue_atoms(residue)
    if normalize_to_la and len(atoms) != 1:
        raise HoldoutPreparationError("selected source metal residue does not contain one atom")
    for atom in atoms:
        clone.add_atom(clone_atom(atom))
    return clone


def write_selected_normalized_structure(
    source: gemmi.Structure,
    model: gemmi.Model,
    chain: gemmi.Chain,
    metal_selector: AtomSelector,
    pqq_selector: ResidueSelector,
    output: Path,
) -> tuple[
    dict[tuple[str, int, str, str, str, str], tuple[float, float, float]],
    dict[str, Any],
]:
    selected = gemmi.Structure()
    selected.name = source.name
    selected.cell = source.cell
    selected.spacegroup_hm = source.spacegroup_hm
    out_model = gemmi.Model(1)
    out_chain = gemmi.Chain("A")
    source_heavy: dict[
        tuple[str, int, str, str, str, str], tuple[float, float, float]
    ] = {}
    excluded: list[dict[str, Any]] = []
    altloc_choices: list[dict[str, Any]] = []
    retained_residues = 0
    retained_source_hydrogens = 0

    for residue in chain:
        selector = residue_selector(chain, residue)
        is_metal = selector == metal_selector.residue
        is_pqq = selector == pqq_selector
        retain = selector.resname in STANDARD_AA or is_metal or is_pqq
        if not retain:
            excluded.append(
                {
                    "selector": selector.label(),
                    "residue_name": selector.resname,
                    "atom_count": len(list(residue)),
                    "atom_names": [atom.name.upper().strip() for atom in residue],
                    "classification": (
                        "source_water" if selector.resname in WATER_NAMES else "noncore_heterogen"
                    ),
                }
            )
            continue
        if selector.resname in {"PQQ"} and not is_pqq:
            raise HoldoutPreparationError("an unselected PQQ occurs on selected chain A")
        if residue_has_altloc(residue):
            altloc_choices.append(
                {
                    "selector": selector.label(),
                    **altloc_decision(residue),
                }
            )
        selected_atoms = selected_residue_atoms(residue)
        retained_source_hydrogens += sum(
            atom.element.name.upper() in {"H", "D"} for atom in selected_atoms
        )
        for atom in selected_atoms:
            if atom.element.name.upper() in {"H", "D"}:
                continue
            key = selected_atom_identity(
                chain.name, residue, atom, normalized_metal_residue=is_metal
            )
            if key in source_heavy:
                raise HoldoutPreparationError(f"duplicate retained source heavy identity {key}")
            source_heavy[key] = (atom.pos.x, atom.pos.y, atom.pos.z)
        out_chain.add_residue(
            clone_residue_selected(residue, normalize_to_la=is_metal)
        )
        retained_residues += 1

    if retained_source_hydrogens:
        raise HoldoutPreparationError(
            "retained raw source already contains hydrogens; deterministic v3 protonation requires none"
        )
    out_model.add_chain(out_chain)
    selected.add_model(out_model)
    output.parent.mkdir(parents=True, exist_ok=True)
    selected.write_pdb(str(output))
    normalized_heavy = heavy_map_from_file(output)
    check = compare_heavy_maps(
        source_heavy, normalized_heavy, "retained_raw_source_to_selected_normalized_PDB"
    )
    normalized_structure = gemmi.read_structure(str(output))
    normalized_model = exact_model(normalized_structure, "1")
    normalized_chain = exact_chain(normalized_model, "A")
    observed_names = [residue.name.upper().strip() for residue in normalized_chain]
    unsupported = sorted(
        name for name in set(observed_names) if name not in STANDARD_AA | {"PQQ", "LA"}
    )
    if unsupported or observed_names.count("PQQ") != 1 or observed_names.count("LA") != 1:
        raise HoldoutPreparationError(
            f"selected normalized structure has unsupported identities {unsupported} "
            f"or PQQ/LA counts {observed_names.count('PQQ')}/{observed_names.count('LA')}"
        )
    return source_heavy, {
        "retained_chain": "A",
        "retained_residue_count": retained_residues,
        "retained_heavy_atom_count": len(source_heavy),
        "retained_source_hydrogen_count": 0,
        "excluded_selected_chain_residue_count": len(excluded),
        "excluded_selected_chain_residues": excluded,
        "excluded_other_models": [
            str(item.num) for item in source if item.num != model.num
        ],
        "excluded_other_chains": [
            item.name for item in model if item.name != chain.name
        ],
        "alternate_conformer_selections": altloc_choices,
        "metal_normalization": {
            "source": metal_selector.label(),
            "normalized": (
                f"A:LA{metal_selector.residue.resnum}{metal_selector.residue.icode}/"
                f"{metal_selector.atom}"
            ),
            "native_atom_name_retained": metal_selector.atom,
            "native_element_retained": True,
            "residue_name_only_changed_to_LA": True,
            "coordinate_changed": False,
        },
        "heavy_coordinate_check": check,
    }


def audit_raw_source(row: Mapping[str, str]) -> dict[str, Any]:
    source = Path(row["source_path"]).resolve()
    if sha256_file(source) != row["source_sha256"]:
        raise HoldoutPreparationError(f"{row['pdb_id']} raw source SHA-256 mismatch")
    expected_format = "PDB" if source.suffix.lower() == ".pdb" else "mmCIF"
    if row["source_format"] != expected_format:
        raise HoldoutPreparationError(f"{row['pdb_id']} source format/path disagree")
    structure = gemmi.read_structure(str(source))
    expected_resolution = float(row["resolution_A"])
    observed_resolution = float(structure.resolution)
    assert_close(
        observed_resolution,
        expected_resolution,
        f"{row['pdb_id']} source resolution",
        0.011,
    )
    model = exact_model(structure, row["model"])
    chain = exact_chain(model, row["selected_chain"])
    metal_selector = parse_atom_selector(row["metal_selector"])
    pqq_selector = parse_residue_selector(row["pqq_selector"])
    if metal_selector.residue.chain != "A" or pqq_selector.chain != "A":
        raise HoldoutPreparationError(f"{row['pdb_id']} frozen cofactors are not on chain A")
    metal_residue = resolve_residue(model, metal_selector.residue)
    metal = resolve_atom(model, metal_selector)
    if residue_has_altloc(metal_residue):
        raise HoldoutPreparationError(f"{row['pdb_id']} selected metal has an altloc")
    if len(selected_residue_atoms(metal_residue)) != 1:
        raise HoldoutPreparationError(f"{row['pdb_id']} selected metal residue is not monatomic")
    if metal.element.name.upper() != row["metal_element"].upper():
        raise HoldoutPreparationError(f"{row['pdb_id']} source metal element changed")
    assert_close(
        float(metal.occ), float(row["metal_occupancy"]),
        f"{row['pdb_id']} metal occupancy", 0.001,
    )
    assert_close(
        float(metal.b_iso), float(row["metal_bfactor_A2"]),
        f"{row['pdb_id']} metal B factor", 0.011,
    )
    expected_xyz = [float(value) for value in row["metal_xyz_A"].split(",")]
    observed_xyz = [metal.pos.x, metal.pos.y, metal.pos.z]
    if len(expected_xyz) != 3:
        raise HoldoutPreparationError(f"{row['pdb_id']} malformed metal coordinate")
    for axis, observed, expected in zip("xyz", observed_xyz, expected_xyz, strict=True):
        assert_close(observed, expected, f"{row['pdb_id']} metal {axis}", 0.001)

    pqq_residue = resolve_residue(model, pqq_selector)
    if residue_has_altloc(pqq_residue):
        raise HoldoutPreparationError(f"{row['pdb_id']} selected PQQ has an altloc")
    if any(atom.element.name.upper() in {"H", "D"} for atom in pqq_residue):
        raise HoldoutPreparationError(f"{row['pdb_id']} raw PQQ unexpectedly contains H/D")
    prepared_pqq = prepare_pqq(pqq_residue, DEFAULT_PQQ_MICROSTATE)
    if (
        prepared_pqq.schema_id != row["pqq_schema"]
        or prepared_pqq.schema_id != "pdb_ccd_pqq_v1"
        or prepared_pqq.microstate.formal_charge != -3
        or prepared_pqq.formula != "C14H3N2O8"
    ):
        raise HoldoutPreparationError(f"{row['pdb_id']} PQQ schema/microstate changed")
    pqq_on_chain = [residue for residue in chain if residue.name.upper().strip() == "PQQ"]
    if len(pqq_on_chain) != 1:
        raise HoldoutPreparationError(f"{row['pdb_id']} chain A does not have exactly one PQQ")

    role_selectors, _, include_plus2, partner_geometry = validate_role_selectors(
        model, row, metal
    )
    if row["core_altlocs"] != "none":
        raise HoldoutPreparationError(f"{row['pdb_id']} frozen core-altloc field is not none")
    observed_altlocs = altloc_inventory(model)
    expected_altlocs = [] if row["source_altloc_residues"] == "none" else sorted(
        row["source_altloc_residues"].split(";")
    )
    if observed_altlocs != expected_altlocs:
        raise HoldoutPreparationError(
            f"{row['pdb_id']} source altloc inventory changed: {observed_altlocs}"
        )

    donor_records = validate_donor_ledger(
        row, typed_source_contacts(model, metal_selector, metal)
    )
    included_core = {pqq_selector}
    included_core.update(
        selector
        for role, selector in role_selectors.items()
        if role != "extra_acidic_ligand_homolog" or include_plus2
    )
    unexpected = sorted(
        {
            parse_atom_selector(record["selector"]).residue.label()
            for record in donor_records
            if parse_atom_selector(record["selector"]).residue not in included_core
        }
    )
    if unexpected:
        raise HoldoutPreparationError(
            f"{row['pdb_id']} qualifying donor lies outside frozen core: {unexpected}"
        )
    water_record = validate_water_inventory(
        row, water_inventory(model, row["selected_chain"], metal)
    )
    if water_record["waters_within_3p6_A"] != 0:
        raise HoldoutPreparationError(f"{row['pdb_id']} has a water inside the frozen dry vacancy")
    noncore_contacts = validate_noncore_contact(
        row,
        noncore_heterogen_contacts(model, metal_selector, pqq_selector, metal),
    )
    if row["pdb_id"] == "4MAE":
        if [record["selector"] for record in noncore_contacts] != ["A:15P603/OXT"]:
            raise HoldoutPreparationError("4MAE 15P603/OXT exclusion is not explicit")
        if "15P603" not in row["dry_policy"] or "do not fill" not in row["dry_policy"]:
            raise HoldoutPreparationError("4MAE dry vacancy policy changed")
    elif noncore_contacts:
        raise HoldoutPreparationError(f"{row['pdb_id']} has an unexpected noncore ligand")

    return {
        "source": source,
        "structure": structure,
        "model": model,
        "chain": chain,
        "metal_selector": metal_selector,
        "metal": metal,
        "pqq_selector": pqq_selector,
        "role_selectors": role_selectors,
        "include_plus2": include_plus2,
        "partner_geometry": partner_geometry,
        "donor_records": donor_records,
        "water_record": water_record,
        "noncore_contacts": noncore_contacts,
        "source_altlocs": observed_altlocs,
        "source_metal": {
            "selector": metal_selector.label(),
            "element": metal.element.name.title(),
            "occupancy": float(metal.occ),
            "b_iso_A2": float(metal.b_iso),
            "xyz_A": [round(value, 6) for value in observed_xyz],
        },
        "pqq": prepared_pqq.metadata(),
    }


def validate_protonation_manifest(
    manifest_path: Path,
    normalized: Path,
    protonated: Path,
    calibration_pins: Mapping[str, Any],
) -> dict[str, Any]:
    manifest = read_object(manifest_path)
    subprotocol = calibration_pins["protonation_subprotocol"]
    policy = manifest.get("nonstandard_definition_policy")
    if (
        manifest.get("protocol_id") != "pdbfixer_standard_residue_protonation_v2"
        or manifest.get("experiment_protonation_protocol_id") != subprotocol["id"]
        or manifest.get("canonical_protonator")
        != calibration_pins["canonical_helpers"]["protonate_cif"]
        or manifest.get("experiment_wrapper")
        != calibration_pins["experiment_helpers"]["protonate_standard_only"]
        or manifest.get("source") != file_record(normalized)
        or manifest.get("output") != file_record(protonated)
        or manifest.get("ph") != 7.0
        or manifest.get("add_missing_residues") is not False
        or manifest.get("repaired_missing_atom_count") != 0
        or manifest.get("repaired_missing_terminal_atom_count") != 0
        or manifest.get("cofactor_protonation") != "not_assigned"
        or manifest.get("software") != calibration_pins["preparation_runtime"]["packages"]
        or not isinstance(policy, dict)
        or policy.get("policy_id") != subprotocol["nonstandard_definition_policy_id"]
        or policy.get("observed_nonstandard_residue_counts") != {"LA": 1, "PQQ": 1}
        or policy.get("output_nonstandard_hydrogen_counts") != {"LA": 0, "PQQ": 0}
        or policy.get("python_random_seed") != subprotocol["python_random_seed"]
        or policy.get("openmm_platform") != "CPU"
        or policy.get("openmm_CPU_Threads") != 1
        or policy.get("forcefield") is not None
        or policy.get("scientific_chemistry_changed") is not False
    ):
        raise HoldoutPreparationError("standard-AA protonation manifest violates v3")
    return manifest


def fixed_key(selector: ResidueSelector, role: str) -> base.ResidueKey:
    return fixed.residue_key(
        {
            "chain": selector.chain,
            "resname": selector.resname,
            "resnum": selector.resnum,
            "icode": selector.icode,
        },
        role,
    )


def contact_record(contact: base.Contact) -> dict[str, Any]:
    record = contact.metadata()
    record.pop("included_in_qm_at_3p3A", None)
    record["included_in_fixed_qm_core"] = True
    return record


def validate_fragment_heavy_atoms_against_raw(
    fragments: Sequence[Mapping[str, Any]], raw_model: gemmi.Model
) -> dict[str, Any]:
    fragment_records: list[dict[str, Any]] = []
    maximum = 0.0
    total = 0
    for fragment in fragments:
        selector = parse_residue_selector(str(fragment.get("id", "")))
        source_residue = resolve_residue(raw_model, selector)
        source_atoms = {
            atom.name.upper().strip(): atom
            for atom in selected_residue_atoms(source_residue)
            if atom.element.name.upper() not in {"H", "D"}
        }
        records = fragment.get("atom_records")
        if not isinstance(records, list):
            raise HoldoutPreparationError(f"fragment {selector.label()} lacks atom provenance")
        checked_names: list[str] = []
        for record in records:
            if not isinstance(record, dict):
                raise HoldoutPreparationError(f"fragment {selector.label()} has malformed provenance")
            element = str(record.get("element", "")).upper()
            if element == "H":
                continue
            if record.get("origin") != "source_heavy_atom":
                raise HoldoutPreparationError(
                    f"fragment {selector.label()} has a non-source heavy atom"
                )
            atom_name = str(record.get("name", "")).upper()
            source_atom = source_atoms.get(atom_name)
            xyz = record.get("xyz_A")
            if source_atom is None or not isinstance(xyz, list) or len(xyz) != 3:
                raise HoldoutPreparationError(
                    f"fragment {selector.label()} heavy atom {atom_name} is absent from raw source"
                )
            if source_atom.element.name.upper() != element:
                raise HoldoutPreparationError(
                    f"fragment {selector.label()} heavy element changed for {atom_name}"
                )
            displacement = math.dist(
                (source_atom.pos.x, source_atom.pos.y, source_atom.pos.z),
                tuple(float(value) for value in xyz),
            )
            if displacement > 1.0e-5:
                raise HoldoutPreparationError(
                    f"fragment {selector.label()}/{atom_name} moved from raw source by "
                    f"{displacement:.8f} A"
                )
            maximum = max(maximum, displacement)
            total += 1
            checked_names.append(atom_name)
        if len(checked_names) != len(set(checked_names)):
            raise HoldoutPreparationError(
                f"fragment {selector.label()} duplicates a raw heavy atom"
            )
        fragment_records.append(
            {
                "fragment": selector.label(),
                "raw_heavy_atom_names": checked_names,
                "raw_heavy_atom_count": len(checked_names),
            }
        )
    return {
        "all_final_fragment_heavy_atoms_match_raw_source": True,
        "checked_heavy_atom_count": total,
        "maximum_displacement_A": maximum,
        "allowed_tolerance_A": 1.0e-5,
        "fragments": fragment_records,
    }


def write_holdout_pair(
    protonated_path: Path,
    target_dir: Path,
    stem: str,
    row: Mapping[str, str],
    raw: Mapping[str, Any],
    normalization_manifest: Path,
    protonation_manifest: Path,
    heavy_check_path: Path,
    gate: Mapping[str, Any],
    holdout_pins_path: Path,
    holdout_pins: Mapping[str, Any],
    calibration_pins_path: Path,
    calibration_pins: Mapping[str, Any],
) -> Path:
    structure = gemmi.read_structure(str(protonated_path))
    model = exact_model(structure, "1")
    exact_chain(model, "A")
    if any(
        residue.name.upper().strip() in WATER_NAMES
        for chain in model
        for residue in chain
    ):
        raise HoldoutPreparationError("a source water survived into holdout carving")
    unsupported = sorted(
        residue_selector(chain, residue).label()
        for chain in model
        for residue in chain
        if residue.name.upper().strip() not in STANDARD_AA | {"PQQ", "LA"}
    )
    if unsupported:
        raise HoldoutPreparationError(f"noncore heterogens survived selection: {unsupported}")

    source_metal_selector: AtomSelector = raw["metal_selector"]
    normalized_metal = AtomSelector(
        ResidueSelector("A", "LA", source_metal_selector.residue.resnum, source_metal_selector.residue.icode),
        source_metal_selector.atom,
    )
    sites = base.select_metal_sites(
        model,
        site_chain="A",
        site_resnum=normalized_metal.residue.resnum,
        site_icode=normalized_metal.residue.icode,
        site_resname="LA",
        site_atom=normalized_metal.atom,
    )
    site = sites[0]
    prepared_pqq = base._prepare_nearby_pqq(model, site.atom.pos, DEFAULT_PQQ_MICROSTATE)
    base._reject_unsupported_nearby_species(model, site, prepared_pqq)
    pqq_selector: ResidueSelector = raw["pqq_selector"]
    pqq_key = fixed_key(pqq_selector, "PQQ")
    if set(prepared_pqq) != {pqq_key}:
        raise HoldoutPreparationError(
            f"{row['pdb_id']} protonated input does not contain exactly selected PQQ"
        )
    pqq_residue, prepared = prepared_pqq[pqq_key]
    if any(atom.element.name.upper() in {"H", "D"} for atom in pqq_residue):
        raise HoldoutPreparationError("PDBFixer supplied PQQ hydrogen")
    if (
        prepared.schema_id != "pdb_ccd_pqq_v1"
        or prepared.microstate.microstate_id != DEFAULT_PQQ_MICROSTATE
        or prepared.microstate.formal_charge != -3
        or prepared.formula != "C14H3N2O8"
    ):
        raise HoldoutPreparationError("holdout PQQ chemistry differs from frozen PQQ(3-)")

    role_selectors: Mapping[str, ResidueSelector] = raw["role_selectors"]
    role_keys = {role: fixed_key(selector, role) for role, selector in role_selectors.items()}
    role_residues = {
        role: fixed.resolve_residue(model, key, role) for role, key in role_keys.items()
    }
    partner_geometry = fixed.asp_partner_geometry(
        role_residues["catalytic_aspartate"],
        role_residues["catalytic_asp_cationic_partner"],
    )
    if partner_geometry != raw["partner_geometry"]:
        raise HoldoutPreparationError("protonation changed frozen Asp/partner geometry")

    qualified, _, nearby_sulfur, nearby_untyped = base._classify_contacts(
        model, site, prepared_pqq
    )
    observed_donors = [
        {
            "selector": AtomSelector(
                ResidueSelector(
                    contact.residue.chain,
                    contact.residue.resname,
                    contact.residue.resnum,
                    contact.residue.icode.strip(),
                ),
                contact.atom_name,
            ),
            "element": contact.element,
            "donor_type": contact.donor_type,
            "distance_A": contact.distance_A,
        }
        for contact in qualified
    ]
    validate_donor_ledger(row, observed_donors)
    include_plus2 = bool(raw["include_plus2"])
    included_roles = {
        "anchor_glutamate",
        "anchor_asparagine",
        "catalytic_aspartate",
        "catalytic_asp_cationic_partner",
    }
    if include_plus2:
        included_roles.add("extra_acidic_ligand_homolog")
    included_keys = {pqq_key, *(role_keys[role] for role in included_roles)}
    unexpected = sorted(
        contact.residue.label() for contact in qualified if contact.residue not in included_keys
    )
    if unexpected:
        raise HoldoutPreparationError(f"protonated holdout has off-core direct donor: {unexpected}")

    pqq_records = fixed.pqq_atom_records(pqq_residue, prepared)
    fragments: list[dict[str, Any]] = [
        {
            "kind": "fixed_core_pqq",
            "role": "pqq_cofactor",
            "id": pqq_key.label(),
            "formal_charge": -3,
            "atom_count": len(prepared.atoms),
            "atom_records": pqq_records,
            "microstate_id": DEFAULT_PQQ_MICROSTATE,
            "selection": "exact_frozen_complete_CCD_cofactor",
        }
    ]
    persistent_atoms: list[tuple[str, float, float, float]] = list(prepared.atoms)
    excluded_plus2: dict[str, Any] | None = None
    for role in fixed.ROLE_ORDER:
        if role == "extra_acidic_ligand_homolog" and not include_plus2:
            excluded_plus2 = {
                "role": role,
                "id": role_keys[role].label(),
                "resname": role_keys[role].resname,
                "included": False,
                "reason": "nonacidic_homolog_has_no_XoxF_extra_carboxylate",
            }
            continue
        if role == "catalytic_asp_cationic_partner":
            atoms, fragment = fixed.cationic_sidechain_fragment(
                role_residues[role], role_keys[role], site.atom.pos
            )
        else:
            atoms, fragment = fixed.canonical_sidechain_fragment(
                role_residues[role], role_keys[role], role, site.atom.pos
            )
        persistent_atoms.extend(atoms)
        fragments.append(fragment)

    scaffold_charge = sum(int(fragment["formal_charge"]) for fragment in fragments)
    final_fragment_heavy_check = validate_fragment_heavy_atoms_against_raw(
        fragments, raw["model"]
    )
    charges = {"La": scaffold_charge + 3, "Ca": scaffold_charge + 2}
    expected_charges = {"La": -2, "Ca": -3} if include_plus2 else {"La": -1, "Ca": -2}
    if charges != expected_charges:
        raise HoldoutPreparationError(
            f"{row['pdb_id']} charge ledger {charges} differs from {expected_charges}"
        )
    la_atoms = [("La", site.atom.pos.x, site.atom.pos.y, site.atom.pos.z), *persistent_atoms]
    ca_atoms = [("Ca", site.atom.pos.x, site.atom.pos.y, site.atom.pos.z), *persistent_atoms]
    if la_atoms[1:] != ca_atoms[1:]:
        raise HoldoutPreparationError("in-memory La/Ca nonmetal coordinates differ")
    base._validate_singlet("La", la_atoms, charges["La"])
    base._validate_singlet("Ca", ca_atoms, charges["Ca"])

    outputs: dict[str, dict[str, str]] = {}
    tasks: list[dict[str, Any]] = []
    for label, atoms in (("La", la_atoms), ("Ca", ca_atoms)):
        xyz_path = target_dir / f"{stem}_{label}_qm.xyz"
        input_path = target_dir / f"sp_{stem}_{label}.inp"
        output_path = target_dir / f"sp_{stem}_{label}.out"
        fixed.write_xyz(xyz_path, stem=stem, label=label, atoms=atoms, charge=charges[label])
        fixed.write_orca_input(
            input_path,
            xyz_name=xyz_path.name,
            charge=charges[label],
            stem=stem,
            label=label,
        )
        outputs[f"{label}_xyz"] = file_record(xyz_path)
        outputs[f"{label}_input"] = file_record(input_path)
        tasks.append(
            {
                "task_id": label,
                "input": {"path": input_path.name, "sha256": sha256_file(input_path)},
                "xyz": {"path": xyz_path.name, "sha256": sha256_file(xyz_path)},
                "output_path": output_path.name,
            }
        )
    la_lines = (target_dir / f"{stem}_La_qm.xyz").read_text().splitlines()
    ca_lines = (target_dir / f"{stem}_Ca_qm.xyz").read_text().splitlines()
    if (
        len(la_lines) < 4
        or len(ca_lines) < 4
        or la_lines[2].split()[0].upper() != "LA"
        or ca_lines[2].split()[0].upper() != "CA"
        or la_lines[2].split()[1:] != ca_lines[2].split()[1:]
        or la_lines[3:] != ca_lines[3:]
    ):
        raise HoldoutPreparationError("serialized La/Ca nonmetal payloads are not byte-identical")
    for label in ("La", "Ca"):
        text = (target_dir / f"sp_{stem}_{label}.inp").read_text()
        if EXPECTED_ORCA_DIRECTIVE not in text or "%pal" in text.lower():
            raise HoldoutPreparationError("ORCA template differs from frozen electronic model")
    if list(target_dir.glob("*.out")):
        raise HoldoutPreparationError("holdout preparation directory unexpectedly contains ORCA output")

    special_limitation = (
        "A:15P603/OXT at 2.747 A was deliberately removed without replacement; "
        "this is a dry fixed-coordinate structural-transfer test with a ligand vacancy, "
        "not the intact crystallographic first shell."
        if row["pdb_id"] == "4MAE"
        else None
    )
    manifest = {
        "schema_version": HOLDOUT_CARVE_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "ready_for_orca_after_passing_calibration_gate",
        "stem": stem,
        "panel_id": row["pdb_id"],
        "holdout": {
            "role": row["holdout_role"],
            "biological_class": row["biological_class"],
            "primary_test_member": row["holdout_role"] == "primary",
            "secondary_geometry_check_only": row["holdout_role"] == "secondary",
            "special_interpretation_limitation": special_limitation,
        },
        "calibration_gate": gate,
        "experiment_provenance": {
            "holdout_pins": file_record(holdout_pins_path),
            "calibration_implementation_pins": file_record(calibration_pins_path),
            "selector_authority": holdout_pins["holdout_spec"],
            "selector_audit": holdout_pins["holdout_selector_audit"],
            "preparer": holdout_pins["holdout_implementation"]["prepare_holdouts"],
        },
        "raw_source": {"path": str(raw["source"]), "sha256": row["source_sha256"]},
        "source_selection": {
            "model": "1",
            "chain": "A",
            "selection_fallback_allowed": False,
            "source_metal": raw["source_metal"],
            "pqq_selector": pqq_selector.label(),
            "pqq_schema": prepared.schema_id,
        },
        "normalization_manifest": file_record(normalization_manifest),
        "protonation_manifest": file_record(protonation_manifest),
        "heavy_coordinate_check": file_record(heavy_check_path),
        "selected_site": {
            "normalized_selector": normalized_metal.label(),
            "xyz_A": [round(site.atom.pos.x, 6), round(site.atom.pos.y, 6), round(site.atom.pos.z, 6)],
        },
        "coordination_policy": {
            "policy_id": COORDINATION_POLICY_ID,
            "cutoff_A": COORDINATION_CUTOFF_A,
            "minimum_coordination_number": 6,
            "maximum_direct_nitrogen": MAX_DIRECT_N,
            "alternate_conformer_policy": ALTLOC_POLICY_ID,
            "fragment_membership_is_distance_independent": True,
        },
        "coordination": {
            "coordination_number": len(qualified),
            "oxygen_count": sum(contact.element.upper() == "O" for contact in qualified),
            "direct_nitrogen_count": sum(contact.element.upper() == "N" for contact in qualified),
            "source_typed_direct_donors_in_frozen_tsv_order": raw["donor_records"],
            "protonated_typed_direct_donors_ordered_by_distance": [
                contact_record(contact) for contact in qualified
            ],
            "nearby_sulfur": [contact.metadata() for contact in nearby_sulfur],
            "nearby_untyped_ON": [contact.metadata() for contact in nearby_untyped],
            "unexpected_direct_donors": [],
        },
        "fixed_core": {
            "policy_id": CORE_POLICY_ID,
            "requested_roles": {
                role: selector.label() for role, selector in role_selectors.items()
            },
            "included_roles": sorted(included_roles),
            "excluded_nonacidic_Dplus2_homolog": excluded_plus2,
            "catalytic_Asp_partner_geometry": partner_geometry,
            "water_policy": "dry_exclude_all_source_and_synthetic_waters",
            "noncore_heterogen_policy": "exclude_all_without_replacement",
            "source_water_inventory": raw["water_record"],
            "excluded_noncore_direct_ligands": raw["noncore_contacts"],
            "synthetic_water_count": 0,
            "replacement_ligand_count": 0,
            "point_charge_embedding": False,
            "geometry_relaxation": False,
            "final_fragment_heavy_coordinate_check": final_fragment_heavy_check,
        },
        "qm_fragments": fragments,
        "pqq": {
            "present": True,
            "schema_id": prepared.schema_id,
            "microstate_id": prepared.microstate.microstate_id,
            "formal_charge": prepared.microstate.formal_charge,
            "multiplicity": prepared.microstate.multiplicity,
            "formula": prepared.formula,
        },
        "charge_ledger": {
            "fragments": fragments,
            "scaffold_formal_charge": scaffold_charge,
            "expected_total_charges": expected_charges,
            "La_total": charges["La"],
            "Ca_total": charges["Ca"],
            "multiplicity": 1,
        },
        "paired_arm_invariant": {
            "nonmetal_coordinates_byte_identical": True,
            "arms_differ_only_in_metal_identity_charge_and_electron_count": True,
        },
        "electronic_structure": {
            "method_id": METHOD_ID,
            "method": "r2SCAN-3c",
            "basis_policy": "native_r2scan3c_def2_mTZVPP_with_native_ECP",
            "solvation": "CPCM(Water)",
            "grid": "DefGrid3",
            "single_point_only": True,
        },
        "execution_policy": {
            "id": EXECUTION_POLICY_ID,
            "task_runner": calibration_pins["canonical_helpers"]["run_orca_task_manifest"],
            "runtime_renderer": calibration_pins["canonical_helpers"]["render_orca_runtime_input"],
            "orca_executable": calibration_pins["orca_runtime"]["executable"],
            "omp_threads_per_rank": 1,
        },
        "tasks": tasks,
        "outputs": outputs,
    }
    manifest_path = target_dir / f"{stem}_carve_manifest.json"
    write_json_atomic(manifest_path, manifest)
    return manifest_path


def prepare_one(
    row: Mapping[str, str],
    index: int,
    out_root: Path,
    gate: Mapping[str, Any],
    holdout_pins_path: Path,
    holdout_pins: Mapping[str, Any],
    calibration_pins_path: Path,
    calibration_pins: Mapping[str, Any],
) -> dict[str, Any]:
    raw = audit_raw_source(row)
    pdb_id = row["pdb_id"]
    stem = f"pmdh_fc_holdout_{pdb_id.lower()}"
    target_dir = out_root / f"{index:02d}_{pdb_id}"
    target_dir.mkdir()
    normalized = target_dir / f"{stem}_selected_normalized.pdb"
    protonated = target_dir / f"{stem}_protonated.pdb"
    source_heavy, selection = write_selected_normalized_structure(
        raw["structure"],
        raw["model"],
        raw["chain"],
        raw["metal_selector"],
        raw["pqq_selector"],
        normalized,
    )
    if pdb_id == "4MAE" and not any(
        item["selector"] == "A:15P603" and item["classification"] == "noncore_heterogen"
        for item in selection["excluded_selected_chain_residues"]
    ):
        raise HoldoutPreparationError("4MAE normalization did not explicitly exclude A:15P603")
    normalization_manifest = target_dir / f"{stem}_normalization_manifest.json"
    normalization = {
        "schema_version": HOLDOUT_NORMALIZATION_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "holdout_id": pdb_id,
        "source": {"path": str(raw["source"]), "sha256": row["source_sha256"]},
        "selector_authority": holdout_pins["holdout_spec"],
        "preparer": holdout_pins["holdout_implementation"]["prepare_holdouts"],
        "selection": selection,
        "source_site_audit": {
            "source_metal": raw["source_metal"],
            "pqq_selector": raw["pqq_selector"].label(),
            "pqq_schema": raw["pqq"]["atom_naming_schema"],
            "typed_donors": raw["donor_records"],
            "source_water_inventory": raw["water_record"],
            "source_altloc_residues": raw["source_altlocs"],
            "excluded_noncore_direct_ligands": raw["noncore_contacts"],
        },
        "output": file_record(normalized),
        "dry_policy": row["dry_policy"],
    }
    write_json_atomic(normalization_manifest, normalization)

    protonation_manifest = standard_only.protonate_standard_only(
        normalized,
        protonated,
        pins_path=calibration_pins_path,
        ph=7.0,
    )
    validate_protonation_manifest(
        protonation_manifest, normalized, protonated, calibration_pins
    )
    normalized_heavy = heavy_map_from_file(normalized)
    protonated_heavy = heavy_map_from_file(protonated)
    checks = [
        compare_heavy_maps(
            source_heavy, normalized_heavy, "retained_raw_source_to_selected_normalized_PDB"
        ),
        compare_heavy_maps(
            normalized_heavy, protonated_heavy, "selected_normalized_to_protonated_PDB"
        ),
        compare_heavy_maps(
            source_heavy, protonated_heavy, "retained_raw_source_to_protonated_PDB"
        ),
    ]
    heavy_check_path = target_dir / f"{stem}_heavy_coordinate_check.json"
    heavy_check = {
        "schema_version": HEAVY_CHECK_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "holdout_id": pdb_id,
        "raw_source": {"path": str(raw["source"]), "sha256": row["source_sha256"]},
        "normalized": file_record(normalized),
        "protonated": file_record(protonated),
        "checks": checks,
        "all_retained_source_heavy_atoms_preserved": True,
        "no_retained_heavy_atoms_added_or_removed": True,
        "excluded_atoms_are_limited_to_unselected_models_chains_altlocs_waters_and_noncore_heterogens": True,
        "passes": True,
    }
    write_json_atomic(heavy_check_path, heavy_check)

    carve_manifest = write_holdout_pair(
        protonated,
        target_dir,
        stem,
        row,
        raw,
        normalization_manifest,
        protonation_manifest,
        heavy_check_path,
        gate,
        holdout_pins_path,
        holdout_pins,
        calibration_pins_path,
        calibration_pins,
    )
    manifest = read_object(carve_manifest)
    return {
        "holdout_index": index,
        "pdb_id": pdb_id,
        "holdout_role": row["holdout_role"],
        "biological_class": row["biological_class"],
        "raw_source": {"path": str(raw["source"]), "sha256": row["source_sha256"]},
        "normalized": file_record(normalized),
        "normalization_manifest": file_record(normalization_manifest),
        "protonated": file_record(protonated),
        "protonation_manifest": file_record(protonation_manifest),
        "heavy_coordinate_check": file_record(heavy_check_path),
        "manifest": file_record(carve_manifest),
        "coordination_number": manifest["coordination"]["coordination_number"],
        "charges": manifest["charge_ledger"]["expected_total_charges"],
        "nonmetal_coordinates_byte_identical": True,
        "orca_executed": False,
    }


def prepare_holdouts(
    *,
    calibration_result: Path,
    output: Path,
    holdout_pins_path: Path,
    include_secondary: bool,
) -> Path:
    holdout_pins_path = holdout_pins_path.resolve()
    holdout_pins, calibration_pins, spec_path = verify_holdout_pins(holdout_pins_path)
    calibration_pins_path = Path(
        holdout_pins["calibration_implementation_pins"]["path"]
    ).resolve()
    if file_record(calibration_result) != holdout_pins.get("calibration_result"):
        raise HoldoutPreparationError(
            "supplied calibration result differs from the locked holdout-pins record"
        )
    gate = verify_calibration_gate(
        calibration_result, calibration_pins_path, calibration_pins
    )
    rows = load_holdout_rows(spec_path)
    selected_rows = [row for row in rows if row["holdout_role"] == "primary"]
    if include_secondary:
        selected_rows.extend(row for row in rows if row["holdout_role"] == "secondary")
    if [row["pdb_id"] for row in selected_rows] not in (
        ["1H4I", "4MAE"],
        ["1H4I", "4MAE", "6OC6"],
    ):
        raise HoldoutPreparationError("holdout target selection changed")

    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    incomplete = output / "HOLDOUT_PREPARATION_INCOMPLETE"
    write_text_new(
        incomplete,
        "This directory is incomplete and must not be submitted to ORCA.\n",
    )
    targets: list[dict[str, Any]] = []
    for index, row in enumerate(selected_rows, start=1):
        targets.append(
            prepare_one(
                row,
                index,
                output,
                gate,
                holdout_pins_path,
                holdout_pins,
                calibration_pins_path,
                calibration_pins,
            )
        )
    preparation = {
        "schema_version": HOLDOUT_PREPARATION_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "ready_for_orca_after_passing_calibration_gate",
        "release_scope": (
            "atomic_primary_pair_plus_nonindependent_secondary"
            if include_secondary
            else "atomic_primary_pair"
        ),
        "release_runnable": True,
        "calibration_gate": gate,
        "holdout_pins": file_record(holdout_pins_path),
        "calibration_implementation_pins": file_record(calibration_pins_path),
        "selector_authority": holdout_pins["holdout_spec"],
        "selector_audit": holdout_pins["holdout_selector_audit"],
        "primary_holdouts": ["1H4I", "4MAE"],
        "secondary_holdout": "6OC6",
        "secondary_included": include_secondary,
        "secondary_omission_disposition": (
            None if include_secondary else "secondary-not-run"
        ),
        "target_count": len(targets),
        "task_count": 2 * len(targets),
        "all_retained_source_heavy_coordinates_preserved": True,
        "all_nonmetal_arm_coordinates_byte_identical": True,
        "water_policy": "dry_exclude_all_source_and_synthetic_waters",
        "noncore_heterogen_policy": "exclude_all_without_replacement_including_4MAE_15P",
        "orca_executed": False,
        "targets": targets,
    }
    preparation_path = output / "holdout_preparation.json"
    write_json_atomic(preparation_path, preparation)
    incomplete.unlink()
    return preparation_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--calibration-result",
        type=Path,
        required=True,
        help="immutable result.json from a fully passing 25-member v3 calibration",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="fresh output directory; existing paths are refused",
    )
    parser.add_argument(
        "--holdout-pins",
        type=Path,
        default=HERE / "holdout_implementation_pins.json",
    )
    parser.add_argument(
        "--include-secondary-6oc6",
        action="store_true",
        help="also prepare nonindependent 6OC6 as a secondary geometry check",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = prepare_holdouts(
            calibration_result=args.calibration_result,
            output=args.output,
            holdout_pins_path=args.holdout_pins,
            include_secondary=args.include_secondary_6oc6,
        )
    except (HoldoutPreparationError, fixed.FixedCoreError, base.PQQCarveError, OSError, ValueError) as exc:
        output = args.output.resolve()
        if output.is_dir() and (output / "HOLDOUT_PREPARATION_INCOMPLETE").is_file():
            write_json_atomic(
                output / "HOLDOUT_PREPARATION_ERROR.json",
                {
                    "schema_version": "alchemical_bvs.pqq_fixed_core_holdout_preparation_error.v1",
                    "protocol_id": PROTOCOL_ID,
                    "error_type": type(exc).__name__,
                    "message": str(exc),
                    "valid_for_orca": False,
                },
            )
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
