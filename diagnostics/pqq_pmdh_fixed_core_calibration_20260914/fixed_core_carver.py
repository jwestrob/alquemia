#!/usr/bin/env python3
"""Prepare a dry, chemistry-defined PQQ-MDH La/Ca vertical pair.

This experiment-local carver deliberately does not alter the production
carver.  It reuses the production implementation's PQQ microstate, residue
fragment, alternate-conformer, metal-selection, and electron-parity helpers,
but replaces its distance-defined QM boundary with an explicit homologous
PQQ-MDH core supplied in the committed ``core_map.tsv``.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Any, Mapping, Sequence

import gemmi


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SCRIPTS = PROJECT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import carve_with_pqq as base  # noqa: E402
from coordination_policy import (  # noqa: E402
    ALTLOC_POLICY_ID,
    COORDINATION_CUTOFF_A,
    MAX_DIRECT_N,
    WATER_NAMES,
    selected_residue_atoms,
)
from pqq_microstates import DEFAULT_PQQ_MICROSTATE  # noqa: E402


PROTOCOL_ID = "pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3"
METHOD_ID = "orca_r2scan3c_cpcm_water_native_basis_v2"
CORE_POLICY_ID = "pqq_pmdh_E_N_D_Dplus2acidic_Dnetwork_cation_v1"
COORDINATION_POLICY_ID = "typed_cn6_3p1A_maxN2_calibration_only_v1"
MANIFEST_SCHEMA = "alchemical_bvs.carve_manifest.v1"
MIN_CALIBRATION_CN = 6
MAX_ASP_PARTNER_DISTANCE_A = 3.5

ROLE_ORDER = (
    "anchor_glutamate",
    "anchor_asparagine",
    "catalytic_aspartate",
    "extra_acidic_ligand_homolog",
    "catalytic_asp_cationic_partner",
)
REQUIRED_ROLES = frozenset(ROLE_ORDER)
ROLE_ALLOWED_RESNAMES = {
    "anchor_glutamate": frozenset({"GLU"}),
    "anchor_asparagine": frozenset({"ASN"}),
    "catalytic_aspartate": frozenset({"ASP"}),
    "extra_acidic_ligand_homolog": frozenset({"ASP", "GLU", "ALA", "SER", "THR"}),
    "catalytic_asp_cationic_partner": frozenset({"ARG", "LYS"}),
}
EXPECTED_CHARGE = {
    "anchor_glutamate": -1,
    "anchor_asparagine": 0,
    "catalytic_aspartate": -1,
    "extra_acidic_ligand_homolog": -1,
    "catalytic_asp_cationic_partner": 1,
}
CATION_HEAVY_ATOMS = {
    "ARG": frozenset({"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}),
    "LYS": frozenset({"CB", "CG", "CD", "CE", "NZ"}),
}
CATION_HYDROGEN_PARENT_COUNTS = {
    "ARG": {"CB": 2, "CG": 2, "CD": 2, "NE": 1, "NH1": 2, "NH2": 2},
    "LYS": {"CB": 2, "CG": 2, "CD": 2, "CE": 2, "NZ": 3},
}
CATION_HBOND_ATOMS = {
    "ARG": frozenset({"NE", "NH1", "NH2"}),
    "LYS": frozenset({"NZ"}),
}


class FixedCoreError(RuntimeError):
    """The frozen calibration core cannot be represented exactly."""


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
        raise FixedCoreError(f"cannot read JSON object {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise FixedCoreError(f"JSON root is not an object: {path}")
    return value


def write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def verify_file_record(record: Mapping[str, Any], label: str) -> Path:
    raw_path = record.get("path")
    expected = record.get("sha256")
    if not isinstance(raw_path, str) or not isinstance(expected, str):
        raise FixedCoreError(f"{label} lacks path/SHA-256")
    path = Path(raw_path).resolve()
    if not path.is_file():
        raise FixedCoreError(f"{label} is missing: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise FixedCoreError(
            f"{label} SHA-256 mismatch: expected {expected}, observed {observed}"
        )
    return path


def verify_pins(pins_path: Path) -> dict[str, Any]:
    pins = read_object(pins_path)
    if pins.get("schema_version") != "alchemical_bvs.pqq_fixed_core_pins.v1":
        raise FixedCoreError("unsupported implementation-pins schema")
    if pins.get("protocol_id") != PROTOCOL_ID:
        raise FixedCoreError("implementation pins carry another protocol ID")
    for name in (
        "preregistration",
        "core_map",
        "core_map_audit",
        "reserved_holdout_spec",
        "reserved_holdout_audit",
        "preparation_amendment",
    ):
        record = pins.get(name)
        if not isinstance(record, dict):
            raise FixedCoreError(f"implementation pins lack {name}")
        verify_file_record(record, f"frozen {name}")
    helpers = pins.get("canonical_helpers")
    if not isinstance(helpers, dict):
        raise FixedCoreError("implementation pins lack canonical_helpers")
    for name in (
        "carve_with_pqq",
        "pqq_microstates",
        "coordination_policy",
        "normalize_af3_cif",
        "protonate_cif",
        "run_orca_task_manifest",
        "render_orca_runtime_input",
    ):
        record = helpers.get(name)
        if not isinstance(record, dict):
            raise FixedCoreError(f"implementation pins lack helper {name}")
        verify_file_record(record, f"canonical helper {name}")
    experiment_helpers = pins.get("experiment_helpers")
    if not isinstance(experiment_helpers, dict):
        raise FixedCoreError("implementation pins lack experiment_helpers")
    for name, record in experiment_helpers.items():
        if not isinstance(record, dict):
            raise FixedCoreError(f"implementation pins lack experiment helper {name}")
        verify_file_record(record, f"experiment helper {name}")
    own_record = experiment_helpers.get("fixed_core_carver")
    if not isinstance(own_record, dict) or Path(str(own_record.get("path", ""))).resolve() != Path(__file__).resolve():
        raise FixedCoreError("fixed-core carver is not self-bound in implementation pins")
    aquo = pins.get("aquo_reference")
    if not isinstance(aquo, dict):
        raise FixedCoreError("implementation pins lack aquo_reference")
    aquo_path = verify_file_record(aquo, "pinned aquo reference")
    aquo_payload = read_object(aquo_path)
    if aquo_payload.get("reference_id") != aquo.get("reference_id"):
        raise FixedCoreError("pinned aquo reference ID mismatch")
    observed_delta = aquo_payload.get("delta_E_aquo_hartree")
    expected_delta = aquo.get("delta_E_aquo_hartree")
    if observed_delta != expected_delta:
        raise FixedCoreError("pinned aquo energy difference mismatch")
    subprotocol = pins.get("protonation_subprotocol")
    if (
        not isinstance(subprotocol, dict)
        or subprotocol.get("id")
        != "pdbfixer_standard_only_rng20260914_openmm_cpu_threads1_v1"
        or subprotocol.get("nonstandard_definition_policy_id")
        != "standard_residues_only_no_nonstandard_CCD_hydrogens_v1"
        or subprotocol.get("python_random_seed") != 20260914
        or subprotocol.get("openmm_platform") != "CPU"
        or subprotocol.get("openmm_cpu_threads") != 1
        or subprotocol.get("forcefield") is not None
        or subprotocol.get("allowed_nonstandard_residue_counts")
        != {"LA": 1, "PQQ": 1}
        or subprotocol.get("required_output_nonstandard_hydrogen_counts")
        != {"LA": 0, "PQQ": 0}
    ):
        raise FixedCoreError("protonation subprotocol pins are incomplete or changed")
    return pins


def residue_key(record: Mapping[str, Any], role: str) -> base.ResidueKey:
    chain = record.get("chain")
    resnum = record.get("resnum")
    icode = record.get("icode", "")
    resname = record.get("resname")
    if (
        not isinstance(chain, str)
        or not chain
        or isinstance(resnum, bool)
        or not isinstance(resnum, int)
        or not isinstance(icode, str)
        or not isinstance(resname, str)
        or not resname
    ):
        raise FixedCoreError(f"invalid residue selector for {role}")
    return base.ResidueKey(chain, resnum, icode.strip(), resname.upper().strip())


def resolve_residue(
    model: gemmi.Model, key: base.ResidueKey, role: str
) -> gemmi.Residue:
    matches = [
        residue
        for chain, residue in base._iter_residues(model)
        if base._residue_key(chain, residue) == key
    ]
    if len(matches) != 1:
        raise FixedCoreError(
            f"{role} selector {key.label()} resolved to {len(matches)} residues"
        )
    return matches[0]


def donor_distances(
    residue: gemmi.Residue,
    key: base.ResidueKey,
    metal_pos: gemmi.Position,
) -> list[dict[str, Any]]:
    if key.resname in {"ARG", "LYS"}:
        names = CATION_HBOND_ATOMS[key.resname]
    else:
        names = {
            "ASP": frozenset({"OD1", "OD2"}),
            "GLU": frozenset({"OE1", "OE2"}),
            "ASN": frozenset({"OD1"}),
            "SER": frozenset({"OG"}),
            "THR": frozenset({"OG1"}),
        }.get(key.resname, frozenset())
    records = []
    for atom in selected_residue_atoms(residue):
        atom_name = base._atom_name(atom)
        if atom_name not in names:
            continue
        records.append(
            {
                "atom": atom_name,
                "element": base._element(atom).title(),
                "distance_A": round(metal_pos.dist(atom.pos), 3),
            }
        )
    records.sort(key=lambda item: (item["distance_A"], item["atom"]))
    return records


def atom_record(
    *,
    name: str,
    element: str,
    origin: str,
    xyz: Sequence[float],
    parent_atom: str | None = None,
) -> dict[str, Any]:
    """Return one ordered, coordinate-bearing fragment provenance record."""
    return {
        "name": name,
        "element": element.title(),
        "origin": origin,
        "parent_atom": parent_atom,
        "xyz_A": [round(float(value), 6) for value in xyz],
    }


def pqq_atom_records(
    residue: gemmi.Residue, prepared: Any
) -> list[dict[str, Any]]:
    records = [
        atom_record(
            name=base._atom_name(atom),
            element=base._element(atom),
            origin="source_heavy_atom",
            xyz=(atom.pos.x, atom.pos.y, atom.pos.z),
        )
        for atom in selected_residue_atoms(residue)
        if base._element(atom) != "H"
    ]
    records.extend(
        atom_record(
            name=str(item["name"]),
            element="H",
            origin=str(item["origin"]),
            parent_atom=str(item["parent_atom"]),
            xyz=item["xyz_A"],
        )
        for item in prepared.metadata()["hydrogens"]
    )
    if len(records) != len(prepared.atoms):
        raise FixedCoreError("PQQ atom provenance count differs from prepared PQQ")
    for record, coordinates in zip(records, prepared.atoms, strict=True):
        element, x, y, z = coordinates
        expected_xyz = [round(float(x), 6), round(float(y), 6), round(float(z), 6)]
        if record["element"].upper() != element.upper() or record["xyz_A"] != expected_xyz:
            raise FixedCoreError("PQQ atom provenance order differs from QM coordinates")
    return records


def canonical_sidechain_fragment(
    residue: gemmi.Residue,
    key: base.ResidueKey,
    role: str,
    metal_pos: gemmi.Position,
) -> tuple[list[tuple[str, float, float, float]], dict[str, Any]]:
    heavy_names = base.SIDECHAIN_HEAVY_ATOMS.get(key.resname)
    if heavy_names is None:
        raise FixedCoreError(f"{role} {key.label()} has no canonical fragment model")
    observed: dict[str, gemmi.Atom] = {}
    for atom in selected_residue_atoms(residue):
        if base._element(atom) == "H":
            continue
        name = base._atom_name(atom)
        if name in observed:
            raise FixedCoreError(f"duplicate heavy atom {name} in {role} {key.label()}")
        observed[name] = atom
    missing = sorted(heavy_names - set(observed))
    if missing:
        raise FixedCoreError(f"{role} {key.label()} lacks atoms {missing}")
    attached = base._attached_hydrogens(residue, heavy_names)
    charge, charge_detail = base._sidechain_formal_charge(residue, key, attached)
    expected = EXPECTED_CHARGE[role]
    if charge != expected:
        raise FixedCoreError(
            f"{role} {key.label()} charge is {charge:+d}, expected {expected:+d}"
        )
    atoms = [
        (atom.element.name, atom.pos.x, atom.pos.y, atom.pos.z)
        for atom in selected_residue_atoms(residue)
        if base._element(atom) != "H" and base._atom_name(atom) in heavy_names
    ]
    atoms.extend(
        ("H", atom.pos.x, atom.pos.y, atom.pos.z) for atom, _, _ in attached
    )
    link_h = base._sidechain_link_h(residue, key)
    atoms.append(link_h)
    atom_records = [
        atom_record(
            name=base._atom_name(atom),
            element=base._element(atom),
            origin="source_heavy_atom",
            xyz=(atom.pos.x, atom.pos.y, atom.pos.z),
        )
        for atom in selected_residue_atoms(residue)
        if base._element(atom) != "H" and base._atom_name(atom) in heavy_names
    ]
    atom_records.extend(
        atom_record(
            name=base._atom_name(atom),
            element="H",
            origin="source_protonation_hydrogen",
            parent_atom=parent,
            xyz=(atom.pos.x, atom.pos.y, atom.pos.z),
        )
        for atom, parent, _ in attached
    )
    atom_records.append(
        atom_record(
            name="H_LINK_CB",
            element="H",
            origin="generated_Cbeta_link_cap",
            parent_atom="CB",
            xyz=link_h[1:],
        )
    )
    return atoms, {
        "kind": "fixed_core_protein_sidechain",
        "role": role,
        "id": key.label(),
        "formal_charge": charge,
        "atom_count": len(atoms),
        "atom_records": atom_records,
        "selection": "explicit_homologous_role_independent_of_metal_distance",
        "metal_facing_atom_distances": donor_distances(residue, key, metal_pos),
        **charge_detail,
    }


def cationic_sidechain_fragment(
    residue: gemmi.Residue,
    key: base.ResidueKey,
    metal_pos: gemmi.Position,
) -> tuple[list[tuple[str, float, float, float]], dict[str, Any]]:
    heavy_names = CATION_HEAVY_ATOMS[key.resname]
    selected = selected_residue_atoms(residue)
    observed = {
        base._atom_name(atom): atom
        for atom in selected
        if base._element(atom) != "H"
    }
    missing = sorted(heavy_names - set(observed))
    if missing:
        raise FixedCoreError(f"cationic partner {key.label()} lacks atoms {missing}")
    attached = base._attached_hydrogens(residue, heavy_names)
    parent_counts: dict[str, int] = {}
    for _, parent, _ in attached:
        parent_counts[parent] = parent_counts.get(parent, 0) + 1
    expected_counts = CATION_HYDROGEN_PARENT_COUNTS[key.resname]
    if parent_counts != expected_counts:
        raise FixedCoreError(
            f"cationic partner {key.label()} hydrogen-parent counts are "
            f"{parent_counts}, expected {expected_counts}"
        )
    atoms = [
        (atom.element.name, atom.pos.x, atom.pos.y, atom.pos.z)
        for atom in selected
        if base._element(atom) != "H" and base._atom_name(atom) in heavy_names
    ]
    atoms.extend(
        ("H", atom.pos.x, atom.pos.y, atom.pos.z) for atom, _, _ in attached
    )
    link_h = base._sidechain_link_h(residue, key)
    atoms.append(link_h)
    atom_records = [
        atom_record(
            name=base._atom_name(atom),
            element=base._element(atom),
            origin="source_heavy_atom",
            xyz=(atom.pos.x, atom.pos.y, atom.pos.z),
        )
        for atom in selected
        if base._element(atom) != "H" and base._atom_name(atom) in heavy_names
    ]
    atom_records.extend(
        atom_record(
            name=base._atom_name(atom),
            element="H",
            origin="source_protonation_hydrogen",
            parent_atom=parent,
            xyz=(atom.pos.x, atom.pos.y, atom.pos.z),
        )
        for atom, parent, _ in attached
    )
    atom_records.append(
        atom_record(
            name="H_LINK_CB",
            element="H",
            origin="generated_Cbeta_link_cap",
            parent_atom="CB",
            xyz=link_h[1:],
        )
    )
    return atoms, {
        "kind": "fixed_core_cationic_sidechain",
        "role": "catalytic_asp_cationic_partner",
        "id": key.label(),
        "formal_charge": 1,
        "atom_count": len(atoms),
        "atom_records": atom_records,
        "selection": "explicit_catalytic_Asp_hydrogen_bond_partner",
        "metal_facing_atom_distances": donor_distances(residue, key, metal_pos),
        "attached_source_hydrogen_parents": sorted(
            parent for _, parent, _ in attached
        ),
    }


def asp_partner_geometry(
    catalytic_asp: gemmi.Residue,
    partner: gemmi.Residue,
) -> dict[str, Any]:
    asp_atoms = [
        atom
        for atom in selected_residue_atoms(catalytic_asp)
        if base._atom_name(atom) in {"OD1", "OD2"}
    ]
    partner_names = CATION_HBOND_ATOMS[partner.name.upper().strip()]
    partner_atoms = [
        atom
        for atom in selected_residue_atoms(partner)
        if base._atom_name(atom) in partner_names
    ]
    pairs = sorted(
        (
            asp_atom.pos.dist(partner_atom.pos),
            base._atom_name(asp_atom),
            base._atom_name(partner_atom),
        )
        for asp_atom in asp_atoms
        for partner_atom in partner_atoms
    )
    if not pairs:
        raise FixedCoreError("catalytic Asp/partner geometry has no donor-acceptor pair")
    distance, asp_atom, partner_atom = pairs[0]
    if distance > MAX_ASP_PARTNER_DISTANCE_A:
        raise FixedCoreError(
            f"mapped catalytic Asp partner is {distance:.3f} A away, beyond "
            f"{MAX_ASP_PARTNER_DISTANCE_A:.1f} A"
        )
    return {
        "asp_atom": asp_atom,
        "partner_atom": partner_atom,
        "heavy_atom_distance_A": round(distance, 3),
        "maximum_allowed_A": MAX_ASP_PARTNER_DISTANCE_A,
        "passes": True,
    }


def write_xyz(
    path: Path,
    *,
    stem: str,
    label: str,
    atoms: Sequence[tuple[str, float, float, float]],
    charge: int,
) -> None:
    with path.open("x") as handle:
        handle.write(f"{len(atoms)}\n")
        handle.write(
            f"{stem} {label} (charge={charge}, mult=1; {PROTOCOL_ID})\n"
        )
        for element, x, y, z in atoms:
            handle.write(f"{element:<3s} {x:>14.6f} {y:>14.6f} {z:>14.6f}\n")


def write_orca_input(
    path: Path, *, xyz_name: str, charge: int, stem: str, label: str
) -> None:
    path.write_text(
        f"# {stem} {label}; {PROTOCOL_ID}\n"
        "! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3\n"
        "%maxcore 8000\n"
        f"* xyzfile {charge} 1 {xyz_name}\n"
    )


def contact_metadata(
    contact: base.Contact, included_keys: frozenset[base.ResidueKey]
) -> dict[str, Any]:
    record = contact.metadata()
    record.pop("included_in_qm_at_3p3A", None)
    record["included_in_fixed_qm_core"] = contact.residue in included_keys
    return record


def load_core_map_entry(
    core_map_path: Path, target_id: str, pins: Mapping[str, Any]
) -> dict[str, str]:
    pins_record = pins.get("core_map")
    if not isinstance(pins_record, dict):
        raise FixedCoreError("implementation pins do not freeze core_map.tsv")
    pinned_path = verify_file_record(pins_record, "frozen core map")
    if pinned_path != core_map_path.resolve():
        raise FixedCoreError("requested core map is not the frozen selector authority")
    with core_map_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 25 or len({row.get("id") for row in rows}) != 25:
        raise FixedCoreError("core map is not the frozen 25-member unique panel")
    matches = [row for row in rows if row.get("id") == target_id]
    if len(matches) != 1:
        raise FixedCoreError(
            f"target ID {target_id!r} resolved to {len(matches)} core-map rows"
        )
    return matches[0]


def parse_residue_selector(value: str, role: str) -> dict[str, Any]:
    match = re.fullmatch(r"([^:]+):([A-Za-z0-9_]+?)(-?[0-9]+)", value)
    if match is None:
        raise FixedCoreError(f"invalid {role} selector in core map: {value!r}")
    return {
        "chain": match.group(1),
        "resname": match.group(2).upper(),
        "resnum": int(match.group(3)),
        "icode": "",
    }


def parse_source_metal_selector(value: str) -> dict[str, Any]:
    match = re.fullmatch(r"([^:]+):([A-Za-z0-9_]+?)(-?[0-9]+):([A-Za-z0-9_]+)", value)
    if match is None:
        raise FixedCoreError(f"invalid source metal selector: {value!r}")
    if match.group(2).upper() != "LIG_B":
        raise FixedCoreError("frozen source metal residue is not LIG_B")
    return {
        "chain": match.group(1),
        "resname": "LA",
        "resnum": int(match.group(3)),
        "icode": "",
        "atom": match.group(4).upper(),
    }


def parse_source_pqq_selector(value: str) -> dict[str, Any]:
    selector = parse_residue_selector(value, "source PQQ")
    if selector["resname"] != "LIG_C":
        raise FixedCoreError("frozen source PQQ residue is not LIG_C")
    selector["resname"] = "PQQ"
    return selector


def map_float(row: Mapping[str, str], field: str) -> float | None:
    raw = row.get(field, "")
    if raw == "":
        return None
    try:
        return float(raw)
    except ValueError as exc:
        raise FixedCoreError(f"invalid numeric core-map field {field}={raw!r}") from exc


def assert_distance(observed: float, expected: float | None, label: str) -> None:
    if expected is None or abs(observed - expected) > 0.01:
        raise FixedCoreError(
            f"{label} distance mismatch: observed {observed:.3f}, expected {expected}"
        )


def carve(
    structure_path: Path,
    out_dir: Path,
    stem: str,
    core_map_path: Path,
    target_id: str,
    pins_path: Path,
) -> Path:
    pins = verify_pins(pins_path)
    protonation_subprotocol = pins["protonation_subprotocol"]
    row = load_core_map_entry(core_map_path, target_id, pins)
    source_record = {
        "path": row["source_cif"],
        "sha256": row["source_sha256"],
    }
    verify_file_record(source_record, f"{target_id} source CIF")

    structure_path = structure_path.resolve()
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if not structure_path.is_file():
        raise FixedCoreError(f"protonated structure is missing: {structure_path}")

    protonation_path = out_dir / f"{stem}_protonation_manifest.json"
    protonation = read_object(protonation_path)
    if protonation.get("protocol_id") != "pdbfixer_standard_residue_protonation_v2":
        raise FixedCoreError("protonation manifest has the wrong protocol")
    nonstandard_policy = protonation.get("nonstandard_definition_policy")
    if (
        protonation.get("experiment_protonation_protocol_id")
        != protonation_subprotocol["id"]
        or protonation.get("canonical_protonator")
        != pins["canonical_helpers"]["protonate_cif"]
        or protonation.get("experiment_wrapper")
        != pins["experiment_helpers"]["protonate_standard_only"]
        or protonation.get("ph") != 7.0
        or protonation.get("add_missing_residues") is not False
        or protonation.get("repaired_missing_atom_count") != 0
        or protonation.get("repaired_missing_terminal_atom_count") != 0
        or protonation.get("cofactor_protonation") != "not_assigned"
        or protonation.get("software") != pins["preparation_runtime"]["packages"]
        or not isinstance(nonstandard_policy, dict)
        or nonstandard_policy.get("policy_id")
        != protonation_subprotocol["nonstandard_definition_policy_id"]
        or nonstandard_policy.get("observed_nonstandard_residue_counts")
        != protonation_subprotocol["allowed_nonstandard_residue_counts"]
        or nonstandard_policy.get("output_nonstandard_hydrogen_counts")
        != protonation_subprotocol[
            "required_output_nonstandard_hydrogen_counts"
        ]
        or nonstandard_policy.get("forcefield")
        != protonation_subprotocol["forcefield"]
        or nonstandard_policy.get("python_random_seed")
        != protonation_subprotocol["python_random_seed"]
        or nonstandard_policy.get("openmm_platform")
        != protonation_subprotocol["openmm_platform"]
        or nonstandard_policy.get("openmm_CPU_Threads")
        != protonation_subprotocol["openmm_cpu_threads"]
    ):
        raise FixedCoreError("protonation manifest violates the frozen preparation")
    normalization_path = out_dir / f"{stem}_normalization_manifest.json"
    normalization = read_object(normalization_path)
    if normalization.get("schema_version") != "alchemical_bvs.pqq_fixed_core_normalization.v1":
        raise FixedCoreError("normalization manifest has the wrong schema")
    if normalization.get("normalizer") != pins["canonical_helpers"]["normalize_af3_cif"]:
        raise FixedCoreError("normalization manifest has the wrong implementation pin")
    if normalization.get("source") != source_record:
        raise FixedCoreError("normalization source does not match frozen core_map.tsv")
    normalized_output = normalization.get("output")
    if not isinstance(normalized_output, dict):
        raise FixedCoreError("normalization output provenance is missing")
    verify_file_record(normalized_output, "normalized structure")
    if protonation.get("source") != normalized_output:
        raise FixedCoreError("protonation source does not match normalized structure")
    protonated_output = protonation.get("output")
    if not isinstance(protonated_output, dict):
        raise FixedCoreError("protonation output provenance is missing")
    if Path(str(protonated_output.get("path", ""))).resolve() != structure_path:
        raise FixedCoreError("protonation manifest points to another output")
    if protonated_output.get("sha256") != sha256_file(structure_path):
        raise FixedCoreError("protonated structure SHA-256 mismatch")
    heavy_check_path = out_dir / f"{stem}_heavy_coordinate_check.json"
    heavy_check = read_object(heavy_check_path)
    if (
        heavy_check.get("schema_version") != "alchemical_bvs.heavy_coordinate_check.v1"
        or heavy_check.get("protocol_id") != PROTOCOL_ID
        or heavy_check.get("passes") is not True
        or heavy_check.get("all_source_heavy_atoms_retained") is not True
        or heavy_check.get("no_heavy_atoms_added") is not True
        or heavy_check.get("source") != source_record
        or heavy_check.get("protonated") != protonated_output
    ):
        raise FixedCoreError("heavy-coordinate preservation record failed")

    structure = gemmi.read_structure(str(structure_path))
    structure.setup_entities()
    if len(structure) != 1:
        raise FixedCoreError(f"expected one structural model, found {len(structure)}")
    model = structure[0]

    waters = [
        base._residue_key(chain, residue).label()
        for chain, residue in base._iter_residues(model)
        if residue.name.upper().strip() in WATER_NAMES
    ]
    if waters:
        raise FixedCoreError(
            "dry calibration forbids all source waters: " + ", ".join(waters)
        )

    selector = parse_source_metal_selector(row["metal_selector_source"])
    sites = base.select_metal_sites(
        model,
        site_chain=selector.get("chain"),
        site_resnum=selector.get("resnum"),
        site_icode=selector.get("icode", ""),
        site_resname=selector.get("resname"),
        site_atom=selector.get("atom"),
    )
    site = sites[0]
    metal_pos = site.atom.pos

    prepared_pqq = base._prepare_nearby_pqq(
        model, metal_pos, DEFAULT_PQQ_MICROSTATE
    )
    base._reject_unsupported_nearby_species(model, site, prepared_pqq)
    pqq_selector = parse_source_pqq_selector(row["pqq_selector_source"])
    pqq_key = residue_key(pqq_selector, "PQQ")
    if set(prepared_pqq) != {pqq_key}:
        observed = sorted(key.label() for key in prepared_pqq)
        raise FixedCoreError(
            f"expected exactly PQQ {pqq_key.label()}, observed {observed}"
        )
    pqq_residue, prepared = prepared_pqq[pqq_key]
    if any(base._element(atom) == "H" for atom in selected_residue_atoms(pqq_residue)):
        raise FixedCoreError(
            "source/protonation path supplied a PQQ hydrogen; the frozen PQQ "
            "specification must supply all PQQ hydrogens"
        )
    if prepared.microstate.microstate_id != DEFAULT_PQQ_MICROSTATE:
        raise FixedCoreError("PQQ microstate is not the frozen oxidized PQQ(3-) state")
    if prepared.microstate.formal_charge != -3 or prepared.microstate.multiplicity != 1:
        raise FixedCoreError("PQQ charge/multiplicity differs from frozen PQQ(3-)")
    if prepared.formula != "C14H3N2O8":
        raise FixedCoreError(f"prepared PQQ formula is {prepared.formula}, not C14H3N2O8")
    if prepared.schema_id != row["pqq_schema"]:
        raise FixedCoreError("normalized PQQ schema differs from frozen core map")

    core_raw = {
        "anchor_glutamate": parse_residue_selector(row["glu_selector"], "anchor Glu"),
        "anchor_asparagine": parse_residue_selector(row["asn_selector"], "anchor Asn"),
        "catalytic_aspartate": parse_residue_selector(
            row["catalytic_asp_selector"], "catalytic Asp"
        ),
        "extra_acidic_ligand_homolog": parse_residue_selector(
            row["plus2_selector"], "D+2 homolog"
        ),
        "catalytic_asp_cationic_partner": parse_residue_selector(
            row["cationic_partner_selector"], "cationic partner"
        ),
    }
    keys: dict[str, base.ResidueKey] = {}
    residues: dict[str, gemmi.Residue] = {}
    for role in ROLE_ORDER:
        raw = core_raw[role]
        if not isinstance(raw, dict):
            raise FixedCoreError(f"{role} selector is not an object")
        key = residue_key(raw, role)
        if key.resname not in ROLE_ALLOWED_RESNAMES[role]:
            raise FixedCoreError(f"{role} has invalid residue identity {key.resname}")
        keys[role] = key
        residues[role] = resolve_residue(model, key, role)

    catalytic_key = keys["catalytic_aspartate"]
    plus2_key = keys["extra_acidic_ligand_homolog"]
    if plus2_key.chain != catalytic_key.chain or plus2_key.resnum != catalytic_key.resnum + 2:
        raise FixedCoreError("extra-ligand homolog is not catalytic-Asp + 2")
    partner_key = keys["catalytic_asp_cationic_partner"]
    if partner_key.chain != catalytic_key.chain:
        raise FixedCoreError("catalytic Asp and cationic partner are on different chains")
    partner_geometry = asp_partner_geometry(
        residues["catalytic_aspartate"],
        residues["catalytic_asp_cationic_partner"],
    )
    expected_asp_atom = row["asp_o_atom"]
    expected_partner_atom = row["partner_n_atom"]
    if (
        partner_geometry["asp_atom"] != expected_asp_atom
        or partner_geometry["partner_atom"] != expected_partner_atom
    ):
        raise FixedCoreError("closest Asp/partner atom pair differs from frozen core map")
    assert_distance(
        partner_geometry["heavy_atom_distance_A"],
        map_float(row, "asp_o_partner_n_A"),
        "catalytic Asp/partner",
    )

    role_distance_fields = {
        "anchor_glutamate": "glu_metal_A",
        "anchor_asparagine": "asn_metal_A",
        "catalytic_aspartate": "catalytic_asp_metal_A",
        "extra_acidic_ligand_homolog": "plus2_metal_A",
    }
    for role, field in role_distance_fields.items():
        expected = map_float(row, field)
        if expected is None and role == "extra_acidic_ligand_homolog":
            continue
        distances = donor_distances(residues[role], keys[role], metal_pos)
        if not distances:
            raise FixedCoreError(f"mapped {role} has no expected metal-facing atom")
        assert_distance(distances[0]["distance_A"], expected, role)

    qualified, _, nearby_sulfur, nearby_untyped = base._classify_contacts(
        model, site, prepared_pqq
    )
    n_direct_n = sum(contact.element.upper() == "N" for contact in qualified)
    n_direct_o = sum(contact.element.upper() == "O" for contact in qualified)
    if len(qualified) < MIN_CALIBRATION_CN or n_direct_n > MAX_DIRECT_N:
        raise FixedCoreError(
            f"calibration geometry gate failed: CN={len(qualified)}, N={n_direct_n}"
        )

    include_plus2_raw = row["include_plus2_as_ln_specific_asp"].lower()
    if include_plus2_raw not in {"true", "false"}:
        raise FixedCoreError("invalid D+2 inclusion flag in frozen core map")
    acidic_plus2 = include_plus2_raw == "true"
    if acidic_plus2 != (plus2_key.resname in {"ASP", "GLU"}):
        raise FixedCoreError("D+2 inclusion flag disagrees with residue chemistry")
    included_roles = {
        "anchor_glutamate",
        "anchor_asparagine",
        "catalytic_aspartate",
        "catalytic_asp_cationic_partner",
    }
    if acidic_plus2:
        included_roles.add("extra_acidic_ligand_homolog")
    included_keys = frozenset({pqq_key, *(keys[role] for role in included_roles)})
    unexpected_direct = sorted(
        {
            contact.residue.label()
            for contact in qualified
            if contact.residue not in included_keys
        }
    )
    if unexpected_direct:
        raise FixedCoreError(
            "qualifying donor lies outside the frozen core: "
            + ", ".join(unexpected_direct)
        )
    expected_fixed_selectors = [
        row["glu_selector"],
        row["asn_selector"],
        row["catalytic_asp_selector"],
    ]
    if acidic_plus2:
        expected_fixed_selectors.append(row["plus2_selector"])
    expected_fixed_selectors.append(row["cationic_partner_selector"])
    if row["fixed_protein_core_selectors"].split(";") != expected_fixed_selectors:
        raise FixedCoreError("fixed-protein selector list is internally inconsistent")

    pqq_records = pqq_atom_records(pqq_residue, prepared)
    fragments: list[dict[str, Any]] = [
        {
            "kind": "fixed_core_pqq",
            "role": "pqq_cofactor",
            "id": pqq_key.label(),
            "formal_charge": prepared.microstate.formal_charge,
            "atom_count": len(prepared.atoms),
            "atom_records": pqq_records,
            "microstate_id": prepared.microstate.microstate_id,
            "selection": "explicit_complete_cofactor",
        }
    ]
    persistent_atoms: list[tuple[str, float, float, float]] = list(prepared.atoms)
    excluded_homolog: dict[str, Any] | None = None
    for role in ROLE_ORDER:
        if role == "extra_acidic_ligand_homolog" and not acidic_plus2:
            excluded_homolog = {
                "role": role,
                "id": plus2_key.label(),
                "resname": plus2_key.resname,
                "included": False,
                "reason": "nonacidic_homolog_has_no_XoxF_extra_carboxylate",
                "metal_facing_atom_distances": donor_distances(
                    residues[role], plus2_key, metal_pos
                ),
            }
            continue
        if role == "catalytic_asp_cationic_partner":
            atoms, fragment = cationic_sidechain_fragment(
                residues[role], keys[role], metal_pos
            )
        else:
            atoms, fragment = canonical_sidechain_fragment(
                residues[role], keys[role], role, metal_pos
            )
        persistent_atoms.extend(atoms)
        fragments.append(fragment)

    scaffold_charge = sum(fragment["formal_charge"] for fragment in fragments)
    charges = {"La": scaffold_charge + 3, "Ca": scaffold_charge + 2}
    expected_charges = (
        {"La": -2, "Ca": -3} if acidic_plus2 else {"La": -1, "Ca": -2}
    )
    if charges != expected_charges:
        raise FixedCoreError(
            f"derived La/Ca charges {charges} differ from frozen expectation "
            f"{expected_charges}"
        )
    la_atoms = [("La", metal_pos.x, metal_pos.y, metal_pos.z), *persistent_atoms]
    ca_atoms = [("Ca", metal_pos.x, metal_pos.y, metal_pos.z), *persistent_atoms]
    base._validate_singlet("La", la_atoms, charges["La"])
    base._validate_singlet("Ca", ca_atoms, charges["Ca"])

    for suffix in ("_La_qm.xyz", "_Ca_qm.xyz", ".inp", ".out"):
        if any(out_dir.glob(f"*{suffix}")):
            raise FixedCoreError(
                f"refusing to overwrite existing scientific artifacts in {out_dir}"
            )

    output_records: dict[str, dict[str, str]] = {}
    tasks = []
    for label, atoms in (("La", la_atoms), ("Ca", ca_atoms)):
        xyz_path = out_dir / f"{stem}_{label}_qm.xyz"
        input_path = out_dir / f"sp_{stem}_{label}.inp"
        output_path = out_dir / f"sp_{stem}_{label}.out"
        write_xyz(
            xyz_path,
            stem=stem,
            label=label,
            atoms=atoms,
            charge=charges[label],
        )
        write_orca_input(
            input_path,
            xyz_name=xyz_path.name,
            charge=charges[label],
            stem=stem,
            label=label,
        )
        xyz_record = {"path": xyz_path.name, "sha256": sha256_file(xyz_path)}
        input_record = {"path": input_path.name, "sha256": sha256_file(input_path)}
        output_records[f"{label}_xyz"] = xyz_record
        output_records[f"{label}_input"] = input_record
        tasks.append(
            {
                "task_id": label,
                "input": input_record,
                "xyz": xyz_record,
                "output_path": output_path.name,
            }
        )

    helper_pins = pins["canonical_helpers"]
    aquo_pin = pins["aquo_reference"]
    script_path = Path(__file__).resolve()
    manifest: dict[str, Any] = {
        "schema_version": MANIFEST_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "ready_for_orca",
        "stem": stem,
        "panel_id": target_id,
        "carver": {
            "name": script_path.name,
            "path": str(script_path),
            "implementation_sha256": sha256_file(script_path),
            "canonical_helper": helper_pins["carve_with_pqq"],
            "pqq_spec_path": helper_pins["pqq_microstates"]["path"],
            "pqq_spec_sha256": helper_pins["pqq_microstates"]["sha256"],
            "coordination_policy_path": helper_pins["coordination_policy"]["path"],
            "coordination_policy_sha256": helper_pins["coordination_policy"]["sha256"],
        },
        "experiment_provenance": {
            "core_map": {"path": str(core_map_path.resolve()), "sha256": sha256_file(core_map_path)},
            "implementation_pins": {"path": str(pins_path.resolve()), "sha256": sha256_file(pins_path)},
            "selector_authority": "committed_core_map.tsv",
            "biological_class_not_copied_to_carve_manifest": True,
        },
        "source_structure": {
            "path": str(structure_path),
            "sha256": sha256_file(structure_path),
        },
        "source_cif": source_record,
        "normalization_manifest": {
            "path": str(normalization_path.resolve()),
            "sha256": sha256_file(normalization_path),
        },
        "normalized_source_structure": normalized_output,
        "protonation_manifest": {
            "path": str(protonation_path.resolve()),
            "sha256": sha256_file(protonation_path),
            "protocol_id": protonation.get("protocol_id"),
            "experiment_protonation_protocol_id": protonation.get(
                "experiment_protonation_protocol_id"
            ),
            "experiment_wrapper": protonation.get("experiment_wrapper"),
        },
        "heavy_coordinate_check": {
            "path": str(heavy_check_path.resolve()),
            "sha256": sha256_file(heavy_check_path),
            "passes": True,
        },
        "site_selection": {"mode": "exact_frozen_panel_selector", "requested": selector},
        "selected_site": {
            "chain": site.chain,
            "resnum": site.resnum,
            "icode": site.icode.strip(),
            "resname": site.resname,
            "atom_name": site.atom_name,
            "source_element": site.element,
            "xyz_A": [round(metal_pos.x, 6), round(metal_pos.y, 6), round(metal_pos.z, 6)],
        },
        "coordination_policy": {
            "policy_id": COORDINATION_POLICY_ID,
            "sha256": helper_pins["coordination_policy"]["sha256"],
            "cutoff_A": COORDINATION_CUTOFF_A,
            "minimum_coordination_number": MIN_CALIBRATION_CN,
            "maximum_direct_nitrogen": MAX_DIRECT_N,
            "forced_core_atoms_count_toward_coordination": False,
            "production_CN7_policy_unchanged": True,
            "alternate_conformer_policy": ALTLOC_POLICY_ID,
        },
        "coordination": {
            "coordination_number": len(qualified),
            "oxygen_count": n_direct_o,
            "direct_nitrogen_count": n_direct_n,
            "passes_geometry": True,
            "source_typed_direct_donors_ordered_by_distance": [
                contact_metadata(item, included_keys) for item in qualified
            ],
            "nearby_sulfur": [item.metadata() for item in nearby_sulfur],
            "nearby_untyped_ON": [item.metadata() for item in nearby_untyped],
        },
        "fixed_core": {
            "policy_id": CORE_POLICY_ID,
            "selection_is_distance_independent": True,
            "requested_roles": core_raw,
            "included_roles": sorted(included_roles),
            "included_residue_keys": sorted(key.label() for key in included_keys),
            "excluded_nonacidic_Dplus2_homolog": excluded_homolog,
            "catalytic_Asp_partner_geometry": partner_geometry,
            "unexpected_direct_donors": [],
            "water_policy": "dry_exclude_all_source_and_synthetic_waters",
            "point_charge_embedding": False,
        },
        "qm_fragments": fragments,
        "pqq": {
            "present": True,
            "microstate_id": prepared.microstate.microstate_id,
            "spec_version": prepared.metadata()["spec_version"],
            "formal_charge": prepared.microstate.formal_charge,
            "multiplicity": prepared.microstate.multiplicity,
            "residues": [
                {
                    **prepared.metadata(),
                    "residue": {
                        "chain": pqq_key.chain,
                        "resnum": pqq_key.resnum,
                        "icode": pqq_key.icode,
                        "resname": pqq_key.resname,
                    },
                }
            ],
        },
        "charge_ledger": {
            "fragments": fragments,
            "scaffold_formal_charge": scaffold_charge,
            "expected_total_charges": expected_charges,
            "La_total": charges["La"],
            "Ca_total": charges["Ca"],
            "multiplicity": 1,
        },
        "electronic_structure": {
            "method_id": METHOD_ID,
            "method": "r2SCAN-3c",
            "basis_policy": "native_r2scan3c_def2_mTZVPP_with_native_ECP",
            "solvation": "CPCM(Water)",
            "grid": "DefGrid3",
            "aquo_reference_id": aquo_pin["reference_id"],
            "aquo_reference_status": "experiment_local_hash_pinned",
            "aquo_reference": aquo_pin,
            "score_formula": "((E_Ca_site-E_La_site)-delta_E_aquo)*627.509474",
        },
        "execution_policy": {
            "id": "two_leg_allocation_aware_manifested_orca_v1",
            "runtime_renderer": helper_pins["render_orca_runtime_input"],
            "task_runner": helper_pins["run_orca_task_manifest"],
            "heavy_legs": ["La", "Ca"],
            "apo_leg": False,
            "water_leg": False,
            "omp_threads_per_rank": 1,
        },
        "tasks": tasks,
        "outputs": output_records,
        "assumptions": [
            "vertical La/Ca swap at identical nuclear coordinates",
            "fixed homologous PQQ-MDH core independent of metal distance",
            "conserved catalytic Asp and its +1 Arg/Lys partner are retained together",
            "D+2 homolog is retained only when it supplies an acidic side chain",
            "no source water, synthetic water, point charges, or geometry relaxation",
        ],
    }
    manifest_path = out_dir / f"{stem}_carve_manifest.json"
    if manifest_path.exists():
        raise FixedCoreError(f"refusing to overwrite manifest {manifest_path}")
    write_json_atomic(manifest_path, manifest)
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--stem", required=True)
    parser.add_argument("--core-map", type=Path, required=True)
    parser.add_argument("--target-id", required=True)
    parser.add_argument("--implementation-pins", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest = carve(
            args.structure,
            args.out_dir,
            args.stem,
            args.core_map,
            args.target_id,
            args.implementation_pins,
        )
    except (FixedCoreError, base.PQQCarveError, OSError, ValueError) as exc:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        error_path = args.out_dir / f"{args.stem}_CARVE_ERROR.json"
        write_json_atomic(
            error_path,
            {
                "schema_version": "alchemical_bvs.pqq_fixed_core_carve_error.v1",
                "protocol_id": PROTOCOL_ID,
                "panel_id": args.target_id,
                "error_type": type(exc).__name__,
                "message": str(exc),
            },
        )
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
