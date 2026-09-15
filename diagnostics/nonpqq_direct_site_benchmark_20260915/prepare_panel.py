#!/usr/bin/env python3
"""Prepare and verify the frozen non-PQQ direct-site La/Ca panel.

This script performs no ORCA calculation and submits no scheduler job.  It
validates the preregistered source hashes, protonates each crystallographic
structure once, makes exact-selector vertical La/Ca carves, and fails closed if
the donor set, source heavy-atom coordinates, or paired nonmetal coordinates
depart from the frozen manifest.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import gemmi


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SCRIPTS = PROJECT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from carve_generic import carve  # noqa: E402
from protonate_cif import protonate  # noqa: E402


MANIFEST_PATH = HERE / "panel_manifest.json"
PREPARED = HERE / "prepared"
REPORT_PATH = HERE / "preparation_report.json"
TASK_INVENTORY_PATH = HERE / "prepared_task_inventory.json"
COORDINATE_TOLERANCE_A = 0.002


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atom_key(entry: dict[str, Any]) -> tuple[str, str, int, str]:
    return (
        str(entry["chain"]),
        str(entry["resname"]).upper(),
        int(entry["resnum"]),
        str(entry["atom"]).upper(),
    )


def structure_atoms(path: Path) -> dict[tuple[str, str, int, str], gemmi.Atom]:
    structure = gemmi.read_structure(str(path))
    if not structure:
        raise RuntimeError(f"no model in {path}")
    atoms: dict[tuple[str, str, int, str], gemmi.Atom] = {}
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                key = (
                    chain.name,
                    residue.name.strip().upper(),
                    residue.seqid.num,
                    atom.name.strip().upper(),
                )
                if key in atoms:
                    raise RuntimeError(f"ambiguous atom selector {key} in {path}")
                atoms[key] = atom
    return atoms


def distance(a: gemmi.Atom, b: gemmi.Atom) -> float:
    return a.pos.dist(b.pos)


def verify_source_coordinates(
    source: Path,
    protonated: Path,
    sites: list[dict[str, Any]],
) -> dict[str, Any]:
    source_atoms = structure_atoms(source)
    protonated_atoms = structure_atoms(protonated)
    checked: list[dict[str, Any]] = []
    maximum = 0.0
    selectors: list[dict[str, Any]] = []
    for site in sites:
        selectors.append(site["metal_selector"])
        selectors.extend(site["expected_donors"])
    seen: set[tuple[str, str, int, str]] = set()
    for selector in selectors:
        key = atom_key(selector)
        if key in seen:
            continue
        seen.add(key)
        if key not in source_atoms:
            raise RuntimeError(f"frozen selector absent from source {source}: {key}")
        if key not in protonated_atoms:
            raise RuntimeError(
                f"frozen selector absent after protonation {protonated}: {key}"
            )
        displacement = distance(source_atoms[key], protonated_atoms[key])
        maximum = max(maximum, displacement)
        checked.append({"selector": list(key), "displacement_A": displacement})
    if maximum > COORDINATE_TOLERANCE_A:
        raise RuntimeError(
            f"source heavy-atom displacement {maximum:.6f} A exceeds "
            f"{COORDINATE_TOLERANCE_A:.3f} A"
        )
    return {
        "status": "pass",
        "coordinate_tolerance_A": COORDINATE_TOLERANCE_A,
        "checked_atom_count": len(checked),
        "maximum_displacement_A": maximum,
        "atoms": checked,
    }


def donor_key(entry: dict[str, Any]) -> tuple[str, str, int, str]:
    return (
        str(entry["chain"]),
        str(entry["resname"]).upper(),
        int(entry["resnum"]),
        str(entry["atom"]).upper(),
    )


def verify_donors(
    frozen_site: dict[str, Any], carve_manifest: dict[str, Any]
) -> dict[str, Any]:
    expected = {donor_key(item): item for item in frozen_site["expected_donors"]}
    observed_entries = carve_manifest["coordination"]["typed_donors"]
    observed = {donor_key(item): item for item in observed_entries}
    if set(expected) != set(observed):
        raise RuntimeError(
            f"donor mismatch for {frozen_site['site_key']}: "
            f"missing={sorted(set(expected) - set(observed))}, "
            f"extra={sorted(set(observed) - set(expected))}"
        )
    distance_deltas: dict[str, float] = {}
    for key in sorted(expected):
        delta = abs(float(expected[key]["distance_A"]) - float(observed[key]["distance_A"]))
        distance_deltas["/".join(map(str, key))] = delta
        if delta > COORDINATE_TOLERANCE_A:
            raise RuntimeError(
                f"donor distance changed by {delta:.6f} A for {key}"
            )
    expected_waters = {
        donor_key(item) for item in frozen_site["first_shell_waters"]
    }
    observed_waters = {
        donor_key(item)
        for item in observed_entries
        if str(item["resname"]).upper() in {"HOH", "WAT"}
    }
    if expected_waters != observed_waters:
        raise RuntimeError(
            f"first-shell water mismatch for {frozen_site['site_key']}: "
            f"expected={sorted(expected_waters)}, observed={sorted(observed_waters)}"
        )
    observed_cn = int(carve_manifest["coordination"]["coordination_number"])
    if observed_cn != int(frozen_site["expected_typed_cn"]):
        raise RuntimeError(
            f"CN mismatch for {frozen_site['site_key']}: "
            f"expected {frozen_site['expected_typed_cn']}, observed {observed_cn}"
        )
    return {
        "status": "pass",
        "coordination_number": observed_cn,
        "donor_keys": [list(key) for key in sorted(observed)],
        "first_shell_water_keys": [list(key) for key in sorted(observed_waters)],
        "maximum_distance_delta_A": max(distance_deltas.values(), default=0.0),
    }


def xyz_records(path: Path) -> tuple[str, list[str]]:
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise RuntimeError(f"truncated XYZ: {path}")
    atom_count = int(lines[0])
    records = lines[2:]
    if len(records) != atom_count:
        raise RuntimeError(
            f"XYZ atom count mismatch in {path}: header={atom_count}, records={len(records)}"
        )
    return lines[1], records


def verify_vertical_pair(la_xyz: Path, ca_xyz: Path) -> dict[str, Any]:
    _, la_records = xyz_records(la_xyz)
    _, ca_records = xyz_records(ca_xyz)
    if len(la_records) != len(ca_records):
        raise RuntimeError("La/Ca XYZ atom-count mismatch")
    la_metal = la_records[0].split()
    ca_metal = ca_records[0].split()
    if la_metal[0].upper() != "LA" or ca_metal[0].upper() != "CA":
        raise RuntimeError("vertical pair does not begin with La/Ca")
    if la_metal[1:] != ca_metal[1:]:
        raise RuntimeError("La and Ca coordinates differ")
    if la_records[1:] != ca_records[1:]:
        raise RuntimeError("La/Ca nonmetal atom records are not byte-identical")
    nonmetal_payload = ("\n".join(la_records[1:]) + "\n").encode()
    return {
        "status": "pass",
        "atom_count": len(la_records),
        "nonmetal_atom_count": len(la_records) - 1,
        "metal_xyz_tokens": la_metal[1:],
        "nonmetal_records_sha256": hashlib.sha256(nonmetal_payload).hexdigest(),
        "La_xyz_sha256": sha256_file(la_xyz),
        "Ca_xyz_sha256": sha256_file(ca_xyz),
    }


def validate_static_hashes(panel: dict[str, Any]) -> dict[str, Any]:
    checked: list[dict[str, str]] = []
    for protein in panel["proteins"]:
        path = Path(protein["source"]["path"])
        observed = sha256_file(path)
        expected = protein["source"]["sha256"]
        if observed != expected:
            raise RuntimeError(f"source hash mismatch: {path}")
        checked.append({"path": str(path), "sha256": observed})
    for item in panel["software_sources"].values():
        path = Path(item["path"])
        observed = sha256_file(path)
        if observed != item["sha256"]:
            raise RuntimeError(f"software/reference hash mismatch: {path}")
        checked.append({"path": str(path), "sha256": observed})
    return {"status": "pass", "files": checked}


def main() -> int:
    panel = json.loads(MANIFEST_PATH.read_text())
    if REPORT_PATH.exists() or TASK_INVENTORY_PATH.exists():
        raise RuntimeError(
            "frozen preparation report already exists; refusing to overwrite it"
        )
    static_hashes = validate_static_hashes(panel)
    PREPARED.mkdir(parents=True, exist_ok=True)

    protein_reports: list[dict[str, Any]] = []
    tasks: list[dict[str, Any]] = []
    global_order_observed: list[str] = []
    for protein in panel["proteins"]:
        protein_dir = PREPARED / protein["panel_key"]
        protein_dir.mkdir(exist_ok=True)
        source = Path(protein["source"]["path"])
        protonated = protein_dir / f"{protein['panel_key']}_protonated.pdb"
        protonation_manifest_path = (
            protein_dir / f"{protein['panel_key']}_protonation_manifest.json"
        )
        if protonated.exists() and protonation_manifest_path.exists():
            protonation = json.loads(protonation_manifest_path.read_text())
            if protonation.get("source", {}).get("sha256") != protein["source"]["sha256"]:
                raise RuntimeError(
                    f"existing protonation source mismatch for {protein['panel_key']}"
                )
            if protonation.get("output", {}).get("sha256") != sha256_file(protonated):
                raise RuntimeError(
                    f"existing protonated output hash mismatch for {protein['panel_key']}"
                )
        elif protonated.exists() or protonation_manifest_path.exists():
            raise RuntimeError(
                f"partial protonation output for {protein['panel_key']}; refusing overwrite"
            )
        else:
            protonation = protonate(
                source,
                protonated,
                ph=float(panel["protocol"]["ph"]),
                add_missing_residues=bool(panel["protocol"]["add_missing_residues"]),
            )
        coordinate_check = verify_source_coordinates(
            source, protonated, protein["sites"]
        )
        site_reports: list[dict[str, Any]] = []
        for site in protein["sites"]:
            selector = site["metal_selector"]
            site_dir = protein_dir / site["site_key"]
            stem = f"{protein['panel_key']}_{site['site_key']}"
            carve_manifest_path = site_dir / f"{stem}_carve_manifest.json"
            if carve_manifest_path.exists():
                required = [
                    site_dir / f"{stem}_La_qm.xyz",
                    site_dir / f"{stem}_Ca_qm.xyz",
                    site_dir / f"sp_{stem}_La.inp",
                    site_dir / f"sp_{stem}_Ca.inp",
                ]
                if not all(path.exists() for path in required):
                    raise RuntimeError(
                        f"partial carve output for {protein['panel_key']}/{site['site_key']}"
                    )
            elif site_dir.exists():
                raise RuntimeError(
                    f"unmanifested carve directory for {protein['panel_key']}/{site['site_key']}"
                )
            else:
                carve(
                    protonated,
                    site_dir,
                    stem,
                    site_chain=selector["chain"],
                    site_resnum=int(selector["resnum"]),
                    site_icode=selector["icode"],
                    site_atom=selector["atom"],
                    qm_inclusion_cut=float(panel["protocol"]["qm_inclusion_cutoff_A"]),
                )
            carve_manifest = json.loads(carve_manifest_path.read_text())
            if carve_manifest.get("status") not in {"prepared", "ready_for_orca"}:
                raise RuntimeError(
                    f"carve ineligible for {protein['panel_key']}/{site['site_key']}: "
                    f"{carve_manifest.get('status')} {carve_manifest.get('reason')}"
                )
            donor_check = verify_donors(site, carve_manifest)
            la_xyz = site_dir / f"{stem}_La_qm.xyz"
            ca_xyz = site_dir / f"{stem}_Ca_qm.xyz"
            pair_check = verify_vertical_pair(la_xyz, ca_xyz)
            ordered_key = f"{protein['panel_key']}/{site['site_key']}"
            global_order_observed.append(ordered_key)
            site_reports.append(
                {
                    "site_key": site["site_key"],
                    "ordered_key": ordered_key,
                    "carve_manifest": {
                        "path": str(carve_manifest_path.resolve()),
                        "sha256": sha256_file(carve_manifest_path),
                    },
                    "donor_check": donor_check,
                    "vertical_pair_check": pair_check,
                }
            )
            for metal in ("La", "Ca"):
                input_path = site_dir / f"sp_{stem}_{metal}.inp"
                tasks.append(
                    {
                        "task_key": f"{ordered_key}/{metal}",
                        "panel_key": protein["panel_key"],
                        "site_key": site["site_key"],
                        "metal": metal,
                        "input_path": str(input_path.resolve()),
                        "input_sha256": sha256_file(input_path),
                        "expected_output_path": str(
                            input_path.with_suffix(".out").resolve()
                        ),
                        "submitted": False,
                    }
                )
        protein_reports.append(
            {
                "panel_key": protein["panel_key"],
                "source_path": str(source),
                "source_sha256": sha256_file(source),
                "protonation_manifest": protonation,
                "source_coordinate_check": coordinate_check,
                "sites": site_reports,
            }
        )

    if global_order_observed != panel["global_site_order"]:
        raise RuntimeError(
            f"global site order mismatch: observed={global_order_observed}, "
            f"frozen={panel['global_site_order']}"
        )

    task_inventory = {
        "schema_version": "nonpqq_direct_site_task_inventory.v1",
        "panel_id": panel["panel_id"],
        "status": "prepared_not_submitted",
        "task_count": len(tasks),
        "task_order": [task["task_key"] for task in tasks],
        "tasks": tasks,
    }
    TASK_INVENTORY_PATH.write_text(
        json.dumps(task_inventory, indent=2, sort_keys=True) + "\n"
    )
    report = {
        "schema_version": "nonpqq_direct_site_preparation_report.v1",
        "panel_id": panel["panel_id"],
        "status": "prepared_not_submitted",
        "panel_manifest": {
            "path": str(MANIFEST_PATH.resolve()),
            "sha256": sha256_file(MANIFEST_PATH),
        },
        "preregistration": {
            "path": str((HERE / "PREREGISTRATION.md").resolve()),
            "sha256": sha256_file(HERE / "PREREGISTRATION.md"),
        },
        "static_hash_validation": static_hashes,
        "site_order": global_order_observed,
        "site_count": len(global_order_observed),
        "orca_task_count": len(tasks),
        "orca_tasks_submitted": 0,
        "task_inventory": {
            "path": str(TASK_INVENTORY_PATH.resolve()),
            "sha256": sha256_file(TASK_INVENTORY_PATH),
        },
        "proteins": protein_reports,
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"prepared {len(global_order_observed)} sites and {len(tasks)} La/Ca inputs")
    print(f"report: {REPORT_PATH}")
    print("submitted: 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
