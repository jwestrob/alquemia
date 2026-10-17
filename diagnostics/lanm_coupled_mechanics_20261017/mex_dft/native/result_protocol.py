#!/usr/bin/env python3
"""Shared result-schema, provenance, and protocol-selection utilities.

The JSONL is an append-only event log.  Its current value for a stem is the
last decodable event for that stem -- not the first positive result and not
the event with the lexicographically largest timestamp.  Both result writers
and downstream trackers import this module so that rule cannot drift.

This module deliberately does *not* identify a historical carver by hashing
the carver currently present in ``scripts/``.  A carver hash is authoritative
only when it was captured in the per-workspace carve manifest at carve time.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable


RESULT_SCHEMA_VERSION = "alchemical_bvs.result.v4"
CARVE_MANIFEST_SUFFIX = "_carve_manifest.json"
CAPTURE_DESCRIPTOR_ID = "vertical_aquo_to_frozen_site_conditioned_panel_v2"
ORCA_EXECUTION_SCHEMA = "alchemical_bvs.orca_task_execution.v1"

# These are the references actually used by the historical production score.
# Despite old prose mentioning CN6, both coordinate files contain M + 8 waters
# and the ORCA inputs are single points (there is no Opt keyword).
LEGACY_AQUO_REFERENCE_ID = (
    "aquo_cn8_independent_handbuilt_sp_mixed_r2scan3c_cpcm_v1"
)
LEGACY_AQUO_REFERENCE: dict[str, Any] = {
    "reference_id": LEGACY_AQUO_REFERENCE_ID,
    "status": "complete",
    "coordination_number": 8,
    "geometry_protocol": "independent_handbuilt_Ca_and_La_CN8_single_points",
    "E_Ca_Ha": -1288.921041427522,
    "E_La_Ha": -642.856485875777,
    "delta_E_aquo_hartree": -646.064555551745,
    "method_id": "orca_r2scan3c_cpcm_water_mixed_la_def2tzvp_v1",
    "source": {
        "Ca_xyz": "legacy/qmmm/aquo_panel/Ca/aquo_Ca_cn8.xyz",
        "Ca_xyz_sha256": (
            "9cd7d490eaa3eb0835b62fd305e9b627f8c453103847abad3625bd2870b8d787"
        ),
        "La_xyz": "legacy/qmmm/aquo_test/aquo_La_cn8.xyz",
        "La_xyz_sha256": (
            "b783ac6f73a3c3f6dfdd4c1567fda6e623cf26b99b4cdbbd7f66016bebdaf9e8"
        ),
    },
}

EXPECTED_AQUO_COORDINATION_NUMBER = 8
NATIVE_AQUO_REFERENCE_RELATIVE = Path(
    "reference_inputs/aquo_cn8_native_r2scan3c_v2/aquo_reference.json"
)


class ProtocolSelectionError(ValueError):
    """Raised when a comparable table would silently pool protocols."""


def infer_qm_stem(qm_dir: Path) -> str:
    """Infer the calculation stem, including nested ``--all-sites`` outputs."""
    manifests = sorted(qm_dir.glob(f"*{CARVE_MANIFEST_SUFFIX}"))
    if len(manifests) == 1:
        return manifests[0].name[: -len(CARVE_MANIFEST_SUFFIX)]
    if len(manifests) > 1:
        raise ValueError(
            f"ambiguous QM directory {qm_dir}: multiple carve manifests "
            f"{[path.name for path in manifests]}"
        )
    return qm_dir.name[:-3] if qm_dir.name.endswith("_qm") else qm_dir.name


def iter_workspace_qm_dirs(workspace: Path) -> list[Path]:
    """Return ordinary QM workspaces plus one-level all-site subdirectories."""
    direct = [path for path in workspace.glob("workspaces/*_qm") if path.is_dir()]
    nested = [
        site_dir
        for parent in direct
        for site_dir in parent.glob("site_*")
        if site_dir.is_dir()
    ]
    return sorted({*direct, *nested}, key=lambda path: str(path))


def sha256_file(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError:
        return None
    return digest.hexdigest()


def _stable_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def load_carve_manifest(qm_dir: Path, stem: str) -> tuple[dict[str, Any] | None, Path]:
    path = qm_dir / f"{stem}{CARVE_MANIFEST_SUFFIX}"
    return _read_json(path), path


def _normalise_simple_input(lines: Iterable[str]) -> list[str]:
    values: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not stripped.startswith("!"):
            continue
        tokens = [token.lower() for token in stripped[1:].split()]
        # NoAutostart changes restart behaviour, not the model chemistry.
        tokens = [token for token in tokens if token != "noautostart"]
        values.extend(tokens)
    return sorted(set(values))


def _normalise_basis_block(lines: list[str]) -> list[str]:
    block: list[str] = []
    in_basis = False
    for line in lines:
        stripped = " ".join(line.strip().lower().split())
        if stripped.startswith("%basis"):
            in_basis = True
        if in_basis and stripped:
            block.append(stripped)
        # ORCA accepts both a multiline block and an inline spelling such as
        # ``%basis NewGTO La \"def2-TZVP\" end``.  Treat an ``end`` token on
        # the opening line as the end of that block; otherwise an inline basis
        # declaration would absorb every following percent block.
        if in_basis and (
            stripped == "end"
            or (
                stripped.startswith("%basis")
                and "end" in stripped.split()[1:]
            )
        ):
            break
    return block


def _normalise_physical_blocks(lines: list[str]) -> list[str]:
    """Canonicalize non-resource ORCA blocks that can change the Hamiltonian."""
    excluded = {"%pal", "%maxcore", "%basis"}
    blocks: list[str] = []
    index = 0
    while index < len(lines):
        stripped = " ".join(lines[index].strip().lower().split())
        if not stripped.startswith("%"):
            index += 1
            continue
        keyword = stripped.split()[0]
        if keyword in excluded:
            if keyword in {"%pal", "%basis"}:
                # Resource/basis blocks may be written wholly on their opening
                # line (for example ``%pal nprocs 4 end``).  In that case the
                # block is already closed; scanning for a later bare ``end``
                # would swallow a subsequent physical block such as ``%cpcm``.
                if "end" in stripped.split()[1:]:
                    index += 1
                    continue
                index += 1
                while index < len(lines):
                    if " ".join(lines[index].strip().lower().split()) == "end":
                        index += 1
                        break
                    index += 1
            else:
                index += 1
            continue
        block = [stripped]
        if len(stripped.split()) == 1:
            index += 1
            while index < len(lines):
                value = " ".join(lines[index].strip().lower().split())
                if value:
                    block.append(value)
                index += 1
                if value == "end":
                    break
        else:
            index += 1
        blocks.extend(block)
    return blocks


def effective_basis_block_for_elements(
    parsed_input: dict[str, Any], elements: Iterable[str]
) -> list[str]:
    """Return only element-specific basis/ECP overrides relevant to a geometry."""
    wanted = {str(element).lower() for element in elements}
    relevant: list[str] = []
    for line in parsed_input.get("basis_block", []):
        match = re.match(r"new(?:gto|ecp)\s+(\S+)\s+", line, re.IGNORECASE)
        if match and match.group(1).strip('"').lower() in wanted:
            relevant.append(line)
    return relevant


_XYZFILE_RE = re.compile(
    r"^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+(\S+)", re.IGNORECASE
)


def parse_orca_input(path: Path) -> dict[str, Any]:
    """Extract model chemistry and charge from an ORCA input.

    The full-file hash proves which concrete input was run.  The recipe hash
    intentionally omits filenames, charges, and resource controls so rows made
    with the same physical model share an electronic-structure method ID.
    """
    if not path.is_file():
        return {
            "path": str(path),
            "sha256": None,
            "simple_input": [],
            "basis_block": [],
            "physical_blocks": [],
            "noautostart": False,
            "charge": None,
            "multiplicity": None,
            "xyzfile": None,
        }
    try:
        lines = path.read_text(errors="ignore").splitlines()
    except OSError:
        lines = []
    charge = multiplicity = None
    xyzfile = None
    raw_simple_tokens = {
        token.lower()
        for line in lines
        if line.strip().startswith("!")
        for token in line.strip()[1:].split()
    }
    for line in lines:
        match = _XYZFILE_RE.match(line)
        if match:
            charge = int(match.group(1))
            multiplicity = int(match.group(2))
            xyzfile = match.group(3)
            break
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "simple_input": _normalise_simple_input(lines),
        "basis_block": _normalise_basis_block(lines),
        "physical_blocks": _normalise_physical_blocks(lines),
        "noautostart": "noautostart" in raw_simple_tokens,
        "charge": charge,
        "multiplicity": multiplicity,
        "xyzfile": xyzfile,
    }


def electronic_structure_metadata(
    la_input: dict[str, Any], ca_input: dict[str, Any]
) -> dict[str, Any]:
    recipe = {
        "La_simple_input": la_input["simple_input"],
        "Ca_simple_input": ca_input["simple_input"],
        "La_basis_block": la_input["basis_block"],
        "Ca_basis_block": ca_input["basis_block"],
        "La_physical_blocks": la_input.get("physical_blocks", []),
        "Ca_physical_blocks": ca_input.get("physical_blocks", []),
    }
    recipe_hash = _stable_hash(recipe)
    la_tokens = set(la_input["simple_input"])
    ca_tokens = set(ca_input["simple_input"])
    la_basis = " ".join(la_input["basis_block"])
    ca_basis = " ".join(ca_input["basis_block"])
    exact_simple = {"r2scan-3c", "cpcm(water)", "defgrid3"}
    exact_shared_recipe = (
        la_tokens == ca_tokens == exact_simple
        and la_input.get("physical_blocks", []) == ca_input.get("physical_blocks", []) == []
    )
    shared_model = (
        exact_shared_recipe
    )
    legacy_mixed_basis = (
        shared_model
        and 'newgto la "def2-tzvp" end' in la_basis
        and 'newecp la "def2-ecp" end' in la_basis
        and not ca_basis
    )
    native_r2scan3c = shared_model and not la_basis and not ca_basis
    if legacy_mixed_basis:
        method_id = LEGACY_AQUO_REFERENCE["method_id"]
        basis_policy = "La_def2-ECP_def2-TZVP_override_Ca_native_composite"
    elif native_r2scan3c:
        method_id = "orca_r2scan3c_cpcm_water_native_basis_v2"
        basis_policy = "native_r2SCAN-3c_basis_ECP_policy"
    else:
        method_id = f"orca_recipe_sha256_{recipe_hash[:16]}"
        basis_policy = "unclassified"
    return {
        "method_id": method_id,
        "recipe_sha256": recipe_hash,
        "basis_policy": basis_policy,
        "La_simple_input": la_input["simple_input"],
        "Ca_simple_input": ca_input["simple_input"],
        "La_basis_block": la_input["basis_block"],
        "Ca_basis_block": ca_input["basis_block"],
        "La_physical_blocks": la_input.get("physical_blocks", []),
        "Ca_physical_blocks": ca_input.get("physical_blocks", []),
        "La_input_sha256": la_input["sha256"],
        "Ca_input_sha256": ca_input["sha256"],
    }


_PROGRAM_VERSION_RE = re.compile(r"Program Version\s+([0-9][0-9A-Za-z._-]*)")


def parse_orca_version(path: Path) -> str | None:
    if not path.is_file():
        return None
    try:
        for line in path.read_text(errors="ignore").splitlines():
            match = _PROGRAM_VERSION_RE.search(line)
            if match:
                return match.group(1)
    except OSError:
        return None
    return None


def orca_execution_record_path(output_path: Path) -> Path:
    return output_path.with_name(output_path.name + ".execution.json")


def verify_orca_execution_record(
    *,
    output_path: Path,
    manifest_path: Path,
    task_id: str,
    input_path: Path,
    xyz_path: Path,
    engrad_path: Path | None = None,
) -> list[str]:
    """Verify that a completed output was run from the manifested artifacts."""
    if not output_path.is_file():
        return []
    prefix = task_id
    record_path = orca_execution_record_path(output_path)
    record = _read_json(record_path)
    if record is None:
        return [f"{prefix}_execution_record_missing_or_invalid"]
    reasons: list[str] = []
    if record.get("schema_version") != ORCA_EXECUTION_SCHEMA:
        reasons.append(f"{prefix}_execution_schema_mismatch")
    if record.get("task_id") != task_id:
        reasons.append(f"{prefix}_execution_task_id_mismatch")
    manifest_record = record.get("manifest")
    if not isinstance(manifest_record, dict) or manifest_record.get(
        "sha256"
    ) != sha256_file(manifest_path):
        reasons.append(f"{prefix}_execution_manifest_sha256_mismatch")
    if record.get("returncode") != 0:
        reasons.append(f"{prefix}_execution_returncode_nonzero")
    if record.get("normal_termination") is not True:
        reasons.append(f"{prefix}_execution_normal_termination_false")
    if record.get("scf_converged") is not True:
        reasons.append(f"{prefix}_execution_SCF_converged_false")
    observed_version = parse_orca_version(output_path)
    if not observed_version or record.get("orca_version") != observed_version:
        reasons.append(f"{prefix}_execution_ORCA_version_mismatch")
    executable = record.get("orca_executable")
    if not isinstance(executable, dict) or not all(
        isinstance(executable.get(key), str) and executable.get(key)
        for key in ("path", "sha256")
    ):
        reasons.append(f"{prefix}_orca_executable_provenance_incomplete")

    source_manifest = _read_json(manifest_path)
    if source_manifest is not None:
        execution_policy = source_manifest.get("execution_policy")
        implementation = source_manifest.get("implementation")
        declared_runner = (
            execution_policy.get("task_runner")
            if isinstance(execution_policy, dict)
            else None
        )
        if not isinstance(declared_runner, dict) and isinstance(implementation, dict):
            declared_runner = implementation.get("execution_runner")
        declared_renderer = (
            execution_policy.get("runtime_renderer")
            if isinstance(execution_policy, dict)
            else None
        )

        def identity_matches(observed: Any, declared: Any) -> bool:
            if not isinstance(declared, dict):
                return True
            if not isinstance(observed, dict):
                return False
            declared_path = declared.get("path")
            observed_path = observed.get("path")
            declared_hash = declared.get("sha256")
            observed_hash = observed.get("sha256")
            if not all(
                isinstance(value, str) and value
                for value in (
                    declared_path,
                    observed_path,
                    declared_hash,
                    observed_hash,
                )
            ):
                return False
            return (
                Path(declared_path).resolve() == Path(observed_path).resolve()
                and declared_hash == observed_hash
            )

        if not identity_matches(record.get("task_runner"), declared_runner):
            reasons.append(f"{prefix}_execution_task_runner_mismatch")
        if not identity_matches(record.get("runtime_renderer"), declared_renderer):
            reasons.append(f"{prefix}_execution_runtime_renderer_mismatch")

    artifacts = record.get("artifacts")
    if not isinstance(artifacts, dict):
        return [*reasons, f"{prefix}_execution_artifacts_missing"]
    expected: dict[str, tuple[Path, str | None]] = {
        "template_input": (input_path, sha256_file(input_path)),
        "xyz": (xyz_path, sha256_file(xyz_path)),
        "output": (output_path, sha256_file(output_path)),
    }
    if engrad_path is not None:
        expected["engrad"] = (engrad_path, sha256_file(engrad_path))
    for name, (path, digest) in expected.items():
        artifact = artifacts.get(name)
        if not isinstance(artifact, dict):
            reasons.append(f"{prefix}_execution_{name}_missing")
            continue
        path_value = artifact.get("path")
        if not isinstance(path_value, str) or Path(path_value).resolve() != path.resolve():
            reasons.append(f"{prefix}_execution_{name}_path_mismatch")
        if digest is None or artifact.get("sha256") != digest:
            reasons.append(f"{prefix}_execution_{name}_sha256_mismatch")
    for name in ("runtime_input", "runtime_sidecar"):
        artifact = artifacts.get(name)
        if not isinstance(artifact, dict):
            reasons.append(f"{prefix}_execution_{name}_missing")
            continue
        path_value = artifact.get("path")
        digest = artifact.get("sha256")
        if not isinstance(path_value, str) or not isinstance(digest, str):
            reasons.append(f"{prefix}_execution_{name}_provenance_incomplete")
        elif sha256_file(Path(path_value)) != digest:
            reasons.append(f"{prefix}_execution_{name}_sha256_mismatch")
    return reasons


def parse_xyz_atoms(
    path: Path,
) -> tuple[int | None, list[tuple[str, float, float, float]] | None]:
    """Return a strict element/coordinate list from an XYZ file."""
    if not path.is_file():
        return None, None
    try:
        lines = path.read_text(errors="ignore").splitlines()
        n_atoms = int(lines[0].strip())
        atom_lines = [line for line in lines[2:] if line.strip()]
        if len(atom_lines) != n_atoms:
            return n_atoms, None
        atoms: list[tuple[str, float, float, float]] = []
        for line in atom_lines:
            fields = line.split()
            if len(fields) < 4:
                return n_atoms, None
            coordinates = [float(fields[index]) for index in (1, 2, 3)]
            if not all(math.isfinite(value) for value in coordinates):
                return n_atoms, None
            atoms.append(
                (fields[0].upper(), *(round(value, 7) for value in coordinates))
            )
        return n_atoms, atoms
    except (OSError, ValueError):
        return None, None


def parse_xyz_signature(path: Path) -> tuple[int | None, list[Any] | None]:
    """Return atom count and a metal-agnostic coordinate signature."""
    n_atoms, atoms = parse_xyz_atoms(path)
    if atoms is None:
        return n_atoms, None
    signature: list[Any] = []
    for index, (element, x, y, z) in enumerate(atoms):
        signature.append(["M" if index == 0 else element, x, y, z])
    return n_atoms, signature


def _manifest_retained_water_count(
    manifest: dict[str, Any] | None,
) -> tuple[int | None, list[str]]:
    """Count chemically complete holo-only waters declared by a v2 carver."""
    if manifest is None:
        return None, []
    fragments = manifest.get("qm_fragments")
    if not isinstance(fragments, list):
        return None, ["capture_manifest_qm_fragments_missing"]

    count = 0
    reasons: list[str] = []
    for index, fragment in enumerate(fragments):
        if not isinstance(fragment, dict):
            reasons.append(f"capture_manifest_fragment_{index}_invalid")
            continue
        fragment_type = fragment.get("fragment_type")
        kind = fragment.get("kind")
        if fragment_type != "water" and kind != "inner_shell_water":
            continue
        count += 1
        if fragment.get("atom_count") != 3:
            reasons.append("capture_manifest_water_not_three_atoms")
        if fragment.get("formal_charge") != 0:
            reasons.append("capture_manifest_water_not_neutral")
        included_in = fragment.get("included_in")
        if included_in is not None and (
            not isinstance(included_in, list)
            or set(included_in) != {"La", "Ca"}
        ):
            reasons.append("capture_manifest_water_leg_membership_invalid")

    coordination = manifest.get("coordination")
    if isinstance(coordination, dict) and isinstance(coordination.get("waters"), list):
        if len(coordination["waters"]) != count:
            reasons.append("capture_manifest_water_lists_disagree")
    return count, reasons


def _manifest_fragment_atom_counts(
    manifest: dict[str, Any] | None,
) -> tuple[dict[str, int] | None, list[str]]:
    if manifest is None or not isinstance(manifest.get("qm_fragments"), list):
        return None, []
    totals = {"La": 1, "Ca": 1, "apo": 0}
    reasons: list[str] = []
    for index, fragment in enumerate(manifest["qm_fragments"]):
        if not isinstance(fragment, dict):
            reasons.append(f"capture_manifest_fragment_{index}_invalid")
            continue
        atom_count = fragment.get("atom_count")
        if isinstance(atom_count, bool) or not isinstance(atom_count, int) or atom_count < 1:
            reasons.append(f"capture_manifest_fragment_{index}_atom_count_invalid")
            continue
        kind = fragment.get("kind")
        fragment_type = fragment.get("fragment_type")
        default_legs = (
            {"La", "Ca"}
            if kind == "inner_shell_water" or fragment_type == "water"
            else {"La", "Ca", "apo"}
        )
        included = fragment.get("included_in")
        legs = set(included) if isinstance(included, list) else default_legs
        if not legs <= {"La", "Ca", "apo"}:
            reasons.append(f"capture_manifest_fragment_{index}_leg_membership_invalid")
            continue
        for leg in legs:
            totals[leg] += atom_count
    return totals, reasons


def capture_water_accounting(
    *,
    manifest: dict[str, Any] | None,
    la_atoms: list[tuple[str, float, float, float]] | None,
    ca_atoms: list[tuple[str, float, float, float]] | None,
    apo_atoms: list[tuple[str, float, float, float]] | None,
    bulk_water_atoms: list[tuple[str, float, float, float]] | None,
) -> tuple[dict[str, Any], list[str]]:
    """Prove that holo = apo + metal + an integer number of neutral waters.

    Coordinates, not merely atom counts, are compared.  A versioned manifest
    supplies an independent fragment-level water count; historical unmanifested
    rows can only receive the explicitly labelled geometry-only accounting.
    """
    reasons: list[str] = []
    manifest_count, manifest_reasons = _manifest_retained_water_count(manifest)
    reasons.extend(manifest_reasons)
    fragment_counts, fragment_count_reasons = _manifest_fragment_atom_counts(manifest)
    reasons.extend(fragment_count_reasons)
    geometry_count: int | None = None

    if la_atoms is None:
        reasons.append("capture_La_geometry_missing")
    if ca_atoms is None:
        reasons.append("capture_Ca_geometry_missing")
    if apo_atoms is None:
        reasons.append("capture_apo_geometry_missing")
    if bulk_water_atoms is None:
        reasons.append("capture_bulk_water_geometry_missing")

    if bulk_water_atoms is not None:
        elements = Counter(atom[0] for atom in bulk_water_atoms)
        if len(bulk_water_atoms) != 3 or elements != Counter({"O": 1, "H": 2}):
            reasons.append("capture_bulk_water_not_one_complete_H2O")

    if la_atoms is not None and ca_atoms is not None:
        if not la_atoms or la_atoms[0][0] != "LA":
            reasons.append("capture_La_geometry_first_atom_not_La")
        if not ca_atoms or ca_atoms[0][0] != "CA":
            reasons.append("capture_Ca_geometry_first_atom_not_Ca")
        if la_atoms[1:] != ca_atoms[1:]:
            reasons.append("capture_La_Ca_nonmetal_geometries_differ")

    if la_atoms is not None and la_atoms and apo_atoms is not None:
        holo_counter = Counter(la_atoms[1:])
        apo_counter = Counter(apo_atoms)
        if any(apo_counter[atom] > holo_counter[atom] for atom in apo_counter):
            reasons.append("capture_apo_atoms_not_subset_of_holo")
        else:
            released_atoms = holo_counter - apo_counter
            released_elements: Counter[str] = Counter()
            for (element, *_), count in released_atoms.items():
                released_elements[element] += count
            n_o = released_elements.get("O", 0)
            n_h = released_elements.get("H", 0)
            if set(released_elements) - {"O", "H"}:
                reasons.append("capture_holo_minus_apo_contains_nonwater_elements")
            elif n_h != 2 * n_o or sum(released_elements.values()) != 3 * n_o:
                reasons.append("capture_holo_minus_apo_not_integer_H2O")
            else:
                geometry_count = n_o

    if (
        manifest_count is not None
        and geometry_count is not None
        and manifest_count != geometry_count
    ):
        reasons.append("capture_manifest_geometry_water_count_mismatch")

    if manifest is not None:
        ledger = manifest.get("charge_ledger")
        clusters = ledger.get("clusters") if isinstance(ledger, dict) else None
        actual_counts = {
            "La": len(la_atoms) if la_atoms is not None else None,
            "Ca": len(ca_atoms) if ca_atoms is not None else None,
            "apo": len(apo_atoms) if apo_atoms is not None else None,
        }
        if isinstance(clusters, dict):
            for label, actual in actual_counts.items():
                cluster = clusters.get(label) or clusters.get(label.lower())
                declared = cluster.get("atom_count") if isinstance(cluster, dict) else None
                if (
                    isinstance(declared, int)
                    and not isinstance(declared, bool)
                    and actual is not None
                    and declared != actual
                ):
                    reasons.append(f"capture_manifest_{label}_atom_count_mismatch")
        if fragment_counts is not None:
            for label, actual in actual_counts.items():
                if actual is not None and fragment_counts[label] != actual:
                    reasons.append(
                        f"capture_manifest_{label}_fragment_atom_count_mismatch"
                    )

    accounting_source = (
        "geometry_and_manifest"
        if manifest_count is not None
        else "geometry_only_legacy_unversioned"
    )
    retained = geometry_count if not reasons else None
    return {
        "retained_inner_water_count": retained,
        "retained_inner_water_count_geometry": geometry_count,
        "retained_inner_water_count_manifest": manifest_count,
        "water_accounting_source": accounting_source,
    }, sorted(set(reasons))


def infer_fold_selection(stem: str) -> dict[str, Any]:
    sample_index = None
    sample_source = None
    for pattern, source in (
        (r"_s(\d+)$", "stem_suffix_s"),
        (r"_sample_(\d+)(?:$|_)", "stem_sample"),
    ):
        match = re.search(pattern, stem, re.IGNORECASE)
        if match:
            sample_index = int(match.group(1))
            sample_source = source
            break
    rank = None
    rank_match = re.search(r"_rank_?(\d+)(?:$|_)", stem, re.IGNORECASE)
    if rank_match:
        rank = int(rank_match.group(1))
    return {
        "fold_sample_index": sample_index,
        "fold_sample_index_source": sample_source,
        "structure_rank": rank,
        "structure_rank_source": "stem_rank" if rank is not None else None,
    }


def _manifest_get_hash(value: Any) -> str | None:
    if isinstance(value, dict):
        for key in ("implementation_sha256", "sha256", "module_sha256"):
            candidate = value.get(key)
            if isinstance(candidate, str) and candidate:
                return candidate
    return None


def _manifest_get_name(value: Any) -> str | None:
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        for key in ("name", "path"):
            candidate = value.get(key)
            if isinstance(candidate, str) and candidate:
                return candidate
    return None


def _contains_pqq(qm_dir: Path) -> bool:
    for path in sorted(qm_dir.glob("*_protonated.pdb"))[:2] + sorted(
        qm_dir.glob("*_normalized.pdb")
    )[:2]:
        try:
            for line in path.read_text(errors="ignore").splitlines():
                if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "PQQ":
                    return True
        except OSError:
            continue
    return False


def protocol_metadata(
    qm_dir: Path,
    stem: str,
    electronic_structure: dict[str, Any],
    manifest: dict[str, Any] | None,
    manifest_path: Path,
) -> dict[str, Any]:
    if manifest is not None:
        protocol_id = manifest.get("protocol_id")
        carver = manifest.get("carver")
        policy = manifest.get("coordination_policy")
        carver_sha = _manifest_get_hash(carver)
        policy_sha = _manifest_get_hash(policy)
        if policy_sha is None and isinstance(carver, dict):
            candidate = carver.get("coordination_policy_sha256")
            policy_sha = candidate if isinstance(candidate, str) and candidate else None
        selected_site = manifest.get("selected_site")
        charge_ledger = manifest.get("charge_ledger")
        declared_electronic = manifest.get("electronic_structure")
        pqq_microstate = manifest.get("pqq_microstate")
        if isinstance(pqq_microstate, dict):
            pqq_microstate = pqq_microstate.get("microstate_id") or pqq_microstate.get("id")
        if pqq_microstate is None and isinstance(manifest.get("pqq"), dict):
            pqq_microstate = manifest["pqq"].get("microstate_id")
        source_structure = manifest.get("source_structure")
        if not isinstance(source_structure, dict):
            source_structure = manifest.get("structure")
        protonation_manifest = manifest.get("protonation_manifest")
        coordination = manifest.get("coordination")
        if not isinstance(coordination, dict):
            coordination = manifest.get("qualification")
        source_structure_complete = (
            isinstance(source_structure, dict)
            and isinstance(source_structure.get("path"), str)
            and bool(source_structure.get("path"))
            and isinstance(source_structure.get("sha256"), str)
            and bool(source_structure.get("sha256"))
        )
        electronic_complete = (
            isinstance(declared_electronic, dict)
            and isinstance(declared_electronic.get("method_id"), str)
            and bool(declared_electronic.get("method_id"))
            and isinstance(declared_electronic.get("basis_policy"), str)
            and bool(declared_electronic.get("basis_policy"))
            and isinstance(declared_electronic.get("aquo_reference_id"), str)
            and bool(declared_electronic.get("aquo_reference_id"))
        )
        required_present = all(
            (
                isinstance(protocol_id, str) and bool(protocol_id),
                bool(carver_sha),
                bool(policy_sha),
                isinstance(selected_site, dict) and bool(selected_site),
                isinstance(charge_ledger, dict) and bool(charge_ledger),
                source_structure_complete,
                electronic_complete,
                isinstance(coordination, dict) and bool(coordination),
            )
        )
        provenance_status = "versioned" if required_present else "versioned_incomplete"
        requested_aquo = (
            declared_electronic.get("aquo_reference_id")
            if isinstance(declared_electronic, dict)
            else None
        )
        aquo_status = (
            declared_electronic.get("aquo_reference_status")
            if isinstance(declared_electronic, dict)
            else None
        )
        values = {
            "protocol_id": protocol_id,
            "protocol_provenance_status": provenance_status,
            "carver_name": _manifest_get_name(carver),
            "carver_version": manifest.get("carver_version")
            or (carver.get("version") if isinstance(carver, dict) else None),
            "carver_sha256": carver_sha,
            "coordination_policy_version": manifest.get("coordination_policy_version")
            or (policy.get("version") if isinstance(policy, dict) else None),
            "coordination_policy_sha256": policy_sha,
            "carve_manifest_schema_version": manifest.get("schema_version"),
            "carve_manifest_path": str(manifest_path),
            "carve_manifest_sha256": sha256_file(manifest_path),
            "selected_site": selected_site,
            "manifest_source_structure": source_structure,
            "manifest_protonation": protonation_manifest,
            "manifest_coordination": coordination,
            "manifest_charge_ledger": charge_ledger,
            "manifest_electronic_structure": declared_electronic,
            "pqq_microstate_id": pqq_microstate,
            "requested_aquo_reference_id": requested_aquo,
            "declared_aquo_reference_status": aquo_status,
        }
    else:
        family = "pqq" if _contains_pqq(qm_dir) else "generic_or_unknown"
        method_id = electronic_structure["method_id"]
        protocol_id = f"legacy_unversioned_{family}_{method_id}"
        values = {
            "protocol_id": protocol_id,
            "protocol_provenance_status": "legacy_unversioned",
            "carver_name": family,
            "carver_version": None,
            "carver_sha256": None,
            "coordination_policy_version": None,
            "coordination_policy_sha256": None,
            "carve_manifest_schema_version": None,
            "carve_manifest_path": None,
            "carve_manifest_sha256": None,
            "selected_site": None,
            "manifest_source_structure": None,
            "manifest_protonation": None,
            "manifest_coordination": None,
            "manifest_charge_ledger": None,
            "manifest_electronic_structure": None,
            "pqq_microstate_id": None,
            "requested_aquo_reference_id": (
                LEGACY_AQUO_REFERENCE_ID
                if method_id == LEGACY_AQUO_REFERENCE["method_id"]
                else None
            ),
            "declared_aquo_reference_status": "legacy_inferred_from_input_recipe",
        }

    comparison_recipe = {
        "protocol_id": values["protocol_id"],
        "electronic_structure_method_id": electronic_structure["method_id"],
        "electronic_structure_recipe_sha256": electronic_structure["recipe_sha256"],
        "carver_sha256": values["carver_sha256"],
        "coordination_policy_sha256": values["coordination_policy_sha256"],
        "pqq_microstate_id": values["pqq_microstate_id"],
    }
    values["comparison_recipe"] = comparison_recipe
    return values


def load_aquo_references(workspace: Path) -> dict[str, dict[str, Any]]:
    references = {LEGACY_AQUO_REFERENCE_ID: dict(LEGACY_AQUO_REFERENCE)}
    # The corrected native reference has one canonical location.  Additional
    # references may be added beneath reference_inputs without changing this
    # reader, but malformed/ambiguous records are never accepted implicitly.
    candidates = [workspace / NATIVE_AQUO_REFERENCE_RELATIVE]
    reference_root = workspace / "reference_inputs"
    if reference_root.is_dir():
        candidates.extend(sorted(reference_root.glob("*/aquo_reference.json")))
    for path in dict.fromkeys(candidates):
        record = _read_json(path)
        if record is None:
            continue
        ref_id = record.get("reference_id") or record.get("id")
        if not isinstance(ref_id, str) or not ref_id:
            continue
        record = dict(record)
        record["registry_path"] = str(path)
        record["registry_sha256"] = sha256_file(path)
        if ref_id in references:
            existing = references[ref_id]
            conflict_paths = existing.get("registry_conflict_paths")
            if not isinstance(conflict_paths, list):
                conflict_paths = [existing.get("registry_path", "built_in_legacy_record")]
            references[ref_id] = {
                "reference_id": ref_id,
                "status": "registry_conflict",
                "registry_conflict_paths": [*conflict_paths, str(path)],
            }
        else:
            references[ref_id] = record
    return references


def aquo_reference_registry_fingerprint(
    references: dict[str, dict[str, Any]],
) -> str:
    """Hash the reference inventory so blocked rows become stale only on change."""
    inventory = []
    for reference_id, record in sorted(references.items()):
        inventory.append(
            {
                "reference_id": reference_id,
                "status": record.get("status"),
                "coordination_number": _reference_cn(record),
                "method_id": record.get("method_id"),
                "delta_E_aquo_hartree": record.get("delta_E_aquo_hartree"),
                # Individual aquo energies cancel from the exchange score but
                # are required by the two capture legs.  They must therefore
                # participate in freshness even when their difference is fixed.
                "E_Ca_Ha": record.get("E_Ca_Ha"),
                "E_La_Ha": record.get("E_La_Ha"),
                "energies_hartree": record.get("energies_hartree"),
                "compatible_protocol_ids": record.get("compatible_protocol_ids"),
                "registry_path": record.get("registry_path"),
                "registry_sha256": record.get("registry_sha256"),
                "live_provenance_status": _capture_reference_provenance_status(
                    record
                ),
            }
        )
    return _stable_hash(inventory)


def _reference_cn(reference: dict[str, Any]) -> int | None:
    for key in ("coordination_number", "water_coordination_number", "cn"):
        value = reference.get(key)
        if isinstance(value, int) and not isinstance(value, bool):
            return value
        if isinstance(value, str) and value.isdigit():
            return int(value)
    return None


def finite_real(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    converted = float(value)
    return converted if math.isfinite(converted) else None


def _capture_reference_provenance_status(reference: dict[str, Any]) -> str:
    reference_id = reference.get("reference_id") or reference.get("id")
    if reference_id == LEGACY_AQUO_REFERENCE_ID:
        return "compatible"
    if reference.get("schema_version") != 1:
        return "aquo_finalized_schema_missing"
    registry_path_value = reference.get("registry_path")
    registry_hash = reference.get("registry_sha256")
    if not isinstance(registry_path_value, str) or not isinstance(registry_hash, str):
        return "aquo_registry_provenance_missing"
    registry_path = Path(registry_path_value)
    if sha256_file(registry_path) != registry_hash:
        return "aquo_registry_sha256_mismatch"
    reference_dir = registry_path.parent
    manifest_name = reference.get("reference_manifest")
    manifest_hash = reference.get("reference_manifest_sha256")
    if not isinstance(manifest_name, str) or not isinstance(manifest_hash, str):
        return "aquo_manifest_provenance_missing"
    if sha256_file(reference_dir / manifest_name) != manifest_hash:
        return "aquo_manifest_sha256_mismatch"
    manifest_path = reference_dir / manifest_name
    source_manifest = _read_json(manifest_path)
    if source_manifest is None or source_manifest.get("schema_version") != 1:
        return "aquo_manifest_invalid"
    calculations = reference.get("calculations")
    energies = reference.get("energies_hartree")
    required = {
        "Ca_on_Ca_geometry",
        "Ca_on_La_geometry",
        "La_on_Ca_geometry",
        "La_on_La_geometry",
    }
    if not isinstance(calculations, dict) or set(calculations) != required:
        return "aquo_calculation_provenance_incomplete"
    if not isinstance(energies, dict):
        return "aquo_absolute_energies_missing"
    source_calculations = source_manifest.get("calculations")
    if not isinstance(source_calculations, dict) or set(source_calculations) != required:
        return "aquo_manifest_calculations_incomplete"
    for label in sorted(required):
        calculation = calculations.get(label)
        if not isinstance(calculation, dict):
            return "aquo_calculation_provenance_incomplete"
        source_calculation = source_calculations.get(label)
        if not isinstance(source_calculation, dict):
            return "aquo_manifest_calculations_incomplete"
        for path_key, hash_key in (
            ("input", "input_sha256"),
            ("xyz", "xyz_sha256"),
            ("expected_output", "output_sha256"),
        ):
            path_value = calculation.get(path_key)
            expected_hash = calculation.get(hash_key)
            if not isinstance(path_value, str) or not isinstance(expected_hash, str):
                return "aquo_calculation_provenance_incomplete"
            if source_calculation.get(path_key) != path_value:
                return f"aquo_{label}_{path_key}_manifest_mismatch"
            if path_key != "expected_output" and source_calculation.get(hash_key) != expected_hash:
                return f"aquo_{label}_{path_key}_manifest_hash_mismatch"
            if sha256_file(reference_dir / path_value) != expected_hash:
                return f"aquo_{label}_{path_key}_sha256_mismatch"
        execution_reasons = verify_orca_execution_record(
            output_path=reference_dir / calculation["expected_output"],
            manifest_path=manifest_path,
            task_id=label,
            input_path=reference_dir / calculation["input"],
            xyz_path=reference_dir / calculation["xyz"],
        )
        if execution_reasons:
            return f"aquo_{label}_execution_provenance_failed"
        recorded_energy = finite_real(calculation.get("energy_hartree"))
        panel_energy = finite_real(energies.get(label))
        if recorded_energy is None or panel_energy is None or recorded_energy != panel_energy:
            return "aquo_recorded_energy_mismatch"
    return "compatible"


def resolve_capture_reference_energies(
    reference: dict[str, Any] | None,
) -> tuple[dict[str, Any] | None, str]:
    """Resolve per-ion aquo energies without weakening exchange-reference QC.

    The legacy record contains one energy per ion.  The symmetric v2 record
    contains both ions on both conditioning geometries; its capture legs use
    the per-ion arithmetic means so their difference is exactly the symmetric
    exchange gauge.  A reference that stores only the pre-subtracted gap can
    still score exchange but cannot score capture.
    """
    if reference is None:
        return None, "capture_aquo_reference_unavailable"
    if reference.get("status") != "complete":
        return None, "capture_aquo_reference_incomplete"
    provenance_status = _capture_reference_provenance_status(reference)
    if provenance_status != "compatible":
        return None, provenance_status

    e_ca = finite_real(reference.get("E_Ca_Ha"))
    e_la = finite_real(reference.get("E_La_Ha"))
    policy: str | None = None
    source_keys: dict[str, list[str]] | None = None
    if e_ca is not None and e_la is not None:
        policy = "explicit_per_ion_single_reference_geometry"
        source_keys = {"Ca": ["E_Ca_Ha"], "La": ["E_La_Ha"]}
        by_geometry = {
            "Ca": {"single_reference_geometry": e_ca},
            "La": {"single_reference_geometry": e_la},
        }
    else:
        energies = reference.get("energies_hartree")
        if not isinstance(energies, dict):
            return None, "capture_aquo_absolute_energies_missing"
        ca_keys = ["Ca_on_Ca_geometry", "Ca_on_La_geometry"]
        la_keys = ["La_on_Ca_geometry", "La_on_La_geometry"]
        ca_values = [finite_real(energies.get(key)) for key in ca_keys]
        la_values = [finite_real(energies.get(key)) for key in la_keys]
        if any(value is None for value in ca_values + la_values):
            return None, "capture_aquo_symmetric_energy_panel_incomplete"
        e_ca = sum(value for value in ca_values if value is not None) / 2.0
        e_la = sum(value for value in la_values if value is not None) / 2.0
        policy = "per_ion_mean_over_Ca_and_La_conditioned_geometries"
        source_keys = {"Ca": ca_keys, "La": la_keys}
        by_geometry = {
            "Ca": {
                "Ca_conditioned_geometry": ca_values[0],
                "La_conditioned_geometry": ca_values[1],
            },
            "La": {
                "Ca_conditioned_geometry": la_values[0],
                "La_conditioned_geometry": la_values[1],
            },
        }

    delta = finite_real(reference.get("delta_E_aquo_hartree"))
    if delta is None:
        return None, "capture_aquo_delta_missing"
    assert e_ca is not None and e_la is not None and policy is not None
    residual = (e_ca - e_la) - delta
    if abs(residual) > 1.0e-9:
        return None, "capture_aquo_absolute_energy_delta_mismatch"
    return {
        "E_Ca_aquo_Ha": e_ca,
        "E_La_aquo_Ha": e_la,
        "energy_policy": policy,
        "source_energy_keys": source_keys,
        "per_ion_energies_by_geometry_Ha": by_geometry,
        "per_ion_geometry_span_Ha": {
            ion: max(values.values()) - min(values.values())
            for ion, values in by_geometry.items()
        },
        "delta_consistency_residual_Ha": residual,
    }, "compatible"


def resolve_aquo_reference(
    protocol: dict[str, Any],
    electronic_structure: dict[str, Any],
    references: dict[str, dict[str, Any]],
) -> tuple[dict[str, Any] | None, str]:
    protocol_id = protocol.get("protocol_id")
    requested = protocol.get("requested_aquo_reference_id")

    # Historical scores are accepted only when their concrete ORCA inputs
    # identify the exact mixed-basis method paired with the built-in CN8 record.
    if (
        protocol.get("protocol_provenance_status") == "legacy_unversioned"
        and electronic_structure.get("method_id") == LEGACY_AQUO_REFERENCE["method_id"]
    ):
        requested = LEGACY_AQUO_REFERENCE_ID

    candidates: list[dict[str, Any]] = []
    if requested:
        record = references.get(requested)
        if record is None:
            return None, "aquo_reference_not_found"
        candidates = [record]
    else:
        for record in references.values():
            compatible = record.get("compatible_protocol_ids", [])
            if isinstance(compatible, list) and protocol_id in compatible:
                candidates.append(record)
        if not candidates:
            return None, "aquo_reference_required"
        if len(candidates) > 1:
            return None, "aquo_reference_ambiguous"

    reference = candidates[0]
    reference_id = reference.get("reference_id") or reference.get("id")
    if reference.get("status") == "registry_conflict":
        return reference, "aquo_reference_duplicate_id"
    if reference.get("status") != "complete":
        return reference, "aquo_reference_incomplete"
    cn = _reference_cn(reference)
    if cn is None:
        return reference, "aquo_reference_coordination_number_missing"
    if cn != EXPECTED_AQUO_COORDINATION_NUMBER:
        return reference, f"aquo_reference_CN{cn}_incompatible_with_CN8_protocol"
    # Downstream rows use one canonical spelling even if a registry used an
    # accepted legacy alias such as ``cn``.
    reference["coordination_number"] = cn
    delta = finite_real(reference.get("delta_E_aquo_hartree"))
    if delta is None:
        return reference, "aquo_reference_delta_missing"
    reference["delta_E_aquo_hartree"] = delta

    compatible = reference.get("compatible_protocol_ids")
    if protocol.get("protocol_provenance_status") == "versioned":
        if not isinstance(compatible, list) or not compatible:
            return reference, "aquo_reference_compatible_protocols_missing"
        if protocol_id not in compatible:
            return reference, "aquo_reference_protocol_incompatible"
    elif isinstance(compatible, list) and compatible and protocol_id not in compatible:
        return reference, "aquo_reference_protocol_incompatible"
    reference_method = reference.get("method_id")
    if (
        protocol.get("protocol_provenance_status") == "versioned"
        and not isinstance(reference_method, str)
    ):
        return reference, "aquo_reference_method_missing"
    if reference_method and reference_method != electronic_structure.get("method_id"):
        return reference, "aquo_reference_method_incompatible"
    if reference_id != LEGACY_AQUO_REFERENCE_ID:
        provenance_status = _capture_reference_provenance_status(reference)
        if provenance_status != "compatible":
            return reference, provenance_status
    return reference, "compatible"


def make_comparison_group_id(
    protocol: dict[str, Any],
    electronic_structure: dict[str, Any],
    aquo_reference_id: str | None,
) -> str | None:
    protocol_id = protocol.get("protocol_id")
    if not protocol_id or not aquo_reference_id:
        return None
    recipe = dict(protocol["comparison_recipe"])
    recipe["aquo_reference_id"] = aquo_reference_id
    recipe["capture_descriptor_id"] = CAPTURE_DESCRIPTOR_ID
    return f"{protocol_id}__{_stable_hash(recipe)[:12]}"


def normalise_result_fields(row: dict[str, Any]) -> dict[str, Any]:
    """Return a copy with historical aliases reconciled.

    ``n_atoms`` and ``workspace_path`` are canonical in current rows.  The two
    old names remain as read aliases because add_cn_3A.py and
    analyze_bidentate.py still consume historical JSONL records.
    """
    result = dict(row)
    if result.get("n_atoms") is None and result.get("n_atoms_la") is not None:
        result["n_atoms"] = result["n_atoms_la"]
    if result.get("n_atoms_la") is None and result.get("n_atoms") is not None:
        result["n_atoms_la"] = result["n_atoms"]
    if not result.get("workspace_path") and result.get("cluster_panel_path"):
        result["workspace_path"] = result["cluster_panel_path"]
    if not result.get("cluster_panel_path") and result.get("workspace_path"):
        result["cluster_panel_path"] = result["workspace_path"]
    ddE = result.get("ddE_kcal")
    if isinstance(ddE, bool) or (
        isinstance(ddE, (int, float)) and not math.isfinite(float(ddE))
    ):
        result["ddE_kcal"] = None
        result["class"] = "QC_FAIL (nonfinite score)"
        result["qc_status"] = "fail"
        reasons = result.get("qc_reasons")
        reasons = list(reasons) if isinstance(reasons, list) else []
        result["qc_reasons"] = sorted(set([*reasons, "nonfinite_ddE_score"]))
    return result


def load_latest_rows(jsonl: Path) -> dict[str, dict[str, Any]]:
    """Load the last decodable JSONL event for every stem."""
    latest: dict[str, dict[str, Any]] = {}
    if not jsonl.is_file():
        return latest
    with jsonl.open() as handle:
        for line_number, line in enumerate(handle, 1):
            try:
                raw = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not isinstance(raw, dict) or not raw.get("stem"):
                continue
            row = normalise_result_fields(raw)
            row["_jsonl_line_number"] = line_number
            latest[row["stem"]] = row
    return latest


def comparable_protocol_groups(
    rows: Iterable[dict[str, Any]], *, allow_legacy_unversioned: bool = False
) -> dict[str, list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        group = row.get("comparison_group_id")
        status = row.get("protocol_provenance_status")
        if not group:
            continue
        if status == "legacy_unversioned" and not allow_legacy_unversioned:
            continue
        if status not in ("versioned", "legacy_unversioned"):
            continue
        groups.setdefault(group, []).append(row)
    return groups


def select_comparable_rows(
    rows: Iterable[dict[str, Any]],
    *,
    requested_protocol: str | None = None,
    allow_legacy_unversioned: bool = False,
) -> tuple[list[dict[str, Any]], str]:
    """Select exactly one comparison group or raise.

    ``requested_protocol`` may be either the full comparison-group ID (most
    precise) or a declared protocol ID.  A protocol ID that maps to multiple
    carver/method/reference hashes remains ambiguous and is rejected.
    """
    materialised = [normalise_result_fields(row) for row in rows]
    groups = comparable_protocol_groups(
        materialised, allow_legacy_unversioned=allow_legacy_unversioned
    )
    if requested_protocol:
        matching = {
            group_id: group_rows
            for group_id, group_rows in groups.items()
            if group_id == requested_protocol
            or any(row.get("protocol_id") == requested_protocol for row in group_rows)
        }
        if not matching:
            raise ProtocolSelectionError(
                f"No comparable rows found for protocol {requested_protocol!r}."
            )
        if len(matching) != 1:
            raise ProtocolSelectionError(
                f"Protocol {requested_protocol!r} maps to {len(matching)} distinct "
                "method/carver/aquo groups; specify the full comparison_group_id."
            )
        selected_id = next(iter(matching))
    else:
        if not groups:
            legacy_hint = (
                " Pass --allow-legacy-unversioned and an explicit --protocol-id "
                "to inspect historical rows."
                if not allow_legacy_unversioned
                else ""
            )
            raise ProtocolSelectionError(
                "No versioned, aquo-calibrated comparison group is available."
                + legacy_hint
            )
        if len(groups) != 1:
            listing = ", ".join(sorted(groups))
            raise ProtocolSelectionError(
                f"Refusing to mix {len(groups)} incompatible comparison groups: "
                f"{listing}. Select one with --protocol-id."
            )
        selected_id = next(iter(groups))

    return groups[selected_id], selected_id
