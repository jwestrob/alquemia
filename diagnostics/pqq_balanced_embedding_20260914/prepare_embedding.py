#!/usr/bin/env python3
"""Prepare and verify the preregistered balanced PQQ embedding panel."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import openmm
from openmm import app, unit


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
PREPARED = HERE / "prepared"
TASK_MANIFEST = HERE / "embedded_tasks.json"
PREPARATION = HERE / "preparation.json"

FORCEFIELD = Path(
    "/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/lib/python3.11/"
    "site-packages/openmm/app/data/amber19/protein.ff19SB.xml"
)
FORCEFIELD_SHA256 = "086bfdb1c05e7d5cb1d330e6c39fab5056a95e89492cdbd3044c746db13e9d09"
RUNNER = PROJECT / "scripts/run_orca_task_manifest.py"
RUNNER_SHA256 = "b9278e74e5317e859e8cb05541e1ec57cb335b26dcdf65a83b89675407343b9a"
RENDERER = PROJECT / "scripts/render_orca_runtime_input.py"
RENDERER_SHA256 = "4d75904eca4296abd892958df700063e272f5f029a0ce256bcf1b579768aeadb"

SIDECHAIN_HEAVY_ATOMS = {
    "ASP": {"CB", "CG", "OD1", "OD2"},
    "GLU": {"CB", "CG", "CD", "OE1", "OE2"},
    "ASN": {"CB", "CG", "OD1", "ND2"},
    "GLN": {"CB", "CG", "CD", "OE1", "NE2"},
    "SER": {"CB", "OG"},
    "THR": {"CB", "OG1", "CG2"},
    "TYR": {"CB", "CG", "CD1", "CE1", "CZ", "OH", "CE2", "CD2"},
    "HIS": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
}

MXAF36 = PROJECT / "workspaces/mxaf_1H4I_qm36_qm"
MONO33 = PROJECT / (
    "workspaces/tf_pqq_C5B120_pqq_la_monomer_v1_"
    "d7eb580a1486_rank0_sample_0_qm"
)
DIMER33 = PROJECT / (
    "workspaces/tf_pqq_C5B120_pqq_la_dimer_v1_"
    "cdf8e2570d6c_rank2_sample_2_qm"
)
BOUNDARY = PROJECT / "diagnostics/pqq_boundary_pair_20260914"

CASES = {
    "mxaf": {
        "source": PROJECT / "workspaces/mxaf_qm/1H4I_protonated.pdb",
        "source_sha256": "e91e6d51e0175a5f92e7eb39a446d65f3e869db7d475b643fc8a5eb358127bab",
        "chain": "A",
        "chain_charge": -8,
        "repair": None,
        "radii": {
            "3.3": {
                "manifest": BOUNDARY / "mxaf_1H4I_qm33/mxaf_1H4I_qm33_carve_manifest.json",
                "sha256": "95a97a9b07c87d8ca420ac86c65bb4e8914b3485664c66f3576e14480447d9d0",
            },
            "3.6": {
                "manifest": MXAF36 / "mxaf_1H4I_qm36_carve_manifest.json",
                "sha256": "0eb4944d69b5ab8897903a206259afe692891bab0e8261b15531001a87f68659",
            },
        },
    },
    "c5_monomer": {
        "source": MONO33 / "tf_pqq_C5B120_pqq_la_monomer_v1_d7eb580a1486_rank0_sample_0_protonated.pdb",
        "source_sha256": "1db4fe7872530e97fb3027ce23036096fc6fa5de49e6c9237def991397cc5531",
        "chain": "A",
        "chain_charge": -3,
        "repair": {"residue": "601", "bad": "OXT", "sound": "O"},
        "radii": {
            "3.3": {
                "manifest": MONO33 / "tf_pqq_C5B120_pqq_la_monomer_v1_d7eb580a1486_rank0_sample_0_carve_manifest.json",
                "sha256": "648af876bc24c0cfb10bbcd62233aea47f386b6389e8dcb27bbe1d77440c571e",
            },
            "3.6": {
                "manifest": BOUNDARY / "C5B120_monomer_qm36/C5B120_monomer_qm36_carve_manifest.json",
                "sha256": "04f27c7a910d996b0a3dbf5bf219aaa63c22afa2ec446befe5d12b69c3deb868",
            },
        },
    },
    "c5_dimer": {
        "source": DIMER33 / "tf_pqq_C5B120_pqq_la_dimer_v1_cdf8e2570d6c_rank2_sample_2_protonated.pdb",
        "source_sha256": "579febf0b4fb1aab9198916b55336258c74dccd3d0b01ee54aade650ac8aa170",
        "chain": "B",
        "chain_charge": -3,
        "repair": {"residue": "601", "bad": "O", "sound": "OXT"},
        "radii": {
            "3.3": {
                "manifest": DIMER33 / "tf_pqq_C5B120_pqq_la_dimer_v1_cdf8e2570d6c_rank2_sample_2_carve_manifest.json",
                "sha256": "f3cd9ffd71ef67e62889c9d30231436595e5334f178adf6f4581ead11a212844",
            },
            "3.6": {
                "manifest": BOUNDARY / "C5B120_dimer_qm36/C5B120_dimer_qm36_carve_manifest.json",
                "sha256": "82fa4bb4cc0ad72fd6e785ffeb6e2b86182ea6a30ae23e9c92cc69e79ac06d81",
            },
        },
    },
}

FRAGMENT_RE = re.compile(r"^([^:]+):([A-Z0-9]{3})(-?[0-9]+)([A-Za-z]?)$")
COULOMB_KCAL_A = 332.063713299


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise RuntimeError(f"refusing to replace divergent prepared artifact: {path}")
        return
    temporary = path.with_name("." + path.name + ".tmp")
    temporary.write_bytes(payload)
    temporary.replace(path)


def write_json(path: Path, value: object) -> None:
    write_bytes(path, (json.dumps(value, indent=2, sort_keys=True) + "\n").encode())


def xyz_tuple(position) -> tuple[float, float, float]:
    value = position.value_in_unit(unit.angstrom)
    return float(value.x), float(value.y), float(value.z)


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(a, factor):
    return tuple(x * factor for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def unit_vector(a):
    norm = math.sqrt(dot(a, a))
    if norm < 1.0e-10:
        raise RuntimeError("zero-length geometry vector")
    return scale(a, 1.0 / norm)


def angle(a, center, b) -> float:
    x = unit_vector(subtract(a, center))
    y = unit_vector(subtract(b, center))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot(x, y)))))


def atom_lookup(topology) -> dict[tuple[str, str, str], object]:
    result = {}
    for atom in topology.atoms():
        key = (atom.residue.chain.id, atom.residue.id, atom.name)
        if key in result:
            raise RuntimeError(f"duplicate atom key in isolated chain: {key}")
        result[key] = atom
    return result


def isolate_and_repair(case_name: str, case: dict) -> tuple[object, list, dict]:
    source = case["source"]
    if sha256(source) != case["source_sha256"]:
        raise RuntimeError(f"source hash drift: {case_name}")
    pdb = app.PDBFile(str(source))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.delete([atom for atom in modeller.topology.atoms() if atom.residue.chain.id != case["chain"]])
    positions = [openmm.Vec3(*xyz_tuple(position)) for position in modeller.positions]
    lookup = atom_lookup(modeller.topology)
    repair_record = None
    spec = case["repair"]
    if spec is not None:
        prefix = (case["chain"], spec["residue"])
        atoms = {name: lookup[(*prefix, name)] for name in ("CA", "C", spec["sound"], spec["bad"])}
        coords = {name: tuple(positions[atom.index]) for name, atom in atoms.items()}
        c = coords["C"]
        ca = coords["CA"]
        sound = coords[spec["sound"]]
        bad = coords[spec["bad"]]
        u = unit_vector(subtract(ca, c))
        v = unit_vector(subtract(sound, c))
        parallel = scale(u, dot(v, u))
        reflected = subtract(parallel, subtract(v, parallel))
        sound_length = distance(c, sound)
        repaired = add(c, scale(reflected, sound_length))
        positions[atoms[spec["bad"]].index] = openmm.Vec3(*repaired)
        repair_record = {
            "residue": f"{case['chain']}:ASN{spec['residue']}",
            "repaired_atom": spec["bad"],
            "sound_atom": spec["sound"],
            "original_C_bad_distance_A": distance(c, bad),
            "sound_C_O_distance_A": sound_length,
            "repaired_C_O_distance_A": distance(c, repaired),
            "CA_C_repaired_angle_deg": angle(ca, c, repaired),
            "sound_C_repaired_angle_deg": angle(sound, c, repaired),
            "repaired_xyz_A": list(repaired),
            "rule": "reflection of sound C-O vector about C-CA axis",
        }
        if not 1.15 <= repair_record["repaired_C_O_distance_A"] <= 1.40:
            raise RuntimeError(f"terminal repair failed geometry gate: {case_name}")

    chain_path = PREPARED / case_name / "target_chain.pdb"
    chain_path.parent.mkdir(parents=True, exist_ok=True)
    if not chain_path.exists():
        temporary = chain_path.with_name("." + chain_path.name + ".tmp")
        with temporary.open("w") as handle:
            app.PDBFile.writeFile(
                modeller.topology,
                unit.Quantity(positions, unit.angstrom),
                handle,
                keepIds=True,
            )
        temporary.replace(chain_path)
    reread = app.PDBFile(str(chain_path))
    chains = list(reread.topology.chains())
    if len(chains) != 1 or chains[0].id != case["chain"]:
        raise RuntimeError(f"isolated-chain identity failed: {case_name}")
    if spec is not None:
        check = atom_lookup(reread.topology)
        c = xyz_tuple(reread.positions[check[(case["chain"], spec["residue"], "C")].index])
        for oxygen in ("O", "OXT"):
            o = xyz_tuple(reread.positions[check[(case["chain"], spec["residue"], oxygen)].index])
            if not 1.15 <= distance(c, o) <= 1.40:
                raise RuntimeError(f"serialized terminal geometry failed: {case_name} {oxygen}")
    provenance = {
        "source": {"path": str(source), "sha256": sha256(source)},
        "target_chain": {"path": str(chain_path), "sha256": sha256(chain_path)},
        "terminal_repair": repair_record,
    }
    return reread.topology, list(reread.positions), provenance


def forcefield_charges(topology) -> tuple[list[float], float]:
    forcefield = app.ForceField(str(FORCEFIELD))
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
    )
    nonbonded = next(
        force for force in system.getForces() if isinstance(force, openmm.NonbondedForce)
    )
    charges = [
        float(nonbonded.getParticleParameters(index)[0].value_in_unit(unit.elementary_charge))
        for index in range(system.getNumParticles())
    ]
    if len(charges) != sum(1 for _ in topology.atoms()):
        raise RuntimeError("force-field particle/atom count mismatch")
    return charges, sum(charges)


def find_residue(topology, fragment_id: str):
    match = FRAGMENT_RE.fullmatch(fragment_id)
    if match is None:
        raise RuntimeError(f"cannot parse fragment id: {fragment_id}")
    chain_id, residue_name, residue_id, insertion = match.groups()
    matches = [
        residue
        for residue in topology.residues()
        if residue.chain.id == chain_id
        and residue.name == residue_name
        and residue.id == residue_id
        and (not insertion or getattr(residue, "insertionCode", "") == insertion)
    ]
    if len(matches) != 1:
        raise RuntimeError(f"fragment residue lookup is not unique: {fragment_id}")
    return matches[0]


def parse_xyz(path: Path) -> list[tuple[str, tuple[float, float, float]]]:
    lines = path.read_text().splitlines()
    if int(lines[0]) != len(lines) - 2:
        raise RuntimeError(f"malformed XYZ: {path}")
    atoms = []
    for line in lines[2:]:
        fields = line.split()
        atoms.append((fields[0], tuple(float(value) for value in fields[1:4])))
    return atoms


def build_environment(
    case_name: str,
    case: dict,
    radius: str,
    manifest_path: Path,
    topology,
    positions,
    base_charges: list[float],
    chain_charge: float,
) -> tuple[dict, list[dict]]:
    manifest_info = case["radii"][radius]
    if sha256(manifest_path) != manifest_info["sha256"]:
        raise RuntimeError(f"anchor-manifest hash drift: {case_name} {radius}")
    manifest = read_json(manifest_path)
    if abs(float(manifest["qm_inclusion"]["cutoff_A"]) - float(radius)) > 1.0e-9:
        raise RuntimeError(f"QM radius mismatch: {case_name} {radius}")

    atoms = list(topology.atoms())
    coords = [xyz_tuple(position) for position in positions]
    adjacency = {atom.index: set() for atom in atoms}
    for first, second in topology.bonds():
        adjacency[first.index].add(second.index)
        adjacency[second.index].add(first.index)

    removed: set[int] = set()
    shifts = {atom.index: 0.0 for atom in atoms}
    ledgers = []
    protein_fragments = [
        fragment for fragment in manifest["qm_fragments"] if fragment["kind"].startswith("protein_")
    ]
    if any(fragment["kind"] != "protein_sidechain" for fragment in protein_fragments):
        raise RuntimeError("this frozen panel only supports side-chain protein fragments")
    if len({fragment["id"] for fragment in protein_fragments}) != len(protein_fragments):
        raise RuntimeError("duplicate protein fragment")

    link_hydrogens = []
    for fragment in protein_fragments:
        residue = find_residue(topology, fragment["id"])
        residue_atoms = list(residue.atoms())
        by_name = {atom.name: atom for atom in residue_atoms}
        heavy_names = SIDECHAIN_HEAVY_ATOMS.get(residue.name)
        if heavy_names is None or not heavy_names <= set(by_name):
            raise RuntimeError(f"incomplete fragment source residue: {fragment['id']}")
        source_qm = {by_name[name].index for name in heavy_names}
        for index in list(source_qm):
            source_qm.update(
                neighbor
                for neighbor in adjacency[index]
                if atoms[neighbor].element.symbol.upper() == "H"
            )
        if len(source_qm) != int(fragment["atom_count"]) - 1:
            raise RuntimeError(
                f"topology-derived source atom count disagrees with carve: {fragment['id']}"
            )
        ca = by_name["CA"]
        if ca.index in source_qm:
            raise RuntimeError(f"CA unexpectedly belongs to QM source fragment: {fragment['id']}")
        local_removed = source_qm | {ca.index}
        if removed & local_removed:
            raise RuntimeError(f"overlapping QM/MM cuts: {fragment['id']}")
        removed.update(local_removed)
        recipients = [
            atom
            for atom in residue_atoms
            if atom.name in {"N", "C"} and atom.index in adjacency[ca.index]
        ]
        if sorted(atom.name for atom in recipients) != ["C", "N"]:
            raise RuntimeError(f"missing bonded N/C correction recipients: {fragment['id']}")
        q_residue = sum(base_charges[atom.index] for atom in residue_atoms)
        retained_before_shift = sum(
            base_charges[atom.index]
            for atom in residue_atoms
            if atom.index not in local_removed
        )
        target_retained = q_residue - float(fragment["formal_charge"])
        correction = target_retained - retained_before_shift
        increment = correction / 2.0
        for recipient in recipients:
            shifts[recipient.index] += increment

        cb_xyz = coords[by_name["CB"].index]
        ca_xyz = coords[ca.index]
        link_xyz = add(cb_xyz, scale(unit_vector(subtract(ca_xyz, cb_xyz)), 1.09))
        link_hydrogens.append((fragment["id"], link_xyz))
        ha = by_name.get("HA")
        if ha is None or ha.index in removed:
            raise RuntimeError(f"retained backbone HA missing: {fragment['id']}")
        ledgers.append(
            {
                "fragment_id": fragment["id"],
                "fragment_formal_charge": fragment["formal_charge"],
                "source_qm_atoms": [
                    {
                        "index": index,
                        "name": atoms[index].name,
                        "element": atoms[index].element.symbol,
                        "original_charge_e": base_charges[index],
                    }
                    for index in sorted(source_qm)
                ],
                "removed_MM1_CA": {
                    "index": ca.index,
                    "original_charge_e": base_charges[ca.index],
                },
                "full_residue_ff_charge_e": q_residue,
                "target_retained_PC_charge_e": target_retained,
                "boundary_correction_e": correction,
                "recipients": [
                    {
                        "name": recipient.name,
                        "index": recipient.index,
                        "original_charge_e": base_charges[recipient.index],
                        "increment_e": increment,
                        "final_charge_e": base_charges[recipient.index] + increment,
                    }
                    for recipient in sorted(recipients, key=lambda item: item.name)
                ],
                "retained_HA": {
                    "index": ha.index,
                    "charge_e": base_charges[ha.index],
                    "distance_to_link_H_A": distance(coords[ha.index], link_xyz),
                },
                "link_H_xyz_A": list(link_xyz),
            }
        )

    points = []
    for atom in atoms:
        if atom.index in removed:
            continue
        points.append(
            {
                "index": atom.index,
                "chain": atom.residue.chain.id,
                "residue": atom.residue.id,
                "resname": atom.residue.name,
                "atom": atom.name,
                "element": atom.element.symbol,
                "charge_e": base_charges[atom.index] + shifts[atom.index],
                "xyz_A": list(coords[atom.index]),
            }
        )
    pc_sum = sum(point["charge_e"] for point in points)
    qm_protein_charge = sum(float(fragment["formal_charge"]) for fragment in protein_fragments)
    expected_pc_sum = chain_charge - qm_protein_charge
    if abs(pc_sum - expected_pc_sum) > 1.0e-6:
        raise RuntimeError(f"point-charge closure failed: {case_name} {radius}")

    environment_dir = PREPARED / case_name / f"qm{radius.replace('.', '')}"
    environment_dir.mkdir(parents=True, exist_ok=True)
    pc_path = environment_dir / "pointcharges.pc"
    pc_lines = [str(len(points))]
    pc_lines.extend(
        f"{point['charge_e']:.12f} {point['xyz_A'][0]:.10f} {point['xyz_A'][1]:.10f} {point['xyz_A'][2]:.10f}"
        for point in points
    )
    write_bytes(pc_path, ("\n".join(pc_lines) + "\n").encode())
    serialized_sum = sum(float(line.split()[0]) for line in pc_lines[1:])
    if abs(serialized_sum - expected_pc_sum) > 1.0e-6:
        raise RuntimeError(f"serialized point-charge closure failed: {case_name} {radius}")

    tasks = []
    qm_coordinates = {}
    for metal in ("La", "Ca"):
        anchor_input_record = manifest["outputs"][f"{metal}_input"]
        anchor_xyz_record = manifest["outputs"][f"{metal}_xyz"]
        anchor_input = manifest_path.parent / anchor_input_record["path"]
        anchor_xyz = manifest_path.parent / anchor_xyz_record["path"]
        if sha256(anchor_input) != anchor_input_record["sha256"]:
            raise RuntimeError(f"anchor input hash mismatch: {case_name} {radius} {metal}")
        if sha256(anchor_xyz) != anchor_xyz_record["sha256"]:
            raise RuntimeError(f"anchor XYZ hash mismatch: {case_name} {radius} {metal}")
        xyz_path = environment_dir / f"{metal}.xyz"
        write_bytes(xyz_path, anchor_xyz.read_bytes())
        qm_coordinates[metal] = parse_xyz(xyz_path)
        input_path = environment_dir / f"{metal}.inp"
        text = anchor_input.read_text()
        if "%pointcharges" in text.lower() or "%pal" in text.lower():
            raise RuntimeError(f"anchor input already has embedding/parallel controls: {anchor_input}")
        lines = text.splitlines(keepends=True)
        coordinate_index = next(
            (index for index, line in enumerate(lines) if line.lstrip().lower().startswith("* xyzfile")),
            None,
        )
        if coordinate_index is None:
            raise RuntimeError(f"anchor input lacks xyzfile directive: {anchor_input}")
        fields = lines[coordinate_index].split()
        fields[-1] = xyz_path.name
        lines[coordinate_index] = " ".join(fields) + "\n"
        lines.insert(coordinate_index, '%pointcharges "pointcharges.pc"\n')
        write_bytes(input_path, "".join(lines).encode())
        output_path = environment_dir / f"{metal}.out"
        if output_path.exists() or Path(str(output_path) + ".execution.json").exists():
            raise RuntimeError(f"refusing an existing embedded energy: {output_path}")
        tasks.append(
            {
                "task_id": f"{case_name}_qm{radius.replace('.', '')}_{metal}",
                "input": {
                    "path": str(input_path.relative_to(HERE)),
                    "sha256": sha256(input_path),
                    "anchor": {"path": str(anchor_input), "sha256": sha256(anchor_input)},
                },
                "xyz": {
                    "path": str(xyz_path.relative_to(HERE)),
                    "sha256": sha256(xyz_path),
                    "anchor": {"path": str(anchor_xyz), "sha256": sha256(anchor_xyz)},
                },
                "point_charges": {
                    "path": str(pc_path.relative_to(HERE)),
                    "sha256": sha256(pc_path),
                },
                "output_path": str(output_path.relative_to(HERE)),
            }
        )

    la_atoms = qm_coordinates["La"]
    ca_atoms = qm_coordinates["Ca"]
    if len(la_atoms) != len(ca_atoms):
        raise RuntimeError(f"La/Ca QM atom count mismatch: {case_name} {radius}")
    for index, (la, ca) in enumerate(zip(la_atoms, ca_atoms)):
        if la[1] != ca[1] or (index > 0 and la[0] != ca[0]):
            raise RuntimeError(f"La/Ca coordinates are not vertical: {case_name} {radius}")
    if la_atoms[0][0].upper() != "LA" or ca_atoms[0][0].upper() != "CA":
        raise RuntimeError(f"first QM center is not the swapped metal: {case_name} {radius}")

    all_qm = [coordinate for _, coordinate in la_atoms]
    min_pc_qm = min(
        distance(tuple(point["xyz_A"]), coordinate)
        for point in points
        for coordinate in all_qm
    )
    min_pc_link = min(
        distance(tuple(point["xyz_A"]), link_xyz)
        for point in points
        for _, link_xyz in link_hydrogens
    )
    if min_pc_qm < 1.0 - 1.0e-9:
        raise RuntimeError(f"PC-to-QM distance gate failed: {case_name} {radius}")

    metal_xyz = tuple(float(value) for value in manifest["selected_site"]["xyz_A"])
    potential_bins = {"le_6_A": 0.0, "6_to_12_A": 0.0, "gt_12_A": 0.0}
    for point in points:
        radius_A = distance(tuple(point["xyz_A"]), metal_xyz)
        key = "le_6_A" if radius_A <= 6.0 else "6_to_12_A" if radius_A <= 12.0 else "gt_12_A"
        potential_bins[key] += COULOMB_KCAL_A * point["charge_e"] / radius_A

    environment = {
        "case": case_name,
        "radius_A": float(radius),
        "anchor_manifest": {"path": str(manifest_path), "sha256": sha256(manifest_path)},
        "chain_charge_e": chain_charge,
        "qm_protein_fragment_charge_e": qm_protein_charge,
        "point_charge_count": len(points),
        "point_charge_sum_e": pc_sum,
        "expected_point_charge_sum_e": expected_pc_sum,
        "combined_chain_charge_e": pc_sum + qm_protein_charge,
        "point_charges": {"path": str(pc_path), "sha256": sha256(pc_path)},
        "boundary_ledger": ledgers,
        "removed_atom_indices": sorted(removed),
        "minimum_PC_to_any_QM_center_A": min_pc_qm,
        "minimum_PC_to_link_H_A": min_pc_link,
        "metal_potential_kcal_mol_per_e": {
            **potential_bins,
            "total": sum(potential_bins.values()),
        },
    }
    return environment, tasks


def prepare() -> None:
    for path, expected in (
        (FORCEFIELD, FORCEFIELD_SHA256),
        (RUNNER, RUNNER_SHA256),
        (RENDERER, RENDERER_SHA256),
    ):
        if sha256(path) != expected:
            raise RuntimeError(f"pinned implementation drift: {path}")
    experiment = HERE / "EXPERIMENT.md"
    amendment = HERE / "PREPARATION_AMENDMENT.md"
    environments = {}
    tasks = []
    structures = {}
    for case_name, case in CASES.items():
        topology, positions, provenance = isolate_and_repair(case_name, case)
        charges, observed_chain_charge = forcefield_charges(topology)
        if abs(observed_chain_charge - float(case["chain_charge"])) > 1.0e-6:
            raise RuntimeError(f"ff19SB chain charge mismatch: {case_name}")
        provenance.update(
            {
                "chain_id": case["chain"],
                "atom_count": len(charges),
                "ff19SB_charge_e": observed_chain_charge,
            }
        )
        structures[case_name] = provenance
        environments[case_name] = {}
        for radius in ("3.3", "3.6"):
            manifest_path = case["radii"][radius]["manifest"]
            environment, radius_tasks = build_environment(
                case_name,
                case,
                radius,
                manifest_path,
                topology,
                positions,
                charges,
                observed_chain_charge,
            )
            environments[case_name][radius] = environment
            tasks.extend(radius_tasks)

        low = environments[case_name]["3.3"]
        high = environments[case_name]["3.6"]
        if abs((high["point_charge_sum_e"] - low["point_charge_sum_e"]) - 1.0) > 1.0e-6:
            raise RuntimeError(f"3.3/3.6 PC charge transfer is not +1 e: {case_name}")
        if abs((high["qm_protein_fragment_charge_e"] - low["qm_protein_fragment_charge_e"]) + 1.0) > 1.0e-6:
            raise RuntimeError(f"3.3/3.6 QM charge transfer is not -1 e: {case_name}")
        if abs(high["combined_chain_charge_e"] - low["combined_chain_charge_e"]) > 1.0e-6:
            raise RuntimeError(f"combined charge changes with QM radius: {case_name}")

    task_manifest = {
        "schema_version": "pqq_balanced_embedding.tasks.v1",
        "experiment": {"path": str(experiment), "sha256": sha256(experiment), "git_commit": "fafd9f7"},
        "preparation_amendment": {"path": str(amendment), "sha256": sha256(amendment), "git_commit": "be1e051"},
        "implementation": {
            "preparer": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
            "execution_runner": {"path": str(RUNNER), "sha256": sha256(RUNNER)},
            "runtime_renderer": {"path": str(RENDERER), "sha256": sha256(RENDERER)},
            "openmm_version": openmm.__version__,
            "forcefield": {"path": str(FORCEFIELD), "sha256": sha256(FORCEFIELD)},
        },
        "execution_policy": {
            "task_runner": {"path": str(RUNNER), "sha256": sha256(RUNNER)},
            "runtime_renderer": {"path": str(RENDERER), "sha256": sha256(RENDERER)},
        },
        "tasks": tasks,
    }
    write_json(TASK_MANIFEST, task_manifest)
    preparation = {
        "schema_version": "pqq_balanced_embedding.preparation.v1",
        "status": "PASS",
        "embedded_energies_present_at_preparation": False,
        "required_orca_legs": len(tasks),
        "structures": structures,
        "environments": environments,
        "task_manifest": {"path": str(TASK_MANIFEST), "sha256": sha256(TASK_MANIFEST)},
    }
    write_json(PREPARATION, preparation)
    verify()
    print(json.dumps({"status": "PASS", "orca_legs": len(tasks), "task_manifest": str(TASK_MANIFEST)}, indent=2))


def verify() -> None:
    manifest = read_json(TASK_MANIFEST)
    implementation = manifest["implementation"]
    for record in (
        implementation["preparer"],
        implementation["execution_runner"],
        implementation["runtime_renderer"],
        implementation["forcefield"],
        manifest["experiment"],
        manifest["preparation_amendment"],
    ):
        if sha256(Path(record["path"])) != record["sha256"]:
            raise RuntimeError(f"prepared implementation/artifact drift: {record['path']}")
    for task in manifest["tasks"]:
        for key in ("input", "xyz", "point_charges"):
            record = task[key]
            path = HERE / record["path"]
            if sha256(path) != record["sha256"]:
                raise RuntimeError(f"prepared task artifact drift: {task['task_id']} {key}")
        input_path = HERE / task["input"]["path"]
        if '%pointcharges "pointcharges.pc"' not in input_path.read_text():
            raise RuntimeError(f"point-charge directive missing: {task['task_id']}")
    preparation = read_json(PREPARATION)
    if preparation["status"] != "PASS" or preparation["required_orca_legs"] != 12:
        raise RuntimeError("preparation record is not a 12-leg PASS")
    if sha256(TASK_MANIFEST) != preparation["task_manifest"]["sha256"]:
        raise RuntimeError("task-manifest hash does not match preparation record")
    print(json.dumps({"status": "PASS", "verified_tasks": len(manifest["tasks"])}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "verify"))
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    else:
        verify()


if __name__ == "__main__":
    main()
