#!/usr/bin/env python3
"""Prepare reviewed PLM monomer sites with the unchanged fixed-core-v3 chemistry.

No model selection, residue discovery, ORCA execution, or scheduler submission.
The manifest names every selected model and every core residue explicitly.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

import gemmi

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CAL = ROOT / "diagnostics/pqq_pmdh_fixed_core_calibration_20260914"
PROTOCOL = "pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3"
STANDARD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


fixed = load_module("plm_original_fixed_core", CAL / "fixed_core_carver.py")
protonator = load_module("plm_original_standard_protonator", CAL / "protonate_standard_only.py")
base = fixed.base


def record(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(item):
    if record(item["path"]) != item:
        raise ValueError(f"File hash/path mismatch: {item['path']}")
    return Path(item["path"])


def read(path):
    return json.loads(Path(path).read_text())


def write(path, payload):
    with Path(path).open("x") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def verify_pins(path):
    pins = read(path)
    if pins["protocol_id"] != PROTOCOL or pins["schema_version"] != "plm.adh9.fixed_core_pins.v1":
        raise ValueError("Wrong candidate pins")
    if verify(pins["wrapper"]) != Path(__file__).resolve():
        raise ValueError("Executing wrapper differs from pinned wrapper")
    cp = verify(pins["calibration_implementation_pins"])
    original = fixed.verify_pins(cp)
    release = read(verify(pins["calibration_result"]))
    if release["calibration"]["verdict"] != "CALIBRATABLE":
        raise ValueError("Calibration release does not pass")
    expected = {"minimum_confidence": 0.9, "minimum_CN": 7, "cutoff_A": 3.1,
                "maximum_direct_N": 2, "maximum_Asp_partner_A": 3.5}
    if pins["candidate_gates"] != expected:
        raise ValueError("Candidate gates changed")
    verify(pins["candidate_manifest"])
    for rec in pins["authorities"].values():
        verify(rec)
    import importlib.metadata
    for package, version in original["preparation_runtime"]["packages"].items():
        if importlib.metadata.version(package) != version:
            raise ValueError(f"Preparation package differs: {package}")
    if sys.version.split()[0] != original["preparation_runtime"]["python"]["version"]:
        raise ValueError("Preparation Python version differs")
    return pins, original, cp


def heavy_map(structure):
    result = {}
    if len(structure) != 1:
        raise ValueError("Expected one structural model")
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                if atom.element.name in {"H", "D"}:
                    continue
                if atom.altloc not in {"\0", " ", ""}:
                    raise ValueError("Alternate conformers are unsupported")
                key = (chain.name, residue.name, residue.seqid.num,
                       residue.seqid.icode.strip(), atom.name, atom.element.name)
                if key in result:
                    raise ValueError(f"Duplicate heavy atom {key}")
                result[key] = (atom.pos.x, atom.pos.y, atom.pos.z)
    return result


def compare_heavy(left, right):
    if left.keys() != right.keys():
        raise ValueError("Heavy atoms added, removed, or renamed")
    movement = max(math.dist(left[k], right[k]) for k in left)
    if movement > 0.001:
        raise ValueError(f"Heavy coordinate displacement {movement} exceeds PDB serialization tolerance")
    return {"atom_count": len(left), "maximum_displacement_A": movement,
            "tolerance_A": 0.001, "passes": True}


def role_state(model, target):
    roles = target["roles"]
    if set(roles) != set(fixed.ROLE_ORDER):
        raise ValueError("Exactly five explicit reviewed core roles are required")
    keys = {role: fixed.residue_key(roles[role], role) for role in fixed.ROLE_ORDER}
    for role, key in keys.items():
        if key.resname not in fixed.ROLE_ALLOWED_RESNAMES[role]:
            raise ValueError(f"Unsupported {role}: {key.resname}")
        if key.chain != target["protein_chain"]:
            raise ValueError("Core role is outside the explicitly selected monomer")
    cat = keys["catalytic_aspartate"]
    extra = keys["extra_acidic_ligand_homolog"]
    if extra.resnum != cat.resnum + 2 or extra.icode != cat.icode:
        raise ValueError("Mapped D+2 is not exactly catalytic-Asp+2")
    residues = {role: fixed.resolve_residue(model, key, role) for role, key in keys.items()}
    partner = fixed.asp_partner_geometry(residues["catalytic_aspartate"],
                                         residues["catalytic_asp_cationic_partner"])
    return keys, residues, partner


def site_state(model, target):
    selector = target["metal"]
    site = base.select_metal_sites(model, site_chain=selector["chain"],
        site_resnum=selector["resnum"], site_icode=selector.get("icode", ""),
        site_resname=selector["resname"], site_atom=selector["atom"])[0]
    if site.element.upper() != "LA":
        raise ValueError("Candidate must be an existing La-conditioned model")
    pqq = base._prepare_nearby_pqq(model, site.atom.pos, fixed.DEFAULT_PQQ_MICROSTATE)
    base._reject_unsupported_nearby_species(model, site, pqq)
    pkey = fixed.residue_key(target["pqq"], "PQQ")
    if set(pqq) != {pkey}:
        raise ValueError("Exactly the explicitly selected complete PQQ is required")
    residue, prepared = pqq[pkey]
    if prepared.microstate.formal_charge != -3 or prepared.formula != "C14H3N2O8":
        raise ValueError("PQQ does not match oxidized PQQ(3-) chemistry")
    if prepared.schema_id != "pdb_ccd_pqq_v1":
        raise ValueError("This PLM wrapper requires the reviewed CCD-PQQ source schema")
    contacts, _, sulfur, untyped = base._classify_contacts(model, site, pqq)
    n = sum(c.element.upper() == "N" for c in contacts)
    if len(contacts) < 7 or n > 2:
        raise ValueError(f"Production donor gate failed: CN={len(contacts)}, direct N={n}")
    return site, pkey, residue, prepared, contacts, sulfur, untyped


def prepare_one(target, output, pins, original, calibration_pins):
    source = verify(target["source_cif"])
    summary_path = verify(target["summary_confidence"])
    if summary_path != source.with_name(source.stem + "_summary.json"):
        raise ValueError("Confidence summary does not name this exact ranked source")
    for rec in target["source_provenance"].values():
        verify(rec)
    structure = gemmi.read_structure(str(source))
    heavy = heavy_map(structure)
    model = structure[0]
    chains = [c.name for c in model]
    if len(set(chains)) != len(chains):
        raise ValueError("Duplicate source chains")
    if set(chains) != {target["protein_chain"], target["metal"]["chain"], target["pqq"]["chain"]}:
        raise ValueError("Expected exactly the selected monomer, La, and PQQ chains")
    for chain in model:
        for residue in chain:
            expected = STANDARD if chain.name == target["protein_chain"] else {"LA", "PQQ"}
            if residue.name not in expected:
                raise ValueError(f"Unsupported source residue {residue.name}; no silent exclusion")
            if any(a.element.name in {"H", "D"} for a in residue):
                raise ValueError("Raw source already contains hydrogen")
    summary = read(summary_path)
    confidence = float(summary["chain_pair_iptm"][chains.index(target["protein_chain"])][chains.index(target["metal"]["chain"])] )
    if not math.isfinite(confidence) or confidence < 0.9:
        raise ValueError(f"Protein–La confidence gate failed: {confidence}")
    keys, residues, raw_partner = role_state(model, target)
    site, pkey, pqq_residue, prepared, contacts, sulfur, untyped = site_state(model, target)
    acidic = keys["extra_acidic_ligand_homolog"].resname in {"ASP", "GLU"}
    included = set(fixed.ROLE_ORDER)
    if not acidic:
        included.remove("extra_acidic_ligand_homolog")
    included_keys = {pkey, *(keys[role] for role in included)}
    if any(c.residue not in included_keys for c in contacts):
        raise ValueError("A qualifying direct donor lies outside the fixed homologous core")
    stem = target["case_id"]
    directory = output / stem
    directory.mkdir(exist_ok=False)
    normalized = directory / f"{stem}_normalized.pdb"
    structure.write_pdb(str(normalized))
    normalized_heavy = heavy_map(gemmi.read_structure(str(normalized)))
    first_check = compare_heavy(heavy, normalized_heavy)
    normalization = directory / f"{stem}_normalization_manifest.json"
    write(normalization, {"protocol_id": PROTOCOL, "source": record(source),
        "output": record(normalized), "normalizer": pins["wrapper"],
        "policy": "already_CCD_PQQ_LA_serialization_only_no_identity_or_coordinate_edit",
        "heavy_coordinate_check": first_check})
    protonated = directory / f"{stem}_protonated.pdb"
    protonation = protonator.protonate_standard_only(normalized, protonated,
        pins_path=calibration_pins, ph=7.0)
    pm = read(protonation)
    if pm["repaired_missing_atom_count"] or pm["repaired_missing_terminal_atom_count"]:
        raise ValueError("Zero-heavy-repair policy failed")
    if pm["software"] != original["preparation_runtime"]["packages"]:
        raise ValueError("Protonation runtime differs from calibration")
    protonated_structure = gemmi.read_structure(str(protonated))
    checks = {"raw_to_normalized": first_check,
        "normalized_to_protonated": compare_heavy(normalized_heavy, heavy_map(protonated_structure)),
        "raw_to_protonated": compare_heavy(heavy, heavy_map(protonated_structure))}
    heavy_path = directory / f"{stem}_heavy_coordinate_check.json"
    write(heavy_path, {"schema_version": "alchemical_bvs.heavy_coordinate_check.v1",
        "protocol_id": PROTOCOL, "passes": True, "all_source_heavy_atoms_retained": True,
        "no_heavy_atoms_added": True, "all_retained_source_heavy_atoms_preserved": True,
        "source": record(source), "protonated": record(protonated), "checks": checks})
    model = protonated_structure[0]
    keys, residues, partner = role_state(model, target)
    site, pkey, pqq_residue, prepared, qcontacts, sulfur, untyped = site_state(model, target)
    def signature(values):
        return sorted((c.residue.label(), c.atom_name, c.element) for c in values)
    if signature(contacts) != signature(qcontacts):
        raise ValueError("PDB serialization/protonation changed source donor membership")
    if any(c.residue not in included_keys for c in qcontacts):
        raise ValueError("Protonated qualifying donor lies outside core")
    if any(a.element.name in {"H", "D"} for a in pqq_residue):
        raise ValueError("PDBFixer supplied PQQ hydrogen")
    fragments = [{"kind": "fixed_core_pqq", "role": "pqq_cofactor", "id": pkey.label(),
        "formal_charge": -3, "atom_count": len(prepared.atoms),
        "atom_records": fixed.pqq_atom_records(pqq_residue, prepared),
        "microstate_id": prepared.microstate.microstate_id}]
    atoms = list(prepared.atoms)
    for role in fixed.ROLE_ORDER:
        if role not in included:
            continue
        if role == "catalytic_asp_cationic_partner":
            frag_atoms, fragment = fixed.cationic_sidechain_fragment(residues[role], keys[role], site.atom.pos)
        else:
            frag_atoms, fragment = fixed.canonical_sidechain_fragment(residues[role], keys[role], role, site.atom.pos)
        atoms.extend(frag_atoms)
        fragments.append(fragment)
    scaffold_charge = sum(f["formal_charge"] for f in fragments)
    charges = {"La": scaffold_charge + 3, "Ca": scaffold_charge + 2}
    expected = {"La": -2, "Ca": -3} if acidic else {"La": -1, "Ca": -2}
    if charges != expected:
        raise ValueError("Fixed-core charge ledger differs from calibration")
    tasks, outputs = [], {}
    for metal in ("La", "Ca"):
        metal_atoms = [(metal, site.atom.pos.x, site.atom.pos.y, site.atom.pos.z), *atoms]
        base._validate_singlet(metal, metal_atoms, charges[metal])
        xyz = directory / f"{stem}_{metal}_qm.xyz"
        inp = directory / f"sp_{stem}_{metal}.inp"
        fixed.write_xyz(xyz, stem=stem, label=metal, atoms=metal_atoms, charge=charges[metal])
        fixed.write_orca_input(inp, xyz_name=xyz.name, charge=charges[metal], stem=stem, label=metal)
        xr, ir = record(xyz), record(inp)
        xr["path"], ir["path"] = xyz.name, inp.name
        outputs[metal + "_xyz"], outputs[metal + "_input"] = xr, ir
        tasks.append({"task_id": metal, "input": ir, "xyz": xr, "output_path": inp.with_suffix(".out").name})
    la = (directory / f"{stem}_La_qm.xyz").read_text().splitlines()
    ca = (directory / f"{stem}_Ca_qm.xyz").read_text().splitlines()
    if la[3:] != ca[3:] or la[2].split()[1:] != ca[2].split()[1:]:
        raise ValueError("La/Ca nuclear coordinate invariant failed")
    helpers = original["canonical_helpers"]
    manifest = {"schema_version": "alchemical_bvs.carve_manifest.v1", "protocol_id": PROTOCOL,
        "status": "ready_for_orca", "stem": stem, "target_id": target["target_id"],
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "candidate_manifest": pins["candidate_manifest"], "wrapper": pins["wrapper"],
        "calibration_authority": pins["calibration_result"],
        "calibration_implementation_pins": pins["calibration_implementation_pins"],
        "source_cif": record(source), "source_structure": record(protonated),
        "normalization_manifest": record(normalization), "normalized_source_structure": record(normalized),
        "protonation_manifest": record(protonation),
        "heavy_coordinate_check": {**record(heavy_path), "passes": True},
        "confidence": {"value": confidence, "basis": "protein_ion_chain_pair_iptm", "minimum": 0.9,
                       "source": target["summary_confidence"]},
        "site_selection": {"mode": "explicit_reviewed_candidate_selector", "requested": target["metal"]},
        "selected_site": {"chain": site.chain, "resnum": site.resnum, "icode": site.icode.strip(),
            "resname": site.resname, "atom_name": site.atom_name, "source_element": site.element},
        "coordination_policy": {"policy_id": "typed_cn7_3p1A_maxN2_candidate_v1", "cutoff_A": 3.1,
            "minimum_coordination_number": 7, "maximum_direct_nitrogen": 2,
            "forced_core_atoms_count_toward_coordination": False},
        "coordination": {"coordination_number": len(qcontacts),
            "oxygen_count": sum(c.element.upper() == "O" for c in qcontacts),
            "direct_nitrogen_count": sum(c.element.upper() == "N" for c in qcontacts),
            "passes_geometry": True,
            "source_typed_direct_donors_ordered_by_distance": [fixed.contact_metadata(c, frozenset(included_keys)) for c in contacts],
            "prepared_typed_direct_donors_ordered_by_distance": [fixed.contact_metadata(c, frozenset(included_keys)) for c in qcontacts]},
        "fixed_core": {"policy_id": fixed.CORE_POLICY_ID, "selection_is_distance_independent": True,
            "requested_roles": target["roles"], "included_roles": sorted(included),
            "included_residue_keys": sorted(k.label() for k in included_keys),
            "unexpected_direct_donors": [], "catalytic_Asp_partner_geometry": partner,
            "raw_catalytic_Asp_partner_geometry": raw_partner,
            "water_policy": "dry_exclude_all_source_and_synthetic_waters", "source_water_count": 0,
            "synthetic_water_count": 0, "replacement_ligand_count": 0,
            "point_charge_embedding": False, "geometry_relaxation": False},
        "qm_fragments": fragments,
        "pqq": {"present": True, "microstate_id": prepared.microstate.microstate_id,
            "formal_charge": -3, "multiplicity": 1, "residues": [prepared.metadata()]},
        "charge_ledger": {"fragments": fragments, "scaffold_formal_charge": scaffold_charge,
            "expected_total_charges": expected, "La_total": charges["La"], "Ca_total": charges["Ca"], "multiplicity": 1},
        "paired_arm_invariant": {"nonmetal_coordinates_byte_identical": True},
        "electronic_structure": {"method_id": fixed.METHOD_ID, "method": "r2SCAN-3c",
            "basis_policy": "native_r2scan3c_def2_mTZVPP_with_native_ECP", "solvation": "CPCM(Water)",
            "grid": "DefGrid3", "aquo_reference": original["aquo_reference"],
            "score_formula": "((E_Ca_site-E_La_site)-delta_E_aquo)*627.509474"},
        "execution_policy": {"id": "two_leg_allocation_aware_manifested_orca_v1",
            "task_runner": helpers["run_orca_task_manifest"], "runtime_renderer": helpers["render_orca_runtime_input"],
            "heavy_legs": ["La", "Ca"], "apo_leg": False, "water_leg": False, "omp_threads_per_rank": 1},
        "tasks": tasks, "outputs": outputs}
    path = directory / f"{stem}_carve_manifest.json"
    write(path, manifest)
    return record(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--implementation-pins", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pins, original, cp = verify_pins(args.implementation_pins)
    if record(args.manifest) != pins["candidate_manifest"]:
        raise ValueError("Manifest differs from pinned selector authority")
    spec = read(args.manifest)
    if spec["schema_version"] != "plm.adh9.candidate_manifest.v1":
        raise ValueError("Unsupported input manifest schema")
    output = args.output.resolve()
    if ROOT.resolve() / "workspaces" not in output.parents:
        raise ValueError("All candidate output must be under Alquemia workspaces")
    output.mkdir(exist_ok=False, parents=True)
    cases = []
    ids = [t["case_id"] for t in spec["targets"]]
    if len(set(ids)) != len(ids) or any(not re.fullmatch(r"[A-Za-z0-9_-]+", x) for x in ids):
        raise ValueError("Duplicate or unsafe case IDs")
    for target in spec["targets"]:
        case = {k: target[k] for k in ("case_id", "target_id", "rank", "source_cif")}
        try:
            case["carve_manifest"] = prepare_one(target, output, pins, original, cp)
            case["status"] = "ready_for_orca"
            case["reason"] = None
        except Exception as exc:
            case.update(status="unsupported", reason=f"{type(exc).__name__}: {exc}", carve_manifest=None)
        cases.append(case)
        print(case["case_id"], case["status"], case["reason"], flush=True)
    report = {"schema_version": "plm.adh9.prepared_pairs.v1", "protocol_id": PROTOCOL,
        "candidate_manifest": record(args.manifest), "implementation_pins": record(args.implementation_pins),
        "cases": cases, "prepared_pair_count": sum(c["status"] == "ready_for_orca" for c in cases),
        "unsupported_count": sum(c["status"] == "unsupported" for c in cases), "orca_executed": False}
    write(output / "prepared_pairs.json", report)
    print(output / "prepared_pairs.json")


if __name__ == "__main__":
    main()
