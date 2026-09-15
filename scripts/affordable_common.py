"""Strict, unit-explicit records shared by opt-in discriminator development."""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
from pathlib import Path

HA_TO_KCAL = 627.509474  # exact factor in the released PQQ calibration
BOHR_TO_A = 0.529177210903


class InvalidArtifact(ValueError):
    pass


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def record(path):
    p = Path(path).resolve()
    return {"path": str(p), "sha256": digest(p)}


def verify(rec):
    p = Path(rec["path"])
    if not p.is_file() or digest(p) != rec["sha256"]:
        raise InvalidArtifact(f"missing or changed artifact: {p}")
    return p


def read_json(path):
    return json.loads(Path(path).read_text())


def write_new(path, value):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("x") as f:
        json.dump(value, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")


def cache_key(settings):
    """Hash the complete scientific configuration, never just the site ID."""
    return hashlib.sha256(json.dumps(settings, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def snapshot_implementation(output):
    """Preserve actual live dependencies without claiming their edits as ours."""
    directory=Path(__file__).resolve().parent
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    files=sorted(directory.glob('affordable_*.py'))+[directory/name for name in (
        'run_orca_task_manifest.py','render_orca_runtime_input.py','carve_generic.py',
        'coordination_policy.py','pqq_microstates.py','site_mechanics.py')]
    result={}
    for p in files:
        destination=output/p.name
        shutil.copyfile(p,destination)
        result[p.name]={'live_source':record(p),'preserved_copy':record(destination)}
    write_new(output/'inventory.json',result)
    return record(output/'inventory.json')


def energy(path):
    text = Path(path).read_text()
    if "ORCA TERMINATED NORMALLY" not in text:
        raise InvalidArtifact(f"abnormal termination: {path}")
    if "SCF CONVERGED" not in text or re.search(r"SCF NOT CONVERGED|SCF failed to converge", text, re.I):
        raise InvalidArtifact(f"SCF nonconvergence: {path}")
    hits = re.findall(r"FINAL SINGLE POINT ENERGY\s+([-+\d.EeDd]+)", text)
    if not hits:
        raise InvalidArtifact(f"missing endpoint energy: {path}")
    result = float(hits[-1].replace("D", "E"))
    if not math.isfinite(result):
        raise InvalidArtifact(f"nonfinite energy: {path}")
    return result


def contrast(e_ca_ha, e_la_ha, aquo_gap_ha=None):
    if not all(math.isfinite(x) for x in (e_ca_ha, e_la_ha)):
        raise InvalidArtifact("nonfinite endpoint")
    r = e_ca_ha - e_la_ha
    return {"R_hartree": r, "R_kcal_mol": r * HA_TO_KCAL,
            "S_kcal_mol": None if aquo_gap_ha is None else (r - aquo_gap_ha) * HA_TO_KCAL,
            "reference_status": "unavailable" if aquo_gap_ha is None else "explicit_reporting_gauge"}


def classify_raw(r_kcal, release, protocol):
    if protocol != release["protocol_id"]:
        return "uncalibrated_protocol"
    c = release["calibration"]
    if c["verdict"] != "CALIBRATABLE":
        return "calibration_unavailable"
    if r_kcal <= c["U_max_Ca_kcal_mol"]:
        return "Ca-supported"
    if r_kcal >= c["L_min_La_kcal_mol"]:
        return "Ln-supported"
    return "indeterminate"


def corrected(baseline, delta_u_ca, delta_u_la):
    if delta_u_ca is None or delta_u_la is None:
        return {"status": "environment_unavailable", "S_env_kcal_mol": None,
                "correction_kcal_mol": None, "decision": "uncalibrated_protocol"}
    if not all(math.isfinite(x) for x in (delta_u_ca, delta_u_la)):
        raise InvalidArtifact("nonfinite transfer energy")
    correction = delta_u_ca - delta_u_la
    return {"status": "computed_descriptor", "correction_kcal_mol": correction,
            "S_env_kcal_mol": None if baseline["S_kcal_mol"] is None else baseline["S_kcal_mol"] + correction,
            "decision": "uncalibrated_protocol"}


def xyz(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2 or int(lines[0]) != len(lines) - 2:
        raise InvalidArtifact(f"invalid XYZ length: {path}")
    rows = []
    for line in lines[2:]:
        f = line.split()
        if len(f) != 4:
            raise InvalidArtifact(f"invalid XYZ row: {path}")
        v = tuple(float(t) for t in f[1:])
        if not all(math.isfinite(t) for t in v):
            raise InvalidArtifact("nonfinite coordinates")
        rows.append((f[0], *v))
    return rows


def paired(la, ca, q_la, q_ca):
    a, b = xyz(la), xyz(ca)
    if not a or len(a) != len(b) or a[0][0] != "La" or b[0][0] != "Ca":
        raise InvalidArtifact("invalid metal pair")
    if a[0][1:] != b[0][1:] or a[1:] != b[1:] or q_la - q_ca != 1:
        raise InvalidArtifact("paired coordinate/charge invariant failed")
    if Path(la).read_text().splitlines()[3:] != Path(ca).read_text().splitlines()[3:]:
        raise InvalidArtifact("nonmetal coordinates not byte identical")
    return {"status": "pass", "atom_count": len(a), "charge_difference": 1}
