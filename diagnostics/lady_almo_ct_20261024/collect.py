"""Balanced La/Dy water->acetate double differences from the four ct_run outputs.

C_X = [X(Dy,mix) - X(La,mix)] - [X(Dy,aq) - X(La,aq)], kcal/mol; negative = Dy favoured.
"""
import json, sys
from pathlib import Path

H = 627.509474
RUN = Path(sys.argv[1])
SYS = {("La", "aq"): "aquo8_central_La", ("Dy", "aq"): "aquo8_central_Dy",
       ("La", "mix"): "hydrated_acetate_mono_central_La", ("Dy", "mix"): "hydrated_acetate_mono_central_Dy"}
r = {k: json.loads((RUN / v / "RESULT.json").read_text()) for k, v in SYS.items()}
bad = {SYS[k]: v["status"] for k, v in r.items() if v["status"] != "passed"}
E = {k: v["energies_hartree"] for k, v in r.items()}


def dd(f):
    return H * ((f(E["Dy", "mix"]) - f(E["La", "mix"])) - (f(E["Dy", "aq"]) - f(E["La", "aq"])))


rows = {
    "full_electronic_U": dd(lambda e: e["U"]),
    "full_U_plus_archived_D4": dd(lambda e: e["U"] + e["D4_archived"]),
    "CT_suppressed_P_plus_archived_D4": dd(lambda e: e["P"] + e["D4_archived"]),
    "frozen_F_plus_archived_D4": dd(lambda e: e["F"] + e["D4_archived"]),
    "FRZ": dd(lambda e: e["FRZ"]), "POL": dd(lambda e: e["POL"]), "CT": dd(lambda e: e["CT"]),
    "archive_ORCA_electronic": dd(lambda e: e["archive_electronic"]),
    "backend_offset_PySCF_minus_ORCA": dd(lambda e: e["U"] - e["archive_electronic"]),
}
per = {SYS[k]: {x: H * E[k][x] for x in ("FRZ", "POL", "CT")} for k in SYS}
out = {"not_passed": bad, "double_differences_kcal_mol": rows, "per_endpoint_kcal_mol": per,
       "solver": {SYS[k]: v["solver"] for k, v in r.items()}}
(RUN / "QUARTET_RESULT.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
