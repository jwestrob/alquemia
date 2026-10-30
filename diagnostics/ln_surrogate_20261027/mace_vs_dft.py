"""MACE-OMol-0 (gas, omega-B97M-V/def2-TZVPD-trained) against gas-phase PBE0-D4 single points on the
La/Dy acetate exchange clusters (geometries: lady_ingredient_tests gas/<name>/start.xyz).

Compares exchange energies dE = E(Dy.L) + E(La aq) - E(La.L) - E(Dy aq) + n E(water) for the
combinations used in collect.py, against large-core PBE0 (gas) and 4f-explicit PBE0 (sc_gas) where
available. Charges from composition; multiplicity La 1, Dy 6 (4f9).
"""
import re, sys, json
from pathlib import Path
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

torch.set_num_threads(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
T = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lady_ingredient_tests_20261026")
CKPT = "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/software_v1/MACE-omol-0-extra-large-1024.model"
H = 627.509474; EV = 23.060548

def charge(n):
    if n == "water": return 0
    return 3 - (1 if "_ac" in n else 0)

def mult(n):
    return 6 if n.startswith("Dy") else 1

def dft(test, n):
    p = T / test / n / "job.out"
    if not p.exists(): return None
    t = p.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in t: return None
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    return float(m[-1]) * H if m else None

calc = MACECalculator(model_paths=CKPT, device="cuda" if torch.cuda.is_available() else "cpu", default_dtype="float64", head="omol")
names = [p.name for p in sorted((T / "gas").iterdir()) if (p / "start.xyz").exists() and not p.name.endswith("_ion")]
E = {"mace": {}, "gas": {}, "sc_gas": {}}
for n in names:
    a = read(T / "gas" / n / "start.xyz")
    a.info["charge"] = charge(n); a.info["spin"] = mult(n)
    a.calc = calc
    E["mace"][n] = a.get_potential_energy() * EV
    E["gas"][n] = dft("gas", n); E["sc_gas"][n] = dft("sc_gas", n)
E["sc_gas"]["water"] = E["gas"]["water"]
nw = lambda s: int(re.search(r"(?:aq|_w)(\d+)$", s).group(1))
combos = {
    "acmono best": ("La_acmono_w8", "Dy_acmono_w7", "La_aq9", "Dy_aq8"),
    "acmono CN9": ("La_acmono_w8", "Dy_acmono_w8", "La_aq9", "Dy_aq9"),
    "acmono CN8": ("La_acmono_w7", "Dy_acmono_w7", "La_aq8", "Dy_aq8"),
    "acbi best": ("La_acbi_w7", "Dy_acbi_w7", "La_aq9", "Dy_aq8"),
    "acbi CN9": ("La_acbi_w7", "Dy_acbi_w7", "La_aq9", "Dy_aq9"),
    "acbi CN8": ("La_acbi_w6", "Dy_acbi_w6", "La_aq8", "Dy_aq8"),
    "merged best": ("La_acmono_w8", "Dy_acbi_w7", "La_aq9", "Dy_aq8"),
    # pure hydration-number steps (same metal): M(aq9) -> M(aq8) + water
    "La aq9->aq8+w": None, "Dy aq9->aq8+w": None,
}
out = {}
for k, c in combos.items():
    row = {}
    for m, Em in E.items():
        if c is None:
            M = k[:2]; v = [Em.get(f"{M}_aq8"), Em.get("water"), Em.get(f"{M}_aq9")]
            row[m] = None if None in v else v[0] + v[1] - v[2]
        else:
            laL, dyL, laA, dyA = c
            n = nw(laL) + nw(dyA) - nw(dyL) - nw(laA)
            v = [Em.get(x) for x in (dyL, laA, laL, dyA, "water")]
            row[m] = None if None in v else v[0] + v[1] - v[2] - v[3] + n * v[4]
    out[k] = row
    print(f"{k:16s} " + "  ".join(f"{m}={'   n/a' if x is None else f'{x:+7.2f}'}" for m, x in row.items())
          + (f"   mace-gas={row['mace']-row['gas']:+.2f}" if row['gas'] is not None else ""))
Path(__file__).with_name("MACE_VS_DFT.json").write_text(json.dumps(out, indent=1))
