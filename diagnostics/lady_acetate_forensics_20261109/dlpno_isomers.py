"""DLPNO-CCSD(T) minus PBE0-D4 (gas, fixed geometry) for same-composition acetate isomers and per-water steps.
Isomer differences (bi_wN vs mono_wN) need no free-water energy (the DLPNO water job failed).
Also prints wall times of single-point gradients usable for costing. Writes dlpno_isomers.json."""
import json, re
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]
T = REPO / "workspaces/lady_ingredient_tests_20261026"
OUT = REPO / "workspaces/lady_acetate_forensics_20261109"
H = 627.509474


def e(test, n):
    p = T / test / n / "job.out"
    if not p.exists(): return None
    t = p.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in t: return None
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])


def off(n):
    a, b = e("dlpno", n), e("gas", n)
    return None if None in (a, b) else H * (a - b)


res = {}
for m in ("La", "Dy"):
    o = {s: off(f"{m}_{s}") for s in ("acmono_w7", "acbi_w7", "acmono_w8", "acbi_w6", "aq8", "aq9")}
    iso = o["acbi_w7"] - o["acmono_w7"]
    steps = {k: (None if None in (o[a], o[b]) else o[a] - o[b]) for k, (a, b) in
             {"mono w7->w8": ("acmono_w8", "acmono_w7"), "bi w6->w7": ("acbi_w7", "acbi_w6"), "aq 8->9": ("aq9", "aq8")}.items()}
    res[m] = {"bi7_minus_mono7": iso, "per_water_step_offsets(includes unknown free-water offset)": steps}
    print(f"{m}: DLPNO-PBE0D4 [bi_w7 - mono_w7] = {iso:+.2f} kcal/mol;  per-water-step offsets "
          + " ".join(f"{k}:{'n/a' if v is None else f'{v:+.2f}'}" for k, v in steps.items()))
(OUT / "dlpno_isomers.json").write_text(json.dumps(res, indent=1))
