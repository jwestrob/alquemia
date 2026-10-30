"""Split each selected exchange dG into D4 dispersion, CPCM dielectric, rest of the electronic energy, and thermal (G - E)."""
import json, re
from pathlib import Path
B = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces")
OLD, NEW = B / "lady_relaxed_exchange_20261024/runs", B / "lady_panel_20261028/runs"
H = 627.509474
def last(pat, txt):
    m = re.findall(pat, txt); return float(m[-1]) if m else None
def comp(n):
    t = ((OLD if (OLD / n).exists() else NEW) / n / "job.out").read_text(errors="ignore")
    E = last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    d = last(r"Dispersion correction\s+(-?\d+\.\d+)\s*\n", t)
    c = last(r"CPCM Dielectric\s+:\s+(-?\d+\.\d+)", t)
    G = last(r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)", t)
    return {"disp": d, "cpcm": c, "elec_rest": E - d - c, "thermal": G - E, "E": E}
R = json.load(open(B / "lady_panel_20261028/PANEL_RESULT.json"))["results"]["conformer_best"]
EXP = {"ac": -0.068, "edtax": -4.010, "nta": -1.582, "dtpa": -4.556, "cdta": -5.674, "tmdta": -4.869}
w = comp("water")
print(f"{'lig':6s} {'exp':>6s} {'dG':>7s} | {'E':>7s} {'disp':>6s} {'cpcm':>7s} {'rest':>7s} {'therm':>6s}   (kcal/mol, exchange; dG also has n*water corr)")
for l, v in R.items():
    b = v["best"]; s = {k: comp(b[k]) for k in ("La_L", "Dy_L", "La_aq", "Dy_aq")}
    dd = {k: H * (s["Dy_L"][k] + s["La_aq"][k] - s["La_L"][k] - s["Dy_aq"][k] + b["n"] * w[k]) for k in w}
    print(f"{l:6s} {EXP[l]:+6.2f} {b['dG']:+7.2f} | {dd['E']:+7.2f} {dd['disp']:+6.2f} {dd['cpcm']:+7.2f} {dd['elec_rest']:+7.2f} {dd['thermal']:+6.2f}   n={b['n']:+d}")
