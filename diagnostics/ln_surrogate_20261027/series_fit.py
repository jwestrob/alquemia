"""Does MACE capture the *shape* of the Ln series (the part donor counts cannot predict)?

calc:  dE_M = [E(M.L) - E(La.L)] - [E(M aq) - E(La aq)]   (CN8, kcal; M = Nd, Gd, Dy, Ho, Lu)
exp:   dG_M = -1.364 [log K_M - log K_La]
Tests (kill criteria stated up front):
  (a) curvature  dE_Gd - dE_Lu/2  vs  dG_Gd - dG_Lu/2 ; partial r on charge must be > 0.3 with p < 0.05
  (b) LOO ridge for dG_Dy and for exp curvature: donors+charge vs donors+charge+calc shape features
usage: series_fit.py <workspace>
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import LeaveOneOut, cross_val_predict

KC = 1.36394
W = Path(sys.argv[1])
C = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ln_selectivity_calibration_20261027")
MS = ["Nd", "Gd", "Dy", "Ho", "Lu"]


def energies(d):
    e = {}
    for fn in ("energies.json", "series.json"):
        if (d / fn).exists(): e.update(json.loads((d / fn).read_text()))
    return e


aq = energies(W / "nist/aqua")
p = pd.read_csv(C / "pairs.tsv", sep="\t", dtype={"lig_id": str})
f = pd.read_csv(C / "features.tsv", sep="\t", dtype={"lig_id": str})[["lig_id", "carboxyl", "amine_N", "ether_O", "alcohol", "macrocycle"]]
rows = []
for _, r in p.iterrows():
    d = W / "nist" / r.lig_id
    if not d.is_dir(): continue
    e = energies(d); k = dict(x.split(":") for x in r.logK_series.split(";"))
    row = {"lig_id": r.lig_id, "name": r["name"]}
    if "La_cn8" not in e: continue
    row["q_lig"] = 3 - e["La_cn8"]["charge"]; row["CN_La"] = e["La_cn8"]["CN"]
    for M in MS:
        key = f"{M}_cn8"
        if key in e and key in aq:
            row[f"c_{M}"] = (e[key]["E_kcal"] - e["La_cn8"]["E_kcal"]) - (aq[key]["E_kcal"] - aq["La_cn8"]["E_kcal"])
            row[f"CN_{M}"] = e[key]["CN"]
        if M in k: row[f"x_{M}"] = -KC * (float(k[M]) - float(k["La"]))
    rows.append(row)
X = pd.DataFrame(rows).merge(f, on="lig_id", how="left")
print(f"ligands with structures: {len(X)}; aqua series: {[m for m in MS if f'{m}_cn8' in aq]}")
for M in MS:
    if f"c_{M}" in X: print(f"  {M}: calc n={X[f'c_{M}'].notna().sum()}  CN lost vs La: {(X[f'CN_{M}'] < X.CN_La).sum()}")

Z = X.dropna(subset=["c_Gd", "c_Lu", "x_Gd", "x_Lu"]).copy()
Z["curv_c"] = Z.c_Gd - Z.c_Lu / 2; Z["curv_x"] = Z.x_Gd - Z.x_Lu / 2
Z["slope_c"] = Z.c_Lu; Z["slope_x"] = Z.x_Lu


def resid(a, b):
    A = np.c_[np.ones(len(b)), b]; return a - A @ np.linalg.lstsq(A, a, rcond=None)[0]


print(f"\n(a) shape test, n={len(Z)}")
for nm, c, x in (("curvature (Gd - Lu/2)", "curv_c", "curv_x"), ("slope (Lu-La)", "slope_c", "slope_x")):
    r, pv = pearsonr(Z[c], Z[x]); rp, pp = pearsonr(resid(Z[c].values, Z.q_lig.values), resid(Z[x].values, Z.q_lig.values))
    print(f"  {nm:22s} r {r:+.2f} (p {pv:.2g})   partial|charge r {rp:+.2f} (p {pp:.2g})   sd exp {Z[x].std():.2f} kcal  sd calc {Z[c].std():.2f}")

D = ["carboxyl", "amine_N", "ether_O", "alcohol", "macrocycle", "q_lig"]
Z = Z.dropna(subset=D)


def loo(cols, ycol):
    y = Z[ycol].values
    pr = cross_val_predict(RidgeCV(alphas=np.logspace(-3, 3, 40)), Z[cols].values.astype(float), y, cv=LeaveOneOut())
    return np.sqrt(np.mean((y - pr) ** 2)), pr


print(f"\n(b) LOO RMSE kcal, n={len(Z)}")
calc = ["c_Gd", "c_Lu", "curv_c"]
for ycol in ("x_Dy", "curv_x"):
    if Z[ycol].isna().any(): continue
    base, _ = loo(D, ycol); full, pr = loo(D + calc, ycol); only, _ = loo(calc + ["q_lig"], ycol)
    print(f"  {ycol:7s} sd {Z[ycol].std():.2f} | donors+charge {base:.2f} | calc+charge {only:.2f} | donors+charge+calc {full:.2f}")
    if ycol == "x_Dy": Z["pred_Dy"] = pr
for nm in ("EDTA", "Cyclohexylene", "Triethylenetetra", "Trimethylenedinitrilo", "Diethylenetrinitrilo", "Tetraazacyclododecane"):
    for _, r in Z[Z.name.str.contains(nm)].head(1).iterrows():
        print(f"   {r['name'][:44]:44s} ddG_Dy exp {r.x_Dy:+.2f} pred {r.pred_Dy:+.2f}  curv exp {r.curv_x:+.2f} calc {r.curv_c:+.2f}")
Z.to_csv(W / "series_fit.tsv", sep="\t", index=False)
