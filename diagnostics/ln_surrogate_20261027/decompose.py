"""Option 1: does MACE gas dE_cn8 track ddH (bonding) rather than ddG?  TddS vs charge.
Option 2: low-rank structure of the full Ln series (log K_M - log K_La over 13 metals).
usage: decompose.py <surrogate_nist.tsv (with dE_cn8)>"""
import sys, numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import LeaveOneOut, cross_val_predict
KC = 1.36394
C = "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ln_selectivity_calibration_20261027/"
p = pd.read_csv(C + "pairs.tsv", sep="\t", dtype={"lig_id": str})
f = pd.read_csv(C + "features.tsv", sep="\t", dtype={"lig_id": str})
s = pd.read_csv(sys.argv[1], sep="\t", dtype={"lig_id": str}).drop_duplicates("lig_id")

print("== Option 1: enthalpy/entropy split ==")
h = p.dropna(subset=["ddH", "TddS"]).merge(s[["lig_id", "dE_cn8", "q_lig", "nlig_La_cn8"]], on="lig_id", how="left")
print(f"ligands with ddH & TddS: {len(h)};  with MACE dE too: {h.dE_cn8.notna().sum()}")
print(f"  corr(ddH, TddS) r {pearsonr(h.ddH, h.TddS)[0]:+.2f}  (enthalpy-entropy compensation)")
hm = h.dropna(subset=["dE_cn8"])
for col in ("ddG", "ddH", "TddS"):
    r, pv = pearsonr(hm.dE_cn8, hm[col]); rs = spearmanr(hm.dE_cn8, hm[col])[0]
    rq = pearsonr(hm.q_lig, hm[col])[0]
    print(f"  {col:5s} vs dE_cn8 r {r:+.2f} (p {pv:.2g}, rho {rs:+.2f})   vs charge r {rq:+.2f}   n={len(hm)}")
# partial: dE vs ddH controlling for charge
def resid(a, b):
    A = np.c_[np.ones(len(b)), b]; return a - A @ np.linalg.lstsq(A, a, rcond=None)[0]
print(f"  partial r(dE, ddH | charge) {pearsonr(resid(hm.dE_cn8.values, hm.q_lig.values), resid(hm.ddH.values, hm.q_lig.values))[0]:+.2f}")
print(f"  partial r(dE, TddS| charge) {pearsonr(resid(hm.dE_cn8.values, hm.q_lig.values), resid(hm.TddS.values, hm.q_lig.values))[0]:+.2f}")
for _, r in hm.sort_values("ddG").iterrows():
    print(f"    {r['name'][:42]:42s} ddG {r.ddG:+5.2f} ddH {r.ddH:+5.2f} TddS {r.TddS:+5.2f} dE {r.dE_cn8:+6.1f} q {r.q_lig:.0f}")

print("\n== Option 2: lanthanide-series structure ==")
M = ["La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"]
R = dict(zip(M, [1.216,1.196,1.179,1.163,1.132,1.120,1.107,1.095,1.083,1.072,1.062,1.052,1.042,1.032]))  # Shannon CN9
rows = []
for _, r in p.iterrows():
    k = dict(x.split(":") for x in r.logK_series.split(";"))
    rows.append({m: float(k[m]) - float(k["La"]) if m in k else np.nan for m in M} | {"lig_id": r.lig_id})
S = pd.DataFrame(rows).set_index("lig_id")[M]
full = S.dropna()
print(f"ligands: {len(S)}; complete La..Lu series: {len(full)}")
U, sv, Vt = np.linalg.svd(full.values - full.values.mean(0), full_matrices=False)
ev = sv**2 / (sv**2).sum()
print("  variance explained by components:", " ".join(f"{x:.3f}" for x in ev[:5]))
print("  PC1 loading:", " ".join(f"{m}:{v:+.2f}" for m, v in zip(M, Vt[0])))
print("  PC2 loading:", " ".join(f"{m}:{v:+.2f}" for m, v in zip(M, Vt[1])))
# predict Dy-La from a cheap probe: one other metal (each), LOO linear
y = full["Dy"].values
print("  predict dlogK(Dy-La) from one measured metal (LOO RMSE kcal):")
for m in ["Ce","Nd","Sm","Eu","Gd","Lu"]:
    pr = cross_val_predict(RidgeCV(alphas=np.logspace(-4,2,20)), full[[m]].values, y, cv=LeaveOneOut())
    print(f"    from {m}-La: {KC*np.sqrt(np.mean((y-pr)**2)):.2f}")
print(f"    sd baseline: {KC*y.std():.2f}")
# Gd break: curvature = (Gd-La)/(Lu-La) fraction vs linear-in-radius expectation
full2 = full.copy()
full2["early_slope"] = full["Nd"] / (R["La"] - R["Nd"]); full2["late_slope"] = (full["Lu"] - full["Gd"]) / (R["Gd"] - R["Lu"])
print(f"  median slope (log K per A of contraction) La->Nd {full2.early_slope.median():.1f}, Gd->Lu {full2.late_slope.median():.1f}")
print(f"  fraction of ligands whose selectivity reverses or plateaus after Gd (late slope <= 0): {(full2.late_slope<=0).mean():.2f}")
# Do donor features predict the PC scores better than they predict dlogK directly?
sc = pd.DataFrame(U[:, :2] * sv[:2], index=full.index, columns=["PC1", "PC2"]).reset_index()
fm = f.merge(sc, on="lig_id")
dc = ["carboxyl","amine_N","ether_O","alcohol","macrocycle"]
for t in ("PC1", "PC2"):
    yy = fm[t].values
    pr = cross_val_predict(RidgeCV(alphas=np.logspace(-3,2,30)), fm[dc].values.astype(float), yy, cv=LeaveOneOut())
    print(f"  donors -> {t}: LOO R^2 {1-np.mean((yy-pr)**2)/yy.var():.2f} (n={len(fm)})")
