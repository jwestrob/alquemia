"""CPCM-DFT La/Dy exchange on MACE CN8 structures vs NIST.  dE = [E(Dy.L)-E(La.L)] - [E(Dy aq)-E(La aq)].
Kill: dE(+charge) must beat donors+charge LOO on the same ligands and order CDTA < EDTA < TTHA.
usage: cpcm_fit.py <ws>"""
import re, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import LeaveOneOut, cross_val_predict
KC, H = 1.36394, 627.509474
W = Path(sys.argv[1]); C = W.parent / "ln_selectivity_calibration_20261027"


def E(t):
    p = W / "cpcm" / t / "job.out"
    if not p.exists(): return None, None
    s = p.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in s: return None, None
    tot = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", s)[-1]) * H
    cp = re.findall(r"CPCM Dielectric\s*:\s*(-?\d+\.\d+)", s)
    return tot, (float(cp[-1]) * H if cp else None)


aqL, aqLc = E("aqua_La"); aqD, aqDc = E("aqua_Dy")
p = pd.read_csv(C / "pairs.tsv", sep="\t", dtype={"lig_id": str})
f = pd.read_csv(C / "features.tsv", sep="\t", dtype={"lig_id": str})[["lig_id", "carboxyl", "amine_N", "ether_O", "alcohol", "macrocycle"]]
rows = []
for _, r in p.iterrows():
    (l, lc), (d, dc) = E(f"{r.lig_id}_La"), E(f"{r.lig_id}_Dy")
    if l is None or d is None: continue
    q = int(re.search(r"\* xyzfile (-?\d+)", (W / "cpcm" / f"{r.lig_id}_La" / "job.inp").read_text()).group(1))
    rows.append({"lig_id": r.lig_id, "name": r["name"], "ddG": r.ddG, "ddH": r.ddH, "TddS": r.TddS, "q_lig": 3 - q,
                 "dE_cpcm": (d - l) - (aqD - aqL),
                 "dE_diel": None if None in (dc, lc, aqDc, aqLc) else (dc - lc) - (aqDc - aqLc)})
X = pd.DataFrame(rows).merge(f, on="lig_id", how="left").dropna(subset=["carboxyl"])
X["dE_nodiel"] = X.dE_cpcm - X.dE_diel
print(f"n={len(X)}  sd ddG {X.ddG.std():.2f}")
for c in ("dE_cpcm", "dE_nodiel", "dE_diel"):
    print(f"  r({c:9s}, ddG) {pearsonr(X[c], X.ddG)[0]:+.2f}   range {X[c].min():+.1f}..{X[c].max():+.1f}")
h = X.dropna(subset=["ddH"])
print(f"  n={len(h)} with ddH:  r(dE_cpcm, ddH) {pearsonr(h.dE_cpcm, h.ddH)[0]:+.2f}  r(dE_cpcm, TddS) {pearsonr(h.dE_cpcm, h.TddS)[0]:+.2f}")
y = X.ddG.values
D = ["carboxyl", "amine_N", "ether_O", "alcohol", "macrocycle", "q_lig"]
for lab, cols in (("donors+charge", D), ("dE_cpcm", ["dE_cpcm"]), ("dE_cpcm+charge", ["dE_cpcm", "q_lig"]), ("donors+charge+dE_cpcm", D + ["dE_cpcm"])):
    pr = cross_val_predict(RidgeCV(alphas=np.logspace(-3, 3, 40)), X[cols].values.astype(float), y, cv=LeaveOneOut())
    print(f"  {lab:24s} LOO RMSE {np.sqrt(np.mean((y-pr)**2)):.2f} kcal  r {np.corrcoef(y, pr)[0,1]:+.2f}")
    X[f"pred_{lab}"] = pr
for nm in ("Ethylenedinitrilotetraacetic", "Cyclohexylene", "Triethylenetetra", "Trimethylenedinitrilo", "Diethylenetrinitrilo", "Tetraazacyclododecane-N", "Nitrilotri", "Ethanoic"):
    for _, r in X[X.name.str.contains(nm)].head(1).iterrows():
        print(f"   {r['name'][:40]:40s} ddG {r.ddG:+.2f}  dE_cpcm {r.dE_cpcm:+6.2f}  pred {r['pred_dE_cpcm+charge']:+.2f}")
X.to_csv(W / "cpcm_fit.tsv", sep="\t", index=False)
