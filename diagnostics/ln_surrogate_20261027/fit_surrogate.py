"""Physics-shaped surrogate for La/Dy selectivity on the NIST carboxylate/aminocarboxylate set.

Per ligand, from MACE-OMol gas relaxations (workspace nist/<id>/energies.json):
  dE_cn = [E(Dy.L, CN) - E(La.L, CN)] - [E(Dy aq, CN) - E(La aq, CN)]   (same waters on both sides)
Experiment: dlogK = log K(Dy) - log K(La); ddG = -1.364 dlogK.
Models (ridge, LOO):  donor counts (reference), dE9 + charge (+ donors), etc.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.linear_model import RidgeCV, LinearRegression
from sklearn.model_selection import LeaveOneOut, cross_val_predict

KC = 1.36394
W = Path(sys.argv[1] if len(sys.argv) > 1 else
         "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ln_surrogate_20261027")
F = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ln_selectivity_calibration_20261027/features.tsv")

aq = json.loads((W / "nist/aqua/energies.json").read_text())
d = pd.read_csv(F, sep="\t", dtype={"lig_id": str})
rows = []
for _, r in d.iterrows():
    f = W / "nist" / r.lig_id / "energies.json"
    if not f.exists(): continue
    e = json.loads(f.read_text()); row = {"lig_id": r.lig_id}
    for cn in ("cn8", "cn9"):
        k = [f"{m}_{cn}" for m in ("La", "Dy")]
        if all(x in e for x in k) and all(x in aq for x in k):
            row[f"dE_{cn}"] = (e[k[1]]["E_kcal"] - e[k[0]]["E_kcal"]) - (aq[k[1]]["E_kcal"] - aq[k[0]]["E_kcal"])
            row[f"nlig_La_{cn}"] = e[k[0]]["n_ligand"]; row[f"nlig_Dy_{cn}"] = e[k[1]]["n_ligand"]
            row[f"conv_{cn}"] = e[k[0]]["converged"] and e[k[1]]["converged"]
            row["charge"] = e[k[0]]["charge"]
    rows.append(row)
X = d.merge(pd.DataFrame(rows), on="lig_id")
X = X[X.I <= 0.5].dropna(subset=["dE_cn8"]).copy()
X["q_lig"] = 3 - X.charge
y = X.dlogK.values
loo = LeaveOneOut()

def cv(cols, label):
    Z = X[cols].values.astype(float)
    pred = cross_val_predict(RidgeCV(alphas=np.logspace(-4, 3, 40)), Z, y, cv=loo)
    rmse = np.sqrt(np.mean((y - pred) ** 2))
    print(f"{label:40s} n={len(y)} LOO RMSE {KC*rmse:.2f} kcal  r {np.corrcoef(y, pred)[0,1]:+.2f}")
    return pred

print(f"sd(ddG) = {KC*y.std():.2f} kcal")
donor_cols = ["carboxyl", "amine_N", "ether_O", "alcohol", "macrocycle"]
cv(donor_cols, "donor counts")
cv(["q_lig"], "ligand charge")
p1 = cv(["dE_cn8"], "MACE dE_cn8")
p2 = cv(["dE_cn8", "q_lig"], "MACE dE_cn8 + charge")
p3 = cv(["dE_cn8", "q_lig", "nlig_La_cn8"], "MACE dE_cn8 + charge + denticity")
if "dE_cn8" in X and X.dE_cn8.notna().all():
    pass
cv(donor_cols + ["dE_cn8", "q_lig"], "donors + dE_cn8 + charge")
m = LinearRegression().fit(X[["dE_cn8", "q_lig", "nlig_La_cn8"]], y)
print("fit dlogK =", " + ".join(f"{c:.3f}*{n}" for c, n in zip(m.coef_, ["dE9", "q", "dent"])), f"+ {m.intercept_:.3f}")
X["pred"] = p3
X["res_kcal"] = -KC * (y - p3)
for nm in ("EDTA", "CDTA", "Cyclohexylene", "TTHA", "Triethylenetetra", "TMDTA", "Trimethylenedinitrilo", "DTPA", "NTA", "IDA", "Acetic", "Malonic", "DOTA"):
    for _, r in X[X.name.str.contains(nm)].head(1).iterrows():
        print(f"   {r['name'][:50]:50s} dlogK {r.dlogK:+.2f} pred {r.pred:+.2f}  dE9 {r.dE_cn8:+7.2f} q {r.q_lig} dent {r.nlig_La_cn8}")
X.to_csv(W / "surrogate_nist.tsv", sep="\t", index=False)
