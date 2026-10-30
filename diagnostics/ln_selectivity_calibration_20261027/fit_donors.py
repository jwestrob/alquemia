"""Empirical La/Dy selectivity: dlogK = log K1(Dy) - log K1(La) (25 C, I <= 0.5) against donor-group
counts from RDKit SMARTS. Leave-one-out CV for: mean baseline, logK_La only, donor counts, donors +
logK_La. Writes fit_report.txt and features.tsv."""
import csv, sys
import numpy as np, pandas as pd
from rdkit import Chem
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import LeaveOneOut, cross_val_predict

KC = 1.36394  # kcal per log unit
SMARTS = {
    "carboxyl":   "[CX3](=O)[OX2H1,OX1-]",
    "aminoacid_N": "[NX3;!$(N-C=O)]-[CH2,CH1]-[CX3](=O)[OX2H1,OX1-]",  # glycinate-type N (5-ring chelate)
    "amine_N":    "[NX3;!$(N-C=[O,S,N]);!$(N-a);!$(N-[N,O]);!$(N=*)]",
    "ether_O":    "[OX2;!H;$(O([#6;!$(C=O)])[#6;!$(C=O)])]",
    "alcohol":    "[OX2H1][CX4]",
    "phenol":     "[OX2H1]c",
    "phosphonate": "P(=O)([OX2H1,OX1-])",
    "pyridine_N": "[nX2]",
    "amide_O":    "[OX1]=[CX3][NX3]",
    "carbonyl_O": "[OX1]=[CX3;!$(C-[OX2H1,OX1-]);!$(C-N)]",
    "sulfonate":  "S(=O)(=O)[OX2H1,OX1-]",
    "hydroxamate": "C(=O)N[OX2H1]",
}
PAT = {k: Chem.MolFromSmarts(v) for k, v in SMARTS.items()}

p = pd.read_csv("pairs.tsv", sep="\t", dtype={"lig_id": str})
s = pd.read_csv("smiles.tsv", sep="\t", dtype={"lig_id": str})
d = p.merge(s[s.formula_ok == "yes"], on="lig_id")
d = d[d.I <= 0.5].copy()
F = []
for smi in d.smiles:
    m = Chem.MolFromSmiles(smi)
    f = {k: len(m.GetSubstructMatches(q)) for k, q in PAT.items()}
    ri = m.GetRingInfo()
    f["macrocycle"] = int(any(len(r) >= 9 for r in ri.AtomRings()))
    f["heavy_atoms"] = m.GetNumHeavyAtoms()
    F.append(f)
F = pd.DataFrame(F, index=d.index)
d = pd.concat([d, F], axis=1)
d.to_csv("features.tsv", sep="\t", index=False)
y = d.dlogK.values
donor_cols = [c for c in F.columns if F[c].sum() > 0 and c != "heavy_atoms"]
sets = {"mean": None, "logK_La": ["logK_La"], "donors": donor_cols, "donors+logK_La": donor_cols + ["logK_La"]}
out = [f"n = {len(d)} ligands (I <= 0.5, PubChem SMILES formula-checked); sd(dlogK) = {y.std():.2f} log = {KC*y.std():.2f} kcal"]
loo = LeaveOneOut()
for name, cols in sets.items():
    if cols is None:
        pred = np.array([np.delete(y, i).mean() for i in range(len(y))])
    else:
        X = d[cols].values.astype(float)
        mdl = RidgeCV(alphas=np.logspace(-3, 2, 30))
        pred = cross_val_predict(mdl, X, y, cv=loo)
    rmse = np.sqrt(np.mean((y - pred) ** 2)); r = np.corrcoef(y, pred)[0, 1]
    out.append(f"{name:16s} LOO RMSE {rmse:.2f} log = {KC*rmse:.2f} kcal  r {r:.2f}")
    if name == "donors":
        mdl = RidgeCV(alphas=np.logspace(-3, 2, 30)).fit(d[cols].values.astype(float), y)
        out.append(f"   alpha {mdl.alpha_:.3g}; intercept {mdl.intercept_:+.3f} log; per-group dlogK (kcal, negative ddG = Dy-favoured):")
        for c, w in sorted(zip(cols, mdl.coef_), key=lambda t: -t[1]):
            out.append(f"     {c:12s} {w:+.3f} log  ({-KC*w:+.2f} kcal)  n_ligands_with={int((d[c]>0).sum())}")
        d["pred_donors_LOO"] = pred
out.append("\nlargest LOO residuals (donor model):")
d["res"] = d.dlogK - d.pred_donors_LOO
for _, r in d.reindex(d.res.abs().sort_values(ascending=False).index).head(12).iterrows():
    out.append(f"   {r['name'][:60]:60s} dlogK {r.dlogK:+.2f} pred {r.pred_donors_LOO:+.2f}")
open("fit_report.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
