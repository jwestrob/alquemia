"""Data audit behind the proposed NIST La-Lu predictor (2026-10-29; LANTHANIDE_NIST_PREDICTOR_REVIEW_2026-10-29.md).

Reads the pipeline outputs in workspaces/ln_selectivity_calibration_20261027/ (ln_ML.tsv, pairs.tsv, smiles.tsv,
features.tsv; see README.md) and writes DATA_AUDIT_20261029.txt next to this script. Single thread, a few seconds.
  1. How many ligands each data source would add.
  2. Scatter of dlogK = log K(Dy) - log K(La) between ionic strengths for the same ligand (the noise floor).
  3. How concentrated the donor-count model's leave-one-out error is.
  4. The diamine tetraacetates (same donors, different backbone), proposed as the held-out block.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import LeaveOneOut, cross_val_predict

KC = 1.36394  # kcal/mol per log unit at 298.15 K
HERE = Path(__file__).resolve().parent
W = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parents[1] / "workspaces/ln_selectivity_calibration_20261027"
LN = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
rms = lambda a: float(np.sqrt(np.mean(np.square(a))))
out = []
say = out.append

ml = pd.read_csv(W / "ln_ML.tsv", sep="\t", dtype={"lig_id": str})
k = ml[(ml.qty == "K") & (ml["T"].astype(str) == "25") & ml.metal.isin(LN)]
pairs = pd.read_csv(W / "pairs.tsv", sep="\t", dtype={"lig_id": str})
smi = pd.read_csv(W / "smiles.tsv", sep="\t", dtype={"lig_id": str})
p = pairs.merge(smi[["lig_id", "formula_ok"]], on="lig_id", how="left")
ok = p.formula_ok == "yes"

# 1. Sources. A ligand without La or Dy still fixes the one-shape amplitude if a light and a heavy Ln share conditions.
w = k.pivot_table(index=["lig_id", "I"], columns="metal", values="value", aggfunc="first")
light, heavy = set(LN[:4]), set(LN[7:])
shape_ok = {L for (L, _), row in w.iterrows() if set(row.dropna().index) & light and set(row.dropna().index) & heavy}
say("1. Ligands by source (25 C, 1:1 complexes, NIST SRD 46 v8)")
say(f"   La and Dy at the same ionic strength:            {len(p)}")
say(f"     in the current model (SMILES, I <= 0.5):       {(ok & (p.I <= 0.5)).sum()}")
say(f"     SMILES matched but I > 0.5:                    {(ok & (p.I > 0.5)).sum()}")
say(f"     PubChem lookup failed (no structure yet):      {(~ok).sum()}")
say(f"   a light (La-Nd) and a heavy (Tb-Lu) Ln, same I:  {len(shape_ok)} (+{len(shape_ok - set(p.lig_id))} beyond the La/Dy pairs)")
say(f"   any Ln log K:                                    {k.lig_id.nunique()}")

# 2. Ionic-strength replicates of dlogK
d = k[k.metal.isin(["La", "Dy"])].pivot_table(index=["lig_id", "I"], columns="metal", values="value", aggfunc="first").dropna()
d["dl"] = d.Dy - d.La
rep = [(L, x) for L, x in d.reset_index().groupby("lig_id") if len(x) > 1]
shift = np.array([x.dl.max() - x.dl.min() for _, x in rep])
sd = np.array([x.dl.std(ddof=1) for _, x in rep])
strong = np.array([x.La.max() > 10 for _, x in rep])
say("\n2. dlogK(Dy-La) for ligands listed at more than one ionic strength")
say(f"   {len(rep)} ligands; median shift {np.median(shift):.2f} log, max {shift.max():.2f} log")
say(f"   per-value sd: all {rms(sd):.2f} log = {KC * rms(sd):.2f} kcal;"
    f" log K_La > 10 ({strong.sum()}) {rms(sd[strong]):.2f} log = {KC * rms(sd[strong]):.2f} kcal")
say("   This mixes real ionic-strength dependence with between-study scatter: read it as the floor for pooled data.")

# 3. Donor-count model: same features and ridge as fit_donors.py
f = pd.read_csv(W / "features.tsv", sep="\t", dtype={"lig_id": str})
cols = [c for c in ["carboxyl", "aminoacid_N", "amine_N", "ether_O", "alcohol", "phenol", "phosphonate", "pyridine_N",
                    "amide_O", "carbonyl_O", "sulfonate", "hydroxamate", "macrocycle"] if f[c].sum() > 0]
y = f.dlogK.values
pred = cross_val_predict(RidgeCV(alphas=np.logspace(-3, 2, 30)), f[cols].values.astype(float), y, cv=LeaveOneOut())
res = KC * (y - pred); se = res ** 2; top = np.argsort(-se)[:7]
say("\n3. Donor-count model (fit_donors.py features), leave-one-out")
say(f"   n {len(y)}, sd {KC * y.std():.2f} kcal, RMSE {rms(res):.2f} kcal; the 7 largest residuals carry"
    f" {se[top].sum() / se.sum():.0%} of the squared error; RMSE on the other {len(y) - 7}: {rms(np.delete(res, top)):.2f} kcal")
for i in top:
    say(f"   {res[i]:+.2f} kcal  {f.name.iloc[i][:70]}")
say("   (positive = more Dy-selective than predicted)")

# 4. Diamine-N,N,N',N'-tetraacetates: the same donor set on different backbones
lp = dict(zip(f.lig_id, pred))
fam = p[p.name.str.contains("dinitrilotetraacetic acid", regex=False)].sort_values("dlogK")
say(f"\n4. Diamine-N,N,N',N'-tetraacetates among the La/Dy pairs ({len(fam)}): the proposed held-out block")
say("   lig_id    I  logK_La  logK_Dy  dlogK  donor-model LOO  name")
for _, r in fam.iterrows():
    pr = f"{lp[r.lig_id]:+.2f}" if r.lig_id in lp else "not in model"
    say(f"   {r.lig_id:>6} {r.I:4.1f}  {r.logK_La:7.2f}  {r.logK_Dy:7.2f}  {r.dlogK:+5.2f}  {pr:>15}  {r['name'][:60]}")
span = fam.dlogK.max() - fam.dlogK.min()
say(f"   spread {span:.2f} log = {KC * span:.2f} kcal; {fam.lig_id.isin(lp).sum()} of {len(fam)} are in the current model")
cdta = d.loc["600"].dl if "600" in d.index.get_level_values(0) else None
if cdta is not None:
    say("   CDTA (600) at every ionic strength: " + ", ".join(f"I={i:g} {v:+.2f}" for i, v in cdta.items()))

txt = "\n".join(out) + "\n"
(HERE / "DATA_AUDIT_20261029.txt").write_text(txt)
print(txt)
