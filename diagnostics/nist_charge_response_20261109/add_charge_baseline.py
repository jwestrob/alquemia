"""ADDITION OUTSIDE THE FROZEN DESIGN (lead request received after PLAN.md was frozen at ~01:43 UTC 2026-11-10
and after the frozen fit had run). Specified here at ~01:48 UTC 2026-11-10, before running:

- Mq: one feature, ligand_charge (net formal charge of the fully deprotonated L; the M1 column), via the imported
  30 Oct nested() with the frozen alpha grid, for every target scored under PLAN.md (same primary sets).
- Paired RMSE(Mq) - RMSE(M1) and RMSE(M0) - RMSE(Mq), bootstrap over ligands and families (2,000, seed 42),
  identical resampling to fit_charge_response.py.
- y1 descriptive: raw-unit coefficient per unit of ligand charge (alpha by leave-one-family-out on the whole primary
  set, ligand and family bootstrap at that alpha); family-held-out Mq residual (obs - pred) vs carboxyl count and vs
  neutral donor count = amine_N + ether_O + alcohol + pyridine_N + amide_O + carbonyl_O (phenol, carboxyl,
  phosphonate, hydroxamate counted as anionic as in D_anO; aminoacid_N omitted as a re-count of amine_N; sulfonate
  omitted). Spearman rho with ligand and family bootstrap intervals.
"""
import hashlib, json, sys
sys.dont_write_bytecode = True
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
D30 = REPO / "diagnostics/ln_structure_increment_20261030"
W30 = REPO / "workspaces/ln_structure_increment_20261030"
OUT = REPO / "workspaces/nist_charge_response_20261109"
assert hashlib.sha256((HERE / "PLAN.md").read_bytes()).hexdigest() == \
    "c2c04a1f98dc55a0dc6519bf2a6cef99e1f390a273eb1dfccfbe8bab7d609246"
sys.path.insert(0, str(D30))
import evaluate  # noqa: E402
from evaluate import nested  # noqa: E402
from structure import M1  # noqa: E402
FROZEN = [.01, .1, 1., 10., 100.]
assert evaluate.ALPHAS == FROZEN
NB, SEED = 2000, 42
NEUTRAL = ["amine_N", "ether_O", "alcohol", "pyridine_N", "amide_O", "carbonyl_O"]

prim = pd.read_csv(W30 / "primary.tsv", sep="\t", dtype={"lig_id": str})
feat = pd.read_csv(W30 / "features.tsv", sep="\t", dtype={"lig_id": str})
base = prim[["lig_id", "name", "family", "normalized_identity", "ligand_charge"]].merge(
    feat.drop(columns=["family", "ligand_charge"]), on="lig_id")
base["neutral_donors"] = base[NEUTRAL].sum(1)
P = pd.read_csv(OUT / "la_pairs_all_metals.tsv", sep="\t", dtype={"lig_id": str})
TARGETS = {"y1": ("Ca/2+", "diff"), "y2": ("Ca/2+", "La"), "y3": ("Ca/2+", "X"), "y7": ("Mg/2+", "diff"),
           "y8": ("Ba/2+", "diff"), "y9": ("Lu/3+", "diff"), "y10": ("Cd/2+", "diff"), "y11": ("Y/3+", "diff")}


def dataset(metal, kind):
    s = P[(P.metal == metal) & P.has_features].copy()
    s = s[(s.match == "same_I") | (s.dI.abs() <= 0.1 + 1e-9)]
    s["y_kcal"] = {"diff": s.logK_La - s.logK_X, "La": s.logK_La, "X": s.logK_X}[kind]
    dd = s[["lig_id", "y_kcal", "match", "dI"]].merge(base, on="lig_id", validate="one_to_one")
    return dd.sort_values("lig_id").reset_index(drop=True)


def boot_idx(dd):
    n = len(dd); fam = dd.family.values; F = sorted(set(fam))
    BL = np.random.default_rng(SEED).integers(0, n, (NB, n))
    rf = np.random.default_rng(SEED); fi = {f: np.flatnonzero(fam == f) for f in F}
    BF = [np.concatenate([fi[g] for g in rf.choice(F, len(F), replace=True)]) for _ in range(NB)]
    return BL, BF


def ci(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    return [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]


def rm(p, y, i):
    return np.sqrt(np.mean((p[i] - y[i]) ** 2))


def r2(p, y, i):
    yy = y[i]; return 1 - ((p[i] - yy) ** 2).sum() / ((yy - yy.mean()) ** 2).sum()


RES = {}
for t, (metal, kind) in TARGETS.items():
    dd = dataset(metal, kind); y = dd.y_kcal.values
    pr = {}; al = {}
    for k, cols in (("M0", None), ("Mq", ["ligand_charge"]), ("M1", M1)):
        pr[k], meta = nested(dd, cols); al[k] = [m["alpha"] for m in meta]
    BL, BF = boot_idx(dd); a = np.arange(len(y))
    out = dict(n=len(dd), n_fam=int(dd.family.nunique()), alphas_Mq=al["Mq"])
    for k in ("Mq", "M1"):
        out[k] = dict(rmse=float(rm(pr[k], y, a)), r2=float(r2(pr[k], y, a)),
                      r2_ci_lig=ci([r2(pr[k], y, b) for b in BL]), r2_ci_fam=ci([r2(pr[k], y, b) for b in BF]),
                      rmse_ci_lig=ci([rm(pr[k], y, b) for b in BL]), rmse_ci_fam=ci([rm(pr[k], y, b) for b in BF]))
    for x, z in (("Mq", "M1"), ("M0", "Mq")):
        out[f"{x}-{z}"] = dict(point=float(rm(pr[x], y, a) - rm(pr[z], y, a)),
                               ci_lig=ci([rm(pr[x], y, b) - rm(pr[z], y, b) for b in BL]),
                               ci_fam=ci([rm(pr[x], y, b) - rm(pr[z], y, b) for b in BF]))
    RES[t] = out
    q = out["Mq"]; d1 = out["Mq-M1"]; d0 = out["M0-Mq"]
    print(f"{t:4s} n={out['n']:3d} fam={out['n_fam']:2d} | Mq RMSE {q['rmse']:.3f} R2 {q['r2']:+.3f} "
          f"lig[{q['r2_ci_lig'][0]:+.2f},{q['r2_ci_lig'][1]:+.2f}] fam[{q['r2_ci_fam'][0]:+.2f},{q['r2_ci_fam'][1]:+.2f}] "
          f"| M1 RMSE {out['M1']['rmse']:.3f} | Mq-M1 {d1['point']:+.3f} lig[{d1['ci_lig'][0]:+.3f},{d1['ci_lig'][1]:+.3f}] "
          f"fam[{d1['ci_fam'][0]:+.3f},{d1['ci_fam'][1]:+.3f}] | M0-Mq {d0['point']:+.3f} "
          f"lig[{d0['ci_lig'][0]:+.3f},{d0['ci_lig'][1]:+.3f}] fam[{d0['ci_fam'][0]:+.3f},{d0['ci_fam'][1]:+.3f}] "
          f"| Mq alphas {al['Mq']}")
    if t == "y1":
        Y1 = (dd, pr)

# ---- y1 descriptive ----
dd, pr = Y1
y = dd.y_kcal.values; X = dd[["ligand_charge"]].values.astype(float); fam = dd.family.values


def fit_coef(Xa, ya, a):
    s = StandardScaler().fit(Xa); m = Ridge(alpha=a).fit(s.transform(Xa), ya)
    return m.coef_[0] / s.scale_[0]


sc = []
for a in FROZEN:
    e = np.full(len(y), np.nan)
    for f in sorted(set(fam)):
        te = fam == f; s = StandardScaler().fit(X[~te])
        e[te] = Ridge(alpha=a).fit(s.transform(X[~te]), y[~te]).predict(s.transform(X[te])) - y[te]
    sc.append(np.mean(e ** 2))
a_full = FROZEN[int(np.argmin(sc))]
BL, BF = boot_idx(dd)
c0 = fit_coef(X, y, a_full)
cl = ci([fit_coef(X[b], y[b], a_full) for b in BL]); cf = ci([fit_coef(X[b], y[b], a_full) for b in BF])
ols = np.polyfit(X[:, 0], y, 1)
print(f"\ny1 per unit ligand charge (ridge alpha {a_full}): {c0:+.3f} log units  lig {cl}  fam {cf}  "
      f"(OLS slope {ols[0]:+.3f}, intercept {ols[1]:+.3f}); negative = more La-selective per extra negative charge")
print("y1 by ligand charge:", dd.groupby("ligand_charge").y_kcal.agg(["size", "median", "min", "max"]).round(2).to_dict("index"))
res = y - pr["Mq"]
dd["res_Mq"] = res
for col in ("carboxyl", "neutral_donors"):
    rho = spearmanr(dd[col], res)[0]
    rl = ci([spearmanr(dd[col].values[b], res[b])[0] for b in BL]); rf = ci([spearmanr(dd[col].values[b], res[b])[0] for b in BF])
    RES[f"y1_res_vs_{col}"] = dict(rho=rho, ci_lig=rl, ci_fam=rf)
    print(f"Spearman(held-out Mq residual, {col}) {rho:+.3f} lig {rl} fam {rf}")
    print("   residual by", col, dd.groupby(col).res_Mq.agg(["size", "median", "min", "max"]).round(2).to_dict("index"))
RES["y1_charge_coef"] = dict(alpha=a_full, coef=c0, ci_lig=cl, ci_fam=cf, ols=ols.tolist())
dd[["lig_id", "name", "family", "ligand_charge", "carboxyl", "neutral_donors", "y_kcal", "res_Mq"]].to_csv(
    OUT / "y1_Mq_residuals.tsv", sep="\t", index=False)
json.dump(RES, open(OUT / "addition_Mq.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
