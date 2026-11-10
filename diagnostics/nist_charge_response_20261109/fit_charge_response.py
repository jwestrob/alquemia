"""Frozen design: PLAN.md (same dir; sha256 checked). Login node, single thread, seconds."""
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
PLAN_SHA = "c2c04a1f98dc55a0dc6519bf2a6cef99e1f390a273eb1dfccfbe8bab7d609246"
assert hashlib.sha256((HERE / "PLAN.md").read_bytes()).hexdigest() == PLAN_SHA, "PLAN changed"
sys.path.insert(0, str(D30))
import evaluate  # noqa: E402  (code only; main() not run)
from evaluate import nested  # noqa: E402
from structure import M1  # noqa: E402

FROZEN = [.01, .1, 1., 10., 100.]
assert evaluate.ALPHAS == FROZEN
NB, SEED = 2000, 42
AN_O = ["carboxyl", "phosphonate", "phenol", "hydroxamate"]
TOT = ["carboxyl", "amine_N", "ether_O", "alcohol", "phenol", "phosphonate", "pyridine_N", "amide_O", "carbonyl_O",
       "hydroxamate"]
GROUPS = [c for c in M1 if c not in ("macrocycle", "ligand_charge")]

# ---------------- Step 0 positive control ----------------
prim = pd.read_csv(W30 / "primary.tsv", sep="\t", dtype={"lig_id": str})
feat = pd.read_csv(W30 / "features.tsv", sep="\t", dtype={"lig_id": str})
d0 = prim.merge(feat.drop(columns=["family", "ligand_charge"]), on="lig_id", validate="one_to_one").reset_index(drop=True)
p0, _ = nested(d0, M1)
r0 = float(np.sqrt(np.mean((p0 - d0.y_kcal.values) ** 2)))
print(f"STEP0 M1 La/Dy family-held-out RMSE {r0:.4f} (target 0.854, |diff| {abs(r0 - .854):.4f})")
if abs(r0 - .854) > .005:
    sys.exit("POSITIVE CONTROL FAILED")

base = prim[["lig_id", "name", "family", "normalized_identity", "ligand_charge"]].merge(
    feat.drop(columns=["family", "ligand_charge"]), on="lig_id")
base["D_anO"] = base[AN_O].sum(1); base["D_tot"] = base[TOT].sum(1)
P = pd.read_csv(OUT / "la_pairs_all_metals.tsv", sep="\t", dtype={"lig_id": str})
MODELS = {"M0": None, "M1": M1, "D_anO": ["D_anO"], "D_tot": ["D_tot"]}
TARGETS = {"y1": ("Ca/2+", "diff"), "y2": ("Ca/2+", "La"), "y3": ("Ca/2+", "X"), "y4": ("Sr/2+", "diff"),
           "y5": ("Al/3+", "diff"), "y6": ("Fe/3+", "diff"), "y7": ("Mg/2+", "diff"), "y8": ("Ba/2+", "diff"),
           "y9": ("Lu/3+", "diff"), "y10": ("Cd/2+", "diff"), "y11": ("Y/3+", "diff")}


def dataset(metal, kind, which):
    s = P[(P.metal == metal) & P.has_features].copy()
    if which == "primary":
        s = s[(s.match == "same_I") | (s.dI.abs() <= 0.1 + 1e-9)]
    elif which == "same_I":
        s = s[s.match == "same_I"]
    s["y_kcal"] = {"diff": s.logK_La - s.logK_X, "La": s.logK_La, "X": s.logK_X}[kind]  # log units; nested() reads y_kcal
    dd = s[["lig_id", "y_kcal", "match", "dI", "I_La", "I_X"]].merge(base, on="lig_id", validate="one_to_one")
    return dd.sort_values("lig_id").reset_index(drop=True)


def boot_idx(dd):
    n = len(dd); fam = dd.family.values; F = sorted(set(fam))
    BL = np.random.default_rng(SEED).integers(0, n, (NB, n))
    rf = np.random.default_rng(SEED); fi = {f: np.flatnonzero(fam == f) for f in F}
    BF = [np.concatenate([fi[g] for g in rf.choice(F, len(F), replace=True)]) for _ in range(NB)]
    return BL, BF


def stats(y, P_, idx):
    yy = y[idx]; sst = ((yy - yy.mean()) ** 2).sum()
    out = {}
    for k, p in P_.items():
        sse = ((p[idx] - yy) ** 2).sum()
        out[k] = (np.sqrt(sse / len(idx)), 1 - sse / sst if sst > 0 else np.nan)
    return out


def ci(a):
    a = np.asarray(a); a = a[np.isfinite(a)]
    return [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]


def score(dd, alphas=None, models=MODELS):
    if alphas is not None:
        evaluate.ALPHAS = alphas
    y = dd.y_kcal.values; preds = {}; sel = {}
    for k, cols in models.items():
        preds[k], meta = nested(dd, cols)
        sel[k] = [m["alpha"] for m in meta]
    evaluate.ALPHAS = FROZEN
    pt = stats(y, preds, np.arange(len(y)))
    BL, BF = boot_idx(dd)
    res = dict(n=len(dd), n_fam=int(dd.family.nunique()), sd=float(y.std()), alphas=sel, models={}, contrasts={})
    bs = {"lig": [stats(y, preds, b) for b in BL], "fam": [stats(y, preds, b) for b in BF]}
    for k in models:
        res["models"][k] = dict(rmse=float(pt[k][0]), r2=float(pt[k][1]),
                                **{f"rmse_ci_{u}": ci([s[k][0] for s in bs[u]]) for u in bs},
                                **{f"r2_ci_{u}": ci([s[k][1] for s in bs[u]]) for u in bs})
    for a, bb in [("M0", "M1"), ("M0", "D_anO"), ("M0", "D_tot"), ("D_anO", "M1")]:
        if a in models and bb in models:
            res["contrasts"][f"{a}-{bb}"] = dict(point=float(pt[a][0] - pt[bb][0]),
                                                 **{f"ci_{u}": ci([s[a][0] - s[bb][0] for s in bs[u]]) for u in bs})
    return res, preds


def show(tag, res):
    print(f"  [{tag}] n={res['n']} fam={res['n_fam']} sd={res['sd']:.2f}")
    for k, m in res["models"].items():
        print(f"    {k:6s} RMSE {m['rmse']:.3f} lig[{m['rmse_ci_lig'][0]:.2f},{m['rmse_ci_lig'][1]:.2f}] "
              f"fam[{m['rmse_ci_fam'][0]:.2f},{m['rmse_ci_fam'][1]:.2f}]  R2 {m['r2']:+.3f} "
              f"lig[{m['r2_ci_lig'][0]:+.2f},{m['r2_ci_lig'][1]:+.2f}] fam[{m['r2_ci_fam'][0]:+.2f},{m['r2_ci_fam'][1]:+.2f}]"
              f"  alphas {res['alphas'][k] if k != 'M0' else '-'}")
    for k, c in res["contrasts"].items():
        print(f"    RMSE({k}) {c['point']:+.3f} lig[{c['ci_lig'][0]:+.3f},{c['ci_lig'][1]:+.3f}] "
              f"fam[{c['ci_fam'][0]:+.3f},{c['ci_fam'][1]:+.3f}]")


ALL = {"step0_M1_LaDy_rmse_kcal": r0, "targets": {}}
PRED = {}
for t, (metal, kind) in TARGETS.items():
    dd = dataset(metal, kind, "primary")
    ok = len(dd) >= 25 and dd.family.nunique() >= 4
    print(f"\n== {t} ({kind} {metal}) primary n={len(dd)} families={dd.family.nunique()} -> "
          f"{'SCORED' if ok else 'below minimum, counts only'}")
    ALL["targets"][t] = {"metal": metal, "kind": kind, "primary_n": len(dd), "primary_fam": int(dd.family.nunique()),
                         "scored": ok}
    if not ok:
        for w in ("all",):
            d2 = dataset(metal, kind, w)
            ALL["targets"][t][f"{w}_n"] = len(d2); ALL["targets"][t][f"{w}_fam"] = int(d2.family.nunique())
            print(f"  ({w}: n={len(d2)}, families={d2.family.nunique()})")
        continue
    res, pr = score(dd); show("primary", res); ALL["targets"][t]["primary"] = res; PRED[t] = (dd, pr)
    for w in ("same_I", "all"):
        d2 = dataset(metal, kind, w)
        if len(d2) >= 10 and d2.family.nunique() >= 3:
            r2_, _ = score(d2, models={"M0": None, "M1": M1}); show(w, r2_); ALL["targets"][t][w] = r2_
    edge = any(a in (FROZEN[0], FROZEN[-1]) for a in res["alphas"]["M1"])
    ALL["targets"][t]["M1_grid_edge"] = edge
    if edge:
        rw, _ = score(dd, alphas=[.001] + FROZEN + [1000.], models={"M0": None, "M1": M1})
        show("OUTSIDE FROZEN DESIGN: grid 0.001-1000", rw); ALL["targets"][t]["wide_grid"] = rw

# ---------------- reading ----------------
m = ALL["targets"]["y1"]["primary"]
r2 = m["models"]["M1"]["r2"]; c = m["contrasts"]["M0-M1"]
reading = ("strong scale" if r2 >= .6 else
           "usable" if (r2 >= .3 and c["ci_lig"][0] > 0 and c["ci_fam"][0] > 0) else "not supported")
ALL["reading_y1"] = reading
print(f"\nFROZEN READING y1 (M1, primary): {reading}  (R2 {r2:+.3f}; RMSE(M0)-RMSE(M1) {c['point']:+.3f} "
      f"lig {c['ci_lig']} fam {c['ci_fam']})")

# ---------------- noise floor ----------------
nz = pd.read_csv(OUT / "la_ca_all_conditions.tsv", sep="\t", dtype={"lig_id": str})
for lab, s in (("all La/Ca ligands", nz), ("featured", nz[nz.has_features])):
    g = s.groupby("lig_id").filter(lambda x: x.I.nunique() >= 2).groupby("lig_id")
    ss = sum(((x.y1 - x.y1.mean()) ** 2).sum() for _, x in g); dof = sum(len(x) - 1 for _, x in g)
    rng_ = [x.y1.max() - x.y1.min() for _, x in g]
    print(f"NOISE FLOOR ({lab}): {g.ngroups} ligands, {dof} dof, pooled SD {np.sqrt(ss / dof) if dof else np.nan:.3f} log; "
          f"median range {np.median(rng_) if rng_ else np.nan:.3f}")
    ALL[f"noise_{lab}"] = dict(n_lig=g.ngroups, dof=dof, sd=float(np.sqrt(ss / dof)) if dof else None,
                              median_range=float(np.median(rng_)) if rng_ else None)

# ---------------- descriptive (y1 primary) ----------------
dd, pr = PRED["y1"]
y = dd.y_kcal.values; X = dd[M1].values.astype(float); fam = dd.family.values


def lofo_alpha(Xa, ya, fama):
    sc = []
    for a in FROZEN:
        e = np.full(len(ya), np.nan)
        for f in sorted(set(fama)):
            te = fama == f
            mdl = Ridge(alpha=a).fit(StandardScaler().fit(Xa[~te]).transform(Xa[~te]), ya[~te])
            e[te] = mdl.predict(StandardScaler().fit(Xa[~te]).transform(Xa[te])) - ya[te]
        sc.append(np.mean(e ** 2))
    return FROZEN[int(np.argmin(sc))]


def raw_coef(Xa, ya, a):
    s = StandardScaler().fit(Xa); mdl = Ridge(alpha=a).fit(s.transform(Xa), ya)
    sc_ = np.where(s.scale_ > 0, s.scale_, np.inf)
    return mdl.coef_ / sc_


a_full = lofo_alpha(X, y, fam)
cf = raw_coef(X, y, a_full)
BL, BF = boot_idx(dd)
bl = np.array([raw_coef(X[b], y[b], a_full) for b in BL]); bf = np.array([raw_coef(X[b], y[b], a_full) for b in BF])
print(f"\n== y1 M1 coefficients (log units per group), alpha={a_full}, n={len(y)} ==")
coef_rows = []
for j, k in enumerate(M1):
    nz_ = int((X[:, j] != 0).sum())
    coef_rows.append(dict(feature=k, n_nonzero=nz_, coef=cf[j], lig_lo=np.percentile(bl[:, j], 2.5),
                          lig_hi=np.percentile(bl[:, j], 97.5), fam_lo=np.percentile(bf[:, j], 2.5),
                          fam_hi=np.percentile(bf[:, j], 97.5)))
CR = pd.DataFrame(coef_rows); print(CR.to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
CR.to_csv(OUT / "y1_coefficients.tsv", sep="\t", index=False)

print("\n== y1 by carboxylate count ==")
cc = dd.groupby("carboxyl").y_kcal.agg(["size", "median", "min", "max"]); print(cc.to_string(float_format=lambda v: f"{v:+.2f}"))


def donors(r):
    return ",".join(f"{k}{int(r[k])}" for k in GROUPS if r[k]) + (",macro" if r.macrocycle else "") + f",z{int(r.ligand_charge)}"


dd["donors"] = dd.apply(donors, axis=1); dd["pred_M1"] = pr["M1"]
srt = dd.sort_values("y_kcal")
print("\n== y1 ten largest ==")
print(srt.tail(10)[::-1][["lig_id", "name", "family", "y_kcal", "pred_M1", "donors", "I_La", "dI"]].to_string(index=False))
print("\n== y1 ten smallest ==")
print(srt.head(10)[["lig_id", "name", "family", "y_kcal", "pred_M1", "donors", "I_La", "dI"]].to_string(index=False))
dd.to_csv(OUT / "y1_primary_predictions.tsv", sep="\t", index=False)

ls = pd.read_csv(REPO / "workspaces/ln_lattice_strain_20261109/nist_per_ligand_fits.tsv", sep="\t", dtype={"lig_id": str})
ab = P[(P.metal == "Ca/2+") & ((P.match == "same_I") | (P.dI.abs() <= .1 + 1e-9))][["lig_id", "y_log"]].merge(
    ls[["lig_id", "a2", "b2"]], on="lig_id")
rho = spearmanr(ab.y_log, ab.a2)[0]; rhob = spearmanr(ab.y_log, ab.b2)[0]
bb = np.random.default_rng(SEED).integers(0, len(ab), (NB, len(ab)))
rb = [spearmanr(ab.y_log.values[i], ab.a2.values[i])[0] for i in bb]
print(f"\nSpearman y1 vs lattice-strain a2: rho {rho:+.3f} CI {ci(rb)} n={len(ab)} (all La/Ca primary-rule pairs, "
      f"featured or not); vs b2 rho {rhob:+.3f}")
ALL["spearman_a2"] = dict(rho=rho, ci=ci(rb), n=len(ab), rho_b2=rhob)
ALL["y1_coef_alpha"] = a_full
json.dump(ALL, open(OUT / "summary.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("wrote", OUT)
