"""Step 1, readings A4 and A5 of PLAN.md: does sequence predict profiles in unseen families, and how close must a
neighbour be? Genus-held-out and component-wise scores are added and labelled post hoc."""
import json, sys
import numpy as np, pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, KFold
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
LOG = open(ROOT / "results/first_look_b.log", "w"); SUM = {}
def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
d, raw, C, mean, valid, qc = load()
o = pd.read_csv(ROOT / "data/ORTHOLOGS.tsv", sep="\t"); pos = {a: i for i, a in enumerate(o.accession)}
I616 = np.load(ROOT / "data/identity_616.npy"); E616 = np.load(ROOT / "data/esm2_650M_mean_616.npy")
idx = np.where(qc)[0]; g = np.array([pos[a] for a in d.accession.values[idx]])
Y = mean[idx]; I = I616[np.ix_(g, g)]; E = E616[g]; n = len(idx)
fam = d.family_original.fillna("unspecified").astype(str).values[idx]; gen = d.genus_original.fillna("unspecified").astype(str).values[idx]
F = pd.DataFrame([loop_features(s) for s in d[SEQ].values[idx]])[["n_loops", "acid_per_loop", "n_D9", "n_P2", "net_charge_per_res", "length"]].values.astype(float)
p(f"QC orthologs {n}; median identity to nearest other QC ortholog {np.median(np.sort(I, 1)[:, -2]):.3f}")

def ridge_cv(Xtr, Ytr, Xte, groups, alphas):
    mu, sd = np.nanmean(Xtr, 0), np.nanstd(Xtr, 0) + 1e-9
    Xtr = np.where(np.isnan(Xtr), mu, Xtr); Xte = np.where(np.isnan(Xte), mu, Xte)
    A, B = (Xtr - mu) / sd, (Xte - mu) / sd; ym = Ytr.mean(0)
    ng = len(set(groups)); best = (1e99, alphas[0])
    if ng >= 3:
        for a in alphas:
            e = 0.0
            for tr, va in GroupKFold(n_splits=min(5, ng)).split(A, groups=groups):
                m = Ridge(alpha=a).fit(A[tr], Ytr[tr] - Ytr[tr].mean(0)); e += (((m.predict(A[va]) + Ytr[tr].mean(0)) - Ytr[va]) ** 2).sum()
            if e < best[0]: best = (e, a)
    m = Ridge(alpha=best[1]).fit(A, Ytr - ym)
    return m.predict(B) + ym, best[1]

def predict_all(tr, te):
    out = {}; mu = Y[tr].mean(0)
    out["B0"] = np.tile(mu, (len(te), 1))
    S = I[np.ix_(te, tr)]
    out["B1"] = Y[tr][S.argmax(1)]
    top = np.argsort(-S, 1)[:, :5]; w = np.take_along_axis(S, top, 1); w = w / w.sum(1, keepdims=True)
    out["B2"] = (Y[tr][top] * w[:, :, None]).sum(1)
    out["M1"], a1 = ridge_cv(F[tr], Y[tr], F[te], gen[tr], [0.1, 1, 10, 100])
    out["M2"], a2 = ridge_cv(E[tr], Y[tr], E[te], gen[tr], [1, 10, 100, 1000, 10000])
    return out, mu, S.max(1)

def score(folds, label):
    P = {k: np.zeros_like(Y) for k in ("B0", "B1", "B2", "M1", "M2")}; MU = np.zeros_like(Y); NN = np.zeros(n); tested = np.zeros(n, bool); grp = np.empty(n, object)
    for name, te in folds:
        tr = np.setdiff1d(np.arange(n), te); out, mu, nn = predict_all(tr, te)
        for k in P: P[k][te] = out[k]
        MU[te] = mu; NN[te] = nn; tested[te] = True; grp[te] = name
    names = sorted(set(grp[tested])); res = {}
    se = {k: ((P[k] - Y) ** 2).sum(1) for k in P}
    cos = {k: np.array([np.dot(P[k][i] - MU[i], Y[i] - MU[i]) / (np.linalg.norm(P[k][i] - MU[i]) * np.linalg.norm(Y[i] - MU[i]) + 1e-12) for i in range(n)]) for k in P if k != "B0"}
    rng = np.random.default_rng(42)
    p(f"\n{label}: {len(names)} held-out groups, {int(tested.sum())} orthologs tested")
    for k in ("B1", "B2", "M1", "M2"):
        fm = np.array([se[k][grp == f].mean() for f in names]); f0 = np.array([se["B0"][grp == f].mean() for f in names]); fc = np.array([cos[k][grp == f].mean() for f in names])
        skill = 1 - fm.mean() / f0.mean(); bs = []; bc = []
        for _ in range(2000):
            j = rng.integers(0, len(names), len(names)); bs.append(1 - fm[j].mean() / f0[j].mean()); bc.append(fc[j].mean())
        pooled = 1 - se[k][tested].mean() / se["B0"][tested].mean()
        p(f"  {k}: skill {skill:+.3f} [{np.percentile(bs, 2.5):+.3f}, {np.percentile(bs, 97.5):+.3f}] (groups weighted equally); pooled over orthologs {pooled:+.3f}; "
          f"mean cosine {fc.mean():+.3f} [{np.percentile(bc, 2.5):+.3f}, {np.percentile(bc, 97.5):+.3f}]; RMSE per element {np.sqrt(se[k][tested].mean() / 15):.3f} vs B0 {np.sqrt(se['B0'][tested].mean() / 15):.3f} log10")
        res[k] = dict(skill=float(skill), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)), pooled=float(pooled), cosine=float(fc.mean()),
                      cos_lo=float(np.percentile(bc, 2.5)), cos_hi=float(np.percentile(bc, 97.5)))
    if len(names) <= 12:
        for f in names:
            m = grp == f
            p(f"     {f:22s} n={int(m.sum()):3d} nearest-training identity median {np.median(NN[m]):.2f} | skill " + " ".join(f"{k} {1 - se[k][m].mean() / se['B0'][m].mean():+.2f}" for k in ("B1", "B2", "M1", "M2")))
    return res, P, MU, NN, tested

# A4 primary: leave one family out (NCBI families with >= 3 QC orthologs; 'unspecified' is never a fold)
vc = pd.Series(fam).value_counts(); fams = [f for f in vc.index if vc[f] >= 3 and f != "unspecified"]
SUM["A4_family"], *_ = score([(f, np.where(fam == f)[0]) for f in fams], "A4. Leave one family out")
b1 = SUM["A4_family"]["B1"]
ok = b1["lo"] > 0 and b1["cosine"] >= 0.3
p(f"  reading A4 (B1 skill interval above 0 and mean cosine >= 0.3): {'PASS' if ok else 'FAIL'}"); SUM["A4_pass"] = bool(ok)
# post hoc: leave one genus out (genera with >= 3 QC orthologs)
vg = pd.Series(gen).value_counts(); gens = [x for x in vg.index if vg[x] >= 3 and x != "unspecified"]
SUM["genus"], *_ = score([(x, np.where(gen == x)[0]) for x in gens], "Post hoc. Leave one genus out")
# A5: ten-fold random CV, skill by identity to the nearest training ortholog
kf = KFold(n_splits=10, shuffle=True, random_state=42)
res, P, MU, NN, tested = score([(f"fold{j}", te) for j, (tr, te) in enumerate(kf.split(np.arange(n)))], "A5. Ten random folds")
SUM["A5_overall"] = res
se = {k: ((P[k] - Y) ** 2).sum(1) for k in P}
p("  by identity to the nearest training ortholog:")
SUM["A5_bins"] = {}
for lo, hi in ((0, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 0.95), (0.95, 1.01)):
    m = (NN >= lo) & (NN < hi)
    if m.sum() >= 3:
        line = f"    {lo:.2f}-{hi:.2f}: n={int(m.sum()):3d} | skill " + " ".join(f"{k} {1 - se[k][m].mean() / se['B0'][m].mean():+.2f}" for k in ("B1", "B2", "M1", "M2")) + f" | RMSE per element B0 {np.sqrt(se['B0'][m].mean() / 15):.3f}, B1 {np.sqrt(se['B1'][m].mean() / 15):.3f}, M2 {np.sqrt(se['M2'][m].mean() / 15):.3f}"
        p(line); SUM["A5_bins"][f"{lo}-{hi}"] = dict(n=int(m.sum()), **{k: float(1 - se[k][m].mean() / se["B0"][m].mean()) for k in ("B1", "B2", "M1", "M2")})
# how much of the between-ortholog variance is the tilt, and how well is the tilt alone predicted (post hoc)
co, yoff, _ = components(Y, 3)
for k in ("B1", "M2"):
    cp, yp, _ = components(P[k], 3)
    r = [np.corrcoef(co[:, j], cp[:, j])[0, 1] for j in (1, 2, 3)] + [np.corrcoef(yoff, yp)[0, 1]]
    p(f"  post hoc, random folds, {k}: correlation of predicted with observed linear {r[0]:+.2f}, quadratic {r[1]:+.2f}, cubic {r[2]:+.2f}, yttrium offset {r[3]:+.2f}")
np.savez(ROOT / "results/cv_predictions_random10.npz", **{f"P_{k}": v for k, v in P.items()}, Y=Y, NN=NN)
json.dump(SUM, open(ROOT / "results/first_look_b.json", "w"), indent=1, default=float)
