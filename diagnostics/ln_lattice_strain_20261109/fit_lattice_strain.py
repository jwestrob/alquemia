"""Lattice-strain (Brice / Blundy-Wood) test of the NIST La..Lu two-shape basis.  Frozen design: PLAN.md (same dir).

Login node, single thread, seconds.  Reads only the inputs named in PLAN.md; writes NIST-derived tables to
workspaces/ln_lattice_strain_20261109/ (untracked) and prints the summary used in REPORT.md.
"""
import json, re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "workspaces/ln_lattice_strain_20261109"
OUT.mkdir(exist_ok=True)
KC = 1.3643
LN = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
R9 = np.array([1.216, 1.196, 1.179, 1.163, 1.132, 1.120, 1.107, 1.095, 1.083, 1.072, 1.062, 1.052, 1.042, 1.032])
R8 = np.array([1.160, 1.143, 1.126, 1.109, 1.079, 1.066, 1.053, 1.040, 1.027, 1.015, 1.004, 0.994, 0.985, 0.977])
ZZ = np.array([57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71], float)
BRICE = 4 * np.pi * 6.02214076e23 * 1e9 * 1e-30 / 4184.0  # kcal/mol per GPa*A^3
assert abs(BRICE - 1.809) < 5e-4
NB, SEED = 2000, 42
G_LO, G_HI = 0.300, 3.000
GRID = np.round(np.arange(G_LO, G_HI + 1e-9, 0.001), 3)

# ---- radius check against the mendeleev Shannon table ----------------------------------------------------------
try:
    from mendeleev import element
    bad = []
    for m, a9, a8 in zip(LN, R9, R8):
        rr = {r.coordination: r.ionic_radius / 100 for r in element(m).ionic_radii if r.charge == 3}
        if abs(rr.get("IX", -1) - a9) > 5e-4 or abs(rr.get("VIII", -1) - a8) > 5e-4:
            bad.append(m)
    print(f"radius check vs mendeleev Shannon table: {'all 28 match' if not bad else 'MISMATCH ' + str(bad)}")
except Exception as e:  # pragma: no cover
    print("radius check unavailable:", e)

# ---- NIST series (parsed exactly as series_rank.py) ----------------------------------------------------------------
d = pd.read_csv(REPO / "workspaces/ln_surrogate_20261027/surrogate_nist.tsv", sep="\t", dtype={"lig_id": str})
rows, ids, names = [], [], []
for _, r in d.iterrows():
    k = dict(x.split(":") for x in r.logK_series.split(";"))
    if all(m in k for m in LN):
        rows.append([-KC * (float(k[m]) - float(k["La"])) for m in LN[1:]]); ids.append(r.lig_id); names.append(r["name"])
X = np.array(rows); n, p = X.shape
assert (n, p) == (47, 13), (n, p)
ev = np.linalg.svd(X, compute_uv=False) ** 2; ev /= ev.sum()
print(f"PC1 inline: n={n}, SVD variance {ev[0]:.3f} {ev[1]:.3f} {ev[2]:.3f}")

dr9, dr8, z = R9[1:] - R9[0], R8[1:] - R8[0], ZZ[1:] - 57
iGd = LN.index("Gd") - 1
LIN = {
    "F2": np.c_[dr9, dr9 ** 2],
    "F3": np.c_[dr8, dr8 ** 2],
    "F5": np.c_[z, z ** 2],
    "F6": np.c_[dr9, np.minimum(0.0, dr9 - dr9[iGd])],
    "F7": dr9[:, None],
}
ALL = np.arange(13)
FIT2 = np.array([LN.index(m) - 1 for m in ("Nd", "Dy", "Lu")])
FIT1 = np.array([LN.index("Dy") - 1])
HO2 = np.setdiff1d(ALL, FIT2); HO1 = np.setdiff1d(ALL, FIT1)

LOOV = []
for i in range(n):
    LOOV.append(np.linalg.svd(np.delete(X, i, 0), full_matrices=False)[2])


def linfit(A, y, f):
    c = np.linalg.lstsq(A[f], y[f], rcond=None)[0]
    return A @ c, c


def brice_S(r, r0):
    r = np.asarray(r)[None, :]; r0 = np.asarray(r0)[:, None]
    return BRICE * ((r0 / 2) * (r - r0) ** 2 + (r - r0) ** 3 / 3)


def f4_fit(r, y, f, intercept=False, ref=None, lo=G_LO, hi=G_HI, epos=False):
    """y[f] = [c] + E*(S(r[f]) - S(ref)).  Profile r0 on the frozen grid, then bounded refine +-0.002 A.
    Returns r0, E, c, prediction at all r, profiled-RMS curve over the grid."""
    yf = y[f]

    def prof(r0s):
        g = brice_S(r[f], r0s) - (brice_S([ref], r0s) if ref is not None else 0.0)
        if intercept:
            gc = g - g.mean(1, keepdims=True); yc = yf - yf.mean()
        else:
            gc, yc = g, yf
        gg = (gc * gc).sum(1); gy = gc @ yc
        E = np.where(gg > 0, gy / np.where(gg > 0, gg, 1), 0.0)
        if epos:
            E = np.maximum(E, 0.0)
        sse = (yc @ yc) - 2 * E * gy + E * E * gg
        c = (yf.mean() - E * g.mean(1)) if intercept else np.zeros_like(E)
        return sse, E, c

    grid = GRID[(GRID >= lo - 1e-9) & (GRID <= hi + 1e-9)]
    sse, _, _ = prof(grid)
    j = int(np.argmin(sse)); a, b = max(lo, grid[j] - 0.002), min(hi, grid[j] + 0.002)
    res = minimize_scalar(lambda t: prof(np.array([t]))[0][0], bounds=(a, b), method="bounded",
                          options={"xatol": 1e-7})
    r0 = res.x if res.fun <= sse[j] else grid[j]
    s, E, c = prof(np.array([r0]))
    g_all = brice_S(r, [r0])[0] - (brice_S([ref], [r0])[0, 0] if ref is not None else 0.0)
    return r0, E[0], c[0], c[0] + E[0] * g_all, (grid, np.sqrt(np.maximum(sse, 0) / len(f)))


# ---- per-ligand error matrices for every form x score ---------------------------------------------------------------
forms2 = ["F1", "F2", "F3", "F4", "F4c", "F5", "F6"]
forms1 = ["F7", "F8"]
ERR = {}       # (form, score) -> 47 x 13 error matrix (pred - obs)
PAR = {k: [] for k in ("a2", "b2", "r0", "E", "r0_c", "E_c", "r0_flat_lo", "r0_flat_hi")}
for i in range(n):
    y = X[i]; V = LOOV[i]
    for sc, f in (("S1", ALL), ("S2", FIT2)):
        for fm in forms2:
            if fm == "F1":
                pred, _ = linfit(V[:2].T, y, f)
            elif fm == "F4":
                r0, E, _, pred, (gr, prms) = f4_fit(R9[1:], y, f, ref=R9[0])
                if sc == "S1":
                    PAR["r0"].append(r0); PAR["E"].append(E)
                    ok = gr[prms <= prms.min() + 0.01]
                    PAR["r0_flat_lo"].append(ok.min()); PAR["r0_flat_hi"].append(ok.max())
            elif fm == "F4c":
                r0, E, _, pred, _ = f4_fit(R9[1:], y, f, ref=R9[0], lo=0.8, hi=1.4, epos=True)
                if sc == "S1":
                    PAR["r0_c"].append(r0); PAR["E_c"].append(E)
            else:
                pred, cc = linfit(LIN[fm], y, f)
                if fm == "F2" and sc == "S1":
                    PAR["a2"].append(cc[0]); PAR["b2"].append(cc[1])
            ERR.setdefault((fm, sc), np.zeros((n, p)))[i] = pred - y
    for fm in forms1:
        for sc, f in (("S1", ALL), ("S3", FIT1)):
            A = V[:1].T if fm == "F8" else LIN["F7"]
            pred, _ = linfit(A, y, f)
            ERR.setdefault((fm, sc), np.zeros((n, p)))[i] = pred - y


def regmap(cols):  # series_rank.py loo(), verbatim logic
    idx = [LN[1:].index(c) for c in cols]; err = []
    for i in range(n):
        tr = np.arange(n) != i
        B = np.linalg.lstsq(X[tr][:, idx], X[tr], rcond=None)[0]
        err.append(X[i, idx] @ B - X[i])
    return np.array(err)


ERR[("F1r", "S2")] = regmap(["Nd", "Dy", "Lu"])
ERR[("F8r", "S3")] = regmap(["Dy"])
print(f"PC1 inline: regression map RMS all-13  Dy {np.sqrt((ERR[('F8r','S3')]**2).mean()):.3f}"
      f"  Nd+Dy+Lu {np.sqrt((ERR[('F1r','S2')]**2).mean()):.3f}")

rng = np.random.default_rng(SEED)
BI = rng.integers(0, n, (NB, n))
MASK = {"S1": ALL, "S2": HO2, "S3": HO1}


def ssv(E, cols):
    return (E[:, cols] ** 2).sum(1), len(cols)


def rms_ci(E, cols):
    ss, m = ssv(E, cols)
    pt = np.sqrt(ss.sum() / (n * m)); bs = np.sqrt(ss[BI].sum(1) / (n * m))
    return pt, np.percentile(bs, 2.5), np.percentile(bs, 97.5), bs


summary = {"scores": {}}
print("\n== Scores (kcal/mol): point [95% CI]; S2/S3 primary = held-out metals; 'all13' = series_rank convention ==")
BS = {}
for (fm, sc), E in sorted(ERR.items(), key=lambda t: (t[0][1], t[0][0])):
    pt, lo, hi, bs = rms_ci(E, MASK[sc]); BS[(fm, sc)] = bs
    a13 = np.sqrt((E ** 2).mean())
    summary["scores"][f"{fm}_{sc}"] = {"rms": pt, "lo": lo, "hi": hi, "rms_all13": a13}
    print(f"  {sc} {fm:4s} {pt:.3f} [{lo:.3f}, {hi:.3f}]" + (f"   all13 {a13:.3f}" if sc != "S1" else ""))


def diff(f_a, f_b, sc):
    pa = summary["scores"][f"{f_a}_{sc}"]["rms"]; pb = summary["scores"][f"{f_b}_{sc}"]["rms"]
    dd = BS[(f_a, sc)] - BS[(f_b, sc)]
    return pa - pb, np.percentile(dd, 2.5), np.percentile(dd, 97.5)


# ---- frozen readings ---------------------------------------------------------------------------------------------
s2 = {fm: summary["scores"][f"{fm}_S2"]["rms"] for fm in ("F2", "F3", "F4")}
best = min(s2, key=s2.get)
f1 = summary["scores"]["F1_S2"]["rms"]; f1r = summary["scores"]["F1r_S2"]["rms"]
f1r_ho = np.sqrt((ERR[("F1r", "S2")][:, HO2] ** 2).mean())
dk = diff("F6", "F2", "S1")
read = {
    "best_radius_form": best, "best_S2": s2[best],
    "adequate": bool(s2[best] <= 0.20),
    "equivalent_to_F1": bool(s2[best] <= f1 + 0.03), "F1_S2": f1,
    "equivalent_to_F1r_heldout": bool(s2[best] <= f1r_ho + 0.03), "F1r_S2_heldout": f1r_ho,
    "S1_F6_minus_F2": dk,
    "bend": "smooth" if dk[1] > 0 else ("kink" if dk[2] < 0 else "undetermined"),
}
for a_, b_ in (("F2", "F1"), ("F3", "F1"), ("F4", "F1"), ("F2", "F3"), ("F2", "F5")):
    read[f"S2_{a_}_minus_{b_}"] = diff(a_, b_, "S2")
read["S3_F7_minus_F8"] = diff("F7", "F8", "S3")
summary["readings"] = read
print("\n== Frozen readings ==")
for k_, v in read.items():
    print(f"  {k_}: {v}")

# ---- residual-by-metal for F2 + universal per-metal offset ---------------------------------------------------------
Res = -ERR[("F2", "S1")]  # obs - pred
mu = Res.mean(0); se = Res.std(0, ddof=1) / np.sqrt(n); t = mu / se
share = n * (mu ** 2).sum() / (Res ** 2).sum()
r2 = np.random.default_rng(SEED); stat0 = (mu ** 2).sum(); cnt = 0
for _ in range(NB):
    s = r2.choice([-1.0, 1.0], n)
    cnt += ((s[:, None] * Res).mean(0) ** 2).sum() >= stat0
p_flip = (1 + cnt) / (NB + 1)
print("\n== F2 (CN9) in-sample residual (obs - fit) by metal, kcal/mol, mean +- SE (t) ==")
for m, a_, b_, c_ in zip(LN[1:], mu, se, t):
    print(f"  {m}: {a_:+.3f} +- {b_:.3f} ({c_:+.1f})")
print(f"  share of residual SS in per-metal means {share:.3f}; sign-flip p {p_flip:.4f}; |t|>2 in {int((abs(t)>2).sum())}/13")
off = {}
for sc, f in (("S1", ALL), ("S2", FIT2)):
    E = np.zeros((n, p))
    for i in range(n):
        c = (Res.sum(0) - Res[i]) / (n - 1)
        pred, _ = linfit(LIN["F2"], X[i] - c, f)
        E[i] = pred + c - X[i]
    pt, lo, hi, bs = rms_ci(E, MASK[sc])
    dd = bs - BS[("F2", sc)]
    off[sc] = {"rms": pt, "lo": lo, "hi": hi, "delta_vs_F2": pt - summary["scores"][f"F2_{sc}"]["rms"],
               "delta_lo": np.percentile(dd, 2.5), "delta_hi": np.percentile(dd, 97.5)}
    print(f"  F2 + LOO per-metal offset {sc}: {pt:.3f} [{lo:.3f}, {hi:.3f}]; change vs F2 {off[sc]['delta_vs_F2']:+.3f} "
          f"[{off[sc]['delta_lo']:+.3f}, {off[sc]['delta_hi']:+.3f}]")
summary["per_metal"] = {"mean": dict(zip(LN[1:], mu)), "se": dict(zip(LN[1:], se)), "share": share, "p_flip": p_flip,
                        "offset": off}

# ---- added descriptive checks (DEVIATIONS.md item 3): basis geometry and the shape of the per-metal residual ------
print("\n== Added descriptive: SVD shapes vs span(dr, dr^2) (CN9), all 47 ligands ==")
Vfull = np.linalg.svd(X, full_matrices=False)[2]
Q = np.linalg.qr(LIN["F2"])[0]
for k in range(3):
    v = Vfull[k]; r2_ = 1 - ((v - Q @ (Q.T @ v)) ** 2).sum() / (v ** 2).sum()
    print(f"  shape {k+1}: fraction of its norm^2 inside span(dr, dr^2) = {r2_:.3f}")
cosang = np.linalg.svd(Q.T @ np.linalg.qr(Vfull[:2].T)[0], compute_uv=False)
print(f"  principal angles between span(V1,V2) and span(dr,dr^2): {np.degrees(np.arccos(np.clip(cosang,-1,1))).round(1)} deg")
Xhat = X @ Vfull[:2].T @ Vfull[:2]
fr = 1 - ((Xhat - (Xhat @ Q) @ Q.T) ** 2).sum() / (Xhat ** 2).sum()
print(f"  share of the rank-2 SVD reconstruction (norm^2) inside span(dr, dr^2): {fr:.4f}")
Pperp = np.eye(13) - Q @ Q.T
cand = {"Gd step (Gd..Lu = 1)": (np.arange(13) >= iGd).astype(float), "dr^3": dr9 ** 3,
        "Gd step + dr^3": None}
print("  per-metal mean residual of F2 explained (R^2, through origin) by quadratic-orthogonalised candidates:")
for nm, v in cand.items():
    B = np.c_[Pperp @ cand["Gd step (Gd..Lu = 1)"], Pperp @ (dr9 ** 3)] if v is None else (Pperp @ v)[:, None]
    fit_ = B @ np.linalg.lstsq(B, mu, rcond=None)[0]
    print(f"    {nm:22s} R^2 {1 - ((mu - fit_) ** 2).sum() / (mu ** 2).sum():.3f}")

# ---- F4 parameter distributions ------------------------------------------------------------------------------------
P = pd.DataFrame({"lig_id": ids, "name": names, **{k: np.array(v) for k, v in PAR.items()}})
lo_dr = dr9.min()
P["vertex_dr"] = -P.a2 / (2 * P.b2)
P["monotonic"] = ~((P.vertex_dr > lo_dr) & (P.vertex_dr < 0))
P["boundary"] = (P.r0 < G_LO + 0.005) | (P.r0 > G_HI - 0.005)
P["unphysical"] = (P.r0 < 0.8) | (P.r0 > 1.4) | (P.E <= 0)
P["flat_width"] = P.r0_flat_hi - P.r0_flat_lo
for k_ in ("F1", "F2", "F4", "F6", "F7"):
    P[f"rms_S1_{k_}"] = np.sqrt((ERR[(k_, "S1")] ** 2).mean(1))
print("\n== F4 (Brice, E sign free) over 47 ligands ==")
q = lambda s: f"median {s.median():.3f}, IQR [{s.quantile(.25):.3f}, {s.quantile(.75):.3f}], range [{s.min():.3f}, {s.max():.3f}]"
print(f"  r0 (A): {q(P.r0)}\n  E (GPa): {q(P.E)}")
print(f"  boundary {int(P.boundary.sum())}, unphysical {int(P.unphysical.sum())} "
      f"(r0<0.8 {int((P.r0<0.8).sum())}, r0>1.4 {int((P.r0>1.4).sum())}, E<=0 {int((P.E<=0).sum())})")
print(f"  monotonic ligands {int(P.monotonic.sum())}: unphysical {int((P.unphysical & P.monotonic).sum())}, "
      f"boundary {int((P.boundary & P.monotonic).sum())}; non-monotonic {int((~P.monotonic).sum())}: "
      f"unphysical {int((P.unphysical & ~P.monotonic).sum())}")
print(f"  r0 interval within +0.01 kcal/mol of the profiled minimum: median width {P.flat_width.median():.3f} A "
      f"(monotonic {P.flat_width[P.monotonic].median():.3f}, non-monotonic {P.flat_width[~P.monotonic].median():.3f})")
phys = P[~P.unphysical]
if len(phys):
    print(f"  physical fits n={len(phys)}: r0 {q(phys.r0)}; E {q(phys.E)}")
print(f"  F4c (E>=0, r0 in [0.8,1.4]): r0 at a bound in {int(((P.r0_c<0.8005)|(P.r0_c>1.3995)).sum())}, E=0 in {int((P.E_c<=0).sum())}")
summary["F4"] = {"n_boundary": int(P.boundary.sum()), "n_unphysical": int(P.unphysical.sum()),
                 "n_monotonic": int(P.monotonic.sum()), "n_unphys_monotonic": int((P.unphysical & P.monotonic).sum()),
                 "r0_median": P.r0.median(), "E_median": P.E.median()}
P.to_csv(OUT / "nist_per_ligand_fits.tsv", sep="\t", index=False, float_format="%.5g")

# ---- Spearman with ligand features ---------------------------------------------------------------------------------
F = pd.read_csv(REPO / "workspaces/ln_selectivity_calibration_20261027/features.tsv", sep="\t", dtype={"lig_id": str})
DON = ["carboxyl", "aminoacid_N", "amine_N", "ether_O", "alcohol", "phenol", "phosphonate", "pyridine_N", "amide_O",
       "carbonyl_O", "sulfonate", "hydroxamate"]
F["donor_total"] = F[DON].sum(axis=1)
# charge (q_lig) and denticity (nlig_La_cn8) live in surrogate_nist.tsv (MACE-built complexes), not features.tsv
M = P.merge(F[["lig_id", "donor_total", "carboxyl"]], on="lig_id", how="inner").merge(
    d[["lig_id", "q_lig", "nlig_La_cn8"]].drop_duplicates("lig_id"), on="lig_id", how="left")
print(f"  feature coverage: q_lig {int(M.q_lig.notna().sum())}, nlig_La_cn8 {int(M.nlig_La_cn8.notna().sum())}")
print(f"\n== Spearman rho (bootstrap 95% CI) over ligands in features.tsv: n={len(M)} of {n} ==")
r3 = np.random.default_rng(SEED)
spear = {}
for par, sub in (("a2", M), ("b2", M), ("r0", M[~M.unphysical]), ("E", M[~M.unphysical])):
    for feat in ("donor_total", "q_lig", "carboxyl", "nlig_La_cn8"):
        s = sub[[par, feat]].dropna()
        if len(s) < 8:
            print(f"  {par:3s} vs {feat:12s}: n={len(s)} too few"); continue
        rho = spearmanr(s[par], s[feat])[0]
        bi = r3.integers(0, len(s), (NB, len(s))); v = s.values
        bs = np.array([spearmanr(v[b, 0], v[b, 1])[0] for b in bi])
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        spear[f"{par}~{feat}"] = (rho, lo, hi, len(s))
        print(f"  {par:3s} vs {feat:12s}: rho {rho:+.2f} [{lo:+.2f}, {hi:+.2f}]  n={len(s)}")
summary["spearman"] = spear

# ---- protein / peptide datasets (builder copied from protein_series_shapes.py) ------------------------------------
NS = np.c_[np.zeros(n), -X / KC]  # log K(M) - log K(La), ligands x 14


def basis(Xs, k):
    _, _, Vt = np.linalg.svd(Xs[:, 1:], full_matrices=False)
    return np.vstack([np.zeros(k), Vt[:k].T])


dd_ = pd.read_csv(REPO / "workspaces/ln_selectivity_calibration_20261027/protein_ln_affinities.tsv", sep="\t")
dd_ = dd_[dd_.metal.isin(LN)].copy()
m_ = (dd_.system == "Mex-LanM (agarose-immobilized)") & dd_.log_K.isna()
dd_.loc[m_, "log_K"] = [np.log10(float(re.search(r"/Keq_rel\(La\) = ([0-9.eE+-]+)", t_).group(1))) / 2 if "/Keq_rel(La) =" in t_ else
                        (0.0 if mm == "La" else np.nan) for t_, mm in zip(dd_.loc[m_, "notes"], dd_.loc[m_, "metal"])]
dd_ = dd_[dd_.log_K.notna() & ~dd_.constant_type.isin(["Km", "EC50_or_K_half_activity"])]
key = ["system", "variant", "doi", "method", "pH", "T_C", "constant_type", "site", "competitor", "buffer_ionic_strength"]
SHORT = [LN.index(m) for m in ("La", "Nd", "Dy", "Lu")]
dR9, dR8, Zf = R9 - R9[0], R8 - R8[0], ZZ - 57
PLIN = {"F2": np.c_[dR9, dR9 ** 2], "F3": np.c_[dR8, dR8 ** 2], "F5": np.c_[Zf, Zf ** 2],
        "F6": np.c_[dR9, np.minimum(0.0, dR9 - dR9[LN.index("Gd")])], "F7": dR9[:, None]}
V2 = basis(NS, 2)


def pfit(form, yG, idx, fi):
    """yG = -KC*logK (kcal). Fit c + f on fi, return RMS over idx-minus-fi (or idx if fi is None)."""
    f = np.array(idx if fi is None else fi); ev_ = np.array(idx if fi is None else [j for j in idx if j not in fi])
    if form == "F4":
        _, _, _, pred, _ = f4_fit(R9, yG, f, intercept=True)
    else:
        A = np.c_[np.ones(14), V2 if form == "F1" else PLIN[form]]
        pred, _ = linfit(A, yG, f)
    return np.sqrt(np.mean((yG[ev_] - pred[ev_]) ** 2))


old = pd.read_csv(REPO / "workspaces/lady_series_20261028/protein_series_shapes.tsv", sep="\t")
prow = []
for k_, s in dd_.fillna({c: "-" for c in key}).groupby(key):
    s = s.groupby("metal").log_K.mean()
    if len(s) < 6:
        continue
    y = np.full(14, np.nan); y[[LN.index(x) for x in s.index]] = s.values
    idx = [i for i in range(14) if not np.isnan(y[i])]
    yG = -KC * y
    rec = {"dataset": f"{k_[0][:34]} | {k_[1][:22]}", "doi": k_[2], "n": len(idx), "res2_recomputed": pfit("F1", yG, idx, None)}
    for fm in ("F2", "F3", "F4", "F5", "F6", "F7"):
        rec[f"res_{fm}"] = pfit(fm, yG, idx, None)
    nref = [pfit("F2", -KC * NS[i], idx, None) for i in range(n)]
    rec["nist_F2_med"] = np.median(nref)
    if all(j in idx for j in SHORT) and len(idx) > 4:
        rec["short_F1"] = pfit("F1", yG, idx, SHORT)
        for fm in ("F2", "F3", "F4"):
            rec[f"short_{fm}"] = pfit(fm, yG, idx, SHORT)
    prow.append(rec)
PR = pd.DataFrame(prow).sort_values("n", ascending=False, kind="stable").reset_index(drop=True)
PR = PR.merge(old[["dataset", "doi", "n", "res1", "res2", "short"]], on=["dataset", "doi", "n"], how="left")
PR = PR[(PR.res2_recomputed - PR.res2).abs() < 0.002].reset_index(drop=True)
dev = (PR.res2_recomputed - PR.res2).abs().max()
if len(PR) != 18:
    sys.exit(f"PC2 FAILED: only {len(PR)} of 18 datasets reproduce res2 within 0.002")
print(f"\nPC2: {len(PR)} protein datasets; max |res2 recomputed - stored| = {dev:.4f} kcal/mol")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
cols = ["dataset", "n", "res1", "res2", "res_F2", "res_F3", "res_F4", "res_F5", "res_F6", "res_F7", "nist_F2_med",
        "short", "short_F2", "short_F3", "short_F4"]
print(PR[cols].to_string(index=False, float_format=lambda v: f"{v:.2f}"))
print(f"  median over datasets: res2 {PR.res2.median():.3f}, F2 {PR.res_F2.median():.3f}, F3 {PR.res_F3.median():.3f}, "
      f"F4 {PR.res_F4.median():.3f}, F6 {PR.res_F6.median():.3f}; F2-res2 median {(PR.res_F2-PR.res2).median():+.3f}, "
      f"max {(PR.res_F2-PR.res2).max():+.3f}")
sh = PR.dropna(subset=["short"])
print(f"  shortcut ({len(sh)} datasets): median two-shape {sh.short.median():.3f}, F2 {sh.short_F2.median():.3f}, "
      f"F3 {sh.short_F3.median():.3f}, F4 {sh.short_F4.median():.3f}")
PR.to_csv(OUT / "protein_lattice_strain.tsv", sep="\t", index=False, float_format="%.4f")
summary["protein"] = {"pc2_max_dev": dev, "n": len(PR)}
json.dump(summary, open(OUT / "summary.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("\nwrote", OUT)
