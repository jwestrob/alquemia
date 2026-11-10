"""Confirmation tests of ADDENDUM_confirmation.md on the 24 ligands with 10-13 lanthanides (never fitted before),
and on the protein/peptide datasets with at least 9 lanthanides. Fixed inputs only. Lead agent, Claude Fable."""
import sys, re, json
import numpy as np, pandas as pd
from scipy.stats import binomtest
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa
from run_proteins import build
RES = ROOT / "results"; LOG = open(RES / "confirm.log", "w"); SUM = {}
def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
GAP = 0.0310
step = (np.arange(14) <= 6).astype(float)
RS = R9 + GAP * step
RKC = pd.read_csv(RES / "coordinate_cubic.tsv", sep="\t").rho_KC.values
Y47, _, _ = load_nist(14); YBAR = Y47.mean(0)
Yall, yYall, meta = load_nist(10)
inc = np.isnan(Yall).any(1)
Y = Yall[inc]; yY = yYall[inc]; meta = meta[inc].reset_index(drop=True); n = len(Y)
nN = np.array([int(m.group(1)) if (m := re.search(r"N(\d+)", f)) else 0 for f in meta.formula])
p(f"confirmation ligands: {n}; metals present per ligand: {sorted((~np.isnan(Y)).sum(1).tolist())}; with Y: {int((~np.isnan(yY)).sum())}; "
  f"N atoms: {dict(pd.Series(nN).value_counts().sort_index())}")
p("missing metals:", dict(pd.Series([LN[j] for i in range(n) for j in range(14) if np.isnan(Y[i, j])]).value_counts()))
p("ligands:", "; ".join(f"{a[:30]} (N{b})" for a, b in zip(meta.name, nN)))
Z0 = np.zeros(14)
MODS = {"Q": (poly_basis(R9, 2), Z0), "QS": (poly_basis(RS, 2), Z0), "K": (poly_basis(R9, 3), Z0), "KS": (poly_basis(RS, 3), Z0),
        "KC": (poly_basis(RKC, 3), Z0), "KO": (poly_basis(R9, 3), YBAR), "O": (poly_basis(R9, 2), YBAR)}

def lomo(y, B, c):
    idx = np.where(~np.isnan(y))[0]; e = np.full(14, np.nan)
    for m_ in idx:
        fi = idx[idx != m_]; e[m_] = y[m_] - predict(B, c, y, fi)[0][m_]
    return e
E = {k: np.array([lomo(Y[i], *MODS[k]) for i in range(n)]) for k in MODS}
p("\nC1/C2/C4. Leave-one-metal-out RMS over the 24 ligands, kcal/mol")
for k in MODS:
    a = boot_rms(E[k] ** 2); SUM[f"rms_{k}"] = a; p(f"  {k:3s} {a[0]:.4f} [{a[1]:.4f}, {a[2]:.4f}]")
for a_, b_ in (("K", "KS"), ("K", "KC"), ("KS", "KC"), ("K", "KO"), ("Q", "QS"), ("Q", "O"), ("Q", "K")):
    d = boot_diff(E[a_] ** 2, E[b_] ** 2); SUM[f"{a_}-{b_}"] = d
    per_a = np.sqrt(np.nanmean(E[a_] ** 2, 1)); per_b = np.sqrt(np.nanmean(E[b_] ** 2, 1))
    p(f"  {a_:2s} - {b_:2s}: {d[0]:+.4f} [{d[1]:+.4f}, {d[2]:+.4f}]; second model better for {int((per_b < per_a).sum())}/{n} ligands "
      f"(sign test p = {binomtest(int((per_b < per_a).sum()), n, 0.5).pvalue:.3f})")
    SUM[f"{a_}-{b_}_n_better"] = int((per_b < per_a).sum())
# interior-only (held-out metal not the first or last available)
def interior(Em):
    out = Em.copy()
    for i in range(n):
        idx = np.where(~np.isnan(Y[i]))[0]; out[i, idx[0]] = np.nan; out[i, idx[-1]] = np.nan
    return out
d = boot_diff(interior(E["K"]) ** 2, interior(E["KS"]) ** 2); p(f"  interior metals only, K - KS: {d[0]:+.4f} [{d[1]:+.4f}, {d[2]:+.4f}]"); SUM["K-KS_interior"] = d

# C3. independent estimate of the gap and of the break position (missing data: per-ligand least squares)
grid = np.arange(-0.060, 0.1501, 0.0005)
def rss(rho, rows, deg=3):
    B = poly_basis(rho, deg); t = 0.0
    for y in rows:
        ok = ~np.isnan(y); b = np.linalg.lstsq(B[ok], y[ok], rcond=None)[0]; t += float(((y[ok] - B[ok] @ b) ** 2).sum())
    return t
def per_lig_curve(rows, k):
    return np.array([[rss(R9 + dd * (np.arange(14) <= k), [y]) for dd in grid] for y in rows])   # ligands x grid
p("\nC3. Gap fitted on the 24 ligands alone (cubic per ligand)")
curves = {k: per_lig_curve(Y, k) for k in range(1, 13)}
f0 = rss(R9, Y)
tot = {k: curves[k].sum(0) for k in curves}
for k in range(1, 13):
    j = int(np.argmin(tot[k])); p(f"  break {LN[k]}|{LN[k+1]}: best delta {1000*grid[j]:+6.1f} mA, residual SS {tot[k][j]:.4f} (no break {f0:.4f})")
bestk = min(tot, key=lambda k: tot[k].min()); j6 = int(np.argmin(tot[6]))
rng = np.random.default_rng(42); db = []; wins = 0
for _ in range(2000):
    i = rng.integers(0, n, n)
    tk = {k: curves[k][i].sum(0) for k in curves}
    wins += min(tk, key=lambda k: tk[k].min()) == 6; db.append(grid[int(np.argmin(tk[6]))])
p(f"  Gd|Tb gap on the 24: {1000*grid[j6]:+.1f} mA [{1000*np.percentile(db,2.5):+.1f}, {1000*np.percentile(db,97.5):+.1f}]; "
  f"best break position: {LN[bestk]}|{LN[bestk+1]}; Gd|Tb best in {wins}/2000 resamples")
SUM["gap_24"] = dict(gap=float(grid[j6]), lo=float(np.percentile(db, 2.5)), hi=float(np.percentile(db, 97.5)),
                     best_break=f"{LN[bestk]}|{LN[bestk+1]}", gdtb_wins=int(wins))
for nm, sel in (("O-only", nN == 0), ("N-containing", nN >= 1)):
    if sel.sum() >= 3:
        t = curves[6][sel].sum(0); jj = int(np.argmin(t)); idx = np.where(sel)[0]
        bb = [grid[int(np.argmin(curves[6][rng.choice(idx, len(idx))].sum(0)))] for _ in range(2000)]
        p(f"    {nm} (n={sel.sum()}): gap {1000*grid[jj]:+.1f} mA [{1000*np.percentile(bb,2.5):+.1f}, {1000*np.percentile(bb,97.5):+.1f}]")
        SUM[f"gap_24_{nm}"] = [float(grid[jj]), float(np.percentile(bb, 2.5)), float(np.percentile(bb, 97.5))]

# C5. yttrium
zz = lambda x: (x - X0) / XS
iy = np.where(~np.isnan(yY))[0]
p(f"\nC5. Yttrium on the confirmation ligands (n = {len(iy)}; O-only {int((nN[iy]==0).sum())}, N-containing {int((nN[iy]>=1).sum())})")
EY = {k: [] for k in ("light", "heavy", "gap-free cubic", "twin Dy", "twin Ho", "twin Tb")}
keep = []
for i in iy:
    ok = ~np.isnan(Y[i]); b = np.linalg.lstsq(poly_basis(RS, 3)[ok], Y[i][ok], rcond=None)[0]
    cv = lambda x, b=b: sum(b[k] * zz(x) ** k for k in range(4))
    bk = np.linalg.lstsq(poly_basis(R9, 3)[ok], Y[i][ok], rcond=None)[0]
    EY["light"].append(yY[i] - cv(R9_Y + GAP)); EY["heavy"].append(yY[i] - cv(R9_Y))
    EY["gap-free cubic"].append(yY[i] - sum(bk[k] * zz(R9_Y) ** k for k in range(4)))
    for tw in ("Dy", "Ho", "Tb"):
        EY[f"twin {tw}"].append(yY[i] - Y[i, LN.index(tw)])
    keep.append(i)
EY = {k: np.array(v) for k, v in EY.items()}; cls = nN[iy]
for k, v in EY.items():
    line = f"  {k:15s} all RMS {np.sqrt(np.nanmean(v**2)):.3f}"
    for nm, sel in (("O-only", cls == 0), ("N-containing", cls >= 1)):
        if sel.sum(): line += f" | {nm} RMS {np.sqrt(np.nanmean(v[sel]**2)):.3f} (mean {np.nanmean(v[sel]):+.3f}, n={int(np.sum(sel & ~np.isnan(v)))})"
    p(line); SUM[f"Y_{k}"] = dict(all=float(np.sqrt(np.nanmean(v**2))), O=float(np.sqrt(np.nanmean(v[cls == 0]**2))) if (cls == 0).sum() else None,
                                  N=float(np.sqrt(np.nanmean(v[cls >= 1]**2))) if (cls >= 1).sum() else None)
for nm, sel in (("all", np.ones(len(iy), bool)), ("O-only", cls == 0), ("N-containing", cls >= 1)):
    if sel.sum() >= 3:
        d = boot_diff(EY["heavy"][sel][:, None] ** 2, EY["light"][sel][:, None] ** 2); SUM[f"Y_heavy-light_{nm}"] = d
        p(f"  RMS(heavy) - RMS(light), {nm}: {d[0]:+.3f} [{d[1]:+.3f}, {d[2]:+.3f}] (positive = light branch fits better)")
fr = np.arange(-1.0, 3.01, 0.02)
def best_frac(idx):
    ss = np.zeros(len(fr))
    for i in idx:
        ok = ~np.isnan(Y[i]); b = np.linalg.lstsq(poly_basis(RS, 3)[ok], Y[i][ok], rcond=None)[0]
        ss += (yY[i] - sum(b[k] * zz(R9_Y + fr * GAP) ** k for k in range(4))) ** 2
    return fr[int(np.argmin(ss))]
for nm, sel in (("all", np.ones(len(iy), bool)), ("O-only", cls == 0), ("N-containing", cls >= 1)):
    idx = iy[sel]
    if len(idx) >= 3:
        bs = [best_frac(rng.choice(idx, len(idx))) for _ in range(500)]
        p(f"  best single fraction of the gap for yttrium, {nm}: {best_frac(idx):.2f} [{np.percentile(bs,2.5):.2f}, {np.percentile(bs,97.5):.2f}]")
        SUM[f"Y_best_frac_{nm}"] = [float(best_frac(idx)), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]

# C6. proteins and peptides with at least 9 lanthanides
p("\nC6. Protein and peptide datasets with at least 9 lanthanides: leave-one-metal-out RMS, kcal/mol")
rows = []
for ds in build():
    y = ds["y"]
    if (~np.isnan(y)).sum() < 9: continue
    rec = dict(dataset=ds["dataset"], doi=ds["doi"], n=int((~np.isnan(y)).sum()))
    for k in ("Q", "QS", "K", "KS", "KC"):
        rec[k] = float(np.sqrt(np.nanmean(lomo(y, *MODS[k]) ** 2)))
    # the dataset's own best gap (cubic and quadratic), descriptive
    c = np.array([rss(R9 + dd * step, [y]) for dd in grid]); rec["own_gap_mA_cubic"] = 1000 * grid[int(np.argmin(c))]
    c2 = np.array([rss(R9 + dd * step, [y], 2) for dd in grid]); rec["own_gap_mA_quadratic"] = 1000 * grid[int(np.argmin(c2))]
    rows.append(rec)
R = pd.DataFrame(rows).sort_values("n", ascending=False); R.to_csv(RES / "confirm_proteins.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 250); p(R.round(3).to_string(index=False))
for a_, b_ in (("K", "KS"), ("Q", "QS"), ("K", "KC")):
    d = (R[a_] - R[b_]).values; pos = int((d > 0).sum())
    p(f"  {a_} - {b_}: median {np.median(d):+.4f}; second better in {pos}/{len(d)}; sign test p = {binomtest(pos, len(d), 0.5).pvalue:.3f}")
    SUM[f"protein_{a_}-{b_}"] = dict(median=float(np.median(d)), n_better=pos, n=len(d), p=float(binomtest(pos, len(d), 0.5).pvalue))
json.dump(SUM, open(RES / "confirm_summary.json", "w"), indent=1, default=float)
