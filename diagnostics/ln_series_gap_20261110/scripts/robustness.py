"""Post hoc robustness checks for the Gd|Tb gap. Lead agent, Claude Fable.
A. Other base scales (Shannon CN8, CN6 as recalled, atomic number).
B. Shape of the half-shell feature: step after Gd against the spin-pairing tent (cusp at Gd) and others,
   each as a one-number distortion of the coordinate; fitted on the 47, scored on the 47 (LOLO + LOMO) and on the 24.
C. Protein and peptide datasets: each dataset's own gap with an F-based 95% interval; pooled gap.
D. Interior-only score on the 47."""
import sys, re, json
import numpy as np, pandas as pd
from scipy.stats import f as fdist
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa
from run_proteins import build
RES = ROOT / "results"; LOG = open(RES / "robustness.log", "w"); SUM = {}
def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
Y, yY, meta = load_nist(14); n = Y.shape[0]; A = np.arange(14); S = Y.T @ Y
Yall, _, _ = load_nist(10); Y24 = Yall[np.isnan(Yall).any(1)]
R6 = np.array([1.032, 1.010, 0.990, 0.983, 0.958, 0.947, 0.938, 0.923, 0.912, 0.901, 0.890, 0.880, 0.868, 0.861])  # recalled
step = (np.arange(14) <= 6).astype(float)

def basis_x(x, deg=3):
    z = (x - x.mean()) / x.std(); return np.vstack([z ** k for k in range(deg + 1)]).T
def obj(x, Sm=S, deg=3):
    Qm, _ = np.linalg.qr(basis_x(x, deg)); return float(np.trace(Sm) - np.trace(Qm.T @ Sm @ Qm))
def rss_rows(x, rows, deg=3):
    B = basis_x(x, deg); t = 0.0
    for y in rows:
        ok = ~np.isnan(y); b = np.linalg.lstsq(B[ok], y[ok], rcond=None)[0]; t += float(((y[ok] - B[ok] @ b) ** 2).sum())
    return t

p("A. Gap between Gd and Tb on other base scales (cubic per ligand, 47 ligands)")
for nm, base in (("Shannon CN9", R9), ("Shannon CN8", R8), ("Shannon CN6 (recalled)", R6), ("atomic number", -ZAT.astype(float))):
    sp = np.abs(np.diff(base)).mean()           # mean spacing of adjacent listed metals
    grid = np.arange(-6, 12.01, 0.05) * sp
    f = np.array([obj(base + d * step) for d in grid]); j = int(np.argmin(f)); f0 = obj(base)
    rng = np.random.default_rng(42); db = []
    for _ in range(500):
        i = rng.integers(0, n, n); Sb = Y[i].T @ Y[i]
        db.append(grid[int(np.argmin([obj(base + d * step, Sb) for d in grid[::2]])) * 2])
    gdtb = abs(base[6] - base[7])
    p(f"  {nm:24s}: no-break SS {f0:.3f}; with gap {f[j]:.3f}; gap = {grid[j]/sp:+.2f} mean element steps "
      f"[{np.percentile(db,2.5)/sp:+.2f}, {np.percentile(db,97.5)/sp:+.2f}]"
      + (f" = {1000*grid[j]:+.1f} mA; scale's own Gd-Tb step {1000*gdtb:.1f} mA" if nm.startswith("Shannon") else ""))
    SUM[f"A_{nm}"] = dict(ss0=f0, ss=float(f[j]), gap_steps=float(grid[j] / sp), lo=float(np.percentile(db, 2.5) / sp), hi=float(np.percentile(db, 97.5) / sp))

p("\nB. One-number distortions of the coordinate, rho = r9 + kappa * shape (cubic per ligand)")
tent = T1 / T1.max(); tentz = T1 * (Q4F + 25) / (T1 * (Q4F + 25)).max(); w = T3 / T3.max()
tri = 1 - np.abs(Q4F - 7) / 7.0
shapes = {"step after Gd": step, "step after Eu": (np.arange(14) <= 5).astype(float), "step after Tb": (np.arange(14) <= 7).astype(float),
          "spin-pairing tent (E1)": tent, "tent with (q+25)": tentz, "triangle cusp at Gd": tri, "E3 shape": w}
kg = np.arange(-0.100, 0.1001, 0.0005)
fit = {}
for nm, g in shapes.items():
    f = np.array([obj(R9 + k * g) for k in kg]); j = int(np.argmin(f)); fit[nm] = kg[j]
    p(f"  {nm:24s}: kappa {1000*kg[j]:+6.1f} mA; residual SS {f[j]:.4f} (none: {obj(R9):.4f}; free coordinate: 4.8727)")
    SUM[f"B_{nm}"] = dict(kappa=float(kg[j]), ss=float(f[j]))
# two-shape combinations
def fit2(g1, g2, Sm=S):
    best = (1e9, 0, 0)
    for k1 in kg[::4]:
        for k2 in kg[::4]:
            v = obj(R9 + k1 * g1 + k2 * g2, Sm)
            if v < best[0]: best = (v, k1, k2)
    return best
b2 = fit2(step, tent); p(f"  step after Gd + tent      : kappas {1000*b2[1]:+.1f}, {1000*b2[2]:+.1f} mA; residual SS {b2[0]:.4f}")
b3 = fit2(tent, w); p(f"  tent + E3 shape           : kappas {1000*b3[1]:+.1f}, {1000*b3[2]:+.1f} mA; residual SS {b3[0]:.4f}")
SUM["B_step+tent"] = dict(k_step=float(b2[1]), k_tent=float(b2[2]), ss=float(b2[0])); SUM["B_tent+E3"] = dict(k1=float(b3[1]), k2=float(b3[2]), ss=float(b3[0]))
# LOLO + LOMO on the 47 with kappa refitted without the held-out ligand; and on the 24 with kappa fixed from the 47
E47 = {nm: np.zeros((n, 14)) for nm in list(shapes) + ["none"]}
for i in range(n):
    tr = np.delete(Y, i, 0); St = tr.T @ tr
    ks = {nm: kg[int(np.argmin([obj(R9 + k * g, St) for k in kg[::2]])) * 2] for nm, g in shapes.items()}
    for nm in E47:
        x = R9 if nm == "none" else R9 + ks[nm] * shapes[nm]; B = basis_x(x)
        for m_ in A:
            fi = A[A != m_]; E47[nm][i, m_] = Y[i, m_] - predict(B, np.zeros(14), Y[i], fi)[0][m_]
def lomo24(x):
    B = basis_x(x); out = np.full(Y24.shape, np.nan)
    for i, y in enumerate(Y24):
        idx = np.where(~np.isnan(y))[0]
        for m_ in idx:
            fi = idx[idx != m_]; out[i, m_] = y[m_] - predict(B, np.zeros(14), y, fi)[0][m_]
    return out
E24 = {nm: lomo24(R9 if nm == "none" else R9 + fit[nm] * shapes[nm]) for nm in E47}
p("  cross-validated RMS (kcal/mol): 47 ligands (LOLO + LOMO) | 24 confirmation ligands (kappa fixed from the 47)")
for nm in E47:
    a = boot_rms(E47[nm] ** 2); b = boot_rms(E24[nm] ** 2)
    ai = boot_rms(E47[nm][:, 1:-1] ** 2)
    p(f"    {nm:24s}: {a[0]:.4f} [{a[1]:.4f}, {a[2]:.4f}] (interior {ai[0]:.4f}) | {b[0]:.4f} [{b[1]:.4f}, {b[2]:.4f}]")
    SUM[f"B_cv_{nm}"] = dict(s47=a, s47_interior=ai, s24=b)
for nm in ("spin-pairing tent (E1)", "tent with (q+25)", "triangle cusp at Gd", "step after Eu", "step after Tb"):
    d = boot_diff(E47[nm] ** 2, E47["step after Gd"] ** 2); d2 = boot_diff(E24[nm] ** 2, E24["step after Gd"] ** 2)
    p(f"    {nm} minus step after Gd: 47: {d[0]:+.4f} [{d[1]:+.4f}, {d[2]:+.4f}] | 24: {d2[0]:+.4f} [{d2[1]:+.4f}, {d2[2]:+.4f}]")
    SUM[f"B_diff_{nm}"] = dict(d47=d, d24=d2)
d = boot_diff(E47["none"][:, 1:-1] ** 2, E47["step after Gd"][:, 1:-1] ** 2); p(f"  D. interior 12 metals, 47 ligands: cubic minus cubic with gap {d[0]:+.4f} [{d[1]:+.4f}, {d[2]:+.4f}]"); SUM["D_interior"] = d

p("\nC. Protein and peptide datasets (at least 9 lanthanides): own Gd|Tb gap with F-based 95% interval (cubic; quadratic in brackets)")
grid = np.arange(-0.150, 0.1501, 0.0005)
rows = []; curves3 = []; curves2 = []; names = []
for ds in build():
    y = ds["y"]; m = int((~np.isnan(y)).sum())
    if m < 9: continue
    rec = dict(dataset=ds["dataset"][:44], doi=ds["doi"], n=m)
    for deg, store in ((3, curves3), (2, curves2)):
        c = np.array([rss_rows(R9 + dd * step, [y], deg) for dd in grid]); store.append(c)
        j = int(np.argmin(c)); dof = m - (deg + 2)
        thr = c[j] * (1 + fdist.ppf(0.95, 1, dof) / dof); okk = np.where(c <= thr)[0]
        rec[f"gap{deg}"] = 1000 * grid[j]; rec[f"lo{deg}"] = 1000 * grid[okk[0]]; rec[f"hi{deg}"] = 1000 * grid[okk[-1]]
        rec[f"rms{deg}_nogap"] = np.sqrt(rss_rows(R9, [y], deg) / m); rec[f"rms{deg}_gap31"] = np.sqrt(rss_rows(R9 + 0.031 * step, [y], deg) / m)
    rows.append(rec); names.append(ds["system"])
R = pd.DataFrame(rows); R.to_csv(RES / "protein_gap.tsv", sep="\t", index=False, float_format="%.3f")
pd.set_option("display.width", 260); p(R.round(2).to_string(index=False))
C3 = np.array(curves3); C2 = np.array(curves2); rng = np.random.default_rng(42)
def pooled(C, idx, norm=True):
    Cn = C[idx] / C[idx].min(1, keepdims=True) if norm else C[idx]     # each dataset weighted by its own best residual
    return grid[int(np.argmin(Cn.sum(0)))]
groups = {"all 12": np.arange(len(R)), "lanmodulin (1)": np.array([i for i, s in enumerate(names) if "LanM" in s]),
          "LBT peptides (6)": np.array([i for i, s in enumerate(names) if "LBT" in s]),
          "not LBT (6)": np.array([i for i, s in enumerate(names) if "LBT" not in s])}
for nm, idx in groups.items():
    if len(idx) == 0: continue
    for deg, C in ((3, C3), (2, C2)):
        g0 = pooled(C, idx)
        if len(idx) >= 4:
            bs = [pooled(C, rng.choice(idx, len(idx))) for _ in range(2000)]
            p(f"  pooled gap, {nm:17s}, degree {deg}: {1000*g0:+.1f} mA [{1000*np.percentile(bs,2.5):+.1f}, {1000*np.percentile(bs,97.5):+.1f}] (bootstrap over datasets)")
            SUM[f"C_pooled_{nm}_deg{deg}"] = [float(g0), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
        else:
            p(f"  pooled gap, {nm:17s}, degree {deg}: {1000*g0:+.1f} mA")
# the same F-based interval for each NIST ligand: how many of the 47 + 24 individually exclude zero / include 31?
def own(y, deg=3):
    m = int((~np.isnan(y)).sum()); c = np.array([rss_rows(R9 + dd * step, [y], deg) for dd in grid]); j = int(np.argmin(c)); dof = m - (deg + 2)
    thr = c[j] * (1 + fdist.ppf(0.95, 1, dof) / dof); okk = np.where(c <= thr)[0]
    return grid[j], grid[okk[0]], grid[okk[-1]]
o = np.array([own(y) for y in np.vstack([Y, Y24])])
p(f"  for comparison, the 71 NIST ligands one at a time (cubic): median own gap {1000*np.median(o[:,0]):+.1f} mA; "
  f"interval excludes 0 on the positive side for {int((o[:,1] > 0).sum())}, on the negative side for {int((o[:,2] < 0).sum())}; includes +31 for {int(((o[:,1] <= 0.031) & (o[:,2] >= 0.031)).sum())}")
SUM["C_nist_own"] = dict(median=float(np.median(o[:, 0])), n_pos=int((o[:, 1] > 0).sum()), n_neg=int((o[:, 2] < 0).sum()),
                         n_incl31=int(((o[:, 1] <= 0.031) & (o[:, 2] >= 0.031)).sum()))
json.dump(SUM, open(RES / "robustness_summary.json", "w"), indent=1, default=float)
