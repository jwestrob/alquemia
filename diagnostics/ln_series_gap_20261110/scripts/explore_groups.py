"""Exploratory (post hoc, 47 complete ligands): is the Gd|Tb gap the same for O-only and N-containing ligands,
and which branch does yttrium take?"""
import sys, re, json
import numpy as np
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa
Y, yY, meta = load_nist(14); n = Y.shape[0]
nN = np.array([int(m.group(1)) if (m := re.search(r"N(\d+)", f)) else 0 for f in meta.formula])
g = {"O-only": nN == 0, "N-containing": nN >= 1, "all": np.ones(n, bool)}
k0 = 6; step = (np.arange(14) <= k0).astype(float)
grid = np.arange(-0.060, 0.1501, 0.0005)
def best_delta(Ys, deg=3):
    S = Ys.T @ Ys
    f = np.array([coord_objective_complete(R9 + d * step, S, deg) for d in grid]); j = int(np.argmin(f))
    return grid[j], f[j], coord_objective_complete(R9, S, deg)
out = {}
for nm, sel in g.items():
    Ys = Y[sel]; d, f, f0 = best_delta(Ys)
    rng = np.random.default_rng(42); db = [best_delta(Ys[rng.integers(0, len(Ys), len(Ys))])[0] for _ in range(1000)]
    # scale-free version: every ligand's series divided by its SD across metals (so strong ligands do not dominate)
    Yn = Ys / Ys.std(1, keepdims=True); dn, fn, fn0 = best_delta(Yn)
    dbn = [best_delta(Yn[rng.integers(0, len(Yn), len(Yn))])[0] for _ in range(1000)]
    print(f"{nm:13s} n={sel.sum():2d}: gap {1000*d:+.1f} mA [{1000*np.percentile(db,2.5):+.1f}, {1000*np.percentile(db,97.5):+.1f}]; "
          f"residual SS {f0:.3f} -> {f:.3f} ({100*(1-f/f0):.0f}% lower) | each ligand scaled to unit SD: gap {1000*dn:+.1f} "
          f"[{1000*np.percentile(dbn,2.5):+.1f}, {1000*np.percentile(dbn,97.5):+.1f}], residual SS {100*(1-fn/fn0):.0f}% lower")
    out[nm] = dict(gap=float(d), lo=float(np.percentile(db, 2.5)), hi=float(np.percentile(db, 97.5)), gap_scaled=float(dn),
                   lo_scaled=float(np.percentile(dbn, 2.5)), hi_scaled=float(np.percentile(dbn, 97.5)))
# difference between groups
rng = np.random.default_rng(7); dd = []
YO, YN = Y[g["O-only"]], Y[g["N-containing"]]
for _ in range(1000):
    dd.append(best_delta(YN[rng.integers(0, len(YN), len(YN))])[0] - best_delta(YO[rng.integers(0, len(YO), len(YO))])[0])
print(f"gap(N-containing) - gap(O-only): {1000*(out['N-containing']['gap']-out['O-only']['gap']):+.1f} mA [{1000*np.percentile(dd,2.5):+.1f}, {1000*np.percentile(dd,97.5):+.1f}]")
# per-ligand gap (4 cubic numbers + 1 gap on 14 points) for ligands with a usable slope
per = np.full(n, np.nan)
for i in range(n):
    per[i] = best_delta(Y[i:i+1])[0]
print("per-ligand gap (mA), median [IQR]: O-only %.1f [%.1f, %.1f]; N-containing %.1f [%.1f, %.1f]" % (
    1000*np.median(per[g["O-only"]]), 1000*np.percentile(per[g["O-only"]], 25), 1000*np.percentile(per[g["O-only"]], 75),
    1000*np.median(per[g["N-containing"]]), 1000*np.percentile(per[g["N-containing"]], 25), 1000*np.percentile(per[g["N-containing"]], 75)))

# yttrium: where does it sit relative to the two branches?  Per ligand: cubic in the two-branch coordinate
# (gap fixed at the 47-ligand value), then the coordinate at which the curve equals y(Y).
d_all = out["all"]["gap"]
rho = R9 + d_all * step                      # light branch raised by the gap; heavy branch = Shannon
B = poly_basis(rho, 3)
zz = lambda x: (x - X0) / XS
heavy_twin = R9_Y                            # an ion of Y's radius on the heavy branch
light_twin = R9_Y + d_all                    # the same radius on the light branch
print(f"\nyttrium on the two-branch coordinate (gap {1000*d_all:.1f} mA): heavy-branch twin at {heavy_twin:.4f}, light-branch twin at {light_twin:.4f}")
pos = np.full(n, np.nan); frac = np.full(n, np.nan); slope = np.full(n, np.nan)
xs = np.arange(1.02, 1.17, 0.0002)
for i in range(n):
    if np.isnan(yY[i]): continue
    b = np.linalg.lstsq(B, Y[i], rcond=None)[0]
    f = lambda x: sum(b[k] * zz(x) ** k for k in range(4))
    s = sum(k * b[k] * zz(heavy_twin + d_all / 2) ** (k - 1) for k in range(1, 4)) / XS
    slope[i] = s
    if abs(s) < 5: continue
    v = f(xs) - yY[i]
    # root nearest to the midpoint of the two twins
    sg = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    if len(sg) == 0: continue
    r = xs[sg] + 0.0001; j = int(np.argmin(np.abs(r - (heavy_twin + d_all / 2))))
    pos[i] = r[j]; frac[i] = (r[j] - heavy_twin) / d_all
for nm in ("O-only", "N-containing"):
    sel = g[nm] & ~np.isnan(pos); rng = np.random.default_rng(42)
    med = np.median(frac[sel]); bs = [np.median(frac[sel][rng.integers(0, sel.sum(), sel.sum())]) for _ in range(2000)]
    print(f"  {nm:13s} n={sel.sum():2d}: median position {np.median(pos[sel]):.4f}; fraction of the way from the heavy to the light branch "
          f"{med:.2f} [{np.percentile(bs,2.5):.2f}, {np.percentile(bs,97.5):.2f}]")
    out[f"Y_frac_{nm}"] = dict(median=float(med), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)), n=int(sel.sum()))
for v in sorted(set(nN)):
    sel = (nN == v) & ~np.isnan(pos)
    if sel.sum(): print(f"     N atoms = {v}: median fraction {np.median(frac[sel]):.2f} (n={sel.sum()})  values: " + " ".join(f"{x:.2f}" for x in np.sort(frac[sel])))
from scipy.stats import spearmanr
ok = ~np.isnan(frac); print(f"  Spearman fraction vs N atoms: {spearmanr(nN[ok], frac[ok]).statistic:+.3f} (n={ok.sum()})")
json.dump(out, open(ROOT / "results/groups_explore.json", "w"), indent=1)
