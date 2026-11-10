"""Exploratory (post hoc, the 42 complete ligands that list Y): does yttrium sit on the heavy branch (with the
lanthanides of its size), on the light branch (an ion of its radius that never crossed the gap), or on a
gap-free Shannon curve?  Each hypothesis predicts log K(Y) for a ligand from its 14 lanthanides, with no
yttrium data used.  delta is refitted without the test ligand."""
import sys, re, json
import numpy as np
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa
Y, yY, meta = load_nist(14); n = Y.shape[0]
nN = np.array([int(m.group(1)) if (m := re.search(r"N(\d+)", f)) else 0 for f in meta.formula])
iy = np.where(~np.isnan(yY))[0]
step = (np.arange(14) <= 6).astype(float); grid = np.arange(-0.060, 0.1501, 0.0005)
zz = lambda x: (x - X0) / XS
def curve(b, x):
    return sum(b[k] * zz(x) ** k for k in range(len(b)))
E = {k: [] for k in ("light branch", "heavy branch", "midway", "gap-free cubic", "gap-free quadratic", "twin Dy", "twin Ho", "twin Tb", "mean Dy,Ho")}
deltas = []
for i in iy:
    tr = np.delete(Y, i, 0); S = tr.T @ tr
    f = np.array([coord_objective_complete(R9 + d * step, S, 3) for d in grid]); d = grid[int(np.argmin(f))]; deltas.append(d)
    b = np.linalg.lstsq(poly_basis(R9 + d * step, 3), Y[i], rcond=None)[0]
    E["light branch"].append(yY[i] - curve(b, R9_Y + d)); E["heavy branch"].append(yY[i] - curve(b, R9_Y))
    E["midway"].append(yY[i] - curve(b, R9_Y + d / 2))
    bk = np.linalg.lstsq(poly_basis(R9, 3), Y[i], rcond=None)[0]; E["gap-free cubic"].append(yY[i] - curve(bk, R9_Y))
    bq = np.linalg.lstsq(poly_basis(R9, 2), Y[i], rcond=None)[0]; E["gap-free quadratic"].append(yY[i] - curve(bq, R9_Y))
    for tw in ("Dy", "Ho", "Tb"):
        E[f"twin {tw}"].append(yY[i] - Y[i, LN.index(tw)])
    E["mean Dy,Ho"].append(yY[i] - 0.5 * (Y[i, LN.index("Dy")] + Y[i, LN.index("Ho")]))
E = {k: np.array(v) for k, v in E.items()}
print(f"yttrium, {len(iy)} ligands; gap refitted without the test ligand: {1000*np.min(deltas):.1f}..{1000*np.max(deltas):.1f} mA")
out = {}
for k, v in E.items():
    a = boot_rms(v[:, None] ** 2)
    go = nN[iy] == 0; gn = ~go
    print(f"  {k:19s} RMS {a[0]:.3f} [{a[1]:.3f}, {a[2]:.3f}] | mean error {v.mean():+.3f} | O-only RMS {np.sqrt(np.mean(v[go]**2)):.3f} (mean {v[go].mean():+.3f}, n={go.sum()}) | "
          f"N-containing RMS {np.sqrt(np.mean(v[gn]**2)):.3f} (mean {v[gn].mean():+.3f}, n={gn.sum()})")
    out[k] = dict(rms=a, mean=float(v.mean()), rms_O=float(np.sqrt(np.mean(v[go]**2))), rms_N=float(np.sqrt(np.mean(v[gn]**2))))
for a_, b_ in (("heavy branch", "light branch"), ("gap-free cubic", "light branch"), ("twin Tb", "light branch"), ("twin Dy", "light branch"), ("mean Dy,Ho", "light branch"), ("midway", "light branch")):
    d = boot_diff(E[a_][:, None] ** 2, E[b_][:, None] ** 2); out[f"{a_} - {b_}"] = d
    print(f"  {a_} - {b_}: {d[0]:+.3f} [{d[1]:+.3f}, {d[2]:+.3f}]")
# best fixed fraction of the gap for yttrium (0 = heavy branch, 1 = light branch), all ligands and by class
fr = np.arange(-0.5, 2.51, 0.02)
def best_frac(sel, boot=None):
    idx = iy[sel] if boot is None else iy[sel][boot]
    ss = np.zeros(len(fr))
    for i in idx:
        b = np.linalg.lstsq(poly_basis(R9 + 0.031 * step, 3), Y[i], rcond=None)[0]
        ss += (yY[i] - curve(b, R9_Y + fr * 0.031)) ** 2
    return fr[int(np.argmin(ss))]
rng = np.random.default_rng(42)
for nm, sel in (("all", np.ones(len(iy), bool)), ("O-only", nN[iy] == 0), ("N-containing", nN[iy] >= 1)):
    m = int(sel.sum()); bs = [best_frac(sel, rng.integers(0, m, m)) for _ in range(500)]
    print(f"  best single position for yttrium, {nm:13s} (gap fixed at 31 mA): fraction {best_frac(sel):.2f} [{np.percentile(bs,2.5):.2f}, {np.percentile(bs,97.5):.2f}] of the way from the heavy to the light branch")
    out[f"best_frac_{nm}"] = [float(best_frac(sel)), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
json.dump(out, open(ROOT / "results/yttrium_branch.json", "w"), indent=1, default=float)
