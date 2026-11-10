"""One-number version of the coordinate: Shannon CN9 radius with the light branch (La..Gd) shifted by delta
relative to the heavy branch (Tb..Lu). delta fitted on the 47 complete ligands under the cubic model.
Also: the same with the break placed between every other pair of neighbours (which break position fits best?).
Uses only the 47 complete-series ligands."""
import sys, json
import numpy as np
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa
Y, yY, meta = load_nist(14); n = Y.shape[0]; A = np.arange(14)
S = Y.T @ Y
def obj(delta, k, deg=3, Sm=S):
    rho = R9 + delta * (np.arange(14) <= k)
    return coord_objective_complete(rho, Sm, deg)
grid = np.arange(-0.060, 0.0601, 0.0005)
print("break after | best delta (mA) | objective (cubic)   [no break: %.5f; free coordinate KC: 4.87270]" % obj(0, 6))
best = {}
for k in range(1, 13):
    f = np.array([obj(d, k) for d in grid]); j = int(np.argmin(f)); best[k] = (grid[j], f[j])
    print(f"  {LN[k]:2s}|{LN[k+1]:2s}   {1000*grid[j]:+6.1f}   {f[j]:.5f}")
k0 = 6
d0 = best[k0][0]
# bootstrap interval for delta at the Gd|Tb break, and how often Gd|Tb is the best break position
rng = np.random.default_rng(42); db = []; wins = 0; nb = 1000
for _ in range(nb):
    i = rng.integers(0, n, n); Sb = Y[i].T @ Y[i]
    fk = {k: min(obj(d, k, 3, Sb) for d in grid[::4]) for k in range(1, 13)}
    wins += min(fk, key=fk.get) == k0
    f = np.array([obj(d, k0, 3, Sb) for d in grid]); db.append(grid[int(np.argmin(f))])
print(f"Gd|Tb break: delta = {1000*d0:+.1f} mA, bootstrap 95% [{1000*np.percentile(db,2.5):+.1f}, {1000*np.percentile(db,97.5):+.1f}]; "
      f"Gd|Tb is the best single break position in {wins}/{nb} resamples")
# quadratic-model version for the record
fq = np.array([obj(d, k0, 2) for d in grid]); print(f"quadratic model, Gd|Tb break: delta = {1000*grid[int(np.argmin(fq))]:+.1f} mA; objective {fq.min():.5f} (no break {obj(0,k0,2):.5f})")
# LOLO + LOMO score of KS (delta refitted without the held-out ligand) beside K
E = {m: np.zeros((n, 14)) for m in ("K", "KS")}; dl = []
for i in range(n):
    tr = np.delete(Y, i, 0); St = tr.T @ tr
    f = np.array([obj(d, k0, 3, St) for d in grid]); di = grid[int(np.argmin(f))]; dl.append(di)
    mods = {"K": poly_basis(R9, 3), "KS": poly_basis(R9 + di * (np.arange(14) <= k0), 3)}
    for m, B in mods.items():
        for mm in A:
            fi = A[A != mm]; E[m][i, mm] = Y[i, mm] - predict(B, np.zeros(14), Y[i], fi)[0][mm]
a = boot_rms(E["K"] ** 2); b = boot_rms(E["KS"] ** 2); d = boot_diff(E["K"] ** 2, E["KS"] ** 2)
print(f"S_LOMO: K {a[0]:.4f} [{a[1]:.4f}, {a[2]:.4f}]; KS {b[0]:.4f} [{b[1]:.4f}, {b[2]:.4f}]; K - KS {d[0]:+.4f} [{d[1]:+.4f}, {d[2]:+.4f}]; "
      f"LOLO delta range {1000*min(dl):+.1f}..{1000*max(dl):+.1f} mA")
np.save(ROOT / "results/lomo_KS.npy", E["KS"])
json.dump(dict(delta_A=float(d0), boot_lo=float(np.percentile(db, 2.5)), boot_hi=float(np.percentile(db, 97.5)),
               best_break_wins=int(wins), n_boot=nb, by_break={f"{LN[k]}|{LN[k+1]}": [float(v[0]), float(v[1])] for k, v in best.items()},
               S_LOMO_K=a, S_LOMO_KS=b, K_minus_KS=d), open(ROOT / "results/step_fit.json", "w"), indent=1)
