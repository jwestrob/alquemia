"""Full-data fit of the shared coordinate under the cubic-per-ligand model (KC), with bootstrap intervals.
Uses only the 47 complete-series ligands. Written before the confirmation test (ADDENDUM_confirmation.md)."""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa
Y, yY, meta = load_nist(14); n = Y.shape[0]
rho, _, al = fit_coordinate(Y, offset=False, deg=3, return_all=True)
f = np.array([a[0] for a in al]); d = np.array([np.abs(a[2] - al[0][2]).max() for a in al])
print(f"KC: best objective {f[0]:.5f}; starts at best (1e-6): {(f - f[0] < 1e-6).sum()}/{len(f)}; distinct optima: {(d > 1e-4).sum()};"
      f" objective range {f.min():.5f}..{f.max():.5f}; monotonic {bool(np.all(np.diff(rho) < 0))}")
for a in al[:8]:
    print(f"   start {a[1]:6s} objective {a[0]:.5f} max|d rho| {1000*np.abs(a[2]-al[0][2]).max():.2f} mA")
rng = np.random.default_rng(42); nb = 1000; B = np.zeros((nb, 14))
for b in range(nb):
    i = rng.integers(0, n, n)
    B[b] = fit_coordinate(Y[i], False, deg=3, starts={"cn9": np.zeros(12), "full": (rho[1:-1] - R9[1:-1]) / 0.01})[0]
t = pd.DataFrame(dict(metal=LN, q4f=Q4F, r9=R9, rho_KC=rho, lo=np.percentile(B, 2.5, 0), hi=np.percentile(B, 97.5, 0),
                      d_mA=1000 * (rho - R9), d_lo=1000 * (np.percentile(B, 2.5, 0) - R9), d_hi=1000 * (np.percentile(B, 97.5, 0) - R9)))
t.to_csv(ROOT / "results/coordinate_cubic.tsv", sep="\t", index=False, float_format="%.5f")
np.save(ROOT / "results/bootstrap_rho_KC.npy", B)
print(t.round(4).to_string(index=False))
# local displacement of each metal relative to the line through its two neighbours (in mA), with bootstrap interval
def bump(r):
    out = []
    for j in range(1, 13):
        w = (R9[j - 1] - R9[j]) / (R9[j - 1] - R9[j + 1])
        out.append((r[..., j] - ((1 - w) * r[..., j - 1] + w * r[..., j + 1])) - (R9[j] - ((1 - w) * R9[j - 1] + w * R9[j + 1])))
    return 1000 * np.array(out)
bb = bump(B); b0 = bump(rho)
print("displacement relative to neighbours, beyond Shannon (mA):")
for j in range(12):
    print(f"  {LN[j+1]}: {b0[j]:+.1f} [{np.percentile(bb[j], 2.5):+.1f}, {np.percentile(bb[j], 97.5):+.1f}]")
