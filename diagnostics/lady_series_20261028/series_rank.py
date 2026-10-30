"""How many computed metals pin down the whole Ln series? LOO over NIST ligands with complete La..Lu series."""
import itertools, numpy as np, pandas as pd
KC = 1.3643  # kcal/mol per log unit at 298 K; ddG(M-La) = -KC * dlogK  (negative = M favoured over La)
M = ["La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"]
d = pd.read_csv("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ln_surrogate_20261027/surrogate_nist.tsv", sep="\t", dtype={"lig_id": str})
rows = []
for _, r in d.iterrows():
    k = dict(x.split(":") for x in r.logK_series.split(";"))
    if all(m in k for m in M):
        rows.append([-KC * (float(k[m]) - float(k["La"])) for m in M[1:]] + [r["name"][:28]])
S = pd.DataFrame(rows, columns=M[1:] + ["name"]).set_index("name")
X = S.values; n, p = X.shape
print(f"ligands with complete La..Lu series: {n}")
print("spread across ligands (SD, kcal/mol) of ddG(M-La):", " ".join(f"{m}:{s:.2f}" for m, s in zip(M[1:], X.std(0))))
U, sv, Vt = np.linalg.svd(X, full_matrices=False)  # uncentred: La is the natural zero
ev = sv**2 / (sv**2).sum()
print("uncentred SVD variance:", " ".join(f"{x:.3f}" for x in ev[:5]))
for i in range(3):
    v = Vt[i] * np.sign(Vt[i][-1] or 1)
    print(f"  shape {i+1}:", " ".join(f"{m}:{x:+.2f}" for m, x in zip(M[1:], v)))
def loo(cols):
    idx = [M[1:].index(c) for c in cols]; err = []
    for i in range(n):
        tr = np.arange(n) != i
        B = np.linalg.lstsq(X[tr][:, idx], X[tr], rcond=None)[0]   # through the origin: all-zero site -> all-zero series
        err.append(X[i] - X[i, idx] @ B)
    E = np.array(err)
    return np.sqrt((E**2).mean()), np.abs(E).max(), E
print("\nLOO reconstruction of all 13 metals (relative to La) from k computed metals, kcal/mol")
print(f"  baseline (predict 0 = no selectivity): RMS {np.sqrt((X**2).mean()):.2f}")
for k in (1, 2, 3):
    res = sorted(((loo(c)[0], loo(c)[1], c) for c in itertools.combinations(M[1:], k)))
    print(f"  k={k} best: " + " | ".join(f"{'+'.join(c)} RMS {r:.2f} max {mx:.2f}" for r, mx, c in res[:3]))
for c in (["Dy"], ["Gd", "Lu"], ["Nd", "Dy"], ["Nd", "Gd", "Lu"], ["Sm", "Dy", "Lu"], ["Nd", "Dy", "Lu"]):
    r, mx, E = loo(c)
    print(f"  {'+'.join(c):10s} RMS {r:.2f} max {mx:.2f}   worst ligands: " + ", ".join(f"{S.index[i]} {np.sqrt((E[i]**2).mean()):.2f}" for i in np.argsort(-(E**2).mean(1))[:3]))
# separation factors that matter: Nd/Dy (magnet recycling), La/Ce, Dy/Ho
for c in (["Dy"], ["Gd", "Lu"], ["Nd", "Gd", "Lu"]):
    r, mx, E = loo(c)
    j = lambda m: M[1:].index(m)
    nd_dy = E[:, j("Dy")] - E[:, j("Nd")]
    print(f"  from {'+'.join(c):9s}: error in Dy-Nd selectivity RMS {np.sqrt((nd_dy**2).mean()):.2f} (signal SD {(X[:, j('Dy')] - X[:, j('Nd')]).std():.2f})")
