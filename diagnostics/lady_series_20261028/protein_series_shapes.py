"""Do protein/peptide Ln-series selectivities live in the same low-dimensional space as small-molecule chelates?

Basis: uncentred SVD of NIST 46 log K(M) - log K(La) over the 47 ligands with complete La..Lu series.
Each protein dataset (one paper, one condition) is fitted as  log K(M) = c + sum_k b_k V_k(M)  (V_k(La) = 0, c free),
with k = 1 or 2 shapes. Its residual is compared with NIST ligands fitted the same way on the same metal subset,
using a leave-that-ligand-out basis. "Shortcut" test: fit on {La, Nd, Dy, Lu} only, then predict the other metals.
Energies in kcal/mol (1.3643 x log units).
"""
import re
import numpy as np, pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
KC = 1.3643
LN = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
SHORT = ["La", "Nd", "Dy", "Lu"]

nist = pd.read_csv(REPO / "workspaces/ln_surrogate_20261027/surrogate_nist.tsv", sep="\t", dtype={"lig_id": str})
NS = []
for _, r in nist.iterrows():
    k = dict(x.split(":") for x in r.logK_series.split(";"))
    if all(m in k for m in LN):
        NS.append([float(k[m]) - float(k["La"]) for m in LN])
NS = np.array(NS)  # ligands x 14, La column = 0


def basis(X, k):
    _, _, Vt = np.linalg.svd(X[:, 1:], full_matrices=False)
    return np.vstack([np.zeros(k), Vt[:k].T])  # 14 x k, La row = 0


def fit(y, idx, V, fit_idx=None):
    """Least squares y[idx] = c + V[idx] b; returns fitted values at idx (fit uses fit_idx if given)."""
    f = idx if fit_idx is None else fit_idx
    A = np.c_[np.ones(len(f)), V[f]]
    coef = np.linalg.lstsq(A, y[f], rcond=None)[0]
    return np.c_[np.ones(len(idx)), V[idx]] @ coef


def nist_reference(idx, k, shortcut=False):
    out = []
    for i in range(len(NS)):
        V = basis(np.delete(NS, i, 0), k)
        y = NS[i]
        if shortcut:
            fi = [LN.index(m) for m in SHORT]
            rest = [j for j in idx if j not in fi]
            out.append(KC * np.sqrt(np.mean((y[rest] - fit(y, rest, V, fi)) ** 2)))
        else:
            out.append(KC * np.sqrt(np.mean((y[idx] - fit(y, idx, V)) ** 2)))
    return np.median(out), np.percentile(out, 90)


# ---- protein datasets -------------------------------------------------------------------------------------------
d = pd.read_csv(REPO / "workspaces/ln_selectivity_calibration_20261027/protein_ln_affinities.tsv", sep="\t")
d = d[d.metal.isin(LN)].copy()
# Mex-LanM on agarose: relative Keq of the homogeneous two-ion complex, normalised to La; per-ion log K = log10(ratio)/2
m = (d.system == "Mex-LanM (agarose-immobilized)") & d.log_K.isna()
d.loc[m, "log_K"] = [np.log10(float(re.search(r"/Keq_rel\(La\) = ([0-9.eE+-]+)", n).group(1))) / 2 if "/Keq_rel(La) =" in n else
                     (0.0 if mm == "La" else np.nan) for n, mm in zip(d.loc[m, "notes"], d.loc[m, "metal"])]
d = d[d.log_K.notna() & ~d.constant_type.isin(["Km", "EC50_or_K_half_activity"])]
key = ["system", "variant", "doi", "method", "pH", "T_C", "constant_type", "site", "competitor", "buffer_ionic_strength"]
rows = []
for k_, s in d.fillna({c: "-" for c in key}).groupby(key):
    s = s.groupby("metal").log_K.mean()
    if len(s) < 6:
        continue
    y = np.full(14, np.nan); y[[LN.index(x) for x in s.index]] = s.values
    idx = [i for i in range(14) if not np.isnan(y[i])]
    rec = {"dataset": f"{k_[0][:34]} | {k_[1][:22]}", "doi": k_[2], "n": len(idx),
           "span": KC * (np.nanmax(y) - np.nanmin(y))}
    for k in (1, 2):
        V = basis(NS, k)
        rec[f"res{k}"] = KC * np.sqrt(np.mean((y[idx] - fit(y, idx, V)) ** 2))
        rec[f"nist{k}_med"], rec[f"nist{k}_p90"] = nist_reference(idx, k)
    if all(m in s.index for m in SHORT) and len(idx) > 4:
        V = basis(NS, 2); fi = [LN.index(m) for m in SHORT]; rest = [j for j in idx if j not in fi]
        rec["short"] = KC * np.sqrt(np.mean((y[rest] - fit(y, rest, V, fi)) ** 2))
        rec["nist_short_med"], rec["nist_short_p90"] = nist_reference(idx, 2, shortcut=True)
    rows.append(rec)

R = pd.DataFrame(rows).sort_values("n", ascending=False)
pd.set_option("display.width", 250)
print(f"NIST basis from {len(NS)} ligands. Residual RMS in kcal/mol; NIST columns = median / 90th pct for the same metal subset.\n")
for _, r in R.iterrows():
    line = (f"{r.dataset:60s} n={r.n:2d} span {r.span:4.1f} | 1 shape {r.res1:4.2f} (NIST {r.nist1_med:.2f}/{r.nist1_p90:.2f})"
            f" | 2 shapes {r.res2:4.2f} (NIST {r.nist2_med:.2f}/{r.nist2_p90:.2f})")
    if not np.isnan(r.get("short", np.nan)):
        line += f" | La,Nd,Dy,Lu->rest {r.short:4.2f} (NIST {r.nist_short_med:.2f}/{r.nist_short_p90:.2f})"
    print(line)
out = REPO / "workspaces/lady_series_20261028/protein_series_shapes.tsv"
R.to_csv(out, sep="\t", index=False, float_format="%.3f")
print("\nwrote", out)
