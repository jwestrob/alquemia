"""Shared code for the lanthanide-series coordinate analysis (PLAN.md, frozen 2026-11-10). Lead agent, Claude Fable."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
KC = 1.3643  # kcal/mol per log unit at 298 K
LN = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
Q4F = np.array([0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])          # 4f electron count
ZAT = Q4F + 57
R9 = np.array([1.216, 1.196, 1.179, 1.163, 1.132, 1.120, 1.107, 1.095, 1.083, 1.072, 1.062, 1.052, 1.042, 1.032])
R8 = np.array([1.160, 1.143, 1.126, 1.109, 1.079, 1.066, 1.053, 1.040, 1.027, 1.015, 1.004, 0.994, 0.985, 0.977])
R9_Y, R8_Y = 1.075, 1.019
_n = np.array([0, 0, 1, 3, 6, 10, 15, 21, 15, 10, 6, 3, 1, 0, 0], float)
_m = np.array([0, 0, 9, 21, 21, 9, 0, 0, 0, 9, 21, 21, 9, 0, 0], float)
T1 = (9.0 / 13.0) * _n[Q4F]
T3 = _m[Q4F]
X0, XS = 1.124, 0.06  # centring and scale for numerical conditioning of polynomial bases


def load_nist(min_metals=14):
    d = pd.read_csv(ROOT / "data/surrogate_nist.tsv", sep="\t", dtype={"lig_id": str})
    Y, yy, keep = [], [], []
    for i, r in d.iterrows():
        k = dict(x.split(":") for x in r.logK_series.split(";"))
        if sum(m in k for m in LN) >= min_metals:
            Y.append([-KC * float(k[m]) if m in k else np.nan for m in LN])
            yy.append(-KC * float(k["Y"]) if "Y" in k else np.nan)
            keep.append(i)
    meta = d.loc[keep].reset_index(drop=True)
    meta["n_donor_N"] = meta.aminoacid_N + meta.amine_N + meta.pyridine_N
    return np.array(Y), np.array(yy), meta


def poly_basis(x, deg=2):
    z = (np.asarray(x, float) - X0) / XS
    return np.vstack([z ** k for k in range(deg + 1)]).T


def projector(B):
    Qm, _ = np.linalg.qr(B)
    return Qm @ Qm.T


def coord_objective_complete(rho, S, deg=2):
    Qm, _ = np.linalg.qr(poly_basis(rho, deg))
    return float(np.trace(S) - np.trace(Qm.T @ S @ Qm))


def coord_objective_missing(rho, Y, C):
    """Sum over ligands of the residual SS of (y - c) on [1, rho, rho^2] over the observed metals."""
    B = poly_basis(rho); tot = 0.0
    for y in Y:
        ok = ~np.isnan(y)
        b = B[ok]; v = y[ok] - C[ok]
        beta = np.linalg.lstsq(b, v, rcond=None)[0]
        tot += float(((v - b @ beta) ** 2).sum())
    return tot


def _pin(u, lo=R9[0], hi=R9[-1]):
    return np.concatenate([[lo], R9[1:-1] + 0.01 * u, [hi]])


def start_set(Y_centred_rows, n_jitter=20, seed=42):
    """Starting coordinates (as u = (rho - r9)/0.01 for the 12 free metals)."""
    starts = {"cn9": np.zeros(12)}
    r8s = R9[0] + (R8 - R8[0]) * (R9[-1] - R9[0]) / (R8[-1] - R8[0])
    starts["cn8"] = (r8s[1:-1] - R9[1:-1]) / 0.01
    zs = R9[0] + (ZAT - ZAT[0]) * (R9[-1] - R9[0]) / (ZAT[-1] - ZAT[0])
    starts["linZ"] = (zs[1:-1] - R9[1:-1]) / 0.01
    _, _, Vt = np.linalg.svd(Y_centred_rows, full_matrices=False)
    v = Vt[0]
    if abs(v[-1] - v[0]) > 1e-9:
        ps = R9[0] + (v - v[0]) * (R9[-1] - R9[0]) / (v[-1] - v[0])
        starts["pc1"] = (ps[1:-1] - R9[1:-1]) / 0.01
    rng = np.random.default_rng(seed)
    for j in range(n_jitter):
        starts[f"jit{j}"] = rng.normal(0, 1.0, 12)
    return starts


def fit_coordinate(Y, offset=False, n_jitter=20, seed=42, starts=None, return_all=False, deg=2):
    """Fit rho (ends pinned to Shannon CN9). offset=True fits about the ligand-mean series (model OC).
    Y may contain NaN (then the mean series is the column nanmean and the slower objective is used)."""
    complete = not np.isnan(Y).any()
    c = np.nanmean(Y, axis=0) if offset else np.zeros(Y.shape[1])
    if complete:
        Z = Y - c
        S = Z.T @ Z
        f = lambda u: coord_objective_complete(_pin(u), S, deg)
        rows = Z - Z.mean(1, keepdims=True)
    else:
        f = lambda u: coord_objective_missing(_pin(u), Y, c)
        rows = np.where(np.isnan(Y), 0.0, Y - c)
        rows = rows - rows.mean(1, keepdims=True)
    st = starts if starts is not None else start_set(rows, n_jitter=n_jitter, seed=seed)
    res = []
    for name, u0 in st.items():
        o = minimize(f, u0, method="L-BFGS-B", options=dict(maxiter=5000, ftol=1e-14, gtol=1e-9, maxfun=200000))
        res.append((o.fun, name, _pin(o.x), o.success))
    res.sort(key=lambda t: t[0])
    best = res[0]
    if return_all:
        return best[2], c, res
    return best[2], c


def svd_shapes(Y, offset=False):
    c = Y.mean(0) if offset else np.zeros(Y.shape[1])
    Z = Y - c
    Z = Z - Z.mean(1, keepdims=True)
    _, _, Vt = np.linalg.svd(Z, full_matrices=False)
    return np.c_[np.ones(Y.shape[1]), Vt[0], Vt[1]], c


def predict(B, c, y, fit_idx):
    """Per-ligand least squares of (y - c) on basis B over fit_idx; returns the full predicted series."""
    beta = np.linalg.lstsq(B[fit_idx], (y - c)[fit_idx], rcond=None)[0]
    return B @ beta + c, beta


def boot_rms(E2, n_boot=2000, seed=42):
    """E2: ligands x k squared errors. Returns point RMS and 95% percentile interval, resampling ligands."""
    rng = np.random.default_rng(seed)
    n = E2.shape[0]; per = np.nanmean(E2, axis=1)
    idx = rng.integers(0, n, (n_boot, n))
    b = np.sqrt(per[idx].mean(1))
    return float(np.sqrt(per.mean())), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def boot_diff(E2a, E2b, n_boot=2000, seed=42):
    """RMS(a) - RMS(b), paired over ligands."""
    rng = np.random.default_rng(seed)
    n = E2a.shape[0]; pa = np.nanmean(E2a, axis=1); pb = np.nanmean(E2b, axis=1)
    idx = rng.integers(0, n, (n_boot, n))
    b = np.sqrt(pa[idx].mean(1)) - np.sqrt(pb[idx].mean(1))
    return float(np.sqrt(pa.mean()) - np.sqrt(pb.mean())), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def boot_spearman(x, y, n_boot=2000, seed=42):
    from scipy.stats import spearmanr
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y)); x, y = x[ok], y[ok]
    rng = np.random.default_rng(seed); n = len(x); out = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if len(set(x[i])) < 2 or len(set(y[i])) < 2:
            continue
        out.append(spearmanr(x[i], y[i]).statistic)
    return float(spearmanr(x, y).statistic), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)), int(n)
