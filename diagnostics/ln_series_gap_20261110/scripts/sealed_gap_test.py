"""Ready-to-run test of FROZEN_PREDICTION.md on a table of protein REE profiles.
NOT run on the sealed ortholog tables by the lead agent; Jacob runs it if and when he chooses.

Usage:
    python3 sealed_gap_test.py PROFILES.tsv [--scale log10|ln|kcal] [--sign assoc|dissoc]

PROFILES.tsv: one row per protein; first column an identifier; further columns named by element symbol
(La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu, optionally Y). Missing values may be blank or NA.
--scale: how the values relate to a binding constant. log10 (default): log10 of an association-like quantity
         (larger = tighter); ln: natural log; kcal: already a free energy in kcal/mol.
--sign:  assoc (default): larger values mean tighter binding; dissoc: larger values mean weaker binding
         (for example log Kd). With --scale kcal, negative = tighter is assumed and --sign is ignored.
Proteins with fewer than 10 lanthanides are skipped.
"""
import argparse
import sys
import numpy as np
import pandas as pd

LN = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
R9 = np.array([1.216, 1.196, 1.179, 1.163, 1.132, 1.120, 1.107, 1.095, 1.083, 1.072, 1.062, 1.052, 1.042, 1.032])
R9_Y = 1.075
KC = 1.3643


def basis(x, deg=3):
    z = (np.asarray(x, float) - 1.124) / 0.06
    return np.vstack([z ** k for k in range(deg + 1)]).T


def rss_curve(Y, k, grid, deg=3):
    """ligands x grid matrix of residual sums of squares with a break after metal index k."""
    out = np.zeros((len(Y), len(grid)))
    for g, d in enumerate(grid):
        B = basis(R9 + d * (np.arange(14) <= k), deg)
        for i, y in enumerate(Y):
            ok = ~np.isnan(y); b = np.linalg.lstsq(B[ok], y[ok], rcond=None)[0]
            out[i, g] = ((y[ok] - B[ok] @ b) ** 2).sum()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table"); ap.add_argument("--scale", default="log10", choices=["log10", "ln", "kcal"])
    ap.add_argument("--sign", default="assoc", choices=["assoc", "dissoc"])
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()
    d = pd.read_csv(a.table, sep="\t", na_values=["NA", "na", ""])
    X = d[[m for m in LN if m in d.columns]].reindex(columns=LN).values.astype(float)
    f = {"log10": KC, "ln": 0.5925, "kcal": 1.0}[a.scale]
    sgn = 1.0 if a.scale == "kcal" else (-1.0 if a.sign == "assoc" else 1.0)
    Y = sgn * f * X                                   # kcal/mol, negative = tighter
    yY = sgn * f * d["Y"].values.astype(float) if "Y" in d.columns else np.full(len(d), np.nan)
    keep = (~np.isnan(Y)).sum(1) >= 10
    Y, yY = Y[keep], yY[keep]; n = len(Y)
    print(f"proteins used: {n} of {len(d)}; with Y: {int((~np.isnan(yY)).sum())}")
    grid = np.arange(-0.080, 0.1201, 0.001)
    rng = np.random.default_rng(42)
    print("\n1. Best single break at each neighbour pair (cubic polynomial per protein)")
    curves = {}
    for k in range(1, 13):
        curves[k] = rss_curve(Y, k, grid); t = curves[k].sum(0); j = int(np.argmin(t))
        t0 = curves[k][:, int(np.argmin(np.abs(grid)))].sum()
        print(f"   {LN[k]}|{LN[k+1]}: gap {1000*grid[j]:+6.1f} mA; residual removed {100*(1-t[j]/t0):5.1f}%")
    c = curves[6]; j = int(np.argmin(c.sum(0)))
    bs = [grid[int(np.argmin(c[rng.integers(0, n, n)].sum(0)))] for _ in range(a.boot)]
    print(f"\n2. Gap between Gd and Tb, pooled: {1000*grid[j]:+.1f} mA "
          f"[{1000*np.percentile(bs, 2.5):+.1f}, {1000*np.percentile(bs, 97.5):+.1f}] (bootstrap over proteins)")
    own = grid[np.argmin(c, axis=1)]
    print(f"   per protein: median {1000*np.median(own):+.1f} mA; IQR {1000*np.percentile(own, 25):+.1f} to {1000*np.percentile(own, 75):+.1f}")
    print("   FROZEN PREDICTION P1: pooled gap below +15 mA (small chelators: +31, interval +18.5 to +61.5).")
    print("   A clearly negative gap points to a measurement referenced to a small competing chelator.")
    if (~np.isnan(yY)).sum() >= 5:
        B = basis(R9); zy = (R9_Y - 1.124) / 0.06; an = []; tw = {m: [] for m in LN}
        for y, v in zip(Y, yY):
            if np.isnan(v):
                continue
            ok = ~np.isnan(y); b = np.linalg.lstsq(B[ok], y[ok], rcond=None)[0]
            an.append(v - sum(b[k] * zy ** k for k in range(4)))
            for m_, m in enumerate(LN):
                if not np.isnan(y[m_]):
                    tw[m].append(v - y[m_])
        an = np.array(an); bb = [np.mean(an[rng.integers(0, len(an), len(an))]) for _ in range(a.boot)]
        print(f"\n3. Yttrium against a lanthanide of its own radius (1.075 A) on each protein's curve: mean {an.mean():+.3f} kcal/mol "
              f"[{np.percentile(bb, 2.5):+.3f}, {np.percentile(bb, 97.5):+.3f}] (positive = yttrium binds more weakly); RMS {np.sqrt(np.mean(an**2)):.3f}")
        print("   closest lanthanides by RMS difference: " + ", ".join(f"{m} {np.sqrt(np.mean(np.square(v))):.2f}" for m, v in
                                                                      sorted(tw.items(), key=lambda kv: np.sqrt(np.mean(np.square(kv[1]))) if kv[1] else 9)[:4]))
        print("   FROZEN PREDICTION P2: mean between -0.05 and +0.25 kcal/mol, and Dy or Ho closest (as for oxygen-only chelators).")


if __name__ == "__main__":
    main()
