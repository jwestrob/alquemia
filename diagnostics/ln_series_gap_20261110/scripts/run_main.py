"""Parts A (offset vs coordinate), C (f-electron shapes) and D (yttrium) of PLAN.md. Lead agent, Claude Fable."""
import json
import sys
import time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa

RES = ROOT / "results"; RES.mkdir(exist_ok=True)
LOG = open(RES / "run.log", "w")
SUM = {}


def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()


def fmt(t):
    return f"{t[0]:.4f} [{t[1]:.4f}, {t[2]:.4f}]"


def main():
    t0 = time.time()
    Y, yY, meta = load_nist(14)
    n = Y.shape[0]; A = np.arange(14)
    p(f"ligands with complete series: {n}; with Y: {int((~np.isnan(yY)).sum())}")
    BQ = poly_basis(R9)

    # ---------------- positive control: 9 November numbers ----------------
    D = Y - Y[:, :1]; dr = R9 - R9[0]
    X = np.c_[dr[1:], dr[1:] ** 2]
    E = np.array([D[i, 1:] - X @ np.linalg.lstsq(X, D[i, 1:], rcond=None)[0] for i in range(n)])
    s1 = float(np.sqrt((E ** 2).mean()))
    p(f"positive control, F2 anchored at La: S1 = {s1:.3f} (9 Nov report: 0.164)")
    p("  mean residual by metal:", " ".join(f"{m}:{v:+.2f}" for m, v in zip(LN[1:], E.mean(0))))
    SUM["positive_control"] = dict(S1_F2=s1, mean_resid=dict(zip(LN[1:], np.round(E.mean(0), 3).tolist())))

    # ---------------- full-data shared fits ----------------
    rhoC, _, allC = fit_coordinate(Y, offset=False, return_all=True)
    rhoOC, cbar, allOC = fit_coordinate(Y, offset=True, return_all=True)
    for nm, rho, al in (("C", rhoC, allC), ("OC", rhoOC, allOC)):
        f = np.array([a[0] for a in al]); d = np.array([np.abs(a[2] - al[0][2]).max() for a in al])
        p(f"model {nm}: best objective {f[0]:.5f}; starts within 1e-6 of best: {(f - f[0] < 1e-6).sum()}/{len(f)};"
          f" distinct optima (max |d rho| > 1e-4 A): {(d > 1e-4).sum()}; objective range {f.min():.5f}..{f.max():.5f};"
          f" monotonic: {bool(np.all(np.diff(rho) < 0))}")
        p("   rho - r9 (mA):", " ".join(f"{m}:{1000 * v:+.1f}" for m, v in zip(LN, rho - R9)))
        SUM[f"rho_{nm}"] = dict(zip(LN, np.round(rho, 5).tolist()))
        SUM[f"rho_{nm}_fit"] = dict(best=float(f[0]), n_starts=len(f), n_at_best=int((f - f[0] < 1e-6).sum()),
                                    n_distinct=int((d > 1e-4).sum()), monotonic=bool(np.all(np.diff(rho) < 0)),
                                    worst=float(f.max()))
    HQ = projector(BQ)
    c_off = (np.eye(14) - HQ) @ Y.mean(0)
    p("model O: offset c (kcal/mol):", " ".join(f"{m}:{v:+.3f}" for m, v in zip(LN, c_off)))
    SUM["offset_O"] = dict(zip(LN, np.round(c_off, 4).tolist()))
    c_oc = (np.eye(14) - projector(poly_basis(rhoOC))) @ Y.mean(0)
    p("model OC: offset c (kcal/mol):", " ".join(f"{m}:{v:+.3f}" for m, v in zip(LN, c_oc)))
    SUM["offset_OC"] = dict(zip(LN, np.round(c_oc, 4).tolist()))

    # in-sample residual SS of every 3-number model on the full data (no cross-validation)
    def rss(B, c):
        Hm = projector(B); Z = Y - c
        return float(((Z - Z @ Hm) ** 2).sum())
    BS, _ = svd_shapes(Y, False); BSO, cSO = svd_shapes(Y, True)
    ins = dict(Q=rss(BQ, 0), O=rss(BQ, Y.mean(0)), C=rss(poly_basis(rhoC), 0), OC=rss(poly_basis(rhoOC), Y.mean(0)),
               S=rss(BS, 0), SO=rss(BSO, cSO))
    p("full-data in-sample RMS (kcal/mol):", " ".join(f"{k}:{np.sqrt(v / (n * 14)):.4f}" for k, v in ins.items()))
    SUM["insample_rms_full"] = {k: float(np.sqrt(v / (n * 14))) for k, v in ins.items()}

    # ---------------- leave-one-ligand-out ----------------
    names3 = ["Q", "O", "C", "OC", "S", "SO"]
    namesX = ["K", "Q4", "T", "O+T", "C+T", "OC+T"]
    E_lomo = {k: np.zeros((n, 14)) for k in names3 + namesX}
    E_4 = {k: np.full((n, 14), np.nan) for k in names3}
    E_in = {k: np.zeros((n, 14)) for k in names3}
    F4 = np.array([LN.index(m) for m in ("La", "Nd", "Dy", "Lu")]); R4 = np.array([j for j in A if j not in F4])
    shared = []
    for i in range(n):
        tr = np.delete(Y, i, 0); ybar = tr.mean(0)
        rc, _ = fit_coordinate(tr, offset=False)
        roc, _ = fit_coordinate(tr, offset=True)
        bs, _ = svd_shapes(tr, False); bso, cso = svd_shapes(tr, True)
        shared.append(dict(rc=rc, roc=roc, ybar=ybar))
        Z0 = np.zeros(14)
        mods = {"Q": (BQ, Z0), "O": (BQ, ybar), "C": (poly_basis(rc), Z0), "OC": (poly_basis(roc), ybar),
                "S": (bs, Z0), "SO": (bso, cso),
                "K": (poly_basis(R9, 3), Z0), "Q4": (poly_basis(R9, 4), Z0), "T": (np.c_[BQ, T1, T3], Z0),
                "O+T": (np.c_[BQ, T1, T3], ybar), "C+T": (np.c_[poly_basis(rc), T1, T3], Z0),
                "OC+T": (np.c_[poly_basis(roc), T1, T3], ybar)}
        y = Y[i]
        for k, (B, c) in mods.items():
            for m_ in A:
                fi = A[A != m_]
                E_lomo[k][i, m_] = y[m_] - predict(B, c, y, fi)[0][m_]
            if k in names3:
                E_4[k][i, R4] = (y - predict(B, c, y, F4)[0])[R4]
                E_in[k][i] = y - predict(B, c, y, A)[0]
        if i % 10 == 0:
            p(f"  LOLO {i + 1}/{n}  ({time.time() - t0:.0f} s)")
    np.savez(RES / "lolo_errors.npz", **{f"lomo_{k}": v for k, v in E_lomo.items()},
             **{f"s4_{k}": v for k, v in E_4.items()}, **{f"in_{k}": v for k, v in E_in.items()},
             rho_C_lolo=np.array([s["rc"] for s in shared]), rho_OC_lolo=np.array([s["roc"] for s in shared]))

    p("\nScores, kcal/mol (RMS [95% bootstrap over ligands])")
    SUM["scores"] = {}
    for k in names3 + namesX:
        a = boot_rms(E_lomo[k] ** 2)
        line = f"  {k:5s} S_LOMO {fmt(a)}"
        rec = dict(S_LOMO=a)
        inner = boot_rms(E_lomo[k][:, 1:-1] ** 2); rec["S_LOMO_interior"] = inner
        line += f" | interior 12 metals {fmt(inner)}"
        if k in names3:
            b = boot_rms(E_4[k][:, R4] ** 2); c_ = boot_rms(E_in[k] ** 2)
            rec.update(S_4=b, S_in=c_)
            line += f" | S_4 {fmt(b)} | S_in {fmt(c_)}"
        SUM["scores"][k] = rec; p(line)

    def diff(a, b, which="lomo"):
        Ea, Eb = (E_lomo[a], E_lomo[b]) if which == "lomo" else (E_4[a][:, R4], E_4[b][:, R4])
        return boot_diff(Ea ** 2, Eb ** 2)
    p("\nPaired differences in S_LOMO (first minus second; positive = second is better)")
    SUM["diffs"] = {}
    for a, b in [("Q", "O"), ("Q", "C"), ("Q", "OC"), ("O", "C"), ("C", "S"), ("OC", "SO"), ("O", "OC"), ("C", "OC"),
                 ("Q4", "T"), ("Q", "T"), ("K", "T"), ("Q", "K"), ("O", "O+T"), ("C", "C+T"), ("OC", "OC+T"), ("O", "S"),
                 ("O", "SO"), ("S", "SO")]:
        d = diff(a, b); SUM["diffs"][f"{a}-{b}"] = d
        line = f"  {a:4s} - {b:5s} {d[0]:+.4f} [{d[1]:+.4f}, {d[2]:+.4f}]"
        if a in names3 and b in names3:
            d4 = diff(a, b, "s4"); SUM["diffs"][f"{a}-{b}_S4"] = d4
            line += f"   | S_4: {d4[0]:+.4f} [{d4[1]:+.4f}, {d4[2]:+.4f}]"
        p(line)

    # LOLO stability of the coordinate
    RC = np.array([s["rc"] for s in shared]); ROC = np.array([s["roc"] for s in shared])
    p("\nLOLO spread of rho_C - r9 (mA): max |LOLO - full| per metal:",
      " ".join(f"{m}:{1000 * v:.1f}" for m, v in zip(LN, np.abs(RC - rhoC).max(0))))
    p("LOLO spread of rho_OC - r9 (mA): max |LOLO - full| per metal:",
      " ".join(f"{m}:{1000 * v:.1f}" for m, v in zip(LN, np.abs(ROC - rhoOC).max(0))))

    # bootstrap interval of the coordinate (resample ligands; starts: Shannon CN9 and the full-data optimum)
    rng = np.random.default_rng(42); nb = 1000
    BC = np.zeros((nb, 14)); BOC = np.zeros((nb, 14)); BOFF = np.zeros((nb, 14)); BOFFQ = np.zeros((nb, 14))
    for b in range(nb):
        idx = rng.integers(0, n, n); Yb = Y[idx]
        BC[b] = fit_coordinate(Yb, False, starts={"cn9": np.zeros(12), "full": (rhoC[1:-1] - R9[1:-1]) / 0.01})[0]
        BOC[b] = fit_coordinate(Yb, True, starts={"cn9": np.zeros(12), "full": (rhoOC[1:-1] - R9[1:-1]) / 0.01})[0]
        BOFF[b] = (np.eye(14) - projector(poly_basis(BOC[b]))) @ Yb.mean(0)
        BOFFQ[b] = (np.eye(14) - HQ) @ Yb.mean(0)
    np.savez(RES / "bootstrap_shared.npz", rho_C=BC, rho_OC=BOC, off_OC=BOFF, off_O=BOFFQ)
    tab = pd.DataFrame(dict(metal=LN, q4f=Q4F, r9=R9, r8=R8, rho_C=rhoC, rho_C_lo=np.percentile(BC, 2.5, 0),
                            rho_C_hi=np.percentile(BC, 97.5, 0), rho_OC=rhoOC, rho_OC_lo=np.percentile(BOC, 2.5, 0),
                            rho_OC_hi=np.percentile(BOC, 97.5, 0), off_O=c_off, off_O_lo=np.percentile(BOFFQ, 2.5, 0),
                            off_O_hi=np.percentile(BOFFQ, 97.5, 0), off_OC=c_oc, off_OC_lo=np.percentile(BOFF, 2.5, 0),
                            off_OC_hi=np.percentile(BOFF, 97.5, 0)))
    tab.to_csv(RES / "coordinate.tsv", sep="\t", index=False, float_format="%.5f")
    p("\ncoordinate table (A; offsets kcal/mol):\n" + tab.round(4).to_string(index=False))

    # ---------------- per-metal regression of the Q residual on local slope (descriptive) ----------------
    beta = np.array([np.linalg.lstsq(BQ, Y[i], rcond=None)[0] for i in range(n)])
    EQ = Y - beta @ BQ.T
    z = (R9 - X0) / XS
    slope = (beta[:, 1:2] + 2 * beta[:, 2:3] * z[None, :]) / XS        # kcal/mol/A, ligands x metals
    p("\nper metal: Q residual = c + d * local slope (c kcal/mol, d in mA); t values")
    rows = []
    for j, m in enumerate(LN):
        Xr = np.c_[np.ones(n), slope[:, j]]; co, *_ = np.linalg.lstsq(Xr, EQ[:, j], rcond=None)
        r_ = EQ[:, j] - Xr @ co; s2 = (r_ ** 2).sum() / (n - 2); cov = s2 * np.linalg.inv(Xr.T @ Xr)
        se = np.sqrt(np.diag(cov))
        rows.append(dict(metal=m, c=co[0], c_t=co[0] / se[0], d_mA=1000 * co[1], d_t=co[1] / se[1],
                         mean_resid=EQ[:, j].mean(), sd_resid=EQ[:, j].std(ddof=1),
                         mean_slope=slope[:, j].mean(), frac_slope_pos=float((slope[:, j] > 0).mean())))
    pm = pd.DataFrame(rows); pm.to_csv(RES / "per_metal_regression.tsv", sep="\t", index=False, float_format="%.4f")
    p(pm.round(3).to_string(index=False))

    # ---------------- f-electron shapes ----------------
    BT = np.c_[BQ, T1, T3]
    bt = np.array([np.linalg.lstsq(BT, Y[i], rcond=None)[0] for i in range(n)])
    A1, A3 = bt[:, 3], bt[:, 4]
    P = np.eye(14) - HQ
    e = P @ Y.mean(0); W = P @ np.c_[T1, T3]; Qw, _ = np.linalg.qr(W)
    share = float(((Qw.T @ e) ** 2).sum() / (e ** 2).sum()); pval = float((1 - share) ** 4.5)
    p(f"\nmean Q residual (free intercept): {' '.join(f'{m}:{v:+.3f}' for m, v in zip(LN, e))}")
    p(f"share of the mean Q residual in the plane of the f-electron shapes: {share:.3f}; exact p = {pval:.3g}")
    # same for generic cubic+quartic plane (descriptive comparator)
    Wg = P @ np.c_[z ** 3, z ** 4]; Qg, _ = np.linalg.qr(Wg)
    share_g = float(((Qg.T @ e) ** 2).sum() / (e ** 2).sum())
    Wc = P @ np.c_[z ** 3]; share_c = float(((Wc.T @ e) ** 2).sum() / ((Wc ** 2).sum() * (e ** 2).sum()))
    p(f"  comparator: share in the cubic+quartic plane {share_g:.3f} (p = {(1 - share_g) ** 4.5:.3g}); cubic alone {share_c:.3f}")
    # residual shares of the theory shapes themselves after the quadratic
    p(f"  norm of the theory shapes left after the quadratic: t1 {np.linalg.norm(P @ T1) / np.linalg.norm(T1 - T1.mean()):.3f},"
      f" t3 {np.linalg.norm(P @ T3) / np.linalg.norm(T3 - T3.mean()):.3f} (fraction of centred norm)")
    SUM["f_shapes"] = dict(share_mean_resid=share, p_exact=pval, share_cubic_quartic=share_g, share_cubic=share_c,
                           A1_mean=float(A1.mean()), A1_se=float(A1.std(ddof=1) / np.sqrt(n)),
                           A3_mean=float(A3.mean()), A3_se=float(A3.std(ddof=1) / np.sqrt(n)))
    p(f"  per-ligand amplitudes: A1 mean {A1.mean():+.4f} (SE {A1.std(ddof=1) / np.sqrt(n):.4f}), "
      f"A3 mean {A3.mean():+.4f} (SE {A3.std(ddof=1) / np.sqrt(n):.4f}) kcal/mol per unit shape")
    nN = meta.n_donor_N.values.astype(float)
    for nm, v in (("A1", A1), ("A3", A3)):
        s = boot_spearman(nN, v); SUM[f"spearman_{nm}_N"] = s
        p(f"  Spearman {nm} vs N donors: {s[0]:+.3f} [{s[1]:+.3f}, {s[2]:+.3f}] (n = {s[3]})")
    for feat in ("carboxyl", "q_lig", "nlig_La_cn8"):
        for nm, v in (("A1", A1), ("A3", A3)):
            s = boot_spearman(meta[feat].values, v); SUM[f"spearman_{nm}_{feat}"] = s
            p(f"  (descriptive) Spearman {nm} vs {feat}: {s[0]:+.3f} [{s[1]:+.3f}, {s[2]:+.3f}]")

    # ---------------- yttrium ----------------
    iy = np.where(~np.isnan(yY))[0]; ny = len(iy)
    grid = np.arange(0.90, 1.3501, 0.0005); zg = (grid - X0) / XS
    Gb = np.c_[np.ones(len(grid)), zg, zg ** 2]

    def y_fit(B, c, with_offset, fit_rows, test_row):
        """Shared position (and offset) from fit_rows; prediction error for test_row. Returns err, rhoY, cY."""
        bet = {j: np.linalg.lstsq(B, Y[j] - c, rcond=None)[0] for j in list(fit_rows) + [test_row]}
        Pm = np.array([Gb @ bet[j] for j in fit_rows])            # ligands x grid: predicted value at each position
        Rm = yY[fit_rows][:, None] - Pm
        cy = Rm.mean(0) if with_offset else np.zeros(len(grid))
        ss = ((Rm - cy) ** 2).sum(0); g = int(np.argmin(ss))
        return yY[test_row] - (Gb[g] @ bet[test_row] + cy[g]), grid[g], cy[g]
    EY = {k: np.zeros(ny) for k in ("Q", "O", "C", "OC", "twinHo", "twinDy", "twinEr", "twinTm", "twinBest", "Q@r9Y")}
    posY = {k: [] for k in ("Q", "O", "C", "OC")}
    for t, i in enumerate(iy):
        s = shared[i]; others = [j for j in iy if j != i]
        mods = {"Q": (BQ, np.zeros(14), False), "O": (BQ, s["ybar"], True),
                "C": (poly_basis(s["rc"]), np.zeros(14), False), "OC": (poly_basis(s["roc"]), s["ybar"], True)}
        for k, (B, c, wo) in mods.items():
            er, ry, cy = y_fit(B, c, wo, others, i); EY[k][t] = er; posY[k].append((ry, cy))
        for tw in ("Ho", "Dy", "Er", "Tm"):
            EY[f"twin{tw}"][t] = yY[i] - Y[i, LN.index(tw)]
        best = int(np.argmin([((yY[others] - Y[others, j]) ** 2).sum() for j in range(14)]))
        EY["twinBest"][t] = yY[i] - Y[i, best]
        b_ = np.linalg.lstsq(BQ, Y[i], rcond=None)[0]; zy = (R9_Y - X0) / XS
        EY["Q@r9Y"][t] = yY[i] - (b_[0] + b_[1] * zy + b_[2] * zy ** 2)
    p(f"\nYttrium (n = {ny} ligands), prediction error RMS kcal/mol")
    SUM["yttrium"] = {}
    for k, v in EY.items():
        a = boot_rms(v[:, None] ** 2); SUM["yttrium"][k] = a
        extra = ""
        if k in posY:
            pr = np.array(posY[k]); extra = f"  shared position {pr[:, 0].mean():.4f} A (LOLO range {pr[:, 0].min():.4f}-{pr[:, 0].max():.4f}), offset {pr[:, 1].mean():+.3f}"
        p(f"  {k:9s} {fmt(a)}{extra}")
    for k in ("Q", "O", "C", "OC"):
        for tw in ("twinHo", "twinDy", "twinBest"):
            d = boot_diff(EY[tw][:, None] ** 2, EY[k][:, None] ** 2); SUM["yttrium"][f"{tw}-{k}"] = d
            p(f"  {tw} - {k}: {d[0]:+.4f} [{d[1]:+.4f}, {d[2]:+.4f}]")
    # apparent position per ligand (full-data shared fits)
    perlig = meta[["lig_id", "name", "n_donor_N", "carboxyl", "q_lig", "nlig_La_cn8", "I"]].copy()
    perlig["A1"] = A1; perlig["A3"] = A3
    for k, (B, c, wo) in {"Q": (BQ, np.zeros(14), False), "O": (BQ, Y.mean(0), True),
                          "C": (poly_basis(rhoC), np.zeros(14), False), "OC": (poly_basis(rhoOC), Y.mean(0), True)}.items():
        bet = np.array([np.linalg.lstsq(B, Y[j] - c, rcond=None)[0] for j in range(n)])
        Pm = bet[iy] @ Gb.T; Rm = yY[iy][:, None] - Pm
        cy = Rm.mean(0) if wo else np.zeros(len(grid)); ss = ((Rm - cy) ** 2).sum(0); g = int(np.argmin(ss))
        ry = grid[g]; zy = zg[g]
        sl = (bet[:, 1] + 2 * bet[:, 2] * zy) / XS
        resid = np.full(n, np.nan); resid[iy] = Rm[:, g] - cy[g]
        app = np.where(np.abs(sl) >= 5, ry + resid / sl, np.nan)
        perlig[f"Yresid_{k}"] = resid; perlig[f"slopeY_{k}"] = sl; perlig[f"rhoY_{k}"] = app
        s = boot_spearman(nN, app); SUM[f"yttrium_spearman_{k}"] = dict(rho=s, shared_pos=float(ry), offset=float(cy[g]),
                                                                    n_used=s[3])
        p(f"  model {k}: shared position {ry:.4f} A, offset {cy[g]:+.3f}; apparent position vs N donors: "
          f"Spearman {s[0]:+.3f} [{s[1]:+.3f}, {s[2]:+.3f}] (n = {s[3]}); median apparent position by N donors: "
          + ", ".join(f"N={int(v)}: {np.nanmedian(app[nN == v]):.3f} (n={int(np.sum(~np.isnan(app[nN == v])))})"
                      for v in sorted(set(nN)) if np.sum(~np.isnan(app[nN == v])) > 0))
        # vertex and curvature
        if k in ("Q", "C", "OC"):
            perlig[f"vertex_{k}"] = X0 + XS * (-bet[:, 1] / (2 * bet[:, 2]))
            perlig[f"curv_{k}"] = bet[:, 2] / XS ** 2
    perlig.to_csv(RES / "per_ligand.tsv", sep="\t", index=False, float_format="%.5f")

    json.dump(SUM, open(RES / "summary.json", "w"), indent=1, default=float)
    p(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
