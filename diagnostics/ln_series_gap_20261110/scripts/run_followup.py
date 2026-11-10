"""Follow-up analyses, written after the frozen readings were scored (post hoc; see DEVIATIONS.md).
1. External, parameter-free candidate for the shared offset: Kawabe's hydration stabilisation of the light aqua
   ions, dGh* (Kawabe 1999a Table 4; Geochem. J. 33:249), as transcribed by the literature subagent.
2. Models with a cubic per ligand, with and without the shared parts; rank-3 ceiling.
3. Kink spectrum and the gadolinium break against N atoms.
4. Yttrium: twin scan, anomaly against N atoms (N counted from the molecular formula), corrected slopes.
5. f-electron shapes with Kawabe's (q + 25) scaling; f-electron shapes after the aqua-ion correction.
Lead agent, Claude Fable."""
import json
import re
import sys
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist, spearmanr, rankdata

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa

RES = ROOT / "results"
LOG = open(RES / "followup.log", "w")
SUM = {}


def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()


def f3(t):
    return f"{t[0]:+.4f} [{t[1]:+.4f}, {t[2]:+.4f}]"


def cos_test(a, b, dim):
    c = float(a @ b / np.linalg.norm(a) / np.linalg.norm(b))
    pv = float(beta_dist.sf(c * c, 0.5, (dim - 1) / 2))
    return c, pv


def partial_spearman(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    Z = np.c_[np.ones(len(rz)), rz]
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]; ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(ex @ ey / np.linalg.norm(ex) / np.linalg.norm(ey))


def boot_stat(fn, n, n_boot=2000, seed=42):
    rng = np.random.default_rng(seed); out = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        try:
            v = fn(i)
        except Exception:
            continue
        if np.isfinite(v):
            out.append(v)
    return float(fn(np.arange(n))), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def main():
    Y, yY, meta = load_nist(14)
    n = Y.shape[0]; A = np.arange(14); I14 = np.eye(14)
    meta["N_atoms"] = [int(m.group(1)) if (m := re.search(r"N(\d+)", f)) else 0 for f in meta.formula]
    nN = meta.N_atoms.values.astype(float)
    dent = meta.nlig_La_cn8.values.astype(float)
    p(f"N atoms per ligand: {dict(pd.Series(nN).value_counts().sort_index())}")
    BQ = poly_basis(R9); BK = poly_basis(R9, 3)
    HQ = projector(BQ); HK = projector(BK)
    ybar = Y.mean(0)
    z = (R9 - X0) / XS

    # ---------- 1. Kawabe's aqua-ion hydration stabilisation as the shared offset ----------
    dGh_kJ = np.array([-7.16, -5.27, -4.09, -3.02, -0.67, -0.17, -0.06, 0, 0, 0, 0, 0, 0, 0])   # La..Lu, Pm dropped
    c_k = -dGh_kJ / 4.184            # kcal/mol added to y = -RT ln K (aqua ion more stable -> complex less favoured)
    p("\n1. Kawabe dGh* as an offset (kcal/mol, raw):", " ".join(f"{m}:{v:+.3f}" for m, v in zip(LN, c_k)))
    boots = np.load(RES / "bootstrap_shared.npz")
    for nm, H, dim in (("quadratic", HQ, 11), ("cubic", HK, 10)):
        P = I14 - H
        obs = P @ ybar; pred = P @ c_k
        c, pv = cos_test(obs, pred, dim)
        k = float(obs @ pred / (pred @ pred))
        rng = np.random.default_rng(42); kb = []; cb = []
        for _ in range(2000):
            i = rng.integers(0, n, n); o = P @ Y[i].mean(0)
            kb.append(o @ pred / (pred @ pred)); cb.append(o @ pred / np.linalg.norm(o) / np.linalg.norm(pred))
        p(f"  after the per-ligand {nm}: observed mean residual  ", " ".join(f"{m}:{v:+.3f}" for m, v in zip(LN, obs)))
        p(f"  after the per-ligand {nm}: Kawabe-predicted residual", " ".join(f"{m}:{v:+.3f}" for m, v in zip(LN, pred)))
        p(f"    cosine {c:+.3f} [{np.percentile(cb, 2.5):+.3f}, {np.percentile(cb, 97.5):+.3f}], exact p = {pv:.2g} "
          f"(random direction in {dim} dimensions); amplitude ratio observed/predicted {k:.3f} "
          f"[{np.percentile(kb, 2.5):.3f}, {np.percentile(kb, 97.5):.3f}]; share of observed pattern explained {c * c:.3f}")
        SUM[f"kawabe_{nm}"] = dict(cosine=c, cos_lo=float(np.percentile(cb, 2.5)), cos_hi=float(np.percentile(cb, 97.5)),
                                   p_exact=pv, amp_ratio=k, amp_lo=float(np.percentile(kb, 2.5)),
                                   amp_hi=float(np.percentile(kb, 97.5)), observed=obs.tolist(), predicted=pred.tolist())
    # generic comparators for the same statistic: how often does a smooth hinge look like this?
    P = I14 - HQ; obs = P @ ybar
    comp = {"step at Gd": (np.arange(14) >= 6).astype(float), "cubic": z ** 3, "linear hinge at Gd": np.maximum(R9 - R9[6], 0),
            "linear hinge at Sm": np.maximum(R9 - R9[4], 0), "linear hinge at Tb": np.maximum(R9 - R9[7], 0),
            "t1 (E1 shape)": T1, "t3 (E3 shape)": T3}
    for nm, v in comp.items():
        c, pv = cos_test(obs, P @ v, 11)
        p(f"    comparator {nm:20s}: cosine {c:+.3f} (share {c * c:.3f})")
        SUM.setdefault("comparators_quadratic", {})[nm] = c

    # ---------- 2. model scores with LOLO + LOMO, including cubic-per-ligand families and the Kawabe offset ----------
    names = ["Q", "O", "H", "Hk", "K", "KO", "KH", "KC", "KOC", "S3", "SO3", "T*", "QH+T", "KG"]
    E = {k: np.zeros((n, 14)) for k in names}
    Tz1 = T1 * (Q4F + 25) / 32.0; Tz3 = T3 * (Q4F + 25) / 32.0     # Kawabe's scaling, normalised at Gd
    kfit = []
    for i in range(n):
        tr = np.delete(Y, i, 0); yb = tr.mean(0); Z0 = np.zeros(14)
        # one shared amplitude for the Kawabe offset, fitted on the training ligands
        pr = (I14 - HQ) @ c_k; ob = (I14 - HQ) @ yb; kk = float(ob @ pr / (pr @ pr)); kfit.append(kk)
        rk, _ = fit_coordinate(tr, offset=False, deg=3, n_jitter=6)
        rko, _ = fit_coordinate(tr, offset=True, deg=3, n_jitter=6)
        Zc = tr - tr.mean(1, keepdims=True); _, _, Vt = np.linalg.svd(Zc, full_matrices=False)
        B3 = np.c_[np.ones(14), Vt[0], Vt[1], Vt[2]]
        Zo = tr - yb; Zo = Zo - Zo.mean(1, keepdims=True); _, _, Vo = np.linalg.svd(Zo, full_matrices=False)
        B3o = np.c_[np.ones(14), Vo[0], Vo[1], Vo[2]]
        # shared gadolinium-step amplitude is not fitted; KG = cubic + one free "Gd and heavier" step per ligand
        step = (np.arange(14) >= 6).astype(float)
        mods = {"Q": (BQ, Z0), "O": (BQ, yb), "H": (BQ, c_k), "Hk": (BQ, kk * c_k), "K": (BK, Z0), "KO": (BK, yb),
                "KH": (BK, c_k), "KC": (poly_basis(rk, 3), Z0), "KOC": (poly_basis(rko, 3), yb),
                "S3": (B3, Z0), "SO3": (B3o, yb), "T*": (np.c_[BQ, Tz1, Tz3], Z0), "QH+T": (np.c_[BQ, T1, T3], c_k),
                "KG": (np.c_[BK, step], Z0)}
        y = Y[i]
        for k, (B, c) in mods.items():
            for m_ in A:
                fi = A[A != m_]
                E[k][i, m_] = y[m_] - predict(B, c, y, fi)[0][m_]
    np.savez(RES / "followup_errors.npz", **{f"lomo_{k}": v for k, v in E.items()})
    p(f"\n2. S_LOMO, kcal/mol (post hoc models). Shared Kawabe amplitude fitted on training ligands: "
      f"mean {np.mean(kfit):.3f}, range {np.min(kfit):.3f}-{np.max(kfit):.3f}")
    SUM["scores"] = {}
    for k in names:
        a = boot_rms(E[k] ** 2); b = boot_rms(E[k][:, 1:-1] ** 2); SUM["scores"][k] = dict(S_LOMO=a, interior=b)
        p(f"  {k:5s} {a[0]:.4f} [{a[1]:.4f}, {a[2]:.4f}] | interior 12 metals {b[0]:.4f} [{b[1]:.4f}, {b[2]:.4f}]")
    SUM["diffs"] = {}
    for a, b in [("Q", "H"), ("Q", "Hk"), ("H", "O"), ("Hk", "O"), ("O", "K"), ("K", "KO"), ("K", "KH"), ("K", "KC"), ("KO", "KC"),
                 ("KO", "KOC"), ("K", "S3"), ("KO", "SO3"), ("K", "KG"), ("Q", "T*"), ("H", "QH+T"), ("Q", "K")]:
        d = boot_diff(E[a] ** 2, E[b] ** 2); SUM["diffs"][f"{a}-{b}"] = d
        p(f"  {a:4s} - {b:5s} {f3(d)}")

    # ---------- 3. kink spectrum (model free) ----------
    p("\n3. Kink at each interior metal: y(M) minus the straight line through its two neighbours in radius (kcal/mol)")
    kink = np.zeros((n, 12))
    for j in range(1, 13):
        w = (R9[j - 1] - R9[j]) / (R9[j - 1] - R9[j + 1])
        kink[:, j - 1] = Y[:, j] - ((1 - w) * Y[:, j - 1] + w * Y[:, j + 1])
    g0 = nN == 0; g1 = nN >= 1
    rows = []
    for j in range(12):
        m = LN[j + 1]; v = kink[:, j]
        rows.append(dict(metal=m, mean_all=v.mean(), se_all=v.std(ddof=1) / np.sqrt(n), mean_O_only=v[g0].mean(),
                         se_O_only=v[g0].std(ddof=1) / np.sqrt(g0.sum()), mean_N=v[g1].mean(),
                         se_N=v[g1].std(ddof=1) / np.sqrt(g1.sum()), spearman_N_atoms=spearmanr(nN, v).statistic))
    kt = pd.DataFrame(rows); kt.to_csv(RES / "kink_spectrum.tsv", sep="\t", index=False, float_format="%.4f")
    p(kt.round(3).to_string(index=False))
    gd = kink[:, 5]
    s = boot_spearman(nN, gd); SUM["gd_break_spearman_N_atoms"] = s
    p(f"  Gd break (kink at Gd): all {gd.mean():+.3f} (SE {gd.std(ddof=1) / np.sqrt(n):.3f}); O-only {gd[g0].mean():+.3f} "
      f"(SE {gd[g0].std(ddof=1) / np.sqrt(g0.sum()):.3f}, n={g0.sum()}); N-containing {gd[g1].mean():+.3f} "
      f"(SE {gd[g1].std(ddof=1) / np.sqrt(g1.sum()):.3f}, n={g1.sum()}); Spearman with N atoms {s[0]:+.3f} [{s[1]:+.3f}, {s[2]:+.3f}]")
    ps = boot_stat(lambda i: partial_spearman(nN[i], gd[i], dent[i]), n)
    p(f"  partial Spearman (Gd break, N atoms | denticity): {f3(ps)}")
    SUM["gd_break"] = dict(all=float(gd.mean()), O_only=float(gd[g0].mean()), N=float(gd[g1].mean()),
                           partial_given_denticity=ps)
    # Kawabe-corrected kink at Gd (aqua offset removed)
    Yc = Y - c_k
    kc = Yc[:, 6] - ((1 - (R9[5] - R9[6]) / (R9[5] - R9[7])) * Yc[:, 5] + ((R9[5] - R9[6]) / (R9[5] - R9[7])) * Yc[:, 7])
    p(f"  same after removing the Kawabe offset: all {kc.mean():+.3f}; O-only {kc[g0].mean():+.3f}; N-containing {kc[g1].mean():+.3f}")

    # ---------- 4. f-electron amplitudes: frozen definition, N atoms, partials, scaled shapes, aqua-corrected ----------
    p("\n4. f-electron amplitudes")
    BT = np.c_[BQ, T1, T3]
    bt = np.array([np.linalg.lstsq(BT, Y[i], rcond=None)[0] for i in range(n)])
    A1, A3 = bt[:, 3], bt[:, 4]
    btc = np.array([np.linalg.lstsq(BT, Yc[i], rcond=None)[0] for i in range(n)])
    A1c, A3c = btc[:, 3], btc[:, 4]
    BTz = np.c_[BQ, Tz1, Tz3]; btz = np.array([np.linalg.lstsq(BTz, Y[i], rcond=None)[0] for i in range(n)])
    BK1 = np.c_[BK, T1]; bk1 = np.array([np.linalg.lstsq(BK1, Y[i], rcond=None)[0] for i in range(n)])
    A1k = bk1[:, 4]
    slope_mid = (bt[:, 1]) / XS
    for nm, v in (("A1", A1), ("A3", A3), ("A1 after Kawabe offset", A1c), ("A3 after Kawabe offset", A3c),
                  ("A1 with (q+25) scaling", btz[:, 3]), ("A3 with (q+25) scaling", btz[:, 4]),
                  ("A1 on a cubic baseline", A1k)):
        s = boot_spearman(nN, v)
        ps = boot_stat(lambda i: partial_spearman(nN[i], v[i], dent[i]), n)
        p(f"  {nm:26s}: mean {v.mean():+.4f} (SE {v.std(ddof=1) / np.sqrt(n):.4f}); O-only {v[g0].mean():+.4f}, "
          f"N-containing {v[g1].mean():+.4f}; Spearman with N atoms {s[0]:+.3f} [{s[1]:+.3f}, {s[2]:+.3f}]; "
          f"partial given denticity {f3(ps)}")
        SUM.setdefault("f_amplitudes", {})[nm] = dict(mean=float(v.mean()), se=float(v.std(ddof=1) / np.sqrt(n)),
                                                      O_only=float(v[g0].mean()), N=float(v[g1].mean()), spearman=s,
                                                      partial=ps)
    # share of the mean residual in the f-electron plane after the aqua-ion correction
    for nm, H, dim in (("quadratic", HQ, 11), ("cubic", HK, 10)):
        P = I14 - H; e = P @ Yc.mean(0); W = P @ np.c_[T1, T3]; Qw, _ = np.linalg.qr(W)
        sh = float(((Qw.T @ e) ** 2).sum() / (e ** 2).sum()); pv = float(beta_dist.sf(sh, 1.0, (dim - 2) / 2))
        e0 = P @ ybar; sh0 = float(((Qw.T @ e0) ** 2).sum() / (e0 ** 2).sum()); pv0 = float(beta_dist.sf(sh0, 1.0, (dim - 2) / 2))
        p(f"  share of the mean residual in the f-electron plane, {nm} baseline: raw {sh0:.3f} (p = {pv0:.3g}); "
          f"after the Kawabe offset {sh:.3f} (p = {pv:.3g}); residual norm raw {np.linalg.norm(e0):.3f}, corrected {np.linalg.norm(e):.3f}")
        SUM[f"f_plane_share_{nm}"] = dict(raw=sh0, p_raw=pv0, corrected=sh, p_corrected=pv,
                                          norm_raw=float(np.linalg.norm(e0)), norm_corrected=float(np.linalg.norm(e)),
                                          resid_corrected=e.tolist())
        p("    mean residual after the Kawabe offset:", " ".join(f"{m}:{v:+.3f}" for m, v in zip(LN, e)))

    # ---------- 5. yttrium ----------
    iy = np.where(~np.isnan(yY))[0]; ny = len(iy)
    p(f"\n5. Yttrium ({ny} ligands)")
    tw = {m: float(np.sqrt(np.mean((yY[iy] - Y[iy, j]) ** 2))) for j, m in enumerate(LN)}
    p("  twin scan, RMS of y(Y) - y(M), all:      ", " ".join(f"{m}:{v:.2f}" for m, v in tw.items()))
    for nm, g in (("O-only", g0), ("N-containing", g1)):
        ii = [i for i in iy if g[i]]
        p(f"  twin scan, {nm:13s} (n={len(ii)}):", " ".join(f"{m}:{np.sqrt(np.mean((yY[ii] - Y[ii, j]) ** 2)):.2f}" for j, m in enumerate(LN)))
    SUM["yttrium_twin_rms"] = tw
    # apparent position from each ligand's own smooth curve in Shannon radius (no shared parameters at all):
    out = {}
    for nm, B, deg in (("quadratic", BQ, 2), ("cubic", BK, 3)):
        bet = np.array([np.linalg.lstsq(B, Y[i], rcond=None)[0] for i in range(n)])
        zy = (R9_Y - X0) / XS
        val = sum(bet[:, k] * zy ** k for k in range(deg + 1))
        sl = sum(k * bet[:, k] * zy ** (k - 1) for k in range(1, deg + 1)) / XS
        anom = yY - val                                   # kcal/mol; positive = Y complex weaker than its size twin
        pos = np.where(np.abs(sl) >= 5, R9_Y + anom / sl, np.nan)
        out[nm] = (anom, sl, pos)
        ok = ~np.isnan(pos)
        s = boot_spearman(nN, pos); sa = boot_spearman(nN, anom)
        ps = boot_stat(lambda i: partial_spearman(nN[ok][i], pos[ok][i], dent[ok][i]), int(ok.sum()))
        p(f"  {nm} baseline in Shannon radius: anomaly at r(Y) = 1.075: O-only {np.nanmean(anom[g0]):+.3f} "
          f"(SE {np.nanstd(anom[g0], ddof=1) / np.sqrt(np.sum(~np.isnan(anom[g0]))):.3f}), N-containing {np.nanmean(anom[g1]):+.3f} "
          f"(SE {np.nanstd(anom[g1], ddof=1) / np.sqrt(np.sum(~np.isnan(anom[g1]))):.3f}) kcal/mol; Spearman anomaly vs N atoms {sa[0]:+.3f} [{sa[1]:+.3f}, {sa[2]:+.3f}]")
        p(f"     apparent position (|slope| >= 5; n = {int(ok.sum())}): O-only median {np.nanmedian(pos[g0]):.4f} "
          f"(IQR {np.nanpercentile(pos[g0], 25):.4f}-{np.nanpercentile(pos[g0], 75):.4f}, n={int(np.sum(ok & g0))}); "
          f"N-containing median {np.nanmedian(pos[g1]):.4f} (IQR {np.nanpercentile(pos[g1], 25):.4f}-{np.nanpercentile(pos[g1], 75):.4f}, "
          f"n={int(np.sum(ok & g1))}); Spearman vs N atoms {s[0]:+.3f} [{s[1]:+.3f}, {s[2]:+.3f}]; partial given denticity {f3(ps)}")
        for v in sorted(set(nN)):
            sel = (nN == v) & ok
            if sel.sum():
                p(f"       N atoms = {int(v)}: median position {np.median(pos[sel]):.4f}, mean anomaly {np.nanmean(anom[nN == v]):+.3f} (n pos {int(sel.sum())})")
        SUM[f"yttrium_{nm}"] = dict(anom_O=float(np.nanmean(anom[g0])), anom_N=float(np.nanmean(anom[g1])),
                                    pos_O_median=float(np.nanmedian(pos[g0])), pos_N_median=float(np.nanmedian(pos[g1])),
                                    spearman_pos=s, spearman_anom=sa, partial=ps, n=int(ok.sum()))
        # bootstrap interval for the two medians and their difference
        io = np.where(ok & g0)[0]; inn = np.where(ok & g1)[0]; rng = np.random.default_rng(42); d = []; mo = []; mn = []
        for _ in range(2000):
            a = np.median(pos[rng.choice(io, len(io))]); b = np.median(pos[rng.choice(inn, len(inn))])
            mo.append(a); mn.append(b); d.append(b - a)
        p(f"     medians with 95% intervals: O-only {np.median(pos[io]):.4f} [{np.percentile(mo, 2.5):.4f}, {np.percentile(mo, 97.5):.4f}]; "
          f"N-containing {np.median(pos[inn]):.4f} [{np.percentile(mn, 2.5):.4f}, {np.percentile(mn, 97.5):.4f}]; "
          f"difference {np.median(pos[inn]) - np.median(pos[io]):+.4f} [{np.percentile(d, 2.5):+.4f}, {np.percentile(d, 97.5):+.4f}] A")
        SUM[f"yttrium_{nm}"]["medians"] = dict(O=[float(np.median(pos[io])), float(np.percentile(mo, 2.5)), float(np.percentile(mo, 97.5))],
                                               N=[float(np.median(pos[inn])), float(np.percentile(mn, 2.5)), float(np.percentile(mn, 97.5))],
                                               diff=[float(np.median(pos[inn]) - np.median(pos[io])), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))])
    tab = meta[["lig_id", "name", "formula", "N_atoms", "carboxyl", "alcohol", "nlig_La_cn8", "q_lig"]].copy()
    tab["gd_kink"] = gd; tab["A1"] = A1; tab["A3"] = A3
    tab["Y_anomaly_quadratic"] = out["quadratic"][0]; tab["Y_position_quadratic"] = out["quadratic"][2]
    tab["Y_anomaly_cubic"] = out["cubic"][0]; tab["Y_position_cubic"] = out["cubic"][2]; tab["slope_at_Y"] = out["cubic"][1]
    tab.to_csv(RES / "per_ligand_followup.tsv", sep="\t", index=False, float_format="%.4f")
    # frozen reading 6 and 7b with N atoms in place of the plan's (double-counting) feature
    p("\n  frozen readings recomputed with N atoms:")
    for nm, v in (("A1", A1), ("A3", A3)):
        s = boot_spearman(nN, v); p(f"    reading 6, {nm} vs N atoms: {s[0]:+.3f} [{s[1]:+.3f}, {s[2]:+.3f}]")
        SUM[f"reading6_{nm}_N_atoms"] = s
    pl = pd.read_csv(RES / "per_ligand.tsv", sep="\t")
    for k in ("Q", "C"):
        s = boot_spearman(nN, pl[f"rhoY_{k}"].values); SUM[f"reading7b_{k}_N_atoms"] = s
        p(f"    reading 7b, model {k} apparent position vs N atoms: {s[0]:+.3f} [{s[1]:+.3f}, {s[2]:+.3f}] (n = {s[3]})")
    json.dump(SUM, open(RES / "followup_summary.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
