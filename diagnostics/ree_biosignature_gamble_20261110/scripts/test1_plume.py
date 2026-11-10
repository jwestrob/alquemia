#!/usr/bin/env python3
"""Test 1 (PLAN.md + ADDENDUM_test1_rule.md): removal spectrum of the rare earths in the Deepwater Horizon methane
plume. Inputs are the tidy tables written by parse_shiller.py. Lead agent, Claude Fable, 2026-11-10."""
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fingerprint import LN
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/test1"; LOG = open(ROOT / "results/test1_plume.log", "w")
EL = ["Y"] + LN
rng = np.random.default_rng(20261110)


def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()


J = pd.read_csv(OUT / "dwh_joined.tsv", sep="\t"); C = pd.read_csv(OUT / "dwh_companion.tsv", sep="\t")
may = J[J.sheet.isin(["Early May 2010", "Late May 2010"])].copy()
may["La_Yb"] = may.La / may.Yb

# oxygen baseline from October 2011 (companion table; 800-1,600 m)
b = C[(C.sheet == "October 2011") & C.depth.between(800, 1600) & C.o2.notna()]
slope_o2, icpt_o2 = np.polyfit(b.depth, b.o2, 1)
res_sd = float(np.std(b.o2 - (icpt_o2 + slope_o2 * b.depth), ddof=2))
p(f"oxygen baseline, October 2011, 800-1,600 m: n={len(b)}, O2 = {icpt_o2:.1f} + {slope_o2:.4f} x depth (uM); residual SD {res_sd:.1f} uM")
may["o2_deficit"] = (icpt_o2 + slope_o2 * may.depth) - may.o2


def boot_ratio(P, B, n=5000):
    """f = 1 - mean(P)/mean(B) for every element, with a bootstrap over samples in both sets."""
    f = 1 - P[EL].mean() / B[EL].mean()
    ip = rng.integers(0, len(P), (n, len(P))); ib = rng.integers(0, len(B), (n, len(B)))
    Pa = P[EL].to_numpy(); Ba = B[EL].to_numpy()
    fb = 1 - Pa[ip].mean(1) / Ba[ib].mean(1)                      # n x elements
    return f, fb


def report_r1(deep, sel_col, label, k=5):
    d = deep.sort_values(sel_col)
    P, B = d.iloc[:k], d.iloc[k:]
    f, fb = boot_ratio(P, B)
    lo, hi = np.percentile(fb, [2.5, 97.5], axis=0)
    p(f"\n--- {label}: plume n={len(P)}, background n={len(B)}")
    p("plume samples:", "; ".join(f"{r.station} {r.depth:.0f} m (La {r.La:.1f}, La/Yb {r.La_Yb:.2f}, O2 deficit {r.o2_deficit:+.0f})" for r in P.itertuples()))
    p(f"background: La {B.La.mean():.2f} +- {B.La.std():.2f} pmol/kg; La/Yb {B.La_Yb.mean():.2f} +- {B.La_Yb.std():.2f}; background CV by element: " +
      " ".join(f"{e} {100 * B[e].std() / B[e].mean():.1f}%" for e in EL))
    kl = np.log(1 - fb) / np.log(1 - fb[:, [EL.index('La')]])
    kr = np.log(1 - f) / np.log(1 - f["La"]); klo, khi = np.nanpercentile(kl, [2.5, 97.5], axis=0)
    rows = []
    for i, e in enumerate(EL):
        rows.append(dict(reading=label, element=e, removal=f[e], lo=lo[i], hi=hi[i], k_over_kLa=kr[e], k_lo=klo[i], k_hi=khi[i]))
        p(f"  {e:2s} removal {100 * f[e]:+6.1f}% [{100 * lo[i]:+6.1f}, {100 * hi[i]:+6.1f}]   k/k_La {kr[e]:+.2f} [{klo[i]:+.2f}, {khi[i]:+.2f}]")
    mid = fb[:, [EL.index(e) for e in ("Sm", "Eu", "Gd")]].mean(1); m0 = float(f[["Sm", "Eu", "Gd"]].mean())
    hv = fb[:, [EL.index(e) for e in ("Er", "Tm", "Yb", "Lu")]].mean(1); h0 = float(f[["Er", "Tm", "Yb", "Lu"]].mean())
    p(f"  mean removal Sm-Eu-Gd {100 * m0:+.1f}% [{100 * np.percentile(mid, 2.5):+.1f}, {100 * np.percentile(mid, 97.5):+.1f}];"
      f" Er-Lu {100 * h0:+.1f}% [{100 * np.percentile(hv, 2.5):+.1f}, {100 * np.percentile(hv, 97.5):+.1f}];"
      f" Nd minus La {100 * (f['Nd'] - f['La']):+.1f} points [{100 * np.percentile(fb[:, EL.index('Nd')] - fb[:, EL.index('La')], 2.5):+.1f}, {100 * np.percentile(fb[:, EL.index('Nd')] - fb[:, EL.index('La')], 97.5):+.1f}]")
    return pd.DataFrame(rows), dict(mid=m0, mid_ci=[float(np.percentile(mid, 2.5)), float(np.percentile(mid, 97.5))], heavy=h0,
                                    f={e: float(f[e]) for e in EL}, plume=[f"{r.station}@{r.depth:.0f}" for r in P.itertuples()])


deep = may[may.depth >= 800]
p(f"May 2010 samples at 800 m or deeper: {len(deep)} (early May {int((deep.sheet == 'Early May 2010').sum())}, late May {int((deep.sheet == 'Late May 2010').sum())}); with oxygen: {int(deep.o2.notna().sum())}")
deep.sort_values("La_Yb")[["sheet", "station", "depth", "La", "Yb", "La_Yb", "o2", "o2_deficit", "ba", "pah"] + EL].to_csv(OUT / "may2010_deep_samples.tsv", sep="\t", index=False, float_format="%.4g")
tabs, S = [], {}
for lab, d, col in [("R1 lowest La/Yb, depth >= 800 m", deep, "La_Yb"), ("R1s lowest La, depth >= 800 m", deep, "La"),
                    ("R1s lowest La/Yb, depth >= 1000 m", may[may.depth >= 1000], "La_Yb")]:
    t, s = report_r1(d, col, lab); tabs.append(t); S[lab] = s

# R3: oxygen deficit as the selector
w = may[may.depth.between(900, 1400) & may.o2_deficit.notna()].copy()
p(f"\n--- R3: May 2010, 900-1,400 m, with oxygen: n={len(w)}; oxygen deficit range {w.o2_deficit.min():+.1f} to {w.o2_deficit.max():+.1f} uM (median {w.o2_deficit.median():+.1f})")
p(f"    correlation of La/Yb with oxygen deficit: r = {np.corrcoef(w.La_Yb, w.o2_deficit)[0, 1]:+.2f}; of La: r = {np.corrcoef(w.La, w.o2_deficit)[0, 1]:+.2f}")
X = np.c_[np.ones(len(w)), w.o2_deficit.to_numpy()]; Y = np.log(w[EL].to_numpy())
beta = np.linalg.lstsq(X, Y, rcond=None)[0][1]
ib = rng.integers(0, len(w), (5000, len(w)))
bb = np.array([np.linalg.lstsq(X[i], Y[i], rcond=None)[0][1] for i in ib])
rows = []
for i, e in enumerate(EL):
    r10 = 1 - np.exp(10 * beta[i]); lo, hi = 1 - np.exp(10 * np.percentile(bb[:, i], [97.5, 2.5]))
    rows.append(dict(reading="R3 per 10 uM oxygen deficit", element=e, removal=r10, lo=lo, hi=hi, k_over_kLa=beta[i] / beta[EL.index("La")],
                     k_lo=np.percentile(bb[:, i] / bb[:, EL.index("La")], 2.5), k_hi=np.percentile(bb[:, i] / bb[:, EL.index("La")], 97.5)))
    p(f"  {e:2s} removal per 10 uM deficit {100 * r10:+6.2f}% [{100 * lo:+6.2f}, {100 * hi:+6.2f}]   slope ratio to La {rows[-1]['k_over_kLa']:+.2f} [{rows[-1]['k_lo']:+.2f}, {rows[-1]['k_hi']:+.2f}]")
mid_b = bb[:, [EL.index(e) for e in ("Sm", "Eu", "Gd")]].mean(1); mid0 = beta[[EL.index(e) for e in ("Sm", "Eu", "Gd")]].mean()
p(f"  mean slope Sm-Eu-Gd: {mid0:+.5f} per uM [{np.percentile(mid_b, 2.5):+.5f}, {np.percentile(mid_b, 97.5):+.5f}]  (La slope {beta[EL.index('La')]:+.5f})")
tabs.append(pd.DataFrame(rows)); T = pd.concat(tabs); T.to_csv(OUT / "plume_removal.tsv", sep="\t", index=False, float_format="%.4f")

# decision terms of the addendum
r1 = S["R1 lowest La/Yb, depth >= 800 m"]; f = r1["f"]
nd_ge_la = f["Nd"] >= f["La"]
mid_r1 = r1["mid"] >= 0.05 and r1["mid_ci"][0] > 0
mid_r3 = (mid0 < 0) and (np.percentile(mid_b, 97.5) < 0)
mono = f["La"] > f["Pr"] > f["Nd"]
p("\n=== decision terms ===")
p(f"Nd removal >= La removal (R1): {nd_ge_la}  (La {100 * f['La']:.1f}%, Ce {100 * f['Ce']:.1f}%, Pr {100 * f['Pr']:.1f}%, Nd {100 * f['Nd']:.1f}%)")
p(f"Sm-Gd removal detectable: R1 {mid_r1} (mean {100 * r1['mid']:+.1f}%, interval [{100 * r1['mid_ci'][0]:+.1f}, {100 * r1['mid_ci'][1]:+.1f}]); R3 {mid_r3}")
p(f"removal falls La > Pr > Nd (R1): {mono}")
if nd_ge_la and mid_r1 and mid_r3:
    verdict = "lanmodulin-like"
elif mono and not (mid_r1 and mid_r3) and (mid_r1 == mid_r3):
    verdict = "uptake-like"
else:
    verdict = "unresolved (see terms)"
p("categorical reading of test 1:", verdict)
json.dump(dict(summary=S, verdict=verdict, o2_baseline=dict(icpt=icpt_o2, slope=slope_o2, n=int(len(b)), resid_sd=res_sd),
               r3=dict(n=int(len(w)), mid_slope=float(mid0), mid_ci=[float(np.percentile(mid_b, 2.5)), float(np.percentile(mid_b, 97.5))], la_slope=float(beta[EL.index('La')]))),
          open(OUT / "test1.json", "w"), indent=1)
