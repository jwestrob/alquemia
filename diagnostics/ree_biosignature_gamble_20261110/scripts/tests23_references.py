#!/usr/bin/env python3
"""Tests 2-4 of PLAN.md and the speciation-corrected comparison for test 1.
Inputs: tidy tables in data/reference_tables/ and results/test1/. Lead agent, Claude Fable, 2026-11-10."""
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fingerprint import LN, fingerprint, VARIANTS
ROOT = Path(__file__).resolve().parents[1]; B = Path("/Users/jacob/Documents/Berkeley/Lanthanides")
REF = ROOT / "data/reference_tables"; OUT = ROOT / "results/tests23"; OUT.mkdir(exist_ok=True, parents=True)
LOG = open(ROOT / "results/tests23.log", "w")
rng = np.random.default_rng(20261110)
AW = dict(Y=88.906, La=138.905, Ce=140.116, Pr=140.908, Nd=144.242, Sm=150.36, Eu=151.964, Gd=157.25, Tb=158.925, Dy=162.500,
          Ho=164.930, Er=167.259, Tm=168.934, Yb=173.045, Lu=174.967)
EL = ["Y"] + LN


def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()


def fp(vals, variant="no_Ce_La"):
    v = {e: vals.get(e) for e in LN}
    return fingerprint(v, vals.get("Y"), VARIANTS[variant])


def cls(f):
    """Categories stated in PLAN.md (test 2)."""
    if f is None or "y" not in f or not np.isfinite(f.get("u", np.nan)) or "Y_minus_Ho" not in f:
        return "undefined"
    prot = abs(f["u"]) < 1 and -0.3 <= f["Y_minus_Ho"] <= -0.1
    chel = f["u"] > 1 and abs(f["Y_minus_Ho"]) <= 0.05
    return "protein-type" if prot else "chelator-type" if chel else "neither"


# ------------------------------------------------------------------ 1. free-ion fractions in seawater
K = pd.read_csv(REF / "ree_complexation_constants.tsv", sep="\t", comment="#").set_index("element")
I = 0.7; dav = 0.509 * (np.sqrt(I) / (1 + np.sqrt(I)) - 0.3 * I)
g3, g2, g1 = 10 ** (-9 * dav), 10 ** (-4 * dav), 10 ** (-dav)
LOGA = dict(SO4=np.log10(0.0282 * 0.104), Cl=np.log10(0.546 * 0.666), HCO3=np.log10(1.9e-3 * 0.55))


def log_alpha(cset, log_a_co3, pH=7.9):
    """log10 of the free-ion fraction of every element. cset: 'LB1998' or 'OK2000' carbonate constants."""
    out = {}
    for e in K.index:
        r = K.loc[e]
        s = 1.0
        s += 10 ** (r[f"b1_CO3_{cset}"] + log_a_co3) * g3 / g1
        s += 10 ** (r[f"b2_CO3_{cset}"] + 2 * log_a_co3) * g3 / g1
        s += 10 ** (r["b1_HCO3_LB1998"] + LOGA["HCO3"]) * g3 / g2
        s += 10 ** (r["b1_OH_M1992"] + (pH - 14.0)) * g3 / g2
        s += 10 ** (r["b1_SO4_M1992"] + LOGA["SO4"]) * g3 / g1
        s += 10 ** (r["b1_Cl_Mir1982"] + LOGA["Cl"]) * g3 / g2
        out[e] = -np.log10(s)
    return out


p("=== 1. free-ion fractions in seawater (constants at 25 C, infinite dilution; Davies coefficients at I = 0.7) ===")
p(f"activity coefficients: 3+ {g3:.3f}, 2+ {g2:.3f}, 1+ {g1:.3f}; log activities: SO4 {LOGA['SO4']:.2f}, Cl {LOGA['Cl']:.2f}, HCO3 {LOGA['HCO3']:.2f}")
ALPHA = {}
for cset in ("LB1998", "OK2000"):
    for la in (-5.80, -5.45, -5.10):
        a = log_alpha(cset, la); ALPHA[(cset, la)] = a
        p(f"  {cset} log a(CO3) {la:+.2f}: free fraction La {100 * 10 ** a['La']:.2f}%  Nd {100 * 10 ** a['Nd']:.2f}%  Gd {100 * 10 ** a['Gd']:.3f}%  Lu {100 * 10 ** a['Lu']:.4f}% | "
          "relative to La: " + " ".join(f"{e} {10 ** (a[e] - a['La']):.2f}" for e in ["Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Yb", "Lu", "Y"]))
CENTRAL = ("LB1998", -5.45)
f_alpha = {k: fp(v) for k, v in ALPHA.items()}
for k, f in f_alpha.items():
    p(f"  fingerprint of log(free fraction), {k[0]} {k[1]:+.2f}: step h {f['h']:+.3f} (se {f['h_se']:.3f}), slope per step {f['slope_step']:+.3f}, u {f['u']:+.2f}, yttrium offset {f['y']:+.3f}, Y-Ho {f['Y_minus_Ho']:+.3f}")

# ------------------------------------------------------------------ 2. reference binders
lm = pd.read_csv(B / "lanm_ortholog_profiles_2026-11-10/results/qc_mean_profiles_log10.tsv", sep="\t")
P = {e: float(lm[e].mean()) for e in EL}
p("\n=== 2. lanmodulin mean capture profile (353 orthologs, log10, relative to La) ===")
p("  " + " ".join(f"{e} {P[e] - P['La']:+.3f}" for e in EL))
ox = pd.read_csv(REF / "kawabe1999_kd_e4.tsv", sep="\t", comment="#")
rows = []
for r in ox.itertuples():
    v = {e: np.log10(getattr(r, e)) for e in EL}
    for vn in VARIANTS:
        f = fp(v, vn)
        rows.append(dict(cls=("iron oxide, NaCl solution" if r.NaHCO3_M == 0 else "iron oxide, carbonate solution"), id=r.run, variant=vn, category=cls(f), **f))
for (cset, la), a in ALPHA.items():
    for name, base in (("flat binder in seawater (computed)", {e: 0.0 for e in EL}), ("lanmodulin in seawater (computed)", P)):
        v = {e: base[e] + a[e] for e in EL}
        for vn in VARIANTS:
            f = fp(v, vn); rows.append(dict(cls=name, id=f"{cset} {la:+.2f}", variant=vn, category=cls(f), **f))
    nacl = ox[ox.NaHCO3_M == 0]
    v = {e: float(np.log10(nacl[e]).mean()) + a[e] for e in EL}
    for vn in VARIANTS:
        f = fp(v, vn); rows.append(dict(cls="iron oxide in seawater (computed from NaCl runs)", id=f"{cset} {la:+.2f}", variant=vn, category=cls(f), **f))
f0 = fp(P); rows.append(dict(cls="lanmodulin mean, free ions", id="mean of 353", variant="no_Ce_La", category=cls(f0), **f0))

# ------------------------------------------------------------------ 3. mussels
T = pd.read_csv(REF / "wang2020_tissues.tsv", sep="\t", comment="#")
ELM = [e for e in EL if e != "Tm"]
bg = pd.read_csv(ROOT / "results/test1/may2010_deep_samples.tsv", sep="\t")          # sorted by La/Yb; first five = plume set
sw = bg.iloc[5:][EL].mean()                                                           # deep Gulf background, pmol/kg (stand-in for deep seawater)
sed = T[T.organ == "sediment"]
p("\n=== 3. mussels (Wang et al. 2020) ===")
p("deep-water stand-in (Gulf of Mexico background, pmol/kg): " + " ".join(f"{e} {sw[e]:.2f}" for e in EL))
sel_rows = []
for (sp, site, sym), g in T[T.organ.isin(["gill", "mantle"])].groupby(["species", "site", "symbionts"]):
    for an, a in g.groupby("animal"):
        if set(a.organ) != {"gill", "mantle"}:
            continue
        gi = a[a.organ == "gill"].iloc[0]; ma = a[a.organ == "mantle"].iloc[0]
        ratio = {e: (np.log10(gi[e] / ma[e]) if pd.notna(gi[e]) and pd.notna(ma[e]) else None) for e in ELM}
        tsw = {e: (np.log10((gi[e] / AW[e]) / sw[e]) if pd.notna(gi[e]) else None) for e in ELM}
        for name, v in ((f"gill over mantle, {sym}", ratio), (f"gill over deep water, {sym}", tsw)):
            for vn in VARIANTS:
                f = fp(v, vn)
                if f:
                    rows.append(dict(cls=name, id=f"{sp} {an}", variant=vn, category=cls(f), **f))
        b = float(np.median([gi[e] / ma[e] for e in ("Ho", "Er", "Yb", "Lu")]))
        exc = {e: gi[e] - b * ma[e] for e in ELM if pd.notna(gi[e])}
        E = {e: (exc[e] / AW[e]) / sw[e] for e in exc}
        sel_rows.append(dict(species=sp, site=site, symbionts=sym, animal=an, b=b, La_gill=gi["La"], La_over_mantle=gi["La"] / ma["La"],
                             Y_Ho_gill=gi["Y"] / gi["Ho"], Y_Ho_mantle=ma["Y"] / ma["Ho"],
                             **{f"k_{e}": (E[e] / E["La"] if E.get("La", 0) > 0 else np.nan) for e in ("Pr", "Nd", "Sm", "Eu", "Tb", "Dy")},
                             **{f"gm_{e}": gi[e] / ma[e] for e in ELM if pd.notna(gi[e])}))
S = pd.DataFrame(sel_rows); S.to_csv(OUT / "mussel_gill_selectivity.tsv", sep="\t", index=False, float_format="%.4g")
for sym, g in S.groupby("symbionts"):
    q = lambda c: f"{g[c].median():.3f} [{g[c].quantile(.25):.3f}, {g[c].quantile(.75):.3f}]"
    p(f"  {sym}: n={len(g)} animals; gill La {g.La_gill.median():.0f} ng/g; gill/mantle La {g.La_over_mantle.median():.1f}x; "
      f"gill excess per unit of dissolved element, relative to La: Pr {q('k_Pr')}, Nd {q('k_Nd')}, Sm {q('k_Sm')}, Eu {q('k_Eu')}, Tb {q('k_Tb')}, Dy {q('k_Dy')}")
    p(f"      gill over mantle by element (median): " + " ".join(f"{e} {g['gm_' + e].median():.2f}" for e in ELM if 'gm_' + e in g and g['gm_' + e].notna().any()))
    p(f"      Y/Ho mass ratio: gill {g.Y_Ho_gill.median():.1f}, mantle {g.Y_Ho_mantle.median():.1f}")
p(f"  sediments (n={len(sed)}): Y/Ho mass ratio {np.mean(sed.Y / sed.Ho):.1f}; deep-water stand-in Y/Ho mass ratio {(sw['Y'] * AW['Y']) / (sw['Ho'] * AW['Ho']):.1f}")

R = pd.DataFrame(rows); R.to_csv(OUT / "natural_and_reference_fingerprints.tsv", sep="\t", index=False, float_format="%.4f")
p("\n=== fingerprints (variant without Ce and La): median [min, max] ===")
for c, g in R[R.variant == "no_Ce_La"].groupby("cls", sort=False):
    q = lambda x: (f"{np.nanmedian(x):+.2f} [{np.nanmin(x):+.2f}, {np.nanmax(x):+.2f}]" if np.isfinite(x).any() else "n/a")
    p(f"  {c:50s} n={len(g):2d} | h {q(g.h.values)} | slope/step {q(g.slope_step.values)} | u {q(g.u.values)} | y {q(g.y.values)} | Y-Ho {q(g.Y_minus_Ho.values)} | {dict(g.category.value_counts())}")

# ------------------------------------------------------------------ 4. the plume against the predictions
d = bg.copy(); Pl, Bg = d.iloc[:5], d.iloc[5:]
n = 5000
ip = rng.integers(0, 5, (n, 5)); ib = rng.integers(0, len(Bg), (n, len(Bg)))
fb = 1 - Pl[EL].to_numpy()[ip].mean(1) / Bg[EL].to_numpy()[ib].mean(1)
f = 1 - Pl[EL].mean() / Bg[EL].mean()
heavy = [EL.index(e) for e in ("Er", "Tm", "Yb", "Lu")]
fh = fb[:, heavy].mean(1, keepdims=True); f_h = float(f[["Er", "Tm", "Yb", "Lu"]].mean())
fsel_b = 1 - (1 - fb) / (1 - fh); fsel = 1 - (1 - f) / (1 - f_h)                   # removal beyond the non-selective part seen in Er-Lu
iLa = EL.index("La")
kb = np.log(1 - fsel_b) / np.log(1 - fsel_b[:, [iLa]]); k = np.log(1 - fsel) / np.log(1 - fsel["La"])
p("\n=== 4. plume: rate constants relative to La (selective part: removal beyond the uniform loss seen in Er-Lu) ===")
p(f"  uniform loss in Er-Lu: {100 * f_h:.1f}%")
a = ALPHA[CENTRAL]
oxc = ox[ox.NaHCO3_M > 0]
pred = {"lanmodulin, free ions (plan)": {e: 10 ** (P[e] - P["La"]) for e in EL},
        "lanmodulin x free-ion fraction": {e: 10 ** (P[e] - P["La"] + a[e] - a["La"]) for e in EL},
        "flat binder x free-ion fraction": {e: 10 ** (a[e] - a["La"]) for e in EL},
        "iron oxide, carbonate solution (measured, median of 5 runs)": {e: float(np.median(oxc[e] / oxc["La"])) for e in EL}}
tab = []
for e in ["Pr", "Nd", "Sm", "Eu", "Gd", "Dy", "Yb", "Y"]:
    i = EL.index(e); lo, hi = np.nanpercentile(kb[:, i], [2.5, 97.5])
    row = dict(element=e, observed=float(k[e]), lo=float(lo), hi=float(hi), **{nm: v[e] for nm, v in pred.items()})
    lo_a = min(10 ** (ALPHA[kk][e] - ALPHA[kk]["La"]) for kk in ALPHA); hi_a = max(10 ** (ALPHA[kk][e] - ALPHA[kk]["La"]) for kk in ALPHA)
    row["alpha_ratio_range"] = f"{lo_a:.2f}-{hi_a:.2f}"; row["intrinsic"] = float(k[e]) / pred["flat binder x free-ion fraction"][e]
    row["intrinsic_lo"] = float(lo) / hi_a; row["intrinsic_hi"] = float(hi) / lo_a
    tab.append(row)
    p(f"  {e:2s} observed k/k_La {k[e]:+.2f} [{lo:+.2f}, {hi:+.2f}] | lanmodulin free ions {row['lanmodulin, free ions (plan)']:.2f} | lanmodulin x free fraction {row['lanmodulin x free-ion fraction']:.2f} | "
      f"flat x free fraction {row['flat binder x free-ion fraction']:.2f} (range {row['alpha_ratio_range']}) | iron oxide in carbonate {row['iron oxide, carbonate solution (measured, median of 5 runs)']:.2f} | "
      f"intrinsic selectivity of the sink (observed / free-fraction ratio) {row['intrinsic']:.2f} [{row['intrinsic_lo']:.2f}, {row['intrinsic_hi']:.2f}]")
pd.DataFrame(tab).to_csv(OUT / "plume_vs_predictions.tsv", sep="\t", index=False, float_format="%.4f")
sf = kb[:, EL.index("Pr")] / kb[:, EL.index("Nd")]; sf = sf[np.isfinite(sf) & (kb[:, EL.index("Nd")] > 0)]
p(f"  Pr over Nd rate ratio in the plume: {k['Pr'] / k['Nd']:.1f} [{np.percentile(sf, 2.5):.1f}, {np.percentile(sf, 97.5):.1f}] (bootstrap draws with positive Nd removal: {len(sf)} of {n}); "
  f"free-fraction ratio Pr/Nd {10 ** (a['Pr'] - a['Nd']):.2f}")
p(f"  lanmodulin mean: Nd - Pr {P['Nd'] - P['Pr']:+.3f} log10 (x{10 ** (P['Nd'] - P['Pr']):.2f}); Nd - La {P['Nd'] - P['La']:+.3f} (x{10 ** (P['Nd'] - P['La']):.2f}); Sm - La {P['Sm'] - P['La']:+.3f}; Gd - Sm {P['Gd'] - P['Sm']:+.3f}")
json.dump(dict(alpha_central={e: a[e] for e in EL}, lanmodulin_profile=P, plume_k={e: float(k[e]) for e in EL}, uniform_loss=f_h), open(OUT / "tests23.json", "w"), indent=1)

# ------------------------------------------------------------------ 5. supplementary numbers for REPORT.md
RF = pd.read_csv(ROOT / "results/reference_fingerprints.tsv", sep="\t"); RF = RF[RF.variant == "no_Ce_La"]
lmf = RF[RF.cls == "lanmodulin library"]; cho = RF[(RF.cls == "chelator") & (RF.donors == "O only")]; chn = RF[(RF.cls == "chelator") & (RF.donors == "with N")]
inbox = lambda g: int(((g.u.abs() < 1) & g.Y_minus_Ho.between(-0.3, -0.1)).sum())
p("\n=== 5. supplementary ===")
p(f"  inside the protein-type box: lanmodulins {inbox(lmf)}/{len(lmf)}; oxygen-donor chelators {inbox(cho)}/{len(cho)}; chelators with nitrogen {inbox(chn)}/{len(chn)}")
p(f"  Y - Ho (log10), median [quartiles]: lanmodulins {lmf.Y_minus_Ho.median():+.3f} [{lmf.Y_minus_Ho.quantile(.25):+.3f}, {lmf.Y_minus_Ho.quantile(.75):+.3f}]; "
  f"oxygen-donor chelators {cho.Y_minus_Ho.median():+.3f} [{cho.Y_minus_Ho.quantile(.25):+.3f}, {cho.Y_minus_Ho.quantile(.75):+.3f}] (n {int(cho.Y_minus_Ho.notna().sum())}); "
  f"chelators with nitrogen {chn.Y_minus_Ho.median():+.3f} (n {int(chn.Y_minus_Ho.notna().sum())})")
def auc(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float); a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    return float((np.sum(a[:, None] > b[None, :]) + 0.5 * np.sum(a[:, None] == b[None, :])) / (len(a) * len(b)))
on = R[(R.cls == "iron oxide, NaCl solution") & (R.variant == "no_Ce_La")]; oc = R[(R.cls == "iron oxide, carbonate solution") & (R.variant == "no_Ce_La")]
p(f"  iron oxide in salt solution against lanmodulins, single-pattern separation (AUC): step u {auc(lmf.u, on.u):.3f}; Y - Ho {auc(lmf.Y_minus_Ho, on.Y_minus_Ho):.3f}")
p(f"  iron oxide in carbonate solution against lanmodulins (AUC): step u {auc(oc.u, lmf.u):.3f}; Y - Ho {auc(lmf.Y_minus_Ho, oc.Y_minus_Ho):.3f}")
nacl = ox[ox.NaHCO3_M == 0]
for kk in (("LB1998", -5.80), ("LB1998", -5.45), ("LB1998", -5.10), ("OK2000", -5.45)):
    aa = ALPHA[kk]
    p(f"  iron oxide in seawater, computed from the 16 salt-solution runs ({kk[0]}, log a(CO3) {kk[1]:+.2f}): rate relative to La: " +
      " ".join(f"{e} {float(np.median(nacl[e] / nacl['La'])) * 10 ** (aa[e] - aa['La']):.2f}" for e in ["Pr", "Nd", "Sm", "Eu", "Gd", "Dy", "Yb", "Y"]))
p("  measured carbonate runs, rate relative to La by run: " + "; ".join(f"{r.run} (pH {r.pH}, NaHCO3 {1000 * r.NaHCO3_M:.1f} mM): Nd {r.Nd / r.La:.2f} Sm {r.Sm / r.La:.2f} Gd {r.Gd / r.La:.2f}" for r in oxc.itertuples()))
