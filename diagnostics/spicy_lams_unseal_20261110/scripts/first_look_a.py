"""Step 1, readings A1-A3, A6, A7 of PLAN.md (no sequence-similarity inputs needed)."""
import json, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
LOG = open(ROOT / "results/first_look_a.log", "w"); SUM = {}
def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
d, raw, C, mean, valid, qc = load()
n = int(qc.sum()); idx = np.where(qc)[0]
fam = d.family_original.fillna("NA").astype(str).values; gen = d.genus_original.fillna("NA").astype(str).values
p(f"QC orthologs: {n} of {len(d)}; families with >= 3 QC orthologs: {int((pd.Series(fam[qc]).value_counts() >= 3).sum())}; "
  f"largest: {dict(pd.Series(fam[qc]).value_counts().head(8))}")
G = mean[qc].mean(0)
co, yoff, rss = components(G[None], 3); Q, at = poly_design(3)
p("\nGlobal mean centred profile (log10 units):", " ".join(f"{e}:{v:+.3f}" for e, v in zip(EL, G)))
p("  its residual from a cubic polynomial over the lanthanides:", " ".join(f"{e}:{v:+.3f}" for e, v in zip(LN, G[1:] - (co @ Q.T)[0])),
  f"| yttrium offset {yoff[0]:+.3f}")
SUM["global_mean_profile"] = dict(zip(EL, G.round(4).tolist()))

# ---------- A1 ----------
three = qc & (valid.sum(1) == 3); i3 = np.where(three)[0]
cr, yr, _ = components(C[i3], 3)                      # proteins x 3 x 4 ; proteins x 3
X = dict(linear=cr[:, :, 1], quadratic=cr[:, :, 2], cubic=cr[:, :, 3], yttrium_offset=yr)
def icc(M):
    k = M.shape[1]; gm = M.mean(); msb = k * ((M.mean(1) - gm) ** 2).sum() / (len(M) - 1); msw = ((M - M.mean(1, keepdims=True)) ** 2).sum() / (len(M) * (k - 1))
    return (msb - msw) / (msb + (k - 1) * msw), msb, msw
rng = np.random.default_rng(42)
p(f"\nA1. Intraclass correlation over three replicates ({len(i3)} orthologs)")
SUM["A1"] = {}
for nm, M in X.items():
    v, msb, msw = icc(M); bs = [icc(M[rng.integers(0, len(M), len(M))])[0] for _ in range(2000)]
    rel = (msb - msw) / msb                           # reliability of the three-replicate mean
    sd_true = np.sqrt(max(0, (msb - msw) / 3)); sd_rep = np.sqrt(msw)
    p(f"  {nm:15s} ICC {v:.3f} [{np.percentile(bs, 2.5):.3f}, {np.percentile(bs, 97.5):.3f}]; reliability of the mean {rel:.3f}; "
      f"between-ortholog SD {sd_true:.3f}, replicate SD {sd_rep:.3f} (component units, log10)")
    SUM["A1"][nm] = dict(icc=float(v), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)), reliability_mean=float(rel),
                         sd_between=float(sd_true), sd_replicate=float(sd_rep))
ok = all(SUM["A1"][k]["icc"] >= 0.5 and SUM["A1"][k]["lo"] >= 0.4 for k in ("linear", "quadratic"))
p(f"  reading A1 (linear and quadratic ICC >= 0.5, lower bound >= 0.4): {'PASS' if ok else 'FAIL'}")
SUM["A1"]["pass"] = bool(ok)
# interpretable sizes: fitted Lu-La difference from the linear term, and centre-minus-ends from the quadratic term
lin_span = Q[13, 1] - Q[0, 1]; quad_span = Q[6, 2] - 0.5 * (Q[0, 2] + Q[13, 2])
cm, ym, rs3 = components(mean[idx], 3)
p(f"  in plain terms, across orthologs (means of replicates): Lu-minus-La from the linear term: mean {np.mean(cm[:,1])*lin_span:+.2f}, SD {np.std(cm[:,1])*abs(lin_span):.2f} log10; "
  f"Gd-minus-ends from the quadratic term: mean {np.mean(cm[:,2])*quad_span:+.2f}, SD {np.std(cm[:,2])*abs(quad_span):.2f} log10; "
  f"yttrium offset: mean {ym.mean():+.3f}, SD {ym.std():.3f} log10")

# ---------- A2 ----------
p("\nA2. Smoothness: residual mean square of polynomial fits against replicate noise of the mean profile")
se2 = np.nanvar(C[idx], axis=1, ddof=1) / valid[idx].sum(1)[:, None]      # proteins x 15
noise = float(np.nanmean(se2[:, 1:]))
SUM["A2"] = dict(noise_ms=noise)
for label, prof in (("mean profiles", mean[idx]), ("profiles minus the global mean profile", mean[idx] - G)):
    for deg in (1, 2, 3, 4):
        _, _, rss_ = components(prof, deg); ms = float(np.mean(rss_ / (14 - deg - 1)))
        p(f"  {label:40s} degree {deg}: residual MS {ms:.5f}; replicate noise MS {noise:.5f}; ratio {ms / noise:.2f}; residual RMS {np.sqrt(ms):.3f} log10")
        SUM["A2"][f"{label}|deg{deg}"] = dict(ms=ms, ratio=ms / noise)
r_frozen = SUM["A2"]["mean profiles|deg3"]["ratio"]
p(f"  reading A2 as frozen (cubic on mean profiles, ratio <= 2): {'PASS' if r_frozen <= 2 else 'FAIL'} (ratio {r_frozen:.2f})")
tot = float(np.mean(((mean[idx] - G)[:, 1:] ** 2).sum(1) / 13))
p(f"  variance of ortholog deviations from the global mean profile (per element): {tot:.5f}; share above replicate noise: {max(0, 1 - noise / tot):.2f}")

# ---------- A3 ----------
p("\nA3. Natural differences between orthologs (log10 units; noise-corrected SD; 5th-95th percentile of ortholog means)")
SUM["A3"] = {}
for a, b in (("Lu", "La"), ("Dy", "Nd"), ("Nd", "Pr"), ("Dy", "Tb"), ("Tb", "Gd"), ("Y", "Dy"), ("Y", "Ho")):
    ia, ib = EL.index(a), EL.index(b)
    rep = C[i3][:, :, ia] - C[i3][:, :, ib]; m = rep.mean(1)
    vt = max(0, m.var(ddof=1) - rep.var(axis=1, ddof=1).mean() / 3)
    p(f"  {a}-{b}: mean {m.mean():+.3f}; noise-corrected SD {np.sqrt(vt):.3f}; raw SD {m.std(ddof=1):.3f}; 5th-95th {np.percentile(m, 5):+.3f} to {np.percentile(m, 95):+.3f}; "
      f"replicate SD {np.sqrt(rep.var(axis=1, ddof=1).mean()):.3f}")
    SUM["A3"][f"{a}-{b}"] = dict(mean=float(m.mean()), sd_true=float(np.sqrt(vt)), p5=float(np.percentile(m, 5)), p95=float(np.percentile(m, 95)),
                                 sd_rep=float(np.sqrt(rep.var(axis=1, ddof=1).mean())))

# ---------- A6 ----------
F = pd.DataFrame([loop_features(s) for s in d[SEQ]])
efh = pd.to_numeric(d.EFhands, errors="coerce").values
p(f"\nEF-loop finder against the workbook's EFhands column: agree {int((F.n_loops.values == efh).sum())} of {len(d)}; "
  f"loops found {dict(F.n_loops.value_counts().sort_index())}; workbook {dict(pd.Series(efh).value_counts().sort_index())}")
p(f"  acidic residues per loop (all 600): {dict(F.acid_per_loop.round(2).value_counts().sort_index())}")
feats = ["n_loops", "acid_per_loop", "n_D9", "n_P2", "net_charge_per_res", "length"]
comp = pd.DataFrame(dict(linear=cm[:, 1], quadratic=cm[:, 2], cubic=cm[:, 3], yttrium_offset=ym), index=idx)
Fq = F.iloc[idx].reset_index(drop=True); famq = fam[idx]
big = pd.Series(famq).value_counts(); big = set(big[big >= 3].index)
p("\nA6. Spearman correlation of profile components with sequence features: pooled | within family (family means removed; families with >= 3)")
SUM["A6"] = {}
def within(v, f):
    v = np.asarray(v, float); out = np.full(len(v), np.nan)
    for g in set(f):
        m = (f == g)
        if g in big: out[m] = v[m] - np.nanmean(v[m])
    return out
fams = np.array(sorted(big))
for c in comp.columns:
    for ft in feats:
        x = Fq[ft].values.astype(float); y = comp[c].values
        okk = ~np.isnan(x)
        r0 = spearmanr(x[okk], y[okk]).statistic
        xw, yw = within(x, famq), within(y, famq); ok2 = ~np.isnan(xw) & ~np.isnan(yw)
        rw = spearmanr(xw[ok2], yw[ok2]).statistic if ok2.sum() > 10 and np.nanstd(xw[ok2]) > 0 else np.nan
        bs = []
        for _ in range(500):
            fs = rng.choice(fams, len(fams)); sel = np.concatenate([np.where((famq == g) & ok2)[0] for g in fs])
            if len(sel) > 10 and np.std(xw[sel]) > 0: bs.append(spearmanr(xw[sel], yw[sel]).statistic)
        lo, hi = (np.percentile(bs, 2.5), np.percentile(bs, 97.5)) if bs else (np.nan, np.nan)
        flag = " *" if (not np.isnan(lo)) and (lo > 0 or hi < 0) else ""
        p(f"  {c:15s} vs {ft:19s}: pooled {r0:+.2f} | within family {rw:+.2f} [{lo:+.2f}, {hi:+.2f}]{flag}")
        SUM["A6"][f"{c}|{ft}"] = dict(pooled=float(r0), within=float(rw), lo=float(lo), hi=float(hi))

# ---------- A7 ----------
mt = pd.to_numeric(d["Total_Metal_Replicate_Average"], errors="coerce").values
okk = ~np.isnan(F.acid_per_loop.values) & ~np.isnan(mt)
p(f"\nA7. Total captured metal against acidic residues per loop (all rows with loops, n = {int(okk.sum())}): Spearman {spearmanr(F.acid_per_loop.values[okk], mt[okk]).statistic:+.3f}")
for lo_, hi_ in ((0, 3.5), (3.5, 4.0), (4.0, 4.5), (4.5, 4.9), (4.9, 5.1)):
    m = okk & (F.acid_per_loop.values >= lo_) & (F.acid_per_loop.values < hi_)
    if m.sum(): p(f"   acidic per loop {lo_:.1f}-{hi_:.1f}: n {int(m.sum()):3d}; median total metal {np.median(mt[m]):6.1f} ppb; share passing QC {np.mean(qc[m]):.2f}")
p(f"   by number of loops found: " + "; ".join(f"{k} loops: n {int((F.n_loops.values == k).sum())}, median {np.nanmedian(mt[F.n_loops.values == k]):.1f} ppb, QC {np.mean(qc[F.n_loops.values == k]):.2f}" for k in sorted(set(F.n_loops.values))))
out = pd.concat([d[["workbook_row", "accession", "family_original", "genus_original", "EFhands"]].reset_index(drop=True), F], axis=1)
out["qc"] = qc; out["mtot"] = mt
for c in comp.columns:
    out[c] = np.nan; out.loc[idx, c] = comp[c].values
out.to_csv(ROOT / "results/ortholog_components.tsv", sep="\t", index=False, float_format="%.5f")
json.dump(SUM, open(ROOT / "results/first_look_a.json", "w"), indent=1, default=float)
