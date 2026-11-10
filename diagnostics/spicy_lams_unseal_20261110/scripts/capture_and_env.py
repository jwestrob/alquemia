"""Post hoc (after the frozen first look): (1) a better EF-loop finder; (2) total capture against loop acidity
with clade controls; (3) Jacob's environmental homologs placed against the library: nearest library neighbour,
loop acidity, expected capture class and, where a neighbour is close enough, a predicted profile."""
import json, re, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from Bio import Align
from Bio.Align import substitution_matrices
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
LOG = open(ROOT / "results/capture_and_env.log", "w"); SUM = {}
def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
CAND = re.compile(r"(?=([DN].[DNS].[DNSTG]......[DE]))")
def loops2(seq, thr=7):
    c = []
    for m in CAND.finditer(seq):
        l = m.group(1)
        s = (2 if l[0] == "D" else 1) + (2 if l[2] == "D" else 1) + (2 if l[4] == "D" else 1) + (l[5] == "G") + (l[7] in "ILVM") + (2 if l[11] == "E" else 1) + (l[8] in "DSTNE")
        if s >= thr: c.append((s, m.start(), l))
    out = []
    for s, st, l in sorted(c, key=lambda t: (-t[0], t[1])):
        if all(abs(st - o[0]) >= 12 for o in out): out.append((st, l))
    return sorted(out)
def feats2(seq):
    lp = loops2(seq); ac = [sum(l[i] in "DE" for i in (0, 2, 4, 8, 11)) for _, l in lp]
    return dict(n_loops2=len(lp), acid2=float(np.mean(ac)) if lp else np.nan, acid2_min=min(ac) if lp else np.nan, acid2_list=",".join(map(str, ac)),
                n_D9_2=sum(l[8] == "D" for _, l in lp), n_P2_2=sum(l[1] == "P" for _, l in lp), n_K2_2=sum(l[1] == "K" for _, l in lp),
                spacing=",".join(str(b[0] - a[0]) for a, b in zip(lp, lp[1:])), loops2=";".join(l for _, l in lp))
d, raw, C, mean, valid, qc = load()
F = pd.DataFrame([feats2(s) for s in d[SEQ]]); efh = pd.to_numeric(d.EFhands, errors="coerce").values
p(f"(1) improved loop finder: agrees with the workbook's EFhands for {int((F.n_loops2.values == efh).sum())} of {len(d)} (frozen finder: 192); found {dict(F.n_loops2.value_counts().sort_index())}")
mt = pd.to_numeric(d["Total_Metal_Replicate_Average"], errors="coerce").values
fam = d.family_original.fillna("unspecified").astype(str).values; gen = d.genus_original.fillna("unspecified").astype(str).values
ok = ~np.isnan(F.acid2.values) & ~np.isnan(mt)
p(f"\n(2) total captured metal (ppb) against mean acidic residues per loop, improved finder, n = {int(ok.sum())}: Spearman {spearmanr(F.acid2.values[ok], mt[ok]).statistic:+.3f}")
SUM["capture_bins"] = {}
for lo, hi in ((0, 3.3), (3.3, 3.8), (3.8, 4.2), (4.2, 4.6), (4.6, 4.9), (4.9, 5.1)):
    m = ok & (F.acid2.values >= lo) & (F.acid2.values < hi)
    if m.sum():
        p(f"   {lo:.1f}-{hi:.1f} acidic per loop: n {int(m.sum()):3d}; median capture {np.median(mt[m]):6.1f} ppb (IQR {np.percentile(mt[m], 25):.1f}-{np.percentile(mt[m], 75):.1f}); "
          f"capture >= 20 ppb in {np.mean(mt[m] >= 20):.2f}; families: {dict(pd.Series(fam[m]).value_counts().head(4))}")
        SUM["capture_bins"][f"{lo}-{hi}"] = dict(n=int(m.sum()), median=float(np.median(mt[m])), frac20=float(np.mean(mt[m] >= 20)))
for k in (3, 4):
    m = efh == k; p(f"   workbook EFhands = {k}: n {int(m.sum())}; median capture {np.nanmedian(mt[m]):.1f} ppb; capture >= 20 ppb in {np.nanmean(mt[m] >= 20):.2f}; mean acidic per loop {np.nanmean(F.acid2.values[m]):.2f}")
# clade control: within-genus and within-family rank correlation; and families that contain both kinds
def within(v, g, minn=4):
    out = np.full(len(v), np.nan)
    for x in set(g):
        m = (g == x) & ok
        if m.sum() >= minn and x != "unspecified": out[m] = v[m] - np.mean(v[m])
    return out
lm = np.log10(mt + 1)
for nm, g in (("family", fam), ("genus", gen)):
    a, b = within(F.acid2.values, g), within(lm, g); m = ~np.isnan(a) & ~np.isnan(b) & (np.abs(a) > 0)
    rng = np.random.default_rng(42); gs = np.array(sorted(set(g[m]))); bs = []
    for _ in range(1000):
        pick = rng.choice(gs, len(gs)); sel = np.concatenate([np.where(m & (g == x))[0] for x in pick])
        if len(sel) > 10 and np.std(a[sel]) > 0: bs.append(spearmanr(a[sel], b[sel]).statistic)
    p(f"   within {nm} (group means removed): Spearman {spearmanr(a[m], b[m]).statistic:+.3f} [{np.percentile(bs, 2.5):+.3f}, {np.percentile(bs, 97.5):+.3f}], n = {int(m.sum())}, groups = {len(gs)}")
    SUM[f"capture_within_{nm}"] = dict(rho=float(spearmanr(a[m], b[m]).statistic), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)), n=int(m.sum()))
p("   families holding both low (< 4.0) and high (>= 4.5) acidity members:")
for x in sorted(set(fam)):
    m = (fam == x) & ok; lo_ = m & (F.acid2.values < 4.0); hi_ = m & (F.acid2.values >= 4.5)
    if lo_.sum() >= 2 and hi_.sum() >= 2:
        p(f"     {x:22s}: low n {int(lo_.sum()):2d}, median capture {np.median(mt[lo_]):5.1f} ppb | high n {int(hi_.sum()):3d}, median capture {np.median(mt[hi_]):5.1f} ppb")
# does the tilt still track acidity with the better finder? (QC orthologs)
comp = pd.read_csv(ROOT / "results/ortholog_components.tsv", sep="\t")
q = qc & ok
for c in ("linear", "quadratic", "cubic", "yttrium_offset"):
    line = f"   {c:15s} vs acidic per loop (improved): pooled {spearmanr(F.acid2.values[q], comp[c].values[q]).statistic:+.2f}; vs loops with Asp9: {spearmanr(F.n_D9_2.values[q], comp[c].values[q]).statistic:+.2f}"
    a, b = within(F.acid2.values, gen), within(comp[c].values, gen); m = q & ~np.isnan(a) & ~np.isnan(b)
    line += f"; within genus vs acidity {spearmanr(a[m], b[m]).statistic:+.2f} (n={int(m.sum())})"
    p(line)
out = pd.concat([d[["workbook_row", "accession", "family_original", "genus_original", "EFhands"]].reset_index(drop=True), F], axis=1); out["mtot"] = mt; out["qc"] = qc
out.to_csv(ROOT / "results/loops_improved.tsv", sep="\t", index=False, float_format="%.4f")

# ---------- (3) environmental homologs ----------
p("\n(3) Environmental homologs against the library")
names, seqs = [], []
for l in open("/Users/jacob/Documents/Berkeley/Lanthanides/lanm_la_cofold_2026-11-09/inputs/all.faa"):
    if l.startswith(">"): names.append(l[1:].strip().split()[0]); seqs.append("")
    else: seqs[-1] += l.strip().rstrip("*")
o = pd.read_csv(ROOT / "data/ORTHOLOGS.tsv", sep="\t"); lib = o.mature_sequence.values; lacc = o.accession.values
posn = {a: i for i, a in enumerate(d.accession)}
al = Align.PairwiseAligner(); al.mode = "local"; al.substitution_matrix = substitution_matrices.load("BLOSUM62"); al.open_gap_score = -10; al.extend_gap_score = -0.5
def best(seq):
    res = []
    for j, t in enumerate(lib):
        a = al.align(seq, t)[0]; s1, s2 = a[0], a[1]
        m = sum(1 for x, y in zip(s1, s2) if x == y and x != "-")
        cov = sum(1 for y in s2 if y != "-") / len(t)
        res.append((m / len(s1) * min(1.0, cov), m / len(s1), cov, j))      # identity over aligned columns, discounted by coverage of the library sequence
    return sorted(res, reverse=True)
Yq = mean; rows = []
libF = F.set_index(d.accession.values)
for nm, s in zip(names, seqs):
    f = feats2(s); r = best(s)
    top = r[0]; acc = lacc[top[3]]
    # neighbours with profiles: QC, non-reserved
    nb = [(x[0], lacc[x[3]]) for x in r if lacc[x[3]] in posn and qc[posn[lacc[x[3]]]]][:5]
    rec = dict(name=nm, length=len(s), n_loops=f["n_loops2"], acid_per_loop=f["acid2"], acid_list=f["acid2_list"], K_at_2=f["n_K2_2"], P_at_2=f["n_P2_2"], Asp9=f["n_D9_2"],
               spacing=f["spacing"], nearest_library=acc, nearest_identity=round(top[1], 3), nearest_coverage=round(top[2], 3), nearest_score=round(top[0], 3),
               nearest_family=o.family.values[top[3]], nearest_capture_ppb=(mt[posn[acc]] if acc in posn else np.nan),
               nearest_qc_identity=round(nb[0][0], 3) if nb else np.nan)
    if nb and nb[0][0] >= 0.70:
        w = np.array([x[0] for x in nb]); P = np.array([Yq[posn[x[1]]] for x in nb]); pr = (P * (w / w.sum())[:, None]).sum(0)
        co, yo, _ = components(pr[None], 3)
        rec.update(pred_Lu_minus_La=round(float(pr[EL.index("Lu")] - pr[EL.index("La")]), 3), pred_Dy_minus_Nd=round(float(pr[EL.index("Dy")] - pr[EL.index("Nd")]), 3),
                   pred_Y_minus_Dy=round(float(pr[0] - pr[EL.index("Dy")]), 3))
    rows.append(rec)
E = pd.DataFrame(rows); E.to_csv(ROOT / "results/environmental_homologs.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 100)
show = E[["name", "length", "n_loops", "acid_per_loop", "acid_list", "K_at_2", "Asp9", "nearest_library", "nearest_identity", "nearest_coverage", "nearest_family", "nearest_capture_ppb", "nearest_qc_identity", "pred_Lu_minus_La", "pred_Dy_minus_Nd"]] if "pred_Lu_minus_La" in E else E
p(show[~show.name.str.startswith("omni")].to_string(index=False))
om = E[E.name.str.startswith("omni")]
p(f"\n   Omnitrophota compact EF-hand proteins (n = {len(om)}): loops {dict(om.n_loops.value_counts().sort_index())}; acidic per loop median {om.acid_per_loop.median():.2f} (range {om.acid_per_loop.min():.2f}-{om.acid_per_loop.max():.2f}); "
  f"nearest library identity median {om.nearest_identity.median():.2f} (max {om.nearest_identity.max():.2f}) at coverage median {om.nearest_coverage.median():.2f}; "
  f"nearest library families {dict(om.nearest_family.value_counts().head(4))}; capture of those nearest library members: median {om.nearest_capture_ppb.median():.1f} ppb")
json.dump(SUM, open(ROOT / "results/capture_and_env.json", "w"), indent=1, default=float)
