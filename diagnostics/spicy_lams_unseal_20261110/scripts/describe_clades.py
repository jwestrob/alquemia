"""Post hoc description: which clades sit where on the tilt axis, the extremes, and tilt against loop acidity
within clades (improved loop finder)."""
import sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
LOG = open(ROOT / "results/describe_clades.log", "w")
def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
d, raw, C, mean, valid, qc = load()
L = pd.read_csv(ROOT / "results/loops_improved.tsv", sep="\t")
ix = {e: i for i, e in enumerate(EL)}
t = pd.DataFrame(dict(accession=d.accession, family=d.family_original.fillna("unspecified"), genus=d.genus_original.fillna("unspecified"), efh=d.EFhands, qc=qc,
                      mtot=L.mtot, acid=L.acid2, acid_list=L.acid2_list, K2=L.n_K2_2, P2=L.n_P2_2, D9=L.n_D9_2, loops=L.loops2,
                      Lu_La=mean[:, ix["Lu"]] - mean[:, ix["La"]], Dy_Nd=mean[:, ix["Dy"]] - mean[:, ix["Nd"]], Nd_Pr=mean[:, ix["Nd"]] - mean[:, ix["Pr"]],
                      Dy_Tb=mean[:, ix["Dy"]] - mean[:, ix["Tb"]], Y_Dy=mean[:, ix["Y"]] - mean[:, ix["Dy"]], La_Nd=mean[:, ix["La"]] - mean[:, ix["Nd"]]))
q = t[t.qc].copy()
p("By family (QC orthologs; log10 units; mean and SD):")
g = q.groupby("family").agg(n=("Lu_La", "size"), Lu_La=("Lu_La", "mean"), Lu_La_sd=("Lu_La", "std"), Dy_Nd=("Dy_Nd", "mean"), Dy_Nd_sd=("Dy_Nd", "std"),
                            La_Nd=("La_Nd", "mean"), Y_Dy=("Y_Dy", "mean"), acid=("acid", "mean"), capture=("mtot", "median")).sort_values("n", ascending=False)
p(g[g.n >= 3].round(3).to_string())
p("\nBy genus (QC orthologs, genera with >= 3):")
g2 = q.groupby("genus").agg(n=("Lu_La", "size"), Lu_La=("Lu_La", "mean"), Dy_Nd=("Dy_Nd", "mean"), Dy_Nd_sd=("Dy_Nd", "std"), La_Nd=("La_Nd", "mean"), acid=("acid", "mean"), capture=("mtot", "median")).sort_values("Dy_Nd")
p(g2[g2.n >= 3].round(3).to_string())
p(f"\nShare of between-ortholog variance in Dy-Nd explained by genus (QC, genera >= 3): "
  f"{1 - q[q.genus.map(q.genus.value_counts()) >= 3].groupby('genus').Dy_Nd.apply(lambda v: ((v - v.mean()) ** 2).sum()).sum() / ((q[q.genus.map(q.genus.value_counts()) >= 3].Dy_Nd - q[q.genus.map(q.genus.value_counts()) >= 3].Dy_Nd.mean()) ** 2).sum():.2f}")
p("\nMost Nd-over-Dy selective orthologs (QC):")
p(q.sort_values("Dy_Nd").head(8)[["accession", "family", "genus", "Dy_Nd", "Lu_La", "La_Nd", "mtot", "acid_list", "loops"]].round(3).to_string(index=False))
p("\nLeast Nd-over-Dy selective orthologs (QC):")
p(q.sort_values("Dy_Nd").tail(6)[["accession", "family", "genus", "Dy_Nd", "Lu_La", "La_Nd", "mtot", "acid_list", "loops"]].round(3).to_string(index=False))
p(f"\nAny QC ortholog that prefers heavy over light (Lu-La > 0)? {int((q.Lu_La > 0).sum())} of {len(q)}; that prefers Dy over Nd? {int((q.Dy_Nd > 0).sum())}")
# tilt against acidity: pooled, within family, within genus
def within(df, col, by, minn=4):
    out = pd.Series(np.nan, index=df.index)
    for k, s in df.groupby(by):
        if len(s) >= minn and k != "unspecified": out[s.index] = s[col] - s[col].mean()
    return out
for tgt in ("Lu_La", "Dy_Nd"):
    line = f"{tgt} vs acidic per loop: pooled {spearmanr(q.acid, q[tgt], nan_policy='omit').statistic:+.2f}"
    for by in ("family", "genus"):
        a, b = within(q, "acid", by), within(q, tgt, by); m = a.notna() & b.notna()
        line += f"; within {by} {spearmanr(a[m], b[m]).statistic:+.2f} (n={int(m.sum())})"
    p(line)
t.to_csv(ROOT / "results/ortholog_contrasts.tsv", sep="\t", index=False, float_format="%.4f")
