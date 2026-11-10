"""Step 0 of PLAN.md: reference distributions of the single-pattern fingerprint for proteins and chelators."""
import json, re, sys
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fingerprint import *  # noqa
ROOT = Path(__file__).resolve().parents[1]; B = Path("/Users/jacob/Documents/Berkeley/Lanthanides")
LOG = open(ROOT / "results/step0.log", "w")
def p(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
rows = []
lm = pd.read_csv(B / "lanm_ortholog_profiles_2026-11-10/results/qc_mean_profiles_log10.tsv", sep="\t")
for r in lm.itertuples():
    v = {e: getattr(r, e) for e in LN}
    for vn, ex in VARIANTS.items():
        f = fingerprint(v, r.Y, ex)
        if f: rows.append(dict(cls="lanmodulin library", id=r.id, variant=vn, donors="protein", **f))
pp = pd.read_csv(B / "ln_series_coordinate_2026-11-10/results/public_protein_profiles.tsv", sep="\t")
for r in pp.itertuples():
    v = {e: getattr(r, e) for e in LN}
    for vn, ex in VARIANTS.items():
        f = fingerprint(v, r.Y, ex)
        if f: rows.append(dict(cls=("siderophore DFOB" if "DFOB" in r.id else "public protein or peptide"), id=r.id, variant=vn, donors="protein", **f))
ns = pd.read_csv(B / "ln_series_coordinate_2026-11-10/data/surrogate_nist.tsv", sep="\t", dtype={"lig_id": str})
for r in ns.itertuples():
    k = {a: float(b) for a, b in (x.split(":") for x in r.logK_series.split(";"))}
    if sum(e in k for e in LN) < 10: continue
    nN = int(m.group(1)) if (m := re.search(r"N(\d+)", r.formula)) else 0
    for vn, ex in VARIANTS.items():
        f = fingerprint(k, k.get("Y"), ex)
        if f: rows.append(dict(cls="chelator", id=r.name[:40], variant=vn, donors=("O only" if nN == 0 else "with N"), **f))
T = pd.DataFrame(rows); T.to_csv(ROOT / "results/reference_fingerprints.tsv", sep="\t", index=False, float_format="%.4f")
def auc(a, b):
    a = np.asarray(a)[~np.isnan(a)]; b = np.asarray(b)[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0: return np.nan
    return float((np.sum(a[:, None] > b[None, :]) + 0.5 * np.sum(a[:, None] == b[None, :])) / (len(a) * len(b)))
SUM = {}
for vn in VARIANTS:
    t = T[T.variant == vn]; p(f"\n=== variant {vn} ===")
    for c, g in t.groupby("cls"):
        for d, gg in g.groupby("donors"):
            q = lambda x: f"{np.nanmedian(x):+.3f} [{np.nanpercentile(x, 25):+.3f}, {np.nanpercentile(x, 75):+.3f}]" if np.isfinite(x).any() else "n/a"
            p(f"  {c:26s} {d:8s} n={len(gg):3d} | step h {q(gg.h.values)} | slope per step {q(gg.slope_step.values)} | u {q(gg.u.values)} (defined {int(gg.u.notna().sum())}) | "
              f"yttrium offset {q(gg.y.values) if 'y' in gg else 'n/a'} (n {int(gg.y.notna().sum()) if 'y' in gg else 0}) | fit rms {np.nanmedian(gg.rms):.3f}")
    lmn = t[t.cls == "lanmodulin library"]; ch = t[t.cls == "chelator"]; cho = ch[ch.donors == "O only"]
    a_u = auc(ch.u.values, lmn.u.values); a_h = auc(np.abs(ch.h.values), np.abs(lmn.h.values))
    a_y = auc(ch.y.values, lmn.y.values); a_yo = auc(cho.y.values, lmn.y.values)
    p(f"  single-pattern separation, chelators against lanmodulins (AUC): u {a_u:.3f}; |h| {a_h:.3f}; yttrium offset {a_y:.3f} (oxygen-only chelators {a_yo:.3f})")
    SUM[vn] = dict(auc_u=a_u, auc_abs_h=a_h, auc_y=a_y, auc_y_O=a_yo)
ok = any(v["auc_u"] >= 0.8 or v["auc_y"] >= 0.8 for v in SUM.values())
p(f"\nstop rule of step 0 (AUC >= 0.8 for u or for the yttrium offset in some variant): {'continue' if ok else 'STOP'}")
json.dump(SUM, open(ROOT / "results/step0.json", "w"), indent=1)
