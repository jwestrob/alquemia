"""Part E of PLAN.md: transfer of the NIST-derived shared structure to protein and peptide series.
Dataset builder copied from prior_scripts/protein_series_shapes.py (same grouping key, same filters).
Nothing shared is refitted on the protein data. Lead agent, Claude Fable."""
import json
import re
import sys
import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa

RES = ROOT / "results"


def build():
    d = pd.read_csv(ROOT / "data/protein_ln_affinities.tsv", sep="\t")
    d = d[d.metal.isin(LN)].copy()
    m = (d.system == "Mex-LanM (agarose-immobilized)") & d.log_K.isna()
    d.loc[m, "log_K"] = [np.log10(float(re.search(r"/Keq_rel\(La\) = ([0-9.eE+-]+)", n).group(1))) / 2
                         if "/Keq_rel(La) =" in n else (0.0 if mm == "La" else np.nan)
                         for n, mm in zip(d.loc[m, "notes"], d.loc[m, "metal"])]
    d = d[d.log_K.notna() & ~d.constant_type.isin(["Km", "EC50_or_K_half_activity"])]
    key = ["system", "variant", "doi", "method", "pH", "T_C", "constant_type", "site", "competitor",
           "buffer_ionic_strength"]
    out = []
    for k_, s in d.fillna({c: "-" for c in key}).groupby(key):
        s = s.groupby("metal").log_K.mean()
        if len(s) < 6:
            continue
        y = np.full(14, np.nan); y[[LN.index(x) for x in s.index]] = -KC * s.values
        out.append(dict(dataset=f"{k_[0][:34]} | {k_[1][:22]}", system=k_[0], doi=k_[2], y=y))
    return out


def main():
    tab = pd.read_csv(RES / "coordinate.tsv", sep="\t")
    Y, _, _ = load_nist(14); ybar = Y.mean(0)
    rhoC, rhoOC = tab.rho_C.values, tab.rho_OC.values
    Z0 = np.zeros(14)
    mods = {"Q": (poly_basis(R9), Z0), "O": (poly_basis(R9), ybar), "C": (poly_basis(rhoC), Z0),
            "OC": (poly_basis(rhoOC), ybar), "K": (poly_basis(R9, 3), Z0)}
    rows = []
    for ds in build():
        y = ds["y"]; idx = np.where(~np.isnan(y))[0]
        rec = dict(dataset=ds["dataset"], doi=ds["doi"], n=len(idx), span=float(np.nanmax(y) - np.nanmin(y)))
        for k, (B, c) in mods.items():
            if len(idx) - 1 <= B.shape[1]:
                rec[f"lomo_{k}"] = np.nan; rec[f"in_{k}"] = np.nan; continue
            e = []
            for m_ in idx:
                fi = idx[idx != m_]
                e.append(y[m_] - predict(B, c, y, fi)[0][m_])
            rec[f"lomo_{k}"] = float(np.sqrt(np.mean(np.square(e))))
            rec[f"in_{k}"] = float(np.sqrt(np.mean((y[idx] - predict(B, c, y, idx)[0][idx]) ** 2)))
        rows.append(rec)
    R = pd.DataFrame(rows).sort_values("n", ascending=False)
    R.to_csv(RES / "protein_transfer.tsv", sep="\t", index=False, float_format="%.4f")
    pd.set_option("display.width", 250)
    print(R.round(3).to_string(index=False))
    out = {}
    for k in ("O", "C", "OC", "K"):
        for which in ("lomo", "in"):
            d = (R[f"{which}_Q"] - R[f"{which}_{k}"]).dropna().values
            pos = int((d > 0).sum()); nz = int((d != 0).sum())
            pv = binomtest(pos, nz, 0.5).pvalue if nz else np.nan
            try:
                wv = wilcoxon(d).pvalue
            except Exception:
                wv = np.nan
            out[f"{which}_Q_minus_{k}"] = dict(median=float(np.median(d)), mean=float(np.mean(d)), n_pos=pos, n=nz,
                                               sign_p=float(pv), wilcoxon_p=float(wv),
                                               median_Q=float(np.median(R[f"{which}_Q"].dropna())),
                                               median_model=float(np.median(R[f"{which}_{k}"].dropna())))
            print(f"{which:4s} Q - {k:2s}: median {np.median(d):+.4f}, mean {np.mean(d):+.4f}, better in {pos}/{nz}, "
                  f"sign test p = {pv:.3f}, Wilcoxon p = {wv:.3f}")
    # one dataset per DOI-and-system (first listed = most metals) as a dependence check
    R1 = R.drop_duplicates(subset=["doi"]).copy()
    for k in ("O", "C", "OC"):
        d = (R1["lomo_Q"] - R1[f"lomo_{k}"]).dropna().values
        pos = int((d > 0).sum()); nz = int((d != 0).sum())
        out[f"lomo_Q_minus_{k}_one_per_doi"] = dict(median=float(np.median(d)), n_pos=pos, n=nz,
                                                    sign_p=float(binomtest(pos, nz, 0.5).pvalue))
        print(f"one dataset per DOI, lomo Q - {k}: median {np.median(d):+.4f}, better in {pos}/{nz}, "
              f"sign test p = {binomtest(pos, nz, 0.5).pvalue:.3f}")
    json.dump(out, open(RES / "protein_transfer_summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
