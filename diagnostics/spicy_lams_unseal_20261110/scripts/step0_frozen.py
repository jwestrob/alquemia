"""Step 0: score the frozen predictions P1 and P2 with the hashed script, unchanged, on the QC orthologs."""
import hashlib, subprocess, sys
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
d, raw, C, mean, valid, qc = load()
print(f"rows {len(d)}; replicates valid per protein: {dict(pd.Series(valid.sum(1)).value_counts().sort_index())}; QC set {int(qc.sum())}")
# sanity: my centred mean profile against the workbook's normalized logD mean (centred)
nl = d[[f"nlogD_mean_{e}" for e in EL]].apply(pd.to_numeric, errors="coerce").values
nlc = nl - np.nanmean(nl, axis=1, keepdims=True)
ok = qc & ~np.isnan(nlc).any(1)
r = [np.corrcoef(mean[i], nlc[i])[0, 1] for i in np.where(ok)[0]]
print(f"sanity: per-protein correlation of my centred mean profile with the workbook's normalized logD (centred): median {np.median(r):.4f}, 5th pct {np.percentile(r, 5):.4f}; "
      f"median absolute difference {np.nanmedian(np.abs(mean[ok] - nlc[ok])):.4f} log10 units")
t = pd.DataFrame(mean[qc], columns=EL); t.insert(0, "id", d.accession.values[qc])
t.to_csv(ROOT / "results/qc_mean_profiles_log10.tsv", sep="\t", index=False, float_format="%.5f")
frozen = "/Users/jacob/Documents/Berkeley/Lanthanides/ln_series_coordinate_2026-11-10/scripts/sealed_gap_test.py"
h = hashlib.sha256(open(frozen, "rb").read()).hexdigest(); assert h.startswith("587d1604"), h
out = subprocess.run([sys.executable, "-W", "ignore", frozen, str(ROOT / "results/qc_mean_profiles_log10.tsv"), "--scale", "log10", "--sign", "assoc"],
                     capture_output=True, text=True)
txt = "\n".join(l for l in out.stdout.splitlines() if "RuntimeWarning" not in l)
open(ROOT / "results/step0_frozen_predictions.log", "w").write(f"frozen script sha256 {h}\n" + txt + "\n")
print(txt); print(out.stderr[-500:] if out.returncode else "")
