"""Share of small-chelator selectivity between two metals that comes from the Gd|Tb gap (71 NIST ligands,
cubic polynomial per ligand in the two-branch coordinate, gap fixed at 31 mA). Post hoc; lead agent."""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa
Yall, _, meta = load_nist(10); GAP = 0.031; RS = R9 + GAP * (np.arange(14) <= 6)
zz = lambda x: (x - X0) / XS
rows = []
for a, b in (("La", "Dy"), ("Nd", "Dy"), ("La", "Lu"), ("Eu", "Tb"), ("Gd", "Tb")):
    sh, tot = [], []
    for y in Yall:
        ok = ~np.isnan(y); c = np.linalg.lstsq(poly_basis(RS, 3)[ok], y[ok], rcond=None)[0]
        f = lambda x: sum(c[k] * zz(x) ** k for k in range(4))
        ia, ib = LN.index(a), LN.index(b)
        full = f(R9[ib]) - f(R9[ia] + GAP)
        nog = 0.5 * ((f(R9[ib]) - f(R9[ia])) + (f(R9[ib] + GAP) - f(R9[ia] + GAP)))   # gap closed from either side
        if abs(full) >= 0.5:
            sh.append((full - nog) / full); tot.append(full)
    sh = np.array(sh); tot = np.array(tot)
    rows.append(dict(pair=f"{a}->{b}", n_ligands=len(sh), median_share=np.median(sh), q25=np.percentile(sh, 25), q75=np.percentile(sh, 75),
                     median_abs_ddG=np.median(np.abs(tot)), median_gap_part=np.median(np.abs(tot * sh))))
t = pd.DataFrame(rows); t.to_csv(ROOT / "results/gap_share.tsv", sep="\t", index=False, float_format="%.3f"); print(t.round(3).to_string(index=False))
