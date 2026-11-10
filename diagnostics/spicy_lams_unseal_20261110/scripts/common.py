"""Shared loading and feature code for the lanmodulin ortholog profile analysis (PLAN.md, frozen 2026-11-10)."""
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EL = ["Y", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
LN = EL[1:]
R9 = np.array([1.216, 1.196, 1.179, 1.163, 1.132, 1.120, 1.107, 1.095, 1.083, 1.072, 1.062, 1.052, 1.042, 1.032])
R9_Y = 1.075
KC = 1.3643
SEQ = "Sequence..SP.Removed."


def load():
    d = pd.read_csv(ROOT / "data/nonreserved_profiles.tsv", sep="\t")
    raw = np.stack([d[[f"{e}-{k}" for e in EL]].apply(pd.to_numeric, errors="coerce").values for k in (1, 2, 3)], axis=1)  # proteins x 3 x 15
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.where(raw > 0, np.log10(raw), np.nan)
    valid = ~np.isnan(L).any(axis=2)                       # replicate has all 15 elements positive
    C = L - np.nanmean(L, axis=2, keepdims=True)           # centred replicate profiles
    C[~valid] = np.nan
    mean = np.nanmean(C, axis=1)                           # proteins x 15
    nrep = valid.sum(1)
    mtot = pd.to_numeric(d["Total_Metal_Replicate_Average"], errors="coerce").values
    rsd = pd.to_numeric(d["Total_Metal_Replicate_RSD"], errors="coerce").values
    qc = (mtot >= 20) & (rsd <= 50) & (nrep >= 2)
    return d, raw, C, mean, valid, qc


def poly_design(deg=3, with_y=False):
    """Orthonormal polynomial components on the 14 lanthanides; x runs from -1 (La) to +1 (Lu).
    Returns Q (14 x deg+1) and a function giving the same components at any radius."""
    x = (R9[0] - R9) / (R9[0] - R9[-1]) * 2 - 1
    V = np.vstack([x ** k for k in range(deg + 1)]).T
    Q, Rm = np.linalg.qr(V)
    sgn = np.sign(np.diag(Rm)); Q = Q * sgn; Rm = (Rm.T * sgn).T
    def at(r):
        xx = (R9[0] - np.atleast_1d(r)) / (R9[0] - R9[-1]) * 2 - 1
        return np.vstack([xx ** k for k in range(deg + 1)]).T @ np.linalg.inv(Rm)
    return Q, at


def components(prof15, deg=3):
    """prof15: (..., 15) centred log10 profile (Y first). Returns coefficients (..., deg+1) on the orthonormal
    components fitted to the 14 lanthanides, the yttrium offset (Y minus the fitted curve at 1.075 A), and the
    residual sum of squares over the 14 lanthanides."""
    Q, at = poly_design(deg)
    y = prof15[..., 1:]
    co = y @ Q
    fit = co @ Q.T
    yoff = prof15[..., 0] - (co @ at(R9_Y).T)[..., 0]
    return co, yoff, ((y - fit) ** 2).sum(-1)


LOOP = re.compile(r"(?=(D.[DNS].[DNSTG]......[DE]))")


def ef_loops(seq):
    out = []; last = -99
    for m in LOOP.finditer(seq):
        if m.start() >= last + 12:
            out.append((m.start(), m.group(1))); last = m.start()
    return out


def loop_features(seq):
    lp = ef_loops(seq)
    n = len(lp)
    acid = [sum(l[i] in "DE" for i in (0, 2, 4, 8, 11)) for _, l in lp]
    return dict(n_loops=n, acid_per_loop=float(np.mean(acid)) if n else np.nan, n_D9=sum(l[8] == "D" for _, l in lp),
                n_P2=sum(l[1] == "P" for _, l in lp), acid_total=int(sum(acid)),
                net_charge_per_res=(sum(seq.count(a) for a in "KR") - sum(seq.count(a) for a in "DE")) / len(seq),
                length=len(seq), loops=";".join(f"{s}:{l}" for s, l in lp))
