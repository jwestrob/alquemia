"""MBAR / TI for the La->Dy windows, then dG(T) -> dH, dS for each leg and the exchange cycle.

Leg (aq, ac, edta): dG_alch(T) = G(Dy) - G(La) in that environment. Exchange
La.L + Dy(aq) -> Dy.L + La(aq):  ddG = dG_alch(L) - dG_alch(aq)  (negative = Dy favoured).
dH and dS from a weighted linear fit dG(T) = dH - T dS over the simulated temperatures
(constant dCp neglected; checked against the curvature when 4 T are available).

usage: analyze_fep.py <workspace root>   (expects <leg>/fep/T<T>/s<seed>/lam<l>/rec.npy)
Seeds are pooled per window for MBAR; the per-seed spread is reported as a check on the
MBAR error. The simulation omits the metal-OW tail beyond 1.0 nm (custom force, no long-range
correction); it is added analytically per frame: dU_tail = -4 pi N_OW/V (dC/rc + dB/(3 rc^3)).
Run with the lanm_bench python (pymbar 4).
"""
import json, sys
from pathlib import Path
import numpy as np
from pymbar import timeseries
from pymbar.other_estimators import bar

R_KJ = 0.0083144626
KCAL = 4.184
DISCARD_PS = 200  # extra production discard on top of the 0.5 ns equilibration
RC_A = 10.0
# metal-OW 12-6-4 coefficients (kcal/mol, A) from the parmed add12_6_4 prmtops (OPC)
B_OW = {"La": 655.403423, "Dy": 416.617688}; C_OW = {"La": 165.0, "Dy": 243.0}
N_OW = {"aq": 1256, "ac": 1857, "edta": 2205}


def tail_kcal(V, leg):
    rho = N_OW[leg] / V
    return -4 * np.pi * rho * ((C_OW["Dy"] - C_OW["La"]) / RC_A + (B_OW["Dy"] - B_OW["La"]) / (3 * RC_A ** 3))


def load_windows(d, seeds):
    """{lambda: concatenated records} over the requested seed directories."""
    W = {}
    for sd in seeds:
        for p in (d / sd).glob("lam*/rec.npy"):
            if not (p.parent / "meta.json").exists(): continue
            rec = np.load(p)[DISCARD_PS:]
            if len(rec) < 200: continue
            W.setdefault(float(p.parent.name[3:]), []).append(rec)
    return {k: np.concatenate(v) for k, v in sorted(W.items())}


def leg_T(d, leg, seeds=None):
    seeds = seeds or sorted(x.name for x in d.glob("s*"))
    W = load_windows(d, seeds)
    lams, dus, struct = [], [], {}
    tails = []
    for lam, rec in W.items():
        tails.append(tail_kcal(rec[:, 2], leg).mean())
        g = timeseries.statistical_inefficiency(rec[:, 0])
        idx = timeseries.subsample_correlated_data(rec[:, 0], g=g)
        lams.append(lam); dus.append(rec[idx, 0])
        if lam in (0.0, 1.0):
            cn = rec[:, 3].astype(int)
            struct["La" if lam == 0 else "Dy"] = {
                "CN_mean": float(cn.mean()), "CN_pop": {int(k): float((cn == k).mean()) for k in np.unique(cn)},
                "MO_mean_A": float(np.nanmean(rec[:, 4])),
                "ligand_donors_pop": ({int(k): float((rec[:, 5].astype(int) == k).mean()) for k in np.unique(rec[:, 5].astype(int))}
                                      if rec.shape[1] > 5 else None)}
    if len(lams) < 3 or 0.0 not in lams or 1.0 not in lams or max(np.diff(lams)) > 0.09:
        return None  # incomplete lambda ladder (gaps break MBAR overlap)
    T = float(d.name[1:]); beta = 1.0 / (R_KJ * T)
    lams = np.round(np.array(lams) * 12) / 12  # dir names are rounded to 3 dp; the simulated lambda is i/12
    N = np.array([len(x) for x in dus])
    # U is linear in lambda, so the reduced work between neighbours is beta*dlam*dU/dlam; sum BAR over neighbours
    # (pymbar MBAR on the pooled u_kn failed its normalization check)
    dG, var, ov = 0.0, 0.0, []
    for i in range(len(lams) - 1):
        dl = lams[i + 1] - lams[i]
        wF = beta * dl * dus[i]; wR = -beta * dl * dus[i + 1]
        b = bar(wF, wR)
        dG += b["Delta_f"]; var += b["dDelta_f"] ** 2
        ov.append(float(np.mean(np.minimum(1, np.exp(-wF + b["Delta_f"])))))  # forward acceptance as an overlap proxy
    dG = dG / beta / KCAL; ddG = np.sqrt(var) / beta / KCAL
    ti = np.trapezoid([x.mean() for x in dus], lams) / KCAL
    tail = float(np.trapezoid(tails, lams))
    dG += tail; ti += tail
    min_ov = float(min(ov))
    return {"T": T, "seeds": seeds, "dG_kcal": dG, "dG_err": ddG, "dG_TI_kcal": ti, "tail_kcal": tail, "n_windows": len(lams),
            "n_uncorrelated_min": int(N.min()), "min_neighbor_overlap": min_ov, "structure": struct}


def fit(rows):
    T = np.array([r["T"] for r in rows]); G = np.array([r["dG_kcal"] for r in rows]); s = np.array([max(r["dG_err"], 1e-3) for r in rows])
    if len(T) < 2: return None
    A = np.vstack([np.ones_like(T), -T]).T; W = np.diag(1 / s ** 2)
    cov = np.linalg.inv(A.T @ W @ A); dH, dS = cov @ A.T @ W @ G
    return {"dH_kcal": float(dH), "dH_err": float(np.sqrt(cov[0, 0])), "dS_kcal_per_K": float(dS),
            "TdS_298_kcal": float(298.15 * dS), "TdS_298_err": float(298.15 * np.sqrt(cov[1, 1])),
            "dG_298_kcal": float(dH - 298.15 * dS), "dG_298_err": float(np.sqrt(cov[0, 0] + 298.15 ** 2 * cov[1, 1] - 2 * 298.15 * cov[0, 1])),
            "resid_kcal": [float(x) for x in G - (dH - T * dS)]}


def main():
    root = Path(sys.argv[1])
    res = {}
    for leg in ("aq", "ac", "edta"):
        rows = [x for x in (leg_T(d, leg) for d in sorted((root / leg / "fep").glob("T*"))) if x]
        for r in rows:  # per-seed check
            d = root / leg / "fep" / f"T{r['T']:.0f}"
            per = [leg_T(d, leg, [sd]) for sd in sorted(x.name for x in d.glob("s*"))]
            r["per_seed_dG"] = [x["dG_kcal"] for x in per if x]
        if rows:
            res[leg] = {"per_T": rows, "fit": fit(rows)}
    for leg in ("ac", "edta"):
        if leg in res and "aq" in res:
            Ts = sorted(set(r["T"] for r in res[leg]["per_T"]) & set(r["T"] for r in res["aq"]["per_T"]))
            rows = []
            for T in Ts:
                a = next(r for r in res[leg]["per_T"] if r["T"] == T); b = next(r for r in res["aq"]["per_T"] if r["T"] == T)
                rows.append({"T": T, "dG_kcal": a["dG_kcal"] - b["dG_kcal"], "dG_err": float(np.hypot(a["dG_err"], b["dG_err"]))})
            res[f"exchange_{leg}"] = {"per_T": rows, "fit": fit(rows)}
    (root / "FEP_RESULT.json").write_text(json.dumps(res, indent=1))
    for k, v in res.items():
        print(k, " ".join(f"T{r['T']:.0f}:{r['dG_kcal']:+.2f}±{r['dG_err']:.2f}" for r in v["per_T"]))
        f = v["fit"]
        if f: print(f"   dG298 {f['dG_298_kcal']:+.2f}±{f['dG_298_err']:.2f}  dH {f['dH_kcal']:+.2f}±{f['dH_err']:.2f}  TdS {f['TdS_298_kcal']:+.2f}±{f['TdS_298_err']:.2f}")
        for r in v["per_T"]:
            if "structure" in r: print(f"   T{r['T']:.0f} overlap_min {r['min_neighbor_overlap']:.3f} TI {r['dG_TI_kcal']:+.2f} {json.dumps(r['structure'])}")


if __name__ == "__main__":
    main()
