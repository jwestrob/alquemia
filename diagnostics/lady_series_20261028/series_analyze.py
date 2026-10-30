"""Series extension: La -> M exchange (M = Nd, Dy, Lu) with the careful cluster-continuum protocol, against NIST 46.

Reaction (negative = M favoured):  La.L(H2O)b + M(H2O)c -> M.L(H2O)a + La(H2O)d + n H2O
States: La/Dy aqua, acetate, EDTA from lady_relaxed_exchange_20261024 (conformer-best or primary-only);
La/Dy malonate..TMDTA from lady_panel_20261028; Nd/Lu everything from lady_series_20261028 (primary only, started
from the conformer-best La/Dy analogue). Same state rule as the panel: lowest-G hydration state per metal, water-balanced.
Panel and series states whose coordination number changed in DFT (CN < 3.1 A != donors + waters) are excluded.

Pre-registered (2026-10-28, before any Nd/Lu job ran; see README.md):
  S-H0  per metal, constant offset from acetate: RMS residual over the six other ligands < 0.8 kcal/mol, for each of
        Nd, Dy, Lu.
  S-H2  per metal, offset + b*n_disp calibrated on acetate + EDTA: RMS over the five other ligands < 0.8.
  Mono  aqua-leg prediction: offsets grow more negative with contraction, c(Nd) > c(Dy) > c(Lu).
  Whole series: calibrated calc (Nd, Dy, Lu) - La mapped to all 13 metals with the NIST basis (map fitted without
        that ligand); RMS vs NIST over every ligand x metal < 0.8 kcal/mol.
Added 2026-10-29, before any Nd/Lu job ran (H0 had failed on La/Dy; see ../lady_panel_20261028/README.md):
  TTHA (NIST 678) joins the ligands; exempt from the CN test, as in the panel.
  S-H3  per metal, stretched selectivity calc = alpha*exp + beta, calibrated on acetate + EDTA: RMS (exp units) over
        the other ligands < 0.8. S-H2 and S-H3 scores now exclude both calibration ligands.
  Whole series uses the form the panel's TTHA call selected (PANEL_RESULT.json; per-water -> H2, stretch/constant ->
        H3, since H0 already failed); calibration ligands are excluded from its score. Other forms are printed too.
  Secondary (reported, not a gate): the same forms calibrated on the panel's four fit ligands (ac, EDTA, NTA, DTPA)
        and scored on its four held-out ligands (mal, CDTA, TMDTA, TTHA). For Nd, acetate and EDTA differ by only
        1.2 kcal/mol experimentally, so the two-point H3 slope is poorly conditioned there.
"""
import itertools, json, re, sys
from pathlib import Path
import numpy as np, pandas as pd

D = Path(__file__).resolve().parent
REPO = D.parents[1]
sys.path.insert(0, str(REPO / "diagnostics/lady_relaxed_exchange_20261024"))
import analyze as A  # noqa: E402

PANEL = REPO / "workspaces/lady_panel_20261028/runs"
SERIES = REPO / "workspaces/lady_series_20261028/runs"
H, WC, KC = A.H, A.WATER_CORR, 1.3643
LN = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
METALS = ["Nd", "Dy", "Lu"]
R9 = {"La": 1.216, "Nd": 1.163, "Dy": 1.083, "Lu": 1.032}  # Shannon CN9
LIG_ID = {"edtax": "576", "mal": "3480", "nta": "465", "dtpa": "655", "cdta": "600", "tmdta": "610", "ttha": "678"}
Q_C = {"ac": 2, "edtax": -1, "mal": 1, "nta": 0, "dtpa": -2, "cdta": -1, "tmdta": -1, "ttha": -3}
DONORS = {"aq": 0, "acmono": 1, "acbi": 2, "edtax": 6, "mal": 2, "nta": 4, "dtpa": 8, "cdta": 6, "tmdta": 6, "ttha": 9}
NO_CN_CHECK = {"ttha"}  # no water to lose; 9 vs 10 bound donors are isomers of the same q = 0 state
CAL = {"primary": ("ac", "edtax"), "secondary": ("ac", "edtax", "nta", "dtpa")}  # calibration ligands per metal
HELD_SECONDARY = ("mal", "cdta", "tmdta", "ttha")
FORM_OF_CALL = {"per-water (H2)": "H2", "stretch/constant (H3/H0)": "H3"}


def exp_series():
    """{lig: {M: ddG(M - La) kcal/mol}} from NIST 46 (same sources/ionic strengths as the La/Dy panel)."""
    out = {}
    ac = {}
    for line in (REPO / "workspaces/lady_relaxed_exchange_20261024/src_untrusted/lit/nist46_acetate_ln.txt").read_text().splitlines():
        f = line.split()
        if len(f) >= 8 and f[3] == "[ML]/[M][L]" and f[4] == "K" and f[5] == "25" and f[6] == "0.1":
            ac[f[1].split("/")[0]] = float(f[7])
    out["ac"] = {m: -KC * (ac[m] - ac["La"]) for m in ac if m in LN and m != "La"}
    t = pd.read_csv(REPO / "workspaces/ln_surrogate_20261027/surrogate_nist.tsv", sep="\t", dtype={"lig_id": str})
    for lig, lid in LIG_ID.items():
        k = dict(x.split(":") for x in t[t.lig_id == lid].iloc[0].logK_series.split(";"))
        out[lig] = {m: -KC * (float(k[m]) - float(k["La"])) for m in k if m in LN and m != "La"}
    return out


def nist_complete():
    t = pd.read_csv(REPO / "workspaces/ln_surrogate_20261027/surrogate_nist.tsv", sep="\t", dtype={"lig_id": str})
    ids, X = [], []
    for _, r in t.iterrows():
        k = dict(x.split(":") for x in r.logK_series.split(";"))
        if all(m in k for m in LN):
            ids.append(r.lig_id); X.append([-KC * (float(k[m]) - float(k["La"])) for m in LN[1:]])
    return ids, np.array(X)


def kind(name):
    return "aq" if "_aq" in name else name.split("_")[1]


def load_new(root, metals):
    R, status = {}, {}
    for p in sorted(root.iterdir()):
        if p.name.split("_")[0] not in metals or not (p / "job.inp").exists():
            continue
        r = A.parse(p.name, p)
        if r is None:
            status[p.name] = "no output"; continue
        cn_all = sum(1 for _, x in r.get("donors", []) if x < 3.1)
        exp_cn = DONORS[kind(p.name)] + A.nw(p.name)
        status[p.name] = {"normal": r["normal"], "imag": r["imag"], "CN_3.1A": cn_all, "CN_expected": exp_cn}
        if r["normal"] and r["G"] is not None and (cn_all == exp_cn or kind(p.name) in NO_CN_CHECK):
            R[p.name] = r
        elif r["normal"]:
            print(f"EXCLUDED {p.name}: CN {cn_all} != expected {exp_cn}")
    return R, status


def exchange(ok, lig, M):
    gw = ok["water"]
    pool = (lambda n, m: A.pool_match(n, m, "ac")) if lig == "ac" else (lambda n, m: n.startswith(f"{m}_{lig}_"))
    sel = {}
    for m in ("La", M):
        for fam, p in (("L", [n for n in ok if pool(n, m)]), ("aq", [n for n in ok if n.startswith(f"{m}_aq")])):
            if not p:
                return None
            ref = max(A.nw(x) for x in p)
            sel[(m, fam)] = min(p, key=lambda x: ok[x]["G"] * H + (ref - A.nw(x)) * (gw["G"] * H + WC))
    laL, mL, laaq, maq = sel[("La", "L")], sel[(M, "L")], sel[("La", "aq")], sel[(M, "aq")]
    n = A.nw(laL) + A.nw(maq) - A.nw(mL) - A.nw(laaq)
    d = lambda k: H * (ok[mL][k] + ok[laaq][k] - ok[laL][k] - ok[maq][k] + n * gw[k])
    return {"dG": d("G") + n * WC, "dH": d("H"), "states": [laL, mL, laaq, maq],
            "n_disp": ((A.nw(laaq) - A.nw(laL)) + (A.nw(maq) - A.nw(mL))) / 2}


def calibrate(r, M, EXP, ligs):
    """Per-metal maps calc -> calibrated ddG (experimental units) for H0 (acetate offset), H2 (offset + b*n_disp) and
    H3 (calc = alpha*exp + beta), the two-parameter forms least-squares fitted on `ligs` (exact for two ligands)."""
    off = r["ac"]["dG"] - EXP["ac"][M]
    cal, par = {"H0": lambda l: r[l]["dG"] - off}, {"H0": {"offset": off}}
    L = [l for l in ligs if l in r]
    if len(L) >= 2 and len(L) == len(ligs):
        Xn = np.array([[1.0, r[l]["n_disp"]] for l in L]); yb = np.array([r[l]["dG"] - EXP[l][M] for l in L])
        c0, b = np.linalg.lstsq(Xn, yb, rcond=None)[0]
        Xe = np.array([[EXP[l][M], 1.0] for l in L]); yc = np.array([r[l]["dG"] for l in L])
        al, be = np.linalg.lstsq(Xe, yc, rcond=None)[0]
        cal["H2"] = lambda l: r[l]["dG"] - c0 - b * r[l]["n_disp"]
        par["H2"] = {"c0": c0, "b": b}
        if al > 0.1:  # a stretch must be positive; a near-zero slope cannot be inverted
            cal["H3"] = lambda l: (r[l]["dG"] - be) / al
            par["H3"] = {"alpha": al, "beta": be}
        else:
            print(f"   S-H3 [{'+'.join(ligs)}] undefined for {M}: fitted alpha {al:+.3f}")
    return cal, par


def panel_call(tag):
    f = PANEL.parent / "PANEL_RESULT.json"
    try:
        return json.loads(f.read_text())["results"][tag + "_summary"]["frozen"].get("ttha_call")
    except (OSError, KeyError, ValueError):
        return None


def rms(v):
    return float(np.sqrt(np.mean(np.square(list(v))))) if len(v) else float("nan")


def main():
    EXP = exp_series()
    ids, NX = nist_complete()
    old = {n: r for n, r in ((p.name, A.parse(p.name)) for p in sorted(A.W.iterdir())) if r and r["normal"] and r["G"] is not None}
    old_conf, _ = A.conformer_best(old)
    panel, pstat = load_new(PANEL, ("La", "Dy"))
    series, sstat = load_new(SERIES, ("Nd", "Lu"))
    report = {"status": {**pstat, **sstat}, "results": {}}
    for tag, base in (("conformer_best", old_conf), ("primary_only", old)):
        ok = {**base, **panel, **series}
        res = {M: {l: exchange(ok, l, M) for l in EXP} for M in METALS}
        print(f"\n===== {tag} (old La/Dy states: {tag}; panel and Nd/Lu: primary only) =====")
        summ, cal_all = {}, {}
        for M in METALS:
            r = {l: v for l, v in res[M].items() if v and M in EXP[l]}
            if "ac" not in r:
                print(f"{M}: no acetate exchange yet"); continue
            off = r["ac"]["dG"] - EXP["ac"][M]
            resid = {l: r[l]["dG"] - off - EXP[l][M] for l in r}
            ev = [l for l in r if l != "ac"]
            rms0 = float(np.sqrt(np.mean([resid[l] ** 2 for l in ev]))) if ev else float("nan")
            print(f"\n{M}-La  offset(acetate) {off:+.2f}   S-H0 RMS over {len(ev)} ligands {rms0:.2f} -> {'PASS' if rms0 < 0.8 else 'FAIL'} (< 0.8)")
            for l in r:
                print(f"   {l:6s} calc {r[l]['dG']:+7.2f} exp {EXP[l][M]:+6.2f} resid {resid[l]:+6.2f}  n_disp {r[l]['n_disp']:.1f}  {','.join(r[l]['states'])}")
            cals, pars, sc = {}, {}, {}
            for cs, ligs in CAL.items():
                cals[cs], pars[cs] = calibrate(r, M, EXP, ligs)
                held = [l for l in r if l not in ligs] if cs == "primary" else [l for l in HELD_SECONDARY if l in r]
                for f in ("H2", "H3"):
                    if f not in cals[cs] or not held:
                        continue
                    err = {l: cals[cs][f](l) - EXP[l][M] for l in held}
                    sc[(cs, f)] = {"err": err, "rms": rms(err.values()), "params": pars[cs][f]}
                    pr = "  ".join(f"{k} {v:+.3f}" for k, v in pars[cs][f].items())
                    print(f"   S-{f} [{cs}: cal {'+'.join(ligs)}] ({pr}): RMS over {len(held)} {sc[(cs, f)]['rms']:.2f}"
                          + (f" -> {'PASS' if sc[(cs, f)]['rms'] < 0.8 else 'FAIL'}" if cs == "primary" else " (secondary)")
                          + "   " + " ".join(f"{l} {v:+.2f}" for l, v in err.items()))
            rms2 = sc.get(("primary", "H2"), {}).get("rms", float("nan"))
            rms3 = sc.get(("primary", "H3"), {}).get("rms", float("nan"))
            L = list(r); y = np.array([r[l]["dG"] - EXP[l][M] for l in L])
            for name, Xd in (("constant", np.ones((len(L), 1))), ("+a*q_c^2", np.c_[np.ones(len(L)), [Q_C[l] ** 2 for l in L]]),
                             ("+b*n_disp", np.c_[np.ones(len(L)), [r[l]["n_disp"] for l in L]])):
                if len(L) > Xd.shape[1]:
                    pr = [Xd[i] @ np.linalg.lstsq(Xd[np.arange(len(L)) != i], y[np.arange(len(L)) != i], rcond=None)[0] for i in range(len(L))]
                    print(f"   LOO bias model {name:10s}: RMS {np.sqrt(np.mean((y - np.array(pr)) ** 2)):.2f} (n={len(L)})")
            summ[M] = {"offset": off, "rms_SH0": rms0, "rms_SH2": rms2, "rms_SH3": rms3, "resid": resid,
                       "forms": {f"{cs}/{f}": v for (cs, f), v in sc.items()}}
            cal_all[M] = cals
        if all(M in summ for M in METALS):
            c = [summ[M]["offset"] for M in METALS]
            print(f"\nMono: offsets Nd {c[0]:+.2f}  Dy {c[1]:+.2f}  Lu {c[2]:+.2f}  -> {'monotonic' if c[0] > c[1] > c[2] else 'NOT monotonic'}"
                  f";  per A of contraction: " + "  ".join(f"{M} {summ[M]['offset'] / (R9['La'] - R9[M]):+.1f}" for M in METALS))
            # whole-series reconstruction from calibrated (Nd, Dy, Lu); form chosen by the panel's TTHA call
            call = panel_call(tag); sel = FORM_OF_CALL.get(call)
            print(f"Panel TTHA call ({tag}): {call} -> whole-series form {sel or 'none selected (all reported)'}")
            cols = [LN[1:].index(m) for m in METALS]
            ws = {}
            for cs, ligs in CAL.items():
                for f in ("H0", "H2", "H3"):
                    if f == "H0" and cs != "primary" or not all(f in cal_all[M][cs] for M in METALS):
                        continue
                    errs = {}
                    for l in EXP:
                        if not all(res[M][l] and M in EXP[l] for M in METALS):
                            continue
                        keep = [i for i, x in enumerate(ids) if x != LIG_ID.get(l)]
                        B = np.linalg.lstsq(NX[keep][:, cols], NX[keep], rcond=None)[0]
                        pred = np.array([cal_all[M][cs][f](l) for M in METALS]) @ B
                        errs[l] = rms([pred[LN[1:].index(m)] - EXP[l][m] for m in EXP[l]])
                    excl = ("ac",) if f == "H0" else ligs  # calibration ligands are excluded from the score
                    allrms = rms([v for l, v in errs.items() if l not in excl])
                    key = f"{cs}/{f}"; ws[key] = {"by_ligand": errs, "overall": allrms, "excluded": list(excl)}
                    flag = "  <- panel-selected (pre-registered)" if f == sel and cs == "primary" else ""
                    verdict = "n/a (no held-out ligand)" if allrms != allrms else "PASS" if allrms < 0.8 else "FAIL"
                    print(f"Whole series [{key}] (13 metals) RMS by ligand: " + "  ".join(f"{l} {v:.2f}" for l, v in errs.items())
                          + f"   overall {allrms:.2f} -> {verdict} (< 0.8){flag}")
            summ["whole_series"] = {"panel_call": call, "selected": f"primary/{sel}" if sel else None, "forms": ws}
        report["results"][tag] = {"exchange": res, "summary": summ}
    out = SERIES.parent / "SERIES_RESULT.json"
    out.write_text(json.dumps(report, indent=1, default=float))
    print("\nwrote", out)


if __name__ == "__main__":
    main()
