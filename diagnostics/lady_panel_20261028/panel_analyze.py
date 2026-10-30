"""Constant-bias panel: La->Dy exchange ΔG with the careful cluster-continuum protocol, against NIST 46.

Reaction (negative = Dy favoured):  La.L(H2O)b + Dy(H2O)c -> Dy.L(H2O)a + La(H2O)d + n H2O
Aqua ions, water, acetate and EDTA come from lady_relaxed_exchange_20261024 (conformer-best and primary-only);
malonate, NTA, DTPA, CDTA and TMDTA come from lady_panel_20261028/runs (primary only, no conformer search).
State choice per metal: the lowest-G hydration state, water-balanced (same rule as analyze.py "best").

Hypotheses tested (stated before the panel ran, 2026-10-28):
  H0 constant bias: calc - exp is the same for every ligand. Pass: RMS residual < 0.8 kcal/mol after
     subtracting the acetate-calibrated offset, and CDTA < TMDTA < EDTA in calculated order.
  H1 Born-like bias: residual tracks complex charge (offset + a*q_c^2).
  H2 per-water bias (added 2026-10-28 18:30): offset + b*n_disp.
Interim, 2026-10-29 (acetate, EDTA, NTA, DTPA finished; malonate, CDTA, TMDTA running): H0 already fails, since
EDTA, NTA and DTPA alone give RMS >= 0.86 over the six test ligands. Two post hoc two-parameter forms fit those four:
  H2  bias = -1.642 - 0.2282*n_disp                  (in-sample RMS 0.41)
  H3  calc = 1.2792*exp - 1.956  (selectivity stretched; in-sample RMS 0.48)
FROZEN (below) applies both, plus the acetate-only H0, unchanged to the held-out ligands. Malonate, CDTA and TMDTA
cannot separate H2 from H3 (predictions within 0.5). TTHA (NIST 678, added 2026-10-29) can: it displaces every water
(n_disp 8.5) but is nearly unselective, so H2 predicts calc -4.26 and H3 -2.83 (H0 -2.56). Call: per-water form if TTHA
calc <= -3.90, stretch/constant form if >= -3.20, inconclusive in between (conformer_best tag).
"""
import itertools, json, sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "diagnostics/lady_relaxed_exchange_20261024"))
import analyze as A  # noqa: E402

PANEL = REPO / "workspaces/lady_panel_20261028/runs"
H, WC = A.H, A.WATER_CORR
EXP = {"ac": (-0.068, +2), "edtax": (-4.010, -1), "mal": (-0.655, +1), "nta": (-1.582, 0),
       "dtpa": (-4.556, -2), "cdta": (-5.674, -1), "tmdta": (-4.869, -1),
       "ttha": (-0.682, -3)}  # (NIST ddG Dy-La, complex charge)
N_LIG_DONORS = {"mal": 2, "nta": 4, "dtpa": 8, "cdta": 6, "tmdta": 6, "ttha": 9}
NO_CN_CHECK = {"ttha"}  # no water to lose; 9 vs 10 bound donors are isomers of the same q = 0 state
# Frozen 2026-10-29 from acetate, EDTA, NTA, DTPA (conformer_best); predicted calc from exp and n_disp
FROZEN = {"H0 constant": lambda e, n: e - 1.879, "H2 per-water": lambda e, n: e - 1.642 - 0.2282 * n,
          "H3 stretch": lambda e, n: 1.2792 * e - 1.956}
SLOPE = {"H0 constant": 1.0, "H2 per-water": 1.0, "H3 stretch": 1.2792}  # d(calc)/d(exp): error in exp units = resid/slope


def states_ok(R):
    return {n: r for n, r in R.items() if r and r["normal"] and r["G"] is not None}


def exchange(ok, lig, pool):
    gw = ok["water"]
    aq = {m: [n for n in ok if n.startswith(f"{m}_aq")] for m in ("La", "Dy")}
    cplx = {m: [n for n in ok if pool(n, m)] for m in ("La", "Dy")}
    rows = []
    for laL, dyaq, dyL, laaq in itertools.product(cplx["La"], aq["Dy"], cplx["Dy"], aq["La"]):
        n = A.nw(laL) + A.nw(dyaq) - A.nw(dyL) - A.nw(laaq)
        d = lambda k, c: H * (ok[dyL][k] + ok[laaq][k] - ok[laL][k] - ok[dyaq][k] + n * gw[k]) + (n * WC if c else 0)
        rows.append({"La_L": laL, "Dy_L": dyL, "La_aq": laaq, "Dy_aq": dyaq, "n": n, "dG": d("G", True), "dH": d("H", False)})
    sel = {}
    for m in ("La", "Dy"):
        for fam, p in (("L", cplx[m]), ("aq", aq[m])):
            if not p: return rows, None
            ref = max(A.nw(x) for x in p)
            sel[(m, fam)] = min(p, key=lambda x: ok[x]["G"] * H + (ref - A.nw(x)) * (gw["G"] * H + WC))
    best = next(r for r in rows if (r["La_L"], r["Dy_L"], r["La_aq"], r["Dy_aq"]) ==
                (sel[("La", "L")], sel[("Dy", "L")], sel[("La", "aq")], sel[("Dy", "aq")]))
    return rows, best


def main():
    old_names = [p.name for p in sorted(A.W.iterdir())]
    old = states_ok({n: A.parse(n) for n in old_names})
    old_conf, _ = A.conformer_best(old)
    new_R = {p.name: A.parse(p.name, p) for p in sorted(PANEL.iterdir())}
    new = states_ok(new_R)
    report = {"panel_status": {}, "results": {}}
    for n, r in new_R.items():
        lig, nw = n.split("_")[1], A.nw(n)
        if r is None: report["panel_status"][n] = "no output"; continue
        dn = r.get("donors", [])
        cn_all = sum(1 for _, x in dn if x < 3.1)
        report["panel_status"][n] = {"normal": r["normal"], "opt_converged": r["opt_converged"], "imag": r["imag"],
                                     "CN_3.1A": cn_all, "CN_expected": N_LIG_DONORS[lig] + nw, "donors": dn}
    for n, s in report["panel_status"].items():  # a state whose coordination changed in DFT is not the intended state
        if isinstance(s, dict) and s["CN_3.1A"] != s["CN_expected"] and n in new and n.split("_")[1] not in NO_CN_CHECK:
            print(f"EXCLUDED {n}: CN {s['CN_3.1A']} != expected {s['CN_expected']}"); new.pop(n)
    for tag, base in (("conformer_best", old_conf), ("primary_only", old)):
        ok = {**base, **new}
        res = {}
        for lig in EXP:
            pool = (lambda n, m: A.pool_match(n, m, "ac")) if lig == "ac" else (lambda n, m, l=lig: n.startswith(f"{m}_{l}_"))
            rows, best = exchange(ok, lig, pool)
            if best: res[lig] = {"best": best, "all": rows}
        report["results"][tag] = res
        if "ac" not in res: continue
        off = res["ac"]["best"]["dG"] - EXP["ac"][0]
        print(f"\n== {tag} (acetate/EDTA states: {tag}; panel: primary only) ==  offset from acetate = {off:+.2f}")
        print(f"{'ligand':7s} {'q_c':>4s} {'calc':>7s} {'exp':>7s} {'calc-exp':>9s} {'resid(after offset)':>20s}  states")
        lig_eval = [l for l in res if l != "ac"]
        resid = {}
        for l in res:
            b = res[l]["best"]; e, q = EXP[l]
            resid[l] = b["dG"] - off - e
            print(f"{l:7s} {q:+4d} {b['dG']:+7.2f} {e:+7.2f} {b['dG']-e:+9.2f} {resid[l]:+20.2f}  {b['La_L']},{b['Dy_L']},{b['La_aq']},{b['Dy_aq']}")
        if lig_eval:
            rms = float(np.sqrt(np.mean([resid[l] ** 2 for l in lig_eval])))
            order = [l for l in ("cdta", "tmdta", "edtax") if l in res]
            calc_order = sorted(order, key=lambda l: res[l]["best"]["dG"])
            print(f"H0: RMS residual after acetate offset = {rms:.2f} kcal/mol over {lig_eval}  (pass < 0.8)")
            print(f"    order calc {calc_order} vs exp ['cdta','tmdta','edtax']  -> {'PASS' if calc_order == ['cdta', 'tmdta', 'edtax'] else 'FAIL'}")
            L = list(res); y = np.array([res[l]["best"]["dG"] - EXP[l][0] for l in L]); q2 = np.array([EXP[l][1] ** 2 for l in L], float)
            # H2 (added 2026-10-28 18:30, before any panel job finished): bias grows with first-shell waters displaced
            # (EXPERIMENT_THERMO.md: Dy desolvation too cheap per displaced water at the dH level). n_disp = mean over the
            # two metals of (aqua waters - complex waters) in the chosen states.
            b_ = {l: res[l]["best"] for l in L}
            nd = np.array([((A.nw(b_[l]["La_aq"]) - A.nw(b_[l]["La_L"])) + (A.nw(b_[l]["Dy_aq"]) - A.nw(b_[l]["Dy_L"]))) / 2 for l in L])
            print("    waters displaced per metal:", " ".join(f"{l}:{x:.1f}" for l, x in zip(L, nd)))
            ex = np.array([EXP[l][0] for l in L])
            for name, Xd in (("constant", np.ones((len(L), 1))), ("constant + a*q_c^2", np.c_[np.ones(len(L)), q2]),
                             ("constant + b*n_disp", np.c_[np.ones(len(L)), nd]), ("constant + c*exp", np.c_[np.ones(len(L)), ex])):
                pr = []
                for i in range(len(L)):
                    m = np.arange(len(L)) != i
                    c = np.linalg.lstsq(Xd[m], y[m], rcond=None)[0]; pr.append(Xd[i] @ c)
                print(f"    LOO bias model {name:20s}: RMS {np.sqrt(np.mean((y - np.array(pr)) ** 2)):.2f}  (n={len(L)}, all ligands incl. acetate)")
            # frozen 2026-10-29 forms on the held-out ligands; error in experimental units (what a user of the calibration sees)
            held = [l for l in L if l in ("mal", "cdta", "tmdta", "ttha")]
            frozen = {}
            if held:
                sh = 0.141 if tag == "primary_only" else 0.0  # primary-only fit: same slopes, intercepts +0.141
                print(f"    FROZEN forms (fit 2026-10-29 on ac, edtax, nta, dtpa) on held-out {held}: calibrated error, kcal/mol")
                for f, fn in FROZEN.items():
                    err = {l: (res[l]["best"]["dG"] - fn(EXP[l][0], nd[L.index(l)]) - sh) / SLOPE[f] for l in held}
                    r_ = float(np.sqrt(np.mean(np.square(list(err.values())))))
                    frozen[f] = {"err": err, "rms": r_}
                    print(f"      {f:13s} " + "  ".join(f"{l} {v:+.2f}" for l, v in err.items()) + f"   RMS {r_:.2f} -> {'PASS' if r_ < 0.8 else 'FAIL'} (< 0.8)")
                if "ttha" in res:
                    t = res["ttha"]["best"]["dG"]
                    call = "per-water (H2)" if t <= -3.90 + sh else "stretch/constant (H3/H0)" if t >= -3.20 + sh else "inconclusive"
                    print(f"    TTHA calc {t:+.2f}  (H2 {-4.26 + sh:.2f}, H3 {-2.83 + sh:.2f}, H0 {-2.56 + sh:.2f})  -> {call}")
                    frozen["ttha_call"] = call
            report["results"][tag + "_summary"] = {"offset": off, "residuals": resid, "rms": rms, "calc_order": calc_order,
                                                   "frozen": frozen}
    out = PANEL.parent / "PANEL_RESULT.json"
    out.write_text(json.dumps(report, indent=1, default=float))
    print("\nwrote", out)


if __name__ == "__main__":
    main()
