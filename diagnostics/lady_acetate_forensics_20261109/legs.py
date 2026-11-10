"""Analyses 2-4: leg decomposition + inheritance models (LOO), exchange term split, CPCM radius knob,
DLPNO-CCSD(T) vs PBE0-D4 at fixed geometry, acetate thermo split, FEP legs. Saved outputs only.
Writes workspaces/lady_acetate_forensics_20261109/legs.json.
"""
import json, re, sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "diagnostics/lady_relaxed_exchange_20261024"))
sys.path.insert(0, str(REPO / "diagnostics/lady_panel_20261028"))
import analyze as A  # noqa: E402
import panel_analyze as P  # noqa: E402

OUT = REPO / "workspaces/lady_acetate_forensics_20261109"
T = REPO / "workspaces/lady_ingredient_tests_20261026"
H, WC = A.H, A.WATER_CORR
EXP_HYD, EXP_HYD_H = -66.9, -69.8
EXP = {k: v[0] for k, v in P.EXP.items()}
QC = {k: v[1] for k, v in P.EXP.items()}


def last(pat, t):
    m = re.findall(pat, t); return float(m[-1]) if m else None


def comps(d):
    t = (d / "job.out").read_text(errors="replace")
    E = last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    disp = last(r"Dispersion correction\s+(-?\d+\.\d+)\s*\n", t)
    cp = last(r"CPCM Dielectric\s+:\s+(-?\d+\.\d+)", t)
    G = last(r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)", t)
    Hh = last(r"Total Enthalpy\s+\.\.\.\s+(-?\d+\.\d+)", t)
    return {"E": E, "disp": disp, "cpcm": cp, "rest": E - disp - cp, "thermal": G - E, "G": G, "H": Hh}


def bound(d):
    L = [l.split() for l in (d / "job.xyz").read_text().splitlines()[2:] if l.strip()]
    el = [l[0] for l in L]; X = np.array([[float(v) for v in l[1:4]] for l in L])
    dm = np.linalg.norm(X - X[0], axis=1); nO = nN = nw = 0
    for i in range(1, len(el)):
        if dm[i] >= 2.95 or el[i] not in ("O", "N"): continue
        if el[i] == "N": nN += 1; continue
        nH = sum(1 for j in range(len(el)) if el[j] == "H" and np.linalg.norm(X[j] - X[i]) < 1.25)
        if nH >= 2: nw += 1
        else: nO += 1
    return {"O_lig": nO, "N": nN, "water": nw}


def loo(X, y):
    pr = []
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        c = np.linalg.lstsq(X[m], y[m], rcond=None)[0]; pr.append(X[i] @ c)
    full = np.linalg.lstsq(X, y, rcond=None)[0]
    return float(np.sqrt(np.mean((y - np.array(pr)) ** 2))), [float(v) for v in full], [float(v) for v in y - np.array(pr)]


def main():
    # ---------- reconstruct conformer_best states + panel exactly as panel_analyze, without writing
    old = P.states_ok({p.name: A.parse(p.name) for p in sorted(A.W.iterdir())})
    old_conf, picked = A.conformer_best(old)
    new_R = {p.name: A.parse(p.name, p) for p in sorted(P.PANEL.iterdir())}
    new = P.states_ok(new_R)
    for n, r in list(new.items()):
        lig = n.split("_")[1]
        if lig in P.NO_CN_CHECK: continue
        cn_all = sum(1 for _, x in r.get("donors", []) if x < 3.1)
        if cn_all != P.N_LIG_DONORS[lig] + A.nw(n): new.pop(n)
    ok = {**old_conf, **new}
    def where(n):
        if n in new: return P.PANEL / n
        ch = picked.get(n, {}).get("chosen", "primary")
        return A.W / n if ch == "primary" else A.W.parent / ch.split("/")[0] / n / ch.split("/")[1]
    gw = ok["water"]["G"] * H + WC
    eL = last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", (T / "gas/La_ion/job.out").read_text())
    eD = last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", (T / "gas/Dy_ion/job.out").read_text())
    hyd = (ok["Dy_aq8"]["G"] - eD) * H - (ok["La_aq9"]["G"] - eL) * H + gw
    e_aq = hyd - EXP_HYD
    print(f"aqua leg: calc {hyd:+.2f} exp {EXP_HYD:+.2f} error e_aq = {e_aq:+.2f} (positive = Dy under-hydrated)")
    rows = {}
    wc = comps(A.W / "water")
    for lig in ("ac", "edtax", "nta", "dtpa", "cdta", "tmdta", "mal", "ttha"):
        pool = (lambda n, m: A.pool_match(n, m, "ac")) if lig == "ac" else (lambda n, m, l=lig: n.startswith(f"{m}_{l}_"))
        try:
            _, b = P.exchange(ok, lig, pool)
        except Exception:
            b = None
        if not b:
            print(f"{lig}: unfinished (no complete state set)"); continue
        st = {k: b[k] for k in ("La_L", "Dy_L", "La_aq", "Dy_aq")}
        s = {k: comps(where(v)) for k, v in st.items()}
        dd = {k: H * (s["Dy_L"][k] + s["La_aq"][k] - s["La_L"][k] - s["Dy_aq"][k] + b["n"] * wc[k]) for k in ("E", "disp", "cpcm", "rest", "thermal", "H")}
        bl = {m: bound(where(st[f"{m}_L"])) for m in ("La", "Dy")}
        nd = ((A.nw(st["La_aq"]) - A.nw(st["La_L"])) + (A.nw(st["Dy_aq"]) - A.nw(st["Dy_L"]))) / 2
        err = b["dG"] - EXP[lig]
        rows[lig] = {"states": st, "n": b["n"], "dG": b["dG"], "dH": b["dH"], "exp": EXP[lig], "err": err, "q_c": QC[lig],
                     "n_disp": nd, "n_Ocarb": (bl["La"]["O_lig"] + bl["Dy"]["O_lig"]) / 2, "bound": bl,
                     "complex_leg_calc": b["dG"] + hyd, "complex_leg_exp": EXP[lig] + EXP_HYD, "complex_leg_err": err + e_aq,
                     "terms": dd, "dG_check": dd["E"] + dd["thermal"] + b["n"] * WC}
    print(f"\n{'lig':6s} {'q':>3s} {'calc':>7s} {'exp':>6s} {'err':>6s} | {'cplx_err':>8s} {'aq_err':>6s} | ndisp nOcarb | {'dE':>7s} {'D4':>6s} {'CPCM':>7s} {'rest':>7s} {'G-E':>6s}")
    for l, r in rows.items():
        t = r["terms"]
        print(f"{l:6s} {r['q_c']:+3d} {r['dG']:+7.2f} {r['exp']:+6.2f} {r['err']:+6.2f} | {r['complex_leg_err']:+8.2f} {e_aq:+6.2f} | "
              f"{r['n_disp']:5.1f} {r['n_Ocarb']:6.1f} | {t['E']:+7.2f} {t['disp']:+6.2f} {t['cpcm']:+7.2f} {t['rest']:+7.2f} {t['thermal']:+6.2f}  {r['states']}")
    # ---------- inheritance models
    L = [l for l in rows]
    y = np.array([rows[l]["err"] for l in L]); one = np.ones(len(L))
    nd = np.array([rows[l]["n_disp"] for l in L]); q2 = np.array([9 - rows[l]["q_c"] ** 2 for l in L], float)
    no = np.array([rows[l]["n_Ocarb"] for l in L])
    models = {"M0 const": one[:, None], "M1 a+b*n_disp": np.c_[one, nd], "M2 a+b*(9-q_c^2)": np.c_[one, q2], "M3 a+b*n_Ocarb": np.c_[one, no]}
    inh = {}
    print(f"\ninheritance models on calc-exp (n={len(L)}: {L}); LOO RMS, full-fit coefficients, LOO residuals")
    for k, X in models.items():
        r_, c, res = loo(X, y); inh[k] = {"loo_rms": r_, "coef": c, "loo_resid": dict(zip(L, res))}
        print(f"  {k:18s} LOO {r_:.2f}  coef {np.round(c, 3)}  resid " + " ".join(f"{l}:{v:+.2f}" for l, v in zip(L, res)))
    naq = 8.5
    z1 = -e_aq * nd / naq; z2 = -e_aq * (1 - (9 - q2) / 9)
    for k, z in (("Z1 -e_aq*n_disp/8.5", z1), ("Z2 -e_aq*(1-q_c^2/9)", z2)):
        rms = float(np.sqrt(np.mean((y - z) ** 2)))
        inh[k] = {"rms": rms, "pred": dict(zip(L, map(float, z)))}
        print(f"  {k:22s} (zero-parameter) RMS {rms:.2f}  pred " + " ".join(f"{l}:{v:+.2f}" for l, v in zip(L, z)))
    print("  corr(n_disp, n_Ocarb) = %.2f, corr(n_disp, 9-q^2) = %.2f" % (np.corrcoef(nd, no)[0, 1], np.corrcoef(nd, q2)[0, 1]))
    # ---------- CPCM metal-radius knob (fixed geometry; thermal terms unchanged)
    def e_t(test, n):
        p = T / test / n / "job.out"
        if not p.exists(): return None
        t = p.read_text(errors="replace")
        return last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t) if "ORCA TERMINATED NORMALLY" in t else None
    knob = {}
    st_ac = rows["ac"]["states"]
    for t in ("r_ionic", "r_small", "r_large"):
        dE = {n: (e_t(t, n) - ok[n]["E"]) * H for n in set(st_ac.values()) | {"La_aq9", "Dy_aq8"}}
        h = hyd + dE["Dy_aq8"] - dE["La_aq9"]
        a = rows["ac"]["dG"] + dE[st_ac["Dy_L"]] + dE[st_ac["La_aq"]] - dE[st_ac["La_L"]] - dE[st_ac["Dy_aq"]]
        knob[t] = {"ddG_hyd": h, "ac": a}
        print(f"  radius {t:8s}: ddG_hyd {h:+.2f} (d {h - hyd:+.2f})  acetate {a:+.2f} (d {a - rows['ac']['dG']:+.2f})")
    # ---------- DLPNO-CCSD(T) vs PBE0-D4 gas, fixed geometry, n=0 combinations
    combos = {"ac selected (La mono w8, Dy bi w7 | aq9, aq8)": ("La_acmono_w8", "Dy_acbi_w7", "La_aq9", "Dy_aq8"),
              "mono, La w8 / Dy w7 | aq9 / aq8": ("La_acmono_w8", "Dy_acmono_w7", "La_aq9", "Dy_aq8"),
              "mono CN9 matched | aq9 / aq9": ("La_acmono_w8", "Dy_acmono_w8", "La_aq9", "Dy_aq9"),
              "bi CN9 matched | aq9 / aq9": ("La_acbi_w7", "Dy_acbi_w7", "La_aq9", "Dy_aq9"),
              "water->acetate swap at CN9, metal aq9 (Dy-La on [mono w8 - aq9])": ("La_acmono_w8", "Dy_acmono_w8", "La_aq9", "Dy_aq9")}
    dl = {}
    for k, (laL, dyL, laA, dyA) in combos.items():
        v = {}
        for t in ("dlpno", "gas", "sc_gas"):
            e = [e_t(t, x) for x in (dyL, laA, laL, dyA)]
            v[t] = None if None in e else H * (e[0] + e[1] - e[2] - e[3])
        v["dlpno_minus_pbe0d4"] = None if None in (v["dlpno"], v["gas"]) else v["dlpno"] - v["gas"]
        dl[k] = v
        print(f"  {k:62s} PBE0-D4 gas {v['gas']:+7.2f}  DLPNO {v['dlpno']:+7.2f}  diff {v['dlpno_minus_pbe0d4']:+6.2f}"
              + (f"  (sc_gas {v['sc_gas']:+.2f})" if v['sc_gas'] is not None else ""))
    # D4 share in gas exchange for the selected acetate combo
    # ---------- per-state DLPNO-PBE0 binding differences (Dy-La) for aqua/complex pairs at the same label
    pair = {}
    for lab in ("aq8", "aq9", "acmono_w7", "acmono_w8", "acbi_w6", "acbi_w7"):
        a, b = e_t("dlpno", f"Dy_{lab}"), e_t("dlpno", f"La_{lab}"); c, d = e_t("gas", f"Dy_{lab}"), e_t("gas", f"La_{lab}")
        pair[lab] = None if None in (a, b, c, d) else H * ((a - b) - (c - d))
    print("  (Dy-La) DLPNO-minus-PBE0 per state label (absolute offset cancels in exchange):",
          {k: (None if v is None else round(v, 2)) for k, v in pair.items()})
    # ---------- thermo split for acetate
    ac = rows["ac"]; dH = ac["dH"]; TdS = dH - ac["dG"]
    g = {"exp_2M": {"ddG": -0.17, "ddH": 0.74, "TddS": 0.89}, "model": {"ddG": ac["dG"], "ddH": dH, "TddS": TdS}}
    g["error"] = {k: g["model"][k] - g["exp_2M"][k] for k in ("ddG", "ddH", "TddS")}
    hydH = (ok["Dy_aq8"]["H"] - eD) * H - (ok["La_aq9"]["H"] - eL) * H + ok["water"]["H"] * H
    g["aqua"] = {"ddH_hyd_calc_approx": hydH, "ddH_hyd_exp": EXP_HYD_H, "TddS_hyd_calc_approx": hydH - hyd, "TddS_hyd_exp": -2.74}
    print(f"\nacetate thermo (model n=0): ddG {ac['dG']:+.2f} ddH {dH:+.2f} TddS {TdS:+.2f} | exp 2M: -0.17 +0.74 +0.89 | "
          f"error ddH {g['error']['ddH']:+.2f} TddS {g['error']['TddS']:+.2f} ddG {g['error']['ddG']:+.2f}")
    print(f"  aqua leg: ddH_hyd calc ~{hydH:+.2f} (exp -69.8), TddS_hyd calc ~{hydH - hyd:+.2f} (exp -2.74)")
    # ---------- FEP legs
    F = json.loads((REPO / "workspaces/lady_explicit_solvent_20261026/FEP_RESULT.json").read_text())
    fep = {}
    for leg in ("aq", "ac", "edta"):
        r = F[leg]["per_T"][0]
        fep[leg] = {"with_tail": r["dG_kcal"], "tail": r["tail_kcal"], "no_tail": r["dG_kcal"] - r["tail_kcal"], "structure": r["structure"],
                    "n_uncorr_min": r["n_uncorrelated_min"]}
    for conv in ("with_tail", "no_tail"):
        aqe = fep["aq"][conv] - EXP_HYD
        for lig, e in (("ac", EXP["ac"]), ("edta", EXP["edtax"])):
            ce = fep[lig][conv] - (e + EXP_HYD)
            fep[f"{lig}_{conv}"] = {"aq_err": aqe, "complex_err": ce, "exch_err": ce - aqe}
            print(f"  FEP {lig:4s} {conv:9s}: aqua-leg err {aqe:+.2f}  complex-leg err {ce:+.2f}  exchange err {ce - aqe:+.2f}")
    s = fep["ac"]["structure"]; sa = fep["aq"]["structure"]
    for m in ("La", "Dy"):
        p = s[m]["ligand_donors_pop"]; dent = sum(int(k) * v for k, v in p.items())
        print(f"  FEP {m}: aqua CN {sa[m]['CN_mean']:.2f}; acetate CN {s[m]['CN_mean']:.2f}, denticity {dent:.2f} {p}, "
              f"waters displaced {sa[m]['CN_mean'] - (s[m]['CN_mean'] - dent):.2f}")
    out = {"ddG_hyd_calc": hyd, "e_aq": e_aq, "ligands": rows, "inheritance": inh, "radius_knob": knob, "dlpno": dl,
           "dlpno_pair_offsets": pair, "thermo": g, "fep": fep}
    (OUT / "legs.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
