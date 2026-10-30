"""Decompose the La->Dy exchange at fixed geometry into model-ingredient contributions.

dE_X(combo) = E_X(Dy.L) + E_X(La aq) - E_X(La.L) - E_X(Dy aq) + n E_X(water),  X in
{lc_cpcm (reference Opt energy), gas, sc_cpcm, sc_gas, dlpno}. Contributions:
  CPCM        = lc_cpcm - gas
  4f explicit = sc_cpcm - lc_cpcm   (and sc_gas - gas)
  functional  = dlpno - gas         (CCSD(T) vs PBE0-D4, both gas, large core)
Also: differential hydration of the bare ions, Dy - La, cluster-continuum at each
model's preferred CN (La 9, Dy 8), against experiment (about -67 kcal/mol, Marcus; recalled).
"""
import json, re, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "diagnostics/lady_relaxed_exchange_20261024"))
import analyze as A

T = REPO / "workspaces/lady_ingredient_tests_20261026"
H = A.H


def e(test, name):
    p = T / test / name / "job.out"
    if not p.exists(): return None
    t = p.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in t: return None
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    return float(m[-1]) if m else None


def s2(test, name):
    p = T / test / name / "job.out"
    m = re.findall(r"Expectation value of <S\*\*2>\s*:\s*(\d+\.\d+)", p.read_text(errors="replace")) if p.exists() else []
    return float(m[-1]) if m else None


def main():
    names = [p.name for p in sorted(A.W.iterdir())]
    ok = {n: r for n, r in ((n, A.parse(n)) for n in names) if r and r["normal"] and r["G"] is not None}
    ref, _ = A.conformer_best(ok)
    E = {"lc_cpcm": {n: r["E"] for n, r in ref.items()}}
    for t in ("gas", "sc_cpcm", "sc_gas", "dlpno", "r_ionic", "r_small", "r_large"):
        E[t] = {n: e(t, n) for n in list(ref) + ["water"]}
    E["sc_cpcm"]["water"] = E["lc_cpcm"]["water"]      # no metal: identical model
    E["sc_gas"]["water"] = E["gas"]["water"]
    for t in ("r_ionic", "r_small", "r_large"):   # water has no metal: radius change does not apply
        E[t]["water"] = E["lc_cpcm"]["water"]
    combos = {
        "ac merged best (La mono w8, Dy bi w7 | La aq9, Dy aq8)": ("La_acmono_w8", "Dy_acbi_w7", "La_aq9", "Dy_aq8"),
        "acmono best (La w8, Dy w7 | La aq9, Dy aq8)": ("La_acmono_w8", "Dy_acmono_w7", "La_aq9", "Dy_aq8"),
        "acmono CN9 matched": ("La_acmono_w8", "Dy_acmono_w8", "La_aq9", "Dy_aq9"),
        "acmono CN8 matched": ("La_acmono_w7", "Dy_acmono_w7", "La_aq8", "Dy_aq8"),
        "acbi best (La w7, Dy w7 | La aq9, Dy aq8)": ("La_acbi_w7", "Dy_acbi_w7", "La_aq9", "Dy_aq8"),
        "acbi CN9 matched": ("La_acbi_w7", "Dy_acbi_w7", "La_aq9", "Dy_aq9"),
        "acbi CN8 matched": ("La_acbi_w6", "Dy_acbi_w6", "La_aq8", "Dy_aq8"),
    }
    out = {}
    for label, (laL, dyL, laA, dyA) in combos.items():
        n = A.nw(laL) + A.nw(dyA) - A.nw(dyL) - A.nw(laA)
        row = {}
        for t, Et in E.items():
            v = [Et.get(x) for x in (dyL, laA, laL, dyA, "water")]
            row[t] = None if None in v else H * (v[0] + v[1] - v[2] - v[3] + n * v[4])
        g = lambda a, b: None if row[a] is None or row[b] is None else row[a] - row[b]
        row["CPCM contribution"] = g("lc_cpcm", "gas")
        row["4f explicit (CPCM)"] = g("sc_cpcm", "lc_cpcm")
        row["4f explicit (gas)"] = g("sc_gas", "gas")
        row["CCSD(T) - PBE0 (gas)"] = g("dlpno", "gas")
        for t in ("r_ionic", "r_small", "r_large"):
            row[f"CPCM metal radius {t}"] = g(t, "lc_cpcm")
        out[label] = row
    # differential hydration of the bare ions, Dy - La, at La CN9 / Dy CN8
    hyd = {}
    gw = ok["water"]["G"] * H + A.WATER_CORR
    for core, t_ion in (("large", "gas"), ("small", "sc_gas")):
        eL, eD = e(t_ion, "La_ion"), e(t_ion, "Dy_ion")
        if None in (eL, eD): continue
        if core == "large":
            gL, gD = ref["La_aq9"]["G"] * H, ref["Dy_aq8"]["G"] * H
        else:  # small core: lc thermal corrections on sc_cpcm electronic energies
            if None in (E["sc_cpcm"]["La_aq9"], E["sc_cpcm"]["Dy_aq8"]): continue
            gL = (E["sc_cpcm"]["La_aq9"] + ref["La_aq9"]["G"] - ref["La_aq9"]["E"]) * H
            gD = (E["sc_cpcm"]["Dy_aq8"] + ref["Dy_aq8"]["G"] - ref["Dy_aq8"]["E"]) * H
        hyd[core] = (gD - H * eD) - (gL - H * eL) + gw   # one extra water stays in the liquid for Dy CN8
    res = {"exchange_dE_kcal": out, "ddG_hyd_Dy_minus_La_kcal": hyd, "ddG_hyd_expt_approx": -67,
           "S2_Dy_sc": {n: s2("sc_cpcm", n) for n in ref if n.startswith("Dy")}}
    (T / "RESULT.json").write_text(json.dumps(res, indent=1))
    for k, r in out.items():
        print(f"{k:45s} " + "  ".join(f"{c}={'   n/a' if v is None else f'{v:+6.2f}'}" for c, v in r.items()))
    print("ddG_hyd(Dy-La):", {k: round(v, 2) for k, v in hyd.items()}, "expt about -67")
    print("Dy sc <S^2> (6-et ideal 8.75):", res["S2_Dy_sc"])


if __name__ == "__main__":
    main()
