"""La->Dy exchange free energies from relaxed PBE0-D4/CPCM Opt+Freq runs.

Reaction for ligand L (L = acetate mono/bi, EDTA):
  La.L(H2O)b + Dy(H2O)c  ->  Dy.L(H2O)a + La(H2O)d + n H2O,   n = b + c - a - d
dG = sum G_orca + n*(RT ln(24.46) + RT ln(55.34))   (1 atm->1 M for each
released water, then water at 55.34 M; the 1 M conversion of the four solutes cancels).
Negative = Dy favoured. Reports every CN combination plus the
lowest-G hydration state per metal/ligand ("best").
EDTA: "edta" = xtb-wrapped starts (came out partly unwrapped, not hexadentate);
"edtax" = starts from COD 7236231 Na[Dy(EDTA)(H2O)3].5H2O (La by metal swap).
"""
import itertools, json, re, sys
from pathlib import Path
import numpy as np

H = 627.509474
RT = 0.0019872 * 298.15
WATER_CORR = RT * (np.log(24.46) + np.log(55.34))
W = Path(__file__).resolve().parents[2] / "workspaces/lady_relaxed_exchange_20261024/runs"


def parse(name, d=None):
    d = d or W / name
    out = d / "job.out"
    if not out.exists():
        return None
    t = out.read_text(errors="replace")
    r = {"normal": "ORCA TERMINATED NORMALLY" in t,
         "opt_converged": "OPTIMIZATION RUN DONE" in t or "THE OPTIMIZATION HAS CONVERGED" in t}
    for key, pat in (("E", r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"),
                     ("G", r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)"),
                     ("H", r"Total Enthalpy\s+\.\.\.\s+(-?\d+\.\d+)")):
        m = re.findall(pat, t)
        r[key] = float(m[-1]) if m else None
    freqs = re.findall(r"^\s*\d+:\s+(-?\d+\.\d+) cm\*\*-1", t, re.M)
    r["imag"] = [float(f) for f in freqs if float(f) < 0]
    xyz = d / "job.xyz"
    if xyz.exists() and name != "water":
        L = [l.split() for l in xyz.read_text().splitlines()[2:] if l.strip()]
        X = np.array([[float(v) for v in l[1:4]] for l in L])
        d = [(L[i][0], float(np.linalg.norm(X[i] - X[0]))) for i in range(1, len(L)) if L[i][0] in ("O", "N")]
        r["donors"] = sorted([(e, round(x, 3)) for e, x in d if x < 3.1], key=lambda z: z[1])
    return r


def cn(r):
    return sum(1 for _, x in r.get("donors", []) if x < 2.95)


def conformer_best(ok):
    """Replace each system by its lowest-G DFT structure among the primary run and
    the CREST-conformer re-optimizations (conf/, conf_v2/), keeping only structures
    with the primary's coordination number and no imaginary mode below -30 cm-1."""
    best, picked = {}, {}
    for n, r in ok.items():
        cands = [("primary", r)]
        for c in sorted(W.parent.glob(f"conf*/{n}/dft_c*")):
            rc = parse(n, c)
            if rc and rc["normal"] and rc["G"] is not None:
                cands.append((f"{c.parent.parent.name}/{c.name}", rc))
        good = [(k, x) for k, x in cands if n == "water" or (cn(x) == cn(r) and min(x["imag"], default=0) > -30)]
        k, x = min(good or cands, key=lambda kx: kx[1]["G"])
        best[n] = x
        picked[n] = {"chosen": k, "n_candidates": len(cands), "n_valid": len(good),
                     "dG_vs_primary": H * (x["G"] - r["G"])}
    return best, picked


def pool_match(n, m, lig):
    """'ac' = all acetate binding modes (mono and bi are one species; the
    lowest water-balanced state per metal is used)."""
    if lig == "ac":
        return n.startswith(f"{m}_acmono_") or n.startswith(f"{m}_acbi_")
    return n.startswith(f"{m}_{lig}_")


def nw(name):
    m = re.search(r"(?:aq|_w)(\d+)$", name)
    return int(m.group(1))


def main():
    names = [p.name for p in sorted(W.iterdir())]
    R = {n: parse(n) for n in names}
    ok = {n: r for n, r in R.items() if r and r["normal"] and r["G"] is not None}
    conf = "--conformers" in sys.argv
    if conf:
        ok, picked = conformer_best(ok)
    report = {"status": {n: (None if r is None else {k: r[k] for k in ("normal", "opt_converged", "imag")}) for n, r in R.items()},
              "donors": {n: r.get("donors") for n, r in ok.items()}}
    if "water" not in ok:
        print(json.dumps(report, indent=1)); return
    gw = ok["water"]
    aq = {m: [n for n in ok if n.startswith(f"{m}_aq")] for m in ("La", "Dy")}
    rows = []
    for lig in ("acmono", "acbi", "ac", "edta", "edtax"):
        cplx = {m: [n for n in ok if pool_match(n, m, lig)] for m in ("La", "Dy")}
        for laL, dyaq, dyL, laaq in itertools.product(cplx["La"], aq["Dy"], cplx["Dy"], aq["La"]):
            n = nw(laL) + nw(dyaq) - nw(dyL) - nw(laaq)
            def d(key, corr):
                return H * (ok[dyL][key] + ok[laaq][key] - ok[laL][key] - ok[dyaq][key] + n * gw[key]) + (n * WATER_CORR if corr else 0)
            rows.append({"ligand": lig, "La_L": laL, "Dy_L": dyL, "La_aq": laaq, "Dy_aq": dyaq, "n_water_released": n,
                         "dG": d("G", True), "dH": d("H", False), "dE": d("E", False)})
    best = {}
    for lig in ("acmono", "acbi", "ac", "edta", "edtax"):
        sel = {}
        for m in ("La", "Dy"):
            for fam, pool in (("L", [n for n in ok if pool_match(n, m, lig)]), ("aq", aq[m])):
                if not pool: continue
                # lowest G per hydration family, water-balanced to the largest hydrate
                ref = max(nw(x) for x in pool)
                sel[(m, fam)] = min(pool, key=lambda x: ok[x]["G"] * H + (ref - nw(x)) * (gw["G"] * H + WATER_CORR))
        if len(sel) == 4:
            r = next(r for r in rows if r["ligand"] == lig and r["La_L"] == sel[("La", "L")] and r["Dy_L"] == sel[("Dy", "L")]
                     and r["La_aq"] == sel[("La", "aq")] and r["Dy_aq"] == sel[("Dy", "aq")])
            best[lig] = r
    report["best_state_exchange_kcal"] = best
    report["all_combinations_kcal"] = rows
    if conf:
        report["conformer_selection"] = picked
    out = W.parent / ("EXCHANGE_RESULT_CONFORMERS.json" if conf else "EXCHANGE_RESULT.json")
    out.write_text(json.dumps(report, indent=1, default=float))
    for lig, r in best.items():
        print(f"{lig:7s} dG={r['dG']:+7.2f} dH={r['dH']:+7.2f} dE={r['dE']:+7.2f}  states {r['La_L']},{r['Dy_L']},{r['La_aq']},{r['Dy_aq']} n={r['n_water_released']}")
    print(json.dumps(report["status"], indent=0))


if __name__ == "__main__":
    main()
