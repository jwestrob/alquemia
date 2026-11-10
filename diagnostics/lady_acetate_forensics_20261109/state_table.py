"""Step 0 positive control + Analysis 1 (acetate/aqua state table and exchange under state treatments).

Saved outputs only. Imports analyze.py helpers; never calls analyze.main() (it rewrites saved JSON).
Writes workspaces/lady_acetate_forensics_20261109/{state_table.tsv, state_treatments.json}.
"""
import itertools, json, math, re, sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "diagnostics/lady_relaxed_exchange_20261024"))
import analyze as A  # noqa: E402

OUT = REPO / "workspaces/lady_acetate_forensics_20261109"
H, WC, RT = A.H, A.WATER_CORR, A.RT
W = A.W


def last(pat, t):
    m = re.findall(pat, t)
    return float(m[-1]) if m else None


def path_of(name, chosen):
    if chosen == "primary":
        return W / name
    a, b = chosen.split("/")
    return W.parent / a / name / b


def extra(d):
    """CPCM dielectric, D4, ligand O count, CN from a job directory."""
    t = (d / "job.out").read_text(errors="replace")
    r = {"cpcm": last(r"CPCM Dielectric\s+:\s+(-?\d+\.\d+)", t),
         "disp": last(r"Dispersion correction\s+(-?\d+\.\d+)\s*\n", t),
         "wall_s": None}
    m = re.search(r"TOTAL RUN TIME:\s+(\d+) days (\d+) hours (\d+) minutes (\d+) seconds", t)
    if m:
        dd, hh, mm, ss = map(int, m.groups()); r["wall_s"] = ((dd * 24 + hh) * 60 + mm) * 60 + ss
    m = re.search(r"%pal\s+nprocs\s+(\d+)", t)
    r["nprocs"] = int(m.group(1)) if m else None
    r["n_opt_cycles"] = len(re.findall(r"GEOMETRY OPTIMIZATION CYCLE", t))
    xyz = d / "job.xyz"
    if xyz.exists():
        L = [l.split() for l in xyz.read_text().splitlines()[2:] if l.strip()]
        el = [l[0] for l in L]; X = np.array([[float(v) for v in l[1:4]] for l in L])
        r["natoms"] = len(el)
        dm = np.linalg.norm(X - X[0], axis=1)
        nO_lig = nO_w = 0
        for i in range(1, len(el)):
            if el[i] != "O" or dm[i] >= 2.95: continue
            nH = sum(1 for j in range(len(el)) if el[j] == "H" and np.linalg.norm(X[j] - X[i]) < 1.25)
            if nH >= 2: nO_w += 1
            else: nO_lig += 1
        r["bound_water_O"] = nO_w; r["bound_lig_O"] = nO_lig
        r["far_water_lt3p1"] = sorted(round(float(x), 3) for i, x in enumerate(dm) if el[i] == "O" and 2.95 <= x < 3.4)
    return r


def label_cn(name):
    nw = A.nw(name)
    if "_aq" in name: return nw
    return nw + (2 if "acbi" in name else 1)


def main():
    names = [p.name for p in sorted(W.iterdir())]
    R = {n: A.parse(n) for n in names}
    ok = {n: r for n, r in R.items() if r and r["normal"] and r["G"] is not None}
    best, picked = A.conformer_best(ok)
    gw = best["water"]["G"] * H + WC
    # ---------------- Step 0
    sel = {}
    pool = {("La", "L"): [n for n in best if A.pool_match(n, "La", "ac")], ("Dy", "L"): [n for n in best if A.pool_match(n, "Dy", "ac")],
            ("La", "aq"): [n for n in best if n.startswith("La_aq")], ("Dy", "aq"): [n for n in best if n.startswith("Dy_aq")]}
    for k, p in pool.items():
        ref = max(A.nw(x) for x in p)
        sel[k] = min(p, key=lambda x: best[x]["G"] * H + (ref - A.nw(x)) * gw)
    def exch(st, key="G", src=best, wcorr=True):
        laL, dyL, laA, dyA = st
        n = A.nw(laL) + A.nw(dyA) - A.nw(dyL) - A.nw(laA)
        v = H * (src[dyL][key] + src[laA][key] - src[laL][key] - src[dyA][key] + n * src["water"][key])
        return v + (n * WC if (wcorr and key == "G") else 0.0), n
    st0 = (sel[("La", "L")], sel[("Dy", "L")], sel[("La", "aq")], sel[("Dy", "aq")])
    ac_dG, n0 = exch(st0)
    # aqua: collect.py formula (gas bare-ion E, large core)
    T = REPO / "workspaces/lady_ingredient_tests_20261026"
    def e_sp(test, name):
        t = (T / test / name / "job.out").read_text(errors="replace")
        return last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    eL, eD = e_sp("gas", "La_ion"), e_sp("gas", "Dy_ion")
    hyd = (best["Dy_aq8"]["G"] * H - H * eD) - (best["La_aq9"]["G"] * H - H * eL) + gw
    hyd_prim = (ok["Dy_aq8"]["G"] * H - H * eD) - (ok["La_aq9"]["G"] * H - H * eL) + gw
    # hydration enthalpy analogue (water released: CPCM H only, no standard-state term) -- indicative only
    hydH = (best["Dy_aq8"]["H"] - eD) * H - (best["La_aq9"]["H"] - eL) * H + best["water"]["H"] * H
    print(f"STEP0 acetate conformer_best lowest-G: dG = {ac_dG:+.3f} (n={n0}) states {st0}  [target -1.95]")
    print(f"STEP0 ddG_hyd(Dy-La) conformer_best = {hyd:+.3f}  primary-only = {hyd_prim:+.3f}  [saved RESULT.json -59.21; brief -59.4]")
    print(f"      ddH_hyd analogue = {hydH:+.2f} (exp -69.8; water enthalpy reference approximate)")
    # ---------------- state table
    rows = []
    for n in sorted(best):
        if n == "water" or not (n.startswith("La_a") or n.startswith("Dy_a")) or "edta" in n: continue
        ch = picked[n]["chosen"]; x = extra(path_of(n, ch)); b = best[n]
        rows.append({"state": n, "chosen": ch, "n_cand": picked[n]["n_candidates"], "n_valid": picked[n]["n_valid"],
                     "dG_conf_vs_primary": round(picked[n]["dG_vs_primary"], 3),
                     "E_Eh": b["E"], "H_Eh": b["H"], "G_Eh": b["G"], "cpcm_Eh": x["cpcm"], "disp_Eh": x["disp"],
                     "Hth_kcal": round(H * (b["H"] - b["E"]), 3), "minusTS_kcal": round(H * (b["G"] - b["H"]), 3),
                     "G_wbal_rel_kcal": None, "CN_lt2p95": A.cn(b), "CN_label": label_cn(n),
                     "bound_lig_O": x.get("bound_lig_O"), "bound_water_O": x.get("bound_water_O"),
                     "O_2p95_3p4": x.get("far_water_lt3p1"), "imag": b["imag"],
                     "selected": n in st0, "wall_s": x["wall_s"], "nprocs": x["nprocs"], "natoms": x.get("natoms"),
                     "opt_cycles": x["n_opt_cycles"]})
    # water-balanced relative G within (metal, family): reference 9 waters (aq) / 8 waters (complex)
    for r in rows:
        ref = 9 if "_aq" in r["state"] else 8
        r["G_wbal_kcal"] = r["G_Eh"] * H + (ref - A.nw(r["state"])) * gw
    for m in ("La", "Dy"):
        for fam in ("aq", "ac"):
            grp = [r for r in rows if r["state"].startswith(f"{m}_{fam}")]
            g0 = min(r["G_wbal_kcal"] for r in grp)
            for r in grp: r["G_wbal_rel_kcal"] = round(r["G_wbal_kcal"] - g0, 3)
    OUT.mkdir(parents=True, exist_ok=True)
    cols = ["state", "selected", "chosen", "n_cand", "n_valid", "dG_conf_vs_primary", "G_wbal_rel_kcal", "CN_label", "CN_lt2p95",
            "bound_lig_O", "bound_water_O", "O_2p95_3p4", "imag", "Hth_kcal", "minusTS_kcal", "E_Eh", "cpcm_Eh", "disp_Eh",
            "H_Eh", "G_Eh", "natoms", "nprocs", "wall_s", "opt_cycles"]
    with open(OUT / "state_table.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows: f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("\nSTATE TABLE (G_wbal_rel: within metal+family, waters balanced with G_water+WATER_CORR)")
    for r in rows:
        print(f"  {r['state']:14s} sel={int(r['selected'])} {r['chosen']:16s} cand={r['n_cand']:2d}/{r['n_valid']:2d} "
              f"Grel={r['G_wbal_rel_kcal']:+6.2f} CN={r['CN_lt2p95']}/{r['CN_label']} ligO={r['bound_lig_O']} wO={r['bound_water_O']} "
              f"near={r['O_2p95_3p4']} imag={r['imag']} -TS={r['minusTS_kcal']:+.2f} Hth={r['Hth_kcal']:+.2f}")
    # ---------------- treatments
    st_names = {m: {fam: [r["state"] for r in rows if r["state"].startswith(f"{m}_{fam}")] for fam in ("aq", "ac")} for m in ("La", "Dy")}
    Gw = {r["state"]: r["G_wbal_kcal"] for r in rows}
    def boltz(names):
        g = np.array([Gw[x] for x in names]); g0 = g.min()
        return g0 - RT * math.log(np.exp(-(g - g0) / RT).sum())
    def pops(names):
        g = np.array([Gw[x] for x in names]); p = np.exp(-(g - g.min()) / RT); p /= p.sum()
        return {x: round(float(v), 3) for x, v in zip(names, p)}
    res = {"step0": {"ac_dG": ac_dG, "n": n0, "states": st0, "ddG_hyd": hyd, "ddG_hyd_primary": hyd_prim, "ddH_hyd_analogue": hydH}}
    # (i) lowest
    lo = lambda names: min(Gw[x] for x in names)
    def ex_from(fn, laL, dyL, laA, dyA):
        return (fn(dyL) - fn(laL)) - (fn(dyA) - fn(laA))
    res["i_lowest"] = ex_from(lo, st_names["La"]["ac"], st_names["Dy"]["ac"], st_names["La"]["aq"], st_names["Dy"]["aq"])
    res["ii_boltzmann"] = ex_from(boltz, st_names["La"]["ac"], st_names["Dy"]["ac"], st_names["La"]["aq"], st_names["Dy"]["aq"])
    res["ii_populations"] = {m: {fam: pops(st_names[m][fam]) for fam in ("aq", "ac")} for m in ("La", "Dy")}
    lab = lambda s: s.split("_", 1)[1]
    common = sorted(set(map(lab, st_names["La"]["ac"])) & set(map(lab, st_names["Dy"]["ac"])))
    cm = {m: [f"{m}_{c}" for c in common] for m in ("La", "Dy")}
    res["iii_matched_common_labels"] = common
    res["iii_matched_lowest"] = ex_from(lo, cm["La"], cm["Dy"], st_names["La"]["aq"], st_names["Dy"]["aq"])
    res["iii_matched_boltzmann"] = ex_from(boltz, cm["La"], cm["Dy"], st_names["La"]["aq"], st_names["Dy"]["aq"])
    # mono-only / bi-only pools
    for mode in ("acmono", "acbi"):
        pl = {m: [x for x in cm[m] if mode in x] for m in ("La", "Dy")}
        res[f"pool_{mode}_lowest"] = ex_from(lo, pl["La"], pl["Dy"], st_names["La"]["aq"], st_names["Dy"]["aq"])
        res[f"pool_{mode}_boltz"] = ex_from(boltz, pl["La"], pl["Dy"], st_names["La"]["aq"], st_names["Dy"]["aq"])
    # (iv) every fixed quadruple
    quads = []
    for laL, dyL, laA, dyA in itertools.product(st_names["La"]["ac"], st_names["Dy"]["ac"], st_names["La"]["aq"], st_names["Dy"]["aq"]):
        g, n = exch((laL, dyL, laA, dyA))
        quads.append({"La_L": laL, "Dy_L": dyL, "La_aq": laA, "Dy_aq": dyA, "n": n, "dG": g,
                      "same_type": lab(laL) == lab(dyL) and laA[3:] == dyA[3:]})
    res["iv_all"] = quads
    same = [q for q in quads if q["same_type"]]
    res["iv_same_type"] = same
    vals = [q["dG"] for q in quads]
    print(f"\nTREATMENTS: (i) lowest {res['i_lowest']:+.2f} | (ii) Boltzmann {res['ii_boltzmann']:+.2f} | (iii) matched-labels {common}: "
          f"lowest {res['iii_matched_lowest']:+.2f}, Boltzmann {res['iii_matched_boltzmann']:+.2f}")
    print(f"  mono-only lowest {res['pool_acmono_lowest']:+.2f} boltz {res['pool_acmono_boltz']:+.2f} | bi-only lowest {res['pool_acbi_lowest']:+.2f} boltz {res['pool_acbi_boltz']:+.2f}")
    print(f"  (iv) all {len(quads)} quadruples: min {min(vals):+.2f} max {max(vals):+.2f}")
    print("  same-type quadruples (identical state labels for La and Dy):")
    for q in same:
        print(f"    {lab(q['La_L']):10s} {q['La_aq'][3:]}: dG {q['dG']:+.2f} (n={q['n']})")
    print("  Boltzmann populations:", json.dumps(res["ii_populations"]))
    # ---------------- per-water reference shift (standard-state convention x state choice)
    scan = []
    for dlt in np.arange(-4, 4.01, 1.0):
        Gs = {s: Gw[s] + (-(9 if "_aq" in s else 8) + A.nw(s)) * 0 for s in Gw}
        # shifting the free energy of a liquid water by dlt changes G_wbal by (ref - nw)*dlt
        Gs = {s: Gw[s] + ((9 if "_aq" in s else 8) - A.nw(s)) * dlt for s in Gw}
        lo2 = lambda names: min(Gs[x] for x in names)
        def bz2(names):
            g = np.array([Gs[x] for x in names]); g0 = g.min(); return g0 - RT * math.log(np.exp(-(g - g0) / RT).sum())
        a = ex_from(lo2, st_names["La"]["ac"], st_names["Dy"]["ac"], st_names["La"]["aq"], st_names["Dy"]["aq"])
        b = ex_from(bz2, st_names["La"]["ac"], st_names["Dy"]["ac"], st_names["La"]["aq"], st_names["Dy"]["aq"])
        selx = {f"{m}_{fam}": min(st_names[m][fam], key=lambda x: Gs[x]) for m in ("La", "Dy") for fam in ("aq", "ac")}
        scan.append({"delta": float(dlt), "lowest": a, "boltz": b, "selected": selx})
        print(f"  water-ref shift {dlt:+.0f}: lowest {a:+.2f} boltz {b:+.2f}  {selx}")
    res["water_ref_scan"] = scan
    # ---------------- selection effect: conformer noise and candidate counts
    noise = {}
    for n in [r["state"] for r in rows]:
        cands = []
        pr = ok[n]
        for c in sorted(W.parent.glob(f"conf*/{n}/dft_c*")):
            rc = A.parse(n, c)
            if rc and rc["normal"] and rc["G"] is not None and A.cn(rc) == A.cn(pr) and min(rc["imag"], default=0) > -30:
                cands.append(rc["G"] * H)
        if min(pr["imag"], default=0) > -30: cands.append(pr["G"] * H)
        cands = np.sort(np.array(cands))
        noise[n] = {"k": int(len(cands)), "spread_kcal": [round(float(x - cands[0]), 2) for x in cands[:8]] if len(cands) else []}
    res["candidates"] = noise
    print("\n  valid candidates per state (G above best, kcal/mol):")
    for n, v in noise.items(): print(f"    {n:14s} k={v['k']:2d} {v['spread_kcal']}")
    (OUT / "state_treatments.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
