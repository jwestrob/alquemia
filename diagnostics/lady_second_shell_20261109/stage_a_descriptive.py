"""Stage A, descriptive breakdown (POST HOC, written after the gate verdict was read; outside the frozen analysis).

The gate verdict comes from second_shell_analyze.py alone. This script only describes where the measured
second-shell shift sits: which free-energy term carries it, and how far each Dy cluster moved from its La start.

Terms (ORCA thermochemistry block, 298.15 K):  G = E_el + ZPE + thermal + (-TS)
  E_el    = "Electronic energy" (final single point incl. CPCM and D4)
  ZPE     = "Zero point energy"
  thermal = "Total Enthalpy" - E_el - ZPE
  -TS     = "Final Gibbs free energy" - "Total Enthalpy"
Dy-minus-La difference per term; the bare-ion energies enter the E_el term only.
Shift per term = second-shell pair (seed s) minus first-shell CN9 pair (La_aq9 conf/dft_c3, Dy_aq9 primary).

usage: stage_a_descriptive.py  -> table on stdout + workspaces/lady_second_shell_20261109/STAGE_A_DESCRIPTIVE.json
Saved outputs only; login node; single thread.
"""
import json, re
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2] if "__file__" in globals() else Path.cwd()
H = 627.509474
OLD = REPO / "workspaces/lady_relaxed_exchange_20261024"
NEW = REPO / "workspaces/lady_second_shell_20261109"
ION = REPO / "workspaces/lady_ingredient_tests_20261026/gas"
RUNS = {
    "La_aq9": OLD / "conf/La_aq9/dft_c3",
    "Dy_aq9": OLD / "runs/Dy_aq9",
    "La_ss_s1": OLD / "run_ss/La_ssaq_w21_s1",
    "La_ss_s2": OLD / "run_ss/La_ssaq_w21_s2",
    "Dy_ss_s1": NEW / "runs/Dy_ssaq_w21_s1",
    "Dy_ss_s2": NEW / "runs/Dy_ssaq_w21_s2",
}
PAT = {
    "E_el": r"Electronic energy\s+\.\.\.\s+(-?\d+\.\d+)",
    "ZPE": r"Zero point energy\s+\.\.\.\s+(-?\d+\.\d+)",
    "Hh": r"Total Enthalpy\s+\.\.\.\s+(-?\d+\.\d+)",
    "G": r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)",
}


def last(pat, t):
    m = re.findall(pat, t)
    if not m:
        raise SystemExit(f"pattern not found: {pat}")
    return float(m[-1])


def thermo(d):
    t = (d / "job.out").read_text(errors="replace")
    assert "ORCA TERMINATED NORMALLY" in t, d
    v = {k: last(p, t) for k, p in PAT.items()}
    v["thermal"] = v["Hh"] - v["E_el"] - v["ZPE"]
    v["mTS"] = v["G"] - v["Hh"]
    v["cycles"] = len(re.findall(r"GEOMETRY OPTIMIZATION CYCLE", t))
    fr = [float(f) for f in re.findall(r"^\s*\d+:\s+(-?\d+\.\d+) cm\*\*-1", t, re.M)]
    nz = sorted(f for f in fr if abs(f) > 1e-6)
    v["lowest_modes_cm1"] = [round(f, 1) for f in nz[:6]]
    v["n_modes_below_100"] = int(sum(1 for f in nz if f < 100))
    return v


def ion(m):
    t = (ION / f"{m}_ion/job.out").read_text(errors="replace")
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])


def xyz(d):
    L = [l.split() for l in (d / "job.xyz").read_text().splitlines()[2:] if l.strip()]
    return [l[0] for l in L], np.array([[float(x) for x in l[1:4]] for l in L])


def kabsch_rmsd(P, Q):
    P0, Q0 = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(P0.T @ Q0)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return float(np.sqrt(((P0 @ R - Q0) ** 2).sum() / len(P)))


def geom(d):
    el, X = xyz(d)
    r = np.linalg.norm(X - X[0], axis=1)
    o = [i for i, e in enumerate(el) if e == "O"]
    first = [i for i in o if r[i] < 2.95]
    second = [i for i in o if r[i] >= 2.95]
    oo = set()
    for a in range(len(o)):
        for b in range(a + 1, len(o)):
            if np.linalg.norm(X[o[a]] - X[o[b]]) < 3.2:
                oo.add((o[a], o[b]))
    # second-shell waters by number of first-shell O neighbours within 3.2 A
    donors = [sum(1 for f in first if np.linalg.norm(X[s] - X[f]) < 3.2) for s in second]
    return {"CN": len(first), "first_mean": float(r[first].mean()), "first_min": float(r[first].min()),
            "first_max": float(r[first].max()), "n_second": len(second),
            "second_mean": float(r[second].mean()) if second else None,
            "second_min": float(r[second].min()) if second else None,
            "second_max": float(r[second].max()) if second else None,
            "second_O_by_n_first_shell_neighbours": {str(k): donors.count(k) for k in sorted(set(donors))},
            "first_idx": first, "second_idx": second, "oo_edges": sorted(oo), "el": el, "X": X}


def main():
    T = {k: thermo(d) for k, d in RUNS.items()}
    eL, eD = ion("La"), ion("Dy")
    terms = ("E_el", "ZPE", "thermal", "mTS", "G")

    def pair(dy, la):
        out = {}
        for k in terms:
            v = (T[dy][k] - T[la][k]) * H
            if k in ("E_el", "G"):
                v -= (eD - eL) * H
            out[k] = v
        return out

    ref = pair("Dy_aq9", "La_aq9")
    res = {"note": "post hoc descriptive breakdown; the gate verdict is GATE.json alone",
           "first_shell_CN9": ref, "seeds": {}, "runs": {}}
    for k, v in T.items():
        res["runs"][k] = {a: (round(b, 8) if isinstance(b, float) else b) for a, b in v.items()}
    G = {k: geom(d) for k, d in RUNS.items()}
    for s in (1, 2):
        ss = pair(f"Dy_ss_s{s}", f"La_ss_s{s}")
        la, dy = G[f"La_ss_s{s}"], G[f"Dy_ss_s{s}"]
        assert la["el"][1:] == dy["el"][1:]
        heavy = [i for i, e in enumerate(la["el"]) if e != "H"]
        ela, edy = set(la["oo_edges"]), set(dy["oo_edges"])
        res["seeds"][f"s{s}"] = {
            "second_shell_pair": ss,
            "shift_by_term": {k: ss[k] - ref[k] for k in terms},
            "rmsd_Dy_final_vs_La_final_all_atoms_A": kabsch_rmsd(la["X"], dy["X"]),
            "rmsd_Dy_final_vs_La_final_heavy_A": kabsch_rmsd(la["X"][heavy], dy["X"][heavy]),
            "OO_contacts_below_3p2A": {"La": len(ela), "Dy": len(edy), "shared": len(ela & edy),
                                       "lost_in_Dy": len(ela - edy), "new_in_Dy": len(edy - ela)},
            "first_shell_membership_same": la["first_idx"] == dy["first_idx"],
            "opt_cycles_Dy": T[f"Dy_ss_s{s}"]["cycles"],
        }
    for k, g in G.items():
        res["runs"][k]["geometry"] = {a: b for a, b in g.items() if a not in ("first_idx", "second_idx", "oo_edges", "el", "X")}
    # seed-to-seed differences of each metal's cluster (what the pairing has to cancel)
    res["seed_difference_s2_minus_s1"] = {
        m: {k: (T[f"{m}_ss_s2"][k] - T[f"{m}_ss_s1"][k]) * H for k in terms} for m in ("La", "Dy")}

    def row(name, d):
        return f"{name:<34}" + "".join(f"{d[k]:>10.2f}" for k in terms)

    print("Dy-minus-La hydration difference by term (kcal/mol); experiment -66.9 for G")
    print(f"{'':<34}" + "".join(f"{k:>10}" for k in ("E_el", "ZPE", "thermal", "-TS", "G")))
    print(row("first shell, CN9 pair", ref))
    for s in (1, 2):
        print(row(f"second shell, seed {s}", res["seeds"][f"s{s}"]["second_shell_pair"]))
    for s in (1, 2):
        print(row(f"shift, seed {s}", res["seeds"][f"s{s}"]["shift_by_term"]))
    m = {k: np.mean([res["seeds"][f"s{s}"]["shift_by_term"][k] for s in (1, 2)]) for k in terms}
    res["mean_shift_by_term"] = {k: float(v) for k, v in m.items()}
    print(row("shift, mean of seeds", m))
    print()
    print("seed 2 minus seed 1, same metal (kcal/mol)")
    for mtl in ("La", "Dy"):
        print(row(f"{mtl} + 21 waters", res["seed_difference_s2_minus_s1"][mtl]))
    print()
    print("geometry")
    for k in RUNS:
        g = res["runs"][k]["geometry"]
        sec = "" if g["second_mean"] is None else f"  second shell n={g['n_second']} mean {g['second_mean']:.3f} [{g['second_min']:.3f}, {g['second_max']:.3f}]  by first-shell neighbours {g['second_O_by_n_first_shell_neighbours']}"
        print(f"{k:<10} CN {g['CN']}  M-O mean {g['first_mean']:.3f} [{g['first_min']:.3f}, {g['first_max']:.3f}]{sec}  cycles {res['runs'][k]['cycles']}  modes<100 cm-1 {res['runs'][k]['n_modes_below_100']}  lowest {res['runs'][k]['lowest_modes_cm1'][:3]}")
    print()
    for s in (1, 2):
        d = res["seeds"][f"s{s}"]
        print(f"seed {s}: Dy final vs La final RMSD {d['rmsd_Dy_final_vs_La_final_all_atoms_A']:.3f} A all atoms, {d['rmsd_Dy_final_vs_La_final_heavy_A']:.3f} A heavy; "
              f"O-O contacts {d['OO_contacts_below_3p2A']}; first-shell membership same: {d['first_shell_membership_same']}; Dy cycles {d['opt_cycles_Dy']}")
    (NEW / "STAGE_A_DESCRIPTIVE.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
