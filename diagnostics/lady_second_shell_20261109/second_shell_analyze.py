"""Stage A gate (PLAN.md, frozen). Saved outputs only; login node, single thread.

ddG_hyd,ss(s) = [G(Dy_ssaq_w21_s) - E(Dy3+)] - [G(La_ssaq_w21_s) - E(La3+)]   (kcal/mol; 21 waters each, no water term)
Reference by final Dy CN: 9 -> first-shell La_aq9/Dy_aq9 (-57.54); 8 -> La_aq9/Dy_aq8 + 1 water (-59.21, collect.py).
G and parsing from analyze.py (A.parse); bare ions from lady_ingredient_tests_20261026/gas/{La,Dy}_ion.
usage: second_shell_analyze.py            -> gate JSON on stdout + workspaces/.../GATE.json
       second_shell_analyze.py --selftest -> exercises the gate logic on synthetic shifts (no files written)
"""
import json, re, sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "diagnostics/lady_relaxed_exchange_20261024"))
import analyze as A  # noqa: E402

H, WC = A.H, A.WATER_CORR
LA = REPO / "workspaces/lady_relaxed_exchange_20261024/run_ss"
DY = REPO / "workspaces/lady_second_shell_20261109/runs"
ION = REPO / "workspaces/lady_ingredient_tests_20261026/gas"
OUT = REPO / "workspaces/lady_second_shell_20261109/GATE.json"
REF_EXPECT = {9: -57.54, 8: -59.21}
GATE, DISAGREE = -4.0, 1.0


def ion_E(m):
    t = (ION / f"{m}_ion/job.out").read_text(errors="replace")
    assert "ORCA TERMINATED NORMALLY" in t
    return float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)[-1])


def references(eL, eD):
    old = {p.name: A.parse(p.name) for p in sorted(A.W.iterdir())}
    ok = {n: r for n, r in old.items() if r and r["normal"] and r["G"] is not None}
    best, picked = A.conformer_best(ok)
    gw = best["water"]["G"] * H + WC
    r9 = (best["Dy_aq9"]["G"] - eD) * H - (best["La_aq9"]["G"] - eL) * H
    r8 = (best["Dy_aq8"]["G"] - eD) * H - (best["La_aq9"]["G"] - eL) * H + gw
    ref = {9: r9, 8: r8}
    for k, v in REF_EXPECT.items():
        if abs(ref[k] - v) > 0.01:
            raise SystemExit(json.dumps({"verdict": "STOP", "reason": f"reference CN{k} {ref[k]:.3f} != {v}"}))
    return ref, {k: picked[k]["chosen"] for k in ("La_aq9", "Dy_aq9", "Dy_aq8")}


def state(d, metal):
    """Admissibility per PLAN.md. Returns dict with ok flag and reasons."""
    r = A.parse(d.name, d)
    if r is None:
        return {"ok": False, "reasons": ["no job.out"]}
    why = []
    if not r["normal"]: why.append("not terminated normally")
    if not r["opt_converged"]: why.append("Opt not converged")
    if r["G"] is None: why.append("no Gibbs energy (Freq missing)")
    if min(r["imag"], default=0) <= -30: why.append(f"imaginary mode {min(r['imag'])} cm-1")
    L = [l.split() for l in (d / "job.xyz").read_text().splitlines()[2:] if l.strip()] if (d / "job.xyz").exists() else []
    cn, amb = None, []
    if L:
        X = np.array([[float(v) for v in l[1:4]] for l in L]); el = [l[0] for l in L]
        assert el[0] == metal, (d, el[0])
        dist = np.linalg.norm(X - X[0], axis=1)
        cn = int(sum(1 for e, x in zip(el, dist) if e == "O" and x < 2.95))
        amb = [round(float(x), 3) for e, x in zip(el, dist) if e == "O" and 2.95 <= x < 3.10]
        allowed = (9,) if metal == "La" else (8, 9)
        if cn not in allowed: why.append(f"CN {cn} not in {allowed}")
        if amb: why.append(f"O in 2.95-3.10 A band: {amb}")
    else:
        why.append("no job.xyz")
    return {"ok": not why, "reasons": why, "G": r["G"], "CN": cn, "O_2p95_3p10": amb, "imag": r["imag"]}


def gate(shifts):
    """shifts: list of two floats, or None for an inadmissible seed."""
    if len(shifts) != 2 or any(s is None for s in shifts):
        return "INCONCLUSIVE", "a seed is missing or inadmissible"
    m, d = float(np.mean(shifts)), abs(shifts[0] - shifts[1])
    both_pass, both_fail = all(s <= GATE for s in shifts), all(s > GATE for s in shifts)
    if m <= GATE and (d <= DISAGREE or both_pass):
        return "PASS", f"mean {m:+.2f} <= {GATE}, seed difference {d:.2f}"
    if m > GATE and (d <= DISAGREE or both_fail):
        return "FAIL", f"mean {m:+.2f} > {GATE}, seed difference {d:.2f}"
    return "INCONCLUSIVE", f"seeds differ by {d:.2f} > {DISAGREE} and straddle {GATE} (mean {m:+.2f})"


def selftest():
    cases = {(-5.0, -4.5): "PASS", (-6.0, -4.2): "PASS", (-3.0, -3.5): "FAIL", (-1.0, -2.5): "FAIL",
             (-4.5, -3.2): "INCONCLUSIVE", (-5.5, -3.0): "INCONCLUSIVE", (-4.0, -4.0): "PASS", (None, -6.0): "INCONCLUSIVE"}
    res = {str(k): gate(list(k))[0] for k in cases}
    bad = {k: v for k, v in res.items() if v != cases[eval(k)]}
    print(json.dumps({"selftest": res, "failures": bad}, indent=1))
    return not bad


def main():
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)
    eL, eD = ion_E("La"), ion_E("Dy")
    ref, ref_states = references(eL, eD)
    seeds, shifts = {}, []
    for s in (1, 2):
        la, dy = state(LA / f"La_ssaq_w21_s{s}", "La"), state(DY / f"Dy_ssaq_w21_s{s}", "Dy")
        rec = {"La": la, "Dy": dy}
        if la["ok"] and dy["ok"]:
            v = (dy["G"] - eD) * H - (la["G"] - eL) * H
            rec.update(ddG_hyd_ss=v, reference=ref[dy["CN"]], reference_CN=dy["CN"], shift=v - ref[dy["CN"]])
            shifts.append(rec["shift"])
        else:
            shifts.append(None)
        seeds[f"s{s}"] = rec
    verdict, why = gate(shifts)
    pending = [k for k, v in seeds.items() if v["Dy"]["reasons"] == ["no job.out"]]
    if pending:
        verdict, why = "PENDING", f"no Dy output yet for {pending}; gate not evaluated"
    out = {"plan": "diagnostics/lady_second_shell_20261109/PLAN.md", "references": {f"CN{k}": v for k, v in ref.items()},
           "reference_states": ref_states, "experiment": -66.9, "seeds": seeds,
           "mean_shift": None if None in shifts else float(np.mean(shifts)), "verdict": verdict, "why": why,
           "stage_B_released": verdict == "PASS"}
    txt = json.dumps(out, indent=1, default=float)
    print(txt)
    if verdict != "PENDING":
        OUT.write_text(txt)


if __name__ == "__main__":
    main()
