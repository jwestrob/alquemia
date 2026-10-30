"""Second-hydration-shell clusters for the La/Dy exchange ladder.

Core = lowest-G first-shell structure per state (analyze.py --conformers selection).
Twelve second-shell waters are placed as H-bond acceptors: candidate O sites at
2.80 A from a first-shell water O (or a non-coordinated acetate O, as donor-side
sites along its lone-pair direction) along each O-H bond, chosen by farthest-point
sampling from a seeded random start; each water's H atoms point away from the metal
with a random twist. The whole cluster is GFN2-xTB/ALPB(water) pre-optimized; DFT
(orca_job.sh, run_ss/<name>/job.inp) relaxes it freely. Run inside SLURM: xtb jobs
run in parallel, one thread each.

Naming: <M>_ss<lig>_w<total waters>_s<seed>, e.g. La_ssaq_w21_s1 (9 + 12).
"""
import json, os, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
WS = REPO / "workspaces/lady_relaxed_exchange_20261024"
XTB = "/home/jwestrob/miniconda3/envs/lanm_bench/bin/xtb"
N2 = 12
SEEDS = (1, 2, 3)
STATES = ["La_aq8", "La_aq9", "Dy_aq8", "Dy_aq9",
          "La_acmono_w7", "La_acmono_w8", "Dy_acmono_w7", "Dy_acmono_w8",
          "La_acbi_w6", "La_acbi_w7", "Dy_acbi_w6", "Dy_acbi_w7"]
CHARGE = {"aq": 3, "acmono": 2, "acbi": 2}


def read_xyz(p):
    L = [l.split() for l in Path(p).read_text().splitlines()[2:] if l.strip()]
    return [r[0] for r in L], np.array([[float(x) for x in r[1:4]] for r in L])


def write_xyz(p, el, X, c=""):
    Path(p).write_text(f"{len(el)}\n{c}\n" + "".join(f"{e:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n" for e, (x, y, z) in zip(el, X)))


def core_path(state):
    sel = json.loads((WS / "EXCHANGE_RESULT_CONFORMERS.json").read_text())["conformer_selection"][state]["chosen"]
    return WS / "runs" / state / "job.xyz" if sel == "primary" else WS / sel.split("/")[0] / state / sel.split("/")[1] / "job.xyz"


def sites(el, X):
    M = X[0]; out = []
    for i, e in enumerate(el):
        if e != "O" or i == 0: continue
        Hs = [j for j, f in enumerate(el) if f == "H" and np.linalg.norm(X[j] - X[i]) < 1.2]
        if len(Hs) == 2:  # water: accept along each O-H
            for h in Hs:
                u = (X[h] - X[i]) / np.linalg.norm(X[h] - X[i]); out.append((X[i] + 2.80 * u, X[i]))
        elif not Hs and np.linalg.norm(X[i] - M) > 3.0:  # free carboxylate O: donor site away from metal
            u = X[i] - M; u /= np.linalg.norm(u); out.append((X[i] + 2.80 * u, X[i]))
    return out


def place(el, X, seed):
    rng = np.random.default_rng(seed); M = X[0]
    heavy = np.array([x for e, x in zip(el, X) if e != "H"]); hyd = np.array([x for e, x in zip(el, X) if e == "H"])
    cand = [(p, a) for p, a in sites(el, X) if np.min(np.linalg.norm(heavy - p, axis=1)) > 2.5
            and np.min(np.linalg.norm(hyd - p, axis=1)) > 1.5]
    chosen = [cand[rng.integers(len(cand))]]
    while len(chosen) < N2:
        rest = [c for c in cand if all(np.linalg.norm(c[0] - d[0]) > 2.4 for d in chosen)]
        if not rest:  # exhausted: radial fill at 4.6 A from the metal
            v = rng.normal(size=3); v /= np.linalg.norm(v); rest = [(M + 4.6 * v, M)]
        dists = [min(np.linalg.norm(c[0] - d[0]) for d in chosen) for c in rest]
        chosen.append(rest[int(np.argmax(np.array(dists) + 0.3 * rng.random(len(dists))))])
    el2, X2 = list(el), [x for x in X]
    for O, anchor in chosen:
        out = O - M; out /= np.linalg.norm(out)
        t = rng.normal(size=3); t -= t.dot(out) * out; t /= np.linalg.norm(t)
        a = np.deg2rad(104.52 / 2)
        for s in (1, -1):
            el2.append("H"); X2.append(O + 0.9572 * (np.cos(a) * out + s * np.sin(a) * t))
        el2.insert(len(el2) - 2, "O"); X2.insert(len(X2) - 2, O)
    return el2, np.array(X2)


def build(args):
    state, seed = args
    m, rest = state.split("_", 1)
    lig = "aq" if rest.startswith("aq") else rest.split("_")[0]
    nw1 = int(rest[2:]) if lig == "aq" else int(rest.split("_w")[1])
    name = f"{m}_ss{lig}_w{nw1 + N2}_s{seed}"
    el, X = read_xyz(core_path(state))
    el2, X2 = place(el, X, 1000 * seed + sum(map(ord, state)))
    q = CHARGE[lig]
    # v2: unconstrained GFN2 pulled second-shell waters into the first shell (CN 8 -> 9-10),
    # so the core (metal + first shell, atoms 1..len(el)) is frozen; if CN still changes,
    # the raw placement is kept and DFT relaxes it.
    def cnum(e, Y):
        return sum(1 for i in range(1, len(e)) if e[i] in ("O", "N") and np.linalg.norm(Y[i] - Y[0]) < 2.95)
    cn0 = cnum(el2, X2)
    with tempfile.TemporaryDirectory() as d:
        write_xyz(f"{d}/in.xyz", el2, X2)
        Path(f"{d}/fix.inp").write_text(f"$fix\n atoms: 1-{len(el)}\n$end\n")
        p = subprocess.run([XTB, "in.xyz", "--opt", "--gfn", "2", "--alpb", "water", "--chrg", str(q), "-P", "1", "--input", "fix.inp"],
                           cwd=d, capture_output=True, text=True, env={**os.environ, "OMP_NUM_THREADS": "1"})
        ok = p.returncode == 0 and Path(f"{d}/xtbopt.xyz").exists()
        if ok:
            e3, X3 = read_xyz(f"{d}/xtbopt.xyz")
            ok = cnum(e3, X3) == cn0
            if ok: el2, X2 = e3, X3
    R = WS / "run_ss" / name
    write_xyz(R / "start.xyz", el2, X2, f"charge {q} core {state} seed {seed} xtb_ok {ok}")
    d1 = sorted(np.linalg.norm(X2[1:] - X2[0], axis=1)[[i - 1 for i in range(1, len(el2)) if el2[i] == "O"]])
    return name, {"core": str(core_path(state)), "charge": q, "natoms": len(el2), "xtb_ok": ok,
                  "first_shell_O_lt_2p95": int(sum(x < 2.95 for x in d1)), "nearest_O": [round(float(x), 3) for x in d1[:10]]}


def main():
    jobs = [(s, k) for s in STATES for k in SEEDS]
    if "--names" in sys.argv:  # print run names only (for input generation before the build job runs)
        for s, k in jobs:
            m, rest = s.split("_", 1); lig = "aq" if rest.startswith("aq") else rest.split("_")[0]
            nw1 = int(rest[2:]) if lig == "aq" else int(rest.split("_w")[1])
            print(f"{m}_ss{lig}_w{nw1 + N2}_s{k} {CHARGE[lig]} {s}")
        return
    with ProcessPoolExecutor(int(os.environ.get("SLURM_CPUS_ON_NODE", "1"))) as ex:
        meta = dict(ex.map(build, jobs))
    (WS / "run_ss" / "SYSTEMS.json").write_text(json.dumps(meta, indent=1))
    for k, v in meta.items(): print(k, v["xtb_ok"], v["first_shell_O_lt_2p95"], v["natoms"])


if __name__ == "__main__":
    main()
