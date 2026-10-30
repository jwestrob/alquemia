"""Build per-metal starting geometries for the relaxed La/Dy exchange ladder.

Ligands are rigidly translated so donor O sit at a metal-specific start radius
(La 2.56 A, Dy 2.42 A); final geometries come from PBE0-D4/CPCM optimization,
so these only need to be sane. EDTA is wrapped with restrained GFN2-xTB.
"""
import json, subprocess, sys, tempfile
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
REF = REPO / "workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1"
OUT = REPO / "workspaces/lady_relaxed_exchange_20261024/start"
XTB = "/home/jwestrob/miniconda3/envs/lanm_bench/bin/xtb"
R0 = {"La": 2.56, "Dy": 2.42}
R_N = {"La": 2.80, "Dy": 2.62}


def read_xyz(p):
    L = [l.split() for l in Path(p).read_text().splitlines()[2:] if l.strip()]
    return [r[0] for r in L], np.array([[float(x) for x in r[1:4]] for r in L])


def write_xyz(p, el, X, comment=""):
    Path(p).write_text(f"{len(el)}\n{comment}\n" + "".join(
        f"{e:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n" for e, (x, y, z) in zip(el, X)))


def molecules(el, X):
    """Connected components excluding the metal (index 0)."""
    n = len(el); idx = list(range(1, n)); seen, comps = set(), []
    for s in idx:
        if s in seen: continue
        stack, comp = [s], []
        while stack:
            i = stack.pop()
            if i in seen: continue
            seen.add(i); comp.append(i)
            for j in idx:
                if j not in seen and np.linalg.norm(X[i] - X[j]) < 1.75:
                    stack.append(j)
        comps.append(sorted(comp))
    return comps


def rescale(el, X, metal):
    X = X.copy(); M = X[0]
    for comp in molecules(el, X):
        donors = [i for i in comp if el[i] == "O" and np.linalg.norm(X[i] - M) < 3.0]
        if not donors: continue
        c = X[donors].mean(0) - M
        r = np.mean([np.linalg.norm(X[i] - M) for i in donors])
        X[comp] += c / np.linalg.norm(c) * (R0[metal] - r)
    return X


def add_water(el, X, metal):
    M = X[0]; k = np.arange(2000) + 0.5
    phi = np.arccos(1 - 2 * k / 2000); th = np.pi * (1 + 5 ** 0.5) * k
    dirs = np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)
    best = max(dirs, key=lambda d: np.min(np.linalg.norm(X[1:] - (M + R0[metal] * d), axis=1)))
    O = M + R0[metal] * best
    perp = np.cross(best, [0, 0, 1] if abs(best[2]) < 0.9 else [1, 0, 0]); perp /= np.linalg.norm(perp)
    a = np.deg2rad(104.52 / 2); H1 = O + 0.9572 * (np.cos(a) * best + np.sin(a) * perp)
    H2 = O + 0.9572 * (np.cos(a) * best - np.sin(a) * perp)
    return el + ["O", "H", "H"], np.vstack([X, O, H1, H2])


def xtb_opt(el, X, charge, constraints=None, tag="x"):
    with tempfile.TemporaryDirectory() as d:
        write_xyz(f"{d}/in.xyz", el, X)
        args = [XTB, "in.xyz", "--opt", "--gfn", "2", "--alpb", "water", "--chrg", str(charge), "-P", "1"]
        if constraints:
            Path(f"{d}/c.inp").write_text("$constrain\n force constant=2.0\n" + "".join(
                f" distance: {i+1}, {j+1}, {r:.3f}\n" for i, j, r in constraints) + "$end\n")
            args += ["--input", "c.inp"]
        p = subprocess.run(args, cwd=d, capture_output=True, text=True)
        if p.returncode != 0 or not Path(f"{d}/xtbopt.xyz").exists():
            raise RuntimeError(f"xtb failed for {tag}: {p.stdout[-2000:]}")
        return read_xyz(f"{d}/xtbopt.xyz")


def edta(metal, nwater):
    sdf = (REPO / "workspaces/lady_relaxed_exchange_20261024/src_untrusted/edta.sdf").read_text().splitlines()
    na, nb = int(sdf[3][:3]), int(sdf[3][3:6])
    el = [l.split()[3] for l in sdf[4:4 + na]]
    X = np.array([[float(v) for v in l.split()[:3]] for l in sdf[4:4 + na]])
    bonds = [(int(l[:3]) - 1, int(l[3:6]) - 1) for l in sdf[4 + na:4 + na + nb]]
    acidic = {b if el[a] == "O" else a for a, b in bonds if {el[a], el[b]} == {"O", "H"}}
    keep = [i for i in range(na) if i not in acidic]
    nbr = {i: set() for i in range(na)}
    for a, b in bonds: nbr[a].add(b); nbr[b].add(a)
    N = [i for i in keep if el[i] == "N"]
    carbC = [i for i in keep if el[i] == "C" and sum(el[j] == "O" for j in nbr[i]) == 2]
    donorO = [sorted([j for j in nbr[c] if el[j] == "O"])[0] for c in carbC]
    new = {old: k + 1 for k, old in enumerate(keep)}
    el2 = [metal] + [el[i] for i in keep]
    M = X[N + carbC].mean(0)
    X2 = np.vstack([M, X[keep]])
    cons = [(0, new[i], R_N[metal]) for i in N] + [(0, new[o], R0[metal]) for o in donorO]
    el2, X2 = xtb_opt(el2, X2, -1, cons, f"edta_wrap_{metal}")
    for _ in range(nwater):
        el2, X2 = add_water(el2, X2, metal)
    return xtb_opt(el2, X2, 3 - 4, cons, f"edta_{metal}_w{nwater}")  # restraints kept; DFT relaxes freely


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    systems = {}
    for m in ("La", "Dy"):
        el, X = read_xyz(REF / "prepared/xyz/aquo8_central_La.xyz"); el[0] = m
        systems[f"{m}_aq8"] = (el, rescale(el, X, m), 3)
        el, X = read_xyz(REF / "prepared/xyz/aquo9_central_La.xyz"); el[0] = m
        systems[f"{m}_aq9"] = (el, rescale(el, X, m), 3)
        for fam, base in (("acmono", 7), ("acbi", 6)):
            name = "mono" if fam == "acmono" else "bi"
            el, X = read_xyz(REF / f"hydrated_mixed_v1/prepared_v2/xyz/hydrated_acetate_{name}_central_La.xyz"); el[0] = m
            X = rescale(el, X, m)
            systems[f"{m}_{fam}_w{base}"] = (el, X, 2)
            systems[f"{m}_{fam}_w{base+1}"] = (*add_water(el, X, m), 2)
        for nw in (2, 3):
            try:
                el, X = edta(m, nw)
            except RuntimeError:  # GFN2 SCF can diverge for Dy-EDTA; start from La-EDTA, DFT relaxes it
                el, X = systems[f"La_edta_w{nw}"][0][:], systems[f"La_edta_w{nw}"][1]; el[0] = "Dy"
                print("Dy EDTA fallback to La start", nw)
            systems[f"{m}_edta_w{nw}"] = (el, X, -1)
    systems["water"] = (["O", "H", "H"], np.array([[0, 0, 0], [0.7572, 0.5865, 0], [-0.7572, 0.5865, 0]]), 0)
    meta = {}
    for k, (el, X, q) in systems.items():
        write_xyz(OUT / f"{k}.xyz", el, X, f"charge {q}")
        d = sorted(np.linalg.norm(X[1:] - X[0], axis=1)[[i - 1 for i in range(1, len(el)) if el[i] in "ON"]])[:10] if k != "water" else []
        meta[k] = {"charge": q, "natoms": len(el), "nearest_donor_A": [round(float(x), 3) for x in d]}
        print(k, meta[k])
    (OUT / "SYSTEMS.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
