"""Independent check of the almo.py constrained minimum.

Method A: almo.solve_almo (analytic gradient, Grassmann L-BFGS, polar retraction).
Method B: scipy L-BFGS on an occupied->virtual mixing parametrization of each
fragment block, C_A = C_A^occ + C_A^vir X_A (isolated-fragment orbitals), with
density built here from scratch, D = 2 C (C^T S C)^-1 C^T, and a central
finite-difference gradient of the backend energy only. Nothing from almo.py
(derivative, retraction, preconditioner, normalization) is used by B.

Pass criterion, declared before running: |P_A - P_B| <= 1e-6 Eh for each case,
and the ordering U <= P <= F. Also: CT -> ~0 when the dimer is separated by 10 A.
"""
import json, os, sys, time
from pathlib import Path
import numpy as np
from scipy import optimize

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "diagnostics/lady_almo_implementation_20261023"))
from pyscf import gto, scf, dft, lib
import almo

lib.num_threads(int(os.environ.get("SLURM_CPUS_ON_NODE", "4")))

WATER_DIMER = [("O", (-1.551007, -0.114520, 0.000000)), ("H", (-1.934259, 0.762503, 0.000000)),
               ("H", (-0.599677, 0.040712, 0.000000)), ("O", (1.350625, 0.111469, 0.000000)),
               ("H", (1.680398, -0.373741, -0.758561)), ("H", (1.680398, -0.373741, 0.758561))]
MG_WATER = [("Mg", (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, 2.05)),
            ("H", (0.0, 0.763, 2.64)), ("H", (0.0, -0.763, 2.64))]


def make_mf(atoms, charge, xc):
    mol = gto.M(atom=atoms, basis="6-31g*", charge=charge, verbose=0, max_memory=20000)
    if xc == "HF":
        mf = scf.RHF(mol)
    else:
        mf = dft.RKS(mol); mf.xc = xc; mf.grids.atom_grid = (75, 302); mf.grids.prune = None
    mf.conv_tol, mf.conv_tol_grad = 1e-12, 1e-8
    return mf


def run_case(name, atoms, frags, charges, xc):
    mf = make_mf(atoms, sum(charges), xc)
    U = mf.kernel()
    fmfs, R = [], 0.0
    for idx, q in zip(frags, charges):
        f = make_mf([atoms[i] for i in idx], q, xc); R += f.kernel(); fmfs.append(f)
    sl = mf.mol.aoslice_by_atom()
    aos = [[j for i in idx for j in range(sl[i, 2], sl[i, 3])] for idx in frags]
    nocc = [f.mol.nelectron // 2 for f in fmfs]

    # Method A
    tA = time.time()
    solA = almo.solve_almo(mf, aos, nocc, fragment_mfs=fmfs, gradient_tolerance=1e-8,
                           energy_tolerance=1e-12, max_iterations=2000, max_evaluations=10000)
    tA = time.time() - tA

    # Method B (independent)
    S = mf.get_ovlp(); nao = S.shape[0]
    blocks = []
    for f, a, n in zip(fmfs, aos, nocc):
        blocks.append((np.array(a), f.mo_coeff[:, :n], f.mo_coeff[:, n:]))
    sizes = [cv.shape[1] * co.shape[1] for _, co, cv in blocks]

    def coeffs(x):
        cols, k = [], 0
        for (a, co, cv), m in zip(blocks, sizes):
            X = x[k:k + m].reshape(cv.shape[1], co.shape[1]); k += m
            c = np.zeros((nao, co.shape[1])); c[a] = co + cv @ X; cols.append(c)
        return np.hstack(cols)

    def energy(x):
        C = coeffs(x)
        D = 2 * C @ np.linalg.solve(C.T @ S @ C, C.T)
        return mf.energy_tot(dm=D)

    h = 1e-5
    nev = [0]
    def fun(x):
        e = energy(x); g = np.empty_like(x)
        for i in range(len(x)):
            xp = x.copy(); xp[i] += h; xm = x.copy(); xm[i] -= h
            g[i] = (energy(xp) - energy(xm)) / (2 * h)
        nev[0] += 1 + 2 * len(x)
        return e, g
    x0 = np.zeros(sum(sizes))
    F_B = energy(x0)
    tB = time.time()
    solB = optimize.minimize(fun, x0, jac=True, method="L-BFGS-B",
                             options={"maxiter": 5000, "ftol": 1e-15, "gtol": 1e-8, "maxcor": 30})
    tB = time.time() - tB
    P_A, P_B, F_A = solA.energy, float(solB.fun), solA.initial_energy
    out = {"case": name, "xc": xc, "nao": nao, "nparam_B": len(x0),
           "U": U, "R": R, "F_A": F_A, "F_B": F_B, "P_A": P_A, "P_B": P_B,
           "P_A_minus_P_B": P_A - P_B, "F_A_minus_F_B": F_A - F_B,
           "CT_kcal": (U - P_A) * 627.509474, "POL_kcal": (P_A - F_A) * 627.509474,
           "A": {"status": solA.status, "iterations": solA.iterations, "grad": solA.gradient_norm, "s": tA},
           "B": {"message": str(solB.message), "nit": solB.nit, "energy_evals": nev[0], "s": tB},
           "pass": bool(abs(P_A - P_B) <= 1e-6 and U <= P_A + 1e-9 and P_A <= F_A + 1e-9)}
    print(json.dumps(out, default=float), flush=True)
    return out


def main():
    outdir = Path(sys.argv[1]); outdir.mkdir(parents=True, exist_ok=True)
    far = WATER_DIMER[:3] + [(s, (x + 10.0, y, z)) for s, (x, y, z) in WATER_DIMER[3:]]
    cases = [("water_dimer", WATER_DIMER, [[0, 1, 2], [3, 4, 5]], [0, 0], "HF"),
             ("water_dimer", WATER_DIMER, [[0, 1, 2], [3, 4, 5]], [0, 0], "PBE0"),
             ("water_dimer_far10A", far, [[0, 1, 2], [3, 4, 5]], [0, 0], "PBE0"),
             ("Mg2+_water", MG_WATER, [[0], [1, 2, 3]], [2, 0], "HF"),
             ("Mg2+_water", MG_WATER, [[0], [1, 2, 3]], [2, 0], "PBE0")]
    results = [run_case(*c) for c in cases]
    far_ct = abs(results[2]["CT_kcal"])
    summary = {"cases": results, "far_dimer_CT_kcal": far_ct,
               "all_pass": all(r["pass"] for r in results) and far_ct < 1e-3}
    (outdir / "RESULT.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print("ALL_PASS", summary["all_pass"])


if __name__ == "__main__":
    main()
