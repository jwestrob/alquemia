"""One lambda window of the La3+ -> Dy3+ alchemical transformation, Li-Merz 12-6-4 / OPC.

The metal's LJ and C4 terms are removed from the prmtop system and replaced by one
CustomNonbondedForce over metal x partner pairs:
    U(lam) = (1-lam) U_La + lam U_Dy,   U_M = A_M/r^12 - B_M/r^6 - C_M/r^4
with A, B, C read from the La and Dy 12-6-4 prmtops (parmed add12_6_4). Charge (+3), mass
and all other terms are shared, so dU/dlam = U_Dy - U_La is the only quantity MBAR needs.
Recorded every 1 ps: dU/dlam, potential energy, volume, metal first-shell CN (O/N < 3.2 A),
<M-donor> and the number of non-water donors in the first shell.

usage: run_window.py <la.prmtop> <dy.prmtop> <rst7> <lambda> <T_K> <prod_ns> <outdir> [seed]
"""
import json, os, sys, time
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import app, unit
import parmed as pmd

KJ = 4.184
CUT_NM = 1.0
SHELL_A = 3.2  # first-shell cutoff for both metals (RDF minima are ~3.1-3.3 A for La/Dy-O)


def pair_coeffs(prmtop):
    """{partner LJ type name: (A, B, C)} for the metal (atom 0) in OpenMM units."""
    p = pmd.load_file(prmtop); r = p.parm_data
    nt, idx = r["POINTERS"][1], r["NONBONDED_PARM_INDEX"]
    mt = p.atoms[0].nb_idx
    tname = {a.nb_idx: a.type for a in p.atoms}
    out = {}
    for j in range(1, nt + 1):
        if j == mt: continue
        k = idx[nt * (mt - 1) + j - 1] - 1
        out[tname[j]] = (r["LENNARD_JONES_ACOEF"][k] * KJ * 1e-12, r["LENNARD_JONES_BCOEF"][k] * KJ * 1e-6,
                         r["LENNARD_JONES_CCOEF"][k] * KJ * 1e-4)
    return p, out


def build(la_top, dy_top, rst7):
    p, cla = pair_coeffs(la_top)
    _, cdy = pair_coeffs(dy_top)
    prm = app.AmberPrmtopFile(la_top)
    inp = app.AmberInpcrdFile(rst7)
    kw = dict(nonbondedMethod=app.PME, nonbondedCutoff=CUT_NM * unit.nanometer, constraints=app.HBonds, rigidWater=True)
    ref = prm.createSystem(**kw)
    sysm = prm.createSystem(**kw)
    nb = [f for f in sysm.getForces() if isinstance(f, mm.NonbondedForce)][0]
    q, s, e = nb.getParticleParameters(0)
    nb.setParticleParameters(0, q, s, 0.0)
    nb.setEwaldErrorTolerance(0.0005)
    c4 = [i for i, f in enumerate(sysm.getForces()) if isinstance(f, mm.CustomNonbondedForce)]
    for i in sorted(c4, reverse=True):  # parser's metal C4 force (all other C4 entries are zero)
        sysm.removeForce(i)
    f = mm.CustomNonbondedForce(
        "(1-lambda)*(a0/r^12 - b0/r^6 - c0/r^4) + lambda*(a1/r^12 - b1/r^6 - c1/r^4);"
        "a0=la1+la2; b0=lb1+lb2; c0=lc1+lc2; a1=da1+da2; b1=db1+db2; c1=dc1+dc2")
    for n in ("la", "lb", "lc", "da", "db", "dc"):
        f.addPerParticleParameter(n)
    f.addGlobalParameter("lambda", 0.0)
    f.addEnergyParameterDerivative("lambda")
    partners = []
    for a in p.atoms:
        if a.idx == 0:
            f.addParticle([0.0] * 6); continue
        A0 = cla.get(a.type, (0, 0, 0)); A1 = cdy.get(a.type, (0, 0, 0))
        f.addParticle(list(A0) + list(A1))
        if any(A0) or any(A1): partners.append(a.idx)
    for i in range(nb.getNumExceptions()):
        a, b, *_ = nb.getExceptionParameters(i)
        f.addExclusion(a, b)
    f.addInteractionGroup([0], partners)
    f.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    f.setCutoffDistance(CUT_NM)
    f.setUseLongRangeCorrection(False)
    f.setForceGroup(1)
    sysm.addForce(f)
    oxy = [a.idx for a in p.atoms if a.element in (7, 8)]
    lig = np.array([a.residue.name not in ("WAT", "HOH") for a in p.atoms])[oxy]
    return ref, sysm, prm, inp, oxy, lig


def platform():
    if os.environ.get("OMM_PLATFORM"):
        return mm.Platform.getPlatformByName(os.environ["OMM_PLATFORM"]), os.environ["OMM_PLATFORM"]
    for name in ("OpenCL",):  # CUDA: UNSUPPORTED_PTX_VERSION with this build/driver
        try:
            return mm.Platform.getPlatformByName(name), name
        except Exception:
            pass
    return mm.Platform.getPlatformByName("CPU"), "CPU"


def main():
    la_top, dy_top, rst7, lam, T, prod_ns, out = sys.argv[1:8]
    seed = int(sys.argv[8]) if len(sys.argv) > 8 else 0
    lam, T, prod_ns = float(lam), float(T), float(prod_ns)
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    ref, sysm, prm, inp, oxy, lig = build(la_top, dy_top, rst7)
    if (out / "meta.json").exists():
        return  # finished window (allows resubmission of a partly done task list)
    plat, pname = platform()
    props = {"Precision": "mixed"} if pname in ("CUDA", "OpenCL") else {}
    if pname == "CPU":
        props["Threads"] = os.environ.get("OMM_THREADS", "1")
    if pname == "OpenCL" and os.environ.get("OMM_DEVICE"):
        props["DeviceIndex"] = os.environ["OMM_DEVICE"]
    # correctness check: modified system at lambda=0 must reproduce the La 12-6-4 prmtop energy
    chk = {}
    for tag, s in (("ref", ref), ("mod", sysm)):
        c = mm.Context(s, mm.VerletIntegrator(0.001), plat, props)
        c.setPeriodicBoxVectors(*inp.boxVectors); c.setPositions(inp.positions)
        if tag == "mod": c.setParameter("lambda", 0.0)
        chk[tag] = c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        del c
    integ = mm.LangevinMiddleIntegrator(T * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    integ.setRandomNumberSeed(seed)
    baro = mm.MonteCarloBarostat(1.0 * unit.bar, T * unit.kelvin, 25)
    baro.setRandomNumberSeed(seed)
    sysm.addForce(baro)
    ctx = mm.Context(sysm, integ, plat, props)
    ctx.setPeriodicBoxVectors(*inp.boxVectors); ctx.setPositions(inp.positions)
    ctx.setParameter("lambda", lam)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 2000)
    ctx.setVelocitiesToTemperature(T * unit.kelvin, seed)
    t0 = time.time()
    integ.step(250_000)  # 0.5 ns NPT equilibration
    t_eq = time.time() - t0
    n = int(round(prod_ns * 1000))
    rec = np.zeros((n, 6))
    t0 = time.time()
    for k in range(n):
        integ.step(500)
        st = ctx.getState(getEnergy=True, getParameterDerivatives=True, getPositions=True)
        X = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        box = np.diag(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom))
        d = X[oxy] - X[0]; d -= box * np.round(d / box); r = np.sqrt((d * d).sum(1))
        sh = r[r < SHELL_A]
        nlig = int((r[lig] < SHELL_A).sum())
        rec[k] = (st.getEnergyParameterDerivatives()["lambda"], st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                  np.prod(box), len(sh), sh.mean() if len(sh) else np.nan, nlig)
        if (k + 1) % 500 == 0:
            np.save(out / "rec.npy", rec[:k + 1])
    t_prod = time.time() - t0
    np.save(out / "rec.npy", rec)
    meta = {"lambda": lam, "T": T, "prod_ns": prod_ns, "seed": seed, "platform": pname,
            "E_check_kJ": chk, "E_check_diff_kJ": chk["mod"] - chk["ref"],
            "ns_per_day": prod_ns / (t_prod / 86400), "eq_s": t_eq, "prod_s": t_prod,
            "columns": ["dU_dlam_kJ", "Epot_kJ", "V_A3", "CN_lt3.2", "mean_MO_A", "ligand_donors_lt3.2"]}
    (out / "meta.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
