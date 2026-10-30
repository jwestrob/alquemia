"""One quartet endpoint: full (U), frozen (F), polarized (P) electronic energies.

Same PySCF backend/grid/basis/ECP as the 2026-10-23 baselines. U comes from the
saved baseline orbitals (re-evaluated once on this backend), or from a short
SCF polish restarted from those orbitals (--polish). F and P use the existing
AO-block ALMO solver warm-started from the saved isolated-fragment orbitals.
"""
import argparse, json, os, sys, time
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OLD = REPO / "diagnostics/lady_almo_implementation_20261023"
WS = REPO / "workspaces/lady_almo_implementation_20261023"
sys.path.insert(0, str(OLD))
os.environ.setdefault("LADY_REPO_ROOT", str(REPO))

BASELINES = {"aquo8_central_La": "baseline_aquo_La_v1",
             "aquo8_central_Dy": "baseline_aquo_Dy_1e_v2",
             "hydrated_acetate_mono_central_La": "baseline_mixed_La_v1",
             "hydrated_acetate_mono_central_Dy": "baseline_mixed_Dy_1e_v2"}
METAL_FRAG = {"La": "fragment_La_v1", "Dy": "fragment_Dy_1e_v2"}
SHELL_FRAG = {"aquo": "fragment_aquo_shell_v1", "mixed": "fragment_mixed_shell_v1"}


def dump(path, obj):
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, default=float) + "\n")
    tmp.replace(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True, choices=list(BASELINES))
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--polish", action="store_true",
                    help="restart ordinary SCF from saved baseline orbitals")
    ap.add_argument("--gradient-tolerance", type=float, default=1e-6)
    ap.add_argument("--energy-tolerance", type=float, default=1e-9)
    a = ap.parse_args()
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    threads = int(os.environ.get("SLURM_CPUS_ON_NODE", "1"))
    mem_mb = int(os.environ.get("SLURM_MEM_PER_NODE", str(4000 * threads)))
    max_memory = int(0.85 * mem_mb)

    from pyscf import lib
    lib.num_threads(threads)
    import archive, almo
    from run_system import configure_mean_field

    t0 = time.monotonic()
    res = {"system": a.system, "job": os.environ.get("SLURM_JOB_ID"), "threads": threads,
           "max_memory_mb": max_memory, "status": "started"}
    dump(out / "RESULT.json", res)

    mol, info = archive.build_molecule(a.system, out / "full_pyscf.log", max_memory)
    mf = configure_mean_field(mol, out / "full.chk")
    mf.max_memory = max_memory
    base = np.load(WS / BASELINES[a.system] / "result/baseline_state.npz")
    saved_E = float(base["e_tot"])

    if a.polish:
        mf.conv_tol_grad = 2e-8  # PySCF's post-convergence extra diagonalization can raise |g|
        hist = []
        mf.callback = lambda env: hist.append({k: float(env[k]) for k in ("cycle", "e_tot", "norm_gorb") if k in env})
        U = float(mf.kernel(dm0=base["density"]))
        dm = mf.make_rdm1()
        g = float(np.linalg.norm(mf.get_grad(mf.mo_coeff, mf.mo_occ, mf.get_fock(dm=dm))))
        res["U_source"] = {"mode": "polish_from_saved", "converged": bool(mf.converged),
                           "cycles": hist, "final_gradient": g,
                           "change_from_saved": U - saved_E}
        if not mf.converged or g > 1e-7:
            res["status"] = "failed_polish"
            dump(out / "RESULT.json", res)
            return 2
    else:
        dm = base["density"]
        U = float(mf.energy_tot(dm=dm))
        res["U_source"] = {"mode": "re-evaluated_saved_density", "saved_energy": saved_E,
                           "replay_minus_saved": U - saved_E}
    res["timing_U_s"] = time.monotonic() - t0
    dump(out / "RESULT.json", res)

    # fragments: saved isolated SCF orbitals, own basis (no ghosts)
    shell_key = "aquo" if a.system.startswith("aquo") else "mixed"
    frag_dirs = [WS / METAL_FRAG[info["metal"]], WS / SHELL_FRAG[shell_key]]
    frag_mfs, R = [], 0.0
    for i, d in enumerate(frag_dirs):
        fmol, finfo = archive.build_molecule(a.system, out / f"frag{i}_pyscf.log", max_memory, fragment=i)
        fmf = configure_mean_field(fmol, out / f"frag{i}.chk")
        z = np.load(d / "result/fragment_state.npz")
        if not np.allclose(fmf.get_ovlp(), z["overlap"], atol=1e-11):
            raise ValueError(f"fragment {i} AO overlap differs from saved state")
        fmf.mo_coeff, fmf.mo_occ, fmf.converged = z["mo_coeff"], z["mo_occ"], True
        frag_mfs.append(fmf)
        R += float(z["e_tot"])
    slices = mol.aoslice_by_atom()
    frag_aos = [[j for i in f["atom_indices"] for j in range(int(slices[i, 2]), int(slices[i, 3]))]
                for f in info["fragments"]]
    frag_nocc = [f["nelectron"] // 2 for f in info["fragments"]]

    hist_path = out / "almo_history.jsonl"
    hist_path.unlink(missing_ok=True)
    def progress(row):
        row = {k: v for k, v in row.items() if k != "invariants"}
        with hist_path.open("a") as s:
            s.write(json.dumps(row, default=float) + "\n")

    t1 = time.monotonic()
    sol = almo.solve_almo(mf, frag_aos, frag_nocc, fragment_mfs=frag_mfs,
                          max_iterations=1000, max_evaluations=5000,
                          gradient_tolerance=a.gradient_tolerance,
                          energy_tolerance=a.energy_tolerance, callback=progress)
    np.savez_compressed(out / "polarized_state.npz", coefficients=sol.coefficients, density=sol.density)
    F, P = float(sol.initial_energy), float(sol.energy)
    res.update({
        "status": "passed" if sol.converged and P <= F + 1e-8 and U <= P + 1e-8 else "failed_almo",
        "solver": {"status": sol.status, "message": sol.message, "iterations": sol.iterations,
                   "evaluations": sol.evaluations, "gradient_norm": sol.gradient_norm,
                   "elapsed_s": time.monotonic() - t1},
        "energies_hartree": {"R": R, "F": F, "P": P, "U": U,
                             "FRZ": F - R, "POL": P - F, "CT": U - P,
                             "D4_archived": info["archive_energies"]["d4_printed_hartree"],
                             "archive_electronic": info["archive_energies"]["electronic_hartree"]},
        "elapsed_s": time.monotonic() - t0})
    dump(out / "RESULT.json", res)
    return 0 if res["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
