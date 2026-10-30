"""Counterpoise bound on the BSSE part of AO-block ALMO CT.

E_A(full basis): fragment A's SCF with the partner's basis functions as ghost
atoms (no nuclei, electrons or ECP). BSSE_A = E_A(full basis) - E_A(own basis) <= 0.
CT_cp = CT - (BSSE_metal + BSSE_shell)   (Khaliullin et al. 2007 convention).
"""
import argparse, json, os, sys, time
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "diagnostics/lady_almo_implementation_20261023"))
os.environ.setdefault("LADY_REPO_ROOT", str(REPO))
from ct_run import WS, METAL_FRAG, SHELL_FRAG


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--fragment", type=int, required=True, choices=(0, 1))
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    threads = int(os.environ.get("SLURM_CPUS_ON_NODE", "1"))
    from pyscf import gto, lib
    lib.num_threads(threads)
    import archive
    from run_system import configure_mean_field
    info = archive.load_system(a.system)
    frag = info["fragments"][a.fragment]
    real = set(frag["atom_indices"])
    atoms, basis = [], {}
    for i, (sym, xyz) in enumerate(info["atoms"]):
        label = sym if i in real else f"GHOST-{sym}"
        atoms.append((label, xyz))
        basis[label] = info["basis"][sym]
    mol = gto.Mole()
    mol.atom, mol.unit, mol.basis = atoms, "Angstrom", basis
    mol.ecp = {s: v for s, v in info["ecp"].items() if any(atoms[i][0] == s for i in real)}
    mol.charge, mol.spin, mol.cart, mol.symmetry = frag["charge"], 0, False, False
    mol.max_memory = int(0.85 * 4000 * threads)
    mol.verbose, mol.output = 4, str(out / "cp_pyscf.log")
    mol.build()
    if mol.nelectron != frag["nelectron"]:
        raise ValueError("electron count mismatch")
    mf = configure_mean_field(mol, out / "cp.chk")
    mf.max_memory = mol.max_memory
    # initial density: saved isolated-fragment density embedded in the real-atom AO block
    shell = "aquo" if a.system.startswith("aquo") else "mixed"
    src = WS / (METAL_FRAG[info["metal"]] if a.fragment == 0 else SHELL_FRAG[shell]) / "result/fragment_state.npz"
    z = np.load(src)
    sl = mol.aoslice_by_atom()
    idx = [j for i in sorted(real) for j in range(sl[i, 2], sl[i, 3])]
    dm0 = np.zeros((mol.nao_nr(),) * 2); dm0[np.ix_(idx, idx)] = z["density"]
    t = time.monotonic()
    e = float(mf.kernel(dm0=dm0))
    res = {"system": a.system, "fragment": a.fragment, "converged": bool(mf.converged),
           "E_full_basis": e, "E_own_basis": float(z["e_tot"]), "BSSE_hartree": e - float(z["e_tot"]),
           "nao": mol.nao_nr(), "elapsed_s": time.monotonic() - t}
    (out / "RESULT.json").write_text(json.dumps(res, indent=2) + "\n")
    print(res)


if __name__ == "__main__":
    main()
