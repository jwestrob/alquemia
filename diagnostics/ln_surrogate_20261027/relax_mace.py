"""MACE-OMol-0 gas-phase relaxation of La and Dy forms of each built complex (Dy by metal swap of the
relaxed La structure). Reads <dir>/La_cn*.xyz, writes <dir>/{La,Dy}_cn*_mace.xyz and <dir>/energies.json
(kcal/mol, plus coordination summary: donors within 3.1 A and how many are water O).

usage: relax_mace.py <dir> [threads]
"""
import json, sys
from pathlib import Path
import numpy as np
import torch
from ase.io import read
from ase.optimize import BFGS
from mace.calculators import MACECalculator

CKPT = "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/software_v1/MACE-omol-0-extra-large-1024.model"
EV = 23.060548
d = Path(sys.argv[1]); torch.set_num_threads(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
calc = MACECalculator(model_paths=CKPT, device="cpu", default_dtype="float64", head="omol")


def coord(a):
    i = [k for k, s in enumerate(a.get_chemical_symbols()) if s in ("La", "Dy")][0]
    r = np.linalg.norm(a.positions - a.positions[i], axis=1)
    sym = a.get_chemical_symbols()
    don = [k for k in range(len(a)) if sym[k] in ("O", "N") and r[k] < 3.1]
    # water O: O bonded to two H and nothing else
    def is_water(k):
        rr = np.linalg.norm(a.positions - a.positions[k], axis=1)
        nb = [j for j in range(len(a)) if j != k and j != i and rr[j] < 1.25]
        return sym[k] == "O" and len(nb) == 2 and all(sym[j] == "H" for j in nb)
    nw = sum(is_water(k) for k in don)
    return {"CN": len(don), "n_water": int(nw), "n_ligand": len(don) - int(nw), "M_donor_mean_A": float(r[don].mean()) if don else None}


res = json.loads((d / "energies.json").read_text()) if (d / "energies.json").exists() else {}
for f in sorted(d.glob("La_cn?.xyz")):
    cn = f.stem.split("_")[1]
    q = int(open(f).readlines()[1].split("charge=")[1].split()[0])
    for metal, mult in (("La", 1), ("Dy", 6)):
        key = f"{metal}_{cn}"
        if key in res: continue
        src = d / f"La_{cn}_mace.xyz" if metal == "Dy" else f
        if not src.exists(): continue
        a = read(src)
        if metal == "Dy":
            a.symbols[[k for k, s in enumerate(a.get_chemical_symbols()) if s == "La"][0]] = "Dy"
        a.info = {"charge": q, "spin": mult}
        a.calc = calc
        opt = BFGS(a, logfile=None, maxstep=0.1)
        conv = opt.run(fmax=0.03, steps=600)
        a.write(d / f"{metal}_{cn}_mace.xyz", format="xyz", comment=f"charge={q} mult={mult}")
        res[key] = {"E_kcal": a.get_potential_energy() * EV, "charge": q, "converged": bool(conv), "steps": opt.nsteps, **coord(a)}
        (d / "energies.json").write_text(json.dumps(res, indent=1))
        print(key, res[key], flush=True)
