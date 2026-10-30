"""Lanthanide-series scan: relax M.L for M in Nd, Gd, Ho, Lu (La, Dy already in energies.json) by metal
swap of the MACE-relaxed La structure, same CN. Target is the *shape* of the series (Gd break / PC2): a
within-ligand second difference in which the charge-linear gas artefact largely cancels.

usage: relax_series.py <dir> [threads]  -> <dir>/series.json {"<M>_cn<N>": {E_kcal, converged, CN, ...}}
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
MULT = {"Nd": 4, "Gd": 8, "Ho": 5, "Lu": 1}
d = Path(sys.argv[1]); torch.set_num_threads(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
calc = MACECalculator(model_paths=CKPT, device="cpu", default_dtype="float64", head="omol")


def coord(a, mi):
    r = np.linalg.norm(a.positions - a.positions[mi], axis=1)
    sym = a.get_chemical_symbols()
    don = [k for k in range(len(a)) if sym[k] in ("O", "N") and r[k] < 3.1]
    return {"CN": len(don), "M_donor_mean_A": float(r[don].mean()) if don else None}


res = json.loads((d / "series.json").read_text()) if (d / "series.json").exists() else {}
for f in sorted(d.glob("La_cn?_mace.xyz")):
    cn = f.stem.split("_")[1]
    q = int(open(f).readlines()[1].split("charge=")[1].split()[0])
    for metal, mult in MULT.items():
        key = f"{metal}_{cn}"
        if key in res: continue
        a = read(f)
        mi = [k for k, s in enumerate(a.get_chemical_symbols()) if s == "La"][0]
        a.symbols[mi] = metal
        a.info = {"charge": q, "spin": mult}
        a.calc = calc
        opt = BFGS(a, logfile=None, maxstep=0.1)
        conv = opt.run(fmax=0.03, steps=600)
        a.write(d / f"{metal}_{cn}_mace.xyz", format="xyz", comment=f"charge={q} mult={mult}")
        res[key] = {"E_kcal": a.get_potential_energy() * EV, "charge": q, "converged": bool(conv), "steps": opt.nsteps, **coord(a, mi)}
        (d / "series.json").write_text(json.dumps(res, indent=1))
        print(key, res[key], flush=True)
