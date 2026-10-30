"""PBE0-D4/def2-TZVP (archived large-core ECP for La, Dy; 4f in core, closed shell), CPCM(water) single
points on the MACE CN8 La and Dy structures of every NIST ligand plus the aqua ion.
Writes <ws>/cpcm/<lig>_<M>/{job.inp,start.xyz} and <ws>/cpcm/tasks.txt"""
import json, sys
from pathlib import Path
W = Path(sys.argv[1]); R = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs")
BAS = {m: (R / f"diagnostics/lady_relaxed_exchange_20261024/basis_{m}.inc").read_text().rstrip() for m in ("La", "Dy")}
HEAD = "! PBE0 D4 def2-TZVP def2/J RIJCOSX DefGrid3 TightSCF CPCM(Water)"
tasks = []
for d in sorted((W / "nist").iterdir()):
    e = json.loads((d / "energies.json").read_text()) if (d / "energies.json").exists() else {}
    for M in ("La", "Dy"):
        x = d / f"{M}_cn8_mace.xyz"
        if f"{M}_cn8" not in e or not x.exists(): continue
        o = W / "cpcm" / f"{d.name}_{M}"; o.mkdir(parents=True, exist_ok=True)
        lines = x.read_text().splitlines(); lines[1] = ""
        (o / "start.xyz").write_text("\n".join(lines) + "\n")
        (o / "job.inp").write_text("\n".join([HEAD, "%maxcore 3000", "%pal nprocs @NPROCS@ end", "%scf\n MaxIter 1000\nend",
                                              BAS[M], f"* xyzfile {e[f'{M}_cn8']['charge']} 1 start.xyz\n"]))
        tasks.append((len(lines), o.name))
tasks.sort(reverse=True)  # largest first
(W / "cpcm" / "tasks.txt").write_text("\n".join(t[1] for t in tasks) + "\n")
print(len(tasks), "inputs")
