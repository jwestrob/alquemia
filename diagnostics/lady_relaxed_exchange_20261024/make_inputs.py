"""ORCA PBE0-D4/def2-TZVP(+archived La/Dy ECP basis)/CPCM(water) Opt+Freq inputs."""
import json, shutil, sys
from pathlib import Path
D = Path(__file__).resolve().parent
W = D.parents[1] / "workspaces/lady_relaxed_exchange_20261024"
HEAD = """! PBE0 D4 def2-TZVP def2/J RIJCOSX DefGrid3 TightSCF TightOpt Freq CPCM(Water)
%maxcore 3000
%pal nprocs @NPROCS@ end
%scf
 Guess PModel
 MaxIter 500
end
%geom
 MaxIter 300
end
"""
meta = json.loads((W / "start/SYSTEMS.json").read_text())
for name, m in meta.items():
    run = W / "runs" / name; run.mkdir(parents=True, exist_ok=True)
    shutil.copy(W / "start" / f"{name}.xyz", run / "start.xyz")
    metal = name.split("_")[0]
    basis = (D / f"basis_{metal}.inc").read_text() if metal in ("La", "Dy") else ""
    (run / "job.inp").write_text(HEAD + basis + f"* xyzfile {m['charge']} 1 start.xyz\n")
    print(run)
