"""Opt+Freq inputs for the panel: byte-identical method block to lady_relaxed_exchange_20261024/make_inputs.py
(PBE0-D4/def2-TZVP, archived large-core La/Dy basis+ECP, CPCM(water), TightOpt Freq, PModel guess)."""
import json, shutil
from pathlib import Path
D = Path(__file__).resolve().parent
OLD = D.parent / "lady_relaxed_exchange_20261024"
W = D.parents[1] / "workspaces/lady_panel_20261028"
src = (OLD / "make_inputs.py").read_text()
HEAD = src.split('HEAD = """', 1)[1].split('"""', 1)[0]
meta = json.loads((W / "start/SYSTEMS.json").read_text())
for name, m in meta.items():
    run = W / "runs" / name; run.mkdir(parents=True, exist_ok=True)
    shutil.copy(W / "start" / f"{name}.xyz", run / "start.xyz")
    metal = name.split("_")[0]
    (run / "job.inp").write_text(HEAD + (OLD / f"basis_{metal}.inc").read_text() + f"* xyzfile {m['charge']} 1 start.xyz\n")
print(len(meta), "inputs;", HEAD.splitlines()[0])
