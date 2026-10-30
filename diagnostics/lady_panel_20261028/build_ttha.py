"""TTHA (NIST 678) La/Dy states for the panel, added 2026-10-29 as the test between the two bias forms (README).

TTHA is decadentate and leaves no water on either metal (q = 0), so it displaces the most water of any NIST ligand
(n_disp 8.5) while being nearly unselective (ddG Dy-La -0.68 kcal/mol). Starts are the MACE "cn8" builds from
ln_surrogate_20261027/nist/678: 9 ligand donors within 3.0 A, no water. DFT may close the 10th donor on La.
Method block byte-identical to the panel (make_inputs.py). Writes start/ and runs/ under workspaces/lady_panel_20261028.
"""
import json, sys
from pathlib import Path

D = Path(__file__).resolve().parent
REPO = D.parents[1]
sys.path.insert(0, str(D)); sys.path.insert(0, str(REPO / "diagnostics/lady_relaxed_exchange_20261024"))
import build_panel as B  # noqa: E402
from build import write_xyz  # noqa: E402

OLD = REPO / "diagnostics/lady_relaxed_exchange_20261024"
W = REPO / "workspaces/lady_panel_20261028"
HEAD = (OLD / "make_inputs.py").read_text().split('HEAD = """', 1)[1].split('"""', 1)[0]
Q = -3

meta = json.loads((W / "start/SYSTEMS.json").read_text())
for m in ("La", "Dy"):
    name = f"{m}_ttha_w0"
    el, X = B.load(B.NIST / "678" / f"{m}_cn8_mace.xyz", m)
    s = B.summary(el, X)
    ok = s["n_water_total"] == 0 and s["n_ligand_donors"] >= 9
    print(f"{name:14s} q={Q:+d} natoms={len(el)} {'OK ' if ok else 'BAD'} {s}")
    write_xyz(W / "start" / f"{name}.xyz", el, X, f"charge {Q}")
    meta[name] = {"charge": Q, "natoms": len(el), "nist_id": "678", "source": f"{m}_cn8_mace.xyz", "edit": None,
                  "start_ok": ok, **s}
    run = W / "runs" / name; run.mkdir(parents=True, exist_ok=True)
    (run / "start.xyz").write_text((W / "start" / f"{name}.xyz").read_text())
    (run / "job.inp").write_text(HEAD + (OLD / f"basis_{m}.inc").read_text() + f"* xyzfile {Q} 1 start.xyz\n")
(W / "start/SYSTEMS.json").write_text(json.dumps(meta, indent=1))
