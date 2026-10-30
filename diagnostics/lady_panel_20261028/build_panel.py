"""Starting geometries for the constant-bias panel (same protocol as lady_relaxed_exchange_20261024).

Source: MACE-OMol gas-relaxed Architector complexes from ln_surrogate_20261027/nist/<id>/ (La and Dy
relaxed separately; DFT-CPCM re-optimises everything). Hydration states:
  malonate (3480, complex +1): w6, w7   (w6 = w7 minus its most distant first-shell water)
  NTA      (465,  complex  0): w4, w5
  DTPA     (655,  complex -2): w1       (q = 1 across the series)
  CDTA     (600,  complex -1): w2, w3   (w3 = w2 + one water in the most open direction)
  TMDTA    (610,  complex -1): w2, w3
Metal is written first (analyze.py measures donor distances from atom 0).
Writes workspaces/lady_panel_20261028/start/{<name>.xyz, SYSTEMS.json}.
"""
import json, sys
from pathlib import Path
import numpy as np
from ase.io import read

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "diagnostics/lady_relaxed_exchange_20261024"))
from build import add_water, write_xyz  # noqa: E402  (same water placement as the acetate/EDTA ladder)

NIST = REPO / "workspaces/ln_surrogate_20261027/nist"
OUT = REPO / "workspaces/lady_panel_20261028/start"
# lig: (nist id, complex charge, {state: (source cn file, edit)})
PANEL = {
    "mal":   ("3480", 1,  {7: ("cn9", None), 6: ("cn9", "drop")}),
    "nta":   ("465",  0,  {4: ("cn8", None), 5: ("cn9", None)}),
    "dtpa":  ("655",  -2, {1: ("cn9", None)}),
    "cdta":  ("600",  -1, {2: ("cn8", None), 3: ("cn8", "add")}),
    "tmdta": ("610",  -1, {2: ("cn8", None), 3: ("cn8", "add")}),
}
EXPECT_DONORS = {"mal": 2, "nta": 4, "dtpa": 8, "cdta": 6, "tmdta": 6}


def load(path, metal):
    a = read(path); el = a.get_chemical_symbols(); X = a.positions.copy()
    mi = el.index(metal); order = [mi] + [i for i in range(len(el)) if i != mi]
    return [el[i] for i in order], X[order]


def waters(el, X):
    """First-shell and other water molecules: O bonded to exactly two H and nothing else (metal excluded)."""
    out = []
    for i, e in enumerate(el):
        if e != "O" or i == 0: continue
        nb = [j for j in range(1, len(el)) if j != i and np.linalg.norm(X[j] - X[i]) < 1.25]
        if len(nb) == 2 and all(el[j] == "H" for j in nb): out.append([i] + nb)
    return out


def summary(el, X):
    r = np.linalg.norm(X - X[0], axis=1)
    W = waters(el, X); wO = {w[0] for w in W}
    don = [i for i in range(1, len(el)) if el[i] in ("O", "N") and r[i] < 3.0]
    return {"n_water_total": len(W), "n_water_bound": sum(1 for i in don if i in wO),
            "n_ligand_donors": sum(1 for i in don if i not in wO),
            "donor_A": sorted(round(float(r[i]), 2) for i in don)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = {}
    for lig, (nid, q, states) in PANEL.items():
        for m in ("La", "Dy"):
            for nw, (src, edit) in states.items():
                el, X = load(NIST / nid / f"{m}_{src}_mace.xyz", m)
                if edit == "drop":
                    r = np.linalg.norm(X - X[0], axis=1)
                    w = max(waters(el, X), key=lambda w: r[w[0]])
                    keep = [i for i in range(len(el)) if i not in w]
                    el, X = [el[i] for i in keep], X[keep]
                elif edit == "add":
                    el, X = add_water(el, X, m)
                s = summary(el, X); name = f"{m}_{lig}_w{nw}"
                ok = s["n_water_total"] == nw == s["n_water_bound"] and s["n_ligand_donors"] == EXPECT_DONORS[lig]
                print(f"{name:14s} q={q:+d} natoms={len(el)} {'OK ' if ok else 'BAD'} {s}")
                write_xyz(OUT / f"{name}.xyz", el, X, f"charge {q}")
                meta[name] = {"charge": q, "natoms": len(el), "nist_id": nid, "source": f"{m}_{src}_mace.xyz",
                              "edit": edit, "start_ok": ok, **s}
    (OUT / "SYSTEMS.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
