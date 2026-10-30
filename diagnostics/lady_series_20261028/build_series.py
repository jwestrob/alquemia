"""Build Nd and Lu Opt+Freq inputs for the series extension of the La/Dy panel.

Each Nd state starts from its La analogue and each Lu state from its Dy analogue (same ligand, same water count),
with the metal swapped. Geometry priority: conformer-best DFT structure (aqua/acetate/EDTA), else the finished panel
geometry if it kept its coordination number, else the panel MACE start. Method block is byte-identical to the La/Dy
runs; only the metal's basis/ECP/AuxJ include differs (basis_Nd.inc, basis_Lu.inc: lcecp-1-TZVP, pinned AutoAux J).
Run this after the La/Dy panel has finished so the panel geometries are final.

2026-10-29 (before any Nd/Lu job ran): TTHA (ttha_w0) added; like the panel, it is exempt from the CN-match test
(no water to lose; 9 vs 10 bound donors are isomers of the same q = 0 state). If its La/Dy analogue did not finish,
the start is the metal-specific MACE build (ln_surrogate_20261027/nist/678/<M>_cn8_mace.xyz). Dropped, because
their La/Dy analogue was never run (tier 2) or lost a water in DFT, so the smaller Nd/Lu ion has no state to test:
Nd_tmdta_w3 (La_tmdta_w3 went to CN 8), Lu_cdta_w3, Lu_tmdta_w3, Lu_mal_w7 (Dy used w2/w2/w6 only).
"""
import re, sys
from pathlib import Path

D = Path(__file__).resolve().parent
REPO = D.parents[1]
OLD = REPO / "diagnostics/lady_relaxed_exchange_20261024"
sys.path.insert(0, str(OLD))
import analyze as A  # noqa: E402

PANEL = REPO / "workspaces/lady_panel_20261028/runs"
OUT = REPO / "workspaces/lady_series_20261028/runs"
HEAD = (OLD / "make_inputs.py").read_text().split('HEAD = """', 1)[1].split('"""', 1)[0]
NIST = REPO / "workspaces/ln_surrogate_20261027/nist"
N_LIG_DONORS = {"mal": 2, "nta": 4, "dtpa": 8, "cdta": 6, "tmdta": 6, "ttha": 9}
NO_CN_CHECK = {"ttha"}
PANEL_STATES = ["mal_w6", "mal_w7", "nta_w4", "nta_w5", "dtpa_w1", "cdta_w2", "cdta_w3", "tmdta_w2", "tmdta_w3", "ttha_w0"]
DROP = {"Nd": {"tmdta_w3"}, "Lu": {"cdta_w3", "tmdta_w3", "mal_w7"}}
TTHA_MACE = {"ttha_w0": "678"}  # fallback start: metal-specific MACE build if the La/Dy analogue is unfinished
TARGETS = {
    "Nd": ("La", ["aq8", "aq9", "acmono_w7", "acmono_w8", "acbi_w6", "acbi_w7", "acbi_w8", "edtax_w2", "edtax_w3"] + PANEL_STATES),
    "Lu": ("Dy", ["aq8", "aq9", "acmono_w7", "acmono_w8", "acbi_w6", "acbi_w7", "edtax_w2", "edtax_w3"] + PANEL_STATES),
}
TARGETS = {m: (src, [s for s in sts if s not in DROP[m]]) for m, (src, sts) in TARGETS.items()}


def charge_of(inp):
    return int(re.search(r"^\*\s*xyzfile\s+(-?\d+)\s+1", inp.read_text(), re.M).group(1))


def source(src_name, metal):
    """Return (xyz path, charge, provenance) for the La/Dy analogue."""
    old = A.W / src_name
    if old.exists():
        r = A.parse(src_name)
        if r and r["normal"]:
            _, picked = A.conformer_best({src_name: r})
            ch = picked[src_name]["chosen"]
            xyz = old / "job.xyz" if ch == "primary" else A.W.parent / ch.split("/")[0] / src_name / ch.split("/")[1] / "job.xyz"
            return xyz, charge_of(old / "job.inp"), f"old run, conformer-best ({ch})"
        return old / "job.xyz", charge_of(old / "job.inp"), "old run, NOT normal-terminated (last geometry)"
    p = PANEL / src_name
    lig, nw = src_name.split("_")[1], A.nw(src_name)
    r = A.parse(src_name, p)
    if r and r["normal"] and (lig in NO_CN_CHECK or sum(1 for _, x in r.get("donors", []) if x < 3.1) == N_LIG_DONORS[lig] + nw):
        return p / "job.xyz", charge_of(p / "job.inp"), "panel run, final geometry"
    st = src_name.split("_", 1)[1]
    if st in TTHA_MACE:
        return NIST / TTHA_MACE[st] / f"{metal}_cn8_mace.xyz", charge_of(p / "job.inp"), "metal-specific MACE build (analogue unfinished)"
    return p / "start.xyz", charge_of(p / "job.inp"), "panel MACE start (analogue unfinished, failed or CN changed)"


def main():
    n = 0
    for metal, sts in DROP.items():  # stale dry-run inputs of dropped states would otherwise be submitted
        for st in sts:
            run = OUT / f"{metal}_{st}"
            if run.exists() and not (run / "job.out").exists():
                (OUT.parent / "dropped").mkdir(exist_ok=True); run.rename(OUT.parent / "dropped" / run.name)
                print(f"moved unrun dropped state {run.name} to dropped/")
    for metal, (src_metal, states) in TARGETS.items():
        for st in states:
            name, src_name = f"{metal}_{st}", f"{src_metal}_{st}"
            xyz, q, prov = source(src_name, metal)
            lines = xyz.read_text().splitlines()
            nat = int(lines[0]); atoms = lines[2:2 + nat]
            assert atoms[0].split()[0] in (src_metal, metal), (src_name, atoms[0])
            atoms[0] = atoms[0].replace(src_metal, metal, 1)
            run = OUT / name; run.mkdir(parents=True, exist_ok=True)
            (run / "start.xyz").write_text(f"{nat}\n{name} from {src_name} ({prov}), metal swapped\n" + "\n".join(atoms) + "\n")
            (run / "job.inp").write_text(HEAD + (D / f"basis_{metal}.inc").read_text() + f"* xyzfile {q} 1 start.xyz\n")
            (run / "SOURCE.txt").write_text(f"{src_name}\t{xyz}\t{prov}\tcharge {q}\n")
            print(f"{name:16s} q={q:+d}  <- {src_name:15s} {prov}")
            n += 1
    print(n, "inputs written to", OUT)


if __name__ == "__main__":
    main()
