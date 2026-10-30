"""Ingredient tests: which part of the model carries the La/Dy exchange error?

Fixed geometries = conformer-best PBE0-D4/CPCM structures from
lady_relaxed_exchange_20261024 (analyze.py --conformers). Single points only:

  gas      PBE0-D4, archived large-core ECP46/55 basis, no CPCM   -> CPCM contribution
  sc_cpcm  PBE0-D4, ORCA def2-TZVP small-core def2-ECP (4f explicit; Dy UKS S=5/2), CPCM
  sc_gas   same, no CPCM                                           -> 4f-in-core effect
  dlpno    DLPNO-CCSD(T) TightPNO, large-core basis, no CPCM       -> functional error vs gas

The reference ("lc_cpcm") is the existing Opt+Freq final energy at the same geometry.
Bare La3+/Dy3+ (gas, CPCM, both cores) give the differential hydration check.
"""
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "workspaces/lady_relaxed_exchange_20261024"
OUT = REPO / "workspaces/lady_ingredient_tests_20261026"
BASIS = {m: (REPO / f"diagnostics/lady_relaxed_exchange_20261024/basis_{m}.inc").read_text() for m in ("La", "Dy")}
STATES = ["La_aq8", "La_aq9", "Dy_aq8", "Dy_aq9", "La_acmono_w7", "La_acmono_w8", "Dy_acmono_w7", "Dy_acmono_w8",
          "La_acbi_w6", "La_acbi_w7", "Dy_acbi_w6", "Dy_acbi_w7"]
HEAD = {
    "gas": "! PBE0 D4 def2-TZVP def2/J RIJCOSX DefGrid3 TightSCF",
    "sc_cpcm": "! UKS PBE0 D4 def2-TZVP def2/J RIJCOSX DefGrid3 TightSCF SlowConv CPCM(Water)",
    "sc_gas": "! UKS PBE0 D4 def2-TZVP def2/J RIJCOSX DefGrid3 TightSCF SlowConv",
    "dlpno": "! DLPNO-CCSD(T) TightPNO def2-TZVP AutoAux TightSCF",
    "lc_cpcm_ion": "! PBE0 D4 def2-TZVP def2/J RIJCOSX DefGrid3 TightSCF CPCM(Water)",
}


def geom(state):
    sel = json.loads((SRC / "EXCHANGE_RESULT_CONFORMERS.json").read_text())["conformer_selection"][state]["chosen"]
    return SRC / "runs" / state / "job.xyz" if sel == "primary" else SRC / sel.split("/")[0] / state / sel.split("/")[1] / "job.xyz"


def write(name, test, xyz_text, charge, metal):
    small = test.startswith("sc")
    mult = 6 if (small and metal == "Dy") else 1
    d = OUT / test / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "start.xyz").write_text(xyz_text)
    body = [HEAD[test], "%maxcore 3000", "%pal nprocs @NPROCS@ end", "%scf\n MaxIter 1000\nend"]
    if metal and not small:
        body.append(BASIS[metal].rstrip())
    body.append(f"* xyzfile {charge} {mult} start.xyz\n")
    (d / "job.inp").write_text("\n".join(body))
    return f"{test}/{name}"


def main():
    runs = []
    for s in STATES:
        m = s.split("_")[0]
        q = 3 if "_aq" in s else 2
        x = geom(s).read_text()
        for t in ("gas", "sc_cpcm", "sc_gas", "dlpno"):
            runs.append(write(s, t, x, q, m))
    w = (SRC / "runs/water/job.xyz").read_text()
    for t in ("gas", "dlpno"):
        runs.append(write("water", t, w, 0, None))
    for m in ("La", "Dy"):
        ion = f"1\nbare {m}3+\n{m} 0.0 0.0 0.0\n"
        for t in ("gas", "sc_gas", "sc_cpcm", "lc_cpcm_ion"):
            runs.append(write(f"{m}_ion", t, ion, 3, m))
    (OUT / "RUNS.txt").write_text("\n".join(runs) + "\n")
    print(len(runs), "inputs")


if __name__ == "__main__":
    main()
