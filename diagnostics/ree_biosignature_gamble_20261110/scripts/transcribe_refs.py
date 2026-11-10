#!/usr/bin/env python3
"""Transcribe published tables (text layer of open-access PDFs; untrusted data, numbers only) into tidy files.
Run with `python3 -I`. Sources:
  Kawabe et al. 1999, Geochem. J. 33, 167-179 (Tables 1-a,b,c): Kd between Fe(-Mn) oxyhydroxide and NaCl solution.
  Kawabe et al. 1999, Geochem. J. 33, 181-197 (Table 1): the same with NaHCO3 in solution.
  Ohta & Kawabe 2000, Geochem. J. 34, 455-473 (Table 1): infinite-dilution complexation constants at 25 C.
Usage: python3 -I transcribe_refs.py LIT_DIR OUT_DIR"""
import csv, os, re, sys
lit, out = sys.argv[1:3]
os.makedirs(out, exist_ok=True)
EL = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y"]
NUM = re.compile(r"\d+\.\d+|\d+")


def block(path, start_pat, n_cols, end_pat=r"^\s*S[ceC]\s"):
    """Rows 'El v1 v2 ...' after the line matching start_pat, until the Sc row."""
    lines = open(path, encoding="utf8", errors="ignore").read().splitlines()
    i0 = next(i for i, l in enumerate(lines) if re.search(start_pat, l))
    tab, meta = {}, {}
    for l in lines[i0:i0 + 60]:
        t = l.strip().replace("'", " ")
        m = re.match(r"^(La|Ce|Pr|Nd|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Y)\s+(.*)$", t)
        if m and m.group(1) not in tab:
            v = [float(x) for x in NUM.findall(m.group(2))]
            assert len(v) == n_cols, (path, m.group(1), v)
            tab[m.group(1)] = v
        for key, pat in (("pH", r"^pH\*+\s+(.*)$"), ("NaCl_M", r"^NaC[lI1][\s(].*?M\)\.?\s+(.*)$"), ("NaHCO3_M", r"^NaHCO3 \(M\)\s+(.*)$"), ("T_C", r"^T \(.C\)\*+\s+(.*)$")):
            mm = re.match(pat, t)
            if mm and key not in meta:
                meta[key] = [float(x) for x in NUM.findall(mm.group(1))]
    assert set(tab) == set(EL), (path, sorted(set(EL) - set(tab)))
    return tab, meta


a = os.path.join(lit, "kawabe1999a/kawabe1999_gj33_167.txt"); b = os.path.join(lit, "kawabe1999b/kawabe1999_gj33_181.txt")
runs = []
for start, n, names in [(r"Table 1-a\. Experimental", 6, ["MnFe-1", "MnFe-2", "Fe1-1", "Fe1-2", "Fe1-3", "Fe1-4"]),
                        (r"Table 1-b\. Experimental", 6, [f"Fe2-{i}" for i in range(1, 7)]),
                        (r"Table 1-c\. Experimental", 4, [f"Fe2-{i}" for i in range(7, 11)])]:
    tab, meta = block(a, start, n)
    for j, nm in enumerate(names):
        runs.append(dict(run=nm, source="Kawabe1999a", solid=("Fe-Mn oxyhydroxide" if nm.startswith("Mn") else "Fe oxyhydroxide"),
                         pH=meta["pH"][j], NaCl_M=meta["NaCl_M"][j], NaHCO3_M=0.0, **{e: tab[e][j] for e in EL}))
tab, meta = block(b, r"Table 1\. Experimental", 6)
for j in range(1, 6):                                            # run 1 repeats Fe1-1 of the first paper
    runs.append(dict(run=f"FeCO3-{j + 1}", source="Kawabe1999b", solid="Fe oxyhydroxide", pH=meta["pH"][j], NaCl_M=meta["NaCl_M"][j],
                     NaHCO3_M=meta["NaHCO3_M"][j], **{e: tab[e][j] for e in EL}))
assert all(abs(tab[e][0] - r) < 1e-9 for e, r in zip(EL, [next(x for x in runs if x["run"] == "Fe1-1")[e] for e in EL])), "run 1 of 1999b should equal Fe1-1"
with open(os.path.join(out, "kawabe1999_kd_e4.tsv"), "w", newline="") as fh:
    fh.write("# Kd/1e4 (mole fraction in precipitate per mol/L in solution). Kawabe et al. 1999 Geochem. J. 33:167 (Tables 1-a,b,c) and 33:181 (Table 1).\n")
    w = csv.DictWriter(fh, fieldnames=list(runs[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(runs)
print("oxide runs:", len(runs), [r["run"] for r in runs])
for r in runs:
    print(f"  {r['run']:8s} pH {r['pH']:.2f} NaCl {r['NaCl_M']:.3f} NaHCO3 {r['NaHCO3_M']:.4f} | " + " ".join(f"{r[e]:g}" for e in EL))

# complexation constants
c = os.path.join(lit, "ohta2000b/ohta2000_gj34_455.txt")
lines = open(c, encoding="utf8", errors="ignore").read().splitlines()
i0 = next(i for i, l in enumerate(lines) if "Table 1. Infinite dilution stability constants" in l)
cols = ["b2_CO3_OK2000", "b1_CO3_OK2000", "b2_CO3_LB1998", "b1_CO3_LB1998", "b1_HCO3_LB1998", "b1_OH_M1992", "b1_SO4_M1992", "b1_Cl_Mir1982"]
K = []
for l in lines[i0:i0 + 30]:
    m = re.match(r"^\s*(La|Ce|Pr|Nd|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Y)\s+(.*)$", l)
    if m:
        v = [float(x) for x in NUM.findall(m.group(2))]
        assert len(v) == 8, (m.group(1), v)
        K.append(dict(element=m.group(1), **dict(zip(cols, v))))
assert [k["element"] for k in K] == EL
with open(os.path.join(out, "ree_complexation_constants.tsv"), "w", newline="") as fh:
    fh.write("# log10 K at infinite dilution, 25 C. Table 1 of Ohta & Kawabe 2000 Geochem. J. 34:455. OK2000 = Ohta & Kawabe 2000a; LB1998 = Liu & Byrne 1998; M1992 = Millero 1992; Mir1982 = Mironov et al. 1982. Y values for OH, SO4, Cl were calculated by the authors after Byrne & Lee 1993.\n")
    w = csv.DictWriter(fh, fieldnames=["element"] + cols, delimiter="\t"); w.writeheader(); w.writerows(K)
for k in K:
    print(" ", k)
