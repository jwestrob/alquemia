#!/usr/bin/env python3
"""Transcribe Tables S-2 to S-5, S-7 and S-8 of Wang et al. 2020 (Geochem. Persp. Lett. 14, 26-30, Supplementary
Information; open access) from the PDF text layer. Untrusted data, numbers only; run with `python3 -I`.
Usage: python3 -I transcribe_mussels.py SI.txt OUT.tsv"""
import csv, re, sys
src, out = sys.argv[1:3]
EL = ["Y", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Yb", "Lu"]
SPEC = {"S-2": ("Gigantidas platifrons", "Site F", "methanotrophs"), "S-3": ("Bathymodiolus aduloides", "Site F", "thiotrophs"),
        "S-4": ("Gigantidas haimaensis", "Haima", "methanotrophs"), "S-5": ("Calyptogena marissinica", "Haima", "thiotrophs"),
        "S-7": ("Mytilus edulis", "France coast", "none"), "S-8": ("sediment", None, None)}
ORG = {"shells": "shell", "gills": "gill", "gill": "gill", "mantles": "mantle", "mantle": "mantle", "feet": "foot", "haslets": "haslet"}
rows, table, organ, place = [], None, None, None
for line in open(src, encoding="utf8", errors="ignore"):
    t = line.strip()
    m = re.match(r"^Table (S-\d)", t)
    if m:
        table, organ, place = m.group(1), None, None
        continue
    if table not in SPEC:
        continue
    low = t.lower()
    if low in ORG:
        organ = ORG[low]; continue
    if table == "S-7" and low in ("ouessant", "chausey", "gouville/mer"):
        place = t; continue
    if table == "S-8" and ("sediments" in low or low.startswith("bcr2")):
        place = t.split(",")[0]; continue
    tok = re.split(r"\s{2,}", t)
    if len(tok) < 15:
        continue
    label, vals = tok[0], tok[1:]
    def num(x):
        try:
            return float(x)
        except ValueError:
            return None
    v = [num(x) for x in vals]
    if len(v) not in (16,) or sum(x is not None for x in v[:14]) < 12:
        continue
    if table == "S-8":
        if place is None or place.startswith("BCR2"):
            continue
        rows.append(dict(table=table, species="sediment", site=place, symbionts="", animal=label, organ="sediment", unit="ug/g", **dict(zip(EL, v[:14]))))
        continue
    mm = re.match(r"^([A-Za-z]+\d+),\s*(\w+)$", label)
    if not mm:
        continue
    sp, site, sym = SPEC[table]
    org = organ if table != "S-7" else {"S": "shell", "ST": "soft tissue"}[mm.group(2)]
    rows.append(dict(table=table, species=sp, site=(place if table == "S-7" else site), symbionts=sym, animal=mm.group(1), organ=org, unit="ng/g", **dict(zip(EL, v[:14]))))
with open(out, "w", newline="") as fh:
    fh.write("# Wang, Barrat, Bayon, Chauvaud & Feng 2020, Geochem. Persp. Lett. 14:26-30, Supplementary Tables S-2..S-5, S-7, S-8. Tm was not reported; blank = not analysed.\n")
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
from collections import Counter
for k, n in sorted(Counter((r["table"], r["species"], r["site"], r["organ"]) for r in rows).items()):
    print(n, k)
print("rows with a missing element:", [(r["animal"], r["organ"], [e for e in EL if r[e] is None]) for r in rows if any(r[e] is None for e in EL)])
