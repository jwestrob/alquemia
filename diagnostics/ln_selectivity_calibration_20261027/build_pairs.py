"""Per ligand: log K1(Dy) - log K1(La) at 25 C, same ionic strength (prefer I=0.1, else the
I with both values closest to 0.1), plus dH/dS for both metals at that I when present.
ddG(Dy-La) = -RT ln10 * dlogK (kcal/mol, negative = Dy favoured). Writes pairs.tsv."""
import csv
from collections import defaultdict
RTLN10 = 1.36394  # kcal/mol at 298.15 K
d = defaultdict(lambda: defaultdict(dict)); meta = {}
for r in csv.DictReader(open("ln_ML.tsv"), delimiter="\t"):
    if r["T"] != "25": continue
    d[r["lig_id"]][(r["qty"], float(r["I"]))][r["metal"]] = float(r["value"])
    meta[r["lig_id"]] = (r["formula"], r["prot"], r["name"], r["tags"])
rows = []
for L, q in d.items():
    Is = [I for (k, I), m in q.items() if k == "K" and "La" in m and "Dy" in m]
    if not Is: continue
    I = min(Is, key=lambda x: abs(x - 0.1))
    k = q[("K", I)]
    h = q.get(("H", I), {}); s = q.get(("S", I), {})
    dH = h["Dy"] - h["La"] if "La" in h and "Dy" in h else None
    dS = s["Dy"] - s["La"] if "La" in s and "Dy" in s else None  # cal/mol/K
    rows.append(dict(lig_id=L, formula=meta[L][0], prot=meta[L][1], name=meta[L][2], tags=meta[L][3], I=I,
                     logK_La=k["La"], logK_Dy=k["Dy"], dlogK=round(k["Dy"] - k["La"], 3),
                     ddG=round(-RTLN10 * (k["Dy"] - k["La"]), 3), ddH=dH, TddS=None if dS is None else round(298.15 * dS / 1000, 3),
                     logK_series=";".join(f"{m}:{k[m]}" for m in ["La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Y"] if m in k)))
rows.sort(key=lambda r: r["dlogK"])
w = csv.DictWriter(open("pairs.tsv", "w"), fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(len(rows), "ligands;", sum(r["ddH"] is not None for r in rows), "with dH for both")
