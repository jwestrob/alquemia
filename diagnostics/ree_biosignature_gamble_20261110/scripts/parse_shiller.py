#!/usr/bin/env python3
"""Parse the two public Deepwater Horizon workbooks (Shiller; GRIIDC) into tidy tables.
Downloaded files are untrusted data: run with `python3 -I`, values only. Lead agent, Claude Fable, 2026-11-10.
Usage: python3 -I parse_shiller.py REE.xlsx TRACE.xlsx OUT_DIR"""
import csv, datetime, os, re, sys
import openpyxl

ree_path, tm_path, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
EL = ["Y", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]


def num(x):
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip()
    try:
        return float(s)
    except ValueError:
        return None


def norm_station(s):
    return re.sub(r"[\s\-_]", "", str(s)).upper() if s is not None else None


# ---- rare earth sheet
ws = openpyxl.load_workbook(ree_path, data_only=True).worksheets[0]
rows = list(ws.iter_rows(values_only=True))
hdr = [str(c).strip() if c is not None else "" for c in rows[1]]
col = {e: hdr.index(e) for e in EL}
c_station, c_date, c_lat, c_lon, c_depth = 1, 2, 3, 4, 5
c_id = next(i for i, h in enumerate(hdr) if h.startswith("Sample ID"))
print("REE header:", [h for h in hdr if h][:30])
print("non-empty cells right of the sample ID column:", sum(1 for r in rows[2:] for c in r[c_id + 1:] if c is not None))
ree, cruise, last_station, last_date, notes = [], None, None, None, []
for i, r in enumerate(rows[3:], start=4):
    if r[0] is not None and str(r[0]).strip():
        txt = str(r[0]).strip()
        if all(r[col[e]] is None for e in EL):
            cruise = txt; last_station = None; last_date = None
            print(f"row {i}: cruise header -> {txt!r}")
            continue
        notes.append((i, txt))
    vals = {e: num(r[col[e]]) for e in EL}
    if all(v is None for v in vals.values()):
        if any(c is not None for c in r):
            notes.append((i, " | ".join(str(c) for c in r if c is not None)[:200]))
        continue
    if r[c_station] is not None:
        last_station = str(r[c_station]).strip()
    if r[c_date] is not None:
        last_date = r[c_date]
    d = last_date.date().isoformat() if isinstance(last_date, datetime.datetime) else str(last_date)
    raw = {e: r[col[e]] for e in EL}
    flag = ";".join(f"{e}={raw[e]!r}" for e in EL if raw[e] is not None and vals[e] is None)
    ree.append(dict(row=i, cruise=cruise, station=last_station, station_key=norm_station(last_station), date=d,
                    lat=num(r[c_lat]), lon=num(r[c_lon]), depth=num(r[c_depth]), sample_id=r[c_id], flag=flag, **vals))
print("REE samples:", len(ree))
for n in notes:
    print("note row", n)
from collections import Counter
print(Counter(x["cruise"] for x in ree))
with open(os.path.join(out, "dwh_ree.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(ree[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(ree)

# ---- companion sheets
wb = openpyxl.load_workbook(tm_path, data_only=True)
tm = []
for ws in wb.worksheets[:4]:
    rows = list(ws.iter_rows(values_only=True))
    hdr = [re.sub(r"\s+", " ", str(c)).strip() if c is not None else "" for c in rows[1]]
    def find(name):
        for j, h in enumerate(hdr):
            if h.lower().startswith(name.lower()):
                return j
        return None
    idx = {k: find(k) for k in ["Station", "Date", "Latitude", "Longitude", "Bottom", "Sample", "Temp", "Sal", "O2", "NO3", "PO4", "Si", "Ba", "Fe", "Mn", "Total PAH"]}
    print(ws.title, {k: v for k, v in idx.items()})
    last_station = last_date = None
    for i, r in enumerate(rows[2:], start=3):
        if r[idx["Sample"]] is None or num(r[idx["Sample"]]) is None:
            continue
        if r[idx["Station"]] is not None:
            last_station = str(r[idx["Station"]]).strip()
        if r[idx["Date"]] is not None:
            last_date = r[idx["Date"]]
        d = last_date.date().isoformat() if isinstance(last_date, datetime.datetime) else str(last_date)
        rec = dict(sheet=ws.title, row=i, station=last_station, station_key=norm_station(last_station), date=d, depth=num(r[idx["Sample"]]))
        for k, key in [("Temp", "temp"), ("Sal", "sal"), ("O2", "o2"), ("NO3", "no3"), ("PO4", "po4"), ("Si", "si"), ("Ba", "ba"), ("Fe", "fe"), ("Mn", "mn"), ("Total PAH", "pah")]:
            rec[key] = num(r[idx[k]]) if idx[k] is not None else None
        tm.append(rec)
print("companion samples:", len(tm), Counter(x["sheet"] for x in tm))
with open(os.path.join(out, "dwh_companion.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(tm[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(tm)

# ---- join on cruise + station + depth (station names differ between the two workbooks: St7 / WS07B, GIP4 / GIP-4)
def join_key(sheet, station):
    k = norm_station(station) or ""
    if sheet == "Late May 2010":
        m = re.match(r"(?:ST|WS)0*(\d+)[A-Z]?(\d*)$", k)      # St2-1 -> 2/1 ; WS02-1 -> 2/1 ; WS07B -> 7 ; WS34B -> 34
        if m:
            k = m.group(1) + ("/" + m.group(2) if m.group(2) else "")
    return k


for x in ree:
    c = (x["cruise"] or "").lower()
    x["sheet"] = ("Early May 2010" if "early may" in c else "Late May 2010" if "late may" in c else
                  "October 2010" if "2010" in c else "October 2011" if "2011" in c else None)
by_station = {}
for t in tm:
    by_station.setdefault((t["sheet"], join_key(t["sheet"], t["station"])), []).append(t)
n_hit = 0
for x in ree:
    cands = by_station.get((x["sheet"], join_key(x["sheet"], x["station"])), [])
    best = min(cands, key=lambda t: abs(t["depth"] - x["depth"])) if cands else None
    ok = best is not None and abs(best["depth"] - x["depth"]) <= 1.5
    x["n_match"] = int(ok)
    x["companion_station"] = best["station"] if ok else None
    for k in ("temp", "sal", "o2", "no3", "po4", "si", "ba", "fe", "mn", "pah"):
        x[k] = best[k] if ok else None
    n_hit += ok
print(f"joined: {n_hit} of {len(ree)} rare-earth samples matched a companion row;",
      sorted(Counter((x['sheet'], x['n_match']) for x in ree).items()))
print("unmatched:", [(x["sheet"], x["station"], x["depth"]) for x in ree if x["n_match"] == 0])
print("matched with oxygen:", sorted(Counter(x['sheet'] for x in ree if x['o2'] is not None).items()))
with open(os.path.join(out, "dwh_joined.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(ree[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(ree)
