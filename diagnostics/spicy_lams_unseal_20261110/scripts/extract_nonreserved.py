"""Unseal SpyCI-LAMBS outcomes for the 600 non-reserved orthologs (Jacob, 10 Nov 2026).
The 16 reserved-panel rows are skipped by workbook row index before any of their outcome cells is stored or printed.
Outputs: nonreserved_profiles.tsv (metadata, raw desorption, total metal, normalized logD mean/SD),
reserved_metadata_only.tsv (identifiers, sequence, taxonomy, total-metal QC; no ion-resolved outcome),
EXTRACTION.json (counts and hashes). Prints counts only."""
import csv, hashlib, json, sys
from pathlib import Path
import openpyxl
ROOT = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs")
WB = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/spicy_lams/41589_2026_2176_MOESM3_ESM.xlsx")
OUT = ROOT / "workspaces/spicy_lams_unseal_20261110"; OUT.mkdir(parents=True, exist_ok=True)
sha = hashlib.sha256(WB.read_bytes()).hexdigest()
assert sha == "391933234dcf5d3862de958da88afd46fe0f6a8bcb8e3fb1a22640bf976a99f1", sha
sel = json.loads((ROOT / "spicy_lams_matched_homolog/selection_manifest.json").read_text())
reserved_acc = {m["accession"] for p in sel["pairs"] for m in p["members"]}
assert len(reserved_acc) == 16
orth = list(csv.DictReader((ROOT / "workspaces/spicy_lams_inventory_20260923/inventory_v2/ORTHOLOGS.tsv").open(), delimiter="\t"))
row_acc = {int(r["workbook_row"]): r["accession"] for r in orth}
reserved_rows = {i for i, a in row_acc.items() if a in reserved_acc}
assert len(reserved_rows) == 16 and len(row_acc) == 616
wb = openpyxl.load_workbook(WB, read_only=True, data_only=True); ws = wb["Data"]
hdr = None; kept = []; held = []
for ridx, row in enumerate(ws.iter_rows(min_col=1, max_col=101, values_only=True), start=1):   # CX:CZ never requested
    if ridx == 1:
        hdr = list(row); continue
    if ridx not in row_acc:
        continue
    if ridx in reserved_rows:
        held.append(list(row[:19]) + list(row[64:71]))     # A:S and BM:BS only; outcome cells discarded unread
        continue
    kept.append(list(row))
assert len(kept) == 600 and len(held) == 16, (len(kept), len(held))
# unique column names: replicate columns already unique; normalized logD mean/SD share names -> suffix them
names = hdr[:71] + [f"nlogD_mean_{h}" for h in hdr[71:86]] + [f"nlogD_sd_{h}" for h in hdr[86:101]]
names = ["workbook_row", "accession"] + names
with (OUT / "nonreserved_profiles.tsv").open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(names)
    i = 0
    for ridx in sorted(row_acc):
        if ridx in reserved_rows: continue
        w.writerow([ridx, row_acc[ridx]] + ["" if v is None else v for v in kept[i]]); i += 1
with (OUT / "reserved_metadata_only.tsv").open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["workbook_row", "accession"] + hdr[:19] + hdr[64:71])
    for ridx, r in zip(sorted(reserved_rows), held):
        w.writerow([ridx, row_acc[ridx]] + ["" if v is None else v for v in r])
info = dict(workbook=str(WB), workbook_sha256=sha, rows_total=616, rows_unsealed=600, rows_still_sealed=16,
            reserved_workbook_rows=sorted(reserved_rows), columns_unsealed="T:BL, BT:CW (plus A:S, BM:BS metadata)",
            columns_never_read="CX:CZ",
            outputs={n: hashlib.sha256((OUT / n).read_bytes()).hexdigest() for n in ("nonreserved_profiles.tsv", "reserved_metadata_only.tsv")})
json.dump(info, (OUT / "EXTRACTION.json").open("w"), indent=1)
print(json.dumps({k: v for k, v in info.items() if k != "outputs"}, indent=1)); print("columns written:", len(names))
