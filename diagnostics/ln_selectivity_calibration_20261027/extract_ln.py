"""Extract Ln(III) ML stability data (log K, dH, dS) for every ligand in NIST SRD 46 v8.

Ligand records (F8547_CODE1.DAT) are 256-byte rows: a float32 pointer to the ligand's first row in
F8558_VALUES.DAT, then "<code>\t<formula> <HnL> <name>\t<class tags>". VALUES rows are 15 bytes
(ligand id u16, metal u8, equilibrium u16, condition u8, value char[7], ...), decoded with the
code tables as in diagnostics/lady_relaxed_exchange_20261024/tools/nist46.py.

Output: ln_ML.tsv, one row per (ligand, metal, quantity, T, I) for the 15 Ln(III) + Y(III),
equilibrium ML = [ML]/[M][L] only.
"""
import re, struct, sys
from pathlib import Path

D = Path(sys.argv[1] if len(sys.argv) > 1 else "nist46_dat")
LN = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y"]


def lines(fn):
    raw = (D / fn).read_bytes().decode("latin1")
    return raw.split("\r\n") if "\r\n" in raw else raw.split("\r")


code2 = {}
for L in lines("F8548_CODE2.DAT")[1:]:
    m = re.match(r"\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)", L)
    if m: code2[int(m.group(1))] = m.group(4)
code3 = {}
for L in lines("F8549_CODE3.DAT")[1:]:
    m = re.match(r"\s*(\d+)\s+(\S+)\s+(\S+)", L)
    if m: code3[int(m.group(1))] = (m.group(2), m.group(3))
raw4 = (D / "F8550_CODE4.DAT").read_bytes().decode("latin1")
code4 = {int(i): v.strip() for i, v in re.findall(r"(\d+) ([KHS\*] [\dtv\*]+ [\d\.\*]+)", raw4)}

v = (D / "F8558_VALUES.DAT").read_bytes()
recs = []
for p in range(15, len(v) - 14, 15):
    r = v[p:p + 15]
    recs.append((int.from_bytes(r[0:2], "little"), r[2], int.from_bytes(r[3:5], "little"), r[5], r[6:13].decode("latin1").strip()))

b = (D / "F8547_CODE1.DAT").read_bytes()
lig = {}
for k in range(1, len(b) // 256):
    row = b[k * 256:(k + 1) * 256]
    ptr = int(struct.unpack("<f", row[:4])[0])
    txt = row[4:].decode("latin1").replace("\x00", " ").strip()
    m = re.match(r"(\d{3,4})\t?(\S+) (\S+) (.*?)(?:\t(.*))?$", txt)
    if not m or not (0 <= ptr < len(recs)): continue
    lig[recs[ptr][0]] = (m.group(1), m.group(2), m.group(3), m.group(4).strip(), (m.group(5) or "").strip())

out = ["lig_id\tcode\tformula\tprot\tname\ttags\tmetal\tqty\tT\tI\tvalue"]
for L, M, E, C, val in recs:
    met = code2.get(M, "?").split("/")[0]
    if met not in LN or code2.get(M, "").split("/")[-1] != "3+": continue
    if code3.get(E, ("",))[0] != "ML": continue
    c = code4.get(C, "").split()
    if len(c) != 3: continue
    try: x = float(val.strip("()<>~"))
    except ValueError: continue
    info = lig.get(L, ("?", "?", "?", "?", "?"))
    out.append("\t".join(map(str, (L, *info, met, c[0], c[1], c[2], x))))
Path("ln_ML.tsv").write_text("\n".join(out) + "\n")
print(len(lig), "ligands;", len(out) - 1, "Ln ML rows", file=sys.stderr)
