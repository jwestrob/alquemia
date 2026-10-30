"""PubChem SMILES for each ligand in pairs.tsv (by name variants), checked against the NIST formula
(neutral, fully protonated HnL form). Writes smiles.tsv: lig_id, query, smiles, formula_ok."""
import csv, re, time, urllib.parse, urllib.request
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def nist_formula(f):  # 'C10H16N2O8' with explicit 1s -> RDKit style
    return re.sub(r"([A-Z][a-z]?)1(?![0-9])", r"\1", f)

def query(name):
    u = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + urllib.parse.quote(name, safe="") + "/property/IsomericSMILES/TXT"
    try:
        return urllib.request.urlopen(u, timeout=20).read().decode().split("\n")[0].strip() or None
    except Exception:
        return None

def variants(name):
    out = [name]
    for m in re.findall(r"\(([^()]+)\)", name): out.append(m)          # parenthetical common name
    base = re.sub(r"\s*\([^()]*\)\s*$", "", name).strip()
    if base != name: out.append(base)
    out.append(re.sub(r"^(DL|D|L|meso|rac|trans|cis)-", "", base))
    return list(dict.fromkeys(o for o in out if o))

rows = list(csv.DictReader(open("pairs.tsv"), delimiter="\t"))
seen = {}
try:
    for r in csv.DictReader(open("smiles.tsv"), delimiter="\t"): seen[r["lig_id"]] = r
except FileNotFoundError:
    pass
w = open("smiles.tsv", "a")
if not seen: w.write("lig_id\tquery\tsmiles\tformula_ok\n")
for r in rows:
    if r["lig_id"] in seen: continue
    best = ("", "", "no")
    for q in variants(r["name"]):
        s = query(q); time.sleep(0.25)
        if not s: continue
        m = Chem.MolFromSmiles(s)
        ok = m is not None and CalcMolFormula(m) == nist_formula(r["formula"])
        if ok: best = (q, s, "yes"); break
        if not best[1]: best = (q, s, "no")
    w.write(f"{r['lig_id']}\t{best[0]}\t{best[1]}\t{best[2]}\n"); w.flush()
