"""NIST SRD46 La-vs-M matched 1:1 pairs (ML=[ML]/[M][L], 25 C, K) for a charge-response scale.

Login node, single thread, seconds. Reads only the NIST binary tables plus the 30 Oct primary/features
tables (pinned). Writes NIST-derived tables to workspaces/nist_charge_response_20261109/ (untracked).
Never opens workspaces/ln_structure_increment_20261030/exclusions.tsv.
"""
import hashlib, json, re, struct, sys
sys.dont_write_bytecode = True
from collections import defaultdict
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ND = REPO / "workspaces/ln_selectivity_calibration_20261027/nist46_dat"
OLD = REPO / "workspaces/ln_selectivity_calibration_20261027"
W30 = REPO / "workspaces/ln_structure_increment_20261030"
D30 = REPO / "diagnostics/ln_structure_increment_20261030"
OUT = REPO / "workspaces/nist_charge_response_20261109"
OUT.mkdir(exist_ok=True)
PLAIN = re.compile(r"[+-]?\d*\.?\d+")
EXCL = "workspaces/ln_structure_increment_20261030/exclusions.tsv"

freeze = json.loads((D30 / "FREEZE.json").read_text())
for p, sha in freeze["pins"].items():
    if p == EXCL:
        continue
    assert hashlib.sha256((REPO / p).read_bytes()).hexdigest() == sha, p
print(f"FREEZE pins verified {len(freeze['pins']) - 1}/{len(freeze['pins'])} (exclusions.tsv not opened)")


def lines(fn):
    raw = (ND / fn).read_bytes().decode("latin1")
    return raw.split("\r\n") if "\r\n" in raw else raw.split("\r")


metals = {}
for L in lines("F8548_CODE2.DAT")[1:]:
    m = re.match(r"\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)", L)
    if m:
        metals[int(m[1])] = m[4]
eqs = {}
for L in lines("F8549_CODE3.DAT")[1:]:
    m = re.match(r"\s*(\d+)\s+(\S+)\s+(\S+)", L)
    if m:
        eqs[int(m[1])] = (m[2], m[3])
conds = {int(i): v.split() for i, v in re.findall(r"(\d+) ([KHS\*] [\dtv\*]+ [\d\.\*]+)",
                                                  (ND / "F8550_CODE4.DAT").read_bytes().decode("latin1"))}
notes = {}
for L in lines("F8552_CODE6.DAT"):
    m = re.match(r"\s*(\d+)\s+(\S+)\s+(.*)", L)
    if m:
        notes[int(m[1])] = m[2] + ": " + m[3].strip()
v = (ND / "F8558_VALUES.DAT").read_bytes()
recs = [(int.from_bytes(v[p:p + 2], "little"), v[p + 2], int.from_bytes(v[p + 3:p + 5], "little"), v[p + 5],
         v[p + 6:p + 13].decode("latin1").strip(), v[p + 13], v[p + 14], p) for p in range(15, len(v) - 14, 15)]
b = (ND / "F8547_CODE1.DAT").read_bytes()
lig = {}
for k in range(1, len(b) // 256):
    row = b[k * 256:(k + 1) * 256]
    ptr = int(struct.unpack("<f", row[:4])[0])
    txt = row[4:].decode("latin1").replace("\x00", " ").strip()
    m = re.match(r"(\d{3,4})\t?(\S+) (\S+) (.*?)(?:\t(.*))?$", txt)
    if not m or not (0 <= ptr < len(recs)):
        continue
    lig[str(recs[ptr][0])] = dict(code=m[1], formula=m[2], prot=m[3], name=m[4].strip())

# ---- Positive control A: historical build_pairs rule (extract_ln eq-species ML, strip()<>~, same I) ----
hist = defaultdict(lambda: defaultdict(dict))
for L, M, E, C, val, _, _, _ in recs:
    met = metals.get(M, "?")
    if met not in ("La/3+", "Dy/3+") or eqs.get(E, ("",))[0] != "ML":
        continue
    c = conds.get(C, [])
    if len(c) != 3 or c[0] != "K" or c[1] != "25":
        continue
    try:
        x = float(val.strip("()<>~"))
    except ValueError:
        continue
    hist[L][float(c[2])][met[:2]] = x
n_hist = sum(any({"La", "Dy"} <= set(m) for m in q.values()) for q in hist.values())
print(f"CONTROL A: historical-rule La/Dy pairs = {n_hist} (pairs.tsv: 171)")

# ---- strict curate-rule raw table: ML=[ML]/[M][L], K, 25 C ----
raw = defaultdict(list)  # (lig, metal, cond) -> rows
for L, M, E, C, val, prec, ann, off in recs:
    if eqs.get(E) != ("ML", "[ML]/[M][L]"):
        continue
    c = conds.get(C, [])
    if len(c) != 3 or c[0] != "K" or c[1] != "25":
        continue
    raw[(str(L), metals.get(M, "?"), C)].append(dict(value_raw=val, ann=ann, I=float(c[2]), off=off))
by_lm = defaultdict(dict)  # (lig, metal) -> {cond: rows}
for (L, M, C), rr in raw.items():
    by_lm[(L, M)][C] = rr


def ok1(rr):
    return len(rr) == 1 and PLAIN.fullmatch(rr[0]["value_raw"]) is not None


def pair(L, X, ref="La/3+"):
    """Return (record, reason). Same-condition first (curate rules); nearest-I only when no shared condition."""
    a, x = by_lm.get((L, ref), {}), by_lm.get((L, X), {})
    if not a or not x:
        return None, None
    shared = sorted(set(a) & set(x))
    if shared:
        adm = [c for c in shared if ok1(a[c]) and ok1(x[c]) and a[c][0]["ann"] == x[c][0]["ann"]]
        if not adm:
            why = ("multiple values within condition" if any(len(a[c]) != 1 or len(x[c]) != 1 for c in shared) else
                   "flagged/censored constant" if any(not (ok1(a[c]) and ok1(x[c])) for c in shared) else
                   "annotation mismatch")
            return None, why
        c = min(adm, key=lambda c: (abs(a[c][0]["I"] - .1), a[c][0]["I"], c))
        ra, rx = a[c][0], x[c][0]
        return dict(match="same_I", cond_La=c, cond_X=c, I_La=ra["I"], I_X=rx["I"], dI=0.0, logK_La=float(ra["value_raw"]),
                    logK_X=float(rx["value_raw"]), ann_La=ra["ann"], ann_X=rx["ann"], ann_equal=True,
                    n_shared_adm=len(adm), off_La=ra["off"], off_X=rx["off"]), None
    ca = [c for c in a if ok1(a[c])]; cx = [c for c in x if ok1(x[c])]
    if not ca or not cx:
        return None, "no single plain value (cross-I)"
    c1, c2 = min(((p, q) for p in ca for q in cx),
                 key=lambda t: (abs(a[t[0]][0]["I"] - x[t[1]][0]["I"]),
                                max(abs(a[t[0]][0]["I"] - .1), abs(x[t[1]][0]["I"] - .1)), t))
    ra, rx = a[c1][0], x[c2][0]
    return dict(match="nearest_I", cond_La=c1, cond_X=c2, I_La=ra["I"], I_X=rx["I"], dI=rx["I"] - ra["I"],
                logK_La=float(ra["value_raw"]), logK_X=float(rx["value_raw"]), ann_La=ra["ann"], ann_X=rx["ann"],
                ann_equal=ra["ann"] == rx["ann"], n_shared_adm=0, off_La=ra["off"], off_X=rx["off"]), None


# ---- Positive control B: strict rule on the 115 reproduces primary.tsv exactly ----
prim = pd.read_csv(W30 / "primary.tsv", sep="\t", dtype={"lig_id": str})
feat = pd.read_csv(W30 / "features.tsv", sep="\t", dtype={"lig_id": str})
bad = 0
for _, r in prim.iterrows():
    rec, _ = pair(r.lig_id, "Dy/3+")
    if (rec is None or rec["match"] != "same_I" or f"nist46_{r.lig_id}_C{rec['cond_La']}" != r.record_id
            or rec["logK_La"] != r.logK_La or rec["logK_X"] != r.logK_Dy):
        bad += 1
print(f"CONTROL B: strict-rule La/Dy on the 115 primary ligands: {len(prim) - bad}/{len(prim)} identical records")
assert bad == 0
allLaDy = [L for (L, M) in by_lm if M == "Dy/3+" and pair(L, "Dy/3+")[0] is not None]
print(f"  strict-rule La/Dy pairs over the whole database: {len(allLaDy)} ligands "
      f"(same-I {sum(pair(L, 'Dy/3+')[0]['match'] == 'same_I' for L in allLaDy)})")

# ---- all partner metals ----
la_ligs = sorted({L for (L, M) in by_lm if M == "La/3+"})
print(f"ligands with any strict-rule La 1:1 K at 25 C: {len(la_ligs)}")
F115 = set(feat.lig_id)
famof = dict(zip(prim.lig_id, prim.family))
rows, why = [], defaultdict(lambda: defaultdict(int))
for X in sorted({M for (_, M) in by_lm if M != "La/3+"}):
    for L in la_ligs:
        rec, reason = pair(L, X)
        if reason:
            why[X][reason] += 1
        if rec:
            info = lig.get(L, {})
            rows.append(dict(lig_id=L, name=info.get("name", "?"), prot=info.get("prot", "?"), metal=X, **rec,
                             y_log=rec["logK_La"] - rec["logK_X"], has_features=L in F115, family=famof.get(L, "")))
P = pd.DataFrame(rows)
P.to_csv(OUT / "la_pairs_all_metals.tsv", sep="\t", index=False)
cnt = (P.groupby("metal").agg(n=("lig_id", "size"), same_I=("match", lambda s: (s == "same_I").sum()),
                              feat=("has_features", "sum"),
                              feat_same_I=("has_features", lambda s: (s & (P.loc[s.index, "match"] == "same_I")).sum()),
                              fam=("family", lambda s: s[s != ""].nunique()))
       .sort_values("n", ascending=False))
cnt["excluded_shared"] = [sum(why[m].values()) for m in cnt.index]
cnt.to_csv(OUT / "pair_counts.tsv", sep="\t")
print("\nPartner metals with >=15 La-paired ligands:")
print(cnt[cnt.n >= 15].to_string())
for m in ("Ca/2+", "Sr/2+", "Al/3+", "Fe/3+"):
    print(m, "shared-condition exclusions:", dict(why[m]))

# ---- La/Ca noise-floor table: every admissible same-condition pair, all I ----
nz = []
for L in la_ligs:
    a, x = by_lm.get((L, "La/3+"), {}), by_lm.get((L, "Ca/2+"), {})
    for c in sorted(set(a) & set(x)):
        if ok1(a[c]) and ok1(x[c]) and a[c][0]["ann"] == x[c][0]["ann"]:
            nz.append(dict(lig_id=L, cond=c, I=a[c][0]["I"], y1=float(a[c][0]["value_raw"]) - float(x[c][0]["value_raw"]),
                           has_features=L in F115))
pd.DataFrame(nz).to_csv(OUT / "la_ca_all_conditions.tsv", sep="\t", index=False)
smi = pd.read_csv(OLD / "smiles.tsv", sep="\t", dtype=str)
ok_smi = set(smi.lig_id[smi.formula_ok == "yes"])
ca = P[P.metal == "Ca/2+"]
print(f"\nLa/Ca: {len(ca)} ligands; with 30 Oct feature vectors {ca.has_features.sum()}; formula-checked SMILES outside "
      f"the 115: {len(set(ca.lig_id) & (ok_smi - F115))}")
