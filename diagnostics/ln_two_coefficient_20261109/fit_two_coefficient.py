"""Two-coefficient radius-quadratic selectivity model vs M1.  Frozen design: PLAN.md (same dir).

Login node, single thread, seconds.  Writes NIST-derived tables to workspaces/ln_two_coefficient_20261109/ (untracked).
"""
import hashlib, json, re, sys
sys.dont_write_bytecode = True
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[2]
D30 = REPO / "diagnostics/ln_structure_increment_20261030"
W30 = REPO / "workspaces/ln_structure_increment_20261030"
OLD = REPO / "workspaces/ln_selectivity_calibration_20261027"
OUT = REPO / "workspaces/ln_two_coefficient_20261109"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(D30))
from evaluate import nested, ALPHAS, ridge_fit  # noqa: E402  (code only; main() not run)
from structure import M1  # noqa: E402

KC = 0.00198720425864083 * 298.15 * np.log(10)
LN = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
R9 = dict(zip(LN, [1.216, 1.196, 1.179, 1.163, 1.132, 1.120, 1.107, 1.095, 1.083, 1.072, 1.062, 1.052, 1.042, 1.032]))
NL = LN[1:]  # 13 non-La metals
DR = np.array([R9[m] - R9["La"] for m in NL]); DRDY = R9["Dy"] - R9["La"]
PHI1, PHI2 = DR / DRDY, (DR / DRDY) ** 2
IDY = NL.index("Dy")
assert abs(PHI1[IDY] - 1) < 1e-12 and abs(PHI2[IDY] - 1) < 1e-12
LAMS = [2 * a for a in ALPHAS]
NB, SEED = 2000, 42
EXCL = "workspaces/ln_structure_increment_20261030/exclusions.tsv"  # reserved-out-of-caution: never opened/hashed

# ---------------- pins (exclusions.tsv deliberately skipped) ----------------
freeze = json.loads((D30 / "FREEZE.json").read_text())
for p, sha in freeze["pins"].items():
    if p == EXCL:
        continue
    assert hashlib.sha256((REPO / p).read_bytes()).hexdigest() == sha, p
print(f"FREEZE pins verified: {len(freeze['pins']) - 1} of {len(freeze['pins'])} (exclusions.tsv not opened)")

# ---------------- Step 0: M1 positive control ----------------
d = pd.read_csv(W30 / "primary.tsv", sep="\t", dtype={"lig_id": str})
f = pd.read_csv(W30 / "features.tsv", sep="\t", dtype={"lig_id": str})
d = d.merge(f.drop(columns=["family", "ligand_charge"]), on="lig_id", validate="one_to_one")
assert len(d) == freeze["n_primary"] == 115 and d.lig_id.is_unique
assert np.allclose(d.y_kcal, -freeze["kcal_per_log_unit"] * (d.logK_Dy - d.logK_La), atol=1e-12)
d = d.reset_index(drop=True)
pm1, meta_m1 = nested(d, M1)
rm1 = float(np.sqrt(np.mean((pm1 - d.y_kcal.values) ** 2)))
print(f"STEP0 M1 family-held-out RMSE {rm1:.4f} (target 0.854; |diff| {abs(rm1-0.854):.4f}); alphas "
      f"{[m['alpha'] for m in meta_m1]}")
if abs(rm1 - 0.854) > 0.005:
    sys.exit("POSITIVE CONTROL FAILED")
fam = d.family.values; FAMS = sorted(set(fam)); LIDS = d.lig_id.values; NLIG = len(d)
X = d[M1].values.astype(float); YDY = d.y_kcal.values

# ---------------- multi-metal extraction at each primary condition ----------------
keep = set(LIDS)
nd = OLD / "nist46_dat"
lines = lambda fn: (nd / fn).read_bytes().decode("latin1").splitlines()
metals = {}
for line in lines("F8548_CODE2.DAT")[1:]:
    m = re.match(r"\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)", line)
    if m:
        metals[int(m[1])] = m[4]
eqs = {}
for line in lines("F8549_CODE3.DAT")[1:]:
    m = re.match(r"\s*(\d+)\s+(\S+)\s+(\S+)", line)
    if m:
        eqs[int(m[1])] = (m[2], m[3])
conds = {int(i): v.split() for i, v in re.findall(r"(\d+) ([KHS\*] [\dtv\*]+ [\d\.\*]+)",
                                                  (nd / "F8550_CODE4.DAT").read_bytes().decode("latin1"))}
WANT = {m + "/3+" for m in LN}
v = (nd / "F8558_VALUES.DAT").read_bytes()
raw = defaultdict(list)
pm_seen = 0
for off in range(15, len(v) - 14, 15):
    r = v[off:off + 15]
    lid = str(int.from_bytes(r[:2], "little"))
    if lid not in keep:  # filter to the 115 admitted ids before parsing any value
        continue
    mt = metals.get(r[2], "")
    if mt == "Pm/3+":
        pm_seen += 1
    if mt not in WANT or eqs.get(int.from_bytes(r[3:5], "little")) != ("ML", "[ML]/[M][L]"):
        continue
    c = conds.get(r[5], [])
    if len(c) != 3 or c[:2] != ["K", "25"]:
        continue
    raw[(lid, r[5], mt[:2])].append(dict(value_raw=r[6:13].decode("latin1").strip(), annotation_code=r[14], off=off))
print(f"Pm rows for the 115 ids: {pm_seen} (excluded, no frozen radius)")
long, drops = [], defaultdict(int)
series = {}
for _, row in d.iterrows():
    cond = int(row.record_id.rsplit("_C", 1)[1])
    got = {}
    for m in LN:
        vv = raw.get((row.lig_id, cond, m), [])
        if not vv:
            continue
        if len(vv) != 1:
            drops["multiple_values"] += 1; continue
        x = vv[0]
        if not re.fullmatch(r"[+-]?\d*\.?\d+", x["value_raw"]):
            drops["censored_or_flagged"] += 1; continue
        if x["annotation_code"] != int(row.annotation_code):
            drops["annotation_mismatch"] += 1; continue
        got[m] = float(x["value_raw"])
    assert got.get("La") == row.logK_La and got.get("Dy") == row.logK_Dy, row.lig_id
    series[row.lig_id] = got
    for m in NL:
        if m in got:
            long.append(dict(lig_id=row.lig_id, family=row.family, metal=m, logK_M=got[m], logK_La=got["La"],
                             y_kcal=-KC * (got[m] - got["La"])))
L = pd.DataFrame(long)
L.to_csv(OUT / "multi_metal_points.tsv", sep="\t", index=False)
nmet = L.groupby("lig_id").size()
complete = [l for l in LIDS if nmet.get(l, 0) == 13]
print(f"admitted non-La points {len(L)} over {L.lig_id.nunique()} ligands; drops {dict(drops)}; "
      f"metals/ligand (incl. La) median {int(nmet.median()) + 1}; Dy-only ligands {(nmet == 1).sum()}; "
      f"complete 13-metal series {len(complete)}")
print("  ligands per metal:", L.groupby("metal").size().reindex(NL).to_dict())
lig_ix = {l: i for i, l in enumerate(LIDS)}
PL = np.array([lig_ix[l] for l in L.lig_id]); PM = np.array([NL.index(m) for m in L.metal]); PY = L.y_kcal.values

# ---------------- R2-plus count ----------------
fv = set(pd.read_csv(OLD / "features.tsv", sep="\t", usecols=["lig_id"], dtype=str).lig_id) | set(f.lig_id)
pr = set(pd.read_csv(OLD / "pairs.tsv", sep="\t", usecols=["lig_id"], dtype=str).lig_id)
plus = sorted(fv - pr)
print(f"R2-plus: ligands with a feature vector {len(fv)}; lacking La or Dy (not in pairs.tsv) {len(plus)} -> "
      f"{'skipped' if not plus else plus}")


# ---------------- R2 model ----------------
def fit_r2(tr_lig, use_u, dy_only, lam):
    sc = StandardScaler().fit(X[tr_lig])
    Z = sc.transform(X)
    if dy_only:
        sel = np.isin(PL, tr_lig) & (PM == IDY)
    else:
        sel = np.isin(PL, tr_lig)
    pl, pm, y = PL[sel], PM[sel], PY[sel]
    p = X.shape[1]
    head = np.eye(13)[pm] if use_u else np.c_[PHI1[pm], PHI2[pm]]
    A = np.c_[head, Z[pl] * PHI1[pm][:, None], Z[pl] * PHI2[pm][:, None]]
    k0 = head.shape[1]
    pen = np.zeros((2 * p, A.shape[1])); pen[:, k0:] = np.sqrt(lam) * np.eye(2 * p)
    th = np.linalg.lstsq(np.r_[A, pen], np.r_[y, np.zeros(2 * p)], rcond=None)[0]
    al, be = th[k0:k0 + p], th[k0 + p:]

    def predict(lig):  # -> (n, 13) predictions, plus a, b (Dy-scaled)
        z = Z[lig]
        if use_u:
            a, b, u = z @ al, z @ be, th[:13]
        else:
            a, b, u = th[0] + z @ al, th[1] + z @ be, np.zeros(13)
        return a[:, None] * PHI1 + b[:, None] * PHI2 + u, a, b
    return predict, th


def nested_r2(use_u, dy_only):
    pred = np.full((NLIG, 13), np.nan); ab = np.full((NLIG, 2), np.nan); meta = []
    for F in FAMS:
        te = np.flatnonzero(fam == F); trm = fam != F
        scores = []
        for lam in LAMS:
            inner = np.full(NLIG, np.nan)
            for H in sorted(set(fam[trm])):
                va = np.flatnonzero(trm & (fam == H)); tr = np.flatnonzero(trm & (fam != H))
                pr_, _ = fit_r2(tr, use_u, dy_only, lam)
                inner[va] = pr_(va)[0][:, IDY]
            scores.append(float(np.mean((inner[trm] - YDY[trm]) ** 2)))
        lam = LAMS[int(np.argmin(scores))]
        pr_, th = fit_r2(np.flatnonzero(trm), use_u, dy_only, lam)
        P, a, b = pr_(te); pred[te] = P; ab[te, 0] = a; ab[te, 1] = b
        meta.append(dict(fold=F, lam=lam, inner=scores, u=(th[:13].tolist() if use_u else None)))
    return pred, ab, meta


VAR = {"R2-dy": (True, True), "R2-all": (True, False), "R2-all-nou": (False, False)}
P = {}; AB = {}; META = {}
for k, (uu, dy) in VAR.items():
    P[k], AB[k], META[k] = nested_r2(uu, dy)
    print(f"{k}: selected lambda per fold {[m['lam'] for m in META[k]]}")

# ---------------- scores + paired bootstrap ----------------
rng = np.random.default_rng(SEED); BL = rng.integers(0, NLIG, (NB, NLIG))
rngf = np.random.default_rng(SEED)
fidx = {F: np.flatnonzero(fam == F) for F in FAMS}
BF = [np.concatenate([fidx[g] for g in rngf.choice(FAMS, len(FAMS), replace=True)]) for _ in range(NB)]


def rmse_ci_pair(e_ref, e_var, idx_l, idx_f):
    """e_*: per-unit squared-error sums and counts. Returns rmse_ref, rmse_var, diff, ligand CI, family CI."""
    (s0, n0), (s1, n1) = e_ref, e_var
    r0 = np.sqrt(s0.sum() / n0.sum()); r1 = np.sqrt(s1.sum() / n1.sum())
    out = [r0, r1, r0 - r1]
    for B in (idx_l, idx_f):
        dd = np.array([np.sqrt(s0[b].sum() / n0[b].sum()) - np.sqrt(s1[b].sum() / n1[b].sum()) for b in B])
        out.append((np.percentile(dd, 2.5), np.percentile(dd, 97.5)))
    return out


res = {"step0_M1_rmse": rm1, "n_points": int(len(L)), "drops": dict(drops), "n_complete": len(complete),
       "R2_plus_n": len(plus), "primary": {}, "secondary": {}}
e1 = ((pm1 - YDY) ** 2, np.ones(NLIG))
print("\n== PRIMARY: family-held-out RMSE of Dy-La, n=115 (kcal/mol) ==")
print(f"  M1 {rm1:.3f}")
for k in VAR:
    ev = ((P[k][:, IDY] - YDY) ** 2, np.ones(NLIG))
    r0, r1, df, cl, cf = rmse_ci_pair(e1, ev, BL, BF)
    res["primary"][k] = dict(rmse=r1, m1_minus=df, ci_lig=cl, ci_fam=cf, mae=float(np.abs(P[k][:, IDY] - YDY).mean()))
    print(f"  {k:11s} {r1:.3f}  M1-variant {df:+.3f}  lig CI [{cl[0]:+.3f}, {cl[1]:+.3f}]  fam CI [{cf[0]:+.3f}, {cf[1]:+.3f}]")
print(f"  R2-dy check |RMSE - 0.854| = {abs(res['primary']['R2-dy']['rmse'] - 0.854):.4f}; "
      f"max |pred R2-dy - pred M1| = {np.abs(P['R2-dy'][:, IDY] - pm1).max():.2e}")

r2a = res["primary"]["R2-all"]
reading = ("Strong" if r2a["rmse"] <= 0.60 else
           "Supported" if (r2a["m1_minus"] >= 0.05 and r2a["ci_lig"][0] > 0) else "Not supported")
if r2a["rmse"] <= 0.60 and r2a["m1_minus"] >= 0.05 and r2a["ci_lig"][0] > 0:
    reading = "Strong (and Supported)"
res["reading"] = reading
print(f"  FROZEN READING (R2-all): {reading}")

# secondary: complete-series held-out ligands
ci = np.array([lig_ix[l] for l in complete])
Yser = np.full((NLIG, 13), np.nan)
for l in LIDS:
    g = series[l]
    for j, m in enumerate(NL):
        if m in g:
            Yser[lig_ix[l], j] = -KC * (g[m] - g["La"])
svd_pred = np.full((NLIG, 13), np.nan)
for F in FAMS:
    te = np.flatnonzero(fam == F)
    trc = [i for i in ci if fam[i] != F]
    V1 = np.linalg.svd(Yser[trc], full_matrices=False)[2][0]
    svd_pred[te] = pm1[te, None] * V1[None, :] / V1[IDY]
fam_c = fam[ci]
rngs = np.random.default_rng(SEED); BLs = rngs.integers(0, len(ci), (NB, len(ci)))
rngsf = np.random.default_rng(SEED); FC = sorted(set(fam_c))
fci = {F: np.flatnonzero(fam_c == F) for F in FC}
BFs = [np.concatenate([fci[g] for g in rngsf.choice(FC, len(FC), replace=True)]) for _ in range(NB)]
es = (((svd_pred[ci] - Yser[ci]) ** 2).sum(1), np.full(len(ci), 13.0))
print(f"\n== SECONDARY: complete-series held-out ligands n={len(ci)} in {len(FC)} families, RMSE over 13 metals ==")
print(f"  M1->SVD shape 1 {np.sqrt(es[0].sum() / es[1].sum()):.3f}")
res["secondary"]["M1_svd1"] = float(np.sqrt(es[0].sum() / es[1].sum()))
for k in VAR:
    ev = (((P[k][ci] - Yser[ci]) ** 2).sum(1), np.full(len(ci), 13.0))
    r0, r1, df, cl, cf = rmse_ci_pair(es, ev, BLs, BFs)
    res["secondary"][k] = dict(rmse=r1, m1svd_minus=df, ci_lig=cl, ci_fam=cf)
    flag = "  (Dy-only training: other metals not identified; extrapolation is arbitrary)" if k == "R2-dy" else ""
    print(f"  {k:11s} {r1:.3f}  M1svd-variant {df:+.3f}  lig CI [{cl[0]:+.3f}, {cl[1]:+.3f}]  "
          f"fam CI [{cf[0]:+.3f}, {cf[1]:+.3f}]{flag}")
for k in ("M1svd", "R2-all"):
    Pk = svd_pred if k == "M1svd" else P[k]
    print(f"  {k} Dy-only RMSE on these {len(ci)}: {np.sqrt(np.mean((Pk[ci, IDY] - Yser[ci, IDY]) ** 2)):.3f}")

# ---------------- descriptive ----------------
print("\n== Per-family RMSE (Dy-La) ==")
pf = []
for F in FAMS:
    i = fidx[F]
    pf.append(dict(family=F, n=len(i), M1=np.sqrt(np.mean((pm1[i] - YDY[i]) ** 2)),
                   R2_all=np.sqrt(np.mean((P["R2-all"][i, IDY] - YDY[i]) ** 2)),
                   R2_all_nou=np.sqrt(np.mean((P["R2-all-nou"][i, IDY] - YDY[i]) ** 2)),
                   n_points=int(np.isin(PL, i).sum())))
PF = pd.DataFrame(pf); print(PF.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
PF.to_csv(OUT / "per_family.tsv", sep="\t", index=False)

print("\n== Named ligands: residual obs - pred (Dy-La), kcal/mol ==")
NAMED = {"TTHA": None, "CDTA": "600", "TMDTA": None, "DTPA-BMA": "666", "oxolane-2,3,4,5-tetracarboxylate": "3687",
         "DTPA": "655", "EDTA": "576"}
named = []
for nm, lid in NAMED.items():
    if lid is None or lid not in lig_ix:
        print(f"  {nm}: not in the 115 (excluded upstream)"); continue
    i = lig_ix[lid]
    named.append(dict(name=nm, lig_id=lid, obs=YDY[i], res_M1=YDY[i] - pm1[i], res_R2all=YDY[i] - P["R2-all"][i, IDY],
                      res_R2nou=YDY[i] - P["R2-all-nou"][i, IDY], n_metals=int(nmet.get(lid, 0)) + 1))
    print(f"  {nm:34s} obs {YDY[i]:+.2f}  M1 {YDY[i]-pm1[i]:+.2f}  R2-all {YDY[i]-P['R2-all'][i, IDY]:+.2f}  "
          f"R2-all-nou {YDY[i]-P['R2-all-nou'][i, IDY]:+.2f}  metals {int(nmet.get(lid, 0)) + 1}")

print("\n== Held-out prediction of per-ligand fitted a, b (Å units) from features, ligands with >=4 metals ==")
ab_rows = []
for l in LIDS:
    g = series[l]
    if len(g) < 4:
        continue
    ms = [m for m in NL if m in g]
    dr = np.array([R9[m] - R9["La"] for m in ms]); y = np.array([-KC * (g[m] - g["La"]) for m in ms])
    c = np.linalg.lstsq(np.c_[dr, dr ** 2], y, rcond=None)[0]
    ab_rows.append(dict(lig_id=l, family=fam[lig_ix[l]], n_metals=len(g), a=c[0], b=c[1]))
AT = pd.DataFrame(ab_rows)
dsub = d.set_index("lig_id").loc[AT.lig_id].reset_index()
abres = {}
for par in ("a", "b"):
    dd_ = dsub.copy(); dd_["y_kcal"] = AT[par].values
    pr_, mt_ = nested(dd_, M1)
    ss = ((AT[par].values - pr_) ** 2).sum(); st = ((AT[par].values - AT[par].mean()) ** 2).sum()
    from scipy.stats import spearmanr
    abres[par] = dict(n=len(AT), n_fam=int(AT.family.nunique()), R2=1 - ss / st, spearman=spearmanr(AT[par], pr_)[0],
                      alphas=[m["alpha"] for m in mt_])
    AT[f"{par}_pred"] = pr_
    print(f"  {par}: n={len(AT)} ligands / {AT.family.nunique()} families, held-out R2 {1 - ss / st:+.3f}, "
          f"Spearman(obs,pred) {abres[par]['spearman']:+.2f}, alphas {abres[par]['alphas']}")
AT.to_csv(OUT / "per_ligand_ab.tsv", sep="\t", index=False)
res["ab"] = abres

# per-metal offsets (fold-wise) for R2-all
U = np.array([m["u"] for m in META["R2-all"]])
print("\n== R2-all u_M across outer folds (median [min, max]) ==")
print("  " + "  ".join(f"{m} {np.median(U[:, j]):+.2f}[{U[:, j].min():+.2f},{U[:, j].max():+.2f}]" for j, m in enumerate(NL)))

pred_tab = d[["lig_id", "name", "family"]].copy(); pred_tab["obs_DyLa"] = YDY; pred_tab["M1"] = pm1
for k in VAR:
    pred_tab[k] = P[k][:, IDY]; pred_tab[f"{k}_a_dy"] = AB[k][:, 0]; pred_tab[f"{k}_b_dy"] = AB[k][:, 1]
pred_tab["n_metals"] = [int(nmet.get(l, 0)) + 1 for l in LIDS]
pred_tab.to_csv(OUT / "predictions.tsv", sep="\t", index=False)
res["per_family"] = pf; res["named"] = named; res["meta"] = META
json.dump(res, open(OUT / "summary.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("\nwrote", OUT)
