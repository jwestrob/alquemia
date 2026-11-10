#!/usr/bin/env python
"""Charge line on every NIST La/X pair. Lead-agent addition, 2026-11-10 UTC.

Outside the subagent's frozen PLAN.md (sha256 c2c04a1f...). Design fixed here before running.

Question
    On every 1:1 pair NIST lists for La and a comparison metal X (not only the ligands that have
    donor-feature vectors), how does log K(La) - log K(X) depend on the net charge of the fully
    deprotonated ligand?

Charge
    NIST writes the neutral form as HnL, so z_L = -n. An explicit suffix (L/-, HL/-, L/3-) gives the
    charge of the listed form, from which the n protons are removed.

Fixed choices
    - pairs measured at the same ionic strength only (match == same_I);
    - ordinary least squares of y on q = -z_L (q = number of negative charges on the ligand);
    - 2,000 ligand bootstraps, seed 42, for slope, intercept and Spearman rho;
    - leave-one-out RMSE of the line against the leave-one-out mean;
    - donor class from name keywords (O-only, N-containing, other), checked against features.tsv on
      the ligands that have feature vectors; class-held-out prediction in both directions;
    - lanthanide-specific lines, y = log K(X) - log K(Ca), from ligands whose La/Ca and La/X pairs
      share the La condition.

Inputs  workspaces/nist_charge_response_20261109/la_pairs_all_metals.tsv
        workspaces/ln_selectivity_calibration_20261027/features.tsv
Outputs workspaces/nist_charge_response_20261109/charge_lines.tsv
        workspaces/nist_charge_response_20261109/charge_line_ligands.tsv
        workspaces/nist_charge_response_20261109/charge_line.log
"""
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

REPO = "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs"
W = f"{REPO}/workspaces/nist_charge_response_20261109"
FEAT = f"{REPO}/workspaces/ln_selectivity_calibration_20261027/features.tsv"
SEED, NBOOT = 42, 2000

N_KEYS = ["nitroso", "amin", "imino", "nitrilo", "aza", "azol", "pyrid", "quinolin", "phenanthrolin",
          "histid", "glycin", "alanin", "hydrazin", "oxime", "piperidin", "pyrrol", "purin", "pyrimid"]
OTHER_KEYS = ["cyan", "thio", "thia", "mercapto", "sulfid", "fluorid", "chlorid", "bromid", "iodid",
              "selen", "phosphin", "arsen"]
LN = ["Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y"]
EARTHS = ["Ca", "Sr", "Ba", "Mg"]


def charge(prot):
    m = re.match(r"^(H(\d*))?L(?:/(\d*)([+-]))?$", str(prot))
    if not m:
        return np.nan
    n_h = 0 if m.group(1) is None else (1 if m.group(2) == "" else int(m.group(2)))
    q_form = 0
    if m.group(4):
        k = 1 if m.group(3) == "" else int(m.group(3))
        q_form = k if m.group(4) == "+" else -k
    return q_form - n_h


def donor_class(name):
    s = str(name).lower()
    if any(k in s for k in OTHER_KEYS):
        return "other"
    if any(k in s for k in N_KEYS):
        return "N_containing"
    return "O_only"


def fit(q, y, rng):
    q, y = np.asarray(q, float), np.asarray(y, float)
    n = len(y)
    out = {"n": n}
    if n < 8 or len(set(q)) < 3:
        return out
    b1, b0 = np.polyfit(q, y, 1)
    res = y - (b0 + b1 * q)
    rho = stats.spearmanr(q, y).statistic
    bs, bi, br = [], [], []
    for _ in range(NBOOT):
        i = rng.integers(0, n, n)
        if len(set(q[i])) < 2:
            continue
        s1, s0 = np.polyfit(q[i], y[i], 1)
        bs.append(s1)
        bi.append(s0)
        br.append(stats.spearmanr(q[i], y[i]).statistic)
    # leave-one-out
    e_line, e_mean = [], []
    for j in range(n):
        m = np.ones(n, bool)
        m[j] = False
        if len(set(q[m])) < 2:
            continue
        s1, s0 = np.polyfit(q[m], y[m], 1)
        e_line.append(y[j] - (s0 + s1 * q[j]))
        e_mean.append(y[j] - y[m].mean())
    pc = lambda v: np.percentile(v, [2.5, 97.5])
    out.update({
        "slope": b1, "slope_lo": pc(bs)[0], "slope_hi": pc(bs)[1],
        "intercept": b0, "intercept_lo": pc(bi)[0], "intercept_hi": pc(bi)[1],
        "rho": rho, "rho_lo": pc(br)[0], "rho_hi": pc(br)[1],
        "sd_y": y.std(ddof=1), "resid_sd": res.std(ddof=2), "r2": 1 - res.var() / y.var(),
        "loo_rmse_line": float(np.sqrt(np.mean(np.square(e_line)))),
        "loo_rmse_mean": float(np.sqrt(np.mean(np.square(e_mean)))),
        "q_min": q.min(), "q_max": q.max(),
    })
    return out


def main():
    log = open(f"{W}/charge_line.log", "w")

    def p(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        log.write(s + "\n")
        log.flush()

    rng = np.random.default_rng(SEED)
    d = pd.read_csv(f"{W}/la_pairs_all_metals.tsv", sep="\t")
    d["zL"] = d.prot.map(charge)
    d["q"] = -d.zL
    d["donor_class"] = d.name.map(donor_class)
    d["el"] = d.metal.str.split("/").str[0]
    p("rows", len(d), "unparsed charge", int(d.zL.isna().sum()))

    # keyword check against the feature table
    f = pd.read_csv(FEAT, sep="\t")
    f["n_true"] = (f[["aminoacid_N", "amine_N", "pyridine_N"]].sum(axis=1) > 0)
    f["cls"] = f.name.map(donor_class)
    chk = f[f.cls != "other"]
    agree = ((chk.cls == "N_containing") == chk.n_true).mean()
    p(f"keyword check on {len(chk)} featured ligands (class other excluded): agreement {agree:.3f}")
    bad = chk[(chk.cls == "N_containing") != chk.n_true]
    for _, r in bad.iterrows():
        p("  disagreement:", r.lig_id, r["name"], "keyword", r.cls, "features N", bool(r.n_true))

    rows = []
    same = d[d.match == "same_I"]
    for el in EARTHS + ["Mn", "Zn", "Cd", "Lu", "Y"]:
        s = same[same.el == el]
        for cls, sub in [("all", s), ("O_only", s[s.donor_class == "O_only"]),
                         ("N_containing", s[s.donor_class == "N_containing"])]:
            r = fit(sub.q, sub.y_log, rng)
            r.update({"target": f"La-{el}", "donor_class": cls})
            rows.append(r)

    # class-held-out, La-Ca
    ca = same[same.el == "Ca"]
    for train, test in [("N_containing", "O_only"), ("O_only", "N_containing")]:
        a, b = ca[ca.donor_class == train], ca[ca.donor_class == test]
        if len(a) >= 8 and len(b) >= 8 and a.q.nunique() >= 3:
            s1, s0 = np.polyfit(a.q, a.y_log, 1)
            e_line = b.y_log - (s0 + s1 * b.q)
            e_mean = b.y_log - a.y_log.mean()
            p(f"class-held-out La-Ca: train {train} (n={len(a)}) -> test {test} (n={len(b)}): "
              f"RMSE line {np.sqrt((e_line**2).mean()):.2f}, RMSE train-mean {np.sqrt((e_mean**2).mean()):.2f}, "
              f"mean error line {e_line.mean():+.2f}")

    # lanthanide-specific lines against Ca
    ca_i = ca.set_index("lig_id")
    for el in ["La"] + LN:
        if el == "La":
            t = ca.assign(y_x=ca.y_log)
        else:
            x = same[same.el == el].set_index("lig_id")
            j = ca_i.join(x[["cond_La", "logK_X"]].rename(columns={"cond_La": "cond_La_x", "logK_X": "logK_Ln"}), how="inner")
            j = j[j.cond_La == j.cond_La_x]
            t = j.assign(y_x=j.logK_Ln - j.logK_X).reset_index()
        for cls, sub in [("all", t), ("O_only", t[t.donor_class == "O_only"])]:
            r = fit(sub.q, sub.y_x, rng)
            r.update({"target": f"{el}-Ca", "donor_class": cls})
            rows.append(r)

    out = pd.DataFrame(rows)
    cols = ["target", "donor_class", "n", "slope", "slope_lo", "slope_hi", "intercept", "intercept_lo",
            "intercept_hi", "rho", "rho_lo", "rho_hi", "sd_y", "resid_sd", "r2", "loo_rmse_line",
            "loo_rmse_mean", "q_min", "q_max"]
    out = out[[c for c in cols if c in out.columns]]
    out.to_csv(f"{W}/charge_lines.tsv", sep="\t", index=False, float_format="%.3f")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    p(out.round(2).to_string(index=False))

    lig = ca[["lig_id", "name", "prot", "zL", "q", "donor_class", "has_features", "family", "I_La",
              "logK_La", "logK_X", "y_log"]].sort_values(["donor_class", "q", "y_log"])
    lig.to_csv(f"{W}/charge_line_ligands.tsv", sep="\t", index=False)
    p()
    p("La-Ca same-I by donor class and q (count, median, min, max):")
    p(ca.groupby(["donor_class", "q"]).y_log.agg(["count", "median", "min", "max"]).round(2).to_string())
    log.close()


if __name__ == "__main__":
    sys.exit(main())
