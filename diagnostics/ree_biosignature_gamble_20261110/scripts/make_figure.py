#!/usr/bin/env python3
"""Figure for REPORT.md: (A) what the methane plume and symbiont-bearing gills removed, against three binder
predictions; (B) the two fingerprint numbers in pure solution and in carbonate-bearing water."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA, GREY, LGREY = "#2a78d6", "#eb6834", "#1baf7a", "#52514e", "#a9a8a2"
plt.rcParams.update({"font.family": "Helvetica Neue", "font.size": 10.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.linewidth": 0.8, "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})
PV = pd.read_csv(ROOT / "results/tests23/plume_vs_predictions.tsv", sep="\t").set_index("element")
MS = pd.read_csv(ROOT / "results/tests23/mussel_gill_selectivity.tsv", sep="\t"); MS = MS[MS.symbionts == "methanotrophs"]
FP = pd.read_csv(ROOT / "results/tests23/natural_and_reference_fingerprints.tsv", sep="\t"); FP = FP[FP.variant == "no_Ce_La"]
RF = pd.read_csv(ROOT / "results/reference_fingerprints.tsv", sep="\t"); RF = RF[RF.variant == "no_Ce_La"]

fig, (ax, bx) = plt.subplots(1, 2, figsize=(13.2, 5.6), gridspec_kw=dict(width_ratios=[1.12, 1], wspace=0.24))

# ---------------- panel A
els = ["La", "Pr", "Nd", "Sm", "Eu", "Gd", "Dy", "Yb"]; x = np.arange(len(els))
def series(col):
    return np.array([1.0] + [PV.loc[e, col] for e in els[1:]])
lines = [("iron oxide, carbonate solution (measured, median of 5 runs)", GREY, "-", "iron oxide in carbonate water (measured)"),
         ("lanmodulin x free-ion fraction", AQUA, "-", "lanmodulin capture x free-ion fraction"),
         ("flat binder x free-ion fraction", LGREY, (0, (4, 2)), "free-ion availability alone")]
for col, c, ls, lab in lines:
    ax.plot(x, series(col), color=c, lw=2, ls=ls, label=lab, zorder=2, solid_capstyle="round")
obs = series("observed"); lo = np.array([1.0] + [PV.loc[e, "lo"] for e in els[1:]]); hi = np.array([1.0] + [PV.loc[e, "hi"] for e in els[1:]])
ax.errorbar(x[1:], obs[1:], yerr=[obs[1:] - lo[1:], hi[1:] - obs[1:]], fmt="none", ecolor=BLUE, elinewidth=1.6, capsize=3, zorder=4)
ax.plot(x, obs, color=BLUE, lw=2, marker="o", ms=8, mec=SURF, mew=1.6, label="methane plume (95% interval)", zorder=5)
gm = {"Pr": "k_Pr", "Nd": "k_Nd", "Sm": "k_Sm", "Eu": "k_Eu", "Dy": "k_Dy"}
gx = [els.index(e) for e in gm]; gmed = [MS[c].median() for c in gm.values()]
gq1 = [MS[c].quantile(.25) for c in gm.values()]; gq3 = [MS[c].quantile(.75) for c in gm.values()]
ax.errorbar(np.array(gx) + 0.13, gmed, yerr=[np.array(gmed) - gq1, np.array(gq3) - gmed], fmt="none", ecolor=ORANGE, elinewidth=1.6, capsize=3, zorder=4)
ax.plot(np.array(gx) + 0.13, gmed, color=ORANGE, lw=0, marker="D", ms=7.5, mec=SURF, mew=1.4, label="gills with methane-eating symbionts (14 mussels; quartiles)", zorder=5)
ax.axhline(0, color=INK2, lw=0.8, zorder=1)
ax.set_xticks(x); ax.set_xticklabels(els); ax.set_ylim(-0.08, 1.08); ax.set_xlim(-0.35, len(els) - 0.55)
ax.set_ylabel("removal rate relative to lanthanum"); ax.set_xlabel("element (large ions at left, small at right)")
ax.grid(axis="y", color=GRID, lw=0.8, zorder=0); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
lab = dict(color=INK, fontsize=9.5, va="center")
ax.text(4.12, 0.475, "iron oxide", ha="left", **lab)
ax.text(3.05, 0.262, "lanmodulin", ha="left", **lab)
ax.text(3.05, 0.125, "free ions only", ha="left", **lab)
ax.text(0.80, 0.215, "methane\nplume", ha="right", linespacing=1.05, **lab)
ax.text(0.98, 0.105, "mussel gills", ha="right", **lab)
ax.set_title("A  Biology removed lanthanum to neodymium, then stopped", loc="left", fontsize=12, color=INK, pad=10, fontweight="bold")
ax.legend(loc="upper right", frameon=False, fontsize=8.6, labelcolor=INK2, handlelength=2.2, borderaxespad=0.0, bbox_to_anchor=(1.0, 1.0))

# ---------------- panel B
bx.add_patch(Rectangle((-1, -0.3), 2, 0.2, fc="#efeee9", ec=INK2, lw=0.8, ls=(0, (3, 2)), zorder=0))
bx.text(-1.10, -0.20, "values predicted\nfor protein sites", fontsize=8.8, color=INK2, va="center", ha="right")
lm = RF[RF.cls == "lanmodulin library"]
bx.scatter(lm.u, lm.Y_minus_Ho, s=7, color=LGREY, alpha=0.55, lw=0, zorder=1, label="353 lanmodulins, pure solution")
ch = RF[(RF.cls == "chelator") & (RF.donors == "O only")]
bx.scatter(ch.u, ch.Y_minus_Ho, s=30, marker="s", fc="none", ec=GREY, lw=1.0, zorder=2, label="oxygen-donor chelators, pure solution")
on = FP[FP.cls == "iron oxide, NaCl solution"]
bx.scatter(on.u, on.Y_minus_Ho, s=40, marker="^", fc="none", ec=ORANGE, lw=1.3, zorder=3, label="iron oxide, salt solution (16 runs)")
oc = FP[FP.cls == "iron oxide, carbonate solution"]
bx.scatter(oc.u, oc.Y_minus_Ho, s=62, marker="^", color=ORANGE, ec=SURF, lw=1.2, zorder=5, label="iron oxide, carbonate solution (5 runs)")
ls_ = FP[FP.cls == "lanmodulin in seawater (computed)"]
bx.scatter(ls_.u, ls_.Y_minus_Ho, s=58, marker="o", color=AQUA, ec=SURF, lw=1.2, zorder=5, label="lanmodulin in seawater (computed, 6 settings)")
gl = FP[FP.cls == "gill over deep water, methanotrophs"]
bx.scatter(gl.u, gl.Y_minus_Ho, s=58, marker="D", color=BLUE, ec=SURF, lw=1.2, zorder=5, label="gills with methane-eating symbionts (4 mussels)")
ar = dict(arrowstyle="-|>", color=INK2, lw=1.0, shrinkA=4, shrinkB=6, mutation_scale=11)
bx.annotate("", xy=(float(ls_.u.median()), float(ls_.Y_minus_Ho.median())), xytext=(float(lm.u.median()), float(lm.Y_minus_Ho.median())), arrowprops=ar, zorder=4)
bx.annotate("", xy=(float(oc.u.median()), float(oc.Y_minus_Ho.median())), xytext=(float(on.u.median()), float(on.Y_minus_Ho.median())), arrowprops=ar, zorder=4)
bx.text(-0.02, 0.05, "lanmodulin moves out", fontsize=9.2, color=INK, ha="right", va="center")
bx.text(-0.55, -0.455, "iron oxide moves in", fontsize=9.2, color=INK, ha="left", va="center")
bx.axhline(0, color=INK2, lw=0.6, zorder=0); bx.axvline(0, color=INK2, lw=0.6, zorder=0)
bx.set_xlim(-3.1, 3.4); bx.set_ylim(-0.6, 0.44)
bx.set_xlabel("step between Gd and Tb (element spacings)"); bx.set_ylabel("yttrium relative to holmium (log10)")
bx.grid(color=GRID, lw=0.8, zorder=0); bx.set_axisbelow(True)
for s in ("top", "right"): bx.spines[s].set_visible(False)
bx.set_title("B  Carbonate in the water swaps the two marks", loc="left", fontsize=12, color=INK, pad=10, fontweight="bold")
bx.legend(loc="upper left", frameon=False, fontsize=8.4, labelcolor=INK2, borderaxespad=0.1, handletextpad=0.4, ncol=1)
fig.savefig(ROOT / "figures/biosignature_gamble.png", dpi=200, bbox_inches="tight")
print("chelators off panel B:", int(((ch.u < -3.1) | (ch.u > 3.4) | ch.u.isna()).sum()), "of", len(ch), "| lanmodulins off:", int(((lm.u < -3.1) | (lm.u > 3.4)).sum()),
      "| gills:", gl[["u", "Y_minus_Ho"]].round(2).values.tolist())
