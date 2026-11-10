"""Summary figure for the Gd|Tb gap analysis. Reads saved results; recomputes only the per-break scan on the
24 confirmation ligands. Lead agent, Claude Fable."""
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from lib import *  # noqa

RES = ROOT / "results"; FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3de"
C47, C24, CPR = "#2a78d6", "#eb6834", "#1baf7a"      # fitted ligands, confirmation ligands, proteins and peptides
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.edgecolor": "#b9b8b2", "axes.linewidth": 0.8,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "axes.spines.top": False, "axes.spines.right": False, "xtick.major.size": 0, "ytick.major.size": 0})

step_fit = json.load(open(RES / "step_fit.json")); conf = json.load(open(RES / "confirm_summary.json"))
rob = json.load(open(RES / "robustness_summary.json")); fol = json.load(open(RES / "followup_summary.json"))
ybr = json.load(open(RES / "yttrium_branch.json"))
coord = pd.read_csv(RES / "coordinate_cubic.tsv", sep="\t"); boot = np.load(RES / "bootstrap_rho_KC.npy")
GAP = 0.031

fig = plt.figure(figsize=(14.2, 8.4))
gt = fig.add_gridspec(1, 2, wspace=0.16, left=0.055, right=0.985, top=0.815, bottom=0.535)
gb = fig.add_gridspec(1, 3, wspace=0.70, left=0.125, right=0.975, top=0.385, bottom=0.085)
axA = fig.add_subplot(gt[0, 0]); axB = fig.add_subplot(gt[0, 1])
axC = fig.add_subplot(gb[0, 0]); axD = fig.add_subplot(gb[0, 1]); axE = fig.add_subplot(gb[0, 2])

def style(ax, grid="y"):
    ax.grid(axis=grid, color=GRID, linewidth=0.8); ax.set_axisbelow(True)

def title(ax, letter, text, sub=None):
    ax.text(0, 1.115 if sub else 1.045, f"{letter}  {text}", transform=ax.transAxes, fontsize=11.5, fontweight="bold", color=INK, va="bottom")
    if sub:
        ax.text(0, 1.035, sub, transform=ax.transAxes, fontsize=9.3, color=INK2, va="bottom")

# ---------------- A: the fitted coordinate ----------------
x = np.arange(14); d = 1000 * (coord.rho_KC.values - R9)
rng = np.random.default_rng(1)
for b in boot[rng.choice(len(boot), 150, replace=False)]:
    db = 1000 * (b - R9)
    axA.plot(x[:7], db[:7], color=C47, lw=0.6, alpha=0.10); axA.plot(x[7:], db[7:], color=C47, lw=0.6, alpha=0.10)
# one-number model in the same gauge (ends pinned): light branch raised by the gap, then the straight line through the ends removed
sm = GAP * (np.arange(14) <= 6) - GAP * (R9 - R9[-1]) / (R9[0] - R9[-1])
axA.plot(x[:7], 1000 * sm[:7], color="#8d8c86", lw=2, solid_capstyle="round"); axA.plot(x[7:], 1000 * sm[7:], color="#8d8c86", lw=2, solid_capstyle="round")
axA.plot(x, d, ls="none", marker="o", ms=8.5, mfc=C47, mec=SURF, mew=1.6, zorder=5)
axA.axhline(0, color="#b9b8b2", lw=0.8)
axA.set_xticks(x); axA.set_xticklabels(LN); axA.set_ylabel("fitted position minus Shannon radius (mÅ)")
axA.set_ylim(-36, 36); style(axA)
axA.annotate("", xy=(6.5, -9.3), xytext=(6.5, 17.2), arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.0))
axA.text(6.72, 24.5, "extra gap between Gd and Tb:\n31 mÅ, about 2.2 element steps", fontsize=9.3, color=INK, va="center")
axA.text(1.0, 20.5, "La to Gd", fontsize=9.5, color=INK2); axA.text(10.0, -24.5, "Tb to Lu", fontsize=9.5, color=INK2)
axA.set_xlabel("points: 12 free positions fitted to 47 ligands; thin lines: 150 bootstrap refits; grey: the one-number gap model", fontsize=8.8)
title(axA, "A", "One size coordinate fits every ligand, and it has a gap",
      "Each ligand's series is a cubic in a shared coordinate; the ends are pinned to the Shannon nine-coordinate radii")

# ---------------- B: where the break belongs ----------------
pairs = [f"{LN[k]}|{LN[k+1]}" for k in range(1, 13)]
f0 = 6.14088; red47 = np.array([100 * (1 - step_fit["by_break"][p_][1] / f0) for p_ in pairs])
Yall, _, _ = load_nist(10); Y24 = Yall[np.isnan(Yall).any(1)]
def rss(rho, rows):
    B = poly_basis(rho, 3); t = 0.0
    for y in rows:
        ok = ~np.isnan(y); b = np.linalg.lstsq(B[ok], y[ok], rcond=None)[0]; t += float(((y[ok] - B[ok] @ b) ** 2).sum())
    return t
grid = np.arange(-0.060, 0.1501, 0.001); g0 = rss(R9, Y24)
red24 = np.array([100 * (1 - min(rss(R9 + dd * (np.arange(14) <= k), Y24) for dd in grid) / g0) for k in range(1, 13)])
xb = np.arange(12); w = 0.34
axB.bar(xb - w / 2 - 0.02, red47, w, color=C47, linewidth=0); axB.bar(xb + w / 2 + 0.02, red24, w, color=C24, linewidth=0)
axB.set_xticks(xb); axB.set_xticklabels([p_.replace("|", "·") for p_ in pairs], fontsize=9)
axB.set_ylabel("residual removed by one break here (%)"); style(axB); axB.set_ylim(0, 27)
axB.text(5, red24[5] + 0.9, "Gd·Tb", ha="center", fontsize=9.5, fontweight="bold", color=INK)
axB.annotate("La alone: Ce is absent\nfrom 22 of these 24 series", xy=(0.3, red24[0] - 0.4), xytext=(1.25, 21.3), fontsize=8.3, color=INK2,
             arrowprops=dict(arrowstyle="-", color="#b9b8b2", lw=0.8), va="center")
title(axB, "B", "The break belongs between Gd and Tb in both sets",
      "Best single break at each of the twelve neighbour pairs; cubic per ligand")

# ---------------- C: out-of-sample error ----------------
rowsC = [("plain cubic", fol["scores"]["K"]["S_LOMO"], conf["rms_K"]),
         ("cubic + gap\n(1 shared number)", step_fit["S_LOMO_KS"], conf["rms_KS"]),
         ("cubic, free coordinate\n(12 shared numbers)", fol["scores"]["KC"]["S_LOMO"], conf["rms_KC"])]
for i, (lab, a, b) in enumerate(rowsC):
    y0 = len(rowsC) - 1 - i
    for off, v, c in ((0.14, a, C47), (-0.14, b, C24)):
        axC.plot([v[1], v[2]], [y0 + off] * 2, color=c, lw=2, solid_capstyle="round")
        axC.plot(v[0], y0 + off, marker="o", ms=8.5, mfc=c, mec=SURF, mew=1.6)
axC.set_yticks(range(len(rowsC))); axC.set_yticklabels([r[0] for r in rowsC][::-1], fontsize=9.3)
axC.set_xlabel("leave-one-metal-out error, kcal/mol (95% interval)"); axC.set_xlim(0.07, 0.29); axC.set_ylim(-0.55, 2.55); style(axC, "x")
title(axC, "C", "The gap predicts held-out metals", "Paired gain: 0.020 [0.009, 0.032] (47); 0.075 [0.050, 0.099] (24)")

# ---------------- D: gap by data source ----------------
rowsD = [("47 fitted ligands", [31.0, 18.5, 61.5], C47), ("24 confirmation ligands", [1000 * conf["gap_24"]["gap"], 1000 * conf["gap_24"]["lo"], 1000 * conf["gap_24"]["hi"]], C24),
         ("12 protein and peptide\ndatasets, pooled", [1000 * v for v in rob["C_pooled_all 12_deg3"]], CPR),
         ("lanmodulin, 14 metals\n(one dataset)", [5.0, -3.5, 16.5], CPR)]
for i, (lab, v, c) in enumerate(rowsD):
    y0 = len(rowsD) - 1 - i
    axD.plot([v[1], v[2]], [y0, y0], color=c, lw=2, solid_capstyle="round"); axD.plot(v[0], y0, marker="o", ms=8.5, mfc=c, mec=SURF, mew=1.6)
    axD.text(v[2] + 3, y0, f"{v[0]:+.0f}", va="center", fontsize=9, color=INK)
axD.axvline(0, color="#8d8c86", lw=1.0)
axD.set_yticks(range(len(rowsD))); axD.set_yticklabels([r[0] for r in rowsD][::-1], fontsize=9.3)
axD.set_xlabel("gap between Gd and Tb (mÅ, 95% interval)"); axD.set_xlim(-30, 95); axD.set_ylim(-0.55, 3.55); style(axD, "x")
title(axD, "D", "Proteins and peptides show no gap", "Gap fitted separately to each source")

# ---------------- E: yttrium ----------------
rowsE = [("oxygen-only ligands", ybr["best_frac_O-only"], C47, 0.14), ("oxygen-only ligands", conf["Y_best_frac_O-only"], C24, -0.14),
         ("nitrogen-containing\nligands", ybr["best_frac_N-containing"], C47, 0.0)]
ypos = {"oxygen-only ligands": 1, "nitrogen-containing\nligands": 0}
for lab, v, c, off in rowsE:
    y0 = ypos[lab] + off
    hi = min(v[2], 2.05)
    axE.plot([v[1], hi], [y0, y0], color=c, lw=2, solid_capstyle="round"); axE.plot(v[0], y0, marker="o", ms=8.5, mfc=c, mec=SURF, mew=1.6)
    if v[2] > 2.05:
        axE.annotate("", xy=(2.13, y0), xytext=(2.0, y0), arrowprops=dict(arrowstyle="->", color=c, lw=2))
        axE.text(2.02, y0 + 0.14, f"to {v[2]:.1f}", fontsize=8.3, color=INK2, ha="center")
axE.axvline(0, color="#8d8c86", lw=1.0); axE.axvline(1, color="#8d8c86", lw=1.0)
axE.set_yticks([0, 1]); axE.set_yticklabels(["nitrogen-containing\nligands", "oxygen-only ligands"], fontsize=9.3)
axE.set_xlabel("yttrium's position, fraction of the gap (95%)"); axE.set_xlim(-0.75, 2.2); axE.set_ylim(-0.55, 1.55); style(axE, "x")
title(axE, "E", "Yttrium crosses only with O donors", "0 = with Dy and Ho (its own size); 1 = La–Gd side")

fig.legend(handles=[Line2D([], [], marker="o", ls="none", ms=8.5, mfc=C47, mec=SURF, mew=1.6, label="47 ligands with complete series (fitted)"),
                    Line2D([], [], marker="o", ls="none", ms=8.5, mfc=C24, mec=SURF, mew=1.6, label="24 further ligands (confirmation, untouched until the test was frozen)"),
                    Line2D([], [], marker="o", ls="none", ms=8.5, mfc=CPR, mec=SURF, mew=1.6, label="protein and peptide series")],
           loc="upper left", bbox_to_anchor=(0.045, 0.928), ncol=3, frameon=False, fontsize=9.6, handletextpad=0.3, columnspacing=2.2)
fig.text(0.055, 0.958, "Lanthanide stability constants fall on two branches, split at the half-filled 4f shell", fontsize=15, fontweight="bold", color=INK)
fig.text(0.055, 0.929, "NIST SRD 46 log K1 for La–Lu (energies in kcal/mol at 298 K); analysis of 10 November 2026", fontsize=10, color=INK2)
fig.savefig(FIG / "lanthanide_two_branches.png", dpi=200)
fig.savefig(FIG / "lanthanide_two_branches.pdf")
print("wrote", FIG / "lanthanide_two_branches.png")
print("break scan, % residual removed: 47:", np.round(red47, 1), " 24:", np.round(red24, 1))
