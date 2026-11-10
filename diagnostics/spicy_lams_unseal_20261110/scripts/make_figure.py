"""Summary figure for the first look at the unsealed lanmodulin ortholog profiles. Lead agent, Claude Fable."""
import json, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3de"
C1, C2, C3, CG = "#2a78d6", "#eb6834", "#1baf7a", "#a9a8a2"
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.edgecolor": "#b9b8b2", "axes.linewidth": 0.8, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "savefig.facecolor": SURF, "axes.spines.top": False, "axes.spines.right": False, "xtick.major.size": 0, "ytick.major.size": 0})
d, raw, C, mean, valid, qc = load()
t = pd.read_csv(ROOT / "results/ortholog_contrasts.tsv", sep="\t"); q = t[t.qc]
A = json.load(open(ROOT / "results/first_look_a.json")); B = json.load(open(ROOT / "results/first_look_b.json"))
env = pd.read_csv(ROOT / "results/environmental_homologs.tsv", sep="\t")
fig = plt.figure(figsize=(14.2, 8.6))
gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.20, left=0.075, right=0.985, top=0.835, bottom=0.08)
axA, axB, axC, axD = (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]))
def title(ax, letter, text, sub):
    ax.text(0, 1.115, f"{letter}  {text}", transform=ax.transAxes, fontsize=11.5, fontweight="bold", va="bottom")
    ax.text(0, 1.035, sub, transform=ax.transAxes, fontsize=9.3, color=INK2, va="bottom")
def style(ax, grid="y"):
    ax.grid(axis=grid, color=GRID, linewidth=0.8); ax.set_axisbelow(True)
# A: profiles by genus
order = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y"]; oi = [EL.index(e) for e in order]; x = np.arange(15)
M = mean[qc][:, oi]; gen = d.genus_original.fillna("unspecified").values[qc]; fam = d.family_original.fillna("unspecified").values[qc]
axA.fill_between(x[:14], np.percentile(M, 5, 0)[:14], np.percentile(M, 95, 0)[:14], color=CG, alpha=0.25, linewidth=0)
hA = []
for nm, sel, c in (("Methylobacterium (n=112)", gen == "Methylobacterium", C1), ("Bradyrhizobium (n=126)", gen == "Bradyrhizobium", C2), ("Xanthobacteraceae (n=5)", fam == "Xanthobacteraceae", C3)):
    m = M[sel].mean(0); axA.plot(x[:14], m[:14], color=c, lw=2, solid_capstyle="round"); axA.plot(14, m[14], marker="o", ms=8.5, mfc=c, mec=SURF, mew=1.6)
    hA.append(Line2D([], [], color=c, lw=2, marker="o", ms=7, mfc=c, mec=SURF, label=nm))
axA.legend(handles=hA, loc="upper right", bbox_to_anchor=(0.93, 1.0), frameon=False, fontsize=8.8)
axA.axvline(13.75, color="#b9b8b2", lw=0.8); axA.set_xticks(x); axA.set_xticklabels(order); axA.set_xlim(-0.4, 14.6)
axA.set_ylabel("captured, log10, relative to the protein's mean"); style(axA)
axA.text(0.2, -0.78, "grey band: 5th to 95th percentile of 353 orthologs", fontsize=8.6, color=INK2)
title(axA, "A", "Orthologs differ mainly by one tilt: light against heavy", "Mean capture profiles from an equal-mass mixture at pH 3; yttrium shown apart, at the right")
# B: tilt against loop acidity
jit = np.random.default_rng(1).normal(0, 0.035, len(q))
col = np.where(q.genus == "Methylobacterium", C1, np.where(q.genus == "Bradyrhizobium", C2, np.where(q.family == "Xanthobacteraceae", C3, CG)))
axB.scatter(q.acid + jit, q.Lu_La, s=16, c=col, linewidths=0.4, edgecolors=SURF, alpha=0.9)
axB.axhline(0, color="#8d8c86", lw=1.0); axB.set_xlabel("acidic residues per EF loop (positions 1, 3, 5, 9, 12)"); axB.set_ylabel("Lu minus La, log10"); style(axB)
axB.text(0.02, 0.93, "rank correlation +0.70 pooled, +0.56 within genus", transform=axB.transAxes, fontsize=9.2)
axB.text(0.985, 0.045, "above the line: prefers Lu to La (17 of 353)", transform=axB.transAxes, fontsize=8.6, color=INK2, ha="right")
axB.legend(handles=[Line2D([], [], marker="o", ls="none", ms=7, mfc=c, mec=SURF, label=l) for c, l in ((C1, "Methylobacterium"), (C2, "Bradyrhizobium"), (C3, "Xanthobacteraceae"), (CG, "other"))],
           loc="center left", bbox_to_anchor=(0.02, 0.52), frameon=False, fontsize=8.8, handletextpad=0.2)
title(axB, "B", "More acidic loops give a flatter, heavier-leaning profile", "Each point is one ortholog (improved loop finder; post hoc)")
# C: size of natural differences
names = ["Lu-La", "Dy-Nd", "Y-Dy", "Y-Ho", "Tb-Gd", "Dy-Tb", "Nd-Pr"]; sd = [A["A3"][k]["sd_true"] for k in names]; rep = [A["A3"][k]["sd_rep"] / np.sqrt(3) for k in names]
yy = np.arange(len(names))[::-1]
axC.barh(yy, sd, height=0.5, color=C1, linewidth=0)
for y0, v in zip(yy, sd): axC.text(v + 0.006, y0, f"{v:.3f}", va="center", fontsize=9)
axC.axvline(0.1, color="#8d8c86", lw=1.0); axC.text(0.106, 2.5, "0.1: the bar set in advance\nfor 'something to mine'", fontsize=8.6, color=INK2, va="center")
axC.set_yticks(yy); axC.set_yticklabels([n.replace("-", " minus ") for n in names]); axC.set_xlabel("spread across orthologs, log10 units (SD after removing replicate noise)"); axC.set_xlim(0, 0.42); style(axC, "x")
title(axC, "C", "Neighbouring elements are separated the same way by every ortholog", "Cross-branch contrasts vary; within-branch contrasts and yttrium hardly do")
# D: predictability by identity, with environmental homologs marked
bins = [("0.5-0.7", 0.60), ("0.7-0.85", 0.775), ("0.85-0.95", 0.90), ("0.95-1.01", 0.98)]
for k, c, lab in (("B2", C1, "five nearest library neighbours"), ("M2", C2, "language-model regression")):
    xs = [b[1] for b in bins if b[0] in B["A5_bins"]]; ys = [B["A5_bins"][b[0]][k] for b in bins if b[0] in B["A5_bins"]]
    axD.plot(xs, ys, color=c, lw=2, marker="o", ms=8.5, mfc=c, mec=SURF, mew=1.6, label=lab)
for b in bins:
    if b[0] in B["A5_bins"]: axD.text(b[1], -0.62, f"n={B['A5_bins'][b[0]]['n']}", ha="center", fontsize=8.6, color=INK2)
axD.axhline(0, color="#8d8c86", lw=1.0); axD.set_ylim(-0.75, 0.95); axD.set_xlim(0.28, 1.02)
axD.set_xlabel("identity to the nearest measured library protein"); axD.set_ylabel("prediction skill (1 = perfect, 0 = library mean)"); style(axD)
marks = [("Omnitrophota EF-hand\nfamily (42), median", float(env[env.name.str.startswith("omni")].nearest_identity.median())),
         ("SR-VP homologs", None)]
om = float(env[env.name.str.startswith("omni")].nearest_identity.median())
axD.annotate("Omnitrophota EF-hand\nfamily (42 proteins)", xy=(om, 0.0), xytext=(om, 0.52), ha="center", fontsize=8.8, color=INK, arrowprops=dict(arrowstyle="-", color="#8d8c86", lw=0.9))
sr = env[env.name.str.startswith("srvp")].nearest_identity.values
for v in sr: axD.plot([v, v], [0.0, 0.80], color="#8d8c86", lw=0.9)
axD.text(np.mean(sr[sr > 0.6]), 0.84, "five SR-VP homologs", ha="center", fontsize=8.8, color=INK)
axD.legend(loc="lower right", bbox_to_anchor=(1.0, 0.13), frameon=False, fontsize=8.8)
title(axD, "D", "A measured neighbour at 70% identity or more predicts the profile", "Ten-fold cross-validation within the 353 orthologs; vertical marks: Jacob's environmental homologs")
fig.text(0.06, 0.958, "616 natural lanmodulins: reproducible differences, one main axis, predictable from close relatives", fontsize=15, fontweight="bold")
fig.text(0.06, 0.928, "SpyCI-LAMBS capture profiles (Nat. Chem. Biol. 2026, Supplementary Data 1); 600 rows opened on 10 November 2026, 353 pass quality control; 16 reserved proteins still sealed",
         fontsize=10, color=INK2)
fig.savefig(ROOT / "figures/lanmodulin_ortholog_first_look.png", dpi=200); print("wrote figure")
