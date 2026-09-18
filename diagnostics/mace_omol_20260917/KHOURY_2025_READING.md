# Khoury et al. 2025: calcium-binding domains as lanthanide binders

Read on 2026-09-18 at Jacob's direct request. Main article and44-page SI read;
Figure5 inspected visually to confirm micromolar units. Seven author PDBs
retrieved and inventoried. The initial main-article notes below are retained;
the exact SI evidence and completed benchmark update supersede pending/proposed
wording. [Completed benchmark](KHOURY_BENCHMARK_REPORT.md):all22sites/44endpoints,
A0A7 and HEW5 supporting successes, RTX failure;6/9 grouped comparisons.

## Completed supplementary verification

EuropePMC package PMC12152738 retrieved through its supplementaryFiles API.
ZIP SHA256:4f3e1f74580bfc6e33d724945a9ccd1d62d6002f6dab4dd98bd26c221d863035.
SI PDF SHA256:a0f23f97777602861556d172dad826621d0ae76a54674a2d267c696dea4d605e.
All PDB hashes and downloads are retained in the reading workspace receipt.

TableS6 exact LaKd (95%CI):A0A7 17±2µM,HEW5 5.2±1.7µM,RTX40±4µM.
Ca≈2000,750,250µM are CD folding thresholds, not fittedKd. ITC:25°C,
30µMprotein,50mMMES/50mMNaCl,pH6. CD:25°C,60µMprotein,1mMMES/50mMKCl,pH6.
The supplied models omit the experimental N-terminal MPVP cleavage scar.
Model sequences and all source heavy coordinates were verified before scoring;
no new folding calculation, La geometry, or scar model was invented.

## Initial main-article reading

Paper: *Mining peptides for mining solutions: evaluation of calcium-binding
peptides for rare earth element separations*, Chemical Science 16,15333–15346.
DOI: https://doi.org/10.1039/D5SC02315G
Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC12152738/

## What it adds

Seven constructs: A0A7 (A0A7L4YJY0,S111–S196), HEW5
(A0A6P0HEW5,P14–G125), K3T(VN) and K3T(VV) (A0A1G7K3T7,
V253–N402 and V1–V63), HJH0 (A0A6P0HJH0,G52–R95), CaM(III,IV)
(M76–K148,1CFF reference), and RTX Block V (5CVW reference).
These are construct descriptions from the article, not yet sequence-verified
preparations. Table 1 lengths sometimes include sequence differences/additions;
use actual supplementary construct sequences before modeling.

HEW5, A0A7 and RTX have solution ITC measurements including La at pH6.
Figure5A visually gives approximate La Kd values of5,16 and40 micromolar,
respectively; these are figure readings, not table-derived exact measurements.
The assay fits average multiple sites. Initial low-ratio points were omitted;
for A0A7/RTX the authors explicitly say the first site is not represented.
No assignment of a measured affinity to each predicted metal site is justified.

Calcium ITC curves could not support Kd fits. Calcium comparison instead uses
concentrations inducing conformational change in CD. Reported approximately
60-fold(A0A7/HEW5) and15-fold(RTX) selectivity compares this calcium proxy with
average lanthanide Kd, not a direct fitted La/Ca affinity ratio. Keep as
qualified protein-level evidence and retain the actual conditions/readout.

The seven-construct FRET screen tests Ce,Nd,Dy,Yb at pH7.4, with CFP/EYFP
fusion constructs; it is not a La/Ca dataset. XO competition tests Yb at pH6
and counts high-affinity sites. Failure to compete effectively does not make
HJH0,CaM or RTX calcium-selective negatives.

The study searched soil microbial sequences with domain HMMs and deliberately
excluded more than19000 PQQ-domain hits. It is useful complementary context
for the hillslope paper, not validation of physiological Ln use in those soils.

## Mechanistic implication for our work

ITC supports net entropy-driven binding in the three down-selected domains;
CD shows disordered apo states and metal-induced structural ordering. Authors
attribute the favorable entropy partly to dehydration and hydrophobic burial.
The experiments do not uniquely separate those entropy contributions.

Inference: even an accurate electronic energy of a single bound structure
cannot be assumed to capture coupled hydration/folding contributions. This
strengthens the case for affordable structural/ensemble and hydration
information alongside electronic scoring. It does not establish the cause of
our GGR or alpha failures, validate a local harmonic entropy correction, or
justify routine long MD/FEP. A full unfolding/folding transition exceeds the
scope of our previously tested two-coordinate local response model.

Immobilized A0A7/HEW5 show a different within-Ln preference than solution.
Immobilized measurements also change pH/conditions, so this is not a clean
causal isolation of scaffold mechanics. Within-Ln preferences remain separate
from La/Ca discrimination.

## Initial recommendation (subsequently executed above)

Prioritize the exact HEW5,A0A7,RTX constructs as a proposed non-PQQ challenge
set with explicit protein-level evidence. Retrieve SI sequences, TableS6,
CD titrations, assay conditions and author-supplied AF3 models first. Preserve
all sites as ordered outputs; no best-site selection or multiplying biological
sample size by predicted site count. Predicted Ca-bound geometries are not
experimental coordinates and do not establish the La-bound geometry.

For method development, ask whether affordable structural response and
hydration descriptors add useful information beyond frozen-site energies.
Do not start an expensive new solvent/backend project based on this paper.
The initial reading itself launched no analysis. The subsequently declared
parallel author-domain benchmark is now complete; see the result linked above.

## Source

Local PDF: /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/d5sc02315g.pdf
SHA256: 9a9c22209962a33c1ffa2b6c76b9333d7cc0852d3c3102fc4e963f614407ad5c
