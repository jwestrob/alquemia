# Global MACE development panel: five real structures

Approved under Jacob's active MACE discriminator goal and full discretionary
pilot authorization. Declared before preparing or scoring these new global
states. Preserve baseline, original carves, labels and historical experiments.

## Question and input selection

Does a whole-protein MACE description, with a matched frozen-monopole OBC-II
solvent contribution, recover useful La/Ca ordering across different proteins?
Single-system medium/large agreement cannot answer this. The preceding solver
check passed numerical gates and reduced checkpoint disagreement923.9→106.4
kcal/mol on1H4I; that remaining discrepancy motivates comparison, not rejection
or a successful-classification claim.

Use exactly the five source structures pinned in
../global_electrostatic_20260916/ACCURACY_INPUTS.json:

| Case | Evidence stratum | Expected relative direction | Group |
|---|---|---|---|
| PQQ_1H4I | Canonical PQQ functional metal association | Lower than4MAE | MxaF |
| PQQ_4MAE | Canonical PQQ functional metal association | Higher than1H4I | XoxF |
| GGR_1GLG | Direct same-assay Ca-over-La affinity | Lower than both alpha structures | GGR |
| ALPHA_1F6S | Qualified qualitative La-over-Ca strong-site evidence | Higher thanGGR | bovine alpha |
| ALPHA_6IP9 | Same evidence, different structure | Higher thanGGR | bovine alpha |

All are consumed development cases. The alpha pair is one biological
observation; no independence or blind-test claim. Cross-protein directions are
qualitative development checks, not common-condition quantitative affinities.
No new experimental label is introduced. Primary evidence and caveats are
retained in ../accuracy_strategy_20260915/CHALLENGE_PANEL_CURATION.md.

## Common whole-protein preparation

Use the recorded deposited catalytic/site chainA, not an inferred biological
assembly or every crystallographic copy. Keep every standard-protein atom in
that chain and preserve every existing heavy coordinate. Explicitly inventory
excluded chains and nonprotein source residues. Keep the selected metal at its
frozen source location. Retain the exact frozen PQQ(3−) atom coordinates and
microstate for the two PQQ enzymes. Retain only the frozen site waters:0 GGR,
2 1F6S,3 6IP9,0 for the PQQ pair. Use their exact prepared neutral-water XYZs.
Other crystallographic waters, GGR's remote GAL,6IP9's remote sulfate/glycerol,
and the already excluded4MAE15P adduct remain outside this defined model.
These are explicit state assumptions, not evidence that those species never
matter experimentally. No state is altered in response to the score.

The canonical PQQ sources lack terminalOXT. Complete only that terminal
carboxylate: reflect the existing C→O vector about the C→CA direction in the
local peptide plane, set C–OXT to1.25A (ff19SB carboxylate equilibrium), retain
existing atoms. Require an actual free C terminus and no missing other heavy
atoms. Record the new atom, distance from metal and formal terminal state;
do not bridge a missing residue or model an entire missing tail. This uniform
terminal-completion policy is applied only where needed, before scoring.
Unsupported topology fails explicitly.

Apply the already tested ff19SB protein-H radial bond-length projection to
all protein H atoms, preserving their directions and attachment identities.
Do not move PQQ/water atoms or any source heavy atom. Standard-protein template
assignment supplies the integral protein formal charge from its explicit
protonation/disulfide state. MACE predicts all atomic charges; no ff19SB atom
charges enter the MACE or GB energy. TotalQ=Qprotein+QPQQ+Qmetal, watersneutral,
La+3/Ca+2, singlet parity checked. La/Ca share exact nonmetal atoms/coordinates,
cavity radii, assembly and chemical state. No caps exist in the full system.

## Scores and fixed comparisons

Primary candidate: analytic-kernel MACE-POLAR-1-M vacuum+frozen-monopole OBC-II.
Report vacuum MACE separately. MACE-POLAR-1-L is a declared model-sensitivity
comparison, not a per-protein alternative selected for the better answer.
Use both checkpoint/software pins already tested; no changed weights or
functional. Same OBC-II settings and approximation limits as the preceding
mace_gb_20260916/PLAN.md, with the validated isolated CUDA12 solver build.

Rvac=E_Ca−E_La. Rsolv=Rvac+G_Ca−G_La. No CPCM is included and no baseline
aquo offset/band is inherited. No calibrated absolute decision yet. Report
4MAE−1H4I and each alpha−GGR contrast, larger expected for the La-associated
member. Keep all values, failures and structural replicates; no post hoc
averaging or site selection. The two alpha results remain an ordered vector.

Predeclared usefulness screen: all three qualitative contrasts positive in
the primary candidate, plus numerical admission checks. This is an initial
research screen only, not goal completion or a validated general classifier.
Negative or checkpoint-dependent outcomes guide an explicitly versioned next
model. Do not fit thresholds to these five cases. Full canonical calibration
and independent composition-challenging evidence are later requirements.

## Finite runs, numerical checks and cost

Five cases×two endpoints×two checkpoints=20 primary MACE energy/force calls.
For medium additionally rotate both endpoints of GGR and4MAE using the existing
37-degree[1,2,3] rotation around metal:4 calls. Total24 new MACE calls,
zero DFT calls. Save learned densities, all component energies and forces.
Reuse validated numerical kernel tests and original medium/large core gates.
Full-geometry preparation is deterministic and tested on these actual inputs.
Rotation gates unchanged:0.01kcal/mol energy/contrast,0.001eV/A force;
charge closure1e-5e. Initial geometry identity/charge/template failures block
that case and remain in the denominator.

Then24 native GB reaction-energy/force calls on those exact saved densities,
including the four rotated cases. Reuse the passed native/custom solver gate;
no repeated scientific GB tuning or new model variants. Inspect paired radius,
geometry, charge and rotation invariants. GB fixed-charge forces are not a
combined descriptor gradient; response/relaxation is unavailable.

One A5000/16CPU/64474MiB per serial model job, reuse existing runner/receipts.
Recent9141-atom cost≈117s medium/244s large per pair; the smaller proteins
should be cheaper. Initial expectation: tens of GPU minutes for the entire
panel, not hours per site. Record actual costs and failed attempts. No artificial
project CPU/time budget. No shared environment install, baseline change,
queue reprioritization, production rescore or automatic promotion.

## Pre-scoring water-H preparation correction

Before any global preparation/output, inspection of the five retained alpha
waters found O–H lengths1.163–1.190A, angles105.34–105.87degrees. As with the
protein H defect, these are the archived protonation coordinates, not measured
water H positions. Extend radial H bond projection to these waters using the
installed TIP3P equilibrium O–H length0.9572A; preserve oxygen positions, each
bond direction, molecule identity, protonation, water count and angle. Do not
move water O, rotate water molecules, add waters or change a metal endpoint.
Pin the actual TIP3P parameter source. This supersedes only the instruction
above to retain exact water-H coordinates, before any new score. PQQ coordinates
remain exact. All original baseline water coordinates remain archived unchanged.
The24MACE/24GB task inventory and all score/numerical criteria are unchanged.
