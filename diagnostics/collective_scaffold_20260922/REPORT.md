# Collective scaffold proposals: all 24 feasible, none passes stationarity

**The complete protein maps now generate real collective pocket proposals for all
eight declared sources.** The corrected run returns all 24 donor-matched targets
with preserved chemistry and physical guards. Every search stops at the existing
0.8 Å heavy-movement limit. No endpoint passes the declared stationarity test.
The motions are dominated by a common protein response also present with the
original donor positions; lower parent energy alone is not classifier utility.

This report covers mechanics. Root separately scores the corrected candidates
with the unchanged compact MACE/GFN2 expression under the
[approved round](../scaffold_restart_round_20260922/PLAN.md). No FF energy is added
to a score, and no prediction improvement is inferred here before that scoring.
Production/defaults, old scores, labels and calibrations are unchanged.

## What actually ran

The final protocol is
`source_ff19SB_collective_matched_donor_target_proposal_v2`, executed by job
**1210206** from `workspaces/collective_scaffold_20260922/searches_v3/`.

Independent parent preparation supports all eight source-exact standard-protein
ff19SB Systems, including the previously recorded distant crystal OXT completions.
The new 8 Å rule selects 1,941–2,099 atoms in 135–146 complete residues; each
search frees 1,938–2,094 atoms after fixing its declared donor oxygens. Direct
source peptide bonds define neighbor closure. Metal, PQQ, waters, proton inventory,
exterior atoms, donor target positions and context membership remain fixed.

Each source has exactly three single-start searches: original geometry, actual
Ca-adaptive target and actual La-adaptive target. Sparse projected L-BFGS with
source-bond retraction follows the frozen [execution plan](EXECUTION_PLAN.md).
OpenMM 8.5.1 uses the actual **OpenCL / NVIDIA H200 / double-precision** backend.
All 2,227 energy/force requests complete, and all 24 final contexts reconstruct
from physical parent coordinates and the unchanged cap-offset recipe. Each metal
receives the same nuclear geometry; endpoint charge/multiplicity remain unchanged.

There are no standalone cap motions, H bond-length normalization, inferred missing
atoms, new water/cofactor states, FF metal/PQQ parameters, Hessians or entropy terms.
The large protein parent is a proposal generator only. Its missing cofactor forces
and unsolvated protein mechanics remain limitations of this generator.

## Physical result and what limits it

All final geometries pass source bond, fixed-target, heavy displacement, actual
stereocenter, peptide cis/trans and overlap checks. All rejected physical trials
in the corrected run are heavy-displacement violations; none fails an overlap,
true-chirality or peptide-basin check. The optimizer then exhausts its line search
along the chosen direction. This does **not** prove that no other feasible
collective direction exists or that the pocket is a constrained minimum.

Projected maximum atom gradients remain **7.73–38.13 kcal/mol/Å**, above the
predeclared 0.1 limit. All 24 proposals are admitted by the separate finite,
feasible, nonincreasing-energy rule, retaining their nonstationary status.

| Source | Moving atom at the limit | In compact context? |
|---|---|---|
| 1H4I, all targets | Asp352 OD2 | No |
| 4MAE, all targets | Asp325 OD1 | No |
| Q9Z4J7, all targets | Asp377 OD2 | No |
| Q88JH5, all targets | Asp385 OD2 | No |
| A0A3F2YLY8 Ca-sample1, original/Ca target | Gly358 O | No |
| Same source, La target | Gly196 O | Yes |
| A0A3F2YLY8 Ca-sample3, all targets | Asp353 OD2 | No |
| A0ACD6B9F2 Ca-sample4, original target | Asp203 OD2 | Yes |
| Same source, Ca/La target | Ser513 OG | No |
| A0ACD6B9F2 La-sample4, all targets | Ser513 OG | No |

Several adaptive starting donor oxygens already lie at the archived movement
boundary. They remain fixed and are distinct from the moving atoms that stop
these searches. Thus this is mostly a scaffold-bound limitation, not another
change to the original donor-torsion bound. It still does not identify a unique
missing physical mechanism.

## Matched original-target control

Compare each adaptive target's response vector (final minus its own start) with
the original-target response over the same free protein heavy atoms. Cosines are
**0.97066–0.99893**: the overall movements are very similar. Target-dependent
differences remain, especially for A0ACD6B9F2 Ca-sample4 under the La target.
This is a descriptor of actual displacements, not a fitted accuracy model.

Parent FF work relative to each target's own initial state is roughly
−3,944 to −8,889 kcal/mol. These large common decreases characterize response to
the chosen protein-only potential and starting coordinates; they are not binding
energies or a justified scaffold correction. Actual source bond lengths, including
X–H, stay fixed, so the run does not relax the known H bond-length convention.
No unique energy component is assigned as the cause without a decomposition.

The original-target control is therefore essential: improved scores at the
metal-specific targets must be interpreted alongside generic protein relaxation.
The full parent vectors, target-specific work and numerical residuals are pinned
in the [result](RESULT.json)'s mechanical summary. Parent work is never used to
choose a metal preference or replace a failed compact score.

## Two preserved technical recoveries

1. Job **1210189**: all 24 CUDA Context initializations fail with unsupported PTX;
   **zero** molecular energy/force requests. The existing project OpenCL route
   supplies the explicit [backend recovery](OPENCL_RECOVERY.md), with no install
   or library substitution. Actual double precision is recorded; OpenCL exposes
   no separate `DeterministicForces` toggle, so none is claimed.
2. Job **1210194**: 942 force calls finish, but the first implementation falsely
   preserves signed volumes at achiral CH2 carbons. Concrete early obstructions
   are Lys387 CE in 4MAE and Pro440 CG in Q9Z4J7. Root approved the uniform
   [chemical-identity correction](STEREOCHEMISTRY_RECOVERY.md) before any candidate
   classification. All 24 old attempts remain immutable and **unscored**. The new
   protocol protects graph-verified non-Gly CA and Ile/Thr CB source handedness;
   no other numerical or physical bound changes.

Correcting that defect lets 4MAE and Q9 responses continue to the real heavy
movement limit. It does not establish successful relaxation or discriminatory
gain. No failed or superseded attempt is omitted from cost.

## Cost, checks and usable output

| Job | Status/use | Wall seconds | Allocated core-seconds | Requested GPU-seconds | Force calls |
|---|---|---:|---:|---:|---:|
| 1210189 | CUDA initialization failed | 39 | 1,248 | 39 | 0 |
| 1210194 | Superseded achiral guard | 375 | 12,000 | 375 | 942 |
| 1210206 | Corrected 24 proposals | 509 | 16,288 | 509 | 2,227 |
| **Total** | | **923** | **29,536** | **923** | **3,169** |

The corrected searches themselves total 486.558 seconds; actual force-evaluation
and compressed-record work totals 62.315 seconds. Geometry checks, projection,
retraction and line searches account for much of the rest. Allocation cost
includes idle reserved resources and collection. Local preparation/tests/reporting
are additional unmetered work; downstream root-owned compact scoring is separate.
GPU seconds describe requested allocation, not measured kernel utilization. No
GPU peak-memory measurement was collected; it is unavailable rather than zero.

**Eleven actual-fixture tests pass, zero skips**, including real Jacobian and
bond-retraction checks, all 24 source/target replays, corrupted real input failures,
achiral-site exclusion and protected actual stereocenters. These parser/geometry
tests are distinct from the 24 actually executed molecular searches. Root's
independent candidate-admission tests additionally replay corrected real receipts.

[Commands](COMMANDS.md) validate, collect and reproduce the read-only summary.
Every input, actual force evaluation, output geometry, failure and receipt stays
under `workspaces/`; [COST.json](COST.json) and [RESULT.json](RESULT.json) pin the
compact delivery. Root receives full parent Cartesian coordinates and reconstructed
contexts, without pretending collective motion is an old four-angle `full_q`.

**Recommendation:** assess these corrected finite proposals with the already
authorized matched compact scoring and original-target control. Do not add the
parent FF energy, expand the movement box, claim minima, or promote a classifier
from mechanical energy lowering. Further mechanics changes need an identified
benefit or limitation from that paired result.
