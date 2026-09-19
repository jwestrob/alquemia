# Context-supported physical coordination proposals

## Scope and scientific distinction

Jacob explicitly requested a delegated test that MACE proposes physically
protein-constrained metal-specific structures, and native DFT checks usefulness
before response becomes classifier features. Root requested the exact plan before
new execution; exact five-case proposal was sent for review on 2026-09-19.
Preparation and geometry-only verification can proceed before that review.

Consumed development cases: 1H4I MxaF, 4MAE XoxF, GGR 1GLG, alpha-lactalbumin
1F6S and 6IP9. Two alpha structures represent one biological group. PQQ functional
association and the condition-qualified alpha/GGR affinity comparison are separate
evidence strata. No untouched test set or new independent validation is claimed.

The failed 2026-09-18 coupled path already used physical sidechain torsions and
peptide crankshafts. Its energy was vacuum DFT plus a masked whole-minus-core
MACE scalar, and its path was one initial hybrid-gradient direction. This pilot
instead iteratively minimizes native unmasked OMOL within the completed compact
104–202-atom second-shell contexts. Its final score is ONLY unchanged original-core
r2SCAN-3c/CPCM(Water)/DefGrid3 singlepoint energy. There is no MACE scalar addition,
GB, tangent DFT anchor, entropy, or fitted label-dependent selection.

## Frozen preparation and numerical policy

- Ten endpoints: both Ca/La for the five listed structures.
- PQQ starts from the actual contextual physical-H pilot proposals and their
  completed original-core DFT checks. Alpha/GGR starts from the prior fixed-context
  pilot's actual original cores, including the improved alpha water H.
- Native OMOL100M checkpoint/software are reused exactly from the fixed-context
  pilot. All original core atoms, context membership, charge, protonation, PQQ
  microstate and explicit water identities/counts remain fixed.
- Three metal translations; original acidic/amide donor sidechain chi torsions;
  actual backbone-donor peptide crankshafts. Expected dimensions 10,12,16,11,11.
  Source connectivity identifies neighbors. Original PQQ and Arg, added context,
  scaffold anchors and all explicit waters are frozen. H follows its covalent
  group; synthetic caps follow both real bond anchors and have no independent DOF.
- Every physical heavy atom moves at most 0.20 Å; each angle at most 0.20 rad.
  The physical restriction applies also to moving source atoms outside the scored
  core. Covalent bond lengths must remain within 1e-9 Å of their starting values;
  analytic coordinate Jacobians must match geometry-only checks within 1e-7.
- One fixed starting structure, SLSQP, at most 200 algorithm iterations,
  ftol 1e-9 eV. A fixed starting-energy constant is subtracted during optimization
  for numerical precision, with all reported energies retained in full.
- Select the final feasible, non-energy-increasing iterate. Domain/chemistry
  failures remain unscorable; optimizer nonconvergence remains explicit. Report
  the tangent-cone projected gradient in physical-coordinate units. No claim of
  an unconstrained minimum or basin free energy follows from an accepted proposal.

## Actual useful output required

Ten compact-context searches, ten final original-core native OMOL evaluations,
and up to ten original-core native DFT singlepoints. No DFT gradients/optimizations.
Compare actual before/after PQQ gap, both alpha−GGR gaps and alpha replica spread;
lower energy alone does not establish improvement. Changed geometries do not
inherit original calibrated bands. Baseline/default and historical artifacts stay
unchanged. Preserve unsuccessful cases and denominators.

Expected allocation: one GPU with 16 CPUs/200000 MiB for a few minutes; four
concurrent 16-rank original-core DFT tasks on 64 CPUs/256 GiB, roughly 5–10 minutes
based on recent actual core receipts. This is a finite manifest, not an arbitrary
total compute budget. Actual Slurm costs and all scientific attempts are reported.

Relevant prior results: `mace_site_response_20260918/COUPLED_PATH_REPORT.md`,
`mace_mechanics_20260916/REPORT.md`, and `second_shell_20260919/{REPORT,PQQ_H_REPORT}.md`.
