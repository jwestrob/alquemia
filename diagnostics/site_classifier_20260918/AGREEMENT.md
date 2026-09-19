# Approved: a small task-specific classifier from archived calculations

Jacob approved the immediately preceding proposal with: “I agree. Vamonos.”
The proposal retained DFT contrast as an input, added donor chemistry/denticity,
coordination geometry, solvent exposure and saved local MACE response, fitted
small constrained classifiers, and compared DFT alone / DFT+structure /
DFT+structure+MACE with related proteins held together. PQQ functional association
and direct affinity direction are separate targets. This is an accuracy experiment,
not further runtime optimization or default promotion.

## Frozen first implementation (before any fitted classifier result)

- Reuse the 25 canonical PQQ preparations and three archived crystal transfers.
  1KB0 stays unavailable for matched whole-protein features. Crystals are not
  training rows and never inherit independent status when their sequence was
  already in training. All cases are consumed development, not blind validation.
- Reuse all three GGR and both alpha-lactalbumin prepared structures with matched
  repaired generic-v3 DFT and frozen masked-MACE outputs. These form two biological
  groups, not five independent labels. Their affinity evidence is condition-qualified.
  Aequorin's unmapped site labels, Khoury cross-readout proxies, parvalbumin's
  supporting labels and other unprepared cases cannot enter this direct-site head.
- No new quantum calculations, MACE forwards, structures, protonation, waters,
  optimization, checkpoint training, label changes or PQQ band changes.
- Use the existing lanm_qmmm Python with NumPy/SciPy/Gemmi. CPU-only extraction,
  alignment and tiny convex fits; no GPU/Slurm allocation needed for this step.

## Features and model

Three nested models, fixed before fitting: DFT; DFT+five structural features;
DFT+the same five structural features+two saved MACE response summaries.

The structural features are: number of typed O/N/S contacts within 3.2 A;
fraction of contacts on Asp/Glu sidechain oxygens; excess contacts per donor
residue (a denticity descriptor); mean contact distance; local water-probe
accessible fraction of the selected metal's sphere. Eligible N contacts are
histidine ring N and PQQ N; S contacts are Cys/Met. Backbone amide N is excluded.
These are frozen geometric descriptors, not a changed QM donor inventory.

The accessibility descriptor uses the existing atom radii, a common selected-metal
radius of 1.8 A, probe 1.4 A and 1024 deterministic sphere points in a source-defined
local frame. It measures local steric occlusion, not a solvent-connected cavity
or solvation energy. Both metals use the same source geometry and radius.

MACE summaries: corrected Ca-minus-La native atomic readouts summed over
metal plus <=6 A, and over (6,12] A. Subtract the frozen disconnected-metal
reference at the selected metal only; retain actual full-sum closure. These
are learned descriptors, not physical atomic energies or new free-energy terms.

Fit L2-regularized logistic regression with lambda=1 on the mean weighted
cross-entropy, penalizing standardized feature coefficients but not the intercept.
Fit scaling on training rows only. Give each class equal total weight and each
biological group within a class equal weight; divide a group's weight among its
structures. A zero-variance feature has coefficient fixed to zero. The decision
boundary is logit=0, fixed beforehand; outputs are uncalibrated classification
scores, not affinity or calibrated probabilities. No hyperparameter search,
feature selection, threshold fitting or rescue variants. At most eight feature
coefficients and one intercept. Fit uncertainty is not inferred from structure
confidence or artificial coordinate noise.

## Validation and reporting

Derive sequence from pinned source protein atoms. Join accession duplicates and
sequence-homology connected components at >=50% identical residues relative to
the longer sequence, using fixed global affine alignment (match +2, mismatch -1,
gap open -10, extension -0.5). This threshold is frozen before fitting. Preserve
all transitive groups. Report that this is within-PQQ-family development testing,
not leave-fold-out validation: the existing ledger places all canonical PQQ
proteins in one broad fold family. Never split a cluster to obtain a nicer test.

For each target, hold each complete group out once, with the three model arms
using identical rows/folds. A training fold missing either class is untrainable;
retain its denominator and unavailable predictions. In the direct-site head,
leaving either biological group out is expected to make training impossible.
Do not replace this with structure-level splits or call training accuracy a gain.

Fit final research models on all eligible training groups, and replay the fixed
crystal transfers as explicitly consumed transfer/sequence-replicate checks.
No independent transfer claim for a crystal whose sequence is in training.
Compare PQQ out-of-group correctness, class balance, opposite/inconclusive calls,
and matched changes between the three arms. Show frozen baseline calls alongside
new fitted predictions, with their calibration exposure noted. Report coefficients,
all excluded/untrainable cases, hashes, actual runtime and runnable score commands.

Success requires corrected held-out errors without lost PQQ performance. A larger
margin, a perfect fit on training data, or a single new favorable GGR structure
is not success. The production baseline remains unchanged regardless of outcome.
