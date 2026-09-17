# Test a local MACE descriptor without global charge conditioning

Declared2026-09-17 after the whole-chain canonical failure and disconnected
spectator failure. Jacob's active MACE goal and blanket research authorization
apply. This is a new representation and protocol, not a correction to the
finished tests. All their scores and failure decisions remain immutable.

## Question and fixed model change

Does the useful local chemical representation survive removal of the global
charge feature that caused the verified nonlocal sensitivity?

Use the same pinned100M OMOL checkpoint and exact float64 native network.
Replace only the output of `joint_embedding.embedders['total_charge']` with
an exactly zero vector of its original shape/dtype/device, before the existing
joint linear projection and SiLU. Keep the spin embedding and every learned
parameter unchanged. Do not choose a favorable charge category or checkpoint,
fit a coefficient, remove an unfavorable protein, or alter geometry.

Record actual physical total charge, electron count and singlet multiplicity
exactly as before. The preparation is not chemically neutralized. The separate
model metadata must state that charge information is masked in the learned
representation. Native charge-feature weights remain immutable; a versioned
forward adapter supplies the zero feature. Check every masked call and record
its zero output and unchanged spin/parameter state.

**This is an energy-like learned descriptor, not a molecular quantum energy
for the recorded electronic state.** Retain the network's numeric output scale
and one eV-to-kcal conversion, labelled as model units. Do not attach a binding
free energy, aquo reference, physical relaxation force, electronic compatibility
claim or baseline threshold. A successful prediction test would justify further
descriptor evaluation, not establish a physically accurate potential. Suppressing
the problematic feature could also destroy useful information; failure is valid.

Protocol: `mace_omol_intact_charge_feature_ablation_descriptor_v1`.
Score: D_M=T_bound,M−T_detached,M;
R_mask=(D_Ca−D_La)×23.06054783061903 in kcal-equivalent model units.
T is the explicitly modified network output in its eV-equivalent scale.

## Fixed development inputs and gates

Reuse exactly the five intact preparations in intact_report_v1: GGR1GLG,
alpha1F6S,alpha6IP9,MxaF1H4I,XoxF4MAE. All1718recorded peptide bonds pass the
independent integrity audit. Keep source coordinates, whole-chain membership,
protonation, waters, PQQ state, charge, spin and original30A detached construction.
No geometry relaxation or fresh structures. All cases are consumed development.

1. Inspect the adapter on the pinned checkpoint: all representable input charge
   categories must yield identical joint features at fixed singlet, and those
   features must equal the explicit zero-charge-vector plus original-spin
   projection. This is component arithmetic, not molecular energy evidence.
2. Eight core forwards: the four real1H4I/4MAE La/Ca core states from1200807,
   each with native and exact edge/product execution. Require energy and paired
   differences<=0.01 in the model's kcal-equivalent scale. Both executions use
   the same declared feature mask. Failure blocks whole-chain calls.
3. Fourteen alpha1F6S whole-chain calls: primary,repeat,rigid-rotation for each
   metal/bound-detached state, plus the two farther detached checks, exactly
   following the original qualification geometries. Require numerical changes
   <=0.01. No analytic gradients are requested or claimed.
4. Sixteen primary whole-chain calls: four endpoints for each of the other four
   structures. The20primaries including alpha must pass all three original
   relative directions with margins>0.02: XoxF>MxaF, each alpha>GGR. Two alpha
   forms remain one qualified affinity comparison; PQQ functional class remains
   a separate stratum. No fitted decision band.
5. Four GGR calls use the already pinned sodium-perturbed geometries/states
   from1200823. Require its score shift<=0.1, exactly the earlier consistency
   tolerance. Feature masking makes insensitivity an architectural expectation;
   passing it is not an independent demonstration of physical electrostatics.

Total initial inventory:8core+14alpha+16other+4spectator=42new model forwards.
All original descriptor values remain side by side; their caches cannot serve
modified tasks. Numerical or predictive failures remain visible. No additional
mask, reference category, threshold, geometry or checkpoint variant is included.

## Conditional canonical extension

Only if the above numerical and three-direction development gates all pass,
prepare100new calibration forwards(25proteins×4), reusing the eight exact
modified1H4I/4MAE primaries. Use the existing strict v2 preparation and all
28original inventory rows. 1KB0 remains unsupported and in the denominator;
do not recompute its invalid whole-chain state. Preserve all exposure and
grouping flags, and the declared calibration/transfer rules from the canonical
plan: all25valid and min(La)−max(Ca)>0.02 before any candidate-specific bands.
Do not infer broad affinity discrimination from canonical separation.

The actual protein charges and original training-domain flags remain recorded,
even though the model feature is masked. Neither this ablation nor its numerical
checks place whole proteins inside the checkpoint's training distribution.

## Execution

Use the existing runner, immutable snapshots, receipts and scientifically typed
cache. Distinguish the ablated output component from native OMOL molecular
energy, including in acceptance checks and reports. Preserve partial failures.
OneA5000/16CPU/64474MiB, existing exact batching and allocator setting, scheduler
QOS only. Recent PQQ endpoints take roughly36seconds including loading;
smaller proteins are quicker. Record all preparation, validation, actual GPU/
CPU/memory costs and failed attempts. No project time/CPU budget, newDFT,
solver, training, gradient, trajectory, environment install or production change.
