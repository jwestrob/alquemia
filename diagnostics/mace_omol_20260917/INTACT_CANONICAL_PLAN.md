# Whole-chain OMOL canonical calibration and transfer

Declared 2026-09-17 before preparing or scoring this extension. The five-case
development comparison passed its three relative criteria. It has no absolute
classification bands. This experiment asks whether the same intact coordination
descriptor supports the established PQQ reference panel, across whole-protein
charge and sequence variation.

## Frozen cases and preparation

Use all 28 rows of `workspaces/mace_canonical_20260916/audit_v2/inventory.json`
(SHA256 f51a8641fa50a34087fca9ba76d0e5ca960ac85c1ea5da65a3fa19bed9d5cb21):
25 designated calibration structures, then the consumed 1H4I, 4MAE and 1KB0
crystal transfers. Keep their existing sequence groups, labels, structures and
exposure history. These functional-class labels are not direct affinity labels.

Reuse the exact successful whole-chain 1H4I/4MAE preparations and eight computed
bound/detached endpoints. Prepare the other 26 sources from their archived
protonated coordinates. Select chain A, retain the existing oxidized PQQ3-minus
core cofactor coordinates and zero canonical site waters, and apply the existing
ff19SB radial protein-H repair. Do not rerun protonation or move source heavy atoms.

Extend the existing terminal OXT construction by a general connectivity rule:
only the final protein residue, with carbonyl C bonded only to its same-residue
CA and O, can receive a missing OXT using the existing reflected geometry at
1.25 Angstrom. Internal gaps, unmatched residue templates, missing nearby atoms,
or unexpected nonstandard protein chemistry remain unsupported. Record every
addition and exclusion. This covers the independently identified missing 1KB0
terminal OXT before its score is known. The other 25 sources already contain OXT.
New preparation ID: `omol_canonical_chain_A_ff19sb_H_terminal_OXT_v1`.

Use each actual integral formal charge and multiplicity one, with Ca charge one
less than La. Four calibration proteins are outside OMol25's reported -10..+10
training charge range; retain them explicitly in the experiment, without charge
clipping, neutralization, protonation changes or outcome-based exclusion. Their
model inputs remain within its categorical representable range. Report this
additional extrapolation per case; all whole proteins also exceed training sizes.
Numerical success does not remove either limitation. Preserve unsupported cases
in denominators, with unavailable scores.

## Model, tasks and decisions

Retain the pinned MACE-OMOL-0 100M checkpoint, float64 native energy, and the
qualified exact execution adapters. No solvent, aquo offset, changed functional,
new reference geometry, relaxation, learned classifier or fitted combination.
Use the existing bound/detached construction, moving only the selected metal to
max(other x)+30 Angstrom at unchanged y/z and requiring zero metal neighbor edges.

For each new case evaluate bound/detached x La/Ca: at most 104 new energy calls,
plus eight explicitly verified earlier endpoint reuses. All new inference waits
for both actual stages of EXACT_PRODUCT_PLAN.md and native full equivalence.
Use 1024-edge and 1024-atom product batches. Record the adapter separately from
the unchanged energy expression and validate its exact source hashes.

R_coord = (E_bound,Ca - E_detached,Ca) - (E_bound,La - E_detached,La), converted
from eV to kcal/mol once. Larger is La-like. Canonical calibration requires all
25 scores numerically valid and min(La)-max(Ca)>0.02 kcal/mol. Only then publish
this candidate's research bands: Ca <= max(Ca), La >= min(La), otherwise
inconclusive. Apply them unchanged to all three transfers. Any missing case,
failed separation or incorrect/inconclusive transfer fails the overall canonical
operational gate. Also report charge-range strata without fitting separate bands.

Do not apply PQQ bands to the alpha/GGR affinity comparison. Preserve the three
original whole-chain relative criteria and their completed results separately.
All evidence is consumed development/retrospective evidence. Calibration accuracy
does not establish generalization, and the canonical panel cannot demonstrate
added information beyond motif/composition by itself.

## Execution and outputs

Reuse the existing runner, snapshot, cache and receipt machinery. Candidate
protocol: `mace_omol_intact_canonical_class_v1`. Keep the production baseline
unchanged. Save preparation status, mappings, state checks, all unrounded terms,
per-case domain flags, baseline comparisons, calibration/transfer decisions and
actual cost. No uncomputed correction is zero-filled.

One A5000, 16 CPUs, 64474 MiB host RAM, expandable-segments allocator, existing
scheduler QOS. Measured 9k-atom energy-only calls take about 26 seconds each;
104 calls would be roughly 45 minutes of inference at that size, before model
loading, preparation, validation, collection and failures. This is an estimate,
not a budget. Larger inputs need the qualified product memory adapter. No DFT,
solver, training, trajectory, optimization or force calls are planned here.

Sources for domain limits: [OMol25 paper](https://arxiv.org/html/2505.08762v1)
and [MACE charge embedding explanation](https://github.com/ACEsuit/mace/discussions/1206).
