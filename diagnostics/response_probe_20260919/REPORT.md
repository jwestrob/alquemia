# Archived structural response: radial hybrid improves consumed direct ordering

Completed 2026-09-19. **The PQQ-trained DFT+radial model fixes the ordering of the
consumed alpha/GGR comparison while preserving PQQ fidelity.** Using saved PQQ
weights without refitting, it orders La-associated alpha above Ca-associated GGR
in 4/4 structural/partition comparisons; DFT-only achieves 0/4. Those four rows
represent one biological comparison. PQQ remains 25/25 plus 3/3 consumed crystal
transfers. This is a promising exploratory transfer result, not independent
validation or an absolute affinity classifier. The baseline/default is unchanged.

**Subsequent result:** the [continuation](CONTINUATION_REPORT.md) passes8GY2 but
reduces direct structural robustness to2/6 across threeGGRstructures. The4/4
result below remains the original1GLG/partition observation.

## Fixed experiment and complete denominators

[Plan](PLAN.md) preceded extraction/fitting. Reused 64 native MACE-POLAR-medium
vacuum force endpoints and eight native r2SCAN-3c/CPCM analytic gradients.
32/32 preparations complete: 25 canonical PQQ controls, three consumed crystal
transfers, and four direct preparations representing two biological groups.
No missing cases, molecular calls, GPU allocations, new folds, PLM rescores,
geometry/protonation changes, default changes, or fitted mechanical corrections.

Two rotation-invariant observables use physical source-mapped donors within
3.2 A: mean radial Delta g and tangential fraction, with Delta g = g_Ca - g_La.
Synthetic cap derivatives map to both true source atoms. PQQ N/O typing retains
the complete cofactor. The common cutoff includes legitimate contacts beyond
the old truncated 3.1 A donor list; it never changes the calculated core.
All donor vectors, cap mappings and descriptive cofactor force/torque values are
retained in the feature artifact. [All 32 cases](ALL_CASES.tsv) lists observables
and every class call. Blank direct class/absent DFT-gradient fields are unavailable,
not zero. PQQ classes and direct affinity directions remain separate targets.

## PQQ grouped classification

Same four homology folds (18,5,1,1), class/group-balanced weights, training-only
standardization and fixed ridge lambda=1 as the earlier classifier. All five
predeclared arms are reported; no model or threshold was selected by outcome.

| Fixed arm | Correct held out / 25 | Consumed crystals / 3 |
|---|---:|---:|
| DFT_R | 25 | 3 |
| Radial only | 20 | 1 |
| DFT_R + radial | 25 | 3 |
| Tangential only | 16 | 0 |
| DFT_R + tangential | 24 | 2 |

Earlier DFT+structure also achieved 25/25; adding masked-MACE energy readouts
achieved 24/25. The present native-force features are different observables and
a different MACE Hamiltonian from those masked OMOL energy readouts.

DFT+radial retains fidelity but reduces the minimum signed held-out logit margin
from 0.120387 to 0.035950. For Ca-associated 1H4I, its signed margin falls from
0.355324 to 0.066250; for La-associated 4MAE it rises from 0.489830 to 0.990633.
These uncalibrated logits do not establish confidence or a better overall spread.
DFT+tangential misclassifies q4w6g0-pqq-la_model and the consumed 1KB0 transfer.
Radial alone misclassifies 1H4I and 1KB0; tangential alone misses all three crystals.
All per-case predictions, coefficients and fold memberships are retained.

1H4I and 4MAE share a known homology group with training cases. 1KB0's group was
not computed by the earlier grouped classifier and remains unknown. All three
have already been inspected; none constitutes new blind validation. Composition
already separates the canonical classes, limiting incremental-information claims.

## Direct preparations: interesting, context-dependent radial ordering

| Preparation | Donors | MACE mean radial | DFT mean radial | MACE/DFT donor-vector cosine |
|---|---:|---:|---:|---:|
| alpha 1F6S | 7 | -2.990282 | 7.220856 | 0.245735 |
| alpha 6IP9 | 9 | -1.413574 | 7.031484 | 0.900102 |
| GGR extended | 7 | 6.673648 | 15.238832 | 0.793314 |
| GGR connected | 7 | 7.552458 | 15.076114 | 0.860555 |

Radial units: kcal/mol/A. Both Hamiltonians put both GGR representations above
both alpha structures. GGR extended/connected are partitions of one 1GLG
structure, not distinct biological observations. Alpha structures retain the
older frozen proton/water preparations; newer water experiments are separate.

**This is not independent affinity validation or physical force agreement.**
MACE vacuum and DFT CPCM have different Hamiltonians and opposite absolute radial
signs for alpha. Individual radial signs agree for only 3/7, 5/9, 5/7 and 5/7
donors in the table order. In PQQ, the fitted radial coefficient is negative:
lower load is more La-like, consistent with alpha being below GGR. A common
ordering does not establish a transferable zero, PQQ threshold or absolute sign.
Tangential fractions do not consistently order the direct groups.

GGR connected minus extended radial changes are +0.878810 (MACE) and -0.162718
(DFT) kcal/mol/A; tangential changes are -0.025674 and +0.000184. These are
representation sensitivities, not relaxation energies or qualified affinities.

## Saved PQQ weights improve direct ordering without refitting

[Transfer extension](TRANSFER_PLAN.md) was frozen before these outputs. Use the
saved all-25 PQQ coefficients/scales for DFT-only and DFT+radial, with no new fit.
Report alpha-minus-GGR logit differences; intercept and training mean cancel.
Positive differences give the experimentally supported ordering. No absolute
class, affinity probability or inherited PQQ band is computed across targets.

| Consumed comparison | DFT-only difference | DFT+radial difference |
|---|---:|---:|
| 1F6S - GGR extended | -0.412424 | 0.270990 |
| 1F6S - GGR connected | -0.209415 | 0.514459 |
| 6IP9 - GGR extended | -0.478508 | 0.105432 |
| 6IP9 - GGR connected | -0.275498 | 0.348901 |

This is incremental discriminatory information in one consumed biological
comparison, learned from a separate PQQ target; the direct labels were not fit.
The fixed radial contribution overcomes the unfavorable DFT ordering in each
preparation. Different protocols/targets and prior inspection limit this result;
it needs further controls, not a claim of general affinity prediction.

Reporting correction: initial v1 direct labels were manually reversed. The
source evidence ledger explicitly has alpha=La and GGR=Ca; v2 now looks these up.
All numerical features and PQQ fitted weights are unchanged. V1 is preserved
with `V1_METADATA_SUPERSEDED.md`; it must not supply direct scientific labels.

## Checks and measured cost

Eleven real-fixture tests pass with zero skips in 2.081 s: paired charges/coordinates;
force-to-gradient sign/units; arbitrary rigid transformations; cap Jacobians
checked against geometric derivatives; corrupted-real-input failures; exact
native DFT method/charge identity; complete PQQ/donor inventory; grouped-fit
replay and no held-group leakage; authoritative label lookup; unchanged numerical
features/models after metadata repair; and saved-weight transfer algebra.
These are parser/algebra/archive checks;
no new molecular scientific integration test ran.

All 32 mappings pass the frozen closure tolerance. Maximum cap force/torque
closure residual is 5.2541e-5 (torque, kcal/mol); maximum rigid-observable algebra
error is 5.3291e-15. Standardized DFT_R reproduces the earlier DFT_S logits within
2.4771e-10, demonstrating the common aquo-offset cancellation.

Extraction: 18.966292 wall s / 18.194511 CPU s, peak RSS 823956 KiB.
Fitting: 0.311326 wall s / 0.149978 CPU s, peak RSS 824196 KiB.
These are measured operation-body times, excluding Python imports. The 64 reused
MACE endpoints originally recorded 70.309680 s of summed model evaluation;
that historical figure excludes preparation/startup and is not new cost or a
production runtime estimate. No new GPU time or scheduler allocations.
Saved-weight transfer adds 0.036166 wall s / 0.034277 CPU s, zero new fits.

## Artifacts and recommendation

- Features: `workspaces/response_probe_20260919/features_v2.json`, SHA256
  `b1aa67065243f226d83be9192863aa94f24c3806da34d6840fa8ac6542574ac7`.
- Evaluation: `workspaces/response_probe_20260919/evaluation_v2.json`, SHA256
  `e2b742ed2cef63ea7e7eb77d3d48413a85ef6059e66f5eb8f391dffd6665d4c8`.
- Saved-weight transfer: `workspaces/response_probe_20260919/transfer_v1.json`.
- Source manifest/force/coordinate hashes and model/software settings are in
  features; source code, tests and predeclared plan are preserved under
  `workspaces/response_probe_20260919/implementation_v1` and `implementation_v2`.
- [Runnable commands](COMMANDS.md); protocol
  `mace_polar_archived_donor_response_probe_v1` is opt-in research only.

**Recommendation: pursue the radial hybrid as an opt-in research candidate,
keeping the baseline/default intact.** It supplies a useful new ordering signal
without fitting direct labels. Reject the tested tangential addition; this
result does not justify a larger fitted classifier or a mechanical correction.
