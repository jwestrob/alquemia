# Conditional next direction is locally downhill under both models

19 October 2026. Saved-data job15176 completed the one-direction audit with
zero molecular calls. At the exact 7CCN-La 80-step endpoint, the direction the
unchanged optimizer formula **would construct if continuation were permitted**
has negative nominal derivatives under both whole MACE and admitted whole PBE0,
each with the same fixed-charge ALPB term.

The original optimizer stops at80 accepted steps before constructing another
direction. No direction was actually proposed to an evaluator, no trial
coordinates were formed and no optimization was extended by this audit.

## Saved direction and slopes

| Energy | dE/dα, kcal/mol | dE/ds, kcal/mol/Å |
| --- | ---: | ---: |
| Whole MACE + fixed-charge ALPB | −0.647635 | −5.542160 |
| Whole PBE0 + fixed-charge ALPB | −0.744240 | −6.368858 |

Here `X(α)=X+αd`, with dimensionless α, and `s=α‖d‖₂`. The original eight
history pairs were retained. The fallback was not taken; saved trust remained
0.1Å and the trust scale was1.0. The full direction length is0.116856Å and its
largest atom displacement is0.039884Å. These are saved direction descriptors,
not evaluated finite displacements.

A negative local derivative does not establish finite-step energy lowering,
force accuracy or a stable minimum. The reference did not select the direction.

## Force disagreement remains distributed

The full256 model-minus-reference force-error RMS is3.574975kcal/mol/Å,
with maximum7.660828. The largest residue L2 errors occur at A/9, A/14 and
A/13; all17 peptide residues and the target remain in the result. This does
not identify a unique terminal correction. Atom/residue tables and rankings
are descriptive force summaries, not residue energies or parameter evidence.

Native reference admission is inherited. Only the terminal orbital-gradient
residual passes its printed tolerance; `all_printed_residuals_pass` remains
false. Physical-state and force-accuracy qualifications remain false. Conditional
free termini, the vacuum dianion with post-added solvent, Lu-conditioned7CCN,
the nonconverged endpoint and frozen-f/no-SOC limitations are retained.

## Decision and preservation

Root separately decided to prepare **one 7CCN-La continuation scout**, bounded
to12 additional accepted steps and40 additional distinct evaluations. This
audit itself authorizes no continuation, full campaign or preference claim;
it supplies a local directional observation for that separate decision.

Job15176 used1CPU for2 scheduler seconds: **2 allocated CPU-seconds**, zero GPU
seconds. Molecular calls, model loads, geometry trials and optimizer steps were
all zero. No numerical recomputation was performed while writing this report.

The [lossless result](actual_result_v1/RESULT.json.gz) retains all vectors and
tables; original RESULT SHA256 is
`cc102d11faeb3d7ee7afccd9f772b72a99d2866448fa42e46223f3bbe53bc987`.
[Evidence pins](actual_result_v1/EVIDENCE_PINS.json) bind exact source,
configuration, specification, static review, launch and receipt snapshots.
[Runtime](actual_result_v1/RUNTIME.json) and [scheduler cost](actual_result_v1/COST.json)
are preserved separately. The [result review](actual_result_v1/ROOT_RESULT_REVIEW.json)
is proactive; it does not acknowledge delivery of the queued completion event.
