# Technical cap-comparison recovery after4814/4819

Scientific scope, coordinates, state, physical maps, cap definitions, Jacobians,
optimizer and admission tolerances are unchanged. The immutable4814 preparation
remains partial5/10, with five `Rebased caps changed` rejections and zero model
calls. This separately versioned recovery repairs only a comparison between
different numerical representations.

## Actual diagnostic evidence

One-CPU Slurm4819 completed all ten saved-map comparisons in0.546393s, with zero
energy calls. The original and flipped raw cap coordinates, raw map cap inputs,
cap links and physical cap anchors are exactly unchanged in every source. Old
and rebased computed q0 caps are also bitwise identical for every source.

The five rejected rows already have the same residual in the old map when its
computed caps are compared with the raw XYZ. Actual maximum cap residuals (Å):
Ca1:4.440892098500626e−16; Ca3:2.220446049250313e−16;
Ca4:8.881784197001252e−16; La0:8.881784197001252e−16;
La4:1.1102230246251565e−16. The previous exact predicate reproduces all ten
prepared/rejected statuses. There is no cap movement caused by the rebase.

Result: `workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_caps_v1/execution/RESULT.json`,
SHA256 `b59de647f1a0764836d04a1f5e17187978d5bda6b6e5470c781b1bcd33510603`.
Its config, raw per-cap values and execution receipt remain immutable.

## Narrow predicate change

Replace exact equality between rebased computed caps and original raw XYZ with:

1. Exact equality of original and flipped raw cap coordinates.
2. Exact equality of old and rebased raw map cap coordinates.
3. Exact equality of old and rebased computed q0 caps.

Retain both existing full-coordinate q0 replay checks at1e−12Å, unchanged.
No new tolerance, coordinate rounding, cap replacement or hydrogen repair is
introduced. All other physical/covalent/contact/Jacobian and task checks remain
in place. Each actual admission records the two residuals and exact predicates.
The new config pins this note, preserved4814 inputs/results and actual4819 evidence.

Run the same ten-source preparation in a fresh `asn_relaxation_v2` workspace on
one Slurm CPU before using its outputs. The original failed version is not
rewritten. This recovery performs no model calls; root owns submission and the
separate first-source optimizer scout.
