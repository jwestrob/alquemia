# Affordable discriminator development

## Authorization and ownership

Jacob supplied the full “Alquemia: preserve the working discriminator and test
affordable improvements” implementation assignment in this conversation on
2026-09-15, then said “Sorry, upped your reasoning effort. Continue.”
That direct assignment authorizes the live audit, versioned peptide-amide
repair, regression/comparison machinery, bounded backend implementation and
gradient interface. It preserves the production default and frozen records.
No production rescore, push, queue changes or watcher changes are authorized.

Owned additions: `scripts/affordable_*`, `tests/test_affordable_*`, this
diagnostic directory, and `workspaces/affordable_challenger_20260915/`.
The existing dirty checkout and untracked implementations belong to prior work;
they are dependencies, not changes claimed by this session. The checkout is
kept in place to use its actual pinned dependencies; only explicitly named new
files will be committed. No baseline source is edited.

## Agreed implementation scope

- Reproduce score extraction and published bands from the archived 25-member
  fixed-core PQQ calibration and two crystal transfers. No DFT reruns.
- Prepare a new graph-based peptide-amide protocol for **all six** previously
  prepared non-PQQ sites, preserving their donor selectors, waters,
  protonation and source heavy coordinates. These are preparation integrity
  checks, not new affinity calculations. Test additional peptide units in
  the same real structures for proline, overlapping selections and breaks.
- Retain baseline, repaired baseline and optional environmental records
  separately, with unavailable numerical corrections represented as null.
- Reuse analytic-gradient parsing, provide physical cap mappings and a gated
  response interface. No validated curvature is assumed; no mechanical score,
  entropy, displacement DFT or trajectory is authorized by this record.

## Approved environmental pilot

The asynchronous question proposes APBS/MBIS frozen-distribution transfer on
consumed 1H4I fixed-core and 3.3/3.6 A geometries: six endpoint evaluations,
two reserved retries, eight total. Dielectrics 1/78.54, zero salt, 298.15 K;
fixed radii and common whole-chain physical boundary. Identity tolerance
0.01 kcal/mol; grid/box/rotation tolerance 0.5 kcal/mol; partition-score
tolerance 2 kcal/mol. No pilot output has been inspected. Exact manifests,
backend feasibility, cost evidence and aggregate solver budget must precede
execution.

Jacob subsequently approved the explicitly restated pilot: “go get em tiger.
you have blanket permissions to launch what you like. this work is high
priority.” This authorizes execution of the above scope and routine technical
fixes; the initial eight-endpoint limit remains in force. New biological
comparisons or changed scientific models are not silently added.
