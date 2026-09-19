# Static second-shell quantum effect across GGR replicas

Authorized by Jacob's parallel accuracy-track instruction and root's explicit
2026-09-19 continuation: use the unchanged static complete-second-shell rule on
both archived 2FW0/2FVY structures; report the full two-alpha × three-GGR matrix.

The completed coordination optimization pilot did not improve discrimination.
This continuation does not rescue or alter it. It tests the separate static
second-shell result: alpha−GGR gaps improved by about 3.9 kcal/mol for 1GLG, while
the water-prepared baseline's weakest archived GGR-replica margin is 0.655 kcal/mol.
Does adding the actual local quantum environment improve that weakness consistently?

## Frozen scope before execution

- Cases: exactly GGR 2FW0 and 2FVY, both Ca/La endpoints. Four new expanded-context
  native r2SCAN-3c/CPCM(Water)/DefGrid3 singlepoints. Reuse all four exact archived
  original-core endpoints; no original-core recomputation.
- Eight native unmasked OMOL100M evaluations: original and expanded × both metals
  for each structure. Original coordinates have not been normalized or substituted
  with the prior response/whole-protein preparations. Exact existing matches may
  be reused only with matching input/charge/model provenance.
- Reuse `second_shell_context.parent_state` and `expansion` without changing their
  rule: complete standard protein polar neighbors within 3.5 Å of the original
  direct-donor groups, actual source connectivity/overlap completion, explicit
  retained/cut source mapping, fixed waters and all retained core coordinates.
  The original donor cutoff remains 3.3 Å. Unsupported nearby species remain
  explicit errors; excluded outer waters preserve the original water inventory.
- Same source structures, protonation, original caps and formal-charge policy.
  Overlap completion may replace a cap with actual source atoms, exactly as in
  the original static pilot. Added charge/composition/cavity changes are reported.
  No geometry search, changed reference, new bands, fitted threshold or label change.

## Energy and evaluation

Each endpoint is the complete native continuum-solvated energy of its explicit
core/context. `R=E_Ca−E_La`; environment expansion changes it by
`ΔR=R_expanded−R_original`. Existing alpha 1F6S/6IP9 and GGR 1GLG static core/context
energies are reused. Report all six alpha−GGR differences, before and after, and
GGR representation spread. Retain PQQ's previously measured static 30.57→28.18
kcal/mol gap as a separate development result. Do not attribute a context effect
uniquely to hydrogen bonding, electronic polarization, added charge or cavity.
Native OMOL vacuum is reported separately, not added to CPCM DFT.

All structures are consumed development; three GGR structures are one protein
observation, two alpha structures another. This tests structural robustness and
the direction on known evidence, not independent accuracy or affinity magnitude.

Use one existing GPU allocation for eight small calls and one 64-CPU allocation
for four 16-rank expanded SPs, with 256 GiB host RAM. Prior 1GLG 125-atom context
receipts establish the same size class; exact dimensions and measured costs are
retained in the new manifest/report. No arbitrary aggregate compute cap; finite
declared tasks and ordinary scheduler policies apply. Production stays unchanged.
