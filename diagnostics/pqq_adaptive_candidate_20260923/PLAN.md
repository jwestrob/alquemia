# Opt-in minimal adaptive PQQ source command

Approved implementation scope, 2026-09-23: parent requested a usable candidate
entrypoint for explicitly selected compatible PQQ sources, and approved the
narrow interface after inspection. This delivery runs no new molecular calls.
Jacob's standing overnight authorization remains recorded in
`diagnostics/overnight_discrimination_20260923/CURRENT.md`.

Keep the released source preparation: explicit crystal-chain or Protenix PQQ
normalization; fixed homologous roles; dry PQQ3-; pH7; complete polar context;
paired heavy/H coordinates and La charge = Ca charge + 1, both singlets. No
new normalization, missing-atom repair, microstate or context choice. Source
preparation or fewer than four independent angular modes can be unsupported.

Reuse actual native OMOL origin forces to select the same four physical modes,
scaled SLSQP200, angle bounds +/-0.8 rad and final heavy displacement <=0.8 A.
Both endpoints start at the common source. Cross-score origin/adaptive_Ca/
adaptive_La with native OMOL + native GFN2 ALPB - vacuum. Require the complete
matrix. Native GFN2 MaxIter500 changes the allowed iterations only; preserve
qualification and the known native SCF limitations. No restart/seed selection.

Frozen canonical25 minimal adaptive reference is required, never fitted here.
Report original released score/bands and accommodated candidate score/bands
separately, both mathematical and operational (0.1 kcal/mol origin retention).
Unknown source-domain evidence stays explicit; compatibility is not validation.

Fresh execution, when explicitly invoked later: one 32-CPU, 200000-MiB GPU
allocation; up to four independent eight-rank ORCA tasks concurrently, between
GPU stages. Per supported site: two origin MACE energy+force calls, two bounded
MACE searches, two cross-MACE calls, four origin and eight candidate GFN2 calls.
No DFT. Preparation/dry-run do not execute these calls.

Archive replay is a separate request type limited to validated completed
origin-recovery pools. It cannot stand in for fresh source-to-score integration.
Tests use the two executed recovery sources and released 1H4I request/preparation,
including explicitly corrupted copies for input failures. No fake energies.
