# Field-aware candidate: reference scout queued, ML comparison blocked

26 September 2026 checkpoint. **The public model has the intended atomwise
potential coupling, but its exact trained weights and matching engine interface
are not available in the inspected release.** Do not substitute a different
potential. Continue the small embedded-electronic reference scout; obtain the
identified checkpoint/interface before committing to an ML integration.

## What is executable

Real1H4I qm36 core, Ca−3/La−2 physical singlets, paired fixed coordinates and
9,087 matching protein point charges. A deterministic+10degree rotation of
Thr159's hydroxyl H defines environmentB. The source O–H length1.18537Å is
inherited, not an optimized bond; preserve that limitation. Native r2SCAN-3c
analytic gradients supply six endpoints (four embedded plus two isolated).
La uses native46-core def2-ECP, not an imposed effective singlet for Dy.

Job1218751 waits for the other session's final PQQ collector1217591. It requests
64CPU slots/256GiB and uses four16-rank workers; allocation-derived MaxCore
reserves25% for other memory. No GPU. Collector1218752 writes terminal results
and accounting without requiring this chat or a login watcher. No new molecular
calculation has executed, so all new energy/force comparisons remain unavailable.
Measured new molecular cost at this checkpoint:0 allocatedCPU-s and0GPU-s.
Audit/preparation effort is additional, not a molecular timing measurement.

## What the checks establish

Twelve field-algebra/native-parser/allocation tests pass, plus five source and
cap-map preparation checks. Actual archived native output confirms analytic
core gradients and external-charge gradient artifacts can be parsed. Algebra
tests check response through both potential and field, rigid transformations,
charged-gauge conventions and cache identities. They do **not** qualify an
unavailable ML checkpoint or replace new native directional-force checks.

Full hybrid energy and relaxation remain unsupported: standard protein
parameters alone do not supply metal/PQQ cross repulsion/dispersion or a matched
boundary reference. The queued calculation tests the electronic embedding
component only, with all omitted terms explicitly unavailable. Its response is
neither a binding affinity nor evidence of improved classification.

## Upstream implementation findings

Model code `e2b0aeed27c2822790a6994c9a64327e5f131158` accepts atomwise phi and field;
the README is outdated. Default ASE mappings can silently fill absent phi with
zero, ndarray field updates can bypass caching, and actual spin inputs use
multiplicity. Public engine `4af44eb9a9ce91428dad142cbaa8bdd598e7a09f` lacks the
paper's field/potential coupling. Checkpoint-specific field sign, executable
element support and licensing remain unresolved. Closed-shell nonmetal field
training is not metal-response qualification. Exact sources and an unsent author
question are in [the audit](upstream/AUDIT.md).

## Next decision

**Pursue the small reference test; defer ML integration and full hybrid
relaxation.** Inspect actual six-cell results, then qualify native MM/boundary
directional derivatives and region dependence under the frozen tolerances.
Do not expand to LanM or claim environmental accuracy from code tests alone.
No production change, remote push, external email, old vacuum retry or library
campaign occurred. All old results and unrelated work remain intact.

See [CURRENT.md](CURRENT.md) for exact recovery commands/job ownership,
[PLAN.md](PLAN.md) for energy accounting and frozen gates, and
[preparation report](preparation/REPORT.md) for physical mappings and limitations.
