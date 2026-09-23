# Reusable experimental execution interface, version2

Parent assignment on2026-09-23: remove the hard-coded two-PLM-ID restriction for
new explicit, complete three-source requests, retaining exactly the tested
four-angle/4.3Å envelope/strict-native method. Parent confirmed: “Keep molecular
scientific model/reference identical while versioning execution interface only.
All new execution remains unrun; historical result remains v1.”

Change only the owned execution adapter and focused tests/docs. New finite plans
derive groups, source IDs and counts from explicit validated preparations. Every
group must contain exactly three complete declared La-conditioned sources with
matching source protein/atom/proton/cap/bond/charge/water state. Unsafe identifiers,
mixed source state, repeated sources/groups, unsupported reference or incomplete
triples remain errors. No hidden baseline fallback or reduced-member median.

Preserve v1 plan validation, immutable executed snapshots and the completed
six-source result. New plans use execution-interface v2; physical method,
checkpoint/software, strict native recipe, optimizer, scalar qualification and
actual envelope reference remain unchanged. No new calibration, label, chemistry
or default promotion. Existing v1 result remains historical; validation/recollection
is permitted without changing its protocol or claiming a fresh v2 calculation.

Only preflight the already prepared real07ab triple alone,8344 alone and both
together. Per supported group, declare6originMACE,6bounded searches,≤6crossMACE,
36strictGFN2 cells. Two-group plan doubles those counts. All plans are dry-runs;
no allocation, molecular call, new source preparation or scientific comparison.
Tests use real source hashes, an explicitly named metadata alias, genuinely
different real protein sequences for rejection, and corrupted copies of actual
state/missing-member records. Never fabricate successful scientific output.

The candidate remains experimental. Its actual canonical25+3crystal calibration
retains fidelity, but A0A3Ca3 is an inconclusive probe and full100-triple transfer
is a separate qualification task. This API change supplies no new accuracy
evidence. Execution commands are documented for deliberate finite requests; none
will be submitted in this task.
