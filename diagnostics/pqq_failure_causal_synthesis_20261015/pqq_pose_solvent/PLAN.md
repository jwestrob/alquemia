# Conditional solvent continuation for intact PQQ pose

Implements the already declared [PQQ pose plan](../PQQ_POSE_PLAN.md). This binding
and collection adapter changes no electronic kernel, parameters, geometry,
threshold, tie order or production default. Root owns execution and stage review.

## Finite input and gate

Actual preparation4606 admitted all24 new fixed geometries in
`pose_preparation_v1/execution/PREPARATION.json`. MACE4608 scores that same finite
matrix under `pose_scoring_v1/CONFIG.json` SHA256
`3fadca004a90a13afa372ed5c1454ffc2b3987c3e9833ca15f4c28f9a73a204a`.
The solvent config can be bound only after its actual terminal matrix and receipt
exist, with all24 new calls complete and the recorded continuation decision.

Inside Slurm, before any solvent dispatch, replay the existing pose scorer's
preflight and saved collection algebra without model calls. Require all24 new
cells, both complete common9 pools and an absolute mathematical common3→9 paired
contrast change ≥1kcal/mol in either source. Either sign qualifies. Require the
replayed rows and gate to equal the recorded matrix. An unmet gate stops before
SCC; it never changes the threshold or drops a source.

## Exact reused implementation

Reuse `MATRIX.solvent` from the successful packing/reciprocal/midpoint campaigns.
It accepts arbitrary finite cells already; its scientific code is unchanged.
Twenty-four new coordinates produce48 vacuum/water inputs, evaluated by both
published standalone xTB6.7.1 GFN2 and independent tblite0.7.0: at most96 SCC calls.
The existing qualified binary then makes at most48 water/vacuum-charge ALPB
replays. No MACE, ORCA, DFT, optimizer, library rebuild or new force call occurs.

Preserve existing component-v2 parity, real SCC convergence, exact executable and
parameter pins, water-diagonal component replay, identical water/vacuum cavity,
charge rounding bounds and informative vacuum precision gates. Historical failed
absolute-v1 acceptance remains recorded separately, never relabeled a pass.
The candidate solvent term remains native ALPB evaluated at GFN2 **vacuum** charges.
Water SCC/replay supplies existing qualification; it is not substituted into the
candidate energy. All existing four ALPB terms remain included.

Parameter lineage is copied from the actual midpoint config. SCC consumes its
pinned published standalone parameter file; native ALPB uses the qualified
internal factory. The corrected full ORCA La JSON and three-array policy are
checked as provenance only. Default ORCA GFN2 is never called or used as fallback.

## Baseline and comparison

Reuse the twelve exact original cells (two sources × three original geometries ×
two metals). Their identities, MACE receipts, component-v2 acceptance and native
ALPB raw records are revalidated with the existing `reusable` helper. This repeats
no energy calculation. Preserve and reproduce the actual historical3 pool.

Order is original `origin, adaptive_Ca, adaptive_La`, then each original's
`__PQQ_minus1deg, __PQQ_plus1deg`. Reuse the frozen `choose()` policy and existing
full-pool helpers for both mathematical minima and the separately reported
operational strict improvement below−0.1kcal/mol. Origin remains a separate policy.
No bands, labels, calibration, votes or class inference are introduced.

Report all36 endpoint records, per-metal origin-relative works, each finite pose
response, selected geometries, common3→9 differences split into MACE/solvent/ALPB
components, and target-minus-control changes. Any required failed new cell leaves
that source's expanded common9 pool unavailable; its preserved original3 is
reported separately and is never passed off as the new result. No partial minimum.

## Runtime, interruption and evidence

One standard shared allocation of seven CPU slots; seven independent one-thread
SCC endpoints through the existing executor, then the unchanged serial native
component dispatcher. No memory/time/priority/account override. Existing thread
limits and actual Slurm CPU/affinity checks are retained. Root arms terminal and
health wakes as appropriate.

Execution writes into a new directory under a workspace lock and refuses to
replace it. Existing input failures, raw outputs, partial native results and
per-child resource receipts remain. Nested finalization preserves an execution
receipt even if collection fails. No automatic retry or scientific fallback.
The receipt reports measured stage elapsed/allocation time and actual child call
counts and resource files; queue time and previously reused computation are not
claimed as new-stage cost.

Preparation, source copying, hashes and syntax checks are metadata-only. All
scientific validation, gate replay, SCC/native evaluation and final algebra run
inside root-submitted Slurm. Missing MACE outputs mean conditional readiness,
not a completed solvent stage.
