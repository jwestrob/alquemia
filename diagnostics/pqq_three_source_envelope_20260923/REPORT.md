# Practical three-source envelope inputs are prepared

Both existing unlabeled PLM triples pass the new opt-in preparation and dry-run.
The adapter reconstructs the common4.3Å source-fragment pocket and retains exact
paired coordinates and shared chemical state. No molecular calls, new protonation,
scores or classifications were produced. Older source/scoring APIs and production
remain unchanged.

| Protein | Sources | Atoms per context | Status |
|---|---:|---:|---|
| PQQSEQ_07ab500e3df76b30d71c | 3 | 213 | prepared |
| PQQSEQ_83440678cbbd658047c9 | 3 | 201 | prepared |

Each uses the previously declared AF3 La-conditioned samples0/1/2 and original
source selectors. Archived protonation is reused explicitly. The first
lexicographic source defines physical-state identity; it is neither a label nor
a historical calibration-member requirement. All source atoms, proton parents,
caps, bonds, charges, water inventory and protein identity match within each
triple. Both metal mappings replay q0 within1e-12Å. Original donor anchors are
unchanged; complete fragments are selected by the same4.3Å rule on all sources.

The reference is the actual envelope canonical25 artifact
`REFERENCE_PINNED_v1.json`, SHA256
`4a1d3ae6e83b5b6ace408f729c39390d44f673554e338df932fd38fc971f67d9`.
Its calibration gap is7.0713962491modelkcal/mol. Canonical25 and three crystal
calls are retained, but A0A3Ca3 is inconclusive among six consumed probes. This
supports a named experimental preparation, not broad discrimination validation
or a production release. Full structural-transfer qualification remains separate.
Prior threefold and tenfold references are rejected rather than substituted.

## Concrete delivery and remaining work

New `scripts/pqq_three_source_envelope.py` exposes request, prepare, dry-run and
report using the existing physical validation/discovery/preparation kernels.
Actual preparation/report pins and finite scope are in
[ARTIFACTS.json](ARTIFACTS.json). [Commands](COMMANDS.md) provide a real request
and replay example. The task manifest declares12 fresh paired-origin force calls,
12 bounded MACE searches, at most12 cross-MACE calls and72 strict scalarGFN2
calls. It contains24 actual origin scalar inputs; candidate coordinates remain
unavailable until their future searches. No archived molecular values are used
to advertise a fresh result.

The two preparation-function timings are8.918709388s and8.036662072s
(16.955371460s summed). These measure coordinate discovery/reconstruction and
input staging, excluding imports, request validation, final dry-run and historical
source protonation. No GPU or solver allocation occurred. Future molecular cost
is unmeasured; the finite resource profile is recorded without a speed claim.

Six real-artifact tests pass, including exact discovery replay, paired mappings,
whole-source state, unchanged older preparation replay, strict scalar inputs,
missing/wrong reference rejection and explicitly corrupted request/state/task
copies. The first test run exposed only Python tuple versus JSON-list comparison
in the test; JSON-normalized equality fixes it without changing any scientific
artifact. No scientific integration test ran because molecular execution is held.

The remaining step is a small adapter to the existing warm executors, followed
by a separately declared integration trial if root authorizes it. This module
deliberately exposes no molecular execute operation. Unsupported reference or
incomplete triple emits zero executable tasks, without fallback to baseline.
