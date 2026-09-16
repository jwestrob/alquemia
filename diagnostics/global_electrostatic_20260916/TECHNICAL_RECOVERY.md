# Technical collection and provenance fixes — 2026-09-16

These are software integrity fixes. They change no physical state, numerical
settings, acceptance threshold, energy expression, or approved task count.
No quantum or surface task is rerun for these changes.

## Malformed input handling

- The TABI parser now rejects an empty/truncated CSV header as `InvalidArtifact`
  instead of allowing `StopIteration` to escape. A single malformed native
  output is therefore recorded as a failed task without aborting collection of
  the other tasks.
- The assessor requires exactly 25 collection rows before indexing them by
  task ID. Duplicate rows cannot silently replace a failed row. Its existing
  unique-ID/schedule check also rejects a duplicate that replaces another task
  while preserving the row count.

Tests use archived real TABI outputs and explicitly corrupted copies of the
actual frozen campaign collection. They do not manufacture successful solver
outputs or claim new scientific integration.

## Adapter identity during the running jobs

Jobs **1199956** and **1199959** imported adapter source with SHA256
`50ed8e24668c3777c66f34c888208548332e357769171f248e3bc3dcec0f2454`.
That exact source is preserved in each task manifest's
`implementation.affordable_tabi.py` record. The pinned NanoShaper bridge also
uses the preserved adapter copy.

The old executor records `record(__file__)` when the subprocess finishes. If
the live source was edited after import, its receipt's `executor` field can
therefore identify the filesystem bytes present **at finish**, rather than
the Python code already loaded into that process. For these two jobs, use the
preserved manifest implementation above as the executing adapter identity;
retain the finish-time record as historical filesystem provenance. Existing
receipts are not rewritten, and jobs are not restarted for this parser fix.
The later same-input hardware recovery of1199959 is documented separately in
PERFORMANCE_OBSERVATION.md and the experiment report.

Future imports capture `LOADED_IMPLEMENTATION` immediately and use that record
for the executor field. Future receipts explicitly mark the identity semantics
as `source_identity_captured_at_module_load`.

Collection verifies immutable task implementation snapshots and execution
artifacts. It deliberately does not require the current parser source hash to
equal the executing adapter source hash. A new collector can therefore parse
the original saved outputs, and its own parser identity is recorded separately.
