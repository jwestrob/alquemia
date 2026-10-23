# Journal initialization and terminal reporting recovery

The original local_comparison_v1 is preserved. Its first request failed because
ledger initialization omitted limits before the journal validated request_start.
The journal remained valid with initialize/stop only. Live request mutation then
made its final snapshot disagree; the stale started RESULT and both failures are
retained. Root failure review f26481cf records zero native and solver attempts,
with controller cost unavailable.

This version adds a detached copy of the existing four-vector/112-program limits
before Journal construction. No journal reconciliation is permitted. Normal
snapshot errors still stop execution. Terminal reporting retains the primary
exception, separate stop/snapshot errors and uncertified live counters plus raw
persistence pins. A snapshot error after apparent success changes the result to
failure and removes any stale evaluator snapshot.

MATH.py, SOLVE.py and SPEC.md are byte-identical to v1. Inputs, four probes, old24
actual spans, admitted zero bridge, binary, gates, two solver settings and resource
envelope are unchanged. The cache/output path is fresh. This is a replacement for
the zero-call failed attempt, not an added scientific campaign.

The bounded manufactured control runs the actual public evaluate wrapper and
journal through one durable attempt_intent, then deliberately stops in a fake
helper before any subprocess. It uses46 synthetic slots and one test worker;
there is no successful molecular vector or numerical qualification claim. A
separate injected request_start rejection checks unchanged journal mismatch and
terminal persistence, including apparent-success failure. The first harness-only
event-envelope schema mistake is retained under controls_v1; corrected controls_v2
passes14 checks, zero native/solver calls. Scientific preparation is inherited by
exact metadata/source-delta checks, without replaying the old46 cells or old24
probes. Root release remains absent until exact source review and binding.
