# Fresh source-to-score integration: preflight complete

The two requested execution paths are ready for integration review. Their actual
fresh scientific operation remains untested: no source preparation, energy,
force, search or Slurm submission has run in the fresh execution directories.

The population is all ten declared A0A3F2YLY8 reference structures, one biological
group. The released static path and the opt-in fixed-union adaptive path retain
separate references and independently prepare their sources. Expected maximum
fresh work is 160 GFN2 cells, 60 non-search MACE evaluations and twenty bounded
MACE searches. No DFT, archived molecular result substitution or default change.

A separate, explicitly archived-geometry preflight verified all twenty paired
origin mappings and forty generated rank-one GFN2 input templates. Those files
are outside the actual execution paths and cannot fill a fresh result. Eight
focused tests pass; the actual molecular integration test is explicitly skipped.
The tests also retain all ten rows/strict structural groups when unexecuted,
reject changed precision/call/reuse policies and confirm the private numerical
profile leaves the original optimizer unchanged.

Ready immutable plans and test/source pins: [PREFLIGHT.json](PREFLIGHT.json).
Executable operations: [COMMANDS.md](COMMANDS.md). Scope and timing interpretation:
[PLAN.md](PLAN.md). Plans v1/v2/v3 are retained unexecuted development snapshots; v4
is the intended execution snapshot. The original source/replay interface remains
unchanged. Current implementation deliberately admits this selected compatible
reference group; broader source eligibility has not been established.

Before launch the parent will review this preflight and the full225 verdict.
Each actual job will use 32 CPUs, one H200 and 200000 MiB, sequentially on the
same host, static first. Full job timing includes preparation, model loading and
batching differences; it will not isolate a speedup from accommodation or rank.
Failed preparation, origins, proposals or required matrix cells remain explicit.

## Same-run restart behavior

Completed fresh preparation and scoring manifests can be reused only after exact
plan/source/state/input/reference/implementation checks. Once the released
static runner has started, this adapter performs collection only: that runner
creates its timing and per-source MACE directories exclusively and has no
stage-resume API. No completed or failed molecular task is retried automatically.
A partial preparation without its final manifest remains explicitly unsupported
for automatic restaging. Previous results are preserved; another wrapper call
writes a job-specific collection. This is a technical v4 fix, not archived-energy
reuse or a changed scientific method. See [RESTART_NOTE.md](RESTART_NOTE.md).
