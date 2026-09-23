# Same-run restart policy — technical v4 correction

The parent review identified that v3 would validate an existing fresh preparation
but then restage it unconditionally. No molecular job had run. V3 remains intact;
v4 corrects this orchestration without changing sources, chemistry, optimizers,
solver settings, references or declared call counts.

Every completed fresh preparation now pins its execution plan. Reuse requires
that exact plan, request, all ten source identities, configuration and fresh-source
status. Candidate union preparation also checks its original fixed membership
and state signatures. Static scoring manifests must match the same preparation,
endpoints, recipes, software, implementation hashes and frozen released bands.
Archived mapping-only fixtures are explicitly rejected as same-run preparation.

The released `compact_solvation_scanner.execute` uses exclusive creation for
`execution_started.json`, its timing directory and each MACE result directory.
It cannot automatically resume after it starts. The adapter therefore calls that
runner only if no execution marker/output exists. Any already-started run is
collection-only; the original failures and completed work remain unchanged.
Unfinished preparation/staging without its final manifest is not automatically
regenerated. Manual investigation would be needed before any further execution.
This limitation is explicit rather than advertising general failed-job recovery.

A static lock protects the prepare/execute path. Wrapper collection still runs
following execution failure. If RESULT.json exists, a new RESULT_JOBID.json is
used, preserving the old result. This is within-run artifact reuse, never an
archive value filling a missing fresh molecular cell.

Actual-fixture checks: the real released source-to-score result is recognized as
already started/collection-only; the new archived-coordinate input-only fixture
is recognized as unstarted; an archived preparation cannot pass the fresh-plan
restart guard. No fake energies, execution receipts or successful molecular
results were generated for these tests. Eight tests pass; actual fresh molecular
integration remains explicitly unrun.
