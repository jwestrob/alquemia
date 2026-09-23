# Fresh source-to-score candidate corrects the two difficult A0A3 calls

The opt-in union/adaptive candidate completed fresh preparation and scoring of
all ten A0A3F2YLY8 structures and classified **10/10 correctly**, versus the
released static MACE+GFN2 path's **8 correct / 1 wrong / 1 inconclusive**. It repairs
Ca-conditioned samples1 and3. Both paths have10/10 coverage, with no archived
energy/force substitution. This demonstrates usable execution of the observed
structural-robustness gain on one consumed reference protein, not ten independent
biological validations. Production and the explicit DFT baseline remain unchanged.

## Actual calls and cost

Sequential jobs1211081 then1211082 ran on node-224-2t-8gpu-1, each32CPUs, oneH200
and200000MiB. All160 native GFN2 cells and332 MACE evaluations succeeded. Candidate
work comprises20 fresh origin evaluations,272 search evaluations and20 cross-metal
evaluations; the static path used20 scalar MACE evaluations. All20 bounded searches
were admitted; five finished on a physical/angular boundary. These are finite
candidates selected by their own energy, not unconstrained stationary minima.

|Measured item|Released static|Union/adaptive candidate|
|---|---:|---:|
|Whole allocation wall|451 s|358 s|
|Allocated CPU-seconds|14,432|11,456|
|Requested GPU-seconds|451|358|
|Reported distinct-step CPU-seconds|3,881.556|1,205.387|
|Fresh source preparation|111.257 s|114.933 s|
|GFN2 cells / normal completions|40 /40|120 /120|
|Actual MACE evaluations / complete|20 /20|312 /312|

Total25,888 allocated CPU-seconds and809 requested GPU-seconds, including failures,
staging and in-job collection. No newDFT or molecular retries. Local submission,
tests, repaired collection and reporting are additional unmetered CPU work.
The candidate allocation was93s shorter in this actual two-path comparison.
This includes batching, model initialization and solver concurrency differences;
it does not isolate a speedup from accommodation or ranks, or outperform an
unimplemented equally optimized static scorer.

Released stage wall:105.621s MACE process stages and223.876s GFN2 stages; its
prepared-input score path is332.487s. The static implementation initializes the
model separately for20 endpoints (25.102s initialization/input setup total).
Candidate uses three persistent stage workers (4.015s measured model initialization
in total),83.708s bounded searches and94.150s across the two native solver executor
stages. Stage timings overlap/nest with receipts and must not be added to allocated
costs. Raw job/step accounting, preparation, solver, model and search receipts are
retained in COSTS_v2.json; zero reported RSS on some steps is missing accounting,
not evidence of zero memory use. [Timing details](TIMING_NOTE.md).

## Decisions, reproduction and limits

Larger R is more La-like; R remains Ca−La in model kcal/mol. Each path uses its own
unchanged frozen reference. Candidate bands are Ca_max−405463.71230032056 and
La_min−405456.3869114709, from the complete25 canonical reference. No threshold
was fitted or changed on these ten structures.

|Source conditioning/sample|Released R|Released call|Candidate R|Candidate call|
|---|---:|---|---:|---|
|Ca0|-405457.618084|La-supported|-405454.078006|La-supported|
|Ca1|-405464.654353|Ca-supported|-405450.717111|La-supported|
|Ca2|-405443.765110|La-supported|-405441.319396|La-supported|
|Ca3|-405461.993448|inconclusive|-405455.470283|La-supported|
|Ca4|-405440.471276|La-supported|-405433.438952|La-supported|
|La0|-405458.912508|La-supported|-405452.986231|La-supported|
|La1 (canonical)|-405459.111382|La-supported|-405444.232685|La-supported|
|La2|-405439.987137|La-supported|-405434.258520|La-supported|
|La3|-405456.894855|La-supported|-405448.445819|La-supported|
|La4|-405437.351371|La-supported|-405427.949089|La-supported|

Both paths already give correct strict La4, Ca5, balanced and all four La-triple
aggregate calls. Thus this fixture improves individual-source calls, not those
already-correct aggregate decisions. Noncanonical Ca-fold spread falls24.1831→
22.0313kcal; La-fold spread grows21.5611→25.0371. Do not claim universally narrower
structural spread.

Against the already computed corresponding artifacts, maximum released R change
is5.12e−9kcal and candidate R change1.3912e−5kcal. All ten candidate metal-row
selections remain identical. Every compared native MACE component is identical;
maximum GFN2 component difference is1.39113e−5kcal. Fresh static source coordinates
differ by at most1.78e−15Å (last-bit serialization); candidate origins match exactly,
and final candidate coordinates differ by at most4.67e−14Å. These comparisons use
old outputs only after execution and do not supply fresh scores. All final state,
source mapping and paired-coordinate checks passed.

This explicitly selected compatible ten-reference-source route is implemented and
executed. General unknown sources, different ensemble membership, fewer folds and
three-fold-specific context preparation remain separate qualifications. The
canonical identity is archival, not an invented anchor for an unknown protein.

## Visible technical failure and tests

Released job1211081 has Slurm FAILED/1:0 because its final report adapter matched
`result_JOB__CASE__GFN2.json` component files as well as the one top-level scanner
result. All released preparation and molecular calculations had already succeeded;
the original scanner result records10/10. The collector now selects only the
actual top-level `result_JOBID.json`. Collection-only recovery retained the original
error/logs and immutable v4 snapshot; no endpoint was repeated. Candidate1211082
completed normally. COLLECTION_RECOVERY_v2.json pins the old error, exact actual
scanner result, fixed collector and final collection.

Eleven actual-fixture tests pass, zero skipped in TESTS_final_v2.txt. These include
full fresh integration, common-pool algebra/all60 candidate cells, all160 real
solver completions, input/state/restart checks, the real component-filename
regression, and deliberately corrupted-policy/missing-result handling. Earlier
preflight scientifically unrun statuses remain preserved. No synthetic successful
scientific output was used.

Recommendation: keep the candidate opt-in and usable for this qualified source
domain; retain the released default. Broader scientific/source qualification
continues separately. [Commands](COMMANDS.md), [artifact pins](ARTIFACTS.json).
