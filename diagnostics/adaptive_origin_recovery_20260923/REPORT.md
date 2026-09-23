# Direct-origin adaptive search recovers two correctly classified folds

**Both previously excluded structures now have complete three-candidate pools
and correct Ca-class calls.** Their original contexts and origin energies were
already valid; the exclusion came from requiring older terminal-only proposals
before attempting the newer adaptive search. The new adapter removes that
dependency while reusing the established selector, optimizer and scorer.

| Consumed source | New operational R | ΔR from static origin | Frozen-reference result |
|---|---:|---:|---|
| P12293, Ca-conditioned sample 4 | −405488.97865446034 | +4.56891138770 | Ca-supported, correct |
| Q60AR6, Ca-conditioned sample 1 | −405478.9378397729 | +4.34729991260 | Ca-supported, correct |

R and changes are model kcal/mol. Both mathematical and operational decisions
agree. The unchanged canonical-only minimal-pool reference is
`workspaces/adaptive_minimal_pool_20260923/REFERENCE_v1.json`; no threshold was
fitted to these folds. Both static origins were already correctly Ca-like.
This is **recovered adaptive coverage**, not correction of two wrong static
predictions or new independent biological validation.

Root owns the full225 overlay, preserving the original unavailable rows and
all other source results. The two-source collection is complete, with no
unavailable/failure denominator hidden. The separate Q60AR6 La-sample0 adaptive
iteration-limit failure was not rerun. These are two structural repeats of
existing PQQ Ca-class references, not direct same-assay affinity observations.

The completed parent overlay, `workspaces/adaptive_minimal_pool_20260923/COMPARISON_recovered_v1.json`,
reports **204 correct / 1 wrong / 1 inconclusive / 19 unavailable** among225.
On the same206 available sources the released scorer has202 correct /2 wrong /
2 inconclusive. The two recovered cases add coverage; they do not create a new
accuracy advantage over their already-correct released calls.

## What actually ran

All four searches completed using the unchanged common four-mode selector and
scaled SLSQP200 rule. P12293 used A/209 χ1/χ3/χ2 and A/335 χ2; Q60AR6 used
A/205 χ1/χ3/χ2 and A/331 χ1. Both metals shared each source's physical modes.
All origin coordinates, charges (Ca−2, La−1), singlet states, forces, mappings
and vacuum/ALPB components were verified and reused.

The searches took 18/16 iterations for P12293 Ca/La and 15/16 for Q60AR6.
Q60AR6 Ca's final candidate touches the allowed boundary; it is an admissible
finite candidate, not a claimed unconstrained minimum. No extra starts,
post-hoc clipping, geometry rescue or alternate physical rule was used.

All four cross-metal MACE evaluations and all16 new native GFN2 cells completed.
For each source both metals see the same origin/adaptive_Ca/adaptive_La pool.
P12293 selects adaptive_Ca for both metals; Q60AR6 selects its own metal's
adaptive geometry. Existing mathematical minima and 0.1 kcal/mol operational
origin retention were preserved. Missing required cells would invalidate the
pool; none were missing here.

Native GFN2 uses the reference-compatible primary NoAutostart/native-mixer,
300 K, MaxIter500 policy. The separate experimental continuation policy was
not imported. No new DFT, source folding, hydrogen preparation or water changes.
No production/default change or remote push.

## Measured cost and tests

| Job | Work | Elapsed | Allocation |
|---|---|---:|---:|
|1210389|4 searches +4 cross-MACE|49 s|32 CPUs, one H200,200000 MiB|
|1210399|16 native GFN2 + collection|56 s|64 CPUs,128 GiB, no GPU|

Total **5152 allocated core-seconds and49 allocated GPU-seconds**. The four
searches made63 actual MACE energy/force calls; four cross calls bring the total
to67. All requests and16 GFN2 attempts succeeded. Archived origins and each
candidate's generating-metal MACE result were reused; their historical cost
is not erased. Local setup/reporting is additional and unmetered. Peak CUDA
allocation was3,756,127,232 bytes; maximum GPU-worker host RSS1,725,988 KiB.
These are measured development costs, not a matched production latency claim.

Five focused actual-artifact tests pass after completion, zero skips. Before
execution the final molecular-result check was explicitly skipped. Geometry/
source checks, malformed actual-manifest rejection and score algebra are
separate from the67 MACE/16 GFN2 calculations actually executed.

The immutable run is `workspaces/adaptive_origin_recovery_20260923/run_v2/`.
`run_v1` was preparation-only; no calculation occurred there. The only later
source addition is a pure `score` reporting command; running snapshots remain
unchanged. Exact pins are in [ARTIFACTS.json](ARTIFACTS.json), with runnable
operations in [COMMANDS.md](COMMANDS.md).

Recommend retaining this direct-origin preparation route for the minimal
adaptive challenger: it preserves the physical rule and recovers supported
coverage. Wider integration/promotion remains a separate decision; the released
baseline remains unchanged.
