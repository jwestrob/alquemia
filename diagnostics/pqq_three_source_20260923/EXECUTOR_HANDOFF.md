# Proposed thin executor handoff: the two explicit PLM triples

Status: implementation handoff only, no molecular launch. Root requested this
plan after the preparation API was sealed. Scientific execution waits for the
full 100-triple transfer verdict and root's finite-plan confirmation. This is
command/usability qualification on two already consumed, experimentally unlabeled
proteins; no accuracy test, source selection, PLM cohort rescore or preferred-sign
target. Existing production and tenfold APIs remain unchanged.

## Exact inputs and physical eligibility

Use the two immutable preparations listed in [EXECUTOR_INPUTS.json](EXECUTOR_INPUTS.json):
07ab and 8344, each the already declared AF3 samples 0, 1 and 2. The six contexts
contain 168 or 190 atoms. All twelve paired q0 geometry checks pass the existing
`accommodation_torsion_profiles.geometry_check`, including source-connected cap
checks. This read-only check used no molecular calculator. The old 07ab sample 2
single-fold admission record remains historical evidence; the present result is
for the exact new three-source union state, not a deletion of that old failure.
Actual selection of four independent modes still requires fresh paired forces.

Reuse only the exact completed source preparation, union and maps; no new
protonation. No archived energy or force supplies this proposed integration.
Each group retains its declared lexicographic physical-state anchor and fixed
three-member union. The anchor is neither a label nor a calibration member.

## Small implementation change, existing kernels

Add a separate `pqq_three_source_execution.py` adapter, preserving the sealed
preparation module and previous tenfold execution interface. Its input is a
finite list of explicit prepared triple pins, reference, qualification pins and
agreement. Expose `prepare`, `dry-run`, `execute`, `collect`, `report`; no parallel
workflow framework or new molecular backend.

Adapt the established private `pqq_adaptive_candidate` engine as the tested
`pqq_union_execution.candidate_engine` already does. Supply generic plan checking,
exact prepared-source handoff and mapping lookup instead of its A0A3/historical
canonical guards. Use unchanged `origin_manifest`, `selected_tasks`,
`proposal_manifest`, native warm worker, private precision optimizer,
`candidates_manifest`, scalar low runner and pool collector. Check actual stage
source/method/state membership before every operation. Pin implementation and
SciPy wrapper/kernel identity. Stage directories and completed same-run receipts
remain restartable only under identical pins; incomplete/failed work stays
visible without automatic scientific retries or alternative sources.

The recipe stays native OMOL plus native GFN2(ALPB water minus vacuum),
MaxIter500, one rank per scalar cell, 300 K. Use the existing four-mode selector
from the paired native force pattern, one SLSQP start per metal, 200 iterations,
ftol 1e-8 in the existing Hartree-equivalent objective, angular bounds ±0.8 rad,
maximum heavy displacement 0.8 Å and existing admission tolerances. No gradients
from GFN2, extra modes, new starts, water/proton changes or numerical tuning.

Each source pool is exactly `{origin, adaptive_Ca, adaptive_La}`. Cross-score both
proposals with both metals. Reuse q0 within that run and each proposal maker's
native MACE energy; do not call these reused cells new calculations. The complete
matrix is mandatory. Missing origins, insufficient modes, failed search, invalid
geometry or failed required energy cell leaves that source unavailable. Preserve
origin results separately; do not substitute them as an accommodated score.

## Exact finite scientific scope

|Stage|New calls across six sources|
|---|---:|
|Native MACE paired origin energies and analytic forces|12|
|Native GFN2 origin vacuum/ALPB scalars|24|
|Bounded native MACE searches|12 starts|
|Opposite-metal native MACE at the two final proposals|At most 12|
|Native GFN2 at both proposals, both metals/media|At most 48|
|Total native GFN2|At most 72|
|DFT / folds / protonation|0|

Search forward/gradient evaluation count depends on the frozen optimizer; record
all actual calls, including failed attempts. Twelve searches do not mean twelve
MACE forward evaluations. Exact within-run geometry reuse may lower counts only
when actual pins/state match. No additional starts or rescues are implied.

## Output and frozen reference

Use exactly the prepared threefold reference, SHA256
`e725ea836d4a76bff8476cdedbc12be5d6a8bc5595da092d312f5dcec5c4ac9c`,
`Nikasha_three_La_membership_ftol1e8_canonical25_v1`. Do not recalibrate on PLM.
Retain per-source unrounded origin R, accommodated mathematical and operational
R, component energies, selected candidate identities, endpoint work, physical
boundary flags, termination and failure reason. Here origin R is the same-union
unaccommodated comparator, not the released local-context score.

Use the existing operational 0.1 kcal/mol candidate-selection tolerance. Compute
strict median R across all three declared available source scores, separately
for the mathematical and operational variants; one unavailable source makes the
group summary unavailable. Apply each frozen variant's own bands. Report all
three source scores and spread even when the median is available. All six source
labels and both group labels remain unknown, so output classification is a
prediction with no accuracy numerator, affinity-probability or improvement claim.

## Resources and expected cost

One finite 32-CPU, one-H200, 200000-MiB allocation, with the established warm GPU
stages and up to 32 concurrent one-rank scalar cells. Use an explicit permitted
host and existing MPI/thread setup. The matched ten-source A0A3 candidate took
358 allocated seconds, including 114.933 s fresh preparation, 20 searches/272
search forwards and 120 scalar calls. Its source:
`workspaces/pqq_union_execution_20260923/COSTS_v2.json`.

A simple six-versus-ten task-count scaling of the remaining allocation is about
146 s, only an engineering reference. These PLM contexts are different and some
are larger; expected scale is minutes, not a guaranteed duration or a limit.
Existing union preparation is reused and timed separately. Actual allocation,
step CPU, GPU allocation, model-load, search, scalar and collection costs must
be reported; nested stage timings must not be added as independent allocations.
No molecular submission or executable manifest has been created by this note.
