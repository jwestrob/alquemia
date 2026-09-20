# Three-fold PQQ scoring is now runnable as an opt-in readout

The existing standard scorer can now retain three declared folds and report their
median, rather than making the user's conclusion depend on choosing one fold.
A fresh, real source-to-score test on the prescribed Q9Z4J7 samples passed: all
three preparations and all 18 molecular endpoints completed, and the median
retains the expected Ca-associated call. The prior robustness result is useful
but limited: **native core matched context composite on the all-three-of-four
reference panel**. This integration establishes usability, not superiority over
that simpler method or new biological accuracy.

## What the readout shows

The operation is `affordable_workflow.py standard ensemble`. Membership must
appear in the original source request pinned before execution. It verifies three
distinct source hashes, actual native La, actual protein sequence and numbering,
identical site-role mapping/assembly, prepared core chemistry, checkpoint,
preparation and frozen reference compatibility. Every member is required. Missing
scores produce an unavailable median; no substitute fold or fallback is chosen.
The existing preparation/scoring code, release parameters, default and explicit
DFT reference remain unchanged.

Individual raw values, molecular components, failures, domain flags and differing
context compositions remain visible. The operation reports minimum, maximum,
spread and ordinary median; component medians are explicitly nonadditive. Its
band transfer is named **developmental** because existing single-structure bands
are being applied to a three-fold descriptor. There is no new threshold,
probability, thermal population, occupancy or binding-free-energy claim.

## Fresh integration outcome

The choice was fixed before submission: lexicographically first three
noncanonical La-conditioned Q9Z4J7 source IDs in the pinned reference inventory,
namely seed-1 samples 1, 2 and 3. All were previously consumed development cases.

| Source sample | Composite R, model-kcal/mol | Context atoms |
|---|---:|---:|
| 1 | −405466.96524348116 | 139 |
| 2 | −405469.6848696391 | 127 |
| 3 | −405488.87675400084 | 145 |

The median is **−405469.6848696391, Ca-supported** under unchanged bands. The
**21.91151051968336 model-kcal/mol spread remains prominent**, and all three
retain **unvalidated_input_domain** under the released canonical-source registry.
A median does not remove that structural/model uncertainty. Context membership
varies despite identical core chemical state; each context adds formal charge−1.
No unique physical explanation is inferred from the score spread.

Fresh context coordinates, including all hydrogens, are exactly identical to
the archived source preparations. Maximum score repeat difference is
6.40×10⁻¹⁰ model-kcal/mol. This is a fresh normal standard run, not reused energies
relabelled as an integration result. Complete receipts and original scalar values
are linked in [RESULT.json](RESULT.json).

## Execution and verification

Job1204055 completed six native MACE scalar evaluations and twelve native GFN2
vacuum/ALPB endpoints, zero DFT, with no scientific failures. Source preparation
plus scoring took117.67742664s; total allocation was120s on one H200 and32CPUs,
200000MiB RAM: **3840allocated core-s and120GPU-s**. Folding is excluded. Peak
reported native-worker host RSS was1822884KiB; peak CUDA allocation2789489152bytes.
Scheduler AllocTRES omits GPU entries here; the frozen wrapper and actual CUDA
worker receipts separately record GPU execution. Costs are integration receipts,
not a new matched-hardware comparison with DFT.

Ten targeted real-fixture tests pass, including actual fresh CLI collection,
missing scores, corrupted sequence/role/conditioning, changed bands/score/domain,
rejection of adding membership to an old request after execution, and an extra
actual protein-chain mismatch. The latter was added during review; re-reading the
completed genuine collection required no new molecular calls. Test inputs
with removed or altered fields are explicitly corrupted copies of real artifacts;
no successful scientific output was fabricated. The scoped release/ensemble suite passes23tests
([log](RELEASE_AND_ENSEMBLE_TESTS.txt)). A broader wildcard additionally imported
unrelated `test_pqq_microstates.py`, which could not run because this environment
lacks pytest. That import failure is preserved in REGRESSION_TESTS_FINAL.txt; no
dependency or scientific executable was substituted.

Both separately prepared three-source PLM requests pass read-only schema,
sequence/site-role and native-La validation. This ran **no preparation or molecular
calculation** and supplies no PLM classification. PQQSEQ_07ab500e3df76b30d71c sample2's
prior failed single-fold admission remains included; input identity validation
is not a new admission decision. See [input-only record](PLM_INPUT_VALIDATION.json).

## Recommendation

Use the new operation as an explicit developmental readout when three supported
La-conditioned folds have been declared. Preserve raw spread and incomplete
coverage in interpretation. Do not promote an ensemble default or rescore PLM
from this integration test. Prior native-core parity and the lack of independent
biological validation remain material limits.

[Commands and request schema](COMMANDS.md), [frozen plan](PLAN.md),
[complete request](REQUEST.json), [actual allocation](SACCT.txt), and
[earlier all-triples result](../accommodation_goal_20260920/THREE_FOLD_REPORT.md).
