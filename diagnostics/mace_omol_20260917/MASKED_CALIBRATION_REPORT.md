# Two-call masked MACE research classifier

The canonical25-case calibration passes, as do both valid crystal transfers.
All189full-panel factorization checks pass; maximum difference2.357415596e-8
model kcal vs the predeclared0.01 tolerance. No new inference was needed.

The two-call numeric bands are Ca<=42.346663801320496 and
La>=51.455254788976355. The gap is9.108590987655859. The four-call bands remain
immutable in CHARGE_ABLATION_CANONICAL_REPORT.md. Both numeric recipes agree
to roundoff; neither inherits the ORCA reference, a physical zero or a generic
La/Ca affinity threshold. These are kcal-equivalent learned descriptor units.

The opt-in report accepts --calibration, checks the actual model, software,
source/reference pins, factorization ID, PQQ3minus/27atoms/zero-water state,
chainA assembly, preparation forcefield and singlet policy. Incompatible or
unavailable calibration yields no class; valid raw values remain reportable.
Non-PQQ inputs cannot inherit these PQQ functional-association bands.

Authoritative reference: workspaces/mace_omol_20260917/masked_calibration_v1/reference.json
SHA217a127f7ba8f6cccb59d1c63b2d389e9022449a5feb09f58eba2e3b04cba54d.
Frozen reporting closure: masked_calibration_source_v1/implementation/.
Compact reference/check/score record: MASKED_CALIBRATION_RESULT.json.

1KB0 remains unsupported and in the3-transfer denominator; the all-three
transfer gate is false. None of these cases is newly blind. Canonical composition
already separates this panel. Alpha/GGR is one separate qualified development
comparison, not broad affinity validation. Whole proteins are outside published
OMOL training size, and masking a learned charge feature does not produce a
validated physical potential. Production baseline/default unchanged.

Cost: canonical job1200830 took3546GPU-allocation seconds for100newforwards;
local report115.17wall seconds,115.02CPU seconds,544552KiBRSS. Actual saved bound
pairs across25proteins took49.85685–60.772997model-evaluation seconds(median53.068669),
or53.744777–65.495802worker-wall seconds(median57.267661). These exclude controller,
source preparation and reporting; they are measurements of those exact calls,
not a separately executed two-call production-throughput benchmark.

Actual prepared-input report replays completed:1H4I gives17.57995558550866,Ca;
GGR gives23.98742145552895 with no class and an explicit non-PQQ-scope status.
Both reuse actual endpoints without new MACE/DFT/solver calculations. All five
actual-artifact guard tests passed in675.322seconds, none skipped: unavailable
transfer retention, corrupted-reference rejection, actual PQQ decision replay,
wrong numeric-recipe rejection, non-PQQ scope rejection and native-model
reference rejection. These are parser/algebra/compatibility tests using actual
scientific outputs; no new model evaluations occurred. The research goal
continues; this is a useful candidate milestone, without automatic promotion.

Validation receipt: workspaces/mace_omol_20260917/masked_calibration_validation_v1.json.
Standalone PDF/SVG/PNG figure: masked_calibration_figure_v1/canonical_masked_mace.*
under the same workspace. Source pins, rendering code and cost receipt retained;
no additional scores or fits were produced. Figure visually inspected.
