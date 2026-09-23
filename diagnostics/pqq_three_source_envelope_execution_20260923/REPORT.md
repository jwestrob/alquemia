# Two real PLM triples run successfully through the candidate

**The practical three-source envelope route now completes end to end on both
declared PLM proteins.** All six sources are available. Accommodation substantially
narrows their within-group raw-score ranges, but the candidate predicts
Ca-supported medians for both groups. Their biological labels remain unknown;
this is operational and structural-repeat evidence, not an accuracy claim.

| Group | Same-context origin range | Accommodated range | Accommodated median R | Prediction |
|---|---:|---:|---:|---|
| PQQSEQ_07ab500e3df76b30d71c |40.138064|2.307280|−405465.300734|Ca-supported|
| PQQSEQ_83440678cbbd658047c9 |72.416250|17.230067|−405469.424200|Ca-supported|

Ranges and R are modelkcal/mol. Group results require all three declared members.
These are the same fixed4.3Å pocket/state before and after accommodation, so the
changes are internally comparable. Historical DFT/released scores were neither
subtracted from these values nor altered. No PLM expression/biological joins were
rewritten by this branch.

## All six source results

| Protein suffix / AF3 sample | Origin R | Accommodated R | ΔR | Prediction |
|---|---:|---:|---:|---|
|07ab /0|−405503.182903|−405466.327442|+36.855460|Ca-supported|
|07ab /1|−405527.002750|−405464.020163|+62.982587|Ca-supported|
|07ab /2|−405486.864686|−405465.300734|+21.563952|Ca-supported|
|8344 /0|−405532.899636|−405469.424200|+63.475436|Ca-supported|
|8344 /1|−405482.180027|−405458.593430|+23.586598|inconclusive|
|8344 /2|−405554.596277|−405475.823497|+78.772780|Ca-supported|

Mathematical and operational choices coincide for every source. Each metal
selects its own proposal from the shared three-candidate pool. Selected composite
Ca work ranges−1.150 to−58.944kcal/mol; La work ranges−24.736 to−137.717kcal/mol.
Full native/solvent work components, both metal-row choices, unrounded contrasts
and all matrix cells remain in the result. No prediction was used to choose a
source, geometry, threshold or extra search.

All12 searches converge under the declared finite SLSQP criterion and pass final
physical checks. Seven proposals touch an angular or physical-displacement bound.
They are admitted finite candidates, not unconstrained minima, thermal ensembles
or complete binding free energies. The larger remaining8344 spread and its one
inconclusive source remain visible.

## Scope, validation and cost

Job1211626 completed0:0 on node-224-2t-8gpu-1, with32CPUs, oneH200 and200000MiB
host memory. Actual calls:12 MACE origins +159 search evaluations +12crossMACE
=183, and72 strict nativeGFN2 scalar endpoints; all complete. No molecular
failure, retry, new DFT, fold, protonation or archive energy/force reuse occurred.
The two exact pre-existing source-protonation caches were reused; within-run
origin/own-proposal reuse is explicit.

The outer allocation took251s:8,032allocatedcore-seconds and251requestedGPU-
seconds, including scalar phases while the GPU was reserved. The executor itself
records212.204s. Step accounting reports784.162CPU-seconds; largest batchRSS is
7,256,416KiB. Native worker peakCUDA allocation was6,429,156,864bytes. Three warm
workers each loaded once (1.177/1.193/1.195s); their measured intervals were
5.679s(origin),53.424s(search),8.245s(cross). Twelve searches sum52.097s. Scalar
attempt walls sum2329.623s with31.106s median under concurrent execution. These
nested times are not additional allocations or isolated hardware speed claims.

Separate coordinate-preparation functions total16.955s, excluding historical
protonation, imports/request validation and final dry-run. Additional local
implementation/tests/reporting overhead is unmetered. This is one batched
six-source/two-protein measurement, not guaranteed independent per-protein latency.

Seven final real-artifact tests pass in12.717s, zero skips. They cover exact
source/state/task scope, qualified scalar input and actual printed tolerance/rank,
unchanged force-selected modes using compatible archived warm receipts, direct
score algebra, complete matrices, strict three-member summaries and final physical
admission. An initial preflight test used the older legacy native-force receipt
format and failed at its parser; the test was corrected to an actual compatible
warm-executor fixture. No scientific output or method changed. Preparation-only
tests separately passed six checks.

## Interpretation and next step

The scorer uses its actual4.3Å envelope canonical25 reference
SHA4a1d3ae6e83b5b6ace408f729c39390d44f673554e338df932fd38fc971f67d9.
That calibration retains25canonical+3crystal calls, but A0A3Ca3 remains an
inconclusive benchmark probe and full structural-transfer qualification is a
separate task. The successful PLM run does not erase this limitation or supply
known biological labels. Keep this as an experimental opt-in route; preserve the
released/DFT defaults and hold cohort rescoring.

Authoritative result:
`workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RESULT_1211626.json`,
SHA06e2c0b7727ff0aa53e6085b144a307c2eb6ca71cd0cc8f06112e3b8ec9e4be7.
[Commands and schema](COMMANDS.md) provide the existing preparation, preflight,
collection/report and exact allocation receipt. [Artifact pins](ARTIFACTS.json)
include all final results, tests and accounting. No further molecular work is
required to attach these separately named candidate fields to the existing PLM
protein export; that join is owned by root.
