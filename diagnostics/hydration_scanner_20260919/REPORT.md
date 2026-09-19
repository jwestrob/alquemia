# Context-prepared MACE improves alpha/GGR discrimination

**The existing masked MACE scorer improves from 2/6 to 6/6 expected alpha/GGR
structural orderings after contextual water preparation.** The benefit previously
seen with DFT also works with MACE scoring. These six comparisons represent one
condition-qualified biological comparison, not six independent validation cases.

| Scorer | Original preparation | Matched normalized water geometry | Context-prepared water H |
|---|---:|---:|---:|
| DFT, original scoring cores | 0/6 | 0/6 | 6/6 |
| Masked MACE, original whole proteins | 2/6 | 2/6 | 6/6 |

All DFT results are archived compatible calculations; 12 new MACE calls ran.
The geometry proposer also uses MACE, so this route needs no per-site DFT once
preparation is available. Preparation has not yet been generalized and tested
across arbitrary water-bearing proteins. Production remains unchanged.

## Matched structural comparisons

Positive alpha-minus-GGR is the supported relative direction. Values below are
MACE descriptor model-kcal, not measured binding free energies or DFT energies.

| Alpha | GGR | Normalized starting waters | Context-prepared waters |
|---|---|---:|---:|
| 1F6S | 1GLG | 12.772 | 40.799 |
| 1F6S | 2FW0 | −8.124 | 19.902 |
| 1F6S | 2FVY | −9.216 | 18.811 |
| 6IP9 | 1GLG | 3.850 | 40.958 |
| 6IP9 | 2FW0 | −17.047 | 20.061 |
| 6IP9 | 2FVY | −18.138 | 18.970 |

The two alpha structures' MACE scores differ by 8.923 before contextual
preparation and 0.159 afterward. This supports improved consistency on this
consumed pair; it does not establish robustness for every conformation or protein.
PQQ and GGR contain no retained site waters in this benchmark. Their inputs are
exact identity operations, preserving the existing scorer and 25 canonical PQQ
reference calls plus two supported crystal transfers. No fresh PQQ energies ran.

## Preparation and energy accounting

Transfer the previously generated native OMOL water-H proposals by source identity
into the existing whole-protein inputs. Protein, cofactors, selected metal,
water oxygens, charge, protonation and water inventory remain fixed. No classifier
weights or thresholds were fitted here.

The existing whole-protein water normalization differs slightly from the DFT
water reference (up to 0.0102 Å). Four extra bound MACE calls supply exactly matched
normalized water inputs, separating normalization from contextual orientation.
The substantial effect survives this matching.

The frozen four-state interaction descriptor is:

`R = (E_Ca,bound − E_Ca,detached) − (E_La,bound − E_La,detached)`.

Endpoint-specific water geometry prevents cancellation of identical detached
protein states. The contextual arm therefore runs four endpoints per structure.
Its detached-environment contribution is −2.113/−2.664 model-kcal; it is not
silently assigned zero. The normalized arm has identical paired nonmetal geometry
and uses the previously established exact two-call cancellation. A separate
atom-reference-corrected bound-total contrast is retained for component analysis;
it was not selected as the score after inspecting outcomes. Component closure
errors are below 8e−15 model-kcal.

## Cost and verification

Job **1202084** completed all 12 MACE calls in **97 wall/GPU-allocation seconds**
on one H200 with 16 CPUs: **1,552 allocated CPU-seconds**. Actual model forwards
total 10.558 seconds; process startup and validation overhead remain included in
97 seconds. Peak allocated GPU tensors: 3,125,964,288 bytes. Maximum worker RSS:
1,579,548 KiB. The host memory request was 200,000 MiB under the existing policy.
No new DFT, optimization, or scientific retry ran.

An unlimited scheduler-duration request initially remained pending due to QOS.
The existing seven-day scheduler ceiling was restored without resubmission or
compute. This is a scheduler requirement, not a project compute budget.

The prior MACE proposal cost of 96 GPU-seconds for both structures is separately
recorded and reused. This full validation plus that preparation totals 193
GPU-seconds, but is not a measured production-throughput benchmark: four control
calls are additional, and login preparation/tests were not individually metered.

Four real-fixture tests pass with zero skips, including actual 12-endpoint
integration, source/scaffold preservation, corrupted mapping rejection, archived
sign/unit algebra, and refusal to use missing detached terms without proven
cancellation. Unsubmitted v1 and the preexecution partial collection are preserved.

## Interface and next work

New opt-in protocol: `masked_omol_context_prepared_water_comparison_v1`.
Code: `scripts/hydration_scanner.py` and six dispatch lines in `mace_omol.py`.
Unrounded values and receipts: [RESULT.json](RESULT.json).
Runnable operations: [COMMANDS.md](COMMANDS.md).

No broad validation, physiological occupancy or absolute affinity decision is
claimed. The active goal continues with PQQ-specific context and response, and
integration into the reusable scanner.
