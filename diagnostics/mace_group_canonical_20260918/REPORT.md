# Full PQQ calibration fails for the grouped vacuum model

All 50 new endpoints completed, with zero execution failures. The unchanged
medium typed-group model gives **106/154** expected cross-class orderings and a
class gap of **−213.31636746175354 model kcal**. These are pairwise comparisons
among 25 consumed references, not 154 independent observations. No calibration
bands or transfer classifications are available under the frozen rule.

The preserved fixed-core DFT baseline still separates 25/25 references with its
8.602932 kcal/mol gap. The earlier **whole-chain masked OMOL descriptor** also
separates 25/25, with its own 9.108591 model-kcal gap. That earlier model is not
a carved-core method, and this failed new model does not invalidate its result.
Neither canonical panel establishes broad affinity discrimination.

## Results and interpretation

Protocol remains `intact_POLAR_medium_typed31_group_vacuum_compatibility_v1`.
All 25 cases, labels, preparations, explicit waters and score rules were retained.
Ca p12293 has the largest Ca score (−405225.35800388106); La q88jh0 has the
smallest La score (−405438.6743713428). Their reversed ordering defines the gap.
These raw values retain the large shared metal-energy offset; no physical zero
or old aquo reference was introduced.

All 54 endpoint charge checks pass. Existing rigid-transform qualification and
the three failed GGR grouping checks remain recorded separately. The two cached
crystals have real scores but no classes because calibration failed. 1KB0 stays
unsupported in the three-transfer denominator: its full-protein preparation has
missing peptide connectivity. This is not a failed prediction or a TRO issue.
The supported crystals repeat calibration sequences and are not blind controls.

## Saved-output diagnostic: local information survives

The declared readout analysis used all 25 completed reference pairs and all 34
previous transfer sites, with zero new molecular evaluations or fitted parameters.

| Recorded component | PQQ min(La)−max(Ca), model kcal | Correlation with baseline S | Correlation with full Ca-endpoint charge |
|---|---:|---:|---:|
| Short interaction | +1.822274 | +0.982478 | +0.395237 |
| Electron | −1.141909 | +0.934658 | +0.429463 |
| Electrostatic | −253.537255 | −0.000015 | −0.880624 |
| Total | −213.316367 | +0.323577 | −0.692774 |

The short learned readout separates the PQQ classes, while adding the other
terms loses that separation. This identifies an algebraic source of failure;
correlation does not prove a charge bug or uniquely identify the missing physics.
The unassigned remainder varies by only 6.490154e−8 across the canonical panel
and 3.451714e−8 across the transfer sites. Its correlation is unavailable rather
than calculated from rounding noise. No new component classifier was fitted.
The earlier short-component study already failed alpha/GGR directions, so simply
dropping electrostatics is not an established broad solution.

## Execution and artifacts

Job 1201524: **3,540 GPU-allocation seconds**, **56,640 allocated core-seconds**,
3,719 reported CPU-seconds; 3,098.347655721009 summed inference seconds.
Peak allocated GPU memory 11,639,341,056 bytes; worker RSS 1,815,544 KiB.
New group preparation took 77.917176 wall / 76.674563 CPU seconds. Historical
physical preparations, four cached crystal endpoints and qualification are
additional; local preflight/tests/reporting costs are not fully CPU-profiled.
The failed initial preflight made zero molecular calls and remains archived.

Workspace `workspaces/mace_group_canonical_20260918/` contains immutable
`report_v1/result.json`, `readout_audit_v1/result.json`, `cost_v1/result.json`
and the inspected `figure_v1/comparison.{pdf,svg,png}`. RESULT.json provides
compact values and source pins. COMMANDS.md gives exact replay operations.
Actual parser/algebra/report tests are documented in VALIDATION.md; these are
distinct from the 50 executed GPU calculations. No new DFT or solvent calls.

Recommendation: retain the baseline and the separate masked-PQQ research path.
Do not expand or calibrate this grouped vacuum model by dropping inconvenient
cases. Its successful original crystal pair did not predict full-panel success.
The active MACE goal continues; production defaults and old references are unchanged.
