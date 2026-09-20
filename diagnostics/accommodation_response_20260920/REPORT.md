# Solvent-aware physical response: completed

**Native analytic ALPB gradients work and materially change some donor responses.**
All 60 calls completed: 28 native GFN2 analytic-gradient endpoints and 32 fixed
physical displacement singlepoints. All 12 predeclared derivative/refinement
checks pass. Fourteen exact-source OMOL force arrays were reused. No DFT or MACE
calculation, optimization, training, threshold fit or default change occurred.

## What this establishes

The installed ORCA 6.1.1 explicitly executes native XTB SCF, coordination-number
and ALPB analytic-gradient drivers. Its `.engrad` atom order, source coordinates
and energies agree with the actual endpoint outputs. The resulting gradient is
that of the tested composite

`E = E_OMOL,vacuum + E_GFN2,ALPB − E_GFN2,vacuum`.

Analytic physical-mode derivatives agree with independent energy differences:
maximum error 0.00685356 kcal/mol/rad; maximum 0.005→0.010 rad refinement change
0.01468752. Each passes its frozen tolerance, `max(0.02, 0.005*abs(derivative))`.
The largest primary-to-tight center pair shift across all seven cases is
0.05199304 kcal/mol (2FVY), below the unchanged 0.20 tolerance; largest endpoint
transfer change is 0.05113551, below 0.10. These results support this native
numerical route; they do not establish uniqueness of the electronic solution.

## Scientific information

Whole-group torsions retain information that averaging separate donor radial
forces can hide. In 4MAE, the Glu172 chi1 Ca-minus-La derivative changes from
−1.937 to −19.638 kcal/mol/rad after adding solvent transfer; chi2 changes from
−4.177 to −18.532. The whole angular response vector changes appreciably
(vacuum/composite cosine 0.687). Extra Asp301 chi2 remains similar (+6.775 to
+6.498), so the solvent effect is not a uniform scale factor.

Across GGR's three consumed structural replicas, 12/13 corresponding angular
response signs already agree in vacuum OMOL; all 13 agree with the composite.
The differing Glu205 chi1 changes in 2FW0 from +3.311 to −2.271; the other two
replicas give −6.781 and −6.478. These are one protein group, not three new
biological observations. Earlier radial work used a different MACE potential and
smaller representation, so improved consistency cannot be attributed solely to
the new coordinate mapping or solvent.

The extreme 6IP9 Asp82 chi1 response persists (−106.099→−107.035). Solvent-aware
forces therefore do not cure every strained or poorly represented state. The
PQQ result motivates controlled whole-carboxylate energy profiles, with native
DFT displacement validation, rather than another unconstrained optimization.
There is **no demonstrated accuracy gain, corrected-affinity scalar, entropy or
occupancy estimate** in this diagnostic.

## State and interpretation limits

The seven consumed contexts are 1H4I, 4MAE, 1F6S, 6IP9, 1GLG, 2FW0 and 2FVY.
Both 1F6S and 6IP9 are alpha-lactalbumin; the last three are GGR replicas.
Original charges, source coordinates, water inventories, PQQ and complete polar
contexts are preserved. Metal-specific alpha water H remain distinct and frozen;
the compared donor-coordinate measure is common, not the complete Cartesian
state. Caps follow physical anchors through the existing analytic mapping.
ALPB remains a cluster continuum correction, not the full protein environment.
No exact-current-state native DFT gradients were available for a direct force
accuracy assessment. Older source-verified DFT/POLAR forces are descriptive
comparators only (`archived_response_replay_v1.json`).

## Execution and artifacts

Jobs 1203465 and 1203499 completed in 24 and 200 allocated wall seconds on 64
CPUs: **14,336 allocated core-seconds, zero GPU-seconds**. Summed Slurm CPU usage
is 7,094.714 seconds. Followup batch MaxRSS was 2,865,820 KiB; the gate's zero
MaxRSS record is unavailable, not a measured zero. Local preparation/testing
cost is not included in these allocation receipts. There were no scientific
retries. An initial preparation path error occurred before chemistry. Collection
replay was repaired for <1e-12 cross-host projection roundoff without recomputing
energies; raw gradients, energies and other receipts still compare exactly.

- Frozen plan and diagnosis: [PLAN.md](PLAN.md), [DIAGNOSIS.md](DIAGNOSIS.md).
- Final result: `workspaces/accommodation_response_20260920/result_v2.json`.
- Final analyzer closure: `workspaces/accommodation_response_20260920/analysis_implementation_v1/`.
- All 60 tasks and 14 reused force sources: `prepared_v2/design.json`.
- Original results/receipts: `prepared_v2/{gate,followup}/`.
- Actual scheduling costs: [COSTS_sacct.tsv](COSTS_sacct.tsv).
- Six real-fixture tests pass, zero skips; parser/algebra/mapping replay and the
  actual scientific derivative checks are distinguished above.

`result_v1.json` is retained as an earlier numerical analysis. Its live analyzer
path was subsequently strengthened to require the observed ALPB driver; v2 pins
the immutable final implementation and has identical numerical results.

**Recommendation:** use this verified analytic response to test a small physical
accommodation hypothesis; do not promote forces to an affinity score.
