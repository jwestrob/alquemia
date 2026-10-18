# Mex consistent-grid derivative test: final result

18 October 2026 UTC. **The candidate fails its frozen energy–gradient consistency gate and is closed.** All three native calculations pass their original admission checks. Saved collector 14107 returns `failed_derivative_gate`; its scheduler exit1 reflects the scientific gate, not an incomplete collection.

## Result

The analytic energy slope is compared with the slope from the two displaced energies, using the exact original O97 direction and h=0.001 Å. These are gradient projections; forces have the opposite sign.

| Contribution | Analytic slope | Energy-difference slope | Difference | Allowed magnitude | Result |
| --- | ---: | ---: | ---: | ---: | --- |
| Embedded PBE0 | −3.673645578 | −3.797017470 | −0.123371893 | 0.005379702 | Fail |
| MM ledger | +17.752465327 | +17.752485207 | +0.000019880 | 0.006775249 | Pass |
| Fixed-charge ALPB | −5.530642594 | −5.530641591 | +0.000001003 | 0.005553064 | Pass |
| Total | +8.548177155 | +8.424826145 | −0.123351010 | 0.005854818 | Fail |

All units are kcal/mol/Å. Gate: `abs(FD−analytic) <= 0.005 + 0.0001*max(abs(FD),abs(analytic))`. The exact original comparator, threshold and printed-string energy subtraction are unchanged. All three result slots are available; there was no survivor selection. The collector made zero molecular calls and reused the admitted MM/ALPB terms exactly.

## What this establishes

The absolute PBE0 discrepancy decreases from the historical 0.414201 to 0.123372 kcal/mol/Å, but matching the late-SCF and final COSX grids does not restore consistency under this recipe. The remaining failed contribution is embedded DFT; the MM and solvent implementations pass this particular numerical check.

This is a consistency result, not physical validation of MM or solvent. It does not prove a unique ORCA bug or identify the remaining cause. SCF nonstationarity, electronic-root effects, other integration effects and finite-step uncertainty are not resolved. The new central also uses the matched 12-rank runtime, whereas the historical central used 112 ranks; the improvement cannot be assigned exclusively to the grid change.

The result does not repair the separate LBT scaffold discrepancy: whole-PBE0 plus the same solvent term charges +12.28/+11.72 kcal/mol for the La/Dy deformation, while the hybrid claims a La gain of 28.24. Whole-MACE agrees in direction with whole-DFT. The terminal MM-dominated forces remain the concrete mechanics target; moving the tested partition boundary alone did not remove them. See the [whole-LBT reference](../../../lbt_whole_reference_continuation/FINAL_REPORT.md) and [partition-force comparison](../../../lbt_partition_force_reference_v1/REPORT.md).

## Native evidence and limits

All cells use the exact admitted 209-atom electronic core, 1840 external charges and matching 12 MPI / 48000 MiB / MaxCore 3100 settings. The fixed field direction moves only O97; no new structure, state, parameter or solvent model was introduced. PBE0-D4/def2-TZVP, the La ECP46 representation and 832 explicit electrons remain pinned. Actual effective COSX grids 2 and 3 match in all cells.

| Cell | Job | SCF cycles / stop | Printed residuals passing | Final native energy / Eh |
| --- | --- | --- | --- | ---: |
| Central |13281|46 / Gradient Check|orbital gradient only|−5636.761994353496|
| Plus |13596|46 / Energy Check|energy change and orbital gradient|−5636.762000206363|
| Minus |13597|34 / Energy Check|energy change only|−5636.761988104498|

Every cell retains `all_printed_residuals_pass=false`. Native convergence does not establish a stable electronic state or a force error bound. Early negative-gap warnings remain recorded; they are not evidence of a final negative gap. Central and plus printed COSX/D4 scalar closures pass; minus COSX closure exceeds its printed-rounding envelope while its D4 closure passes. This is a retained diagnostic, not a new post hoc acceptance rule.

Printed-energy quantization contributes a reported slope scale of 3.14e−7 kcal/mol/Å; binary64 operand scale is 5.71e−7. The collector records SCF sensitivity proxies, explicitly not error bounds. No step-halving, global-force, physical-state or La/Dy-affinity qualification follows.

## Decision and preservation

Close this finite candidate without further native calls, solver/grid sweeps, step changes, tolerance changes or optimization. The original failed protocol remains failed; its remaining six qualification cells and eight response cells stay stopped. No production default changes.

All three fresh native results, full central assembly, exact source/config pins, independent reviews and failures are preserved. Root owns integration. The broader mechanics goal remains unmet. A next mechanics task should target the unsupported hybrid scaffold response against the existing whole-peptide reference, with a separately declared finite comparison; this failed numerical test does not authorize an automatic continuation.

## Cost

| Job | Purpose | Elapsed seconds | CPUs | Allocated CPU-seconds |
| --- | --- | ---: | ---: | ---: |
|13281|Central gradient|41475|12|497700|
|13596|Plus energy|39481|12|473772|
|13597|Minus energy|32661|12|391932|
|13279|Preserved failed preflight|10|1|10|
|13280|Successful preflight|10|1|10|
|14107|Saved final comparison|13|1|13|

Molecular allocations total 1,363,404 CPU-seconds; listed support adds 33, for 1,363,437. Allocation cost is not measured process utilization. Prior failed-recipe and saved-diagnosis costs remain separate. No new molecular calls were made by the final collector.

## Reproduction and evidence

- [Execution specification](../EXECUTION_SPEC_V1.md), including frozen settings, gate and finite stop.
- [Full result](RESULT.json), SHA256 `4aa50af0ea9ab26219306aff4b2380ab64f38b01282017f19a7cc98d0d782741`.
- [Configuration](CONFIG.json), SHA256 `ecee028c850f1939d3f02c66de784e697cb40c665793f36e36c1f9ccc0c3ece6`.
- [Submission](SUBMISSION.json), [root review](ROOT_REVIEW.json), [independent review](INDEPENDENT_REVIEW.json), [cost ledger](COSTS.json).

The pinned collector source is unchanged: `1600a0a55635d9c592295cc13cd24d45c2911591889875f4fe06dee363348101`. Its exact Slurm command is in SUBMISSION.json. Do not rerun completed chemistry to reproduce this saved comparison.
