# Tighter native stopping does not remove the donor-response solvent discrepancy

Completed 2026-09-23. **The strict native solver changes these displacement contrasts by at most0.080820 kcal/mol; the solvent-added disagreement with archived DFT remains.** All36 logical cells are available:32 new fresh-SAD nativeGFN2 calls and4 exact strict32 origin reuses. No new MACE, DFT, geometry search, force evaluation, calibration or production change.

## Complete differential comparison

Work is E(q)−E(0) at each metal. The differential is work(Ca)−work(La); positive means the displacement stabilizes La more. The two PLM proteins remain biologically unlabeled. Values are kcal/mol (MACE/composite model energies).

| Source | Angle(rad) | Archived nativeDFT | NativeMACE | Old composite | Strict composite | Strict−old |
|---|---:|---:|---:|---:|---:|---:|
| 4MAE | −0.2 | −1.476170 | −1.685754 | −1.613440 | −1.598234 | +0.015206 |
| 4MAE | +0.2 | +1.087421 | +1.013542 | +0.948262 | +0.973700 | +0.025438 |
| PLM8344 | −0.2 | −32.158738 | −30.738355 | −28.490526 | −28.487349 | +0.003176 |
| PLM8344 | +0.2 | +29.867826 | +28.535393 | +26.957148 | +26.960435 | +0.003287 |
| PLM07ab | −0.2 | −20.354831 | −19.838034 | −18.039220 | −18.119886 | −0.080666 |
| PLM07ab | +0.2 | +20.335087 | +18.916561 | +17.674404 | +17.593584 | −0.080820 |

All12 endpoint-work signs and all6 differential signs agree with archived DFT. The largest endpoint absolute error is0.912009 for MACE,3.151856 for old composite and3.155099 for strict composite. Maximum differential absolute error is1.420383,3.668213 and3.671389 respectively. The strict solvent correction still worsens all four PLM differential errors relative to MACE alone; it slightly improves the 4MAE comparison. Exact values and every endpoint—not only summary extrema—are in [ENDPOINT_WORKS.csv](ENDPOINT_WORKS.csv) and [DIFFERENTIALS.csv](DIFFERENTIALS.csv).

The largest individual scalar shift is−0.081731 kcal/mol for PLM07ab Ca origin vacuum. Its changed origin reference shifts both displacement works similarly; it is not a new donor-geometry effect. All36 component changes, actual SCF cycles and output hashes are retained in [CELL_SHIFTS.csv](CELL_SHIFTS.csv).

## Scientific interpretation

These results exclude the tested stopping change as an explanation of the earlier multi-kcal solvent-added discrepancy on these particular paths. They do not establish that solvent physics is unnecessary, nor identify a unique remaining cause. Native DFT uses r2SCAN-3c/CPCM; the composite uses OMOL plus GFN2(ALPB−vacuum), with finite electronic smearing. Different electronic and continuum models remain genuine differences.

The original experiment froze no general response-accuracy threshold. Accordingly this report assigns no new accuracy pass/fail threshold, force qualification, affinity correction or classifier calibration. The native stopping comparison also uses the independently qualified one-rank scalar execution profile versus the archived eight-rank outputs. The frozen geometries, state and parameter exports are unchanged. No outcome-dependent starts, selected retries or discarded cells were used.

**Recommendation:** retain the tighter scalar protocol's independently demonstrated numerical utility, but do not claim it repairs composite donor-response physics. This contained diagnostic is complete; further calls on the same paths are not warranted by this result.

## Execution, source reuse and checks

Sources are exactly the original4MAE and PLM `PQQSEQ_83440678cbbd658047c9`/`PQQSEQ_07ab500e3df76b30d71c` sample0 contexts, allq=0/−0.2/+0.2 extra-Asp torsions, both metals/media. The four reused4MAE q0 cells have identical XYZ hashes, state, exact strict recipe, fresh-SAD initialization and valid one-rank receipts. Archived MACE and all18 nativeDFT point states/coordinates were checked before execution; their energy-derived works reproduce the pinned original comparison. The new manifest is `workspaces/strict_donor_response_20260923/run_v1/manifest.json`, SHA256 `34a7a11814d42543cf2c0912e21eff71897786f4da1393f675af0de2d9e39193`.

All36 outputs confirm native mixer, fresh SAD, effective TolE1e−10, one rank, matching parameters/electron count/charge and normal SCF completion. SCF took24–115 cycles. Seven outputs exceed at least one printed density threshold; these exceptions remain explicit. The native mixer has its own stopping behavior, so ordinary printed density thresholds do not establish failure of that mixer, but scalar success does not qualify gradients or density convergence.

Six actual-artifact tests pass in2.597s: all18 physical-map replays, exact source recipe/coordinates, actual native receipts/state and raw-output reconstruction of all12 works/all6 differentials. No synthetic scientific data or extra molecular calls were used by tests. The first execution and collection succeeded; there were no failed calls or retries.

Job1211802: **38 allocated seconds×32CPUs=1216 allocated core-seconds**, zeroGPU allocation. Engine execution33.752742s;32 molecular task durations sum615.708289s (median18.206044s, maximum32.784718s). Task-duration sums overlap in time and are not allocated CPU time. Four prior reused calls and local preparation/report CPU are separate historical/local costs. `COSTS.json` retains actual scheduler and executor receipts.

[Commands](COMMANDS.md), [frozen plan](PLAN.md), [summary](SUMMARY.json), [artifact pins](ARTIFACTS.json). The immutable collection SHA256 is `842dc2969beee310b31268f767fdba6b5eba827ee3cf8338c502e62ec91deb78`.
