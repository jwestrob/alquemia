# Classical derivative question closed: quarter-step passes unchanged tolerances

**All four classical components pass the original analytic-versus-finite-difference tolerance at±0.0125degree for both target assignments.** The prior coarse and halfstep failures remain unchanged. This closes the diagnosed classical finite-step question on this real mode; it does not qualify the missing electronic/full hybrid derivative.

Before execution, O(h²) scaling predicted MM-LJ residual.0129420218kcal/mol/radian. Actual residual is.0129405472, below the unchanged.0270116631 tolerance. The halfstep/quarter-step residual ratio is4.000456, as expected for central-difference truncation.

| Component | Actual quarter-step residual, kcal/mol/rad | Original tolerance | Result |
|---|---:|---:|---|
| Retained bonded |4.17e−9|.0200000|pass|
| MM LJ |.01294055|.0270117|pass|
| MM Coulomb |9.30e−6|.0218997|pass|
| RealQM–MM LJ, La |.000790963|.0205188|pass|

Dy has the same MM residuals and essentially identical realQM–MM LJ residual. Half/quarter Richardson discrepancies are−1.97e−6 for MM-LJ and−1.59e−8 for La crossLJ; these remain separate diagnostics rather than new acceptance thresholds. No further step, quantum call or optimization ran.

## Exact execution and receipt status

Job1220331 performed16component energy/force queries: two prescribed configurations, two targetLJ assignments, four components. Actual cost2wall seconds/2allocatedCPU-seconds on1sharedCPU mem0, noGPU; script measured1.847s. All four per-configuration receipts were written immediately, including source IDs, exact coordinates, component energies and full Cartesian forces with explicit units/sign. Complete `RESULT.json` then persisted successfully using the scoped NumPy-safe writer.

A redundant final stdout `json.dumps` failed **after all scientific artifacts were written**, so Slurm reportsFAILED/exit1. This process failure is retained, separate from the passed numerical tests. No molecular rerun occurred. Every saved configuration's energy matches the aggregate, all force arrays are finite1891×3, and the complete saved JSON was rendered byte-exactly in a read-only follow-up.

The immutable executed `QUARTER_STEP.py` remains unchanged. Replacement `QUARTER_STEP_SAFE.py` prints the already serialized file and requires an explicit unused output directory for any future execution. Its `--render` branch was tested on the actual result; its molecular branch was not executed. The shared production writer and other agent's originalA/B assembly remain untouched.

## Artifacts

`workspaces/metal_environment_response_20260926/coupled_scaffold_classical_v1/quarterstep_v1/`:

- `RESULT.json`: actual numerical pass and original comparisons.
- `{La,Dy}_{minus,plus}_receipt.json`: actual component energies/full forces, force sign explicitly retained.
- `{minus,plus}_physical.json`: the two physical configurations.
- `VALIDATED_RECEIPTS.json`: exact energy/force/render checks and pins; no molecular calls.
- `RENDERED.json`: byte-identical rendering recovery.

Scheduler evidence: `QUARTER_ACCOUNTING.txt`. Across this complete classical branch including earlier receipt failure/replay and pair audit, cost is19allocatedCPU-seconds; this quarter-step adds only2. No subsequent compute is running or required to establish this result.

The next scientific question belongs to the complete electronic-plus-classical response, with existing solvent and parameter-transfer limitations unchanged. Passing this classical check does not establish La/Dy affinity discrimination or justify unconstrained whole-protein relaxation.
