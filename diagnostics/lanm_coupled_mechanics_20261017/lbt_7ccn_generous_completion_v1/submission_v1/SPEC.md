# Single generous 7CCN–La completion attempt

## Question and fixed selection

Can the unchanged whole-256-atom MACE plus fixed-FF ALPB objective reach its original full-force convergence criteria when continued from the exact admitted step 172 checkpoint with a larger finite allowance?

This implements root `SCIENCE_DECISION.json` (`e8c4bda02f6714f452251b2c4edd83dbbdebfdbafa81b84cfdef75cf5574cc80`). The sole state is 7CCN–La, charge −2, physical multiplicity 1, original 256 atom IDs and original physical source. The checkpoint is `14cd220d924040caf046060c0b6fe3869941bff4bfeb9b210cc1a1e74a5eeacc`; its canonical JSON digest is `2d9a27465af93f883099a3e9d0562af1ffe98f41ad8d590587cb10dce3a06535`. The starting terminal POINT is `67690a33f7233a3e1b17b1abb687d41717fefd5903e3c04199a5060b1d6b7d81`.

## Admission before any new calls

The one run first calls unchanged completion adapter `551d964b` with its exact actual `CONFIG_RUN.json` (`8a78864e`). That saved loader admits the original and step 92 chain and creates one evaluator with no model initialized. Its original actual preflight and root run decision are re-admitted. The new adapter then admits actual 15458 (`9159ea43`), SEARCH (`4d59031b`), all 193 archived POINT records/189 cache keys, and all 172 prior requests. Every prior key is retained; saved chemistry is not requested again. The original initial point, eight history pairs, saved trust 0.1 Å, exact state/recipe and source-relative physical guard are preserved.

Actual saved audit 15489 (`228ef785`) is the completed audit gate: 81 points, 80 transitions, all 2,828 replay checks passed, zero molecular calls/model loads/optimizer steps/geometry trials. Its exact predecessor, checkpoint, source, recipe, settings and helper bindings must agree. This is admission of recorded audit evidence, not a repeat of its numerical analysis. No separate preflight allocation is required. Within the sole run, raw endpoint gradient/cache/work and physical guard replay still precede lazy native startup, using inherited exact producer arithmetic and gates.

## Allowed change and finite budget

`CARTESIAN_OPTIMIZER_GENEROUS.py` is copied from optimizer `c5092f24`. `OPTIMIZER_DIFF.txt` records only the exact checkpoint entry, old/new limit transition, digest and explanatory messages. Numerical direction, rigid-mode handling, history, trust, Armijo search, physical gates, native evaluator and all tolerances remain unchanged.

The allowance is at most **1,000 additional accepted steps and 2,000 distinct new evaluations**, with cumulative caps **1,172 accepted / 2,172 evaluations** from actual 172/172. Each new composite evaluation attempts at most one whole MACE energy/force and one native fixed-FF ALPB energy/gradient, so each component has a 2,000-call ceiling. One model is loaded lazily and kept warm. No DFT, GFN, MM, other metal, other source, or other search is allowed.

Stop using the original rules: maximum atom force ≤0.10 and atom-vector RMS force ≤0.03 kcal/mol/Å together; either finite allowance; native failure; state departure; or original no-admissible-decreasing-trial termination. Preserve unavailable outcomes and partial native files/checkpoints. No retry, budget extension, reset, parameter change, intermediate DFT, or automatic follow-up follows this attempt.

## Runtime and outputs

Root submits exactly one run on `gpu_h200`, node `node-224-2t-8gpu-1`, with 12 CPUs, one H200 and `standard_expand` defaults, without memory/time overrides. The exact original MACE virtual environment/Python and version/asset gates are inherited. OMP uses 12 threads; other BLAS thread controls remain 1. All numerical imports/replay and molecular work run under Slurm. Login validation is limited to source/JSON/hash/AST inspection.

CLI: `CONTINUE.py --config CONFIG.json --output <fresh run/execution>`. Root owns the launcher. READY records the exact command. The directory must not exist.

START saves successful zero-call chain admission before native startup. RESULT, SEARCH, all POINTs, raw requests and optimizer checkpoints preserve the cumulative history, actual new calls/costs, component work from step 172 and from the original start, full gradients, final force metrics and state/physical outcomes. The result denominator remains one declared state even on failure. Other-metal force availability is unchanged.

## Interpretation

This tests convergence of the declared model objective. Prior energy comparisons support the tested La descent locally but do not qualify force accuracy or establish a physical minimum, affinity or metal preference. Solvent-dominated lowering and opposing component forces remain visible. The inherited PBE0 references retain their vacuum-dianion, frozen-f, native-SCF and unqualified-gradient limitations. No scientific threshold or new accuracy gate is inferred from audit success.
