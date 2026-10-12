# Persistent CPU OMOL proposals: finite scout passed

12 October 2026. **The unchanged four-angle proposal method works on CPU for the
tested202-atom Ca-conditioned fixture.** Five analytic-force replays and both
full SLSQP starts passed. Final energies exactly match archived GPU values;
maximum coordinate differences are below2×10⁻¹⁴Å. Actual cost was33 model
forwards,443.418s elapsed,2766.676 processCPU-s, zero failures and zeroGPU.
No107-case campaign, solvent/DFT calculation or production edit occurred.

## Actual qualification

Source: `a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1`, existing
4.3Å envelope `e05e200a13d024a5`,202 atoms, chargesCa−3/La−2, singlets. The
scout was declared in [SCOUT_PLAN.md](SCOUT_PLAN.md) before any model call.

| Check | Actual result |
|---|---|
| Ca/La origins and saved final geometries, four fresh E/F calls | All energies exactly match GPU; maximum Cartesian force error1.470×10⁻¹³eV/Å |
| Repeat Ca origin after switching states/geometries | Energy identical; CPU repeat force error8.88×10⁻¹⁶eV/Å |
| Paired-force mode selection | Same four IDs, order, criteria and redundancy decisions |
| Ca optimization |14iterations,16function evaluations,15new force calls; final energy identical, coordinate difference1.948×10⁻¹⁴Å; same0.8Å boundary |
| La optimization |14iterations,14function evaluations,13new force calls; final energy identical, coordinate difference1.518×10⁻¹⁴Å; same interior endpoint |

Frozen replay gates were the existing `mace_omol.py` energy0.01kcal/mol and
force0.001eV/Å tolerances. Final proposal energy used the existing
`slsqp_precision.py`0.001kcal/mol gate, plus unchanged optimizer success,
energy nonincrease, physical geometry and boundary admission. Iteration trajectory
identity was not required; both iteration counts nevertheless match the archive.
Actual residual gradients and full-q differences remain in the raw checks.

[SCOUT_RESULT.json](SCOUT_RESULT.json) is the completed result. Its collector
independently verified all33 request/result/force/XYZ pins, actual charge/spin/
head/atom counts, one forward per request, unchanged parameter versions, and
reconstructed all five raw force comparisons. Earlier [RESULT.json](RESULT.json)
is the file-only feasibility inventory, not an execution result. The final
declaration used **no added call/time ceiling**, no retry and no fallback.

## Smallest safe execution path

[cpu_scout.py](cpu_scout.py) adds a persistent CPU worker with the existing
`evaluate(request_pin) -> (result, result_pin)` interface. The parent directly
imports the historical source snapshot and calls its unchanged
`CompletionProposal`/`adaptive_completion.optimize`, selected through the private
`slsqp_precision.engine()` settings. It preserves actual archived q0 reuse and
evaluates every new trial on the CPU worker. Fresh cache keys include the new
manifest, so GPU/historical outputs cannot be overwritten or misidentified.

Keep the two environments separate:

| Role | Exact invocation | Relevant identity |
|---|---|---|
| Optimizer parent | `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python` |NumPy1.26.4, SciPy1.17.1; original SLSQP kernel `4e71d8a9…` |
| Warm CPU model child | `workspaces/mace_hybrid_20260916/software_v1/venv/bin/python` |Torch2.8.0, MACE0.3.16, NumPy2.4.6, ASE3.29.0, e3nn0.4.4 |

The model environment also reports SciPy1.17.1 but its compiled SLSQP kernel is
`ff528471…`, different from the pinned parent binary. A single-process move into
that environment would change the optimizer. No installation is needed.

The worker loads the exact existing checkpoint once (SHA256
`9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`), `omol` head,
float64, nonperiodic native vacuum model, and actual charge/multiplicity inputs.
Weights have gradients disabled and mutation counters remain unchanged; Cartesian
autograd stays enabled. Each request clears only the ASE result cache, retaining
the model; energy and forces within that request share one actual forward.
The energy-only `context_factorial_mace.py` direct no-grad call is not reused for
force computation. No CUDA/srun/Slurm gate is faked or disabled in shared code.

The actual scout ran this parent command (historical execution record; the
completed directory must not be rerun):

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -u \
  diagnostics/orca_parameter_audit_20261012/cpu_proposal_feasibility/cpu_scout.py \
  execute --manifest workspaces/orca_parameter_audit_20261012/cpu_proposal_scout_v1/manifest.json
```

The child command and owned process group are recorded in
`cpu_proposal_scout_v1/CPU_WORKER_START.json`. Parent numerical work uses one
thread while the child waits; child Torch work uses eight while the parent
blocks. Interop and other BLAS thread counts are one. This serial handshake
uses the assigned eight active numerical threads. Both processes have completed
and the worker was reaped; no allocation remains owned by this scout.

## Method preserved

- Common four angular modes from the existing alternating differential/individual
  heavy-normalized two-pass MGS selector, rank tolerance1e-5. No outcome labels
  enter the selector.
- One q=0 start per metal; analytic physical and cap Jacobians; unchanged
  `Kinematics.evaluate`, `project`, `preview`, `MACEProposal.evaluate` and geometry
  checks. Source files match their archived snapshot bytes.
- SLSQP, shifted Hartree-equivalent objective (`eV×0.03674932217934773`), gradient
  converted consistently, ftol1e-8, maxiter200, angular bounds±0.8rad,
  heavy-displacement bound0.8Å. The old `ftol_eV` metadata1e-9 is inherited but
  is not the active precision-stop tolerance.
- Final physical tolerance1e-7Å, angular tolerance1e-12rad, energy increase limit
  1e-7eV, physical-boundary margin1e-5Å, angular-boundary margin0.02rad. Infeasible
  intermediate SLSQP trials remain recorded; final admission is unchanged.

## Measured cost and implication for the missing preparations

| Quantity | Actual |
|---|---:|
| New model forwards |33 =5replays +15Ca +13La |
| Failed calls |0 |
| Full elapsed time |443.418s |
| Worker processCPU time |2746.393s |
| Parent processCPU time |20.283s |
| Force-call time, min / median / max |10.922 /12.420 /12.948s |
| Sum of force-call elapsed time |407.995s |
| Model construction inside warm worker |2.083s; imports/startup also included in full elapsed time |
| Peak worker / parent RSS |6.926 /0.275GiB; separate process maxima |

For context, the archived68 GPU proposals used998 fresh trial forwards, range7–20
per start, median15. Those are historical observations, not predictions for new
inputs. No CPU timing was inferred from the prior energy-only run, from GPU
speed, or from an assumed8→16-thread speedup.

The next contained implementation can use this exact warm transport and optimizer
parent on Curie's frozen preparations. New inputs require actual Ca/La origin
energy/forces first, then the existing `motion_envelope_run.build_tasks` paired
projection/selector and two proposal starts. Preserve native-only origin status;
do not invent missing solvent components to satisfy an old preparation API.
For107 genuinely new prepared sources that means214 origin force calls and at
most214 declared optimizer starts before any downstream scoring, with trial
counts unknown until execution. Existing compatible origins should be reused.

This scout qualifies the tested CPU execution path and preserves one boundary
and one interior proposal. It does not establish throughput for all context
sizes, scientific improvement, force-field accuracy, or permission to silently
replace historical scores. Root owns the expanded-coverage decision. The launcher
remains a finite scout rather than a new production default.
