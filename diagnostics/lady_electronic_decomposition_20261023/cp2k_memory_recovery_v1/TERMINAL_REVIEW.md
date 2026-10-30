# Independent terminal review: La-aquo pilot 21623

**Conclusion: the predeclared native wall-time cap stopped an unfinished first-fragment SCF. No fragment or ALMO scientific result is admitted.** This termination was not the earlier OOM and is not a native electronic-method incompatibility determination. The bounded pilot did not meet its admission requirements; the other three quartet cells remain unreleased. No restart, threshold change or further calculation follows from this review.

## Exact stop mechanism

`RUN_RESULT.json` records `status: failed`, one molecular-parent call, `timed_out: true`, a 7,073-second native cap, and `RuntimeError('Predeclared process cap reached: 7073s')`. The native launcher returned 15 after 7,078.734081486007 seconds. Its log explicitly reports `Terminated (signal 15)`.

The reviewed runner explains this chain exactly. `RUN.py:41–52` catches `subprocess.TimeoutExpired`, marks the command timed out, sends SIGTERM to its own process group, waits up to 15 seconds for cleanup, then raises the recorded error. The extra approximately 5.7 seconds in the native elapsed receipt are termination/cleanup time beyond the timeout trigger, not an authorized SCF continuation. The outer handler saves the failure and returns process exit 1. The saved completion event reports Slurm FAILED, exit 1:0, elapsed 7,087 seconds, 24 CPUs. There is no OOM termination in this recovery's supplied receipts and no native `SCF run NOT converged`, `Leaving inner SCF loop`, or electronic-method abort in the inspected output. A timeout before those native endings does not demonstrate that the SCF would ultimately converge.

`RESULT.json` is **absent**, not passing. This is explained by source control flow: the timeout exception at the native call (`RUN.py:130–131`) skips the subsequent parser invocation (`134–138`) and reaches the failure/finalization handlers. There are only the syntax and native commands in RUN_RESULT. The saved raw output and failed RUN_RESULT remain the terminal authority; no parser was rerun in this review.

## Work actually reached

The same frozen input hash remains `15b07a8ba0a5cb8037847908dd87c6977fa50ed53c8a07ca7489fc7a5d92883b`. Syntax checking exited zero. The runtime grant is 262,144 MiB and the recorded launch remains 12 MPI × 2 OpenMP.

The output has one molecular-guess fragment banner: fragment 1 of 2, La charge +3 and multiplicity 1. Its native readback is eight electrons, four occupied orbitals, and 47 orbital functions. The neutral-atom auxiliary guess was followed by charge rescaling and the actual isolated-ion OT calculation. That ion calculation printed 43 numbered `OT CG` iterations.

The last numbered row, native output line 1673, is:

```text
43 OT CG  0.34E+00  103.9  0.00000040  -30.1574842906  1.46E-12
```

Thus the last printed SCF convergence measure is **4.0e-7**, above the configured **EPS_SCF = 1e-10**. A following line-search energy is also present. The nearly unchanged displayed energy is not evidence that the required density/orbital convergence threshold passed, and the last energy is not an admitted isolated-ion reference.

There are zero native fragment-converged markers; fragment 2 never appears. There is no native fragment sum R, frozen initial energy F, converged block P, fully delocalized U, ALMO optimization row, final FORCE_EVAL energy, or normal program end. No requested cell or constrained electronic state is complete. Partial restart files exist in the execution directory but carry no scientific admission or permission to resume.

## Outstanding requirements and interpretation

Both fragment SCFs, the frozen assembly, block polarization and fully delocalized SCF remain uncompleted/unadmitted. Consequently none of the full-energy archive-agreement, fragment-sum, decomposition-closure, variational-order or quartet-contrast requirements can be evaluated from this attempt. CT's magnitude and its relation to the model discrepancy are **unavailable**. This is a failure to complete the bounded computational diagnostic, not evidence for or against substantial CT or a specific force-field defect.

The worker receipt records 7,086.493903486058 seconds and 170,075.8536836654 allocated CPU-seconds; these are wall-based capacity figures, not measured CPU utilization. Root owns final scheduler/peak-memory and cumulative cost accounting, including failed attempt 21618. All failures and raw files remain preserved.

## Reviewed pins

Task base: `workspaces/lady_electronic_decomposition_20261023/cp2k_quartet_memory_v1/tasks/aquo8_central_La/`.

| Artifact | SHA256 |
|---|---|
| `RUN_RESULT.json` | `536e8fcd4bb727394b3f4df203990b70c349394ea889d22b55a0fe769c67548a` |
| `completion_event.json` as read | `536602e587f8f6df2dd8fba119b55658bedd9abfddc4e85b23e0e9ffbf54f1c8` |
| `execution/endpoint.out` | `329c094e3a97a00c298cafb72eae4935a04267b0183c6c1f3cb9bc52978e84ea` |
| `execution/native_launcher.log` | `51320e4584c86023312d3bf6b96e16929f89e2c43987ede298d3dcceb4a420ab` |
| Diagnostic `cp2k_memory_recovery_v1/RUN.py` | `9244ba6036bf220760e27bca95f1e4eb0e63e18e55e039600c84f3317961d7f8` |

Scope: saved receipts/raw text and source inspection, with metadata hash/record counts only. No scientific numerical analysis, parser execution, native call, submission, cancellation or source/input change.
