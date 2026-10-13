# Hans8FNR La continuation succeeded — 12 October 2026

**The previously missing La endpoint is now available.** Job 2665 converged
and terminated normally; the existing collector accepted its electronic energy,
198-atom gradient and point-charge gradient. Five of the six declared repaired
origins are now available. Dy 2666 is still running; comparison 2669 already waits
for both jobs and must not be duplicated.

This completes one interrupted-orbital continuation under the unchanged
PBE0-D4 embedded recipe. It does not resolve the Hans-source reversal by itself,
qualify the complete hybrid forces, establish an electronic global minimum,
or yield a binding affinity/classification. The fixed-geometry endpoint has
nonzero nuclear gradients; SCF convergence is not structural relaxation.

## Actual result and state

- Source: consumed Hans8FNR, La_A; energy **−5804.400869806574 Hartree**.
- Output reports **SCF CONVERGED AFTER 52 CYCLES** and normal termination.
- LaIII hypothesis, charge −1, physical/effective multiplicity 1; ECP 46,
  816 explicit electrons (408 alpha / 408 beta), matching the declared model.
- Native electronic gradient and point-charge gradient are present. Collector
  records `status=complete`; force/physical-state qualification flags remain
  false, and affinity/classification remain null.
- The pinned interrupted same-state orbitals initialize a new solve;
  `target_energy_reused=false`. Geometry, state, Hamiltonian and admission
  tolerances were not changed. The original interrupted run remains preserved.
- This is PBE0-D4 DFT, independent of the native-GFN2 La shell-parameter repair.

[The pinned result](LA_RESULT_2665.json) preserves collection/output/gradient,
seed, execution and manifest identities. All four output/gradient/receipt hashes
were verified as metadata. The actual archived parser qualification remains
job 2664's eight passed tests; no new scientific test or molecular call ran
while writing this report.

## Independent output review

A second reader confirmed convergence directly in the raw output: final orbital
residual 1.1421e-7 against 2e-6, and energy change −2.8194e-11 Eh against 1e-9.
The final energy matches the native gradient header. The 198-QM and 2008-PC
count headers match the inputs. MORead reports matching geometry/basis and
Aufbau occupation reassignment; input, output, seed and gradient hashes agree.
Small-gap warnings remain, and the final RIJCOSX integration changes the reported
SCF energy by +0.000052490 Eh. No stability/ground-state or complete hybrid-force
qualification is inferred from successful collection.

## Cost and next operation

Scheduler: **32,640 seconds (9h04m), 12 CPUs, 391,680 allocated CPU-seconds**, zero GPU.
Executor interval: 32,631.268078360008 seconds and 391,575.2169403201 allocated
CPU-seconds; the scheduler interval also includes setup/collection and rounding. This is a small reference
calculation under the present allocation, not demonstrated screening throughput.
The shared grant was 48,000 MiB, with 12 MPI ranks and 3,100 MB MaxCore per rank.
These are allocations/limits, not peak utilization or measured peak memory.

Dy 2666 is still RUNNING. Its last read health record showed improving residuals
but no converged output; no Dy energy is accepted. Its terminal/health watchers
2088244/2089708 and the existing comparison watcher 2157957 were verified alive.
No job was cancelled or submitted. Read the eventual Dy collection and then
`workspaces/claude_lanm_20261012/comparison_root_recovery_v1/RECEIPT.json` and
`RESULT.json` after 2669 completes. Both Hans sources must remain separate.
