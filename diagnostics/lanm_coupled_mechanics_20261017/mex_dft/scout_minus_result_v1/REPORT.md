# Mex derivative scout: minus endpoint complete

17 October 2026 UTC. Job 12978, task `MM_O97_La_minus`, completed normally. The unchanged endpoint adapter reports `complete` and DFT `available`. This is one displaced singlepoint; the paired derivative test still awaits plus 12977.

## Actual execution

- Native SCF converged after 32 cycles under the existing energy-based check. ORCA terminated normally with return code 0.
- Final gridx3 exchange recomputation took 1,920.673 s. Native total runtime was 6 h 08 m 43.651 s; scheduler elapsed 22,132 s on 12 CPUs, or 265,584 allocated CPU-seconds.
- Executed work: one DFT singlepoint, one MM configuration/four group queries, and one fixed-charge ALPB evaluation. No central gradient, MACE, GFN, optimization or new geometry search.
- QM state remains charge 0, 832 explicit electrons, 416 alpha/416 beta, La ECP46, effective multiplicity 1. The frozen-f physical-state qualification remains false.

## Admission limits

The native energy residual passes its printed criterion. Printed maximum/RMS density change, orbital gradient, orbital rotation and last DIIS error do not all pass their displayed tolerances. The original admission rule explicitly retains these observations rather than imposing a new all-criteria gate. Energy convergence is not electronic stability or force accuracy. The original pair finite-difference test remains required, followed by its measured-cost review before the remaining six qualification points. The eight-cell response panel still requires the full qualification.

## Quiet postprocessing review

The last output was `Recomputing exchange energy using gridx3` while all 12 owned MPI ranks remained active. CPU activity and unchanged files alone did not prove useful progress. Independent review by Arendt (`mex_health_review`) found compatible archived phase timings of 255.757 s/227.934 s on 112 MPI ranks, compared with the present 12. The 900 s health alert required review, with no declared post-SCF silence cancellation rule. Actual normal termination subsequently resolved the uncertainty; no restart, cancellation or changed setting was applied.

## Preserved evidence and next action

This directory contains byte-identical native output/receipt, inputs, runtime configuration, endpoint and MM/ALPB results, submission, scheduler accounting and a source/hash inventory in `ROOT_REVIEW.json`. The workspace remains the authoritative execution source. Large scratch files and wavefunctions remain there.

Await the actual plus endpoint, bind both terminal results into the original `scout_pair` collector, and run that saved comparison on 1 CPU through Slurm. No paired derivative, metal preference, mechanics accuracy, or production change is claimed here. Whole-peptide DFT 13051 continues separately.
