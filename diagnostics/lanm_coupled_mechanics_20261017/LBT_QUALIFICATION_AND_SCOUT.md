# LBT broad-coordinate qualification and scout

17 October 2026 UTC. **The full-coordinate interface passes its numerical checks, and the first scout lowers previously unresolved strain. It is not converged or independently validated.**

The exact prior 7CCO/La restricted terminal retained a maximum full-atom force of 83.2803 kcal/mol/Å. Twelve accepted full-coordinate steps lowered the unchanged MACE-vacuum + fixed-charge ALPB energy by 63.9031 kcal/mol (MACE −63.5420; solvent −0.3612). Maximum force fell to 12.7867 and RMS atom-vector force to 3.43465 kcal/mol/Å. All topology, chirality, peptide and contact gates passed. The target stayed bound; maximum motion from the source was 0.29917 Å. This scout stopped at its declared 12-step limit, far above the 0.10 maximum/0.03 RMS force convergence requirements.

## Qualification and reuse

The original zero-call preflight12954 failed exact comparison of an old coordinate transformation. Diagnostic12956 isolated only 1.11–1.94e−16 representation differences; every other lineage condition passed. V2 uses the pinned archived transform as authoritative after bounded representation/physical-basis checks and the original whitening gate. It changes no calculator, geometry or scientific admission tolerance. Both the original failure and diagnostic remain preserved.

Preflight12957 passes all four source/metal starts. Qualification12958 passes all eight directional cases: coupled whole-peptide motion and a real peptide-bond stretch, for both metals and sources, at both frozen finite-difference steps. All MACE, ALPB and total derivatives pass. The 32 new composite evaluations used one warm model. Numerical derivative consistency alone does not establish accuracy.

Saved baseline12959 reproduces all 44 cells in the 22-geometry, source-separated old pools and rejects a deliberately omitted required cell. It makes zero molecular calls.

## Measured cost and next decision

| Job | Actual work | Scheduler time | Allocation |
| --- | --- | --- | --- |
| 12957 | Four-state zero-call replay | 17 s | 1 CPU |
| 12958 | 32 MACE + 32 fixed-charge ALPB evaluations | 30 s | 12 CPUs, 1 GPU |
| 12959 | Saved pool replay and missing-cell test | 16 s | 1 CPU |
| 12960 | 12 MACE + 12 fixed-charge ALPB evaluations | 23 s | 12 CPUs, 1 GPU |

Application time for the scout was 21.813 s, including a 6.543 s warm-model startup. There were no new GFN2, MM or DFT evaluations. Scheduler allocations are costs, not measured utilization.

The scout admission/cost gate enables the declared six-search comparison: four LBT states and Mex La/Dy. Each is limited to 80 accepted steps or 200 distinct evaluations, including the already completed LBT scout. Both metals will cross-score the same compatible source-local pools. The next scientific question is whether broader response changes the metal-dependent work and source consistency; extra energy lowering alone is insufficient.

Exact source/config/result/receipt hashes and submitted continuation IDs are in INITIAL_EXECUTION_RECORD.json and workspaces/lanm_coupled_mechanics_20261017/SIX_SEARCH_SUBMISSIONS_v1.json.
