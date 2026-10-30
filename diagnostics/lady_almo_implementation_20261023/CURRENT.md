# ALMO implementation: current checkpoint

> **2026-10-24: superseded by `../lady_almo_ct_20261024/README.md`** — quartet CT done: ΔΔCT = −0.41 (raw) / −0.60 (counterpoise) kcal/mol, borderline vs 0.5 threshold; backend offset cancels to 0.01 kcal/mol.

New-session handoff requested24October2026:
`../../ALMO_SESSION_HANDOFF_2026-10-24.txt`. It records exact reusable assets,
terminal evidence, unresolved gates and the proposed replacement goal. The user
has not set that replacement goal; the handoff turn launches no scientific work.

The user goal is **implement the ALMO constraint as described**. Root owns all
submissions. Use standard_expand; preserve other campaigns and running jobs.
This goal authorizes the published AO-block implementation and its validation,
not fitting, new La/Dy geometries, production changes, commits or pushes.

## Completed evidence

| Job | Result |
| --- | --- |
| 22332 | Core solver tests pass: orbital finite differences, density/support/gauge invariants, one-fragment SCF and constraint-release recovery. |
| 22333 | PySCF PBE water benchmark passes its internal numerical tests. Independent native comparison remains unvalidated. |
| 22343 | Native400Ry benchmark failed at first water SCF iteration300 before any ALMO stage. |
| 22344 | Native600Ry benchmark reached its 3500-second cap during first-water SCF, at iteration246; no ALMO stage reached. |
| 22439 | Strict saved comparison rejected native return codes1/15 before comparing energies. No components admitted. |
| 22334 | Isolated La fragment converged; all native checks pass. |
| 22336 | Aquo ligand shell converged; all native checks pass. |
| 22411 | Isolated Dy fragment converged after explicit initial-guess recovery; all native checks pass. |
| 22337 | Mixed ligand-shell fragment converged; all native checks pass (gradient 4.38e-8). |
| 22429 | Fixed-density XC diagnosis passes; documented ORCA constants explain only 8.5034e-5 Eh of the La gap. |
| 22327 | La-aquo SCF converged but archive parity failed: +0.000825580291 Eh / 0.518059454 kcal/mol. |
| 22409 | Dy-aquo SCF converged but archive parity failed: +0.000956700399 Eh / 0.600338564 kcal/mol. |
| 22329 | La-mixed SCF converged but archive parity failed: +0.000893668327 Eh / 0.560785342 kcal/mol. |
| 22410 | Native mixed-Dy SCF reports convergence, but final gradient1.067060660e-7 fails strict1e-7 admission. |

The actual delivered 22328 and22330 failure events were read and acknowledged
in `baseline_aquo_Dy_v1/ROOT_ACK.json` and
`baseline_mixed_Dy_v1/ROOT_ACK.json` respectively. Other proactively inspected terminal results
have ROOT_REVIEW records; these are separate from acknowledging later delivery.

The failed Dy startup jobs 22328/22330/22335 stopped inside PySCF MINAO before
SCF. Explicit recoveries 22409/22410/22411 use native `1e` initialization and
preserve the original Hamiltonian, inputs, grid and every admission gate.
No failed attempt is overwritten.

Mixed-Dy22410 has now ended. Cycle28 satisfied the original gradient criterion,
but PySCF's final extra diagonalization raised the returned-state gradient to
1.0670606603907637e-7. Installed hf.py extra-cycle logic relaxes its internal
limits and accepts energy OR gradient convergence. The wrapper correctly
retains the original1e-7 criterion and rejects that final state; archive
comparison was not reached. Saved orbitals/checkpoint remain available.


## Finite benchmark round is terminal

The remaining native benchmark22344 and saved comparison22439 have ended.
There are no pending jobs in this declared validation round. Do not duplicate
either run or raise its limit automatically.

The600Ry native process ran3501.126seconds and returned15 after the declared
3500-second cap. Its last completed first-water SCF row is iteration246:
residual7e-8 versus1e-10, energy-76.5599956682Eh and change+8.38e-13Eh.
No second fragment or frozen/polarized/full ALMO states were reached. This is
a timeout, distinct from400Ry's300-iteration convergence failure. No terminal
integrated-charge summary exists for600Ry; the400Ry charge error must not be
assigned to it.

Collector22439 correctly retained native return codes1/15 and rejected the
comparison before parsing energies. It used zero new native calls. Its result
and collection provenance are in `cp2k_saved_comparison_v1/RESULT.json` and
`RESULT.json.collection.json`. The two benchmarks and collector consumed
48,175 allocated CPU-seconds (scheduler elapsed times2514,3507,7seconds).

22416 saved the first fixed-density energy/Fock/gradient, reproducing the La
baseline within 9.09e-13 Eh, then failed serializing NumPy integer metadata keys.
22419 fixed serialization and reused that first evaluation. Its explicit PBE0
component expansion then matched the native energy/Fock, but the custom-name
lookup failed because PySCF lowercases names in its dispersion parser while its
custom-functional registry is case sensitive. The lowercase-name recovery 22429 completed the final parameter variant,
reusing both earlier evaluations. The residual remains +0.000740546245 Eh,
above the unchanged archive gate; no baseline admission follows.
No functional constants, density or scientific gate change in either repair.

Frozen-density ORCA attempts 22432/22433 failed before computation: eight
threads in one Slurm task supplied only one MPI slot. Exact inputs/GBW remain
preserved. A fresh metadata bundle binds the failures for the same two native
checks with eight MPI tasks, one CPU each, 32 GiB, one hour and native cap900s.
Recovery jobs22436/22437 both ended normally; 22436 passes. Direct22437
failed only the parser predicate banning any startup auxiliary V-matrix setup.
Both independent source audits find this predicate incorrect: actual native
property timings record J for direct and SPLITRIJ/COSX for the control.
Collector22438 preserves the original refusal. Saved-only recovery22441
passes, requiring native J versus SPLITRIJ/COSX timing evidence and every
other original gate. Direct minus RIJCOSX = +0.000735320167 Eh
(+0.461420371 kcal/mol); replay minus archive = -0.000045388230 Eh,
inside the unchanged replay gate. No native rerun. The collector dependency
correction is saved separately.
No successful energy evaluation is repeated. See ORCA_FROZEN_SPEC.md.

All submissions have immutable source snapshots, CONFIG/SUBMISSION records,
Slurm logs and actual completion watchers under
`workspaces/lady_almo_implementation_20261023/`.

## Unresolved validation

The fixed per-endpoint archive gate remains **1e-4 Eh**. No real-system ALMO
quartet has been submitted or admitted. Native400Ry failed before any ALMO
stage: first monomer residual1e-8 exceeded1e-10 at300iterations, despite stable
energy. Integrated net charge was -0.0306942524, a separate numerical
representation/integration discrepancy with no established unique cause.
No component agreement is available:600Ry timed out in the same first-monomer
phase, and the strict collection rejected the failed native runs. This does not
demonstrate disagreement between two ALMO implementations: the native ALMO
states were never obtained. No retry, ignored-convergence option or tolerance
change is released.
Read CROSS_CODE_AUDIT.md and PLAN.md for inspected numerical differences and
finite diagnostic scope. No basis/ECP conversion defect has been established.

The core solver remains unchanged at SHA256
`ccaa8a5041373691796958ecb81b41965021cd46cf8c640083e6c39028d4680e`.
Verify actual test source pins before reuse.

## Next actions

1. Read existing actual terminal receipts and results. Preserve and reuse them.
2. The declared fixed-density diagnosis is complete. RIJCOSX/direct J/K and
   documented XC constants explain much of the observed La discrepancy.
   Preserve density, final-COSX-phase and quadrature distinctions; this is not
   an exact additive cross-code closure or automatic baseline admission.
3. The native benchmarks and strict saved comparison are now terminal. Preserve
   their failures; repeated unchanged calculations or longer time limits are not
   supported by these results. Independent native validation remains unresolved.
4. Resolve implementation/native/archive qualification before any scientific CT
   double difference. Keep software correctness distinct from physical accuracy.
5. If admitted, reuse ordinary/fragment states for the four concurrent ALMO runs,
   then use collect_quartet.py. Final reporting must preserve any unresolved limit.

The controller still records an earlier capacity-only blocked status. Capacity
has opened and this completion event resumed authorized work; that old queue
condition is not the present scientific blocker. The API does not expose an
active/resume status mutation. Do not replace or mark the goal complete.
Previous checkpoint text is preserved in immutable source snapshots, including
`xc_fixed_density_diagnosis_v2/source/CURRENT.md`.
