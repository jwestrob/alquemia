# Dy diagnosis: separate numerical stationarity from the physical model

28 September 2026. Five saved failed outputs audited; no new molecular calls.

## What the actual outputs establish

All five molecular outputs lack a qualified converged Dy state. The embedded
195-atom cases end with orbital residuals approximately 187 and 91. The isolated
50-atom native starts instead approach residuals 1.26e-4 and 4.77e-4; isolated
PBE0 reaches 7.22e-5, still above the printed 1e-5 tolerance. The isolated cases
therefore are a different numerical regime from the gross embedded divergence.
This is not evidence that either isolated energy is acceptable.

The smallest overlap eigenvalues are 1.54–1.55e-5 in the embedded systems and
1.16–1.47e-4 in the isolated systems, above the printed 1e-7 cutoff. There is no
observed removed-basis-vector explanation. PModel atomic initialization succeeds;
its neutral atomic state is not the molecular Dy(III) state. No qualified local
spin/4f occupation assignment is available. Do not identify ligand reduction,
incorrect oxidation state or wrong 4f filling solely from these unfinished SCFs.

The two isolated native starts approach different intermediate energies. Their
unconverged difference is not a state splitting, but cautions against declaring a
unique solution from one apparent energy plateau. Snapshot GBWs can be retained
as explicitly unconverged diagnostic seeds, never as validated wavefunctions.

## Two separate questions for the next model decision

1. **Numerical:** is the near-root residual limited by integration/integral
   precision or a difficult electronic surface? ORCA documents that numerical
   errors exceeding requested convergence accuracy prevent convergence. PBE0's
   COSX approximation adds a possible numerical source, but native r2SCAN also
   plateaus, so COSX alone cannot explain all failures. No such causal test has
   yet run. A targeted precision comparison would hold the Hamiltonian, geometry,
   spin and basis fixed and compare actual orbital residuals; not another
   functional sweep, tolerance relaxation, or new biological score.
2. **Representational:** can a declared Dy(III) 4f-in-core model capture the
   relative environmental and mechanical response needed for discrimination?
   Existing ECP55 calculations can establish executability, not the missing
   explicit-f response agreement. Their energies must not be mixed with ECP28
   offsets or La endpoints from a different method.

Primary numerical reference, checked 28 September:
https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html
Sections 2.6.1 and TRAH discuss numerical precision and convergence. This supports
numerical diagnosis as a hypothesis, not a diagnosis of these particular outputs.

## Next finite experiment design, pending representation audit

Reuse the actual small Hans EF3 coordinates and separately mapped source-derived
perturbation. First require one frozen-f analytic-gradient endpoint with the exact
published/pinned ECP and matching orbital/auxiliary basis, explicit electron/core
bookkeeping, physical Dy(III) sextet metadata and distinct effective valence spin.
Do not inherit the old large-core SP results as gradient validation.

Then qualify a matched A/B response and its directional force, keeping geometry,
charge, proton/water inventory and environment fixed within comparisons. Use
same-metal response differences before any cross-metal score. La/Dy selectivity
requires compatible matched La calculations and a balanced cycle; no raw totals,
PQQ thresholds, or sign-selected source. Preserve both Hans source contexts before
claiming within-series transfer. This remains a research candidate with omitted
4f-response/spin-orbit physics, not an established reference oracle.

No new molecular submission has been made by this audit. The frozen-f assessment
must determine actual basis availability and missing force/transfer checks first.

## Reproduce this saved-output audit

From repository root:

```bash
python diagnostics/metal_environment_response_20260926/electronic_diagnosis_20260928/audit_saved.py --output /tmp/nikasha_dy_saved_audit.json
```

Output must not already exist. EXPLICIT_F_EVIDENCE.json pins all five source
outputs and retains macroiterations, including NR steps. Failed iterate energies
are explicitly diagnostic only. No synthetic molecular result is used.
