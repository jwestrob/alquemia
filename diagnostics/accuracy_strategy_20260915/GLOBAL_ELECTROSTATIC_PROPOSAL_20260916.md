# Proposed next accuracy experiment: one global electrostatic model

**2026-09-16 — PROPOSED, not approved or executed.** Jacob asked whether there
is a coherent way forward. This note defines the recommendation. It advances
the direction in STRATEGY.md by identifying an available solver and an explicit
energy expression. It does not change the baseline or authorize calculations.

## Question and model

Does including the actual protein's permanent electrostatics and solvent
boundary improve the existing La/Ca ordering conflict, at practical cost?

Use the existing small quantum cores, but represent the whole recorded protein
chain/assembly in the classical electrostatics and dielectric boundary. For
each metal M, evaluate the frozen-distribution descriptor

    Etilde_M = E_QM,vac,M + C(q_M,Q) + G_RF,D(q_M+Q).
    R_global = Etilde_Ca - Etilde_La.

- E_QM,vac is a new native ORCA r2SCAN-3c gas-phase endpoint, including its
  ordinary native composite-method corrections. No CPCM is present.
- q_M is that endpoint's MBIS charge distribution, with the actual charge sum
  and native La ECP convention checked. No formal-metal-charge substitute.
- Q is one fixed protein charge distribution after the recorded QM boundary
  exclusions. C is the complete core/environment Coulomb interaction.
- D is one physical protein/solvent surface shared by the two endpoints and
  by the two representations in the partition test.
- G_RF is **only** reaction-field/solvation energy, not total PB electrostatic
  energy. In the linear model it is one half of total charge times reaction
  potential. Intracore Coulomb energy is already in QM and is not added again.

The reaction-field expression includes core self response, core/environment
response coupling, and environment self response. The environment-only term
cancels in Ca-minus-La only because Q and D are identical. Other omitted
environment-only bonded/nonbonded terms cancel only when endpoint-independent;
metal-dependent dispersion/repulsion across the QM boundary is not claimed to
be included. This is an electrostatic descriptor, not a complete hybrid binding
free energy or self-consistent QM/PB calculation.

## Difference from the failed transfer model

The old model added a point-charge PB transfer correction to a quantum CPCM
energy, subtracting an isolated-core point-charge PB reference. That reference
changed strongly with the numerical partition and did not exactly represent
the quantum CPCM energy.

The new expression uses vacuum QM plus one global solvent model. It contains
neither CPCM nor the old isolated-core PB counterterm. This removes that
particular mismatch. It does **not** predict success: gas-phase densities also
change, and caps/charge ownership can still make the partition inconsistent.

## Concrete backend and fixed physical choices

Use the archived APBS 3.4.1 distribution's TABI-PB surface solver and NanoShaper.
Both binaries are present under
`workspaces/affordable_challenger_20260915/software/apbs-3.4.1/APBS-3.4.1.Linux/bin/`.
Their standalone files currently lack executable mode; any implementation
should make pinned executable copies in its own workspace, leaving the archive
untouched. Availability and dependencies were inspected; runtime is untested.

TABI solves on the molecular surface rather than a volume grid. This changes
the numerical approach, not the requirement to demonstrate convergence.
Its separate Coulomb output must not be added to the formula above.
[Official solver and energy documentation](https://apbs.readthedocs.io/en/latest/using/input/old/elec/tabi.html).

Retain the earlier environment model's physical settings: dielectric 1 inside,
78.54 outside; 298.15 K; zero salt; 1.4-A solvent probe; the recorded Bondi
H/C/N/O/S radii and common 1.80-A Ca/La radius. Dielectric 1 is consistent with
vacuum QM and unscreened direct Coulomb. No dielectric/radius fitting.
Freeze and record the surface-resolution and solver-convergence sequence
before scientific execution, using the installed interface's actual units.

Reuse the exact archived source coordinates, core membership, protonation,
waters, charges/multiplicities and assembly selections for each named case.
Record full source-to-QM/cap ownership; no force-field charge on a represented
QM source atom and no overlapping charged cap/environment nucleus. Reuse the
audited local exclusion/closure scheme as a declared approximation, judged by
the partition test, not by charge closure alone. The physical cavity uses
source atoms, not cap spheres. Unsupported nearby chemistry remains unscorable.

## Proposed experiment and run count

1. **Physical feasibility:** the already-consumed 1H4I qm33/qm36 partition pair:
   two representations times La/Ca = **four new QM endpoints**. Same actual
   protein and chemistry. Use their charges for surface refinement, rigid
   transformations, component accounting, and isolated-environment reduction
   checks. Reuse the prior numerical acceptance scales: at most 0.5 kcal/mol
   change in the paired contrast under refinement/rigid transforms, and at
   most 2 kcal/mol between the two numerical partitions. Reduction/algebra
   consistency tolerance 0.01 kcal/mol. No new QM calculation is needed for
   these solver checks. Record the exact solver-only task manifest beforehand.
2. **Accuracy trial, conditional on physical feasibility:** GGR 1GLG; alpha
   1F6S and 6IP9; canonical PQQ 1H4I and 4MAE: **five pairs, ten new QM
   endpoints**. Use the repaired generic v3 cores and exact canonical fixed
   PQQ cores. These are consumed development cases, not blind tests.

Total proposed main sequence: **14 new local QM endpoints**, with ordinary
technical retries recorded separately. This specifies the scientific scope;
it is not a CPU-time or wall-time budget. No numerical model change is an
automatic rescue branch. The five-case full-environment preparation must be
reviewable before their jobs run; unavailable cases remain explicit rather
than receiving guessed chemistry.

The performance readout is whether both alpha preparations move above GGR on
R_global while the La-associated PQQ control remains above the Ca-associated
one. Report the continuous margins and changes from archived baseline ordering,
including worsening or partial improvement. Two alpha structures are one
biological observation. These distinct evidence strata are not pooled accuracy.

There is no compatible new aquo reference or calibrated threshold yet; neither
the old PQQ bands nor a universal zero applies. Relative differences between
sites cancel a common aquo offset and suffice for this first accuracy question.
No new reference calculation or threshold fit is part of this proposal.

## Expected utility and limits

The practical target is two local QM endpoints plus classical surface solves
per site, with no full-protein quantum orbitals. The actual solver cost has
not been measured. Record pair cost and solver cost separately on suitable
hardware; no large allocation should sit idle around a serial solver.

This tests a concrete candidate improvement, rather than producing derivatives
without a correction model. It is motivated by omitted protein interactions,
not by a proven diagnosis of alpha or GGR. The main scientific limitation is
the frozen gas-phase density, especially for anionic cores: it does not respond
self-consistently to the protein/reaction field. Relaxation, hydration-state
changes and entropy remain outside this one model.

If physical consistency fails, the proposed energy model has not earned the
accuracy trial. If consistency passes but ordering does not improve, report
that static electrostatic model's failure rather than starting a parameter
sweep. If it improves, the curated mutant comparisons are the next independent
development opportunity; this small consumed panel alone cannot justify
promotion. Baseline/default remains intact in every outcome.
