# Autonomous improvement work — 2026-09-16

Jacob explicitly authorized continued contained scientific experiments after
the global pilot's partition and rotation failures were reported:

> i authorize you to continue working as much as you want- no pressure of
> course- and to try out ideas (in contained manners; preserve the baseline
> for other agents to use) to improve the discriminator.

He then confirmed:

> full permissions. proceed as you like, do whatever you need, as long as you
> like. just as long as you follow the rules you already know.

This grants discretion to choose and execute contained improvement experiments;
it does not authorize changing production defaults, historical experiments,
labels or thresholds. The agent still records each concrete scope and energy
expression before execution. Costs remain measured without CPU/time stopping
budgets. No model is promoted automatically. Existing shared execution and
directory policies, real scientific inputs, and explicit failures still apply.

The frozen global v1 campaign will finish unchanged and retain its failed gate.
Its conditional accuracy stage is not reinterpreted as approved by that gate.
Follow-on work below is a separate development diagnostic under this new grant.

## First concrete experiment

Question: does approximating the quantum density by atomic monopoles, or
freezing that density in vacuum, materially cause the observed partition
dependence? These are separate hypotheses, not established diagnoses.

Use exactly the four already consumed 1H4I qm33/qm36 × La/Ca states from
`workspaces/global_electrostatic_20260916/states_v1/`. Keep every coordinate,
charge/multiplicity, QM atom, cap, permanent environment charge and exclusion
unchanged. Native ORCA 6.1.1 r2SCAN-3c, DefGrid3, same basis/ECP and composite
corrections. No geometry, water, protonation, assembly or dielectric search.

1. Four `orca_vpot` evaluations of the **saved vacuum density** at every actual
   environment point-charge coordinate. Calculate the exact-density Coulomb
   coupling `C_density = sum_j Q_j * phi_QM(r_j)` in atomic units, then convert
   once. Compare with MBIS monopoles in the same atomic-unit convention and
   retain the small separate TABI constant-convention difference. Retain
   source-atom/residue contributions to locate representation errors.
2. Four new native vacuum single points with the same full fixed environment
   inserted through `%pointcharges`, with external charge self-energy disabled
   (ORCA's documented default, explicitly recorded). Request MBIS and preserve
   the wavefunctions for possible later use. These are isolated copies and
   use the existing manifested ORCA runner. Evaluate
   `D_response = E_QM[Q] - E_QM[0] - C_density[rho_0,Q]`.

The embedded total already contains core/environment electrostatics; adding
the direct Coulomb term again would double count it. No CPCM or PB energy is
present in step 2. No old reaction-field energy is silently paired with a new
density and called a validated global score. This experiment yields charge-
representation and electronic-response diagnostics, not a new affinity class.

Report all four endpoint components, their Ca-minus-La differences, and the
qm36-minus-qm33 differences. No fitting, classification threshold, or success
selected from subsets. A residual comparable to the existing 18 kcal/mol
problem would argue against pursuing this approximation without a new physical
idea. A smaller residual is motivation for a separately recorded global test,
not a pass of the failed v1 gate.

Expected compute: four utility calls, likely seconds to minutes based on the
completed ESP job, then four small-core endpoints. The prior four vacuum/MBIS
endpoints used 59,648 allocated CPU-seconds; embedding may change convergence.
Request four 16-rank workers on 64 allocated CPUs, preserving existing MPI fixes.
No GPU, new solvent sweep or full-protein quantum calculation is included.

## Primary interface evidence

- [ORCA external point charges](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html#inclusion-of-point-charges)
  documents the electronic/nuclear coupling and omission of external-external
  interactions unless DoEQ is enabled.
- [ORCA vpot utility](https://www.faccts.de/docs/orca/6.1/manual/contents/utilitiesvisualization/utilities.html#orca-vpot)
  documents potential evaluation from GBW/density at user-supplied Bohr points.

Retrieved source pages are preserved under
`workspaces/density_embedding_20260916/source_inspection/`.
