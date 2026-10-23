# Independent water orientation validation — scheduling amendment

## Scientific question and frozen scope

Does the already selected17605 water model reproduce the four original reserved
aquo8 water-H orientations (minus/plus15 degrees, La/Dy)? The17621 donor fit alters
no represented water-system pair. Independent assessment is recorded in
`paired_water_independent_validation_scope_v1`. This branch advances those four
water slots while donor development continues. It changes execution order only;
it does not relax the full16/24 validation criteria or change any outcome threshold.
Original Stage A declaration, code, configurations and failed historical fits stay
unchanged. No reserved numerical outcome has been read before this declaration.

Bind only actual17605 result6027e319a261eb85171fb06e30aa3e519411f6371a7905eeca0e8a79812c48fc,
its independently reviewed passing six-parameter selection, and exact17612 transfer
e3eeb870aac22f408c8b38a5f050404cfee5a8f4c0468633c07122fbc5e07254.
The latter provides passing consumed water Cartesian checks: both aquo8 center
max-component errors are below9.1e-8 kcal/mol/A. Its donor failures remain explicit.
Do not bind a partial17621 iterate, change water values or rerun consumed checks.

## Inputs and numerical calls

Retain original reserved geometry panel
`eed99b0fb6dbe7b49857d68853fe809481f6b657c928a7e75197d77c61f3432b`
and native24 panel
`75d3fc69703c8a5acca61f547e3fcf5d614cfe6ac653d2f3ecd3f34a93bcaf5a`.
Only IDs `validation_aquo8_orientation_{minus,plus}_{La,Dy}` run. All coordinates,
states, atom ordering and water topology remain exact. Atom types are36,37 and
500/501; no donor interaction recipient is present. Use exact17612 aquo8-center
key/parameter files and its selected numerical evidence. New model inputs differ
only in the original reserved coordinates; no relaxation or new water geometry.

Maximum four model E/F cells: analyze and analytic testgrad for each, eight native
programs total, no new Cartesian finite differences. Maximum four reference
EnGrad calls using unchanged original REFERENCE.py, ECP-PBE0-D4 recipe/assets,
La46/Dy55 cores, effective singlets and physicalLa1/Dy6 state labels. No fragments,
SCF retries, changed reference settings, parameter search or extra orientations.

## Fixed scoring and interpretation

Reuse unchanged StageA SCORE.py1c4bb0fa and rigid-body helperbf65a20a. Supply all16
original metadata slots and explicit unavailable donor numerical records. The
water subset has exactly11 comparisons: two aquo-referenced Dy-minus-La contrasts,
three plus-minus works (La,Dy,paired), and six rigid responses (two configurations
byLa,Dy,paired). Use exact17612 central-aquo model anchors and unchanged admitted
central-aquo reference energies. Every energy absolute error<=1kcal/mol; every
rigid RMS<=2 and maximum component<=5kcal/mol/A. No averages may hide failures.

All4 water cells and all11 comparisons are required for `water_subset_passed`.
The unchanged full16 scorer must remain partial and full44 validation false until
the other12 StageA cells are evaluated and admitted. All24 original reserved slots
remain in the eventual inventory; this water readout never substitutes for donor
or protein validation. No outcome-guided candidate selection or refit is included.
Any later use of these outcomes for development consumes their validation status.

Later combined validation must reuse these exact four reference results. Model
reuse additionally requires exact effective water Hamiltonian equality, including
all represented type/class pairs; inactive donor-key lines cannot affect the water
system. Do not repeat these cells simply to satisfy a different wrapper. Declare
any necessary saved-input collector compatibility change separately before use.

A passing subset qualifies only these water orientations against the provisional
reference surface. Bulk hydration, protein mapping, reference-method uncertainty,
electronic-state populations and the exchange free energy remain unqualified.
A failing or missing result stays visible; it does not release automatic repair.

## Execution and ownership

Singer owns minimal4-cell reference wrapper and metadata preparation. Arendt owns
4-cell model binding/evaluation and exact consumed-water numerical replay. Root
owns combined scoring/collection, source review, release and submissions. Planck
provides independent source review when available. Local shares root<=2,
Singer1,Arendt<=2,Planck1, under16total. No source preparation opens reserved
numerical outcomes or runs native science.

After actual source/preparation review, root may submit four independent4-MPI
reference jobs and one4-CPU/4single-thread model job understandard_expand with
default memory/time, actual-grant reference memory, and completion watchers.
Collector uses1CPU only after all terminal. Record immutable inputs, failures,
actual calls and costs. No changes to running17621, its objective, water invariants,
stopping limits, three terminal collectors or existing production configuration.
