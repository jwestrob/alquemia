# Finite fixed-source one-water energy panel

Declared16October2026 before any energy call on the modeled one-water panel.
Protocol: `lbt3_fixed_source_one_water_composite_v1`. Root owns Slurm execution.
This follows the separately frozen `LBT3_WATER_PLAN.md`; its geometry, selected
rays, contact cutoffs and rejected slots are immutable.

## Question and fixed scope

Does transferring one represented water from the declared4.5Å outer position
to the2.6Å inner position change the La/Dy response in the two original fixed
LBT3 sources? Every state has the same whole peptide plus one neutral H2O:
259atoms, charge−2, target255, water256–258, physical LaIII multiplicity1 or
DyIII sextet6. Original256 coordinates remain exact; no peptide/metal/water
optimization, proton change, new placement, additional water or counterion.

Retain all16 declared slots: source7CCO then7CCN; within source
inner_plane0,inner_plane90,outer_plane0,outer_plane90; within pose La thenDy.
Use only slots geometrically admitted by the pinned actual preparation.
Rejected slots remain unavailable counterfactuals in the full16 denominator.
Never replace them or regenerate placements using energies.

## Exact energy and execution

E = native MACE-OMOL vacuum + standalone GFN2 ALPB(water,gsolv) − GFN2 vacuum.
Inherit the actual12500/12510 lineage: same checkpoint/omol head,float64,
embedding/version/native batch checks, no extra dispersion, no masks; exact
standalone xTB6.7.1 executable, original f-in-core Dy parameters and runtime.
Physical La1/Dy6 remain explicit. GFN2 effective unpaired-electron count0 and
multiplicity1 are a separate model state, not a physical Dy singlet claim.

Use unchanged GFN2 accuracy0.0001,500iterations,electronic temperature300K,
no restart, no point charges; JSON/populations and `--sp` energy-only output.
Water uses ALPB water/gsolv,230grid points, no ion screening. Reuse the existing
energy-only FD admission checks: actual state/settings, normal termination,
convergence, finite energy/populations, charge closure1e−5e,
JSON/stdout energy agreement1e−7hartree and populations agreement1e−7e.
No force or numerical derivative requested. GFN2 effective electron counts
come from the pinned parameter profile for the ACTUAL259 elements.
Units remain23.06054783061903kcal/mol per eV and627.509474 per hartree.

Reuse pinned `LBT3_FD.py::NativeMACE` initialization/evaluation and its
`admit_energy`. A private, recorded AST adaptation replaces ONLY the integer
atom-count256 with259 in `NativeMACE.evaluate` (two occurrences) and
`admit_energy` (one occurrence). Require those exact counts and the pinned
helper hash before creating new private callables. No frozen source or module
is mutated. No old256-atom energy/cache receipt is reused for a259-atom state.

One warm model per allocation. For each cell, native MACE finishes before
launching its two single-thread GFN2 medium workers; no CPU oversubscription
with the MACE thread grant. Record actual component calls, timings, native
model-load cost, versions, runtime pins and complete backend receipts. Maximum
panel cost16MACE+32GFN2 calls; these are ceilings, not claimed execution costs.
No quantum calculation, optimizer or fallback engine.

## Paired scout and remaining gate

`scout`: only7CCO/inner_plane0/La andDy. BOTH must be geometrically admitted
before model loading. If either is rejected, record the unavailable scout and
stop with zero calls; root must review before any alternative scout.

`remaining`: requires the exact prior scout result/config/implementation/plan,
both complete259-atom cells and their actual component/state/coordinate
receipts. Continue the other14 slots unchanged, scoring only admitted ones.
Gate on execution/physical admission, never energy magnitude or preference
sign. Reuse the two scout receipts exactly; never repeat their cells.

Each backend is attempted at most once per declared cell/medium. Failure leaves
that cell unavailable and its raw evidence intact; no restart, altered settings
or replacement geometry. A MACE failure stops later MACE calls in that stage.
An isolated GFN2 failure leaves that cell unavailable; other declared cells may
continue. Failure of either scout cell blocks the remaining stage.

## Reports and finite stopping criterion

Emit JSON and a short Markdown table covering all16 slots and component/call
availability. Per source and orientation report each metal's inner-minus-outer
transfer and H=(Dy_inner−Dy_outer)−(La_inner−La_outer), calculating balanced
component differences before unit conversion. Retain missing constituent cells
as unavailable; never substitute a different orientation or source.

Give both metals the same four-pose menu per source. Report relative energies,
selection and paired contrast of that common menu only when all8 source cells
are available. A partial panel has no complete-four-pose selector. Retain every
individual pose; no pooled source vote or favorable-subset result.

Predeclare1kcal/mol as interpretive resolution for H, not a statistical error
bar or classifier calibration. Report near-zero values and source/orientation
sign disagreements. Complete this fixed panel once, regardless of sign or gap
reduction. No automatic placement, water-count or torsion expansion follows.
This fixed-peptide test cannot rule out coupled water/Asn displacement.

All comparisons retain one explicit water. The outer state is NOT bulk water;
there is no water chemical potential, hydration population or binding affinity.
Never select raw minima across dry/wet inventories or label this a validated
Mex–LBT3 preference profile.

## Interface

`python LBT3_WATER_SCORE.py --stage scout|remaining --config CONFIG --output FRESH`

Config: `repository`, pins `implementation`,`plan`,`water_preparation_result`,
`solvent_result` (actual12510),`fd_result` (actual12521), optional `pins`.
Remaining adds `scout_result`. All geometry replay, energy and reporting run
inside Slurm. Root binds frozen configs and inspects the actual scout gate.
