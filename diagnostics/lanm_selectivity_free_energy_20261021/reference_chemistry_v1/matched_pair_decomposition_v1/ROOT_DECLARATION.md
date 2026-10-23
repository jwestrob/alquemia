# Matched pair versus cooperative water response

## Question and fixed inputs

Does the failed W8→W9 La/Dy coordination contrast arise mainly in individual
metal–water interactions, or in the additional interaction of the full shell?
Use only the four already consumed central W8/W9 La/Dy geometries and the fixed
terminal 17676 model selection. The reference counterpoise diagnostic is complete:
its +0.324723 kcal/mol shift increases the discrepancy from 4.338 to 4.663.
It supplies context, not a replacement target or a corrected decomposition.

## Finite new inventory

1. **Thirty-four quantum dimer single points:** independently extract each
   original water with its metal, retaining the literal coordinates and atom
   order. Eight and nine waters for each of La and Dy give 34 four-atom inputs.
   Do not average, rotate, optimize or replicate one orientation by symmetry.
2. **Two quantum bare-ion single points:** one La and one Dy at the common
   original metal position, using exactly the matching own basis/core/state.
3. **Thirty-four native polarizable-model energy evaluations:** those exact
   four-atom geometries, with the frozen 17676 parameters and native analyze
   route. The isolated model ion has zero interaction energy under the stated
   nonperiodic Hamiltonian; justify this from the implementation. No fitting,
   numerical derivatives or new force qualification is part of this panel.

Reuse all four original complex quantum energies and four selected model
energies. No new water-only or shell-only calls: their terms cancel between
metals at each water count because their coordinates and chemistry match.
No counterpoise fragments, changed core/electronic recipe, extra distances,
reserved validation outcomes, bulk dynamics or protein sampling are included.

## Consequential settings and execution

Reference chemistry remains the original ECP-PBE0-D4/def2-TZVP, def2/J, RIJCOSX,
DefGrid3 recipe with the exact source metal orbital/auxiliary bases and ECPs,
charge +3 and singlet effective state. La/Dy retain their original 46/55-electron
cores. Preserve all remaining SCF settings, stopping metadata and numerical
limitations from the qualified native adapter. Dimer and ion calculations are
electronic single points with the recipe's dispersion treatment, not gradients.

Root submits only after exact input preparation, bounded saved controls and
independent source review. Each reference job gets four MPI ranks under
`standard_expand`, scheduler-default memory/time, allocation-derived MaxCore,
one native invocation and no automatic retry. Schedule the independent cells
in parallel within the current shared limits. The model uses four or fewer
single-thread workers in a separate allocation. Its bound is 34 analyze programs.
Root owns all submissions. Preserve immutable results, unavailable cells and
actual scheduler/application costs. Preparation/review uses at most one local
numerical thread per agent, inside the authorized 16-thread aggregate allowance.

The previous small water reference timings suggest roughly 1.6 allocated
CPU-hours for the 36-reference panel; this is an estimate, not a time cap or
guarantee. Do not start more cells to meet an expected outcome.

## Predeclared arithmetic and interpretation

For either method, define `C=(E_Dy,W9−E_La,W9)−(E_Dy,W8−E_La,W8)`.
Let `S_M,n` be the sum of the literal isolated metal–water dimer energies for
that source. The paired sum of individual interactions is

`P = S_Dy,9 − S_La,9 − S_Dy,8 + S_La,8 − (E_Dy,ion − E_La,ion)`.

The ion difference is necessary because there is one more water at W9. Own-water
energies cancel between metals; they must not be numerically substituted with
zero before proving the coordinate/state correspondence. The cooperative
remainder is `B=C−P`. Report model/reference C, P and B, their signed differences,
and exact closure `error(C)=error(P)+error(B)`. Keep electronic and D4 terms
separate in the quantum result. Missing or inadmissible terms make the affected
complete sum unavailable; no partial-denominator averages or silent exclusions.

Large error in P supports prioritizing the individual-interaction description;
large error in B supports examining cooperative response. Report cancellation
or comparable errors explicitly instead of assigning a unique dominant cause.
No new fit or validation release follows automatically from either outcome.

This is an own-basis fragment decomposition under fixed provisional references.
Its B includes finite-basis nonadditivity and cannot uniquely identify
polarization, charge transfer or a specific missing quantum term. The previous
metal–whole-shell counterpoise cycle is separate: never insert its ghost-ion
energies into this own-basis formula or label a mixed-basis B as corrected truth.
W9 rearranges the other waters, so the contrast is not an isolated ninth-water
binding energy. The full exchange free-energy objective remains unmet.
