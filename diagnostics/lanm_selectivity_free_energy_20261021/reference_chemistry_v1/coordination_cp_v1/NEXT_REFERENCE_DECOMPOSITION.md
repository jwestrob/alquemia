# Matched pair versus cooperative hydration discriminator

**The smallest algebraic test needs four new quantum single points:** one representative exact cluster-extracted M–water dimer per metal and one own-basis bare ion per metal. All seventeen extracted dimers are congruent. A literal sum at every original Cartesian orientation needs 36 single points and avoids assuming numerical rotational invariance. Both implement the same energy decomposition; neither changes the original fitting targets or runs here.

## Exact geometry and available evidence

The consumed W8 and W9 geometries contain eight and nine intact waters. Saved-coordinate analysis finds all six interatomic distances of every extracted M/O/H/H group identical within **1.014×10⁻¹² Å**: M–O 2.500000 Å, O–H 0.957200 Å, H–H 1.513900654528 Å. La/Dy use identical coordinates. Each dimer can therefore be represented by the four unchanged coordinates of W8 water 1, with a recorded rigid mapping to the others.

The previously computed published dimer is incompatible: M–O 2.288157949 Å and O–H 0.997846796 Å. Its energy, correlated correction and displaced derivatives cannot substitute. No compatible own-basis bare-ion result is present in these consumed reference panels. Existing ion results retain ghost ligand/shell bases and cannot substitute for an own-basis ion.

Already available are four real W8/W9 ECP-PBE0-D4 complexes and eight whole-shell CP fragments. The latter shift the paired contrast by +0.324722512 kcal/mol and increase the selected-model discrepancy from −4.338283833 to −4.663006345. They constrain metal–whole-shell basis sharing, but do not decompose it into individual water interactions.

## One matched own-basis energy cycle

For each method separately, let C = (E_Dy9−E_Dy8)−(E_La9−E_La8). For every extracted water i, use its original M/O/H/H coordinates and the same electronic recipe as the full complex. Define

`P = [Σ_i=1..9 E(DyWi) − Σ_i=1..8 E(DyWi) − E(Dy)]`

`    − [Σ_i=1..9 E(LaWi) − Σ_i=1..8 E(LaWi) − E(La)]`

`B = C − P`.

P is the paired contribution of independent metal–water interactions; B is the paired remainder after their removal. The own-water monomer and intact-shell energies cancel between metals at each size because their coordinates, atoms and own-basis recipe match. **The bare-ion difference does not cancel:** the shell-size difference is one. Two admitted own-ion calculations are therefore required. No water-only or shell-only calls are needed for this paired remainder; individual absolute cooperative energies are outside this test.

Use original uncorrected PBE0-D4 consistently for C and the dimers. Retain electronic and D4 parts of P and B separately; isolated-ion D4 is the exact empty-pair/triple sum, zero. Do not subtract own-basis dimer sums from the CP-corrected cluster interaction. B includes every electronic nonadditive contribution and environment-dependent D4 contribution; it does not uniquely identify polarization, charge transfer or another physical mechanism. Finite-basis errors also remain.

Compare P_model−P_reference and B_model−B_reference; their sum must recover the original C discrepancy. This distinguishes pair-dominated, remainder-dominated and mixed errors at these consumed configurations. Agreement is conditional on the provisional quantum method and decomposition, not a free-energy or force qualification. Arendt owns the matching fixed-parameter native-model proposal; no optimization or force dataset is needed for this first question.

## Finite count, cost and numerical choice

| Implementation | New quantum cells | What is established |
| --- | ---: | --- |
| Congruence-reduced own-basis test | 2 dimers + 2 ions = **4** | Algebraic minimum: because 9−8=1, P reduces to the paired single-dimer interaction. Relies on the electronic calculation's numerical rotational invariance. |
| Literal original-orientation own-basis test | 34 dimers + 2 ions = **36** | Same identity at every archived Cartesian orientation, without numerical rotation replication. |
| Every pair in its full cluster ghost basis | 34 MWi + 34 Wi = **68** | Strict pair/remainder decomposition using the existing metal–whole-shell CP basis convention; four ion and four shell ghost terms are reused. |

The four-cell version is the smallest useful discriminator if root explicitly retains its rotation/grid limitation. **For attribution of the exact native-coordinate discrepancy, prefer 36:** molecular congruence does not by itself prove finite RIJCOSX/DefGrid3 energies are numerically identical under rotation. This is not an additional geometry scan; all 34 dimers are exact atom subsets. The 68-cell full-basis option is unnecessary for this first localization and is substantially more costly.

Keep original frozen-f ECP46/55 and metal/light/auxiliary bases, PBE0, RIJCOSX, DefGrid3, VeryTightSCF, PModel and MaxIter300; use energy-only single points. Dimers have 90 AO/18 electrons, bare ions 47 AO/8 electrons, all charge +3/effective singlet, with physical La1/Dy6 metadata retained. New native admission remains mandatory. There are no new full-cluster calls, geometry changes, fitted parameters or automatic retries.

Actual existing dimer EnGrad timings were 39.55/40.30 s at four MPI ranks. At that rate, four new calls cost about **0.18 allocated CPU-hours**, and 36 about **1.6 CPU-hours**; energy-only calls may be cheaper, but these are planning estimates. Existing CP8 required 7,425.52 native-elapsed×rank CPU-seconds (2.06 CPU-hours) with 391/434 AO. Applying its average to 68 full-basis fragments gives roughly 17.5 CPU-hours, a coarse comparison rather than a bound; electron counts differ. Root chooses the finite count before inputs or results.

## Evidence and scope

- Geometry/reference inventory: `reference_chemistry_v1/prepared_v2/PANEL.json`, SHA `d5929526cf70ac957c847a5936ed229e9edf1322764b0cc287a59292e01204f9` under the 20261021 workspace. Its four original complex results are pinned in the coordination CP panel `c109ab96100e49633171a56c5f0b738522a565e088c6623eb257d99965c63a84`.
- Coordinate construction is the frozen original `REFERENCE.py::aquo`, source SHA `eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129`; comparisons above use saved rounded coordinates, not regenerated fixtures.
- CP evidence: `coordination_cp_v1/collection/run_v1/result/RESULT.json` and its unchanged report. Dimer timing evidence: original `pilot_v1/tasks/ecp_pbe0__water_central_{La,Dy}/execution/RESULT.json`, jobs 17255/17256.

One-thread saved-coordinate/scalar analysis only; zero native calls, no executable inputs, no target changes and no held validation outcomes accessed. The four-cell, 36-cell and 68-cell identities are proposals, not computed results.

**Subsequent root decision:** prepare the literal **36-cell own-basis energy-only** branch as `matched_pair_decomposition_v1`: all 34 exact dimers plus two bare ions, no symmetry replication. The model author will use this common geometry inventory for 34 fixed-model analyze calls. This records a preparation decision, not execution or results; separate declaration, frozen inputs, review and root submission follow.
