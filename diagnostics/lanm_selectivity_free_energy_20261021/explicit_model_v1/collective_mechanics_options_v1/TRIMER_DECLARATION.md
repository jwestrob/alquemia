# Finite common-parent-basis three-body discriminator

## Question and fixed selection

Does the missing collective energy/electric response already occur when two intact consumed waters interact simultaneously with La or Dy? Preserve the failed terminal17676 six-parameter model, the original PBE0 reference recipe, all coordinates and prior failures. No fitting, geometry relaxation, new field, EXCHPOL, reserved data or target change.

For each original W8 and W9 geometry, enumerate unordered water pairs. Select the minimum and maximum O–O distance using distance rounded to 1e−9 Å for ties, then the lexicographically smallest one-based water indices. La and Dy parents must have exactly matching coordinate arrays; apply the same pair identities to both metals. This selection uses geometry only. Keep both selections and all four parents even if results are uninformative.

| Parent | Near pair | O–O Å | Far pair | O–O Å |
| --- | --- | ---: | --- | ---: |
| W8 | 1,5 | 2.999312242227 | 1,6 | 4.717427908683 |
| W9 | 1,5 | 2.505029439719 | 1,3 | 4.683748498798 |

## Observables and exact seven-term cycles

For each parent and selected pair `(i,j)`, compute electronic energy and total lab-frame dipole with identical coefficients:

```
X3(M,i,j) = X(MWiWj) − X(MWi) − X(MWj) − X(WiWj)
             + X(M) + X(Wi) + X(Wj)
```

Every quantum term uses the same full original parent orbital and auxiliary basis, at unchanged parent coordinates, with the complementary atoms ghosted. Each charge-neutral combination is origin invariant; do not subtract charged COM-centered printed moments without their exact charge/origin correction. Retain full vectors, parent frame, native origin, total/electronic/nuclear components, ghost populations and actual print bounds. Prefer the existing independent full-precision property route where available. Printing uncertainty is not physical-reference uncertainty.

E3 is the third-order increment in a many-body expansion over intact molecular units. The existing `B=C−P` equals the sum of all metal-containing increments of order three and higher. These eight selected E3 values are samples, not a sum approximation to B. The same algebra applies to dipoles. Original full C/P/B and their failures remain unchanged.

New quantum calls are electronic PBE0 without D4, matching common-basis68. Do not allocate the original shell D4 value to chosen trimers. Report electronic E3 explicitly; if a future dispersion decomposition is wanted it needs its own declared, consistent seven-term data. Model component E3 and total E3 are retained, with the expectation that unchanged pairwise permanent/vdW and bonded terms cancel at printed precision. Do not call all quantum E3 polarization.

## Exact finite inventory and resource proposal

- Eight `metal_two_water` fragments: two metals × two parents × near/far; seven real atoms, charge +3, effective singlet, 28 electrons.
- Eight `two_water` fragments: matching six real atoms, charge 0, singlet, 20 electrons; metal ghost keeps its original basis and has no ECP potential.
- All original parent centers remain: W8 25 centers/391 AOs; W9 28 centers/434 AOs. Original Aux-J counts are 919/937 for W8 La/Dy and 990/1008 for W9 La/Dy. Real metal ECP is 46 electrons for La and 55 for Dy; ghost metal has zero ECP core.
- Reuse exact admitted common-basis metal–water, water-only and ion fragments for the other five terms. No own-basis substitutions, rotational replication or cross-parent ghost-ion cancellation.
- Reuse model78 monomers/dimers/ions. Prepare all 16 literal new model fragments with exact frozen keys, parameter bytes, masses, complete water bonds and local frames. No assumption that equal water-only energies across metals permits an unverified cache substitution; at most 16 model programs.
- Native reference ceiling: **16 single points**, four one-thread MPI ranks per cell, default scheduler memory/time, allocation-derived MaxCore `floor(0.75*granted_MiB/4)`. Same PBE0/RIJCOSX/DefGrid3/VeryTightSCF/PModel/MaxIter300, original ECP/orbital/auxiliary data. Retain native stopping and density/DIIS flags; no silent retry/state change.
- Model ceiling: **16 analyze EM calls**, at most four single-thread workers/four allocated CPUs. `POLAR-PRINT` and `SAVE-UINDUCE` reporting only. Zero new testgrad/Cartesian-FD/force calls. Finite nonnegative induction residual ≤1e−8 Debye; complete site counts and parameter-derived native labels. Use the already qualified native-announced `.uind` readback.

Require all 16 new rows per method and every explicitly reused term for a complete comparison. Missing/failing terms stay unavailable; never form a partial cycle or choose the better pair. Root owns execution and any global failure stop, scheduling, and actual accounting. Preparation is not a native release.

## Prospective numerical and structural controls

Exact selection replay; all parent coordinates unchanged; full ghost/AO/Aux/ECP/electron inventory; literal retained atom mapping; complete water groups and native local frames; original keys/parameter bytes; same selected six-vector; per-cycle atom and pair coefficient cancellation; charge neutrality; finite signed arithmetic; malformed/missing terms fail closed. No raw-admission campaign is repeated during preparation; bind existing independent admissions and source pins. Native run/collector adapters must reuse existing qualified parsers and preserve immutable receipts.

For the model cycle the sum of absolute site coefficients is 28, so six-decimal induced-file printing contributes at most **14e−6 Debye per Cartesian component**. Seven three-decimal total-M rows give 0.0035 Debye per component plus the existing COM arithmetic allowance, retained as a secondary consistency check. Seven eight-decimal energies contribute at most 3.5e−8 kcal/mol from printing alone. Quantum bounds are calculated from actual precision/charges without relaxing previous criteria. No new physical acceptance threshold is invented for this diagnostic; report signed energy and vector discrepancies with their scope.
