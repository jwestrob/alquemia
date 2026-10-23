# Metal–water-shell counterpoise diagnostic

## Result

**Metal–whole-shell basis sharing does not explain the failed coordination-energy
test. Correcting this contribution increases the model discrepancy from
4.338 to 4.663 kcal/mol.** All eight new fragments and four reused complexes
pass the frozen collector. Independent actual-result review passes 399 checks,
including raw energies and separate high-precision cycle arithmetic.

The quantity below is the fixed-geometry, paired coordination contrast
`(Dy W9 − La W9) − (Dy W8 − La W8)`. W9 also rearranges the other eight waters;
this is not an isolated ninth-water binding energy or an exchange free energy.

| Quantity | kcal/mol |
|---|---:|
| Original quantum reference | −0.534207 |
| Metal–whole-shell counterpoise shift | +0.324723 |
| Alternative reference after that correction | −0.209484 |
| Unchanged polarizable model | −4.872491 |
| Model minus original reference | −4.338284 |
| Model minus alternative reference | −4.663006 |

The unchanged original D4 contribution is +0.106638 kcal/mol. The new fragment
calculations contain no dispersion term. The electronic fragment cycle gives
`ΔG = −0.324723 kcal/mol`; the counterpoise shift is `−ΔG`. Own-ion terms cancel
within each metal across sizes, and own-shell terms cancel between metals at
each size. Exact source coordinates, states and basis definitions establish
these cancellations. No individual absolute corrected complex energy is claimed.

## Interpretation and decision

Close this finite eight-fragment diagnostic. Retain the original fitting targets
and both reference results. This particular basis-set error is small relative to
the observed discrepancy and moves it in the opposite direction. It provides
no reason to restart the same parameter fit.

The result does not establish a complete-basis reference, remove water–water
basis sharing within the shell, qualify forces, or settle frozen-core/electronic
state uncertainty. The next proposed mechanics test will separate matched
individual metal–water interactions from the additional response of a full shell;
its inventory and settings must be declared before any execution.

## Execution and provenance

- Native jobs 17701–17708: eight electronic single points, four MPI ranks each,
  `standard_expand`, default scheduler memory/time. All completed successfully.
- Saved collector 17712: one CPU, zero molecular calls; all 8/8 fragments and
  4/4 original complexes admitted, with no partial-cycle substitution.
- Native scheduler cost: **7,428 allocated CPU-seconds**. Collector: **2 allocated
  CPU-seconds**. Total: **7,430 allocated CPU-seconds** (2.064 CPU-hours).
  Native application-derived allocation time is separately 7,425.525 CPU-seconds.
  Historical reused-complex costs are not new computation.
- Frozen native source `93ffa5f9`, panel `c109ab96`, root release `1e26b109`;
  startup review `6aef8947` passed 136 checks.
- Collector source `52b047c7`, configuration `47e857a0`, 26 saved controls and
  independent source review `4092748e`; root collection release `0137b426`.
- Complete collection result SHA256:
  `3647298ae095a95ff1d2149eea2d273432c933cf33bc41e9928abc4a0c811e79`.
- Independent actual review SHA256:
  `e6bb2dae760fa2cb720c8626ca85b0ce2e046dcd6ad1d5127c328771ac0e53fd`;
  zero molecular calls, 3.007 seconds wall and 2.776 seconds process CPU.

Exact results and raw-output provenance are under
`workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_cp_v1/`.
Original failures, fitting targets, held validation outcomes and production
defaults are preserved.
