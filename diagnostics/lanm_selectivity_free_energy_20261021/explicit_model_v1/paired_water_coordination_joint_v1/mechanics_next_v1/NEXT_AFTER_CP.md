# Next after the coordination CP result

## Decision: exact fragment additivity, with no parameter change

Perform one energy-only, fixed-geometry decomposition of the consumed W8/W9 La/Dy contrast into **isolated metal–water pairs** and a **metal-dependent cooperative remainder**. Root selected the literal native-coordinate inventory: 34 extracted dimers plus two own-basis ions for the reference, and 34 matching frozen-model dimer energies. Reuse all four original cluster energies. This tests which response is missing before choosing another functional term or parameter; it is not another fit, H-Thole probe, force campaign or validation test.

The completed whole-shell CP correction changes the reference contrast from −0.534206667 to −0.209484155 kcal/mol. The unchanged model is −4.872490500, so its discrepancy changes from −4.338283833 to −4.663006345. That particular BSSE correction does not explain the discrepancy. The alpha probes also did not provide an admissible repair: none improved the aggregate objective, and the favorable unconstrained projection required invalid/out-of-range parameters. Neither result establishes global model-family impossibility or a unique QM-error cause.

## The matched quantity

Let `E_MN` be the total energy of the exact saved metal M plus N-water cluster; `e_MNi` the energy after retaining only that metal and intact water i at the **same Cartesian coordinates**; and `a_M` the isolated own-basis metal-ion energy. For each method separately define:

```
C = (E_Dy9 − E_La9) − (E_Dy8 − E_La8)
p_M = sum_i=1..9 e_M9i − sum_i=1..8 e_M8i − a_M
P = p_Dy − p_La
B = C − P
error_C = error_P + error_B
```

`P` is the change in summed metal–water two-body interaction contrast. The free-water energies cancel between metals at each size because the water coordinates and Hamiltonian are identical. The own-ion term **does not cancel**: nine minus eight leaves one ion term per metal. `B` is therefore the change in the metal-dependent nonadditive remainder. Metal-free whole-shell energies also cancel in this paired quantity, so no new shell or free-water calculations are required. No individual absolute shell-cooperativity energy is claimed.

Use consistently uncorrected, own-basis PBE0-D4 energies for this identity, with the exact original frozen-f state, ECP/basis and numerical recipe. Retain electronic and D4 contributions separately. The completed full-shell CP contrast is a separate sensitivity result; subtracting raw `P` from CP-corrected `C` would not be a consistently CP-corrected decomposition. Existing ghost-ion energies cannot replace `a_M`. Pair-fragment BSSE and reference/state uncertainty remain unresolved by the old whole-shell CP check.

## Why this distinguishes useful hypotheses

The actual Tinker key uses fixed permanent multipoles, a pair-additive buffered 14–7 vdW term and self-consistent induced dipoles. In `ehal.f` the vdW energy is accumulated independently for each atomic pair (lines 163–210). `induce.f` first forms alpha times the permanent field (lines 232–245), then solves against the field of the other induced dipoles (lines 374–408). Its fixed alpha/Thole prescription can produce cooperative response, but does not guarantee that response matches the electronic calculation. The exact active outputs contain bonded, vdW, permanent-multipole and polarization terms; no separate charge-transfer or exchange-polarization term is present in these inputs.

At these fixed geometries, the metal-dependent contribution of every pair-additive vdW/permanent term is reproduced in the extracted-pair sum. Shared water-only terms cancel between metals. Thus the model's `B` should reduce to its polarization remainder, apart from printed precision. Changing an O-pair radius or epsilon changes `P`, not this nonadditive remainder. O-pair damping can change both. The weak alpha response does not establish that collective water response, charge redistribution, short-range exchange/penetration or nonlinear susceptibility is correct; the QM remainder includes these possible effects plus method/basis artifacts and D4 environment dependence.

Interpret the two signed errors together:

- If `error_P` explains the discrepancy while `error_B` is small, prioritize the isolated one-water interaction near the actual cluster distance. This localizes a pair-level error, **not uniquely an M–O vdW error**: the dimer also contains H interactions and polarization. Resolve its physical origin before refitting.
- If `error_P` is small and `error_B` carries the discrepancy, changing only an additive short-range pair term would compensate for a cooperative error. A subsequent response test should target the collective interaction identified here, rather than adding another radius/epsilon freedom.
- If both are appreciable or cancel, retain a mixed explanation. Do not select the more convenient component or claim a repair from cancellation.

The previously saved model contrast, vdW −5.29009248 plus polarization +0.41760198 kcal/mol, is not already this decomposition: isolated dimers are needed to split pair polarization from cooperative polarization. Component magnitudes alone cannot assign the QM error.

## Exact inventory and reuse

| Item | New work if separately released | Immutable reuse |
|---|---:|---|
| Reference extracted dimers | 34 singlepoint energies: (8+9) waters × 2 metals | Exact deposited-in-panel Cartesian subsets, unchanged electronic recipe |
| Reference own-basis ions | 2 singlepoint energies, La3+ and Dy3+ | Same metal basis/ECP and effective singlet, 8 explicit electrons |
| Frozen-model extracted dimers | 34 `analyze` energies; no `testgrad` | Exact selected17676 keys/prm; types 500/501,36,37,37; inherited water bonds |
| Frozen-model own ions | 0 native calls | Zero analytically: isolated charge has no partner, external field, induced dipole, bonds or pair energy in this vacuum Hamiltonian |
| Cluster energies | 0 calls | Four original reference and four selected17676 model records |

All 17 extracted metal–water geometries are congruent within 1.014e−12 Å in their six pair distances: M–O 2.50 Å, O–H 0.9572 Å, H–H 1.513900654528 Å. A rotationally reduced calculation would need only two representative dimers plus two ions, but would assume numerical rotational equivalence of RIJCOSX/grid energies. Root chose all 34 literal fragments to remove that assumption. No coordinate rotation or averaging is needed. The published water dimer is not reusable: its M–O distance is 2.288157949 Å and O–H distance 0.997846796 Å. W9 also rearranges the other eight waters; `C` is a coordination/arrangement contrast, not a pure ninth-water insertion energy.

The 36-reference/34-model inventory is the minimum literal-coordinate energy panel for this paired decomposition. No new cluster, free-water, whole-shell, geometry optimization, gradients, finite differences, fitting, donor cells or held-out outcomes are included. Full denominators and native admission remain mandatory; an unavailable term leaves the corresponding decomposition unavailable. Root will declare and release the executable panel separately.

## Reusable authorities

- Compact consumed ten-cell data: workspace `explicit_model_v1/paired_water_coordination_joint_v1/metal_alpha_response_v1/preparation_v2/DATA.json`, SHA `e2e58963297984305856210f93cd664f08903485d9413f5f7a3d10ba1a616250`. `rows` supplies original geometries/reference authorities; `baseline_model` supplies exact selected energies/components and `provenance.effective_input` supplies selected keys, not the old `rows.input.key`.
- Frozen17676 result `5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832`; independent selected-ten raw review `fd7ea77167c833c3990d7624ff136042ce03f3ccf9dc4169ba7cff2aa51bb6a2` passes raw admission while retaining scientific failure.
- Exact original four complex raw sources/energies are pinned in `reference_chemistry_v1/coordination_cp_v1/prepared_v4/PANEL.json`, SHA `c109ab96100e49633171a56c5f0b738522a565e088c6623eb257d99965c63a84`.
- Actual CP collection: `reference_chemistry_v1/coordination_cp_v1/collection/run_v1/result/RESULT.json`, SHA `3647298ae095a95ff1d2149eea2d273432c933cf33bc41e9928abc4a0c811e79`.
- Native `ehal.f` SHA `fd5a4a56349e4ad26f2aaae2eb12e8f780b092879e66f0ad88c4d706416d2740`; `induce.f` SHA `759e544790074aa8e7a32a0176c1ef979bceca3787252e58efcffe4a942c0e5f`; alpha/Thole coupling remains documented in the preceding support note.

This note used saved sources/evidence only. No molecular calculations, candidate changes or fresh validation outcomes were used; the original proposal and all failed fits remain intact.
