# Matched fixed-FF force feasibility at four LBT terminal starts

**Feasible with zero new molecular calls for the four own-metal starts.** Raw vacuum MACE forces and exact fixed-FF native ALPB gradients are archived for geometry_11/7CCO/La, geometry_13/7CCO/Dy, geometry_15/7CCN/La and geometry_17/7CCN/Dy. Energy/state and coordinate order are already joined by the12720 inventory. This inspection performs no force assembly, geometry calculation or search.

## Exact model and assembly

The whole256 region has255 peptide atoms plus one trivalent target at index255; charge−2, physical multiplicityLa1/Dy6, same OMOL float64 checkpoint9b64b4fd. Native ff19SB charges total−5 on the peptide, plus formal+3 on the target; all four serialized charge files have SHA256 `7baa18e804fd93edce495c2a1b649738669d448221a2b09b1dc908256e5cfed5`. The free charged termini are the existing conditional preparation, not a new experimental claim. There are no caps, waters, spectator ions or external region. Field-response and MM-ledger energies/gradients are exactly absent by partition; adding a peptide MM vacuum term would double count.

Use the12720 energy unchanged and assemble `g_matched = -F_MACE * 23.06054783061903 + g_ALPB * 627.509474 / 0.529177210903` in kcal/mol/Å. Raw MACE is force in eV/Å; native ALPB is dG/dR in Hartree/Bohr,256 rows×3. The old POINT `GRADIENT.npy` includes self-consistent GFN solvent and must not be reused. No old GFN correction, charge-response force or calibrated offset belongs in the matched gradient.

`LBT_MATCHED_MODEL_V2.py::solvent` archives gradient.dat after the qualified native wrapper zeroes its gradient and calls born%addGradient on the same cavity/fixed charges used by getEnergyParts. “Archived only” describes its prior use, not missing gradient data. The wrapper lineage is native force qualification12566, fixed-FF Mex12685/12709, then unchanged256-atom LBT12715/12720. This establishes a reusable analytic component, not a numerical claim about the new assembled residual.

## Own-state artifact pins

Paths in this table are relative to `workspaces/lanm_coupled_response_20261016/`. Each MACE/POINT directory contains whole.xyz and PHYSICAL.json pinned by the inventory. Each ALPB directory contains START.json, physical.xyz, charges and its RESULT.json.

| Candidate/state | MACE and POINT directory | Raw MACE SHA256 | Matched ALPB directory | gradient.dat SHA256 |
|---|---|---|---|---|
| geometry_11/7CCO/La | `lbt3_basis_scout_v1/execution/7CCO_La_candidate/evaluation_003` | `0c41a148d5296b9813ea32208a749a74ca9493eedb0dbe729c308fc5014bad93` | `lbt_matched_pool_v1/execution/geometry_11_La` | `60e680e2466257ceefd6749df030cc36cbee10cda7fd51e39f0d1cf34d1d8f69` |
| geometry_13/7CCO/Dy | `lbt3_basis_remaining_v1/execution/7CCO_Dy_candidate/evaluation_012` | `bbefd83aafa5d8ff87df4dac2e7d8410a41e5e81bbe44c4eb12ce8f51e0679a8` | `lbt_matched_pool_v1/execution/geometry_13_Dy` | `f676a6897a8025c5829680e23d0cdfa76a70420b40b512f27af1e2056ac8ef22` |
| geometry_15/7CCN/La | `lbt3_basis_remaining_v1/execution/7CCN_La_candidate/evaluation_012` | `1ed06abbaf784e409f686587349a271d1c84acfbeb4c6d2799db0f13d6d89099` | `lbt_matched_pool_v1/execution/geometry_15_La` | `82e222975f30168b8534151e47e3bd61066c247613bc903d8e8faa20d00ea4ba` |
| geometry_17/7CCN/Dy | `lbt3_basis_remaining_v1/execution/7CCN_Dy_candidate/evaluation_013` | `6205ca1528dd5bd1e221591cbabadd9f8212afee04887262ad30e6ca9a7a6968` | `lbt_matched_pool_v1/execution/geometry_17_Dy` | `92eeff2133885db86c10a16918b06687a7d52f40f384d3af21ca89055af9ec91` |

Parent12720 result: `lbt_matched_pool_v1/execution/RESULT.json`, SHA256 `ee867343321293d83c8845aee7ad6b7f5e6e4661a738f1c7c70c015ae1262928`. Exact full POINT/XYZ/PHYSICAL/ALPB result pins are in its four matching rows and `LBT_MATCHED_POOL_INVENTORY.json` (`aa81b392…`). Native binary: `lbt_matched_model_v2/execution/lbt_fixed_solvent`, SHA256 `d7e6d151316ab8d555388921700efe9c4d4ea31bdebc5dd2206ef00c516d548a`; source `LBT_FIXED_SOLVENT.f90` SHA256 `6dc6fca05cb0ca9b16fb9dd56cdb65d2a07d9322a190af5a99f765e8f71411d6`.

## Original-source headroom and stop reasons

The following exact headroom expressions use saved scalar diagnostics only; they are not recomputed geometries or predicted admissible step lengths. A nearly exhausted outward radius still permits inward or tangential motion, subject to the original full guard. No source rebasing or envelope expansion is implied.

| Start | Source displacement headroom, Å | Target headroom, Å | Total torsion headroom, rad | Saved stop |
|---|---|---|---|---|
| 7CCO/La | `0.40 − 0.2252202882832968` | `0.20 − 0.0887399206714831` | `0.20 − 0.19892161967094218` | no_admissible_decreasing_trial after3 accepted |
| 7CCO/Dy | `0.40 − 0.39972827959033214` | `0.20 − 0.1691256795322797` | `0.20 − 0.19093285984860334` | accepted_step_limit after12 accepted |
| 7CCN/La | `0.40 − 0.39994430847323903` | `0.20 − 0.12710936140781537` | `0.20 − 0.17524921987651218` | accepted_step_limit after12 accepted |
| 7CCN/Dy | `0.40 − 0.35785679797286973` | `0.20 − 0.1345240667393522` | `0.20 − 0.19962350552031166` | accepted_step_limit after12 accepted |

7CCO/Dy and7CCN/La have almost exhausted the source displacement radius;7CCO/La and7CCN/Dy approach the torsion limit. All saved boundary flags are false, but the numeric margins are small. 7CCO/La final normal trials were blocked by torsion; two geometrically admitted target-coordinate fallbacks then failed the old-objective Armijo condition. Its stop is therefore not proof that every matched-objective direction is geometrically unavailable. The other three stopped at12 accepted steps.

The exact SEARCH.json pins are in the candidate-arm entries of `lbt3_basis_pool_v1/execution/RESULT.json` (`b3d1b4db…`). They retain full_q, original-source coordinate scaling, offset coordinates, guards, termination and last trial evidence. The ten active indices are `[61,62,63,23,24,8,37,9,38,32]`: target xyz; A7 CA–CB/CB–CG; A3 CA–CB; A11 CA–CB; A3 CB–CG; A11 CB–CG; A10 N–CA. Existing source/target/total-torsion/trial bounds are0.40Å/0.20Å/0.20rad/0.05Å. Keep the original-source metric and guard; reset optimizer history because the energy objective changes.

## Minimal next useful action and limits

A single saved-data Slurm preflight can verify all four exact energy/gradient joins, assemble the matched gradients and project them through the existing current Jacobians. It needs zero MACE, ALPB, GFN2 or MM calls and zero new geometry. Inspect whether the unchanged search rule admits a predicted descent direction before any energy request; do not repeat an outward bound-exhausted direction blindly. No matched residual or feasible new direction has been calculated by this note.

All eight energy cells at these four geometries already exist. Four opposite-state raw MACE records are energy-only (`forces_eV_A=null`); all corresponding ALPB gradients are already archived. Complete paired gradients would require four fixed-state MACE force evaluations, but those are unnecessary to initialize the four own searches. Any later new geometry needs a fresh MACE force/energy and native fixed-FF ALPB call only; old molecular caches must not be labeled with the new objective.

Matched-objective response would address the current sampling asymmetry, not guarantee improved accuracy or a favorable sign. Fixed charge, conditional terminal/occupancy assumptions, source disagreement and absent equilibrium sampling remain. No code, molecular run, new state or submission was created by this feasibility inventory.
