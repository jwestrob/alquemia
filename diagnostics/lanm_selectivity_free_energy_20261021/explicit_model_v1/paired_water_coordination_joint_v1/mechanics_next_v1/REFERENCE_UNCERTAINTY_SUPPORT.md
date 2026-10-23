# Reference support for the consumed W8→W9 contrast

**The −4.338284 kcal/mol model discrepancy is established against the chosen ECP-PBE0-D4 target, but is not resolved against physical reference uncertainty.** No existing correlated, counterpoise, basis-convergence or matched all-electron energy cycle qualifies this exact W8/W9 contrast. Recommend one eight-calculation metal–whole-water-shell counterpoise discriminator. This is a proposal only: no inputs, native calls or submissions were created, and no target changes follow automatically.

## What the archives establish

Define ΔX = X(Dy,W9) − X(La,W9) − X(Dy,W8) + X(La,W8). The four admitted original energies give ΔE = −0.5342066673 kcal/mol; the selected 17676 model gives −4.87249050. W9 also rearranges the original eight waters. This is a fixed-configuration coordination contrast, not a pure ninth-water insertion or hydration free energy.

| Existing evidence | Coverage and limitation |
| --- | --- |
| Four W8/W9 ECP-PBE0-D4 energies and gradients | Correct native inputs/states and raw results are admitted. Their final-COSX contribution to ΔE is only +0.0157687 kcal/mol, so this known stage effect cannot explain the 4.338284 discrepancy. It does not bound other numerical or physical errors. |
| Canonical CCSD(T) versus matched PBE0-D4 counterpoise | One-water paired gap +0.467767 kcal/mol; dry water→amide gap −1.401006 and water→bidentate-acetate gap −0.254895. These are same-basis, frozen-f dimer comparisons, not an error bar transferable to W8/W9. |
| Water dimer two-step force test | Native PBE0 analytic/energy derivatives agree in the sampled direction; paired CC-versus-PBE0 CP slopes still differ by about 1.562 kcal/mol/Å. No W8/W9 energy-cycle qualification follows. |
| All-electron scalar-ZORA water dimers | Paired radial slope differs from ECP by −2.074505 kcal/mol/Å, with a La sign change. Metal/light/auxiliary bases, core representation and Dy spin also change. There is no matched AE fragment energy cycle or AE W8/W9 result; the recovered Dy final-grid update also lacks derivative qualification. |

Finite-basis/BSSE, functional/correlation, frozen-f/core, scalar-relativistic and spin–orbit/state uncertainty therefore remain open. Previously proposed larger-basis and correlated hydrated tests are plans, not completed evidence.

## One finite discriminator: two-fragment M | Wn counterpoise

For each of the exact four original complexes, evaluate **two electronic-only fragments in its full complex basis**:

| New role | W8 | W9 | Charge / effective multiplicity |
| --- | --- | --- | --- |
| Real metal, every water atom ghosted | 8 electrons; 391 AO | 8 electrons; 434 AO | +3 / 1 |
| Real intact Wn shell, metal ghosted | 80 electrons; 391 AO | 90 electrons; 434 AO | 0 / 1 |

Two metals × two shell sizes × two roles = **8 new single points**, with **4 original complex results reused**. The shell is the entire interacting water cluster, not separate water monomers. This tests metal–shell basis sharing; it does not remove intra-shell BSSE or establish a basis-set limit.

The archived La/Dy coordinates and atom order are identical within each n. For each metal, W8/W9 use the same metal basis, auxiliary basis, ECP and effective state. Consequently the own-basis metal energies cancel across n and the own-basis intact-shell energies cancel across metals. Let G(M,n) = E(metal in full basis) + E(intact shell in full basis). Then:

`C_CP,electronic = ΔE_complex,electronic − ΔG`

`C_CP,total = C_original,total − ΔG`.

No extra own-basis fragment calculations are needed for this paired contrast. These cancellations do not supply individually resolved binding energies. Original real-complex ΔD4 = +0.1066383 kcal/mol stays fixed and is reported separately; **no D4 is evaluated on ghost bodies**. The original target remains archived unchanged.

Reuse the actually exercised `canonical_water_cp_v2` → `dft_water_cp_v1` → `donor_cp_v1` conventions: colon-marked ghosts; real metal retains ECP46/55; ghost metal retains its NewGTO and auxiliary basis but removes NewECP and uses DelECP. Require actual zero ghost nuclear charge/electron/ECP contributions, full per-center basis inventory, correct charge/spin and native convergence/energy closure. ORCA documents retaining basis functions while removing ghost nuclei and electrons. [ORCA 6.1 counterpoise documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/counterpoise.html).

Keep the original PBE0, def2-TZVP light basis, pinned metal/auxiliary assets, RIJCOSX, DefGrid3, VeryTightSCF, PModel, MaxIter300 and no-autostart settings. Single points omit EnGrad and ghost bodies omit D4. Preserve actual SCF stopping signals, residuals and final-COSX metadata; no new force-accuracy or physical-state claim. Physical La singlet/Dy sextet metadata remains distinct from the effective closed-shell Hamiltonian.

## Cost and decision

Use four MPI ranks per independent cell, standard_expand with default memory/time; the established four-rank grant used 16,000 MiB and MaxCore3000. Original jobs 17306–17309 took 375, 380, 476 and 469 scheduler seconds: 6,800 allocated CPU-seconds total. Eight calls at those full-gradient rates imply roughly **3.8–4.2 allocated CPU-hours** as a planning estimate, not a measured upper bound. Ghost single points omit gradients and have fewer real electrons, but retain the full 391/434 AO basis and may converge differently. Record actual costs; no automatic retry or expansion.

Report ΔG and both original/CP contrasts, including the resulting model discrepancy. A correction comparable with the current gap would challenge target robustness to this specific basis-sharing effect. A small correction would constrain that explanation only; it would not validate PBE0, remove intra-shell BSSE, qualify forces, or release held fresh validation. Do not fit to or replace the target automatically.

## Exact evidence

Paths below are relative to `workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/`, unless stated otherwise.

- Original `ecp_remaining_v1/tasks/ecp_pbe0__aquo8_central_La/execution/RESULT.json`: `f1be7b98603569899240fed0d80802785f641a346e2a54926cf01e33d3246daf`; W8 Dy: `c1ff6777a8a64dcc7816c0ab1e3314e43d6221d56c50c35d1d0bf3c30b524482`; W9 La: `dae81a20caf7851de3bf000a832108a55fb5fce0ca657b307d40063632c5847c`; W9 Dy: `1fca94594167c125794b8938c578d969d87675c0350b96fdddb2c655195ffd81`.
- Stage audit `ae_water_solver_recovery_v1/paired_comparison_v1/ECP_COORDINATION_STAGE_REVIEW.json`: `8c478786b53ddcf53deda1f74b9b46ce1433a29481c699e34ea50b5846b0e7e1`; AE paired comparison RESULT: `c09e6906d7d6e11e2834179031b0c134601f2079417ad83abe6c58108a346ab8`.
- Dimer water CP `dft_water_cp_v1/final_collection_v1/RESULT.json`: `6465f8c9002c31ee9c79a945152476c91ad740da0a5ddca4e2327d98321bda64`; donor CP `donor_cp_v1/full_collection_v1/RESULT.json`: `c9ca7b157973c513fc808f3caf6976463da7863ebcdaeb2fc19f296fd9b74348`.
- Reusable diagnostic source `dft_water_cp_v1/REFERENCE.py`: `b1f35c46cf4512b6088d5a81183cc6ab1532fd997fee05ab5fd904b3c4595089`; `donor_cp_v1/REFERENCE.py`: `f4baaeddce09cb26551e845f05e7f3c1bfff604a84fc02c666908ada2420dea2`. The former imports the canonical native parser; the latter adds strict per-center basis and ghost checks.
- Joint 17676 RESULT: `5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832`; selected raw/score review: `fd7ea77167c833c3990d7624ff136042ce03f3ccf9dc4169ba7cff2aa51bb6a2`, raw pass with scientific controls false.

Saved-source assessment and scalar bookkeeping only; one local thread, zero native calls. Planck independently checked the geometry/state/basis cancellation joins. The earlier 17670 lost pre-submit review artifact remains a documented provenance defect; nothing here reconstructs it.
