# Independent review of the completed Mex whole-system collection

Reviewed 19 October 2026. **PASS for complete saved collection and provenance. The measured per-metal spin sensitivity prevents treating one aggregate-spin input as a unique prediction.** No independent Mex accuracy or physical-state qualification follows.

This review read source, saved results and scheduler receipts, checked artifact hashes and array dimensions, and copied the reported quantities below. It performed no numerical recomputation, scientific import, model call or submission.

## Actual result and denominator

Slurm collector 14794 completed successfully in 3 scheduler seconds on one CPU, with zero model loads, forwards, other molecular calls or optimizer steps. Its RESULT is:

`workspaces/lanm_coupled_mechanics_20261017/mex_whole_collection_v1/execution/RESULT.json`

SHA-256 `4c0c8b7379a2210f64267f0df94b1df72b54902d7233ec4b9614d4e941349d27`.

The executed source is `dd611ca03982e4f2548a88dabf9477e33aba7592fc2b41d832e60d13edd2add7`, collection config `1536e02656d4887c8828c11c72f31d696b8fd72b4302337195968f0d01c4d6d7`, and unchanged preflight/collection manifest `acdaa4623a85a6e36fa7edc23d43a94fd94bc51ff9465aea69834a41f1b8d8d0`.

All eight ordered cells are complete: two exact job 12970 geometries × La/Dy × the declared aggregate multiplicities 10/15 or 4/9. Each preserves 2,049 atoms, charge +2 and target index 2045. The actual native records retain the three Nd spectators and 167 waters. The collection source requires the exact admitted scout, seven-forward remaining result, unchanged manifest and exact raw/component/total replay for every cell before reporting work. The scout row is reused, not reevaluated. Four ALPB receipts serve both spin hypotheses.

Parent results are scout `5eb4ca1738a19098de9d4788d812aa31ac610b4b32aedcbcc022f62475a6dd8c` and remaining `4df04ed07607c7340a3cac1b7e6c17db87f74e2de8bc943e926f50996fa43fe5`. Both parent pins, the source/config/manifest, and all eight raw and assembled artifact pins and four unique ALPB result pins were independently verified. No survivor filtering or additional call is present. The original failed preflight and closed derivative experiments remain separate preserved evidence.

## Reported deformation work

Values are after minus before, in kcal/mol, at the same fixed coordinate pair.

| Aggregate conditioning | La total | Dy total | Dy minus La work |
| --- | ---: | ---: | ---: |
| Aligned counting, La 10 / Dy 15 | −1937.990241 | −1942.176256 | −4.186016 |
| One-Nd-reversed counting, La 4 / Dy 9 | −1944.455151 | −1948.078345 | −3.623194 |

The stored change from aligned to one-Nd-reversed counting is **−6.464910 for La, −5.902089 for Dy, and +0.562821 for their work difference**. The original greater-than-1 kcal/mol flags are therefore true for each metal and false for the differential. The smaller differential reflects cancellation in this two-hypothesis comparison. It does not establish equivalence of the hypotheses or remove the material per-metal sensitivity.

Both hypotheses lower the candidate energy along the existing deformation. Vacuum whole-MACE work is strongly favorable, while the unchanged ALPB contribution opposes it: +552.925885 for La and +551.835625 for Dy. The old hybrid comparator reports totals −2945.113073 and −2934.389151. The replacement therefore changes the magnitude and metal ordering of this fixed-deformation response; it does not demonstrate that the remaining favorable work is physically correct. Old hybrid work is a comparator, not an independent reference. Removing the partition changes scaffold, boundary, long-range and spectator-ion treatment together, so these differences do not identify a unique faulty term.

The negative Dy-minus-La work means Dy benefits more from this particular motion under the candidate. It is not a binding free energy, exchange preference or source-robust profiler result. No absolute total energies from different Hamiltonians are compared.

## Force changes and scientific limits

Every stored component and total force-change array has 2,049 rows of three components. The source defines these as **force after minus force before**. Reported total-force changes are:

| Conditioning | Metal | Atom-vector RMS, kcal/mol/Å | Maximum atom-vector norm, kcal/mol/Å |
| --- | --- | ---: | ---: |
| Aligned | La | 22.171159 | 106.813121 |
| Aligned | Dy | 21.679724 | 97.277876 |
| One Nd reversed | La | 22.361701 | 109.871263 |
| One Nd reversed | Dy | 22.123908 | 105.753318 |

These demonstrate a substantial change of forces between the two endpoints. They are neither endpoint convergence norms nor reference force errors, and they are not direct force differences between spin hypotheses. Similar summary norms cannot establish similar force directions or local-spin robustness. A single finite deformation does not measure a Hessian or general scaffold stiffness.

The native feasibility gate passed and all fixed cells are available. Aggregate multiplicity categories still do not specify localized Nd moments, identify a reversed ion, determine exchange coupling or establish a physical electronic state. Selecting whichever hypothesis gives a lower inspected energy would discard the declared sensitivity test. The 1 kcal/mol threshold is reporting resolution, not an error bar or accuracy bound.

The finite panel supports feasibility and quantifies model sensitivity. It does not qualify whole-protein accuracy, authorize a unique-spin optimizer or promote a production model. The separate 7CCN reference challenge can test transfer on that peptide deformation; its pending outcome cannot retroactively establish the mixed-ion Mex state.
