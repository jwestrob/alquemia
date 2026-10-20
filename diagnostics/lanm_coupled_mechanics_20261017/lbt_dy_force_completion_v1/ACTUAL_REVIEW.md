# Independent actual-result review

**Job 15720 completes the declared force-coverage test: eight available endpoints and four available matched-metal comparisons. No admission or provenance blocker was found.** The result shows substantial cancellation of model/reference force errors between metals at fixed geometry, with a remaining metal-dependent discrepancy.

## Exact evidence

- RESULT.json: `f42cffcf5717a988ed569db925b840f02dfb8bd3b4c04750b77e33bc231c4790`.
- Executed source: `d77cc40e6cb76f30beb5408c0fafa3d8107c0979f526beb069ee7060ed3a177f`.
- CONFIG.json: `b1655fb8d13b23e036365fc068dfcde2906d83469af693bd94bf526456fb65a6`.
- Inventory: `e2c196d1b5961cc74d50647e0565780928c3c1c103b282d563833f4b7e003957`.
- Completion event at inspection: `240c11f7905bab3356c08260fda9bc7a2b82f03800615bf8a0f7f7ec6c91013b`. It records scheduler COMPLETED, exit 0:0, 12 seconds / 12 CPUs on the qualified host; notification delivery remains marked queued/unconfirmed in that inspected snapshot. Root owns acknowledgement.

The result and eight force reports are under `workspaces/lanm_coupled_mechanics_20261017/lbt_dy_force_completion_v1/run/execution/`. I verified 58 selected artifact pins, plus both parent collectors, and exact saved metadata joins. All eight report sources/native rows match the frozen inventory; atom IDs and 256-row coverage match. Each report directly names its actual raw MACE receipt and task.

| New supplemental raw MACE receipt | SHA256 |
| --- | --- |
| 7CCO_before_Dy | `9980e5781103caaf5d60d211b36661ae9da36f2a910ce5db6529bfb55c57ada1` |
| 7CCN_before_Dy | `20484956b06efc18cc1a7f471a8867293fe7941007805576a024d2ae7e323526` |
| 7CCN_after_Dy | `9e717cb5e0baa17409b842bc624f318bb0ea504e6bceec146b9f3de798d2a249` |

All three raw energies equal their original saved energies, with recorded replay deltas of 0.0 eV. Native input-state records are identical to their originals: all 256 atoms, charge −2, physical multiplicity 6, target 255, one graph, omol head, no charge-prediction claim. Supplemental tasks differ only by the declared task ID and original-task provenance link. The original three energy-only raw records still have null forces.

The producer records admission of all eight saved cells before model loading, exactly three attempted forwards, one model load and NVIDIA H200. Application elapsed time is 10.067909 seconds; model loading took 6.749189 seconds. Five force arrays are reused. New DFT/ALPB/GFN/MM/geometry calls and optimizer steps are all zero. Original energy comparisons and reference caveats are copied exactly.

## Reported force comparison

Values below are copied from saved metrics, not recomputed. Units are kcal/mol/Å; RMS uses full atomic-vector magnitudes.

Individual MACE-minus-PBE0 composite force-error RMS spans 3.511987–3.609827 across eight endpoints. The error in the matched `F_Dy − F_La` vector is smaller:

| Source / geometry | RMS contrast error | Maximum atom contrast error |
| --- | ---: | ---: |
| 7CCO before | 0.501500 | 2.704859 |
| 7CCO after | 0.476265 | 2.259919 |
| 7CCN before | 0.726763 | 6.242802 |
| 7CCN after | 0.563471 | 4.148825 |

This supports a substantial error component shared between metals at each fixed geometry. It does not locate that component anatomically or show that the uncancelled component is harmless. Global RMS cancellation also does not preserve every directional sign: the original 7CCO-after chord has composite Dy-minus-La projections of −0.044940 for MACE and +0.054644 for PBE0. No new force-accuracy threshold was declared.

The original energy/work conclusions are unchanged, including the 7CCN differential work discrepancy. Shared ALPB cancels from model-versus-reference errors and therefore receives no independent validation. Both fixtures are consumed peptide deformations; the 7CCN after geometry is step 80, not the later converged La endpoint.

All eight native rows retain `all_printed_residuals_pass=false`, `physical_state_qualified=false` and `force_consistency_qualified=false`. Physical Dy multiplicity 6 MACE versus frozen-4f9/effective-singlet DFT remains a comparison of distinct approximate electronic representations. Cancellation does not establish a physical state, reliable preference, accurate optimizer, or natural LanM transfer.

## Next bounded interpretation

A saved-only localization by predeclared atom categories can usefully distinguish shared terminal/scaffold errors from residual target/donor discrepancies. Keep all four matched pairs and all atoms, use existing chemical identities rather than error-selected subsets, and describe attribution without fitting or new acceptance gates. This requires no new molecular calls. The separate 15545 reference still addresses stationarity at the converged La geometry.

Review used source text, hashes, exact saved metadata equality and reported metrics only. No vector arithmetic, numerical replay, scientific imports, jobs, source changes or receipt acknowledgement were performed.
