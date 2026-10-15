# A0AC alternative-amide native-reference feasibility

**No exact reusable native DFT vacuum endpoint was found for either current Ca4/La4 origin or either relaxed-flip geometry. No new reference calculation is prepared or submitted.** This is a bounded metadata/hash audit, extending the [earlier Ca4 cache check](../../a0ac_asn_inventory/CACHE_NOTE.md); it is not an exhaustive geometric-equivalence search through historical outputs.

## Exact inputs and cache distinction

The authoritative [completed MACE collection](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_remaining_v1/execution/FINAL_COLLECTION.json), SHA256 `2ed05f5a241f4f1eccbb2a680428a4b032a51449df56ad323e138b83cf32fd82`, supplies both source rows, all original/relaxed XYZ and MACE receipts, final physical admission and proposal pins. Sources are `a0acd6b9f2-pqq-la_model__conditioned_{Ca,La}__seed-1_sample-4`, selection `a0acd6b9f2-pqq-la_model__envelope_2eda778295ffed58`.

Both retain **204 atoms, 14 caps, dry oxidized PQQ3−, Ca charge −4 / La charge −3, singlets**. Original H geometry is retained; the separate two-H projection is not substituted. Exact state signature SHA256 is `92cf7c52d222f9b925eb30346ba5ea701e1b5fa5d3197f70d78299227d20707a` for both. [Preparation v2](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_v2/preparation/PREPARATION.json), SHA256 `6405de39b9d3ffa88aea7f99f895a685ddcd8ec5b1b37529766b72c365d9d57d`, retains each original source/map/state/protonation pin and its distinct rebased map. Original Ca4/La4 map SHAs are `eddbfcd00812399be7f9fbec4c3a642c4b2ae5bcc2da88f0a97dd0415b05182c` / `74a8bdb08f8fdc74bab773d9d5bc85c6116a94a05877d0f44ef91d5b436df0cf`.

The minimal shared-path inputs are these exact collection `matrix[metal][candidate].xyz` pins:

| Source | Metal | Origin SHA256 | relaxed_flip_La SHA256 |
|---|---|---|---|
| Ca4 | Ca | `fc080c643e09252b95dac53b73ae733d2d462c3851c647bd620588c2831b7179` | `2bc7be4df3ecc7f4d5fb5945617831764f46818902478a7e81686e69c6200e4a` |
| Ca4 | La | `a9704aabef3cb2cd389a0a6651282e99855854c16600b244335146b7439fb45c` | `bbf2aa9a03d1ca00f8e4d0bf7f3934763f76a167103d638adca506e3e29a5aa4` |
| La4 | Ca | `ed1569dfe0cd70ba12dbd4708914059cab75be033a77d09968fac3310700ff53` | `6338947b745636f7ee1c96040ea759c7296475ec3bd1503fd662b981543743c3` |
| La4 | La | `a5a12c9650f1e13bb5f980bd9f72a45c96182aa13655fb673dfa427817b7c13e` | `c0898baeb9b16dcd5e80b7f8c62462b20d5730d27b7bb4c8d4737cbadf58a897` |

All 16 current origin/fixed-flip/two-relaxed XYZ hashes are absent from the targeted fold-DFT root/four shard manifests, environment-PQQ prepared_v2, reference-profile prepared_v1/low, ExaA reference_v1 and strict-donor-response run_v1 manifests. The same Ca4 **and** La4 source IDs do occur in [fold-DFT shard1](../../../../workspaces/accommodation_fold_DFT_20260920/run_v1/shard_1/manifest.json), but those are **80-atom CPCM inputs, charges −3/−2**, with different XYZ hashes. The older environment context is a different 161-atom preparation; 180-atom basin/scaffold evidence is also incompatible. Recent MMOL/ExaA native references concern different proteins. No compatible wavefunction was identified through these references.

## Smallest informative reference, if later declared

**Recommendation: finish uniform transfer first; if its outcome warrants independent local-energy validation, start with four exact Ca4 endpoints only.** The same-composition/state energy pair is origin → saved `relaxed_flip_La` **within each metal**. Four energies are necessary to obtain both works and their Ca−La differential; two energies alone test only one metal. Identical coordinates are used for the two exchanged metals at each geometry. No optimization, H repair, context expansion or solvent belongs in that comparison.

This would materially distinguish an independently supported large response from a MACE-specific surface artifact **along that actual saved path**. Both metals are evaluated on the La-generated geometry, avoiding a mixed-path subtraction. It tests the full origin-to-alternative-basin vacuum response; it does not independently validate the incremental common3→6 pool gain or all selected paths. Agreement supports the energy component of the missing-basin explanation, without proving that basin is physically populated. Strong disagreement would weaken that explanation before further reference work. The analogous four La4 cells are a **conditional control expansion**, only if the Ca4 outcome is interpretable and leaves a specific source-comparison question; eight cells are not automatically required. Cross-scoring both relaxed geometries plus origin would require six cells per source and is outside this minimal recommendation.

Reuse the native **ORCA 6.1.1 r2SCAN-3c / DefGrid3 / TightSCF / NoAutostart vacuum** executor and counted collector from [MMOL reference](../../mmol_reference/REFERENCE.py), with actual basis/ECP/electron/convergence checks and preserved energy-only-ending/state warnings. Its current adapter explicitly binds MMOL193/13 and charges −3/−2: it is **not directly runnable** on these 204/14 inputs. A separately pinned source adapter and Slurm preflight would be required; no existing map should be relabeled to satisfy it.

Actual comparable 12-MPI MMOL193 costs were 3683–5497 wall seconds per endpoint; jobs 4603/4604/4637/4638/4639/4640 used respectively 47736/45468/49116/55020/65964/44196 allocated CPU-seconds. [Their final report](../../mmol_reference/FINAL_REPORT.md) retains endpoint identities and receipt links. The [ExaA192 eight-cell reference](../../../exaa_donor_response_20261014/FINAL_REPORT.md) consumed 522696 allocated CPU-seconds, with individual times 3670–6994 seconds. These are cost precedents, not a forecast for A0AC or measured process utilization.

## Execution boundary and interpretation

The completed relaxation gate released matched frozen-charge solvent only; its [plan](../solvent/PLAN.md) explicitly excludes DFT. The [uniform transfer plan](../transfer/PLAN.md) likewise introduces no DFT. Neither canonical success nor existing solvent support silently releases a native-reference branch. Root's current instruction is to finish uniform transfer before deciding whether this independent validation earns a separately declared run. **This note is feasibility only, not a new qualification dependency for the transfer result.**

Agreement would support MACE's local vacuum energy response for these exact finite states; disagreement would motivate a local energy-surface/electronic-state diagnosis. Neither outcome validates frozen-qvac solvent, experimental preference, a unique electronic ground state, or universal physical accuracy. All work here was saved-report/source/metadata reading; no scientific calculation or submission occurred.
