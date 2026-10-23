# Matched donor counterpoise comparison

**PBE0-D4 and canonical CCSD(T) agree on both exchange signs, but their agreement depends on the donor.** The acetamide paired contrast differs by 1.40 kcal/mol. Bidentate acetate differs by 0.255 kcal/mol, narrowly outside the declared 0.25 resolution. This energy comparison does not qualify either method as physical truth or establish force accuracy.

All values below are kcal/mol. The paired quantity is **Dy minus La** for replacing bound water with the specified donor, using the same fixed coordinates and counterpoise (CP) correction to account for fragment basis sharing. Positive values mean this fixed electronic exchange is less favorable for Dy relative to La; negative values mean more favorable. These are not aqueous binding free energies.

| Donor | CCSD(T), CP | PBE0 electronic, CP | PBE0-D4, CP | PBE0-D4 minus CC |
| --- | ---: | ---: | ---: | ---: |
| Acetamide | +6.070361 | +4.333381 | +4.669355 | −1.401006 |
| Bidentate acetate | −5.298568 | −6.103436 | −5.553463 | −0.254895 |

## What this resolves

The retained water-only PBE0-D4-minus-CC paired gap is +0.467767. After balanced water-to-donor exchange, its magnitude grows to 1.401006 for amide and falls to 0.254895 for acetate, with both signed gaps negative. Thus the water method discrepancy does not supply a universal correction for donor exchange. The strict 0.25 test fails for both donors; acetate exceeds it by only 0.004895 and should be described as a close comparison, not a substantial contradiction.

D4 shifts the paired exchange by +0.335974 for amide and +0.549974 for acetate. It improves agreement with CC here, but the electronic-only and D4-added results remain separately reported. D4 was evaluated only for real molecules, with the real free-ligand contribution subtracted in each interaction. CC received no D4.

The recipes retain matched bases, frozen-f ECPs and effective singlet RHF states. Physical Dy spin, spin-orbit effects, basis convergence and the acetate diffuse-basis sensitivity remain unqualified. No geometry optimization, force test or parameter fit is implied. The separately authorized four-parameter water representability scout retains its original scope and criteria.

## Admission and evidence

All 22 new native cells completed: 12 CC and 10 PBE0 fragment calculations. Four existing PBE0 complexes and the existing water comparison were reused. The strict saved-output reviewer passed **510 checks**, binding each result to its actual cell, input, source and native energy. Independent decimal cycle arithmetic agrees with the original collector within 2.21×10⁻¹¹ kcal/mol.

All CC final residuals satisfy the frozen 10⁻⁸ criterion; the largest is 8×10⁻⁹. T1 values range from 0.001373326 to 0.023312701, with the maximum in Dy–acetamide. Raw amplitudes, triples contributions and SCF stop signals remain preserved; no new diagnostic cutoff or physical-accuracy claim is introduced.

Exact paths, all 22 result pins, scheduler receipts, source/helper/basis lineage and per-cell costs are in [FINAL_EVIDENCE.json](FINAL_EVIDENCE.json), SHA256 `3d442f16c25af5c488047899e0edfc9502ce1200d8ef7853720eef9562b04821`.

| Artifact | SHA256 |
| --- | --- |
| Frozen panel | `1a4bca36f2e5f0d0e223b235b0915708a455627444bcb355fc7cbcba2930a701` |
| Native adapter / original collector | `f4baaeddce09cb26551e845f05e7f3c1bfff604a84fc02c666908ada2420dea2` |
| Strict review source | `eef93c84d9a469dd3a15fbe2ed2213140b0e4a51f2c6ddacbbfc4299c76cd379` |
| Full strict review result | `0fdcf05d079a2e846b8d2dccc7ac8d8a764fea4c88c03c88659027a74f8131a1` |
| Full collection result | `c9ca7b157973c513fc808f3caf6976463da7863ebcdaeb2fc19f296fd9b74348` |
| Reused water comparison | `6465f8c9002c31ee9c79a945152476c91ad740da0a5ddca4e2327d98321bda64` |

## Actual incremental cost

All jobs used four CPU ranks and no GPU. Scheduler allocation cost is **60,572 CPU-seconds (16.83 CPU-hours)**: 57,948 for the 12 CC jobs and 2,624 for the 10 PBE0 jobs. Native elapsed time multiplied by ranks totals 60,403.884 CPU-seconds; these are allocation measures, not measured CPU utilization. Previously completed complex and water calculations are excluded from this incremental cost.

The largest scout ran 2,303 scheduler seconds; its measured TotalCPU was 02:25:17 and peak RSS was 16,383,484 KiB, close to the 16,000 MiB grant. It terminated normally with empty stderr. The strict local review took 6.78 seconds on one thread. Collection and this report made zero native calls or submissions.
