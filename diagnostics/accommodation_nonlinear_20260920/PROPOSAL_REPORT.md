# MACE proposals capture accommodation cheaply; accuracy gain is unproven

**All 60 proposals and 120 solvent endpoints completed.** Actual composite energy selected 49 proposals and retained 11 origins. The method reproduces the expensive four-context candidate contrasts closely, but one canonical reference becomes inconclusive under the unchanged bands. No threshold or production default changed.

## Known-reference discrimination

Canonical references: **24 correct, 0 wrong, 1 inconclusive, 0 unavailable out of 25**, compared with 25 correct at the original geometry. Consumed crystals remain **3/3 correct**. Q9Z4J7 is the only reference decision transition: Ca-supported → inconclusive after a +2.164030 model-kcal/mol change.

All canonical Ca/La raw scores remain separated. The gap (minimum La score minus maximum Ca score) grows **5.076366 → 5.844397 model-kcal/mol**. Within-class ranges change from 28.817426 to 28.388790 for Ca and 32.272923 to 29.684206 for La. These are descriptive changes on consumed calibration proteins, not a new accuracy result. The old band boundary no longer covers Q9Z4J7; no new boundary is fitted here.

| Case | Expected | Old call → selected call | ΔR | Old margin → selected margin |
|---|---|---|---:|---:|
| a0a3f2yly8-pqq-la_model | La | La-supported → La-supported | 3.941901 | 0.000000 → 3.941901 |
| a0acd6b9f2-pqq-la_model | La | La-supported → La-supported | -0.702461 | 20.370509 → 19.668049 |
| c5atj3-pqq-la_model | La | La-supported → La-supported | 0.343343 | 32.272923 → 32.616266 |
| c5axv8-pqq-la_model | La | La-supported → La-supported | -0.401375 | 3.333435 → 2.932061 |
| c5b120-pqq-la_model | La | La-supported → La-supported | 0.885534 | 14.091642 → 14.977176 |
| i0jwn7-pqq-la_model | La | La-supported → La-supported | -1.168882 | 25.508286 → 24.339404 |
| mmol_1770-pqq-la_model | La | La-supported → La-supported | 0.602835 | 11.041624 → 11.644460 |
| mmol_2048-pqq-la_model | La | La-supported → La-supported | -1.554598 | 18.114168 → 16.559569 |
| q88jh0-pqq-la_model | La | La-supported → La-supported | 0.000000 | 4.832205 → 4.832205 |
| q89gy2-pqq-la_model | La | La-supported → La-supported | -0.667197 | 17.002837 → 16.335640 |
| q92wy9-pqq-la_model | La | La-supported → La-supported | 1.180770 | 23.276991 → 24.457761 |
| a8r3s4-pqq-la_model | Ca | Ca-supported → Ca-supported | 0.561240 | 13.616007 → 13.054766 |
| atq70401.1-pqq-la_model | Ca | Ca-supported → Ca-supported | 3.855640 | 28.147925 → 24.292284 |
| bbl57595.1-pqq-la_model | Ca | Ca-supported → Ca-supported | 4.573229 | 28.148320 → 23.575091 |
| o24759-pqq-la_model | Ca | Ca-supported → Ca-supported | 3.408366 | 17.569584 → 14.161218 |
| p12293-pqq-la_model | Ca | Ca-supported → Ca-supported | 2.982444 | 28.458808 → 25.476364 |
| p15279-pqq-la_model | Ca | Ca-supported → Ca-supported | 3.038727 | 27.963248 → 24.924522 |
| p16027-pqq-la_model | Ca | Ca-supported → Ca-supported | 2.592666 | 28.817426 → 26.224760 |
| p38539-pqq-la_model | Ca | Ca-supported → Ca-supported | -1.255261 | 9.870018 → 11.125280 |
| q4w6g0-pqq-la_model | Ca | Ca-supported → Ca-supported | 0.705003 | 10.179234 → 9.474232 |
| q60ar6-pqq-la_model | Ca | Ca-supported → Ca-supported | 2.778788 | 15.787933 → 13.009145 |
| q88jh5-pqq-la_model | Ca | Ca-supported → Ca-supported | 2.271595 | 4.922303 → 2.650708 |
| q8gr64-pqq-la_model | Ca | Ca-supported → Ca-supported | 6.275438 | 25.498947 → 19.223509 |
| q9l935-pqq-la_model | Ca | Ca-supported → Ca-supported | -0.138109 | 25.566774 → 25.704884 |
| q9z4j7-pqq-la_model | Ca | Ca-supported → inconclusive | 2.164030 | -0.000000 → -2.164030 |
| 1H4I | Ca | Ca-supported → Ca-supported | -0.119405 | 42.460234 → 42.579638 |
| 4MAE | La | La-supported → La-supported | -0.606237 | 17.444415 → 16.838179 |
| 1KB0 | Ca | Ca-supported → Ca-supported | -0.369260 | 17.473764 → 17.843024 |

Margins are signed distances to the expected class's frozen support boundary, in model kcal/mol; positive means the correct supported side. Raw R and unrounded margins are preserved in [the CSV](PROPOSAL_REFERENCE_TRANSITIONS.csv) and [report data](PROPOSAL_REPORT_DATA.json).

## Unknown PLM predictions

| Source (original sample 0) | Original R | Selected R | ΔR | Call transition |
|---|---:|---:|---:|---|
| PQQSEQ_83440678cbbd658047c9 | -405530.256481 | -405469.310512 | 60.945969 | Ca-supported → Ca-supported |
| PQQSEQ_07ab500e3df76b30d71c | -405497.955553 | -405460.192016 | 37.763537 | Ca-supported → inconclusive |

Both move substantially toward La-like electronic compatibility; neither becomes La-supported. Their biological labels remain unknown. Native DFT validated earlier bounded displacement responses, **not these final proposals or an affinity classification**.

## Comparison with the full composite search

The separate full search used Tight native GFN2 energies/analytic gradients at every step. This proposal experiment uses vacuum MACE during search and the released primary native GFN2 policy only at q0 and the proposal. The comparisons below therefore include both search and numerical-policy differences.

| Context | Proposal selected ΔR | Full-search candidate ΔR | Tight − primary q0 R | Selected primary R − full Tight candidate R |
|---|---:|---:|---:|---:|
| 1H4I | -0.119405 | -0.080880 | 0.001066 | -0.039591 |
| 4MAE | -0.606237 | -0.558278 | -0.025545 | -0.022414 |
| PQQSEQ_83440678cbbd658047c9 | 60.945969 | 60.931457 | 0.000501 | 0.014011 |
| PQQSEQ_07ab500e3df76b30d71c | 37.763537 | 37.540529 | 0.081121 | 0.141887 |

Final candidate contrasts differ by at most **0.141887 model-kcal/mol** on these four consumed development contexts. This is useful evidence that cheap proposals capture the observed accommodation without solvent evaluations at every optimizer step. It does not establish an exact composite minimum: the full search qualified zero of four complete pairs under its unchanged curvature requirements, and this two-geometry method makes no minimum or entropy claim.

## Execution, physical checks and cost

- GPU job **1204169**: 60/60 successful bounded optimizations, 283 actual native MACE energy/force calls, zero failed requests. All fresh q0 MACE energies reproduce archived values exactly. Every proposal is at least 0.02 radian inside the declared ±0.8 bounds; maximum source-heavy displacement is 0.758463 Å. Source/cap geometry, state, active-group and fixed-atom checks passed.
- Solvent job **1204171**: 120/120 actual primary native GFN2 singlepoints converge, both media/state audits pass, zero failed execution receipts. All 120 old q0 quantum outputs are reused from checked matching artifacts. No fresh q0, DFT or extra poses ran. There are no failed scientific attempts hidden by a fallback.
- Scheduler allocation: GPU stage 95 s × 32 CPUs = **3,040 core-seconds and 95 GPU-seconds**; CPU solvent stage 308 s × 64 CPUs = **19,712 core-seconds**. Total **22,752 allocated core-seconds and 95 GPU-seconds**. This includes scheduler allocation overhead; the executor intervals were 86.678916 s and 290.797782 s.
- Incremental batch average: **758.4 allocated core-seconds and 3.167 GPU-seconds per context**, with pre-existing geometry and q0 quantum results. This is not full source-to-score cost or a matched-hardware comparison with the old DFT method. Peak GPU allocation was 6,051,946,496 bytes; worker peak host RSS was 1,879,092 KiB. Scheduler CPU-stage batch MaxRSS was 6,338,132 KiB.
- Six real-artifact preflight tests passed, zero skips, including original reference replay, real-energy selection/sign algebra and missing-member behavior. Actual optimizer and scientific solver execution are reported separately above.

## Artifacts and recommendation

Implementation/checkpoint commit: `d3b69dd`. Full immutable actual result: `workspaces/accommodation_nonlinear_20260920/proposals_v1/final_1204171.json` and `.md`. Endpoint directories retain actual positions, forces, derivatives, traces and selection components; all submission, execution and scheduler receipts are preserved. The full-search comparison pins `collection_final_v1.json`. [Commands](PROPOSAL_COMMANDS.md) reproduce collection without rerunning chemistry.

**Pursue unchanged-rule structural robustness testing; retain the current production scorer.** The proposal strategy is affordable and captures the observed accommodation, with no class-directed movement. It has not improved demonstrated accuracy: the frozen-band canonical result loses one decisive call while raw class separation remains intact. Test all declared noncanonical reference folds with a separately recorded calibration policy, preserving old-band results alongside it; do not fit to those fold outcomes.
