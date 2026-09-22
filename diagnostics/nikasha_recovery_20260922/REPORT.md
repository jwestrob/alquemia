# Keep the original composite: donor proposals did not improve discrimination

**The original context composite is the strongest decisive single-fold method on this consumed reference challenge.** It gives 203 correct calls, two wrong calls and two inconclusives among 207 available sources. Native MACE core gives fewer decisive calls but zero wrong calls. Preserved DFT, native core and the original context composite tie on all 94 available three-fold summaries. The new accommodation proposal method does not improve this result, and its separate canonical-only bands worsen transfer.

This is a structural robustness test on 25 already-consumed protein groups, not 225 independent biological controls or prospective validation. The 25 canonical geometries remain calibration-only. All 225 noncanonical sources remain in denominators; no thresholds, samples or production settings changed during recovery.

## Completed fixed comparisons

Cells are **correct / wrong / inconclusive / unavailable**. The La-conditioned arm is the main scanner use case.

| Method | All225 | La100 | Ca125 | All100 La triples |
|---|---:|---:|---:|---:|
| Preserved DFT | 199/3/6/17 | 96/1/1/2 | 103/2/5/15 | 94/0/0/6 |
| Native MACE core | 199/0/9/17 | 95/0/3/2 | 104/0/6/15 | 94/0/0/6 |
| Native MACE context | 201/6/1/17 | 93/4/1/2 | 108/2/0/15 | 90/2/2/6 |
| Core composite | 189/0/18/18 | 89/0/8/3 | 100/0/10/15 | 87/0/4/9 |
| Original context composite | 203/2/2/18 | 97/0/1/2 | 106/2/1/16 | 94/0/0/6 |
| Proposal, old bands | 201/1/3/20 | 97/0/1/2 | 104/1/2/18 | 94/0/0/6 |
| Proposal, new fixed bands | 197/2/6/20 | 95/0/3/2 | 102/2/3/18 | 92/0/2/6 |

The four triples per protein are correlated subsets, with every member required. DFT, native core and original context composite each have 23/25 proteins with all four triples correct; two groups contain unavailable preparations.

## Common coverage and actual transitions

- all225: on the same 207 available sources, DFT 198/3/6/0; original context composite 203/2/2/0.
- La: on the same 98 available sources, DFT 96/1/1/0; original context composite 97/0/1/0.
- Ca: on the same 109 available sources, DFT 102/2/5/0; original context composite 106/2/1/0.

On La-conditioned sources, the composite makes one additional correct call and removes the one DFT wrong call. It does not establish broader La/Ca affinity accuracy. Raw per-case transitions and all matched subsets are in summary.json and proposal_comparison.json.

| Case | Condition | DFT → original context composite |
|---|---|---|
| a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-0 | Ca | inconclusive → La-supported |
| a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1 | Ca | inconclusive → Ca-supported |
| a0a3f2yly8-pqq-la_model__conditioned_La__seed-1_sample-0 | La | inconclusive → La-supported |
| a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-4 | La | Ca-supported → inconclusive |
| mmol_1770-pqq-la_model__conditioned_Ca__seed-1_sample-1 | Ca | Ca-supported → La-supported |
| q4w6g0-pqq-la_model__conditioned_Ca__seed-1_sample-2 | Ca | inconclusive → Ca-supported |
| q4w6g0-pqq-la_model__conditioned_Ca__seed-1_sample-3 | Ca | inconclusive → Ca-supported |

## Accommodation utility

All 415 eligible optimizers succeeded; 2,046 native MACE energy/force calls completed. There were 822 fresh proposal GFN2 attempts: 820 complete and two failed. Eight identical proposal=q0 endpoints reused the exact prior result. Selection chose 345 proposals and 68 origins; three metal endpoints remained unavailable. No failed proposal fell back to a baseline success.

The proposal protocol has 205 available site pairs. Its old-band result is 201/1/3; the independently frozen new-band result is 197/2/6. Both retain 20 unavailable cases. Changing calibration is therefore material and was not a solution to transfer accuracy.

- La4, same 23 complete protein groups: median within-protein raw-score range 9.083191 → 10.298768 model kcal/mol; range decreased in 13 groups and increased in 10.
- Ca5, same 19 complete protein groups: median within-protein raw-score range 12.128668 → 10.026598 model kcal/mol; range decreased in 9 groups and increased in 10.

The exact missing proposal pairs are the pre-existing A8R3S4 Ca-conditioned sample3 solvent failure, plus new failures for P12293 Ca-conditioned sample4 (La vacuum) and Q60AR6 Ca-conditioned sample1 (La ALPB). The original 17 preparation failures remain separate. Lower energy or narrower geometric spread does not by itself establish better classification.

## Cost and independent electronic reference

- Preserved DFT: 416/416 new endpoints completed; **2,802,688 allocated core-seconds**, no GPU. Archived25 pairs reused.
- Proposal continuation: **200,032 allocated core-seconds and 633 GPU-seconds**, including the two failed solvent attempts and allocation/collection overhead. This is incremental to archived q0/context preparation and not full source-to-score cost or a matched-hardware speed comparison.
- Independent CC job1202429: **12,343,360 allocated core-seconds**, no GPU; only3/6 endpoints complete. All three Ca outputs show MDCI error termination and MPI segmentation faults. No complete Ca/La pair exists, so it cannot adjudicate whether the electronic treatment improves classification. It is too costly for routine use. The failure cause was not isolated by this recovery.
- Original preparation and parser/collection CPU were not separately metered. Recovery itself made zero molecular calls; scheduler costs above are actual Sept20–22 receipts, not estimates.

## Recommendation and reusable artifacts

**Retain the original context composite for this development comparison; retain the production baseline/default.** Reject promotion of the current proposal-selection method. Three-fold aggregation is robust on the available consumed references, with no gain over preserved DFT. The next experiment must provide a transferable accuracy benefit rather than a favorable energy shift or recalibrated canonical separation.

Final reusable proposal selection: `workspaces/nikasha_recovery_20260922/final_selection.json`. Full fixed join: `proposal_comparison.json`; full DFT/MACE join: `final_DFT_MACE_comparison.json`; exact100-triple DFT comparison: `final_DFT_triples/result.json`; compact machine-readable summary and transitions: `summary.json`. All are under that recovery workspace and pin the original immutable runs. Parent separately owns the shared-candidate-pool experiment; no new chemistry belongs to this recovery.
