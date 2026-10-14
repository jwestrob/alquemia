# Target response: independent vacuum DFT supports the large La-directed motion

Read-only review of completed target comparison 4125 and its four endpoint records. **The large metal-dependent energy response is reproduced by r2SCAN-3c on the exact saved motion.** This narrows the case for a MACE-specific local-energy error, while retaining the reference's numerical and chemical-state limitations. No scientific calculation, submission, criterion change or result replacement was performed for this review.

## Sign and scope

The comparator uses each metal's own origin: `W_M = E_M(adaptive_La) − E_M(origin)`, followed by `D = W_Ca − W_La`. Its saved joins retain the same source, geometry path, composition and method within each work. Different metals' absolute DFT totals are not interpreted as affinities.

| Target sample 3, kcal/mol | MACE vacuum | DFT vacuum | DFT − MACE |
|---|---:|---:|---:|
| Ca work | −10.379323 | −9.567475 | +0.811848 |
| La work | −26.782804 | −26.702941 | +0.079862 |
| Ca − La work | +16.403481 | +17.135466 | +0.731985 |

Both metals' energies fall along this path; La's falls more. Positive `D` therefore increases the Ca-minus-La contrast and is **La-directed**, not Ca-directed. The stored differential discrepancy is below the predeclared 1 kcal/mol qualitative resolution floor. That floor is neither a demonstrated uncertainty bound nor a newly imposed acceptance threshold.

These are matched **vacuum** responses. They do not qualify the ALPB correction or constitute a DFT classifier. Both metals traverse the La-generated proposal; the operational scorer instead selects its own lowest admitted common-pool candidate for each metal. This diagnostic therefore must not be substituted for its full composite/minimum-selection response. The four-cell snapshot explicitly leaves the control works and target-minus-control DFT contrast unavailable; it does not imply that every separately finishing control job remains incomplete.

## Numerical and electronic qualification

- Ca origin 4117 ends in TRAH/NR with residual **9.438462e−7**, below its printed **1e−5** threshold. The raw final macro and SCF convergence are at lines 1839/1843 of `sample3_Ca_origin_vacuum/endpoint.out`; all 17 recorded checks pass. Charge −3, singlet and 838 explicit electrons are printed at 1603–1605. Two negative Hessian diagonals and a negative gap remain near the last macro (1830–1831), followed by a small negative-gap warning after convergence (1846). Its 1 Eh runtime level shift is recorded at 1716; numerical admission is not state uniqueness.
- Ca proposal 4118 has 16 passing checks and an **SOSCF energy-convergence** ending (1789–1793). Its final gradient **1.0425e−5** exceeds **1e−5**, and its final rotation also exceeds the printed tolerance (1831–1832). Earlier iterations with smaller gradients do not replace this final result. The existing collector correctly follows the frozen solver-specific rule, but this endpoint is not qualified to the scouts' extra TRAH gradient standard. The excess residual alone supplies no numerical energy-error bound.
- The La pair passes the actual final TRAH gates, as documented in [PAIR_STATE_REVIEW.md](PAIR_STATE_REVIEW.md). Its late gap/level-shift caveats remain. Large density/DIIS summaries are retained throughout. None of these calculations establishes a unique stable electronic state, physical PQQ oxidation state, repeatability or analytic-force accuracy.

## What this narrows, and the next decision

The earlier expanded packing calculation changed the target's operational composite differential only from **+16.691741 to +16.316953 kcal/mol**, while retaining its choices. The admitted crystal chain had a much smaller composite response (**+1.929997 kcal/mol**), under the same conditional dry, oxidized-PQQ local-state convention. Those experiments and the new matched vacuum reference concern different response definitions; their numerical values must not be treated as a single interchangeable series.

Together, they weaken two specific explanations: the selected omitted local packing alone does not remove the modeled response, and a MACE-only exaggeration of this target's local vacuum response is not supported by the independent endpoint comparison. They do **not** prove that the fold is wrong, that the crystal is the uniquely relevant state, or that the declared isolated dry/cofactor state faithfully represents the protein. Shared local-model assumptions, source-dependent strain and omitted environmental/state effects remain possible.

**Next decision:** finish the already-declared control cells and inspect the actual target-minus-control differential with each endpoint's convergence limitations visible. If DFT also reproduces the smaller control response, the evidence favors investigating why this source/state carries unusually large releasable strain, rather than suppressing the MACE response or refitting bands to it. A control disagreement would instead identify a narrower transfer/reference-state question. Neither outcome justifies an automatic new solvent, optimization or functional campaign before the finite matrix is interpreted.

## Evidence pins

Authoritative snapshot: `workspaces/exaa_donor_response_20261014/target_comparison_v1/execution/RESULT.json`, SHA256 **`f17068d7c0432960eecd4c4b6282c7be35d0e625c153afb9d448b986aeaf437a`**. It pins all four collections, raw outputs, receipts and the source manifest. Its comparator is pinned as **`cf41923fe324652e375debd4f403e236318fa9bb13ba87636cbedbed166b4a69`**; work/sign formulas are at lines 53–78.

Ca origin raw output SHA256: `31f1d326645747857450719a0a841a27a9608cffede184186b53544926c6564e`; Ca proposal: `974c01bc8891a5511985cbdb9fcba708a760ab68f632f5a4a75ed3c38040b933`. Raw line references above use their existing `reference_v1/tasks/` directories. Prior evidence is retained in `packing/FINAL_REPORT.md` and `crystal_scoring/RESULT_REVIEW.md` under `diagnostics/exaa_source_packing_20261014/`.
