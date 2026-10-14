# Target La pair: large response reproduced, state caveat retained

Read-only review of completed native jobs 4109/4110 and saved comparison 4111. **No new blocker to the six remaining declared vacuum endpoints was found.** Both La endpoints pass all 17 existing native checks. This extends [the origin review](ORIGIN_STATE_REVIEW.md); no calculation, input change, acceptance change or submission was performed.

## What the pair establishes

The actual Slurm comparison records La origin→adaptive-La work of **−26.70294112086986 kcal/mol** with r2SCAN-3c vacuum, versus **−26.782803599660028 kcal/mol** with MACE vacuum. Its stored DFT-minus-MACE difference is **+0.07986247879016872 kcal/mol**. Thus the large La energy decrease on this exact saved motion is reproduced by the independent vacuum treatment. This is a finite coordinate response, not a DFT optimization or an affinity.

Only two of the declared eight reference cells are available in this comparison. The Ca response, target-minus-control differential and operational shared-pool score remain unavailable. Agreement for La alone does not establish the metal-selective response or classifier benefit.

## Proposal ending and state bookkeeping

Line references here concern `reference_v1/tasks/sample3_La_adaptive_La_vacuum/endpoint.out` under the campaign workspace.

- The last NR/TRAH macro has residual **5.466705e−7**, below the printed **1e−5** criterion (1892; thresholds 1646/1690). SCF convergence after 108 cycles appears at 1896; normal termination at 6326. Final scalar **−6073.139251170734 Eh** at 6121 includes native D4/gCP, rather than substituting the earlier SCF total.
- The output retains molecular charge **−2**, multiplicity **1**, **828 explicit electrons** (1625–1627), and La ECP46 (1457). Equal alpha/beta integrals of 414.000153787722 and total 828.000307575443 appear at 1921–1923; Mulliken charges sum to −2 at 2569. These agree with the declared origin bookkeeping. They do not measure the physical cofactor oxidation state.
- The last occupied orbital is printed at **−0.010045 Eh**, occupation 2; the first unoccupied orbital at **+0.089955 Eh**, occupation 0 (2356–2357). The resulting printed 0.100000 Eh separation has the same level-shift caveat as the origin: a 1 Eh shift was explicitly enabled at 1744, a negative Hessian diagonal and gap −0.000144 appear just before the final macro (1880–1881), and a small negative-gap warning remains after convergence (1899–1900). This table does not establish an unshifted physical gap, validated Aufbau occupation or stable unique ground state.
- Large density/DIIS summary values at 1933–1935 remain recorded. Admission uses the actual ending TRAH criterion; this is not an all-criteria density-convergence or repeatability qualification. No independent stability/occupation/redox check exists for either endpoint.

The agreement therefore supports continuing the unchanged finite diagnostic with the electronic-state limitation explicit. It does not resolve that limitation or justify an automatic solver change. Native proposal runtime was **1 h 48 min 7.423 s** (6327); root retains scheduler and aggregate cost accounting.

## Immutable evidence

Workspace prefix: `workspaces/exaa_donor_response_20261014/`.

| Record | SHA256 |
|---|---|
| `scout_comparison_v1/execution/RESULT.json` | `233b8c100791215f11564dc9ae053c5a897b5a4c4ee6bb6c67adfa77f9dcfb7d` |
| `reference_v1/attempts/4110_execute_sample3_La_adaptive_La_vacuum/COLLECTION.json` | `42d54e3c79fd9c9f42b28b518b5cf89e68c4a704eb2e8fd5cb647c45e7b96c8b` |
| Proposal `endpoint.out` | `e93f387f5adf1c2acfb062ff2b58ead8538dbcc2539408186558373310916e9a` |
| Proposal `endpoint.out.execution.json` | `c9896770aab911ba8ff56350e958a338d765fe8e63c37e197968e847aa6443bf` |

Origin pins are preserved in `ORIGIN_STATE_REVIEW.md`. The comparison retains six missing cells explicitly; no substitute energies or login-side scientific algebra were used for this review.
