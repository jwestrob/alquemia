# Three saved folds do not change either difficult PLM prediction

**Both preselected PLM targets remain Ca-supported in all three saved folds.**
The median equals the original sample0 score for each. This rules out choosing
one of these three folds as the explanation for the persistent calls; it does
not prove the calls are biologically wrong, because neither protein has a verified
La/Ca label. The practical benefit is a clearer development choice: averaging
these saved structures does not resolve the concern, so the separate physical
response work remains relevant to testing structural/model compatibility.

These are isolated developmental predictions. No production scan, original score,
scoring parameter, threshold or default was changed. Both proteins/all three
seed101 samples were fixed before calculation; no sample was selected by score.

## Actual outcomes

Raw R is E_Ca−E_La for native MACE+[nativeGFN2ALPB−vacuum], in model-kcal/mol.
The unchanged context bands transfer to a three-fold median explicitly as
**developmental**. All six sources retain **unvalidated_input_domain**.

| Target | Sample0 R | Sample1 R | Sample2 R | Median | Spread |
|---|---:|---:|---:|---:|---:|
| PQQSEQ_07ab500e3df76b30d71c | −405497.956 | −405520.765 | −405481.404 | −405497.956 | 39.361 |
| PQQSEQ_83440678cbbd658047c9 | −405530.256 | −405482.212 | −405551.390 | −405530.256 | 69.177 |

Every individual call and both medians are Ca-supported. The large spreads are
retained, not interpreted as a calibrated uncertainty or hidden behind the median.
Archived source contacts give the extra-Asp nearest O–La distances across
samples0/1/2 as1.847/1.735/2.120Å for07ab and1.766/2.055/1.700Å for8344. None of these
calculations moved a donor or optimized the structure. This does not establish a
unique cause for the scores or guarantee that physical response will improve them.

All six sources passed the current frozen preparation; **07ab sample2's prior
single-fold admission failure remains explicit** in every exported member record.
The historical admission and current chemistry-preparation gate are different
checks. Passing the latter does not retroactively reverse the former. No failed
member was removed or replaced, and no recipe was patched to obtain a result.

## Exact sample0 comparison

Both fresh sample0 contexts exactly match the previous torsion experiment's q0
contexts, including every coordinate, charge and multiplicity. Native checkpoint,
ORCA executable and actual GFN2 input text also match. Score differences from
those genuine archived origins are0 and6.40×10⁻¹⁰ model-kcal/mol. Individual MACE
and GFN2 differences are retained in [RESULT.json](RESULT.json). No unlike
preparations were forced into a repeat comparison, and no old energies were
substituted for the fresh standard execution.

## Execution and checks

Jobs1204068/1204069 completed12nativeMACE+24nativeGFN2 calls, zeroDFT/folding,
with six fresh source preparations and no endpoint failures. Source-to-score
latencies were134.8345s and133.3738s. Allocations were137s and135s on one H200/
32CPUs/200000MiB each: **8704allocated core-s and272GPU-s total**. Scheduler GRES
accounting omits GPU entries; pinned wrappers and actual CUDA worker receipts
record the one-GPU allocation for each job.

All6source records,12MACE receipts and24GFN2 endpoint records were checked against
their hashes and component algebra. All medians, spreads and unchanged-band
calls recompute exactly. The preceding implementation's23release/ensemble tests
passed; no scorer code changed for this sensitivity experiment. Full traces are
under `workspaces/plm_fold_sensitivity_20260920/`; compact all-member values,
historical admission details, costs and comparison checks are in RESULT.json.

## Recommendation

Do not expect the three-fold median alone to resolve these two predictions. Keep
both as unresolved, out-of-reference-domain predictions while evaluating the
already-running structural-response checks. Do not call this an accuracy failure
or success without biological labels, and do not rescore the production PLM panel
or promote an ensemble default from these two cases.

[Frozen scope](PLAN.md), [commands](COMMANDS.md), [source/plan manifest](MANIFEST.json),
[actual scheduler accounting](SACCT.txt), and
[reference integration/limitations](../pqq_ensemble_20260920/REPORT.md).
