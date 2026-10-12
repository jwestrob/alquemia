# Exact strict208 alternate-fold recomputation

The original strict-panel canonical 25 separation fails after fresh standalone GFN2 reevaluation: the operational gap changes from **+7.3333422580 to -4.5566019908 kcal/mol**. All 25 canonical cases have their complete 12-point candidate pools. The released extrema rule returns **unsupported_class_separation**, with bands unavailable. Fixed previous selections, operational reselection and mathematical reselection give the same overlapping extrema. The original 207-correct/0-wrong/1-inconclusive alternate-fold result therefore has no supported corrected classification under the original calibration rule.

The branch evaluated **2796/2796 fresh primary single points**, with **2796 valid and 0 failed** attempts. The scope is the exact 25 canonical references plus the 208 historically available alternate folds, each at origin/adaptive_Ca/adaptive_La × Ca/La × vacuum/ALPB water. Every Ca and La GFN2 value is fresh. Native MACE values and the original candidate geometry pool are preserved. All 225 originally declared alternate folds are retained, including their 17 original unsupported rows and reasons. All original production labels remain preserved.

## Native baseline and exact scope

Before numerical execution, all 233 original pools were reconstructed exactly from their saved native components; all 208 operational selections and the original canonical 25 band reproduced exactly. The independent audit replayed 932 mathematical/operational metal choices and every numeric score component with zero difference. The native alternate 225 outcome denominator remains 207 correct, 0 wrong, 1 inconclusive and 17 unavailable. The original inconclusive case is c5axv8 conditioned-La seed-1 sample-3.

The strict-panel canonical reference uses 14 Ca and 11 La labels. Its scored coordinate hashes are disjoint from the current 34-reference/39-CN9 matrix. This recomputation calibrates only its own exact 25 cases. The 208 supported alternate folds represent the same 25 parent proteins and inherited labels; they measure computational transfer/robustness rather than independent biological validation.

Nine recovered auxiliary files preserve six original native task manifests, the exact selector implementations and the original REFERENCES.json. All 2796 charge, multiplicity and electron observations are explicit in those manifest records and independently agree with exact-XYZ valence checksums. All states are singlets, and each Ca/La pair retains equal electron and atom counts with La charge one unit higher. Native SCF detail records are available for 300 canonical points; the 2496 alternate native SCF/population audits remain unavailable locally. The independent state audit retains 48 tasks with no explicit manifest temperature field; every fresh standalone point explicitly uses 300 K.

A single La/origin cell (a8r3s4 conditioned-Ca sample-3, both media) has different pool-source and actual native-scored XYZ byte hashes. Both actual task bodies were recovered and verified. Their atom order matches the pool source, and the maximum coordinate-component difference is exactly 1.7763568394002505×10⁻¹⁵ Å. The calculation uses the exact actual native-scored bytes for this one cell, with both source pins preserved. All other input bindings match directly.

## Calibration result and availability

| Panel | Ca maximum case | Ca maximum, kcal/mol | La minimum case | La minimum, kcal/mol | Gap, kcal/mol |
|---|---|---:|---|---:|---:|
| native | q88jh5-pqq-la_model | -405463.7132819406 | c5axv8-pqq-la_model | -405456.3799396825 | +7.3333422580 |
| standalone | q9z4j7-pqq-la_model | -405450.6906886346 | c5axv8-pqq-la_model | -405455.2472906254 | -4.5566019908 |

The corrected Ca maximum moves from q88jh5 to q9z4j7, while c5axv8 remains the La minimum. The q9z4j7 case shifts +13.5465794909 kcal/mol and c5axv8 shifts +1.1326490572 kcal/mol. Both use adaptive_La for the Ca and La evaluations. The canonical extrema now overlap by 4.5566019908 kcal/mol, exceeding the opposite-direction 0.02-kcal/mol separation requirement.

Numerical score availability and classification availability are recorded separately. A complete fresh candidate pool supplies a numerical composite score. A supported calibration band is additionally required to classify that score. The completed strict208 scores remain available, while every corrected classification is unavailable because the canonical band is unsupported. The 17 original unsupported rows retain their separate original-pool status. Primary 500-iteration failures, if present, additionally make the corresponding fresh pool unavailable and remain preserved as raw primary outcomes.

| Profile | Complete fresh alternate pools /208 | Fresh failed pools | Changed alternate selections | Changed canonical selections | Supported corrected classifications |
|---|---:|---:|---:|---:|---:|
| fixed_previous_operational | 208 | 0 | 0 | 0 | 0 |
| reselected_operational | 208 | 0 | 14 | 2 | 0 |
| reselected_mathematical | 208 | 0 | 14 | 2 | 0 |

The report supplies numerical scores, source labels, selection identities and score shifts for every original available fold. The released origin policy is retained exactly: choose the row minimum only if its improvement over origin is strictly greater than 0.1 kcal/mol, otherwise retain origin; exact ties follow origin/adaptive_Ca/adaptive_La list order. The calibration requires every designated canonical case and a La-minimum-minus-Ca-maximum gap greater than 0.02 kcal/mol. The unsupported band is preserved, with no fallback thresholds or alternative classifier fitted.

## Boundary changes from origin to accommodation

The following contrasts use the exact origin pair and the operationally selected pair at each boundary protein. MACE contrast plus solvent-transfer contrast equals total R. The selected-minus-origin component differences close each total score change; their signs are conditional on these contexts and reference cases.

| Case | Backend | Origin MACE | Origin solvent | Origin R | Selected MACE | Selected solvent | Selected R | Selected minus origin R |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| q9z4j7-pqq-la_model | native | -405347.71744580 | -120.03227726 | -405467.74972305 | -405339.19976145 | -125.03750668 | -405464.23726813 | +3.51245493 |
| q9z4j7-pqq-la_model | standalone | -405347.71744580 | -115.05938983 | -405462.77683563 | -405339.19976145 | -111.49092719 | -405450.69068863 | +12.08614699 |
| c5axv8-pqq-la_model | native | -405341.00031262 | -114.77805342 | -405455.77836604 | -405342.14557710 | -114.23436259 | -405456.37993968 | -0.60157364 |
| c5axv8-pqq-la_model | standalone | -405341.00031262 | -112.90874129 | -405453.90905391 | -405342.14557710 | -113.10171353 | -405455.24729063 | -1.33823671 |

All values are in kcal/mol. `canonical_origin_accommodation.tsv` also records the MACE and solvent components of each selected-minus-origin change and their numerical closure. A separately requested cached origin-only ablation is reported in `origin_only_ablation/REPORT.md`; it preserves this primary accommodated calibration failure.

For the corrected q9z4j7/c5axv8 pair, accommodation lowers the La-minus-Ca gap by 13.4243837078 kcal/mol: the MACE contrast changes contribute −9.6629488272 and solvent-transfer contrast changes contribute −3.7614348806. These contributions close the decrease from the corrected origin-only gap +8.8677817170 to the accommodated gap −4.5566019908. At q9z4j7, MACE contributes +8.5176843486 to the selected-minus-origin score and corrected solvent contributes +3.5684626445; the old native solvent contribution was −5.0052294192. At c5axv8, the corrected changes are −1.1452644786 from MACE and −0.1929722362 from solvent. These are conditional score differences at the frozen candidate geometries and do not identify independent chemical truth.

## Comparison with the current positive-gap reference context

The current 34/39 panel retains a positive selected gap +3.4721710483 kcal/mol using its own exact references. Its q9z4j7 and c5axv8 contexts differ substantially from the strict panel. Formal charges remain Ca−3/La−2 and selected candidates remain adaptive_La/adaptive_La, while atom and electron content change:

| Reference | Strict atom count | Current atom count | Strict electrons | Current electrons | Strict S atoms | Current S atoms |
|---|---:|---:|---:|---:|---:|---:|
| q9z4j7 | 145 | 192 | 452 | 588 | 0 | 2 |
| c5axv8 | 172 | 192 | 532 | 592 | 0 | 2 |

The elemental inventories and exact scored hashes are direct source observations. `reference_edge_context_comparison.tsv` lists both panels’ native MACE, native/corrected solvent contrasts, full-precision endpoint energies and metal populations. The current-minus-strict component differences close the composite edge-score difference:

| Edge protein | MACE contrast difference | Corrected solvent contrast difference | Corrected composite score difference |
|---|---:|---:|---:|
| q9z4j7 | -14.49692966 | +7.62316295 | -6.87376671 |
| c5axv8 | -0.73060665 | +1.88561298 | +1.15500633 |

The two panels therefore support different calibrated outcomes on different contexts. This comparison identifies their input/component differences; it assigns a causal contribution to neither an individual atom group nor an independent biological preference.

## Frozen geometric comparator

The predeclared strongest comparator is original_24start evaluated on raw original geometries: 208/208 correct on the exact supported alternate-fold membership. It remains frozen. The accommodated energy protocol chooses among its original three optimized candidates, while this geometric protocol uses the original raw pose; these are distinct end-to-end geometry-selection protocols. The comparison remains the requested end-to-end benchmark. Historical geometric checks at energy-selected poses (205/208) and geometry-selected candidate pools (207/208) are secondary records and supply no replacement primary calls.

The corrected energetic classifier lacks its required canonical separation on this original strict panel. The geometry 208/208 comparator remains available. Accuracy for the corrected energetic method is consequently unavailable under its original rule; neither the complete numeric scores nor the unchanged current 34/39 positive gap substitute for this required calibration.

## Numerical audit and resources

The installed xTB 6.7.1 CLI used GFN2, 300 K, accuracy 0.0001, NoRestart, 500 maximum iterations, singlet UHF 0, original charge, and vacuum/ALPB water. Full stdout TOTAL ENERGY values have 12 decimal places; eight-decimal JSON energy values are retained as supporting output. Every input/result/output hash, explicit state and completed-point SCF residual validates. The resource audit observes at most 4 concurrent jobs × 2 OMP threads (8 requested threads), with BLAS/VECLIB/NUMEXPR capped at 1. The maximum elapsed point time is 4.3113 s and maximum RSS 141.5938 MiB.

A harmless actual process-monitor read and one bounded xTB point passed with approved monitor access before the batch. The frozen runner holds each lifetime semaphore until its owned process group is verified empty. Any process-monitor/launcher fault sets a global STOP, cancels queued work and stops active siblings before releasing slots. Per-point bounds are 180 s wall time, 360 s CPU time and 6 GiB monitored RSS, with core dumps disabled. All primary failures and the 500-iteration policy are preserved; native Ca values and higher-iteration rescues are absent from the fresh primary matrix.

Completed standalone state/SCF checks cover 2796 valid points. Original native canonical SCF criteria are available for 300 points; the retained criterion counts are {"energy": {"pass": 300}, "max_density": {"overrun": 56, "pass": 244}, "rms_density": {"overrun": 1, "pass": 299}}. The other 2496 original native records are explicitly unavailable. Original native convergence declarations and saved scores remain unchanged.

All three printed native criteria pass at 244/300 canonical points. Complete native selected Ca/La vacuum/water criteria pass at 10/25 reference cases; all twelve candidate-pool points pass at 4/25. The selected q9z4j7 La vacuum max-density residual is 1.8264e−5 against 1e−5; selected c5axv8 Ca vacuum is 1.1225e−5 against 1e−5. Their other selected printed criteria pass. The retained strict subset lacks complete canonical coverage and supplies no filtered recalibration; `native_scf_filter_audit.json` preserves the exact case exclusions. These marginal original flags remain explicit despite the fresh standalone checks passing.

The separate [final independent audit](../independent_review/strict208_final_audit.json) verifies all 2796 raw point/output hashes, stdout energies, states, SCF checks and execution receipts; all 233 pool replays and 699 profile records agree exactly. It also reproduces the origin-only ablation, its eight decision changes and every boundary-table closure. No arithmetic, state, scope or resource discrepancies were found.

## Deliverables and reproduction

- `summary.json`, `calibration.json`, `scores.json`: full precision scores, all three selection profiles and unsupported-calibration status.
- `alternate225_diagnostics.tsv`, `per_protein_fold_outcomes.tsv`: all 225 original statuses and 25-parent label lineage, with score/classification availability separated.
- `canonical25_scores.tsv`, `calibration_edges.tsv`: exact reference scores, selected candidates and old/new extrema.
- `singlepoints.tsv`, `endpoint_transfers.tsv`, `validation.json`: all 2796 raw-point receipts, complete transfer comparisons and resource/state checks.
- `reference_edge_context_comparison.tsv`, `reference_edge_context_differences.tsv`: direct comparison against the current positive-gap reference inputs.
- `native_baseline_reconciliation.json`, `coordinate_equivalence_audit.json`, `analysis_source_pins.json`, `command_ledger.json`, `deliverable_manifest.json`: original replay, exact-input exception audit, sources and reproducible commands.
- `runs/*`: preserved full raw stdout/stderr, JSON, coordinates, population output, lifecycle and command/resource receipts.

From this directory, reproduce the lightweight postprocessing with `python3 analyze.py`, `python3 validate.py`, `python3 compare_reference_edges.py`, `python3 compare_origin_accommodation.py`, then `python3 write_report.py`, using one thread in all numerical libraries. Afterwards run `python3 origin_only_ablation/analyze_origin.py` for the separate cached arithmetic ablation. The numerical runner is `python3 run_xtb.py` under the approved eight-thread allocation and process-monitor permission. Existing incomplete primary results require explicit inspection and remain preserved.
