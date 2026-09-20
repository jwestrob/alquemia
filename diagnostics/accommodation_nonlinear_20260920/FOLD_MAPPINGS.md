# All-primary reference donor mappings: complete

**All 208 supported primary folds have compatible physical mappings.** The full
225-fold denominator retains the 17 original preparation failures unchanged.
This is geometry preparation for the frozen fold-transfer experiment, not an
energy calculation or discrimination result.

| Source conditioning | Supported maps | Original unavailable | Total |
|---|---:|---:|---:|
| Ca | 110 | 15 | 125 |
| La, noncanonical | 98 | 2 | 100 |
| All primary folds | 208 | 17 | 225 |

The adapter reused **57** exact compatible mappings and made **151** new ones.
Five of the previous 62 maps belong to canonical replays excluded from this
primary population. No distance warning, score or label chose the cases.

Every supported case exposes terminal anchor-Glu chi3. Its extra homolog position
adds Asp chi2 only when the actual residue is Asp: **93 two-mode cases and 115
Glu-only cases**. Other original donor/scaffold coordinates are inactive. Both
carboxylate oxygens move together; the existing source/cap graph and pair checks
are reused. The maximum q0 arithmetic difference is 8.88e-16 Å, below the existing
1e-12 Å mapping policy; source XYZ files were not rewritten.

## Reusable output

`workspaces/accommodation_nonlinear_20260920/fold_maps_v1/design.json`

`cases` contains 208 supported proposal-compatible records, including exact
Ca/La context XYZ/charge/multiplicity pins, shared physical-map pins, source core
and context preparation, active mode IDs, protein group, root reference ID and
source-metal conditioning. The experimental class is explicitly named
`expected_class_for_later_report_only`; it does not affect mapping.

`rows` contains all 225 original IDs in order, including each unavailable reason.
Model/software/ORCA source pins are included for the proposal adapter; no engine
ran. The independent proposal agent has the completed design.

## Actual execution

Job **1204176** completed normally on64CPUs: 3.968 seconds in the mapping process;
Slurm records5 allocated wall seconds, **320 allocated core-seconds**, and
122.780 actual CPU-seconds. No GPU, ORCA or MACE call. Batch MaxRSS reports zero,
which is unavailable measurement, not measured zero memory use.

The frozen executed source closure is
`workspaces/accommodation_nonlinear_20260920/fold_maps_implementation_v1/`.
The exact wrapper is [run_fold_mappings.sbatch](run_fold_mappings.sbatch), and
the approved population/rule remains [FOLD_TRANSFER_PLAN.md](FOLD_TRANSFER_PLAN.md).
Do not resubmit the completed geometry job. All source preparation failures,
canonical exclusion, existing map files and production behavior remain unchanged.

Tests exercise the real primary population, actual Asp/ALA/SER/THR homolog maps,
and the completed 225-row output with exact original coordinate pins. They call
no molecular model and do not assert unmeasured scientific success.
