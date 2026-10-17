# Mex electronic-response comparison: prepared, awaiting reference qualification

**The eight-cell comparison is prepared and its saved-data tests pass. No new response DFT calculation has run.** The next scientific decision depends on the independent DFT scouts, not further cache repair.

## What this comparison will answer

At the existing four Mex geometries, compare both metals under embedded PBE0 plus the same MM and fixed-charge solvent terms. This replaces the approximate local electronic terms while keeping the scaffold treatment fixed. The comparison can expose electronic-response disagreement; it cannot independently establish that the shared MM stiffness is accurate.

All 2,049 physical atoms, 167 waters, three Nd spectators, electronic states, boundaries and source identities remain unchanged. These are conditional EF3 response tests, not binding affinities or validated preference predictions.

## Actual verification

- Preparation **13067**: all eight tasks prepared; the original seven-geometry/fourteen-cell menu reproduces exactly. Full maps, state metadata, native decks, MM ledger and saved solvent components pass.
- Preflight **13068**: all eight preparations replay successfully. Full derivative qualification and both measured-cost decisions remain explicitly absent.
- Regression **13068**: five controls pass, covering valid saved preflight, missing or incompatible qualification, failure after provisional success, and rejection of a deliberately corrupted in-memory gradient test copy. No molecular observation was fabricated or modified.
- All preparation, diagnosis and test calls together: **zero new molecular evaluations**, 143 summed scheduler seconds, **1,716 allocated CPU-seconds**, no GPU allocation. Actual TotalCPU and batch memory records remain separate in `preflight_v2/SACCT_ALL_PREPARATION.txt`.

## Preserved preparation failures

13059 failed exact terminal replay in a different numerical environment. Restoring the original NumPy environment resolved that check; 13060 first exposed a missing Gemmi import. An isolated, byte-pinned Gemmi/native-library overlay supplies that parser without modifying either existing environment.

13063 then exposed a cache representation error. Actual diagnosis13065 checked all fourteen cells: native energies, complete physical coordinates and full gradients were exactly unchanged. Only six old origin/pilot cache wrappers had different paths. Its deliberate nonzero exit records the reproduced mismatch, not a failed molecular calculation.

The reviewed V2 adapter permits only those diagnosed representations after exact content checks and retains the archived authoritative pins. The eight intended reference cells still require their exact point pins. Numerical tolerances, native acceptance, physical state, source selection and execution gates are unchanged. Original source and every failure remain preserved.

## Restart and next decision

Read `../ACTIVE_CHECKPOINT.json` and the actual terminal receipts before submitting anything. Preparation and tests are complete; **do not rerun them**.

- **12977/12978:** existing Mex La external-charge derivative pair. Admit actual raw endpoints and collect the pair using the original qualification adapter. Actual derivative residuals and measured costs decide the remaining six checks.
- **13051:** existing whole-peptide before-La DFT force/cost scout. Its force comparison tests the large terminal response without inheriting MM stiffness. Further whole-peptide cells require that actual result and cost review.
- Only a passing full eight-point Mex derivative qualification permits the first response reference, `La_terminal__La`, after root records the measured-cost decision. The other seven cells require that actual response scout and a second recorded decision.

Exact successful commands, configuration pins, dependency hashes and watcher records are in `preparation_v4/`, `preflight_v2/`, and their pinned workspace submission records. `RESPONSE_REFERENCE_V2.py` is the tested source (`3acf0887`); `REGRESSION_V2.py` is `5178ecdc`. PQQ production is unchanged. The mechanics goal remains active and unmet.
