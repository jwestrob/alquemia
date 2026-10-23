# Independent review supplement

The held source/data handoff `PILOT_HELD_HANDOFF.json` (SHA256 `bd14c763d35d28250acd3d5b7b420752a404bf04a2ca0f49aaa2be80cbf3140f`) remains unchanged. No fit, molecular call or reserved-validation result was executed/opened in these reviews.

## Saved training and physical interpretation

The independent reference-chemistry agent reviewed `TRAINING.py`, `OBJECTIVE.py` and the actual saved training result `99ea27f8…`. It found the force sign/units, connected-body sums, geometric-centroid torque, category/block weights, paired cancellation,28-cell joins and exact ligand-translation slopes consistent. It found no implementation defect explaining the fixed-La control failure. The0.15Å acetamide-La discrepancy remains a coarse consistency flag, not an analytic-derivative failure claim.

The same agent independently checked all28 parameter files against the exact SI block and retained AMOEBA22 La keys. Water36/37, amide O105/class55 and C104/class53, and acetate403 were correct. Exact audit evidence is in `diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/la_baseline_audit_v1/{REPORT.md,EVIDENCE.json}`; report SHA256 begins `7f6fde7a`. The primary article `reference_chemistry_v1/assets/la_article.txt` lines181–210 distinguishes MP2-optimized structures from the CP-DLPNO-CCSD(T) distance-curve energy fitting target; SI2 lines17–32 explain BSSE. Thus MP2 geometry provenance and the coupled-cluster calibration target are distinct facts, not interchangeable descriptions.

## Held pilot source

The independent experimental-target agent reviewed `PILOT.py` SHA256 `5dc01f56a4ab8040890d055f34dcfc448cb49a3f5f8e6bd680316c68c4bd6c93` against `FIT_SPEC.md`, `FAMILY.json` and the base key. No blocking source issue was found. It confirmed all-seven-parameter plus cell cache identity, reuse of original Dy baselines, no La native evaluation, preservation of the water-stage choice in donor scoring,19/41 search-vector caps with6/8 sensitivity reserves inside25/49, and the conservative1276-program ceiling. Partial polls retain budget-limited status; sensitivity probes cannot replace the selected search candidate. Bound/identifiability flags remain separate from physical controls, and automatic validation/physical release remain false.

The template's `root_release:null` prevents native execution. The pilot has only AST/source review, not runtime qualification. Root has explicitly held it while the La/reference water discrepancy is investigated. A budget-limited best point would not establish failure of the seven-parameter family; no such optimization has started.
