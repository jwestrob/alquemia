# Standalone xTB gives stable forces and useful PQQ ordering

**The installed standalone xTB backend passes the declared numerical and force
checks, while separating all four consumed examples by relative class ordering.**
This is a useful reason to test its own reference panel. It is not a calibrated
four-case classifier, an affinity result or an automatic production promotion.

All224 calls completed successfully: two fixed accuracies on80 pool cells and32
physical displaced cells, with16 analytic-gradient requests included in that
count. No new MACE, DFT, geometry optimization or molecule preparation occurred.
The reported accuracy is always0.02;0.2 is sensitivity, not a competing selectable
result. Every fixed cell is retained.

## Useful signal and its calibration limit

Composite E = archived native OMOL(vac) + standalone GFN2(ALPB−vacuum).
R = E_Ca−E_La, in model kcal/mol. The same five candidate geometries are offered
to both metal rows, with mathematical minima and the separate unchanged0.1kcal/mol
origin-retention rule. The table uses the operational result.

| Source | Known class | Native primary R | Standalone0.02 R | Native adaptive-band transfer |
|---|---|---:|---:|---|
|1H4I|Ca|-405505.420051|-405500.305981|Ca-supported|
|4MAE|La|-405444.304425|-405444.376695|La-supported|
|q88jh5-pqq-la_model|Ca|-405464.122696|-405458.330543|inconclusive|
|a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1|La|-405463.305519|-405455.984060|La-supported|

The largest Ca-class R is Q88JH5; the smallest La-class R is the known problematic
A0A3F2YLY8 Ca-conditioned sample1. Their gap is **2.346483model kcal/mol**.
That is a report of observed ordering, not a fitted decision threshold. The
problematic A0 repeat changes from a native Ca call to the correct relative La
side, while Q88 shifts toward the native decision gap. Old adaptive bands give
three correct and one inconclusive; released static bands give three correct
and the Q88 wrong call. These are explicitly incompatible-reference transfer
checks, not the new backend's final accuracy.

Both metals choose adaptive_La for1H4I andQ88;4MAE andA0 select their own metal's
adaptive proposal. No label entered selection. All four sources and their
repeated structures were already consumed development examples. No calibration
was derived from them, and no broader predictive guarantee follows.

## Numerical and derivative consistency

All4 pool cases pass the frozen0.1kcal/mol cell and0.2kcal/mol contrast gates.
Across all112 accuracy pairs, the largest absolute cell change is
0.0000199824kcal/mol. Largest pooled contrast change is0.0000047240kcal/mol.
At accuracy0.02, **all14 declared derivative quantities pass** (14/14 also pass
at0.2): four individual components, two metal solvent transfers and their
Ca−La difference at each of4MAE/Q88. Largest tight analytic-vs-fine discrepancy
is0.0000760622kcal/mol/radian. No third accuracy, step change, smoothing or retry.

| Source | Quantity | Finest centered derivative | Analytic derivative |
|---|---|---:|---:|
|4MAE|Ca_solvent|-0.20967537|-0.20971555|
|4MAE|La_solvent|-0.42183065|-0.42179476|
|4MAE|Ca_minus_La_solvent|0.21215528|0.21207922|
|q88jh5-pqq-la_model|Ca_solvent|1.22645847|1.22649454|
|q88jh5-pqq-la_model|La_solvent|0.36052116|0.36052572|
|q88jh5-pqq-la_model|Ca_minus_La_solvent|0.86593731|0.86596881|

Units in this table arekcal/mol/radian. Exact source/cap Jacobians map actual
Cartesian analytic gradients; gradient energy, atom order and coordinates pass.
Physical RMS normalization and all component derivatives are retained in the
full comparison. Numerical agreement on two modes does not establish a global
minimum or accuracy of every future gradient.

## Actual backend/state accounting

Unchanged installed xTB6.7.1(edcfbbe), pinned parameter file, fresh/no-restart,
GFN2, exact charges, zero unpaired electrons, electronic300K,500 maximum SCC
iterations. Actual outputs show20–91 iterations, correct even effective-electron
counts, and maximum partial-charge closure error1.2e−7e. No tblite substitution.

Actual ALPB uses water, gsolv[1M gas/solution], P16 kernel, GBOBC Born integration,
230-point surface, H-bond correction and no ionic screening. Its built-in
**solvent parameter temperature is298.15K**, distinct from electronic300K.
Dielectric80.2 and gsolv shift0.001080759698Eh are printed; the state convention
has no added concentration conversion but does retain this solvent parameter
shift. Both accuracy levels share all these settings. Radii/surface and state
are pinned by the binary/parameter version and actual outputs, not equated to
native ORCA by name.

Raw printed totals and Gsolv/Gelec/Gsasa/Ghb/Gshift components remain in every
output; JSON occupations/charges and gradient files are retained. The correction
uses both standalone media consistently. Backend score differences cannot be
attributed uniquely to tighter SCC convergence: different implementations and
solvent terms remain involved. No native calibration is silently inherited.

## Cost, actual tests and technical correction

Job1210445:62s on64CPUs,128GiB requested, **3968 allocated core-seconds, zeroGPU**.
The molecular executor took57.2376s (3663.21core-seconds); collection/startup
occupy the remainder. Individual calls took0.852–7.020s, median1.899s under eight
concurrent eight-thread workers. Batch peakRSS595688KiB is a scheduler reading.
Archived MACE/source computation is reused, not erased from historical cost.
This is measured development throughput, not a matched end-to-end production
speedup claim. Preparation and local reporting cost are additional.

Six real-artifact tests pass in5.548s, zero skips: exact finite input/state pairs,
accuracy-only command differences, source/cap Jacobian check, corrupted actual
state rejection, real completed output parsing and actual selection/derivative
algebra. The molecular-result test was explicitly unrun before execution.

The first collector expected ORCA-style printed charge/spin labels and marked
all224 cells unavailable despite successful executable returns. A parser-only
correction reads actual xTB net-charge/unpaired-electron fields and validates
JSON/state/SCC/solvent output. `collection_v1.json` and the executed snapshot
remain untouched; corrected`collection_v2.json` contains224complete rows. No
molecular call was repeated. This is not224 scientific failures or a hidden
successful-output substitute.

**Recommendation:** pursue standalone0.02 as a separately calibrated backend.
The next authorized test is the unchanged28 known reference/crystal sources,
comparing static origins and the minimal three-candidate pool under their own
canonical25 calibrations. No further optimization is bundled into that test;
225-fold transfer stays unlaunched until calibration is frozen. Production and
native baseline remain unchanged.

Full result: `workspaces/standalone_xtb_20260923/run_v1/COMPARISON_v1.json`.
Final collection: same directory,`collection_v2.json`; costs:`COSTS.json`.
[Commands](COMMANDS.md), [frozen plan](PLAN.md), and [software evidence](SOFTWARE_EVIDENCE.json).
