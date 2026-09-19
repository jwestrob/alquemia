# Complete-context coordination proposals: no discrimination gain

All ten native MACE searches and ten original-core native DFT checks completed.
The three expected orderings remain correct, but every separation narrows and
the alpha-lactalbumin structural spread worsens. This preparation policy did not
improve discrimination or robustness. Baseline/default remain unchanged.

| Expected La-like minus Ca-like, kcal/mol | Exact starting DFT | Proposed geometry DFT |
|---|---:|---:|
| 4MAE XoxF − 1H4I MxaF | 30.071899 | 26.275723 |
| Alpha 1F6S − GGR 1GLG | 10.923654 | 3.820231 |
| Alpha 6IP9 − GGR 1GLG | 14.572699 | 13.974842 |

Alpha replica spread increases from **3.649045 to 10.154612 kcal/mol**. These are
consumed development structures: the two alpha structures are one biological
group, PQQ functional labels differ from the qualified affinity-direction evidence,
and no new biological observations or blind tests were added. Changed coordinates
do not inherit old decision bands or a universal zero.

Nine of ten DFT endpoint energies decrease; alpha 1F6S Ca increases slightly,
by 0.077473 kcal/mol. All four PQQ energies decrease. The strong optimization
success therefore does not establish improved discrimination. Native MACE-core
gaps likewise change 104.027304→103.698799 (PQQ), 8.303640→2.365819 (1F6S/GGR),
and 9.613261→11.776741 (6IP9/GGR); this is a separate vacuum Hamiltonian.

## What was actually tested

Protocol `native_OMOL_context_physical_coordination_proposal_v1` uses the exact
complete 104–202-atom contexts from the previous static second-shell pilot.
PQQ starts from the completed contextual-H preparation; alpha uses the actual
improved-water-H cores. Ten physical mode sets preserve covalent bonds and source
mapping: metal translations, original donor sidechain chi rotations and actual
peptide crankshafts. PQQ/Arg, added neighbors, scaffold and waters remain fixed;
caps move only through their actual bond anchors. Each physical heavy displacement
is limited to 0.20 Å and each angle to 0.20 rad.

Native OMOL100M unmasked vacuum energy supplies an iterative geometry proposal.
The scored energy is **only** unchanged original-core native ORCA 6.1.1
r2SCAN-3c/CPCM(Water)/DefGrid3, original charges and singlet policy. No extra MACE
scalar, masked whole-minus-core term, GB, DFT tangent, entropy or label-directed
selection appears. This is scientifically distinct from the failed 2026-09-18
hybrid fixed path, although the useful physical-coordinate machinery is reused.

All searches report SLSQP convergence in 12–19 iterations. All ten solutions
touch the prescribed displacement boundary (3–8 physical heavy atoms); they are
constrained proposals, not demonstrated unconstrained minima. The tangent-cone
projected gradient diagnostic ranges approximately 6.41e-7–3.58e-5 eV per declared
coordinate unit. No curvature, entropy, thermal population or free-energy
claim follows. All fixed-water, bond, paired composition and analytic mapping
checks pass. The source protein is not assumed mechanically flexible outside
these explicitly supported coordinates.

## Execution and costs

| Job | Actual work | Wall seconds | CPUs | Allocated CPU-seconds | GPU-seconds |
|---|---|---:|---:|---:|---:|
| 1202425 | Import failure; zero model/energy calls | 1 | 16 | 16 | 1 |
| 1202426 | Ten searches, ten final core evaluations | 53 | 16 | 848 | 53 |
| 1202428 | Ten native DFT singlepoints | 464 | 64 | 29696 | 0 |
| **Total** | | | | **30560** | **54** |

Native search worker time 45.046228 s; 250 objective invocations were recorded,
plus ten final core evaluations. ASE can reuse identical positions; these counts
are not a claim of 260 distinct forwards. Peak measured CUDA allocation was
6,107,492,352 bytes. Worker peak host RSS 2,152,332 KiB; DFT Slurm batch MaxRSS
12,458,860 KiB. Preparation/tests/reporting are additional local work, not silently
included in allocation totals. This development experiment uses paired archived
starts; it is not a controlled production speedup measurement.

The first GPU attempt failed because preparation imports required gemmi, absent
from the GPU environment. Moving those imports into CPU preparation fixed it.
Scientific settings/coordinates were unchanged; `prepared_v2` and all failed logs
remain preserved. `prepared_v3` ran. No scientific retry or DFT rescue occurred.

Eight real-fixture tests pass, zero skips (`TESTS.txt`), including actual proposal
mapping and native input/state preservation. Scientific computations above are
separate from parser/geometry checks. Exact results and receipts:

- `workspaces/coordination_preparation_20260919/prepared_v3/collection_1202426.json`
- `workspaces/coordination_preparation_20260919/dft_v1/collection_1202428.json`
- `MACE_SUBMISSION.json`, `MACE_RETRY_SUBMISSION.json`, `DFT_SUBMISSION.json`

## Next useful investigation

Do not continue rescuing this geometry policy or promote it. Preserve its code
and actual negative utility result. The separate static second-shell quantum
pilot modestly improved alpha/GGR gaps by about 3.9 kcal/mol. Testing its unchanged
preparation rule on both archived GGR replicas can determine whether that
environmental effect improves the weakest existing margins or is structure
dependent. This continuation is separately frozen and does not change this report.
The prior environment audit is in `ENVIRONMENT_AUDIT.md`.
