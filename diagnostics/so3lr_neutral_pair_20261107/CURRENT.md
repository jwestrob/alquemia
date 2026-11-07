# Current state — 7 November 2026, 00:37 UTC

The active goal is to continue until there is a good reason to stop. Read
[PLAN.md](PLAN.md) for the protocol declared before these evaluations.

The neutral cis/trans preparation passes the distant-water cancellation check.
Across three starting conformations, two water orientations, two separations,
and the declared ligand motions, the largest energy change is **0.00694
kcal/mol**, below the declared 0.5 threshold. The directional force check also
agrees with the finite difference. These checks establish numerical consistency
and cancellation in these preparations; they do not establish selectivity
accuracy or a binding free energy.

Three full periodic systems are prepared: cis/trans ligands, La/Dy, two K+
counterions, and 2,400 waters, totaling 7,288 atoms. Stereochemistry, charge,
atom correspondence and periodic contacts passed the preparation checks.
These have no added background salt. The proposed experimental electrolyte is
still unverified, so these are development systems.

Completed Slurm jobs, all `standard_expand`:

| Job | Purpose | Result | Allocated CPU seconds |
|---|---|---|---:|
| 52753 | Packmol build and three preparations | Passed | 81 |
| 52797 | Neutral finite-system controls, 128 evaluations | Passed | 1,320 |
| 52830 | First full periodic CPU evaluation | Out of memory; no prediction | 816 |
| 52877 | Same periodic calculation, 24 CPUs/default memory | Both endpoints completed | 2,736 |
| 52948 | Start-0 translation and remote-water controls | Partial; internal time limit reached | 7,968 |
| 52991 | Saved-data and bibliography checks | Passed | 1 |
| 53022 | Independent saved-data replay | Passed | 1 |

Total so far: 12,923 allocated CPU seconds (3.59 CPU hours), 4,806.566 measured
CPU seconds (1.34 CPU hours), zero GPU time. The CPU recovery changed resources only. Its
peak resident memory was about 52.4 GiB; the failed allocation provided
48,000 MiB (46.875 GiB). Raw inputs, outputs, forces and partial charges are preserved under
`workspaces/so3lr_neutral_pair_20261107/`.

The completed periodic controls are small: rigid translation changes the swap
energy by -0.00000160 kcal/mol; removing the first remote water changes it by
-0.001032 kcal/mol. That water is 27.04 A from the closest solute heavy atom.
The metal force-difference vectors change by 0.070/0.159 kcal/mol/A. These are
useful robustness results, with a residual force sensitivity; two water tests
and the other preparations remain to be evaluated on GPU. Job 52948 reached
its declared internal timeout. Its saved results are retained; do not retry it.
The independent replay reproduced both CPU summaries and the finite-system
results exactly, using no new model calls.

The native bibliography now links cis599 to logical identifier440 and its sole
La/Dy source, Charlier, Merciny and Fuger (1987), DOI
10.1016/S0003-2670(00)85691-3. Trans600 uses logical identifier441 and has multiple
source papers. Original methods/tables and the exact source for trans's 1.0 M
row are still needed before experimental admission. See PROVENANCE_STATUS.json.

**Job 52844 is pending for an H200** (`QOSMaxGRESPerUser`), with the full periodic energy/force
controls and measured throughput. It uses the primary L model, float64 and
the same physical charge/spin. The internal bound is 900 seconds. Its completion
watcher is recorded in the workspace's `WATCHERS.json`; do not submit a duplicate.
Other sessions' allocations remain untouched. Scheduler waiting is not evidence
against the model.

Next: assess translation/wrapping consistency, remote-water sensitivity and
actual H200 speed. Only then allocate unrestrained endpoint dynamics within the
total one-H200-day budget. Do not infer a free energy from the current endpoint
energy difference. No new DFT, training, protein simulations or production
baseline changes are part of this continuation.

MD_SCOUT.md and md_scout.yaml define the conditional first 0.25 ps trajectory,
including bounded packing relaxation, state, thermostat and output precision.
It has not been submitted. Its launcher requires a root review record of the
periodic controls and measured throughput, not another user approval.
