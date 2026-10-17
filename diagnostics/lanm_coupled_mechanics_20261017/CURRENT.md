# Coupled La/Dy mechanics: broad response complete, accuracy checks active

17 October 2026 UTC. Jacob activated the saved PLAN.md. Root owns submissions and integration. Goal remains active and unmet; PQQ production is unchanged.

## Main result

**All six broad-coordinate searches and both source-local paired pools are complete.** All six preserve the declared physical state, and none reaches the full-force convergence target. The LBT-source discrepancy in the legacy conditional Mex/LBT exchange decreases from **24.120991 to 6.209455 kcal/mol**. The two conditional values are now **+9.681807 and +3.472352 kcal/mol**, versus +0.982619/+25.103610 before broad response.

This isolates the effect of accessible motion under each unchanged model. It does not qualify preference: LBT still uses whole-peptide MACE, Mex uses local electronic/MM mechanics, and all searches stopped at the declared80-step budget. Read [the six-search report](broad_response/REPORT.md), exact tables/results and BROAD_RESPONSE_CHECKPOINT.json.

| System/source | Metal | Additional work kcal/mol | Max/RMS force kcal/mol/Å | Job |
| --- | --- | ---: | --- | --- |
| LBT7CCO | La | −98.572405 | 2.141/0.921 |12961|
| LBT7CCO | Dy | −104.193162 |5.176/1.143|12967|
| LBT7CCN | La | −137.647592 |9.740/1.988|12968|
| LBT7CCN | Dy | −112.937903 |10.175/1.800|12969|
| Mex8FNS | La | −2945.113073 |11.054/1.418|12970|
| Mex8FNS | Dy | −2942.034643 |21.892/1.772|12971|

All80steps include the successful12-step original LBT scout. Five initial host/compiler-path failures12962–12966 remain preserved; fresh recoveries kept identical configurations. LBT pool12974 evaluates only four missing cross cells; Mex pool12981 evaluates only two. Saved comparator12982 passes all six slots, old-state/component/force replay and missing-menu tests. Source062f's pre-execution LBT-stage schema mismatch is preserved; V2f9cdb6a4 changes only LBT `run` versus Mex `pool` stage dispatch. No numerical failure or repeated chemistry is invented for that static finding.

## What moved

Saved motion12983 uses both metals and no new molecular calls. Mex MM lowering is mostly common to both metals: bonded−724 to−731, LJ−387 to−389, Coulomb−2014 to−2017 kcal/mol, partly offset by solvent. Water-oxygen aligned RMS movement is0.65–0.69Å; outside-QM peptide-heavy RMS is0.18–0.19Å. These global energy changes are not residue-energy decompositions.

The predeclared largest three-residue window is A77–79, only7.42% of outside-QM peptide squared motion; about49 effective residues participate. This supports distributed response and limits what one local stiffness reference could establish. [Motion report](mex_motion/result/REPORT.md); complete topology/HB feasibility review is pending, with no new expanded-reference submission.

## Active independent references

Embedded PBE0 preparation12976 passes two reused origin-gradient checks and constructs the eight declared displaced tasks with zero DFT calls. First external-charge La± pair12977/12978 is RUNNING,12MPI each on node-64-768g-28. Zero-call hostcheck12979 passed before moving only the originally pending node constraints. Same job IDs, original inputs, requested resources and watchers remain. Both actual manifest-aware SCF-health watchers are armed.

Read mex_dft/README.md for the exact saved collector command. When both actual endpoints finish, run its one-CPU scout collector. Actual derivative admission and measured cost must precede the remaining six endpoints. No DFT derivative result exists yet.

## Hybrid LBT failure and next action

Local-electronic/MM CPU preflight12972 passed all eight saved component replays. Its H200 nine-cell scout12975 failed in constructor origin-map replay before any molecular call: `Exact origin mapping replay differs`. Preserve its immutable execution/RESULT.json and passing CPUpreflight. Mechanics agent is preparing a zero-call exact-leaf/ULP diagnostic on the actual failed host; root submits it. No geometry or numerical acceptance tolerance has been changed and no chemistry retry is launched.

The eight-origin matrix, remaining24 derivative checks and complete104-cell saved-pool comparison are staged source only and gated on actual repaired scout admission. The fixed pool includes all26 actual source-local geometries, both metals and both regions; eight origins reused, at most96 new composites, exact old whole-solvent reuse. Never select a favorable subset or silently substitute an available baseline for a failed region.

## Ownership, resources and restart

Root: all launches, integration, status, commits/pushes. Mechanics agent: origin-map diagnostic/recovery. Harvey: independent review and bounded stiffness-region feasibility. Leibniz: paired-pool/remaining derivative adapters prepared; no submissions. Exact jobs/configs/watcher commands are in workspace SUBMISSION.json and INITIAL_EXECUTION_RECORD.json.

All numerical work uses Slurm. Standard policy currently12CPU/job120CPU/user; H200 partition separately limits one concurrent job/GPU per user. Allocation-derived native resources, no memory/time overrides or priority changes. Preserve failures, source inventories, charge/spin conventions and original pools. No source-position springs, rejected clamped-core branch, automatic optimizer extension, or production promotion. Further response budgets require a concrete diagnosis and a separately recorded finite scope.
