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

The predeclared largest three-residue window is A77–79, only7.42% of outside-QM peptide squared motion; about49 effective residues participate. This supports distributed response and limits what one local stiffness reference could establish. [Motion report](mex_motion/result/REPORT.md); topology feasibility review is complete in mex_motion/FEASIBILITY.md. A connected A77–95 expansion is possible, but its new hydrogen-bond neighborhood is unqualified and the motion remains distributed. No expanded DFT is launched; one local segment would not validate global MM lowering.

## Active independent references

Embedded PBE0 preparation12976 passes two reused origin-gradient checks and constructs the eight declared displaced tasks with zero DFT calls. First external-charge La± pair12977/12978 is RUNNING,12MPI each on node-64-768g-28. Zero-call hostcheck12979 passed before moving only the originally pending node constraints. Same job IDs, original inputs, requested resources and watchers remain. Both actual manifest-aware SCF-health watchers are armed.

Read mex_dft/README.md for the exact saved collector command. When both actual endpoints finish, run its one-CPU scout collector. Actual derivative admission and measured cost must precede the remaining six endpoints. No DFT derivative result exists yet.

## Hybrid LBT: first force gate passed

Actual-host diagnostic12984 isolated derived cap-Jacobian roundoff only (all coordinates exact). Versioned evaluator686eba84/helper9990c5b3 permit only bounded derived-coefficient roundoff; physical, state, energy and FD gates remain unchanged. Actual zero-call preflight12985 passes four maps,12 negative controls and eight native replays. Original12975 failure remains immutable.

Nine-cell scout12986 completes and passes every component/total derivative in both directions at both step sizes. Largest total residual is1.30947e−5 kcal/mol/Å against the unchanged minimum0.005 tolerance. Actual45.04 application seconds:9MACE/18GFN/9MM/8newALPB, one oldALPB reused. This qualifies only the declared baseline7CCO/La directions, not physical accuracy or all states.

Eight-origin matrix12987 is complete:7new composites/14GFN/7MM/0newALPB,49.67applicationseconds, all8states available. Zero-call12988 passes parent replay. Remaining24derivatives12989 are COMPLETE: all three cases (baselineDy/expandedLa/expandedDy) pass both directions and both step sizes;24MACE/48GFN/24MM/16newALPB,142.09applicationseconds. Directional force consistency is qualified; electronic accuracy is separate.

Complete-pool zero-call12990 passes26geometries/104requiredcells and52exact priorALPB receipts. Actual12992 stopped on the59thnew composite: 7CCN/baseline/La at saved cartesian_03 has isolated-core GFN2 SCC failure after500iterations; its embedded component passes. The104-cell denominator is66complete/1failed/37unattempted. Both7CCO menus are complete, but the full comparison is unavailable. New optimization stays gated; no favorable subset or alternate region is substituted. Actual360.37applicationseconds,59MACE/118GFN/59MM/0newALPB; originalfailure preserved. Read lbt_hybrid_pool/result_v2/REPORT.md and ROOT_REVIEW.json.

The last isolated iterations oscillate substantially in charge and energy; this is not a rounding-threshold failure. Root is diagnosing a finite same-Hamiltonian numerical repair, with unchanged state, temperature,500iteration ceiling and convergence gates. No retry is launched yet. The baseline-first response rule was recorded before poolresult inspection, and its source-only wrapper passed independent review; it remains unbound/unsubmitted until fullparent admission. Exact pins/commands are in HYBRID_EXECUTION_CHECKPOINT.json and workspace receipts.

## Ownership, resources and restart

Root: all launches, integration, status, commits/pushes. Mechanics agent: isolatedGFN failure diagnosis. Harvey: exact partial-result/cache reuse review. Leibniz: unbound finalcomparison adapter; fullparent gate retained. Response wrapper is independently reviewed but gated. No agent submissions. Exact jobs/configs/watcher commands are in workspace SUBMISSION.json and INITIAL_EXECUTION_RECORD.json.

All numerical work uses Slurm. Standard policy currently12CPU/job120CPU/user; H200 partition separately limits one concurrent job/GPU per user. Allocation-derived native resources, no memory/time overrides or priority changes. Preserve failures, source inventories, charge/spin conventions and original pools. No source-position springs, rejected clamped-core branch, automatic optimizer extension, or production promotion. Further response budgets require a concrete diagnosis and a separately recorded finite scope.
