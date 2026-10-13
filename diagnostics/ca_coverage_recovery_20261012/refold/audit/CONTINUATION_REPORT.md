# Remaining two-Ca source audits: ten of ten unique PQQ sites

Actual C5AXV8 folding 2912 / audit 2925 and Q88JH5 folding 2913 / audit 2931 are complete. All ten prescribed seed 1 models have one schema-valid PQQ, two native Ca ions and exactly one Ca strictly within 4.0 Å of a PQQ heavy atom. Together with Q9Z4J7 audit 2884, all 15 declared models pass this source-geometry rule. This is not full protonation, context, spectator or scoring admission. Every sample remains in the ledger; no ranking selection or repeat fold occurred.

| Protein | Sample | Active Ca chain | Active nearest PQQ (Å) | Typed CN | Spectator nearest PQQ (Å) | Ca–Ca (Å) |
|---|---:|---|---:|---:|---:|---:|
| C5AXV8 | 0 | B | 2.463773 | 7 | 25.074189 | 31.125944 |
| C5AXV8 | 1 | D | 2.474940 | 8 | 24.935333 | 31.145722 |
| C5AXV8 | 2 | B | 2.478665 | 7 | 25.027490 | 31.055170 |
| C5AXV8 | 3 | B | 2.495634 | 7 | 25.079433 | 31.049918 |
| C5AXV8 | 4 | B | 2.509520 | 7 | 25.122580 | 31.134302 |
| Q88JH5 | 0 | D | 2.425659 | 7 | 25.072189 | 30.628448 |
| Q88JH5 | 1 | D | 2.405293 | 7 | 25.060325 | 30.703631 |
| Q88JH5 | 2 | D | 2.419649 | 6 | 25.022341 | 30.631237 |
| Q88JH5 | 3 | D | 2.408854 | 6 | 25.091140 | 30.640888 |
| Q88JH5 | 4 | D | 2.445305 | 6 | 25.064338 | 30.691867 |

All 4543 C5 atoms and all 4922 Q88 atoms, occupancies and coordinates are preserved exactly in every full normalized mmCIF, including both native Ca ions. Each PQQ contains 24 source heavy atoms and passes the existing complete schema/connectivity check. Each source still has zero hydrogens: the helper validation microstate is not output protonation. Only the validated PQQ residue name changes.

C5 spectators have the six typed donors Asp36 O/OD1, Thr39 O/OG1, Asp42 O/OD1. Q88 spectators have Asp53 O/OD1, Thr56 O/OG1, Asp59 O/OD1. This prediction pattern supports the declared occupancy hypothesis but does not experimentally establish a second Ca site in ExaF/PedE; the independent crystal evidence applies to ExaA.

**Per-model selectors are mandatory.** C5 sample 1 has active D / spectator B; the other C5 samples have active B / spectator D. Every Q88 sample has active D / spectator B. Q9 sample 3 separately has active D / spectator B. Root and Curie received the actual raw/full-source pins and selectors. Neither ion identity nor PQQ carboxyl oxygen identity may be inferred from input order.

## Actual costs and reproducibility

| Target | Fold wrapper wall (s) | Fold child CPU (s) | Child peak RSS (KiB) | Audit wall (s) | Audit CPU (s) |
|---|---:|---:|---:|---:|---:|
| C5AXV8 | 117.305290 | 212.315653 | 10495036 | 1.177023 | 1.061214 |
| Q88JH5 | 117.879817 | 211.272870 | 10520544 | 1.253245 | 1.151876 |

The raw inference logs report 69.72 s model / 1.12 s extraction for C5 (70.83 s combined, as logged) and 69.47 s / 1.23 s for Q88 (70.70 s combined). Fold jobs each used twelve CPUs/one H200; audits used one actual Slurm CPU. Same frozen AUDIT.py, RUN.sbatch, scientific helper pins and distance rules as Q9 were verified from both manifests. Only target/model/output metadata and the explicit continuation-scope plan changed. No molecular code ran on login.

Artifacts: `workspaces/ca_coverage_recovery_20261012/refold/audit_c5axv8_v1/execution/FINAL_AUDIT.json` and `audit_q88jh5_v1/execution/FINAL_AUDIT.json`. Exact pins, source IDs, model maps and copied terminal evidence are in CONTINUATION_LEDGER.json and evidence/{c5axv8,q88jh5}/.

Root owns downstream integration and completion acknowledgements; Curie owns the separately declared spectator-aware preparation. These checks do not alter the original 125 outcome rows, infer new biological class calls, fit new bands or promote a production default.
