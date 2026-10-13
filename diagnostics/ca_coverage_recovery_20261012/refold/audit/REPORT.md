# Q9 two-Ca geometry audit: all five source sites identified

Actual inference job 2881 produced all five prescribed seed 1 samples. Actual one-thread Slurm audit 2884 finds one schema-valid PQQ and exactly one PQQ-associated Ca in every sample under the frozen strict <4.0 Å rule. Both native Ca ions and every one of the 4849 source atoms, occupancies and coordinates are preserved exactly in each full normalized mmCIF. The only source edit is the schema-validated PQQ residue name; all five sources retain zero hydrogens. No model ranking or sample selection occurred.

| Sample | Active chain | Active–PQQ nearest (Å) | Active typed CN | Spectator–PQQ nearest (Å) | Ca–Ca (Å) |
|---|---|---:|---:|---:|---:|
| 0 | B | 2.421077 | 6 | 25.039368 | 30.660849 |
| 1 | B | 2.421406 | 7 | 25.145292 | 30.777659 |
| 2 | B | 2.400389 | 7 | 25.097887 | 30.754069 |
| 3 | D | 2.403414 | 7 | 25.197612 | 31.174870 |
| 4 | B | 2.420167 | 6 | 25.052314 | 30.695966 |

Every spectator has typed CN 6 with Asp45 O/OD1, Thr48 O/OG1 and Asp51 O/OD1, matching the structural-site donor identities mapped in the prior 1FLG audit. Active sites have the expected PQQ donors and Glu213/Asn300, with Asp350 and Glu213 denticity varying under the unchanged 3.1 Å typed-contact definition. Full contact lists, exact ions, source IDs and model/normalized-source hashes remain in the per-sample audit and terminal collection.

**Sample 3 uses Ca chain D as active and chain B as spectator.** The other four reverse that assignment. Downstream preparation must use actual per-row selectors; input ligand order does not fix physical occupancy.

## Scope and costs

This is source geometry qualification only. No hydrogens are added to outputs, no context/spectator boundary or full chemical-state admission has run, and no MACE, GFN2 or scoring call occurs. PQQ helper metadata includes its validation microstate and internal proposed H coordinates; those are explicitly not source protonation and do not enter output structures. Original 125 outcomes and corrected score bands stay frozen.

Inference receipt: workspaces/ca_coverage_recovery_20261012/refold/run_v3_h200/folds/q9z4j7/EXECUTION.json. Five outputs; wall 118.328887s, childCPU 212.992071s, child peakRSS 10528928KiB. The actual log records 69.79 s model inference plus 1.19 s output extraction (70.98 s combined); the wrapper wall includes surrounding container work. Audit 2884: wall 1.226249s, CPU 1.111978s, peakRSS 72488KiB. No repeated fold or local scientific execution occurred.

Authoritative result: `workspaces/ca_coverage_recovery_20261012/refold/audit_q9_v1/execution/FINAL_AUDIT.json`, SHA256 082762f924d51b7fb263f6c28c71ffc28abc224d0fd7f396b869c1223a7a222b. Root owns the already-declared remaining protein folds and spectator-aware preparation; no downstream result is claimed here.
