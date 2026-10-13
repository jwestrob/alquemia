# Ca reference coverage recovery — 12 October 2026

Jacob authorized recovery of the eighteen unavailable Ca-conditioned models.
Keep the original 125-source result and frozen 25-reference calibration intact.
No production promotion or new La/Dy research is part of this recovery.

## Completed recovery

Both interrupted MxaF models are recovered and classify correctly as Ca under
all three policies. Native job 2869 and corrected-solvent job 2873 completed;
all cached-component replay and new component gates passed. The separate recovery
overlay gives **107 correct / 1 wrong / 1 inconclusive / 16 unavailable**. The
original result remains 105 / 1 / 1 / 18. See mxaf/RESULT.md and its exact pins.

MMOL_2048 XoxF4-1 Ca sample0 had an intact Cys127–128 disulfide; one existing
Cys128 HB2 was too far from CB to receive a source-graph parent. Diagnostic 2864
and repair 2870 completed. Only that H moved 0.287 A to the ff19SB 1.09 A C-H
length. Heavy atoms, original core, protons/waters and gates are unchanged.
The 199-atom envelope and frozen state pass; both real controls replay.
Native job 2875 completed all six matrix cells under the original proposal policy.
Corrected-solvent job 2880 completed and passed component gates. The recovered
MMOL source is **inconclusive under accommodation**, while its static origin is
correctly La-like under its own frozen bands. Do not tune this result away. The
combined recovery overlay is **107 correct / 1 wrong / 2 inconclusive / 15 unavailable**;
all 125 primary rows remain intact. See disulfide/solvent_v1/FINAL_COLLECTION.json
and FINAL_COLLECTION_PRIMARY_RECOVERY_OVERLAY.json in the matching workspace.

## Fifteen models need a distinct source-occupancy hypothesis

Audit 2866 found valid PQQ in every source, with Ca 25.0–25.2 A away. No metadata
repair is justified. Five models each belong to C5AXV8 ExaF, Q88JH5 PedE and
Q9Z4J7 ExaA. Raw heavy coordinates and normalization are unchanged.

The ExaA 1FLG crystal provides independent evidence: the exact Q9 Asp45/Thr48/
Asp51 positions map to crystal Asp11/Thr14/Asp17 in both chains, where they
coordinate the real second, structural Ca site. The other Ca occupies the PQQ
site. Thus the original single-Ca fold placed its ion at a real alternative site.
See refold/REPORT.md and https://www.rcsb.org/structure/1FLG.

Successful preparation job 2879 on H200 verified all three exact AF3 cached inputs
and mapped the crystal. Each new input preserves every original field and adds
only a second Ca ligand. This occupancy is a homology-based hypothesis for PedE
and ExaF, not experimental proof for them. No new MSA search. AF3 3.0.4 and actual
weights are pinned, seed1, five samples, ten recycles; original fold runtime and
recycles are unverified, so this is not a controlled one-variable fold comparison.

**Q9 scout 2881 and source geometry audit 2884 COMPLETED**: all five sources
have exactly one PQQ-associated Ca (2.400–2.421 A) and one retained structural Ca
(~25 A from PQQ). All 4849 source atoms/coordinates/occupancies are preserved.
Sample3 uses active chain D/spectator B; the other four use active B/spectator D.
Use per-row selectors. These are source-geometry results, not scores.

The same frozen AF3 continuation is complete: C5AXV8 fold2912/audit2925 and
Q88JH5 fold2913/audit2931 each retain all five models and pass the same unique
PQQ-Ca geometry rule. Thus **all fifteen newly generated sources have a PQQ-bound
Ca plus the preserved structural Ca**. This resolves source placement, not final
chemical preparation or classification. C5 sample1 selects activeD (other fourB);
Q88 selects activeD in all five. Always consume the actual selectors.

Each fold used12CPUs/oneH200 and took about118s in its wrapper; no new MSA search,
memory/time request or parameter change. GPU QOS imposed12GiB hostRAM. All outputs,
terminal receipts, completion wakes and root acknowledgements remain in the
matching refold workspace. Do not rerun any of these completed fifteen folds.

## All fifteen preparations admitted; scoring in progress

Spectator-aware preparation2969/2985 completed **15/15** under the same new source
adapter. The previously scored one-Ca P12293 sample0 exactly reproduces its old
core/context coordinates, maps and state. Both native Ca and all source heavy
coordinates are retained exactly. All paired local contexts have192atoms and
match their own frozen canonical state and fragment membership. Spectators are
23.74–24.21A from the nearest retained heavy atom, outside the unchanged4.3A
boundary. Their exclusion from the local Hamiltonian is explicit; this is not a
global electrostatic model. See spectator_preparation/FINAL_REPORT.md and actual
FINAL_PREPARATION.json in spectator_preparation_q9_v2 and _c5q88_v1 workspaces.

The existing native CPU and corrected-solvent kernels now score these real inputs:

| Batch | Native jobs (7CPUs each) | Collector | Corrected solvent (12CPUs) | Completion wake PID |
|---|---|---|---|---|
| Q9/ExaA, five sources | 2996,2997 | 2998 afterany | 2999 afterok collector | 3714564 |
| C5/ExaF + Q88/PedE, ten sources | 3018,3019,3020 | 3021 afterany | 3022 afterok collector | 3731153 |

Q9 native preflight and both real force/energy qualification fixtures pass in
both lanes. Model searches are running; no new composite classifications exist
yet. The same bounded four-mode search, common candidate pool, corrected GFN2
components and frozen bands apply. No source ranking, threshold refit or new
biological labels. Normal scheduler/user CPU limits may queue a lane; no bypass.

Exact submission/dependency/wake commands are in the matching
`two_ca_scoring_q9_v1` and `two_ca_scoring_c5q88_v1` workspaces. Both batch overlays
preserve primary125 plus the earlier three recoveries. The final cumulative
report must combine both actual5+10 results; neither batch alone is the complete
fifteen-source outcome. Gauss owns the separate finite merger; root integrates.

## Runtime failures retained

Job 2877 stopped before scientific work because the node-local container launcher
has a different hash from login. Technical v3 records each actual node launcher
while retaining exact shared container/weights/input pins. Job 2878 then exposed
an Apptainer 1.1.3 user-namespace startup failure on CPU node5; schema never ran.
H200 Apptainer 1.5.4 ran all actual schema/reference checks successfully as 2879.
No privilege flags, installations or policy changes were made. Original failed
outputs remain. Unrun disulfide native_v1 metadata had an ID suffix bug; native_v2
corrected it before any molecular call. No completed chemistry was repeated.

## Resume precisely

Root owns integration and submissions. Curie/Heisenberg preparation and refold
reports are complete; Gauss finishes the small cumulative result merger.
Read actual Slurm status and these receipts before acting:

- Completed disulfide/solvent_v1/FINAL_COLLECTION.json and its primary/recovery
  overlay retain the new inconclusive result. No further MMOL optimization.
- workspaces/ca_coverage_recovery_20261012/refold/audit_q9_v1/execution/FINAL_AUDIT.json
  after2884. Read every sample and distinguish implementation errors from
  source geometry failures; this audit is not full endpoint preparation.
- All three fold EXECUTION.json and audit FINAL_AUDIT.json are complete. Use
  their exact per-row source/active/spectator selectors; no fold/audit reruns.
- Do not repeat completed preparation2969/2985. Inspect the exact native lane
  receipts and the dependent NATIVE_COLLECTION.json in both scoring workspaces.
- Read corrected-solvent outputs under two_ca_solvent_q9_v1 and
  two_ca_solvent_c5q88_v1 when their jobs finish. Retain unavailable cells and
  actual costs. Combine the two disjoint new-source batches on the same immutable
  primary125/prior-three-recovery ledger; no automatic production promotion.

Existing LanM 2665/2666 and collector2669 remain preserved. Readiness probes stay
off; owned completion wakes are armed. Avoid stale/mismatched node-48-256g-18.
All science uses Slurm. Scoped commits are backed up to the existing main branch;
check live git history for the latest push. No push of unrelated dirty files.
