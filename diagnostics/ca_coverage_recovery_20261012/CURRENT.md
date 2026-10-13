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

**Q9 five-model scout 2881 COMPLETED on H200, 12 CPUs/one GPU**: all five
raw outputs exist, 125 scheduler seconds. Inference success is not admission.
**Five-source geometry audit 2884 runs on one CPU**, watcher3556764, using
refold/audit_q9_v1/manifest.json; no new energies or protonation.
Actual scheduler allocation is 12 GiB host RAM; no memory or time request was
added. Manifest and submission/watch receipts are under refold/run_v3_h200/.
Both full-source ions will remain explicit. Curie's spectator adapter is PLAN
ONLY, under spectator_preparation/; do not score until exactly one PQQ-associated
Ca is found and the spectator is demonstrably outside the retained local context.
All generated samples and failures remain visible; no score-based selection.
Other two five-model batches require an executable scout and a root submission.

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

Root owns integration and current jobs; agents finish scoped result reports.
Read actual Slurm status and these receipts before acting:

- Completed disulfide/solvent_v1/FINAL_COLLECTION.json and its primary/recovery
  overlay retain the new inconclusive result. No further MMOL optimization.
- workspaces/ca_coverage_recovery_20261012/refold/audit_q9_v1/execution/FINAL_AUDIT.json
  after2884. Read every sample and distinguish implementation errors from
  source geometry failures; this audit is not full endpoint preparation.
- For a successful Q9 scout, run source/occupancy admission under Slurm before
  implementing the pinned spectator-aware protonation path and scoring it.
  Preserve all samples. Use the same frozen corrected bands, without refitting.

Existing LanM 2665/2666 and collector2669 remain preserved. Readiness probes stay
off; owned completion wakes are armed. Avoid stale/mismatched node-48-256g-18.
All science uses Slurm. Scoped commits are backed up to the existing main branch;
check live git history for the latest push. No push of unrelated dirty files.
