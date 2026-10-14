# Reciprocal placement changes the saved donor response in both directions

Independent read-only review of completed preparation and MACE job 4130. **All eight new cells are available, and both directions pass the predeclared descriptive test.** The source-dependent fixed metal position materially contributes to the different donor-motion responses within this local vacuum model.

| Host | Original Ca−La donor work | Transplanted Ca−La donor work | Change |
|---|---:|---:|---:|
| Sample 0 control, receiving sample 3 metal point | +5.261497 | +8.634833 | +3.373337 |
| Sample 3 target, receiving sample 0 metal point | +16.403481 | +7.074878 | −9.328603 |

All values are kcal/mol from the actual saved collection. Work means energy at the saved `adaptive_La` arrangement minus that metal's energy at its corresponding origin. Positive Ca−La work means La benefits more from that donor motion. The target-minus-control contrast changes from **+11.141984 to −1.559956 kcal/mol**. The response does not merely shrink uniformly: the control response grows while the target response falls.

## What changed physically

One proper alignment of all 24 named PQQ heavy atoms and its inverse moved only the metal point by **0.3643115 Å** in each host. The PQQ fit RMSD is **0.0224387 Å**; determinant is one. All nonmetal coordinate lines, source membership, proton inventory, dry oxidized-PQQ assumption, caps, charges and spin inputs remain unchanged. The same transferred point is used for both host donor arrangements and both metal identities.

Both hosts pass the existing covalent/cap, displacement and severe-overlap gates. There are no new severe overlaps. The target's saved donor arrangement still reaches its original approximately 0.8 Å maximum heavy-displacement boundary; the transferred metal itself remains inside that bound. This is a source-conditioned fixed point transfer, **not an optimized metal position or a newly relaxed protein**.

## Warranted conclusion and next decision

The intervention supports a causal role for fixed metal placement in this model's source-sensitive donor response. Together with the independent DFT reproduction of the original source contrast, it makes a targeted test of the frozen-metal constraint more compelling than indiscriminate context expansion or suppression of the MACE response.

It does not establish which placement is physically correct, prove a folding error, identify a unique metal–PQQ interaction, or demonstrate biological specificity. Changing the metal point also changes its interactions with all retained host atoms. The new transplanted surfaces have not been checked by DFT, and neither source becomes an equilibrium state by passing geometry gates. The two interventions are development evidence for one consumed protein.

The declared next step—evaluate the qualified frozen-charge solvent contribution at these exact same eight cells—is justified. It will test whether the reciprocal effect survives the current composite endpoint energy. No new optimizer, calibration band, class assignment or fitted correction follows from this vacuum result alone.

## Measured scope and records

The existing eight original cells were reused. The worker records **eight new forwards, zero failed requests**, 105.7171 s worker wall time, 597.2930 CPU-s, and 6,896,616 KiB peak RSS on seven CPU threads. The whole scoring-stage receipt reports 118.3786 s, including its surrounding work. No GPU, solvent, DFT or optimization call occurred in this stage; these measurements are not end-to-end site latency.

Authoritative result: `workspaces/exaa_metal_position_20261014/scoring_v1/execution/FINAL_COLLECTION.json`, SHA256 `39d8b38c2480aa5c5ab91606edcf7738cf6934418f12d9c8d6e12dcecd03e25e`.

It pins the scoring config (`0874c1088e9eeb6ea29fb0408916b7c1926b5ca90396772392c2657ef4ba2154`), actual preparation (`f067efbfcbb8c75cc4a4923d9369c128173c824edc4a80545471677ddec00369`), all native receipts and the original inventory. This reviewer read existing outputs only; no scientific algebra, molecular calls, refitting or submissions were performed on the login host.
