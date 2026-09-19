# Frozen response-model continuation

Parent-approved 2026-09-19, before new outputs. Preserve all existing models,
weights, thresholds, water states, coordinates and defaults. No folds/rescores.

1. GGR 2FW0/2FVY: reuse existing `stage_b_prepared_v1/*/alpha_caps/`
   preparation_manifest.json. Each 58-atom core uses the same alpha-cap policy
   as GGR_extended/alpha in the response archive, original frozen protonation,
   no waters, Ca charge -1/La0. Run four native MACE-POLAR-medium vacuum analytic
   force endpoints, unchanged model/software/kernel from the fitted PQQ panel.
   Existing stageB DFT SPs use exact same coordinates/Hamiltonian but NormalSCF;
   initial direct response centers use TightSCF. This difference is reported
   before execution; TightSCF repeats require root coordination, not silent reuse.
2. 8GY2/O05542: use already prepared canonical fixed-core v3 pair, 71 atoms,
   Ca -2/La -1, no water, oxidized PQQ3-. Run two native MACE core forces and two
   native r2SCAN-3c/CPCM/DefGrid3 SPs with its unchanged canonical input policy.
   Heme/other chains remain excluded by the established core policy. This is
   not a validated heme-free whole-protein preparation. Evidence is Ca-associated
   enzyme/structural transfer, not measured direct La/Ca affinity. Source max
   sequence identity43.7068% gives a new group under the previous50% rule.

Total initial manifest:6MACE gradients+2DFTSP. Existing A5000/16CPU/64474MiB runner
for MACE;16MPI ranks per DFT endpoint, concurrent endpoints on standard/memory.
Historical core MACE model cost~1.1s/endpoint excluding process/startup. Historical
GGR58atom SP~2–3min/endpoint on16ranks;8GY271atom baseline expected same order.
Measure actual cost, retain failures, no numerical gradients/geometry changes.

Use saved all25PQQ DFT-only and DFT+radial models without refit. Report both heads
on8GY2 at fixed logit0, alongside separate canonical baselineS/released bands.
Keep structural association target distinct from affinity. No adaptive threshold.
For all three GGR structures against both alphas report alpha-minus-GGR logit
differences only, plus original connected1GLG partition separately. The two new
GGR structures are robustness tests of the same protein, not new biological groups.
Reuse exactly compatible original1GLG/alpha rows; never pair newwater-stateDFT
energies with old force arrays. No absolute cross-target affinity classifications.

Essential checks: pinned inputs/method/software, exact pair geometry/charge,
source atom/cap mapping, nativeCUDAdevice, force and score algebra, actual DFT
normal termination and execution receipts. Use existing qualified worker and
runner, not a new potential/backend. A generic explicit-preparation scoring
operation will expose the same frozen mapping/model for opt-in scanner use.

Root approved the four TightSCF GGR SPs before new outputs. Final finite inventory:
6MACE gradient calls and6DFTSP calls (four GGR TightSCF, two canonical8GY2 inputs).
DFT executes up to four endpoints concurrently,16MPI ranks each. Original GGR
NormalSCF results remain separately archived; no broad numerical revalidation.
