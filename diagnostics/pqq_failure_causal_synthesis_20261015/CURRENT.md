# Current causal investigation

The active user goal remains: "keep working until we get a conclusion and figure
out why our predictions are messed up in these cases."

The cavity exclusion branch is complete and closed; all callbacks through4587
were acknowledged. No repeat chemistry or wider cavity sweep is pending.
The nine-source causal synthesis in review/REPORT.md distinguishes established
solvent/calibration effects from source-specific unresolved mechanisms.

Current bounded tasks:

- Saved ExaA MACE force projection onto explicitly defined PQQ internal/rigid
  motions. Eight real cached arrays; no DFT gradients exist for these endpoints.
  Heisenberg prepared pqq_force_inventory/. Root submitted the saved-only
  analysis as4598 (one shared CPU), now complete8/8. Exact command
  and config870b906b are in workspaces/pqq_failure_causal_synthesis_20261015/
  pqq_forces_v1/SUBMISSION.json. Strong quinone-shortening loads are shared;
  a distinct rigid-PQQ pose response persists across sources. Read
  pqq_force_inventory/RESULT_REPORT.md. No new molecular calls; DFT forces remain
  unavailable. Actual cost2wallseconds/1CPU; algebra checks are not force-model
  qualification. Actual queued wake has now been received and ROOT_ACK written.
- Exact MMOL1770 Ca0 six-cell native vacuum reference, common3 pool under both
  metals. Read MMOL_REFERENCE_PLAN.md. Curie owns mmol_inventory/, Gauss owns
  mmol_reference/. Source1 is inventory only. Preflight4601 passed all six
  real-input checks and actual archived parser/collection regressions. Two native
  origin scouts4603(Ca)/4604(La) submitted,12MPI each, with terminal and health
  wakes. Read mmol_reference_v1/SCOUT_LAUNCH.json and each task'sSUBMISSION.json.
  Four proposal cells are prepared but await actual origin state/cost review.
- Finite intact-PQQ pose experiment: read PQQ_POSE_PLAN.md. Curie prepares24new
  MACE cells with +/-1degree cofactor rotations across BOTH sources and ALL three
  donor candidates; twelve baselines reused. Heisenberg reviews mapping/geometry.
  Geometry preparation4606 completed24/24 (33wallseconds/one shared CPU); reviewed source code
  and source-union gates are pinned in pose_preparation_v1/manifest.json SHAeb2763a1.
  Terminal wake3352635 is armed; exact command in that workspace/SUBMISSION.json.
  Both sources pass frozen internal/paired/exterior-contact gates. No pose energy
  call had run at preparation. Actual queued preparation wake is now acknowledged.
  Independently reviewed warm scorer4608 is RUNNING on node-64-768g-8,7CPUs;
  actual preflight passed and model loading began. Exact config3fadca00 and
  submission/watcher3372210 are in pose_scoring_v1/SUBMISSION.json. This executes
  the24 new vacuum cells and reuses12 exact original cells. No ideal bond target,
  state change or optimizer. The conditional solvent adapter in pqq_pose_solvent/
  is now independently reviewed and ready; execution requires the complete actual
  MACE matrix and declared response gate. Existing corrected standalone/component
  kernels are reused, not rebuilt. No solvent job has been submitted.
- Curie reviews existing same-protein source evidence and prior experiments for
  the separate XoxF1A0ACD6B9F2 Ca4/Asn290 ambiguity in a0ac_asn_inventory/.
  This is a bounded read-only investigation, with no new molecular task declared.

All scientific preparation/analysis/calculation runs on Slurm. Login-host work
is metadata, reading, editing and submission. Other users/sessions' jobs and the
production discriminator remain untouched. Current partial plans are not results.
