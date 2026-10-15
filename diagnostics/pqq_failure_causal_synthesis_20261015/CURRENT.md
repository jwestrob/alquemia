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
  qualification. ROOT_REVIEW records proactive collection; ACK awaits actual wake.
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
  This branch is not yet submitted. No ideal bond target, state change or optimizer.

All scientific preparation/analysis/calculation runs on Slurm. Login-host work
is metadata, reading, editing and submission. Other users/sessions' jobs and the
production discriminator remain untouched. Current partial plans are not results.
