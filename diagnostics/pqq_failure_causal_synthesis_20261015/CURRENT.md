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
  origin scouts4603(Ca)/4604(La) complete as energy diagnostics (3978/3789s,
  12MPI each). Actual state/basis/coordinate checks pass; density/state stability
  remains unqualified with energy-only endings/warnings retained. Root reviewed
  both in mmol_reference_v1/SCOUT_REVIEW.json, then submitted four unchanged saved
  proposals4637–4640 (12MPI each), comparator4641 afterany. All terminal/health
  wakes recorded in taskSUBMISSION.json andPROPOSAL_LAUNCH.json; no repeat origins.
- Finite intact-PQQ pose experiment: read PQQ_POSE_PLAN.md. Curie prepares24new
  MACE cells with +/-1degree cofactor rotations across BOTH sources and ALL three
  donor candidates; twelve baselines reused. Heisenberg reviews mapping/geometry.
  Geometry preparation4606 completed24/24 (33wallseconds/one shared CPU); reviewed source code
  and source-union gates are pinned in pose_preparation_v1/manifest.json SHAeb2763a1.
  Terminal wake3352635 is armed; exact command in that workspace/SUBMISSION.json.
  Both sources pass frozen internal/paired/exterior-contact gates. No pose energy
  call had run at preparation. Actual queued preparation wake is now acknowledged.
  Independently reviewed warm scorer4608 completed24/24 on node-64-768g-8,7CPUs.
  Control paired change−0.7132293300,target+0.9782094679kcal/mol: neither reaches
  the frozen1kcal solvent gate. Close this finite angle test; prepared solvent
  adapter stays unrun. Read pqq_pose/REPORT.md and RESULT.json. Config3fadca00 and
  submission/watcher3372210 are in pose_scoring_v1/SUBMISSION.json. Twelve exact
  original cells reused; actual426wallseconds/2982allocatedCPU-seconds,zeroGPU.
  No ideal bond target, state change, optimizer or calibrated classification.
- XoxF1A0ACD6B9F2 inventory complete: Ca4 and La4 have a distinct terminal-amide
  orientation, nearby ND2/distant OD1; all other eight origins show the opposite.
  Both Asn modes already converge inside their local domain. Correctly calledLa4
  remains a counterexample; this is not evidence of a missing donor or a bad fold.
  ASN_ROTAMER_PLAN.md declares a separate whole-amide180degree origin flip across
  allten sources, source-graph/covalent/fullsourceH+heavy checks first. Curie owns
  a0ac_asn_rotamer preparation, Heisenberg mapping review, Gauss count-aware204atom
  MACE scoring. Preparation4614 completed20/20 on1CPU (about5s), manifestae418440
  in asn_rotamer_preparation_v1. Root reviewed actual covalent/maps/fullsource
  contacts: omitted-source nearest heavy4.48–5.27A,H2.74–3.90A; local shortH
  contacts remain represented in the204atom scoring context, not repaired.
  Warm MACE4617 completed20/20,60baselinesreused,363wallseconds/7CPU. OnlyCa4
  changes: La selectsflip,pairedvacuumcontrast+1.460981253638387kcal/mol; other9
  selections unchanged. Both policies agree. Actual frozen1kcal stage gate passes.
  Read a0ac_asn_scoring/REPORT.md andRESULT.json. Matchedsolvent4623 completes
  all20new/60reused,80SCC+40nativecalls,475wallseconds/7CPU. Ca4commonpoolcontrast
  +10.974999511265196kcal/mol,La4+.6730578923597932,other8zero. Both sample4 La
  solvent gains are about9.5kcal; different vacuum headroom explains most of
  their different final gains. Read a0ac_asn_solvent/{REPORT,INTERPRETATION}.md.
  No class/band change. The old0.8
  domain is not relaxed or retrospectively passed. Root owns submissions.
  Saved-geometry4626 completes10/10: Ca4 AsnO-metal4.595→2.562A,N2.446→4.718A;
  correctlycalledLa4 showsanalogouschange. TargetincludedHD22-Gly292H1.382A
  remainsanunrelaxedcontact; read a0ac_asn_geometry_response/REPORT.md. Noenergycalls.

All scientific preparation/analysis/calculation runs on Slurm. Login-host work
is metadata, reading, editing and submission. Other users/sessions' jobs and the
production discriminator remain untouched. Current partial plans are not results.

Next finite research work:

- Native geometry/charge factorial4813 is complete:8 new off-diagonal calls,
 8 verified saved diagonals, all4 source/metal factorials available. No SCC,
 MACE or DFT. Both La solvent benefits persist at either saved charge vector;
 replacing charges opposes the geometry benefit. Large mixed terms are retained.
 Two Ca mixed-term aggregate precision flags remain false. Read
 a0ac_asn_solvent_factorial/ and asn_solvent_factorial_v1/ROOT_REVIEW.json.
 Actual cost4 scheduler seconds/1CPU, zeroGPU. No repeat or new score follows.
- Bounded flipped-origin preparation4821 passes10/10 sources,20 prepared tasks,
 zero energy calls. Original4814 remains partial5/10. Actual saved-cap4819 proves
 raw inputs/anchors/links and old/new computed caps identical; only pre-existing
 arithmetic residuals up to8.9e−16Å caused the raw-versus-computed exact check.
 Reviewed v2 uses three exact like-for-like comparisons and keeps both existing
 1e−12Å q0 checks. No molecular inputs/modes/tolerances changed. Root review is
 asn_relaxation_v2/ROOT_REVIEW.json. Read a0ac_asn_relaxation/PREPARATION_REPORT.md.
 Curie finishes the separate warm scout wrapper and Gauss reviews it; no
 optimizer submitted yet. Next: bind the actual v2 prepared first Ca0 pair,
 submit seven-CPU warm scout with watcher, inspect its actual geometry/energy
 admission and costs before remaining18 searches. Later solvent requires credible
 complete common6 pools and absolute common6−common4 change ≥1kcal/mol.
- Four MMOL proposal jobs4637–4640 are actually RUNNING; comparator4641 pending.
 All exact commands and completion/health watchers are in PROPOSAL_LAUNCH.json.
 No duplicate origin or proposal submissions.

Heisenberg owns focused independent review; root owns launch/score integration.
