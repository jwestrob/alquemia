# Current causal investigation

The active user goal remains: "keep working until we get a conclusion and figure
out why our predictions are messed up in these cases."

The cavity exclusion branch is complete and closed; all callbacks through4587
were acknowledged. No repeat chemistry or wider cavity sweep is pending.
The nine-source causal synthesis in review/REPORT.md distinguishes established
solvent/calibration effects from source-specific unresolved mechanisms.

## Latest result — 15 October 2026, 04:35 UTC

Solvent4971 completed all40 new endpoints with80 old cells reused. The alternative
amide response survives the complete frozen-qvac composite: common3→6 paired
changes are +20.883598kcal/mol for Ca-conditioned4 and +19.290964 for La-conditioned4.
Other eight sources keep their old scores and selected geometries. Both selectors
agree. Each new sample4 minimum uses relaxed_flip_Ca for Ca and relaxed_flip_La
for La. This supports missed local orientation as a source of the inconsistency;
it does not yet establish a calibrated classifier improvement.

All90 saved algebra checks pass (largest residual6.65e-11kcal/mol). Actual160SCC
calls and80native component replays took941.196s internally; scheduler942s on7CPUs.
No new MACE/DFT/search occurred. Historical failed absolute-v1 parity is preserved;
the declared component-v2 and native replay gates determine admission. Read
asn_relaxation_solvent_v1/execution/FINAL_COLLECTION.json SHA256
64668c7ec758da42d34b039855bdc370b39d06815db47a00dfaf69c41a69368d.

Next: the same alternative-state rule must retain canonical Ca/La reference
fidelity. One-CPU geometry-only rebase5016 is submitted, watcher3932797, exact
manifest4d7a83d… and command in asn_canonical_relaxation_v1/SUBMISSION.json.
It cannot call molecular engines. Gauss owns the separate finite canonical search
adapter; Heisenberg owns saved solvent/range/old-band reporting. Curie owns a
separate fixed-heavy-geometry Asn-H bond-length sensitivity preparation. The latter
uses archived ff19SB targets to test a known generic H-builder artifact, without
changing the completed matched experiment. Root owns all submissions. Production
and thresholds remain unchanged; no new quantum campaign is launched.

### Preceding stage — 04:15 UTC

Full Asn searches4844 and contact review4910 are now complete. All20 searches,
40 final endpoint cells and ten common6 pools are available. Relative to the
fixed-flip common4 pool, vacuum paired contrasts change by +16.527363kcal/mol at
Ca4 and +11.463790 at La4; the other eight remain unchanged. Both mathematical
and operational selectors agree. Selected new proposals are inside the domain;
four other-source proposals retain boundary flags. See
a0ac_asn_relaxation/FINAL_MACE_REPORT.md and contact_review/FINAL_REPORT.md.

The specific HD22–Gly292H contact lengthens to2.17–2.28A in bothsample4 sources,
with oxygen-facing Asn retained. Other included H–H contacts remain short;
this is admitted finite geometry/contact relief, not an equilibrium claim.
The frozen sign-blind solvent gate passes. Actual preflight4968 replayed80prior
cells and20geometry/40final-frame mappings successfully with zero molecular
calls (61s/7CPUs). The unchanged matched solvent stage4971 is RUNNING on7CPUs,
watcher3871230. Configf9b8b4a36c1d1944251307bad16679f088b2a67e71e68428e9e86e03036d4ebe;
exact command/review/wake in asn_relaxation_solvent_v1/SUBMISSION.json.
This finite stage adds40endpoints across allten sources, reuses80old cells,
and adds no MACE/DFT/search or class/calibration change. Read its actual terminal
FINAL_COLLECTION.json and EXECUTION.json before choosing the next experiment.

MMOL's six native reference endpoints and comparison4641 are complete. Both DFT
and MACE select Ca/adaptive_Ca and La/adaptive_La. The common-pool differential
relaxation is −5.221985 versus −4.004941 kcal/mol: the Ca-directed response persists
and is stronger in DFT. This closes the finite MACE-only-artifact test, without
establishing a unique cause or a better classification. The Ca-path and full-pool
differences exceed the declared1kcal interpretive resolution. Density residuals
and nonuniform inverted frontier ordering remain electronic-state limitations.
See mmol_reference/FINAL_REPORT.md; authoritative result SHA256
cf0216391de0680b50b746c9799edd5365c23cc6ee1b9855e39844173903ab61.
Six molecular allocations total307500 allocatedCPU-s; no new solvent or forces.

Canonical fixed-amide geometry4936 completes25/25 sources and50/50 paired
geometries, including an exact real A0AC regression,11seconds/1CPU and zero
energy calls. Read a0ac_asn_relaxation/canonical_geometry/REPORT.md. Preparation
success does not authorize automatic canonical scoring or change calibration.

No exact compatible204atom A0AC Ca4 vacuumDFT cache was identified in the relevant
inventories; read a0ac_asn_inventory/CACHE_NOTE.md. No new DFT is submitted.

The chronological preparation notes below describe the stages preceding these
completed results. They do not authorize repeat submissions.

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
 Warm scout4831 completed the first Ca0 pair and both opposite-metal cells.
 Both searches converge without boundary flags; source/cap/bond/Jacobian and
 full-source severe-contact gates pass. Original common3/4/6 selections stay
 identical, paired change zero. Actual30 new forwards,570 scheduler seconds on
 seven CPUs. Read a0ac_asn_relaxation/SCOUT_REPORT.md and actual ROOT_REVIEW.json.
 Root released the exact remaining18 searches as4844 (RUNNING on7CPUs), watcher
 3586320. Confige1692709 and exact command are in asn_relaxation_remaining_v1/
 SUBMISSION.json. No solvent/DFT expansion or settings change. Complete common6
 and final named-contact analysis remain pending. Real4831 saved-contact check4859
 completed on one CPU/one second; configf7b1ddd2 in asn_relaxation_contact_scout_v1.
 All four tests pass, including actual nine unexecuted rows; zero energy calls.
 Ca0 contact analysis retains short H contacts in its unselected flipped proposals;
 it does not establish sample4 contact relief. Final40-geometry analysis binds
 only terminal4844 evidence. It is now queued as4910 afterany4844 on one CPU;
 watcher3688523 and exact launch are in asn_relaxation_contact_followup_v1/
 SUBMISSION.json. Reviewed AFTER.sbatch verifies qualified4859 source hashes,
 preserves partial/failed outcomes and records binding failures. Do not duplicate.
 See a0ac_asn_relaxation/contact_review/READY.md for exact binding command.
 Canonical25 source-graph inventory4847 completed on1CPU/one second. All25
 neutral existing ASN anchors pass identity/H-parent/bond/map/cap checks; some
 terminal modes are absent from their original four-mode selectors. Exact
 manifest8ba06c39 in asn_canonical_scope_v1; see canonical_scope/REPORT.md.
 Geometry4936 subsequently passed all25 sources/50 paired fixed flips, including
 exact A0AC regression; all new canonical energies remain unrun. Its manifest
 8271e29a and full pins are in asn_canonical_geometry_v2 and canonical_geometry/.
 Later solvent requires credible
 complete common6 pools and absolute common6−common4 change ≥1kcal/mol.
 The independently reviewed conditional adapter is ready at
 a0ac_asn_relaxation/solvent/READY.md: exact metadata binding, separate Slurm
 preflight-only command, then root-reviewed execution. It reuses80 prior cells
 and evaluates at most40 relaxed endpoints with unchanged qualified solvation.
 No solvent config is bound or job submitted while4844/contact results are
 incomplete; static review is not scientific qualification.
- MMOL4637–4640 and comparator4641 are all complete. Actual4638 completion and
 health callbacks are acknowledged. The last cell4639 (La/adaptive_Ca) ended with
 TRAH/NR residual5.899231e−6 below the printed1e−5 criterion, but large density
 residuals. Its final frontier ordering is ordinary; inverted ordering persists
 in both origins and other proposals. Electronic stability/same-state continuity
 remain unqualified. No post-hoc gate change or rerun. The complete response is
 reported above and in mmol_reference/FINAL_REPORT.md.
 The completed Ca proposal audit distinguishes4637 (ordinary final occupied/empty
 ordering) from4638 (inverted ordering and negative diagonal-Hessian warnings).
 Both meet their actual printed orbital-gradient criterion, with larger density
 changes retained. Thus electronic-state uncertainty is not uniform across the
 geometries; neither a passed native stop nor agreement alone proves a stable
 ground-state response. No new stability campaign or acceptance change follows.
 All exact commands and completion/health watchers are in PROPOSAL_LAUNCH.json.
 No duplicate origin or proposal submissions.

Heisenberg owns focused independent review; root owns launch/score integration.
