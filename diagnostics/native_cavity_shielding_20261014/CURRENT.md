# Origin cavity scout complete — job 4562

Jacob explicitly instructed execution. The finite12-update scope and numerical
checks are frozen in PLAN.md; no new permission gate applies.

Actual preparation4539 completed7s/1CPU on node-64-768g-16. All four source/metal
origins pass: local192 atoms, capped union9474, cap-removed9460, paired geometry
and exact charge/source maps preserved. Three scout rigid/permutation inputs are
prepared. No new SCC/MACE/DFT or energy call ran in preparation.

Thin-caller build4543 completed4s/1CPU against the unchanged qualified native
libraries. Worker4562 completed12/12 updates in136s with four shared CPUs;
all scout and per-endpoint gates pass. No duplicate chemistry ran. Read REPORT.md
and RESULT.json (saved summary4570). Ca-minus-La transfer is+5.031545/+4.936869
kcal/mol at control/ambiguous origins: primarily a common shift. Raw source
difference−0.094676 does not explain the ambiguity; its aggregate precision flag
is false under the frozen conservative bound ceiling. No tolerance repair/rerun.

Current config: workspaces/native_cavity_shielding_20261014/native_v1/CONFIG.json,
SHA f6620e7651f55b0df07c3dbc90de5db21468961da2e924928fc33f9158c594d7.
Exact submission and watcher command are in native_v1/SUBMISSION.json;
watcherPID3249783 queues event native-cavity-shielding-v1 to the actual root thread.
Inspect native_v1/execution/{SCOUT_QUALIFICATION,FINAL_COLLECTION,EXECUTION}.json.
Do not resubmit successful chemistry after a delayed callback.
Preparation4539/build4543 callbacks also require acknowledgement, not reruns.

No production score or calibration changes. Caps remain explicit model limitations;
this is a component exclusion sensitivity, not a source-only physical cavity.
The origin-shift explanation is closed. The separately declared next question is
whether burial changes accommodation work in the existing three-geometry pools:
../native_cavity_accommodation_20261014/PLAN.md. Eight saved nonorigin cells need
16native updates, zero electronic reruns; preparation/runner review is in progress.
