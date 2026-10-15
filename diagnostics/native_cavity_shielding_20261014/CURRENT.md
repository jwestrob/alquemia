# Native cavity scout running — job 4562

Jacob explicitly instructed execution. The finite12-update scope and numerical
checks are frozen in PLAN.md; no new permission gate applies.

Actual preparation4539 completed7s/1CPU on node-64-768g-16. All four source/metal
origins pass: local192 atoms, capped union9474, cap-removed9460, paired geometry
and exact charge/source maps preserved. Three scout rigid/permutation inputs are
prepared. No new SCC/MACE/DFT or energy call ran in preparation.

Thin-caller build4543 completed4s/1CPU against the unchanged qualified native
libraries. The independently reviewed runner is submitted as4562 and was verified
RUNNING on node-64-768g-16, four shared CPUs with16000MiB scheduler-granted RAM.
The qualified engine is serial; slots provide memory, not four-way computation.
The job runs six sample0Ca numerical controls, then six further updates only if
the frozen gates pass. There is no separate first-two scout or duplicate chemistry.

Current config: workspaces/native_cavity_shielding_20261014/native_v1/CONFIG.json,
SHA f6620e7651f55b0df07c3dbc90de5db21468961da2e924928fc33f9158c594d7.
Exact submission and watcher command are in native_v1/SUBMISSION.json;
watcherPID3249783 queues event native-cavity-shielding-v1 to the actual root thread.
Inspect native_v1/execution/{SCOUT_QUALIFICATION,FINAL_COLLECTION,EXECUTION}.json
when available. Do not resubmit successful chemistry after a delayed callback.
Preparation4539/build4543 callbacks also require acknowledgement, not reruns.

No production score or calibration changes. Caps remain explicit model limitations;
this is a component exclusion sensitivity, not a source-only physical cavity.
