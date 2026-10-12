# 03_source_interpretation — collect the matched six origins first

Session UUID: 7e19f112-d807-461c-a981-5561c419c00f
Canonical repository: /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
Owned worktree: /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/claude_lanm_20261012/worktree_03
Owned output: /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/claude_lanm_20261012/source_interpretation_03_v1
Owned code/docs: diagnostics/claude_lanm_20261012/03_source_interpretation/ within your worktree.
Read COMMON_RULES.md. This dispatch supersedes old task03 execution scope and obsolete task02 ownership/status prose.

## First deliverable: one cheap, automatic six-cell comparison

Do this before a broad interpretation or mechanism survey. Root owns the currently running Hans8FNR La2665 and Dy2666 (12CPUs each). Read actual scheduler records and their actual manifests; never resubmit, cancel, modify or restart those jobs. PQQ/Ca-source work is separately owned and must remain untouched.

Use the exact canonical implementation:
/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/metal_environment_response_20260926/repaired_exchange_20260928/COMPARE.py

Create one NEW versioned comparison config by cloning:
/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/metal_environment_response_20260926/repaired_exchange_20260928/COMPARISON_RECOVERY_v1.json

Change ONLY the manifest/collection paths in the two Hans8FNR cells:
- La manifest: /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/claude_lanm_20261012/recovery_02_v1/La_singleton_origin_v1/manifest.json
- La collection: /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/claude_lanm_20261012/recovery_02_v1/La_singleton_origin_v1/COLLECTION_2665.json
- Dy manifest: /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/claude_lanm_20261012/recovery_02_v1/Dy_singleton_origin_v1/manifest.json
- Dy collection: /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/claude_lanm_20261012/recovery_02_v1/Dy_singleton_origin_v1/COLLECTION_2666.json

Keep all four admitted Hans8DQ2/Mex8FNS electronic origins, all six existing matching classical terms, source IDs, task IDs, source ordering, signs and physical limitations exactly. You may add separately named analysis-role/provenance metadata documenting the source declaration and interrupted historical attempts; preserve the original failure records. Retain each source independently. No refit, geometric change, new labels, averaging across Hans sources or reference offset.

Pin actual config, canonical COMPARE.py, scientific imported source dependencies, existing input manifests/classical/accepted collection files, and the two new endpoint manifests. New collections do not exist yet: record expected paths until actual execution, then retain actual hashes. Do not fabricate future pins or copy/relabel accepted scientific data. Use exact absolute paths because worktree workspaces are intentionally absent.

Create a task-local thin Slurm wrapper with an explicit real SLURM_JOB_ID guard BEFORE importing/running scientific code. Existing COMPARE.py remains untouched. Scientific imports, COMPARE.py execution, parser tests and scientific table generation are Slurm-only. Stdlib metadata/config edits and hash inspection are allowed on login; no numpy/import probes there. The historical16-thread login permission is expired. Root's actual Slurm2664 test receipt passed all8 corrected collector cases; do not rerun those tests as task03 work.

Submit ONE1CPU analysis/collection job on standard, normal priority, with afterany:2665:2666. Current QOS is12CPU/job,60CPU/user;4GB/CPU default RAM. Do not add memory/time/account/QOS flags or rely on whole-node resources. Both currently reviewed node-64-768g-1/-2 are valid choices if still healthy; use the actual scheduler. This finite job runs the existing comparison once, records scheduler accounting and actual output/error, and may export a small source-level table from that output. No ORCA, orca_2json, MACE, GFN2, OpenMM, geometry preparation/optimization, extra reference calculations, or new scientific branch. No call to old September sbatch scripts carrying obsolete memory/date settings.

Use afterany so failed endpoints remain explicit unavailable rows in the declared six-cell table. The endpoint wrappers already collect on EXIT. If a collection is absent after a hard failure, the comparison's missing-cell behavior must remain unavailable; do not invent an energy or silently substitute an old failed run. A parser/pin failure must remain a visible failed analysis with a diagnostic receipt. Do not repair scientific acceptance by weakening checks. Future collector fallback/recovery is root's decision, not part of this one-job dispatch.

Before yielding, arm the existing actual-job watcher on your comparison job:
/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/slurm_wakeup_20260927/watch_and_queue.py
with --jobs ACTUAL_JOB_ID, --thread 01a0a63a-ee36-7483-b726-1f7f5b7f75cc,
--codex /home/jwestrob/.nvm/versions/node/v22.23.1/bin/codex,
--receipt a fresh absolute receipt under your output, and a unique task03 event name.
Record exact command, PID and receipt. Root already watches2665/2666; do not create duplicate endpoint watches or a cluster-readiness monitor. Reconcile ambiguous submission with actual scheduler records before any retry.

## Required child work and stopping point

The installed delegation tool is Task. Invoke reviewer for the finite config/runner/ownership/Slurm guard and analyst for bounded archived mechanism/state caveats. They are read-only, may not calculate or submit. Do not wait for new science using model polling; read their finished results, write CHECKPOINT.md and STATUS.json, commit only owned task-local code/docs, and end with status waiting_slurm after submission/watcher. Include actual instructions_read, session UUID, config/source pins, job/dependency, source denominator, missing outcomes, exact next action, and child records. A human final under200words is enough. No push, merges, global changes, environment installs or production promotion.

If scheduling is blocked, report blocked_scheduler with actual error, keep runnable artifacts and stop. Do not replace queue access with login computation. A running job or successful LLM process is not an accepted scientific result.

## On root-authorized resume after comparison finishes

Read the actual six-cell output and raw receipts before conclusions. Deliver INTERPRETATION.md and the source-level table: both Hans-minus-Mex electronic contrasts, matching classical contributions, finite sums separately, complete/missing denominator, initialization/convergence assumptions, conditional EF3 occupancy, spectator identity, source waters, and region differences. These are consumed development controls, not blind predictions or binding free energies.

Read existing repaired compact/full-region, hydrogen-repair, Glu91 and classical-bridge reports for context. Distinguish geometric response from changed water inventories, spectator ions, and membership. Select at most two plausible mechanisms and one discriminating next experiment, with predicted distinguishing outcomes (no invented numbers), actual cost basis and rejection criteria, in EXPERIMENT_REQUEST.json. If either endpoint failed, report the failure and what can actually be concluded. No coupled-response or other next experiment may be scheduled until root interprets all actual evidence.
