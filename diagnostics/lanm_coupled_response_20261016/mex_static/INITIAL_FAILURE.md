# Initial static assembly: stopped before molecular calls

Job 12673 failed during reuse validation: `Saved field raw replay differs:
mapped_physical_gradient`. The saved raw inputs and numerical results remain
intact. No new MACE, ALPB, GFN2 or MM evaluation was attempted.

The job used one scheduler second on 12 CPUs and one GPU; its application
reported 1.269 seconds. This is a software/replay failure, with no scientific
inference about the model. A saved-only diagnostic will identify the actual
difference before any changed validation policy or rerun. Original source,
config and traceback remain in `workspaces/lanm_coupled_response_20261016/
mex_static_v1/`; no old calculation will be repeated to resolve it.

An earlier review concern about distinct CPU/GPU interpreter paths was disproved
by the actual equal Python binary hashes. No source/runtime change was made for
that concern; it was not this job's failure.
