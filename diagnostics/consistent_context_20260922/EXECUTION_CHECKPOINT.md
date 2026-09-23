# Full transfer execution checkpoint

Historical running checkpoint. The completed result is in [REPORT.md](REPORT.md).

Frozen full scope: 225 transfer sources, 208 supported and 17 original unsupported.
There are 188 fresh MACE and 378 fresh GFN2 endpoints; 228 and 454 respective
components are reused with exact state and coordinate checks, including every
successful pilot component. The same canonical reference remains fixed.

MACE job 1210038 completed all 188 new endpoints in 23.248 s worker wall, including
1.281 s model load. Its allocation is one H200, 16 CPUs and 200000 MiB. The original
32-CPU job 1210032 was cancelled while pending, before any molecular work, to use
the 24 then-free H200 CPUs without changing other agents' jobs. The existing
warm worker accepts the declared CPU count. No scientific manifest, GPU/RAM,
checkpoint or method changed. The cancelled original collector 1210034 likewise
never started; its replacement is 1210039.

Native solvent job 1210033 uses 64 CPUs and 128 GiB on node-128-512g-8gpu-1, with
no GPU reservation. Collector 1210039 depends on both molecular jobs with afterany,
so failures remain visible. It writes collection_final_v1.json, comparison_v1.json
and autocollect_execution.json under the primary225_v1 workspace. It performs no
molecular calculations, parameter adjustments or score-dependent substitutions.

Exact commands, original/corrected submission receipts and resource changes are
recorded under workspaces/consistent_context_20260922/primary225_v1/. Full results
are pending at this checkpoint; no accuracy gain or production promotion claimed.
