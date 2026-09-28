# Small isolated Dy: no converged endpoint

Both explicit starts on the archived 50-atom Hans EF3 core failed to produce an
admissible native r2SCAN-3c sextet endpoint. PModel terminated with explicit SCF
nonconvergence after 858 printed cycles. HCore remained above tolerance after a
near-root plateau and was stopped after the live-health alert. No intermediate
energy was accepted. This establishes that the large protein field is not
necessary for this numerical difficulty; changing both size and environment does
not identify a unique cause or prove which electronic state is correct.

Evidence: workspaces/metal_environment_response_20260926/dy_small_guess_v3/
contains endpoint outputs, ROOT_ACK_HEALTH_STOP.json, FINAL_COLLECTION.json,
health_receipt.json and completion_event.json. The generic collector says
PModel invalid and HCore missing because cancellation prevented its final
execution receipt. Neither means the endpoint was unattempted.

Job 1219790 was cancelled after 7010 seconds on 24 allocated CPUs; collector
1219791 completed in one second on one CPU: 168241 allocated CPU-seconds,
zero GPUs. Earlier startup-only failures/probes remain separately recorded and
are not included in that molecular-job subtotal. Allocated time is not utilization.

Next decision: one isolated PBE0 capability test, declared separately in
DY_PBE0_CAPABILITY_PLAN.md. No iteration increase, forced singlet, affinity
claim, large-region restart or production change follows from these failures.
