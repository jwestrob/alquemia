# Initial DFT scout cost context

Read-only audit by Harvey, root reviewed17OctoberUTC. Current12977/12978 remain unchanged.

Reused original Mex references are La_A/Dy_A in workspaces/metal_environment_response_20260926/lady_repaired_origins_Mex8FNS_v1, scheduler1220336. They used112MPI each with MaxCore10800MB; two endpoints shared224CPUs/3094843MiB. La33cycles/SCF4390.360s/total1h28m53.155s; Dy43cycles/SCF5348.757s/total1h39m57.052s. Allocation6002s/224CPUs,1344448allocatedCPU-s.

Current scouts have12MPI each,MaxCore3100MB from48000MiB. Their firsttwo electronic iterations are1185.9/547.7s(plus) and1179.6/544.8s(minus), versus196.8/122.7s for oldLa. This is actual progress under fewer ranks/different host and memory, not evidence of a stall or a reliable completion prediction. Remaining six derivative endpoints stay gated on actual admission and measured cost.

Old wavefunctions exist at La_A/endpoint.runtime.gbw (127992317bytes,SHA c2fbbe1c9258ddb1e8405064771c35622dd399511e95954dbc14459ab330c434) and Dy_A/endpoint.runtime.gbw (127993533bytes,SHA32cfebc4940b0ccd40876fa4e15dead71cf46dee5734c2720323eb76d89cf97f). Current La coreXYZ is byte-identical, method/basis/state identical (3945basis/832explicit electrons/q0/effectiveM1); template removes EnGrad and externalPC line833 moves. Thus a future origin-orbital seed is available if a concrete measured repair is needed. It is not a converged displaced-state restart; GBW internals were not parsed. Current NoAutostart/GuessPModel intentionally does not use it. Any future seeded comparison must record initialization and retain native/density/state gates. No change or new seeded calculation is made by this audit.
