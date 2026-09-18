# Bounded response improves both alpha margins; neither direction is corrected

2026-09-18. The declared fixed0.20A metal-response model passes all native
energy/gradient checks. It reduces both alpha-lactalbumin versus GGR ordering
errors, but **0/2 primary directional comparisons pass**, unchanged from the
baseline. These are two structural observations of one alpha biological group,
not two independent validation proteins. All cases are consumed development.
Production/default, historical references and earlier experiments are unchanged.

| Alpha structure | Original alpha−GGR margin | Bounded-response margin | Margin improvement | Correct direction? |
|---|---:|---:|---:|---|
| 1F6S | −14.918692 | −8.300315 | +6.618377 | No |
| 6IP9 | −17.309145 | −6.073599 | +11.235545 | No |

Units are kcal/mol; comparator is GGR extended under the shared core policy.
Positive margin is the declared La-like alpha direction. These differences
cancel the common aquo offset; there is no new absolute reference or calibrated
classification. Improved margins must not be reported as improved accuracy.

## Actual conditional responses

| State | Movement,A | Actual DFT+J change,kcal/mol | Prediction minus actual,kcal/mol | Check |
|---|---:|---:|---:|---|
| Alpha1F6S Ca | 0.200000 | −13.224273 | +0.846736 | Pass |
| Alpha1F6S La | 0.200000 | −21.032509 | +1.565895 | Pass |
| Alpha6IP9 Ca | 0.200000 | −12.646982 | +0.691237 | Pass |
| Alpha6IP9 La | 0.200000 | −25.072386 | +0.761848 | Pass |
| GGR extended Ca | 0.057355 | −0.353519 | +0.001181 | Exact archive reuse/pass |
| GGR extended La | 0.081585 | −1.543378 | +0.037965 | Exact archive reuse/pass |
| GGR connected Ca | 0.048823 | −0.257748 | −0.002598 | Exact archive reuse/pass |
| GGR connected La | 0.084584 | −1.623583 | +0.022352 | Exact archive reuse/pass |

Actual Ca−La corrections: alpha1F6S +7.80823649314,alpha6IP9 +12.42540410768,
GGRextended +1.18985906045,GGRconnected +1.36583559959kcal/mol. The response
partition difference is0.17597653913kcal/mol, below the frozen2kcal requirement.
That checks the correction's representation sensitivity, not equality of the
underlying baseline scores across different carves.

All28 endpoint checks pass (16new alpha,12reused GGR). Alpha tangent-gradient
norms are0.472–5.162kcal/mol/A; all radial signs are consistent with the sphere
constraint. Energy errors0.691–1.566kcal/mol pass the original nonlinear-energy
scaled tolerances1.366–3.462. Full vectors, tolerances and components are retained.
The nonzero outward driving gradients mean these alpha points are **not free
stationary minima**. Neither the radius nor the acceptance criteria were changed.
All donor inventories, scaffold/cap coordinates, charges, waters and protonation
remain fixed. All eight coarse/fine matrices retain positive curvature.

The energy is the evaluated change in native r2SCAN-3c/CPCM core plus full-minus-
core short MACE interaction J at the prescribed bounded point. The inexpensive
MACE+GB curvature chooses that point; actual DFT energies/analytic gradients
independently check it. No numerical DFT gradients, quantum optimization,
conductor/GK correction, entropy term or additional solvation energy was added.
Long-range scaffold electrostatics and scaffold movement remain omitted; original
source-H limitations remain. This is a conditional response descriptor, not
binding free energy or a proof of physical equilibration.

## What this establishes

Metal motion carries a real, substantial differential signal in these alpha
structures: the La endpoint benefits more than Ca, consistently moving the
relative score toward the experimental direction. DFT confirms that signal at
the independently proposed MACE points. Small native discrepancies and the GGR
partition pass support using this response component in further development.

It does not establish that mechanical response alone explains the remaining
failure. The protein exterior stays fixed, the same structure is used for both
ions, and the endpoint model omits other binding/folding terms. The new paper
benchmark's RTX failure is a separate intact masked-MACE result; this pilot did
not compute a response correction for RTX, PQQ, parvalbumin or the other domains.
No blind-validation claim or general performance improvement follows.

**Recommendation:** retain baseline; keep this validated response component as
research support. Further development should test an independently declared
physical addition capable of addressing the remaining6–8kcal margins, rather
than expand this radius or refit a classification threshold to these cases.
No larger pilot or production rescore was launched in this step.

## Execution, costs and reproducibility

New calls:4native analytic DFT +8short MACE; zero failed scientific executions.
Reuse:4GGR native DFT +8short checks and all preceding curvature grids/centers.
Jobs1201383 and1201384 completed. Quantum used64CPU for217wall seconds;
short used one A5000/16CPU/64474MiB for81wall seconds. New total:
**15184allocatedcore-seconds,11930.002reportedCPU-seconds,81GPU-allocation-seconds**.
The allocations overlapped. This is incremental development cost, not an
end-to-end production benchmark; reuse does not make initial curvature free.

All three preparation attempts together recorded11.708221wall/10.754155CPU
seconds. Read-only report recorded2.016501wall/1.606895CPU seconds. Other local
housekeeping/tests were not fully profiled, not zero. Parent grid/GGR development
cost remains separately preserved:5730GPU-s,127328allocatedcore-s,
31251.232reportedCPU-s. No matched production affordability claim is established.

V1/V2 preparation remained unexecuted after frozen GPU preflight exposed an
unavailable preparation-only import. V3 performs full archived receipt checks
in CPU preparation/report and artifact/prediction checks in inference. Exact
scientific predictions/coordinates were preserved. No scientific retry occurred.
Four real-fixture tests pass in3.124s, none skipped: sphere/KKT/energy algebra and exact
reuse; corrupted real model rejection; unchanged native inputs/scaffold and
cache corruption rejection; actual DFT energy/gradient/sign/margin replay.

Protocol:`DFT_anchored_MACE_GB_bounded_metal_response_v1`.
Workspace:`workspaces/mace_bounded_response_20260918`.
Primary result:`report_v1/result.json`, locally pinned reporter.
Preparation:`prepared_v3/preparation.json`; costs:`cost_v1/result.json`.
Submission receipts hold exact argv and immutable manifests.
[Executable operations](BOUNDED_METAL_RESPONSE_COMMANDS.md).
