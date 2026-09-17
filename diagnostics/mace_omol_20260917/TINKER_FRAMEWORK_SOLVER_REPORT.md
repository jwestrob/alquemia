# Native coupled polarization solver: all controls pass

Job1201015 completed the declared twelve framework energy evaluations.
**All30 checks pass.** Native AMOEBA2018 electrostatics with coupled GK is
fast and numerically usable on these real protein/water frameworks. This is
implementation evidence; the metal is absent from these explicitly labeled
framework controls, so no La/Ca score or predictive improvement is established.
Production baseline/default remain unchanged.

| Framework | Atoms | Frozen source sites | Tightening change, kcal/mol | Rigid change, kcal/mol |
|---|---:|---:|---:|---:|
| GGR1GLG |4697|57|−0.001325844498|+0.000000069955|
| Alpha1F6S |1931|51|−0.000467879863|−0.000000016678|
| Alpha6IP9 |1897|54|−0.000570844031|+0.000000006897|

The four states per structure are all-polarizable, source-frozen, source-frozen
at tighter convergence, and a rigid transform of the latter. Physical atoms,
hydrogens, waters, disulfides, coordinates and permanent multipoles are preserved.
The adapter restores the source mask after native initialization, preserving
polarizability and damping parameters. All four induced-dipole arrays are
exactly zero on frozen sites. Standard solves took11iterations; tight solves16.
Every reported residual satisfies its1e−5 or1e−7Debye target. Each requested
energy produced one native SCF summary; no hidden extra energy requests.

Total energy equals native `em+ep+es` within the frozen1e−8kcal tolerance; every
disabled component is zero. Native `es` includes nonpolar cavity/dispersion
and must not be called pure reaction-field energy. Native GK solvent dielectric
is source-pinned78.3, internal1. All per-atom radii/scales and solvent constants
were exported and verified before energies. No grid, trajectories or optimization.

## Costs and validation

The12native kernels took **17.122408376seconds** summed wall time; complete
native subprocesses took26.024873793seconds and458.458163CPU-seconds.
GGR kernels took2.81–3.14seconds; alpha0.56–0.76seconds on64CPU node-64-768g-15.
The complete allocation lasted37seconds:2368allocated core-seconds,
467.119actualCPU-seconds, zeroGPU. Runtime child peak RSS414640KiB; scheduler
sampled batch MaxRSS224424KiB, a separately measured quantity.

Preparation took13.657238085wall-seconds,4.852405314parentCPU-seconds plus
7.911318childCPU-seconds for12parameter-only preflights. One-file frontend
compile took0.860898882wall/0.559939CPU-seconds. Earlier native-library build
and parameter export costs remain in TINKER_CAPABILITY_REPORT.md. Full hybrid
production-pair cost remains unmeasured; these are development controls.

Five tests pass in8.005seconds on the actual completed outputs: exact source
geometry/masks, damaged parameter inventory/mask, parameter-only output rejected
as energy evidence, and real energy parsing with damaged completion rejection.
The integration test was explicitly skipped while jobs were still incomplete,
then passed after all twelve outputs existed. No fabricated scientific output.

Protocol `Tinker26_2_AMOEBA2018_GK_framework_solver_control_v1`.
Workspace `workspaces/mace_omol_20260917/tinker_framework_solver_v1/`;
manifest SHA256 `ee7cadb5d9be64a70b0faf920f4f48213d4a9830982431e686afc3962347e2ed`.
Native source/library/frontend hashes, receipts, settings, unrounded energies,
response arrays, checks and costs are linked in TINKER_FRAMEWORK_SOLVER_RESULT.json.

## Judgment and next step

Numerical credibility: **passes these framework controls**. Scientific utility:
**untested**. Cost: **practical for this solver component**, full scorer unmeasured.
Retain baseline; pursue the hybrid. Next check whether the saved QM charge
representation reproduces electric fields that drive environmental induction.
Potential accuracy alone cannot answer that. Full source damping, covalent
boundary accounting, matched subtraction and predictive/partition checks remain.
