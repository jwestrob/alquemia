# Conditional endpoint dynamics scout

Declared before any dynamics. This implements the next part of PLAN.md; it is
conditional on the periodic checks and real-box H200 benchmark. No dynamics
job is released by writing this file.

Question: can the primary L potential run stable, unrestrained dynamics of this
neutral system at useful speed, and does either prepared complex immediately
lose its bound arrangement? A short trajectory answers feasibility and reveals
gross preparation problems. It cannot establish an equilibrium population,
selectivity, or the accuracy of a force field.

Use start 0, endpoint A, as the first scout. Preserve its exact initial file.
Use the pinned native SO3LR command-line implementation, standard engine,
float64 model and state, highest matmul precision, charge 0, five unpaired
electrons, head 5, L-2606, real cutoff 12 A, dispersion damping 2 A and PME
(smearing 4 A, grid spacing 1 A, interpolation order 4). Do not constrain
fragment charges or atoms. Keep the 42 A cell fixed.

Preparation under that same potential: at most 200 native FIRE steps
(20 cycles of 10), default native FIRE timestep parameters, with an early stop
at maximum force 0.2 eV/A. This only relieves packing strain. Preserve the
optimization trajectory and report the final force; do not call it a physical
minimum or extend it automatically.

Then run 1,000 NVT steps at 298.15 K, timestep 0.25 fs: **0.25 ps**. Use the
native Nose–Hoover thermostat, chain length 3, two integration substeps,
`nhc_thermo=100` (nominal damping time 25 fs). Seed 54100. Write full-system,
float64 HDF5 positions/momenta/partial charges every 10 steps, retain the native
restart and complete log. The YAML boundary uses femtoseconds; the native loader
converts once to picoseconds. No restraints hold a metal on its ligand.

Bound this first dynamics scout by an internal one-hour timeout, including
compilation and preparation. Count its entire H200 allocation against the
one-GPU-day total, alongside the static benchmark and failures.

Inspect the preparation and trajectory before allocating the other starts:
finite positions/energies, temperature, covalent geometry and stereochemistry,
metal–donor distances, metal–ligand centroid separation, denticity, hydration,
water identity changes, proton transfers and counterion locations. Report
distances as well as counts; a counting cutoff is not an enforced coordination
number. Keep native full precision; any future cross-evaluated energies must
come from saved device arrays, not rounded log text.

The scout's measured time, rather than an assumed speed, will determine whether
the six independent endpoint preparations fit the remaining budget. Do not
extend this scout to chase equilibration. A failed initialization or loss of
binding is evidence to diagnose; neither should be silently turned into an
affinity result. No midpoint or alchemical ladder is automatically released.

7 November preparation amendment, before dynamics: the initial observations
found zero first-shell waters in every original box. HYDRATION_PREPARATION.md
declares one geometry-only repair; job53157 completed it and verified two/three/
two initial waters per metal across the three starts. The scout now uses
`hydrated_preparation_v1/start0/A.extxyz`. Ligand coordinates, cell, composition,
state, model, minimization bound and all dynamics settings remain as declared.
The original periodic benchmark retains its original inputs and results.
Analyze the eventual native trajectory with `trajectory.py --md-dir ...
--preparation workspaces/so3lr_neutral_pair_20261107/hydrated_preparation_v1/RESULT.json
--output-dir <fresh-workspace-directory>` in a one-CPU Slurm allocation.
