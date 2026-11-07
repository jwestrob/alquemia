# Neutral cis/trans Tauri continuation

**Continue the bounded feasibility test.** The checks completed so far support
energy cancellation in the intended neutral system. Chemical accuracy and
solvent dynamics remain untested. The next executable step is the already
queued H200 periodic benchmark; it measures speed and finishes the controls
before any trajectory is released.

The donor-count/charge model remains the established predictive baseline.
This continuation has not produced a new selectivity prediction. The missing
qualified quantum reference from the earlier screen remains unresolved; it
does not establish that Tauri failed.

## What these tests establish

All tests use the primary, unchanged L-2606 model. The proposed reaction swaps
La and Dy between trans600 and cis599 while preserving composition, total
charge, total spin and all coordinates. This is useful because donor counts
cannot distinguish these stereoisomers. It also avoids subtracting two systems
with different global charge/spin conditioning.

| Test | Why it is useful | Result so far |
|---|---|---|
| Add the same distant water to both endpoints | Detects a spurious change in the swap energy from remote surroundings | Maximum 0.00694 kcal/mol across three neutral preparations, two distances/orientations and the declared ligand motions; passes the frozen 0.5 criterion |
| Compare forces with an energy finite difference | Checks that the computed response matches the energy derivative on a real preparation | Residual 0.000000035 kcal/mol/A versus the declared 0.003 tolerance |
| Translate and wrap the entire periodic box | Represents the same physical configuration; exposes numerical dependence on its position in the cell | Swap-energy change -0.00000160 kcal/mol; maximum endpoint-energy change 0.000220 kcal/mol |
| Remove a remote water identically from both endpoints | Tests sensitivity to bulk surroundings in the actual full box | First water: swap-energy change -0.001032 kcal/mol at a closest solute distance of 27.04 A; two waters and the other starts remain pending |
| Measure the real box on H200 | Determines whether useful sampling fits the budget | Job 52844 pending, `QOSMaxGRESPerUser`; no GPU speed claim yet |

Energy cancellation does not eliminate force sensitivity. In the neutral finite
systems, the largest change in a force-difference component is 0.193 kcal/mol/A,
and the trans-metal vector changes by up to 0.173. Removing the first remote
water in the full periodic box changes the two metal force-difference vectors
by 0.070 and 0.159 kcal/mol/A. These responses are reported separately from the
energy criterion. Removing water also changes physical electrostatics, so a
response is not automatically a potential defect.

The current, unrelaxed start-0 endpoint gap is -59.411 kcal/mol. It is retained
for numerical reproduction and **is not an estimate of the experimental
selectivity contrast**. No comparison or fitted correction to experiment is made.

## System, settings and experimental status

Three preparations contain the exact curated cis/trans stereochemical graphs,
three archived initial ligand conformations, La/Dy, two K+ counterions and
2,400 waters: 7,288 atoms in a periodic 42 A cube. The ligand states are L4-,
with neutral tertiary amines. Preparation verifies relative stereochemistry,
covalent distances, exact A/B correspondence, charge and periodic contacts.
The cached restrained-MMFF poses supply starting coordinates only; they are not
measured structures or demonstrated minima under Tauri.

These boxes have **no added background salt**. They test feasibility and do not
represent the unverified 1.0 M experimental medium. Following CODE1's first-record
pointer recovers the correct bibliography identifiers: native599 maps to440,
and native600 to441. The cis La/Dy bibliography identifies
[Charlier, Merciny and Fuger, 1987](https://doi.org/10.1016/S0003-2670(00)85691-3).
Trans has multiple cited studies, including
[Merciny and Fuger, 1984](https://doi.org/10.1016/S0003-2670(00)84511-0).
Publisher metadata was retrieved; the original methods and tables were not.
The exact electrolyte and the source of trans's 1.0 M row remain to be verified.
Neither the primary pair nor the previously declared fallback is promoted to
an admitted common-medium target by this source recovery. No third pair is searched.

The potential and environment remain pinned to SO3LR70707f5 and GLP2294f4f,
Python3.12, JAX0.10.2, L-2606 weights SHA256
`0b5998bfd45929f344dd55174f542a4cfde6582458749477fcb6ad9e1d1d3d64`.
Use float64, x64 enabled before loading, highest matmul precision, charge0,
multiplicity6/five unpaired electrons, head5. Periodic settings: short cutoff5 A,
long cutoff12 A, dispersion damping2 A, PME smearing4 A, mesh spacing1 A and
interpolation order4. Native state propagation, finite float64 forces and
charge conservation are checked. Raw arrays retain full precision.

## Execution and cost

7 November preparation update: the new initial-observation analysis found no
water oxygens within3.5 A of either metal in any original box. That limits their
use as distinct hydrated starts. A single declared preparation repair adds
two/three/two first-shell waters per metal, respectively, using rigidly rotated
archived aquo-water coordinates and geometric clearance. It preserves all ligand
coordinates and repacks the remaining bulk water, keeping the total at2400.
These waters are fixed during packing only; dynamics will remain unrestrained.
The completed periodic results above refer to the original boxes. Their inputs
are preserved, including the pending benchmark's inputs.

Job53157 completed the repair and initial observations. All ligand bond and
stereochemistry checks pass; the closest chosen insertion clearance is1.881 A,
above the declared1.6 A minimum. Main-cutoff ligand O/N contacts at3.3 A are
trans/cis=5/2,5/5,3/4 across the starts, alongside2/3/2 waters per metal. These
are distance counts, not a forced coordination number or proof of binding.
The native model must still relax these preparations. Before that dynamics
scout, six additional fixed-geometry checks on the repaired start0 are declared
in HYDRATION_PREPARATION.md, using the same0.5-kcal/mol criterion. They are
unsubmitted and count within the existing one-GPU-day envelope.

The trajectory analyzer covers ligand integrity, stereochemistry, hydration,
coordination, metal–ligand separation, counterions and proton-transfer candidates.
It has been verified on real initial structures. Native HDF5 parsing and log
timing analysis await the first actual trajectory; no synthetic fixture is
represented as native-path validation.

All numerical work used Slurm `standard_expand`, with default scheduler memory
and time. The first 12-CPU periodic attempt exhausted its 48,000-MiB (46.875-GiB)
allocation. A resource-only 24-CPU recovery completed unchanged calculations,
peaking near52.4 GiB. The subsequent CPU controls reached their declared
330-second internal limit; their saved partial results remain usable. No
automatic CPU retry follows.

| Job | Work | Status | Allocated CPU seconds |
|---|---|---|---:|
| 52753 | Prepare three systems; build pinned Packmol20.15.0 | Completed | 81 |
| 52797 | 128 neutral finite-system evaluations | Completed | 1,320 |
| 52830 | First periodic attempt | Failed: memory | 816 |
| 52877 | Resource-only periodic recovery, two endpoints | Completed | 2,736 |
| 52948 | Periodic controls, four saved new predictions | Partial: internal timeout | 7,968 |
| 52991 | Saved-data analysis and bibliography checks | Completed | 1 |
| 53022 | Independent saved-data replay | Completed | 1 |
| 53123 | Analyze all six original endpoints | Completed | 7 |
| 53157 | Hydration repair and initial observations | Completed | 51 |
| 53167 | Repaired observation replay | Completed | 2 |
| 53205 | Combined saved-data and24-record observation replay | Completed | 5 |
| 52844 | H200 periodic controls and throughput | Pending | 0 so far |

This phase has used **12,988 allocated CPU seconds (3.61 CPU hours)** and
4,839.253 measured CPU seconds (1.34 CPU hours), with zero GPU time, zero DFT
and zero MD. These costs exclude the earlier, separately reported Tauri screen
and charged spectator test. Partial calculations without a saved complete
prediction are not counted as completed model calls.

The independent replay exactly reproduced the finite-system summary and both
CPU periodic summaries from the real saved arrays. It verifies hashes, element
changes, coordinate correspondence, precision and state records; it makes no
new molecular calls and uses no synthetic fixtures. Reproduce with:

```bash
bash diagnostics/so3lr_neutral_pair_20261107/reproduce.sh
```

The same command now also reproduces all24 original/repaired initial-observation
records; job53205 passed both this comparison and the earlier potential-table
checks exactly. HYDRATED_INPUTS.json pins the repaired input set separately.

Inputs, raw arrays, actual executed sources and complete logs are retained under
`workspaces/so3lr_neutral_pair_20261107/`; compact numerical results are in
CONTROL_SUMMARY.json and LOCAL_RESULTS.json. The active GPU completion watcher
is recorded in the workspace's WATCHERS.json. Preserve that job, other sessions'
executors, historical calculations, M1 and the unopened reserved outcomes.

## Next action and stopping point

**Run the queued H200 benchmark and review its full periodic controls.** If
the implementation is consistent and its cost is practical, the next contained
test is the unrestrained 0.25 ps endpoint scout specified in MD_SCOUT.md. Its
purpose is to expose unstable preparations, gross coordination problems and
actual integration cost. It cannot measure affinity or demonstrate sampling
convergence. The six endpoint starts remain bounded by one H200-day total,
including benchmarks and failures. No alchemical ladder or protein campaign is
released automatically.

Current decision: **continued feasibility testing is justified; selectivity
accuracy remains unresolved**. The live wait is a shared GPU allocation limit,
not a scientific rejection of Tauri.
