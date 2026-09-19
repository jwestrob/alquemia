# Bound-water radial response: all four states pass

**DFT-anchored MACE predicts the tested water response accurately, including the
previously difficult 6IP9 arrangement.** This supports pursuing inexpensive
bound-water mechanics for the occupancy model. It does not yet provide occupancy
probabilities or a new classification result. Baseline/default/PQQ are unchanged.

## What ran

Jacob authorized with “proceed”; [frozen scope](AGREEMENT.md).
Use completed, consumed alpha-lactalbumin structures 1F6S and 6IP9, each with Ca/La:

- 1F6S 11: variable waters A211/A212, frozen outer A209; 70 atoms.
- 6IP9 110: variable waters A310/A322, frozen outer A329/A353; 73 atoms.
- Translate both variable waters rigidly along their own original metal–O vectors
 by a common ±0.05Å. Keep protein, metal, outer waters, orientation, internal water
 geometry, charge and protonation fixed. Ca/La have the same displacement mapping.
- Eight native r2SCAN-3c/CPCM(Water)/DefGrid3/TightSCF analytic energy-gradient
 evaluations; sixteen native unmasked MACE-OMOL0-100M float64 vacuum calls,
 including ±0.025Å refinement. Four center energies/gradients are reused.

The input O−H lengths are the earlier normalized 0.957162848Å geometry; inherited
`groups.OH_lengths_A` fields describe the original unnormalized source and are
not the current coordinates. The pinned XYZs and actual gradient coordinates
are authoritative. Water shapes are unchanged during these displacements.

New protocol: `native_r2scan3c_cpcm_anchored_mace_water_radial_motion_v1`.
No optimizer, numerical DFT gradient/Hessian, water insertion/deletion, new
reference, threshold fit, or production rescore was run.

## Prediction and result

For physical displacement h, g is dE/dh, with forces carrying the opposite sign:

```text
ΔE_pred(h) = h g_DFT(0) + E_MACE(h) − E_MACE(0) − h g_MACE(0)
g_pred(h)  = g_DFT(0) + g_MACE(h) − g_MACE(0)
```

DFT supplies the target energy and starting gradient. Vacuum MACE supplies the
nonlinear change. This is a local approximation to DFT/CPCM response, not a
self-consistent solvent model. No corrections are counted twice.

| State | DFT curvature | MACE curvature | Relative error | Largest energy error |
|---|---:|---:|---:|---:|
|1F6S 11 Ca|100.241|102.779|2.53%|0.00402|
|1F6S 11 La|236.801|238.182|0.58%|0.00464|
|6IP9 110 Ca|114.404|111.917|2.17%|0.00423|
|6IP9 110 La|228.653|230.245|0.70%|0.00707|

Curvature units:kcal/mol/Å²; energy errors:kcal/mol. These are secant curvatures
along the declared collective coordinate, not eigenvalues of a full Hessian.
All are positive. All four states pass all five frozen criteria: anchored energy,
even curvature, MACE refinement, analytic-gradient consistency and anchored
projected gradient. Largest projected-gradient error:0.377465kcal/mol/Å.
Coarse/fine even-energy disagreement:0.000301–0.000723kcal/mol.

| Structure | Motion | Actual Δ(E_Ca−E_La) | Predicted | Error |
|---|---|---:|---:|---:|
|1F6S|inward0.05Å|−1.10306|−1.10537|−0.00231|
|1F6S|outward0.05Å|+0.76166|+0.76686|+0.00521|
|6IP9|inward0.05Å|−0.62464|−0.62578|−0.00114|
|6IP9|outward0.05Å|+0.33901|+0.32996|−0.00905|

All entries:kcal/mol. The common aquo offset cancels in these displacement
contrasts. These expanded-context changes are not calibrated baseline decisions.
[Machine result](RESULT_v2.json), [table/CSV](export_v3/TABLE.md),
[figure](export_v3/energy_changes.svg).

## What this teaches us

The poor starting MACE force agreement for 6IP9 A322 did not prevent accurate
**collective radial force changes after a DFT anchor**. This distinction is
useful: cheap curvature can contribute even when its raw force is imperfect.
Because two waters move together, this test does not isolate A322's individual
response or exclude cancellation of errors between waters.

La is about twice as stiff as Ca along this coordinate in both tested structures.
Metal-dependent stiffness therefore exists in this local direction. Its entropy
contribution cannot be inferred from one coordinate, and neither its sign nor
its size establishes a discrimination improvement. These are four metal states
from two structures of one consumed biological group, not four independent tests.

**Recommendation: pursue the anchored water-response component.** Next extend
the cheap model to coupled physical translations/rotations of the retained
waters, preserving a common coordinate measure and all cross-water terms.
Check stable basins and nonzero gradients before estimating a partition integral;
validate representative soft/coupled directions with native DFT before using
free-energy terms. No next-stage calculations have been launched here.
Water internal vibration/solvent contributions and site-exclusion rules also
remain needed for an occupancy model. Unsupported terms stay unavailable.

This experiment is numerically credible for its tested coordinate, informative
about local mechanics, and inexpensive on the MACE side. A complete occupancy
model's accuracy and production cost remain unmeasured. Numerical entropy,
relaxation and occupancy outputs remain null (`response_model_not_validated`).

## Cost, checks and recovery

- DFT 1201953: 1044 s on 64 CPUs, four 16-rank workers;66,816allocated core-seconds.
 Eight endpoint wall times 472.793–563.049 s, median 514.629 s. Slurm peak batch-step
 RSS 12,116,232 KiB; 256 GiB requested. Executor elapsed 1041.064 s is separately retained.
- MACE 1201958: 136 s on one RTX A5000/16 CPUs;2,176allocated core-seconds.
 Peak worker RSS 1,553,904 KiB; peak allocated device memory 1.85724 GiB.
- Startup 1201954: 1 s/16 CPUs/one GPU; zero inference calls. Machine-precision vector
 recomputation differed by 2.22e−16 across NumPy builds. Exact input coordinates
 were preserved and a 1e−14 roundoff-only check fixed the startup; scientific
 tolerances were unchanged. The GPU environment dry-run then passed.
- Total: 69,008 allocated core-seconds and 137 allocated GPU-seconds including startup.
 Preparation, tests and plotting CPU were not individually metered. These are
 development-validation costs, not a measured production occupancy score cost.

The first preparation failed DFT path containment before submission. Identical
XYZ bytes were copied into the existing runner's supported directory layout;
no DFT result was repeated. Transient PMIX warnings in the second DFT batch did
not stop execution; all eight outputs/receipts confirm SCF convergence, normal
termination, complete gradients and return code 0.

**56 distinct real-fixture tests pass, zero skips:**26 earlier hydration tests,
seven motion tests and 23 archived-baseline/preparation tests. Scientific
integration is the actual eight DFT/sixteen MACE evaluations above. Tests cover
rigid geometry, coordinate/gradient mapping, force sign, paired displacements,
corrupted real manifests, technical cache invalidation, actual contrast algebra,
missing occupancy outputs and released PQQ energy/band regression.

The final result uses an immutable analyzer snapshot; its scientific contents
equal the initial RESULT.json. Final exports include the reused center explicitly.
Older exports remain retained; no energy or acceptance criterion changed.
[Receipts](RECEIPTS.json), [accounting](ACCOUNTING.txt), [commands](COMMANDS.md).
Next runnable operation is the documented collection/comparison replay; it
requires no new scientific calculation. Native orientation comparator jobs
1201824/1201825 are separate prior work.

Native orientation comparator1201824 has completed all four1F6S optimizations.
They converge normally but fail the frozen-coordinate tolerance (drift2.26e-5 to
7.46e-4Å); one also exceeds the rigid-water tolerance. Actual final energies and
gradient traces are retained in orientation_1f6s_v1/collected_opt_v2.json, without
promoting them to qualified exact-geometry minima. No rerun or tolerance change.
Comparator1201825 for6IP9 remains live; preserve it. These are separate older
runs and do not affect the successful exact-coordinate motion checks.
Job1201824 took10,572s on64CPUs in the earlier allocation; its676,608allocated
core-seconds are separate from the motion-pilot total above.
