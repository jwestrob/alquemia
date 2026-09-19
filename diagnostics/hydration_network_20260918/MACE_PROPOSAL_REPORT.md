# Native MACE passes the water-orientation proposal check

Job1201847 completed16 actual energy/analytic-force evaluations. No new DFT ran.
All predeclared proposal checks pass on the consumed alpha-lactalbumin geometry
sample, one biological group:

- Eight/eight first-to-third DFT energy-change directions reproduced.
- Four/four source-versus-radial starting orientation rankings reproduced.
- Median rotational-gradient cosine0.991329 across16 geometries.
- Mean absolute error in the eight local energy changes0.424043kcal/mol;
  maximum0.875042. The four larger seed-switch errors range1.69–9.54kcal/mol:
  this supports proposal directions, not substituting MACE absolute energies.

The unchanged native vacuum OMOL model is empirically useful for these water
motions despite its solvent mismatch with target CPCM-DFT. No masked model,
charge fit, water deletion or biological threshold was selected from the results.
No occupancy, affinity classification or generalization claim follows yet.

Actual cost:137s elapsed allocation,16CPUs (2,192allocated core-seconds), oneGPU
(137allocated GPU-seconds). Median per-geometry inference0.811408s including
first-call overhead; model/process startup was repeated by the existing runner.
Peak framework allocated GPU memory1.99124GiB; max per-worker hostRSS1,568,696KiB.
Requested host memory64474MiB, the existing one-eighth resource policy.

Results/receipts:
`workspaces/hydration_network_20260918/mace_proposal_v2/collection_job_1201847.json`.
Plan and pre-inference geometry-drift finding: [MACE_PROPOSAL_PLAN.md](MACE_PROPOSAL_PLAN.md).
The source/radial DFT starts shift R by24.824611/31.161834kcal/mol in1F6S/6IP9;
restoring neighbors plus tighter SCF changes source-orientation R by only
1.049171/1.533610. Those are completed SCF/gradient checkpoints of still-running
optimizations, not final minima. See [initial checkpoint results](INITIAL_CHECKPOINT_RESULT.json).

In both structures the previously problematic water (A211/A310) has a generated
H only1.852/1.880A from the metal in the normalized source seed. Its water
bisector points partly toward the metal. The label-independent radial seed puts
water H farther away and lowers both endpoint energies, especially La. This is
a concrete preparation/geometry mechanism worth testing to completion; it does
not prove the final occupancy or isolate all electrostatic contributions.

Next: use MACE to optimize the same two starts with exact rigid-water rotations,
then score each metal's independently selected MACE minimum with native DFT.
Keep ongoing native DFT searches as the direct comparator. No production change.
