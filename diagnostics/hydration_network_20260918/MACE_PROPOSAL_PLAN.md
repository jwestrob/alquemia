# MACE as a water-orientation proposal model

Contained extension under Jacob's approved water-modeling work and standing
autonomy. No production change or additional quantum calculation. The purpose
is to find an affordable way to generate water configurations, with DFT still
judging their energies. This is not a new direct La/Ca classifier.

## Frozen experiment, before inference

Use each of the eight running full-network DFT searches: two structures, two
metals, two orientation starts. Copy the existing ORCA output into immutable
checkpoint artifacts and extract the first and third completed analytic-gradient
evaluations. Thus16 actual geometries, all from already executed DFT work. This
selection is independent of favorable energy, metal direction or classification.
The final optimized endpoint is not needed, and these iterates are not minima.
Printed coordinates have1e-6A resolution. Verify fixed atoms and rigid waters.

Evaluate unchanged native MACE-OMOL0-100M energies and analytic forces, float64,
actual net charge (La−1/Ca−2), multiplicity1, no categorical-charge masking,
no added dispersion, no dielectric or invented solvent correction. Use the
existing pinned checkpoint/environment and manifest executor. Request one GPU,
16CPU and64474MiB host RAM on node-128-512g-8gpu-1 (one eighth of host resources).
There are16 inference calls, no training or MACE optimization in this check.

Compare only differences within fixed metal/composition. Native charge-readout
constant offsets cancel within these comparisons. MACE's vacuum potential is
being tested empirically as a proposal proxy for the target CPCM-DFT surface;
its gradient is not called the gradient of that DFT Hamiltonian.

## Predeclared utility checks

- Correct sign for at least7/8 first-to-third DFT energy changes, counting only
  |DFT change|>=0.5kcal/mol; if fewer than8 are resolved, require>=87.5% and report
  the denominator. A zero MACE change is not a correct sign.
- Correct source/radial starting-seed ordering for all four metal/structure
  comparisons whose |DFT difference|>=0.5kcal/mol; report unresolved pairs.
- Median cosine>=0.5 between MACE and DFT rotational derivatives, over the16
  geometries. For each physical water, the derivative for an infinitesimal
  rotation is sum_H[(r_H−r_O) cross grad_H]. Concatenate all variable waters.
  Norms below1e-6kcal/mol/radian are undefined and explicitly excluded/reported.
- Report absolute errors in electronic energy changes, all torque norms/cosines,
  runtime/memory and disagreements. These are proposal-utility checks, not a
  claimed2kcal/mol free-energy certification or biological accuracy gate.

If it fails, do not use it to discard water states or replace DFT orientation
energies. If it passes, a separate recorded proposal trial can test whether its
suggestions actually lower target DFT energies. This is a development sample
from one biological group, with no new affinity labels or predictive validation.

## Pre-inference preparation finding

The first preparation stopped before inference: native redundant-internal DFT
optimization lets nominally frozen coordinates drift by up to6.3e-5A at the
third iterate. Four source-seed traces exceed the predeclared2e-5A geometry
check. Preserve that failed check; do not relax the criterion or call those
iterates exact water-only rotations. No completed optimized result is available.

The proposal comparison still uses exactly the same actual geometry for each
MACE/DFT pair, retaining every selected case and reporting frozen-coordinate and
water-internal errors. Its torque comparisons are well-defined physical water
rotations at those geometries. First-to-third energy changes include the small
other-atom movements and are labelled accordingly. This revision occurs before
any inference and does not change the utility thresholds. The original DFT final
constraint gate remains unchanged. `mace_proposal_v1` is an incomplete unexecuted
preparation; preparev2 preserves the audit and continues the matched proxy test.
