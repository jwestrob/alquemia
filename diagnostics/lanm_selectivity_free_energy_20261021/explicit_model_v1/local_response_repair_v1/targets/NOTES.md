# Consumed targets for a local-response repair

`TARGET_INVENTORY.json` is a saved-data catalog, not a new loss or fit. Its
`artifacts` table holds absolute paths and SHA-256 hashes; each group identifies
the relevant JSON fields, literal coefficients, units and dependencies. The
model comparator is the frozen six-parameter selection from job17676. All four
target groups below have already been inspected and are development evidence.

| Group | Existing observations | Acceptance status |
| --- | --- | --- |
| Original joint water development | Ten energy checks and fifteen rigid-force fixtures on ten cells | Preserve every original energy error≤1 kcal/mol and force RMS≤2/max≤5 kcal/mol/Å gate |
| Pair/collective dipoles | 34 literal pair cycles; C/P/B for four parent clusters | Descriptive; no physical-error threshold or prospective loss assigned |
| Trimer nonadditivity | Eight E3 values and eight three-component μ3 vectors | Descriptive; same full-parent-basis seven-term cycles |
| Conditional susceptibility | Four state χzz values, two derived increments and their paired difference | Numerical consistency passed; physical-error threshold/loss remains undefined |

## Units and signs

Energy errors are model minus reference. Paired quantities are Dy minus La.
The original W8/W9 contrast is DyW9−LaW9−DyW8+LaW8. Radial work is displaced
minus central. Original reference forces are negative energy gradients. The
force helper uses ion net forces and each intact water's net force and centroid
torque divided by1Å, with its existing category/body/component weighting.

Original energy/force targets use the uncorrected ECP-PBE0-D4 surface. Trimer E3
uses electronic PBE0 without D4 in the full original parent basis. These are
different declared quantities; CP, CC, or trimer values must not silently replace
an original energy/gradient target. The small canonical-control failure and
remaining reference-method uncertainty are not removed by this inventory.

Dipoles use Debye in the original laboratory frames. Quantum moments include
nuclear and electronic parts with the recorded common-origin conversion;
effective metal ECP nuclear charge is11 and ghost charge is0. Individual pair
cycles are μ(MWi)−μ(M)−μ(Wi), with all three quantum terms in the same parent
basis. P sums these pairs, C=μ(MWn)−μ(M)−μ(Wn), and B=C−P. Preserve the N-fold
ion subtraction and the original water indices; no rotational replication or
substitution of an incompatible isolated-dimer geometry is implied.

The unchanged model's permanent moments cancel in these fixed-coordinate
cycles, allowing its admitted induced sums to represent the cycle's total
response. A new response functional must provide its own **energy-conjugate
total moment**. Do not assume that an old induced-array convention remains valid
after changing dual-field or environment-dependent terms.

Field data use laboratory z, ±0.001 and ±0.0005 atomic units; χ is in a0³.
χE=−[E(+s)+E(−s)−2E0]/s² and χμ=[μz(+s)−μz(−s)]/(2s).
The four states are W9 far(1,3) trimer and matching M+water1 dimer, La/Dy.
Their increments describe adding a second water; they are not α3. The existing
0.5 a0³ state,1 a0³ increment/paired, and1e−4 e·a0 odd-slope limits are numerical
consistency gates, not physical-equivalence margins. Printing bounds exclude
SCF, finite-basis and reference-method uncertainty.

## Avoid counting derived quantities as independent observations

- The ten energy checks have exact linear rank6. Paired radial works equal
  Dy−La individual works and displaced-minus-central balanced exchanges.
  All ten acceptance checks still remain required.
- Paired force vectors derive from the two individual metal vectors. Body
  projections and atomic residual diagnostics reuse the same atomic force sets.
- P is the sum of pair dipoles; C=P+B. Vector norms, xyz components and metal
  differences reuse those same moments. Near-zero W8 net vectors reflect
  symmetry cancellation and do not establish negligible pair response.
- The trimer cycles share monomer/dimer terms. Eight sampled E3 values cannot
  be summed into the full-shell B. μ3 and E3 are distinct derivative/value
  observations, but remeasuring μ3 through an odd field slope adds no new target.
- The28 susceptibility table rows are four states measured by two estimators
  at two steps plus derived increments and paired contrasts. They do not form
  28 independent physical targets. The four zero-field moments also reuse
  existing pair/trimer component data.

## Excluded validation

The four fresh W9 radial cells remain reserved and held after failed joint
development. Twelve donor-orientation and eight protein Stage B slots remain
reserved and excluded. Their numerical outcomes were not opened here. The four
earlier aquo8 orientation outcomes were already inspected: they are historical
validation/regression evidence, not fresh holdouts, and are also excluded from
this repair inventory. Exact IDs and authority records are in `exclusions`.

Inventory construction used one local thread,0.120s wall/0.112CPU-s and zero
native calls. It performed only saved-value extraction and coefficient sums;
raw admission was inherited from the pinned reviews. An initial metadata-key
error stopped before writing the inventory and is retained in
`BUILD_FAILURE_v1.json`; no scientific source or threshold changed.
