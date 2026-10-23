# Isolated CPU metal-pair response: numerical qualification

## Question and scope

Does the opt-in direct/mutual cubic-Thole split implement a consistent energy,
induction solve, and Cartesian gradient while preserving the original engine
when disabled and preserving every unrelated pair kernel when enabled?
This is a numerical implementation test, not a fit or a physical-accuracy gate.
Root owns the isolated build, source audit, and one bounded execution. The
production executable and all completed results remain unchanged.

The proposed keyword is `POLAR-PAIR-RESPONSE typeA typeB delta_direct delta_mutual`.
Only unordered type pairs `{500,36}` and `{501,36}` are admitted initially.
At their existing `POLPAIR` coefficient gamma0, use
`gammaD = gamma0*(1+delta_direct)` and
`gammaM = gamma0*(1+delta_mutual)` with the original cubic radial exponent.
Finite deltas must exceed -1. Zero changes must route to the old kernel;
unmatched pairs, including metal–H and all water–water pairs, do so regardless
of whether a metal is present elsewhere. Do not set global `use_tholed`.
Intrinsic polarizabilities, local multipoles, exclusion/group scales, and all
other energy terms remain fixed. The implementation audit must cover both
native induced fields (`uind`/`uinp`), direct and mutual field builders,
preconditioner, polarization energy, Cartesian force, and local-frame torque.
The simple stationary-energy schematic is not a substitute for the native
dual-field/scaling derivation.

Initial supported execution is nonperiodic CPU MUTUAL polarization without
neighbor lists. Reject alternate global Thole, charge penetration, EXCHPOL,
Ewald, continuum solvent, OPT/TCG, and other unsupported routes with a nonzero
pair response. The source audit must verify this and the symmetric selector;
no new feature family is introduced by the harness.

## Frozen numerical probes and consumed fixtures

Same fixed probes for both metals: D-only `(0.05,0)`, M-only `(0,-0.05)`, and
combined `(0.05,-0.05)`. These are modest algorithm tests chosen before output;
no direction, magnitude, weight, or threshold may be selected from efficacy.
Also test absent keyword and explicit `(0,0)` records for both allowed pairs.
The pure-water fixture still receives these metal-pair records, so a record's
mere presence cannot silently switch global water physics.

Use nine exact frozen-17676 consumed fixtures, with unchanged typed graphs,
parameters, and local-frame definitions:

* `ion_La`, `ion_Dy` from the admitted 78-cell model manifest.
* `aquo9_central_La__near__two_water`: one intact six-atom water pair.
* `aquo9_central_{La,Dy}__dimer_01`: two four-atom metal–water dimers.
* `aquo9_central_{La,Dy}__{near,far}__metal_two_water`: four seven-atom trimers.

The existing geometry-only choices are near waters 1/5 and far waters 1/3;
do not replace symmetry-suppressed fixtures. Source manifests are
`saved_redistribution_counterpart_v1/recovery_preparation_v1/MANIFEST.json`
(`141048e19aed8337ead156d2aa9025fc546f8c2252cb64c23cb44f30a2239e92`)
and `collective_mechanics_options_v1/preparation_v1/MANIFEST.json`
(`a8fdf0c1ed724657bd642862d6af52632854e48ad8006f68cb0d87c2aef6f34c`).
No reference energies, reserved outcomes, or physical scores are inputs.

## One finite execution and fixed acceptance

1. On all nine fixtures, run original CPU analytic `testgrad system.xyz Y N`,
   then the isolated prototype with keyword absent and with explicit zeros.
   Require each prototype energy difference ≤1e-4 kcal/mol and every force
   component difference ≤1e-3 kcal/mol/Å versus the original. Preserve exact
   printed equality as an additional observation. Original analytic gradients
   are not available for every fragment, so these nine calls are counted.
2. Combined probe on all nine fixtures: bare-ion and pure-water E/F must retain
   those same identity tolerances. On the two dimers and four trimers use
   `testgrad system.xyz Y Y 0.00001`, with all metal/O/H coordinates active.
   This is the existing **total separation** 1e-5 Å, or ±5e-6 Å per component.
   Require all Cartesian analytic/numerical component differences ≤1e-3
   kcal/mol/Å. It exercises local-frame torque redistribution as well as
   metal and oxygen translation.
3. On both near trimers, run D-only and M-only with the same full Cartesian
   test. Their combined calculation from step 2 is reused. This separates
   direct and mutual implementation paths; no arbitrary minimum physical
   response or prescribed response sign is an acceptance criterion.
4. Supplement the combined near-trimer checks with two fixed paths per metal:
   (a) translate metal by `-u*q/2` and all three atoms of water 1 by `+u*q/2`,
   where u is its original metal→O unit vector; (b) rotate the intact first
   water around its O about normalized lab axis (1,2,3), keeping other sites
   fixed. For rotation let L be the norm of the original Cartesian rotational
   tangent and angle=q/L. Both q tangents have unit Cartesian norm, so slopes
   have units kcal/mol/Å. Evaluate energies at q=±1e-4 and ±5e-5 Å; compare
   `[E(+h)-E(-h)]/(2h)` with `gradient·tangent = -force·tangent`.
   Require error ≤1e-3 plus the known eight-decimal energy/gradient printing
   bound, and agreement between the two slopes within 1e-3 plus their two
   printing bounds. Energy rounding contributes `1e-8/(2h)`; the analytic
   projection contributes `5e-9*sum(abs(tangent))`. Use 17-digit displaced
   coordinates. Do not infer a physical error bar from these bounds.

Every program must exit normally and yield finite, complete arrays and every
reported final induction residual must be finite, nonnegative, and ≤1e-8 Debye.
Retain original polarization tolerance and iteration settings. Stop on the
first native, parser, identity, or derivative failure; preserve all nine
fixture slots and all unattempted program slots. Do not auto-retry or loosen
a setting. Actual parent/child CPU time, wall time, invocation counts, internal
FD counts, executable hashes, arguments, and raw outputs are retained.

The whole finite panel is **56 molecular program invocations**: 40 `testgrad`
(27 identity + 9 combined + 4 separate) and 16 directional `analyze E` calls.
Ten of the testgrad programs use Cartesian finite differences, totaling
**384 internal displaced energy evaluations** (216 combined +168 separate).
The 16 directional energy calls are separate invocations, not part of 384.
Serial, one CPU and one numerical thread is sufficient for these tiny cells.
There are no force-reference evaluations, new quantum calculations, fits,
iterations, or automatic continuation.

## Exact unrelated-pair fallback

Total water forces in a metal-containing system may change through the
induction solve; requiring them to remain unchanged would be incorrect.
The implementation audit instead binds the actual symmetric pair selector
and proves that unmatched pair direct/mutual damping factors **and their
radial derivatives** use the original branch. Inspect all O–O, O–H, H–H,
metal–H, and intrawater/excluded pair classes represented in each trimer,
both argument orders, and verify zero perturbations take the identical old
path. A small callable selector/kernel readback may support that source
proof; it must not evaluate an extra molecular panel or replace the full
gradient checks. Root binds this audit with the source and build pins in the
single runtime configuration. Absent/zero and pure-water numerical identity
provide complementary evidence, not a substitute for the in-system proof.

The prototype already prints native `PAIR RESPONSE TYPES/GAMMA0/D/M:` or
`PAIR RESPONSE ABSENT TYPES:` during parameter setup. Every enabled program
must provide exactly the two declared records. Actual represented metal/O
pairs must resolve gamma0 to their last frozen POLPAIR record and gammaD/M
to the prescribed profile (absolute/relative numerical allowance 1e-14);
absent types must remain absent. Absent-key programs must print neither
record. This verifies keyword recognition/application without requiring any
minimum physical energy response. Together with the bound selector audit it
requires no additional executable or molecular calls.

Passing establishes only this isolated CPU route at these small consumed
fixtures and perturbations. It does not qualify other engines, reciprocal
electrostatics, periodic salt, production dynamics, protein transfer, or
improvement of the model's failed physical comparisons.
