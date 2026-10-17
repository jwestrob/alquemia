# Mex full-coordinate evaluator qualification

Declared 17 October 2026 UTC before new molecular calls. Root owns submission.

## Question

Does the existing coherent Mex electronic/MM/solvent energy and its analytic
physical gradient remain consistent when previously fixed protein, boundary,
water and spectator coordinates move? This is an interface qualification, not
proof of accurate mechanics or a new affinity result.

## Exact source and energy

Reuse the complete Mex common5 result from job 12731. The two starts are the
actual own-metal terminal POINT.json artifacts of La job 12727 and Dy job 12729.
`--stage bind` reads these receipts and writes immutable SHA-256 pins. Constructor
replays the archived source, original pilot, both searches and common5 references
without any new energy evaluation. The old q7 geometry guards apply only to
admitting those historical points. New evaluations accept all 2,049 physical
coordinates; archived q=[] means no reduced coordinate representation exists.

The energy remains exactly:

```
MACE vacuum(Q209) + GFN2 embedded(Q209;1840 charges) - GFN2 vacuum(Q209)
+ existing MM interaction ledger(R2049) + fixed FF charge ALPB(R2049)
```

All physical coordinates propagate to electronic atoms, all external charge
coordinates, both anchors of all four link atoms and their Jacobians, MM and
ALPB. The frozen preparation's `spectator_occupancy` is immutable identity/state
provenance; current positions are the physical-coordinate and mapping artifacts.
No environment charge is dynamically polarized. The target is La(III) or
Dy(III), with three fixed-identity Nd(III) spectators; spectator coordinates are
free in the new contract. Whole charge is +2, electronic-region charge zero.
MACE's audited `spin` input is **multiplicity** (1 for La, 6 for Dy). Original
spin-independent GFN2 uses effective multiplicity 1 / UHF 0, with 612 effective
valence electrons, separately retaining physical Dy multiplicity 6. This is an
inherited approximation, not a new spin validation.

Keep the exact installed standalone GFN2 executable, published parameter file,
runtime libraries, fixed-charge native ALPB binary, MACE weights/software and
MM ledgers pinned. No ORCA default parameter substitution is involved. The cache
key covers every physical coordinate, full model/protocol identity, state and
complete environment. No cache keyed only on the electronic XYZ is permitted.

## Finite matrix and predeclared acceptance

At each metal's actual own terminal:

1. Read/replay the exact cached energy and gradient twice (zero new calls).
2. Re-evaluate it once with cache bypass, retaining both receipts. These two
   total central calls are deliberate implementation/force checks, not two new
   independent molecular configurations. Every component work must agree within
   1e-6 kcal/mol and the full gradient within 1e-6 kcal/mol/angstrom.
3. Check four sign-independent, source-defined directions using central
   differences at h=0.0002 and h=0.0001 angstrom:
   - `core_backbone`: A/84/C along its C-to-O direction, exercising newly free
     peptide internal motion;
   - `boundary_anchor`: A/83/CA perpendicular to its C-to-CA vector in the
     CA-to-N plane, exercising link geometry and omitted-anchor derivatives;
   - `external_water`: rigid translation of the first lexicographically named
     complete, wholly external O/H/H water, directed away from the target;
   - `spectator`: first lexicographically named Nd spectator along its
     target-to-spectator vector.

Each full Cartesian direction has unit Euclidean norm. Water moves its three
atoms together, with each atom's displacement h/sqrt(3). Direction selection is
frozen before scoring either sign. This declares 4 directions x 2 metals x
2 step sizes x 2 signs = 32 displaced composite evaluations, plus the two fresh
central replays: **34 composite evaluations maximum**, no optimization or DFT.
Each composite uses one local MACE call, two local GFN2 calls, one complete MM
configuration/four component queries, and one fixed-charge ALPB call. MACE stays
warm; the existing pinned MM interpreter is a sequential per-cell subprocess.
The two single-thread GFN2 calculations may run concurrently within allocation.

For every direction, metal, step size and component (including their total):

```
abs(finite_difference - analytic_projection)
<= 0.005 kcal/mol/angstrom + 0.0001 * max(abs(FD), abs(analytic))
```

This is substantially tighter than the prior 0.2 + 0.05 relative test and resolves
the intended 0.10 maximum/0.03 atom-vector-RMS force convergence criteria. Both
step sizes must pass independently; no choosing the favorable step. Preserve
residuals and all unavailable slots. A failed check blocks full-coordinate Mex
optimization pending a concrete numerical or coupling diagnosis. It does not
permit relaxing these criteria after inspecting response signs.

## Checkpoint and optimizer contract

```
ev = MexEvaluator(pinned_config, fresh_output_directory)
state_ev = ev.bound("La")
point = state_ev.evaluate(ev.starts["La"])
work = state_ev.difference(new_point, old_point)
```

Returned values use kcal/mol and dE/dR in kcal/mol/angstrom; forces are -gradient.
`new_evaluations` counts actual total evaluations. Each request gets an explicit
Cartesian receipt referencing immutable component artifacts, full state, complete
coordinates and the complete identity hash. `restore_cartesian(receipt)` validates
this contract without reapplying q7 mode or source-displacement restrictions.
Optimization requires the separate full physical connectivity/stereochemistry/
contact and bound-state occupancy guards; the evaluator's domain guard only checks
finite arrays and nonsingular cap anchors. A dissociation decision is not hidden
inside the electronic evaluator.

## Commands and execution

`bind` only hashes existing metadata and may run on login. `preflight` and
`qualification` run with the pinned MACE Python under Slurm; root creates the
allocation and completion wake. No scientific work is authorized on login.

```
python diagnostics/lanm_coupled_mechanics_20261017/mex/MEX_CARTESIAN.py --stage bind --config /fresh/path/CONFIG.json
PINNED_MACE_PYTHON diagnostics/lanm_coupled_mechanics_20261017/mex/MEX_CARTESIAN.py --stage preflight --config /fresh/path/CONFIG.json --output /fresh/path/preflight
PINNED_MACE_PYTHON diagnostics/lanm_coupled_mechanics_20261017/mex/MEX_CARTESIAN.py --stage qualification --config /fresh/path/CONFIG.json --output /fresh/path/qualification
```

The real executable is recovered from the admitted old MACE config; no environment
upgrade. New qualification requires one allocated GPU and at least two CPUs, with
scheduler thread limits respected. Inspect successful zero-call preflight before
qualification. No automatic downstream optimization is launched by this script.
