# Solvent-consistent physical accommodation response

Jacob authorized the overnight goal of finding transferable accommodation
information beyond one frozen-structure score. Parent explicitly approved seven
consumed contexts, 28 native GFN2 analytic-gradient calls, reuse of 14 exact OMOL
force arrays, and small matched physical-displacement checks. No label fit,
relaxation/entropy scalar, new biological case or production change is authorized
by this pilot. This plan precedes all new scientific outputs.

## Question and distinction from prior work

Does the solvent contribution materially change metal-dependent forces on whole
donor torsions/peptides, where vacuum model proposals and atomwise radial averages
were inconsistent? Different Hamiltonians guided and judged the old coordination
proposals; all ten touched the displacement boundary. Earlier GGR radial means
also mixed changing donor identity and mutually opposing atomic loads. Retain
physical per-mode derivatives without fitting their weights or averaging away
their signs.

Use the existing compact composite at each exact source endpoint:

`E_M(q) = E_OMOL,vac(x_M(q)) + E_GFN2,ALPB(x_M(q)) − E_GFN2,vac(x_M(q))`.

`dE_M/dq = J_M(q)^T [g_OMOL + g_GFN2,ALPB − g_GFN2,vac]`.

The reported differential derivative is `dE_Ca/dq − dE_La/dq`. Translational
components have units kcal/mol/Å; torsional components kcal/mol/radian. A force
has the opposite sign. Keep vacuum and solvent terms separately. This is the
derivative of the composite descriptor, not a DFT gradient or binding free
energy. Electronic smearing remains 300 K; no nuclear thermal interpretation.

## Exact sources and fixed coordinates

Contexts, in order: 1H4I, 4MAE, 1F6S, 6IP9, 1GLG, 2FW0, 2FVY from
`workspaces/compact_solvation_20260920/full_v1/manifest.json`. Use both Ca/La
source coordinates, charges, singlets, atom order, complete original context,
water inventories and exact native OMOL0-100M checkpoint. All14 force arrays
already exist; no GPU recomputation is planned.

Reuse `coordination_preparation_context.geometry` and `mace_site_kinematics`:
metal translations, original donor sidechain chi rotations and actual peptide
crankshafts, with caps following both source bond anchors. PQQ, added neighbors,
scaffold and water geometry remain fixed. Rebuild mappings around the *current
compact* source coordinates, not the different PQQ H-prepared optimization starts.
No independent fragment translation or cap coordinate is introduced. Alpha's
endpoint-specific prepared water H stay distinct and frozen; the common donor
coordinates define the comparison, not a false identical-water Cartesian point.

## Finite manifest: 60 new low-level evaluations

1. **Analytic support gate:** four native GFN2 EnGrad calls for 1H4I,
   Ca/La × vacuum/ALPB(Water). Use the native mixer, installed parameters,
   `Convergence Tight` (actual TolE ≤1e−8 Eh), NoAutostart and unchanged300K.
   Require actual converged analytic gradient artifacts and no numerical-gradient
   fallback. Verify total-energy agreement with original SPs and the previously
   frozen ≤0.10 kcal/mol endpoint-transfer / ≤0.20 pair-difference gates.
2. **Remaining centers:**24 identical-policy EnGrad calls for the other six sites.
3. **Derivative verification:**32 GFN2 single points. Before gradients are known,
   fix two geometrically supported modes:1H4I A/Glu177 chi3 (whole carboxylate
   rotation around CG–CD), and1GLG A/Gln140 peptide crankshaft. Evaluate each at
   ±0.005 and±0.010 radians, both metals and both media. These are energy-based
   checks of the analytic low-level derivatives, not numerical DFT gradients.
   All protein covalent bond lengths and frozen atoms must remain intact; the
   maximum physical heavy displacement must remain≤0.05 Å. Unsupported mapping
   means explicit failure, not score-dependent replacement of a mode.

The support gate runs first. Remaining centers and displacement tasks launch
only if native analytic support and same-state energy accounting pass. No new
DFT, MACE, optimizer, Hessian, trajectory or model training is included.

## Checks fixed before outputs

- Exact source hashes, atom order, charge/spin/electron counts, native GFN2
  parameter export, ALPB selection and unchanged300K. Charge closure≤5e−4 e,
  max atomic charge≤4 e remain the existing diagnostics.
- Gradient energy matches the same output within1e−8 Eh and gradient coordinates
  within1e−6 Å; reject numerical differentiation in output or input.
- Each requested physical mapping passes existing analytic Jacobian/bond/frozen
  atom checks. Ca/La use identical mode identities; synthetic caps are derived.
- At h=0.005 rad, central-energy derivative must agree with the analytic projected
  derivative to `max(0.02 kcal/mol/rad, 0.005*abs(derivative))`; the two h estimates
  must agree to the same tolerance. The absolute floor exceeds the conservative
  1e−8-Eh energy stopping noise divided by h, while the0.5% relative bound checks
  meaningful large forces. Apply to each medium and the ALPB−vacuum difference.
- Report all failed calculations and checks; do not silently retain vacuum
  forces when solvent gradients fail. No numerical corrected affinity appears.

Compare archived DFT gradients only when exact coordinates, atom order and state
match. The older POLAR/DFT force table is separately replayed as historical
mechanism evidence; it is not a validation of the new full-context composite.

## Interpretation and resources

Report the complete per-mode Ca/La and difference vectors, solvent-induced sign
changes, mode-specific replica sensitivity, derivative-check errors and actual
cost. No fitted classifier or scalar force norm is selected using labels.
Useful solvent response would justify a separately frozen accommodation test;
this pilot alone cannot establish improved biological accuracy.

Existing ORCA6.1.1 runner,64 CPUs/128GiB, eight concurrent eight-rank tasks, no GPU.
Prior primary context calls averaged roughly10–20 s each;60 calls predict minutes
of allocation, with actual gradient overhead measured. No project compute/time
budget stop condition; existing scheduler/convergence policies remain.

The official [6.1 changelog](https://www.faccts.de/docs/orca/6.1/manual/contents/appendix/detailedchangelog.html#native-xtb-methods)
documents native GFN2 energies/gradients and native ALPB. The
[native method section](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html#native-gfn-xtb-and-gfn2-xtb)
documents mixer, parameter and smearing controls. No local completed native GFN2
EnGrad was found, so actual combined analytic support remains the first gate.
