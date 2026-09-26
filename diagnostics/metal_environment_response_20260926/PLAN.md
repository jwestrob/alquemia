# Declared experiment and interaction ledger

Version 1, 26 September 2026; before new endpoint energies. The new handoff
authorizes this work. No historical whole-protein vacuum job is restarted.

## Scientific question and initial scope

Can a field-aware local model reproduce the **electronic component** of the
metal-dependent response to a real environmental displacement? Begin with
consumed 1H4I, complete archived qm36 donor/PQQ core. Source identities, boundary
redistribution and hydrogen inventory must come from one consistent archived
preparation. Do not mix the older chain charge −8 with the newer −9 preparation.

Environment A is the archived fixed-charge protein. Environment B rotates HG1
of the nearest intact exterior hydroxyl (THR159) +10 degrees about CB→OG1.
Selection is by proximity to core heavy atoms, not computed response. Preserve
heavy atoms, proton count, charges, O–H length and CB–O–H angle. Require no new
nonbonded H/heavy distance below 1.2 Å. This is a modeled perturbation, not an
experimental conformation or population. Preparation must confirm these gates.

Native reference: ORCA 6.1.1, native r2SCAN-3c, DefGrid3, TightSCF, EnGrad,
NoAutostart; native D4/gCP and La def2-ECP46, no substituted basis or composite
correction. Ca and La physical singlets, actual source charges recorded per
endpoint and checked against output electron count. No continuum term in this
finite-field diagnostic. Point charges in e/Å; potential zero at infinity.

Six initial endpoints: Ca/A, Ca/B, La/A, La/B, Ca/isolated, La/isolated.
Corresponding isolated endpoints are small-core diagnostics, not whole-protein
vacuum subtraction. Compatible previous endpoints may be reused only with exact
method/input/receipt correspondence. No molecular energy has been read for B.

## Energy ledger and explicit limitation

For this diagnostic evaluate

`E_el,M(r;R,Q) = min_rho [E_r2SCAN-3c,M(rho,r) + V_nuc-PC(r,R,Q) + V_e-PC(rho,R,Q)]`.

The local native energy includes its electronic/nuclear, D4 and gCP terms.
ORCA `%pointcharges` supplies the electronic and nuclear electrostatic coupling
exactly once. `%method DoEQ false end` excludes external–external Coulomb energy.
No separate charge–charge approximation is added to this embedded energy.
The local density responds self-consistently to the finite permanent field.
Movement of fixed MM charges is not self-consistent MM polarization.

The complete protein hybrid energy would additionally need MM internal energy,
boundary bonded terms and exclusions/1–4 conventions, and cross-region
repulsion/dispersion. Available standard-protein parameters do not establish
metal/PQQ cross terms or a matched capped-local subtraction. These are **missing
for a full hybrid**, not zero. We do not claim that metal-dependent missing cross
terms cancel, nor call this component test a complete environmental response,
mechanical force model, affinity or discriminator. No optimization uses it.

`delta_env,el = [E_el,La(B)-E_el,La(A)] - [E_el,Ca(B)-E_el,Ca(A)]`.

External-only terms would cancel in this double difference if truly identical
between metals; this algebra does not supply the missing metal-dependent terms.
Report each metal's response and the double difference in Hartree and kcal/mol,
converting once using the existing factor 627.509474 kcal/mol/Hartree.

Core `.engrad` and external `.pcgrad` are derivatives of this declared component.
They must be combined through actual cap/boundary maps before physical-coordinate
force claims. Synthetic cap coordinates are not independent protein motions.
Until checked, combined-force status is unavailable. The broader additive QM/MM
and MACEPOL engine paths remain separately gated by complete interaction coverage.

## Numerical gates, frozen before results

- Repeat energy difference ≤1e−7 Hartree; repeat component-gradient maximum
  difference ≤1e−6 Hartree/bohr. Native output must establish analytic gradients,
  intended ECP/electron count and converged SCF.
- Joint rigid-transform energy difference ≤1e−5 Hartree and maximum transformed
  gradient residual ≤1e−4 Hartree/bohr. These accommodate finite DFT quadrature;
  report residuals and rerun no chemical interpretation through a failed gate.
- Directional central differences at 0.001 and 0.0005 Å: error ≤0.05
  kcal/mol/Å plus 0.5% of the projected derivative; step-to-step agreement must
  meet the same limit. Initial checks will be declared as a finite follow-on,
  after native gradient availability, including an MM and mapped boundary mode.
- Constant-potential convention: adding c to phi changes energy by Q_core*c;
  no geometry force from a spatially constant shift. This is an algebra/interface
  test unless the actual backend exposes the scalar potential input.
- Candidate response accuracy target: absolute double-difference error ≤0.5
  kcal/mol and per-metal response error ≤0.5 kcal/mol; projection error ≤1
  kcal/mol/Å. Reference numerical uncertainty must be <0.05 kcal/mol and <0.05
  kcal/mol/Å respectively. A near-zero contrast within uncertainty has no robust
  sign. These are development targets, not evidence of biological accuracy.

## Execution gates and next stages

Pin and audit actual MACEPOL-EF model and interface before integration. Missing
paper checkpoint/engine inputs means no substitute checkpoint and no claimed ML
comparison. Continue the independent embedded-reference component harness.

Run only behind the other session's PQQ discriminator dependency chain; refresh
the queue immediately before submission. Use CPU allocations for native DFT,
actual MPI slots and allocation-derived memory. Record all attempts. No automatic
CPU job expansion while the scout is incomplete, no GPU reservation for absent
weights, no new functional sweep or automatic biological classification.

4MAE/non-PQQ expansion, expanded-region qualification, physical derivative
checks and then LanM EF3 each require the preceding relevant gates. Missing
full hybrid terms preclude coupled relaxation. All unused matrix cells remain
explicitly unavailable, including MACE cells if weights are absent.

## Primary implementation references

- [ORCA point charges](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html#inclusion-of-point-charges)
- [ORCA multiscale](https://www.faccts.de/docs/orca/6.1/manual/contents/multiscalesimulations/qmmm-molecules.html)

Retrieved exact pages are preserved under the matching workspace/reference_docs.
Upstream paper/code/weight findings are maintained separately in `upstream/`.
