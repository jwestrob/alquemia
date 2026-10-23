# Finite metal-triggered exchange-polarization discriminator

## Question and frozen interpretation

Can a small anisotropic, short-range change in water's induced response produce a useful **collective** response beyond the six fitted radius/damping/epsilon directions? This is a directional functional-form diagnostic, with no optimization, new target, candidate selection or release to sampling. Preparation and saved controls use no molecular calls. Root owns native release.

The admitted 17676 point remains scientifically failed. Its literal fixed-geometry coordination error is −4.338284 kcal/mol: pair contribution +1.263316 and remainder −5.601600. The model remainder is polarization (+2.797568); its pair-additive vdW contribution cancels from the remainder. This does not identify the quantum remainder uniquely as polarization. The completed whole-shell CP shift does not rescue the full error; common-full-cluster-basis references can change attribution. Alpha probes and the consumed donor fit also failed their scientific gates. No held-out geometry or reference is used here.

## One supported mechanism

Use the existing CPU Tinker **EXCHPOL** term, preserving all selected 17676 geometry, radius, epsilon, POLPAIR, intrinsic POLARIZE, multipoles, bonded terms and water-water parameters. No engine edits. `kexpol` reads class-indexed spring/size/exponent/activity records; `kpolar` activates EXCHPOL from nonzero springs; `alterpol` builds and inverts anisotropic susceptibility matrices; `induce` solves the resulting mutual response; `epolar1` calls `dexpol` for radial and angular tensor derivatives. The energy analysis route remains `analyze → analysis → epolar3 → induce`.

For SI water class 34, set a passive positive spring, size 1; class 35 H has zero spring and size, passive. Both metal classes 500/501 have zero spring, size 1 and are active. All three kinds use native exponent **b=4 Å⁻¹**. The pair activation rule is logical OR of activity flags: M–O changes water O; M–H is zero; O–O/H pairs skip. The metal's own tensor stays identity. No O-only or metal-only polarizability replacement occurs.

Select native `EXCHANGE-POLAR S2U`. For equal exponents b, its exact overlap is

`t = b*r/2; S²(r) = [(1+t+t²/3)*exp(−t)]²`.

Two fixed amplitudes, shared across La and Dy: **q=0.01 and 0.05**. Render the O spring as `k=q/S²(2.50 Å)`. Thus each single-metal water-O matrix is `I + q*S²(r)/S²(2.50)*u*uᵀ`; at 2.50 Å its radial inverse susceptibility increases by q (dimensionless), and radial effective alpha is divided by 1+q. Tangential susceptibility and intrinsic alpha stay unchanged. b is a diagnostic native overlap exponent, not a fitted or published La/Dy parameter. No exponent/amplitude selection after outcomes. Only positive semidefinite additions; all source-reconstructed site tensors must be positive definite. This does not prove positive definiteness of the full mutual-interaction operator.

The overlap model and anisotropic variable-polarizability mechanism are described by [Chung et al., JPCB 126, 7579–7594 (2022)](https://dasher.wustl.edu/ponder/papers/jpcb-126-7579-22.pdf). That work uses a HIPPO context, so this restricted AMOEBA perturbation is a hypothesis, not transfer of a validated parameterization. The [native keyword documentation](https://dasher.wustl.edu/tinker/distribution/doc/sphinx/tinker/_build/html/text/keywords.html) defines active/passive EXCHPOL records.

## Exact finite inventory and order

1. **Numerical pre-gates, 4 cells / 8 programs:** zero-q W9 La and Dy (analyze + analytic testgrad), E/F replay against admitted selected outputs under the inherited engine tolerances below; one exact metal-deleted 27-atom W9 water shell at q0 and q.05 (each analyze + analytic testgrad), E/F equality under the inherited engine tolerances below. Waters retain exact positions, types and bonds. No isolated-shell reference target is introduced.
2. **Positive-q derivative gate, 4 cells / 8 programs:** q.05 W8/W9 × La/Dy, analyze + `testgrad system.xyz Y Y 0.00001`. Full analytic/numerical Cartesian error ≤1e−3 kcal/mol/Å. Exactly `6*(25+25+28+28)=636` internal displaced-energy evaluations. All four must finish and pass before step 3.
3. **Remaining fixed matrix, 84 cells / 100 programs:** both amplitudes on all ten consumed water fixtures (analytic E/F), and all 34 literal matched dimers (energy only). The four q.05 cluster rows from step 2 are reused exactly. No rotational replication or new geometry.

Total **92 distinct program-input cells, at most 116 native programs, 636 internal FD energies**, four single-thread workers in one four-CPU Slurm task, default scheduler memory/time. Thread environments all 1. Fixed denominator retained on failure; first native/parse/admission failure prevents new program starts, drains siblings and preserves receipts. No retries or continuation. `POLAR-PRINT` is appended solely to record induction convergence: every printed final residual must be finite, nonnegative and ≤1e−8 Debye (native output rounds to 1e−10). Existing polar-eps/iteration controls remain unchanged.

## Saved reporting

For each amplitude retain all ten original scores/25 controls without retuning, per-cell components, C/P/B and component decompositions, and changes from the frozen baseline. C is Dy9−La9−Dy8+La8; P uses all literal dimer sums with the existing analytic zero isolated-model-ion convention; B=C−P. Report `(response(q=.05)/.05)−(response(q=.01)/.01)` as scale dependence, not extrapolated improvement. No best amplitude or candidate is selected. Report actual program/internal-FD count, wall/aggregate child CPU and full missing slots. Root's common-basis reference analysis remains separate.

## Limits and compatibility

This term also changes isolated M–water dimers. It is not a B-only correction or explicit coordination-count model. Native source-reconstructed off-diagonal tensors are recorded; the native executable does not print those tensors. Native residuals and Cartesian checks are independent runtime evidence. All represented sites must have strictly positive intrinsic alpha because the native derivative implementation divides spring by alpha even for inactive sites. No whole-protein zero-alpha support is assumed.

Metal-absent water retains its existing physical Hamiltonian, including the SI/BIO Urey–Bradley distinction. A later BIO mapping would target only water class 73 and metal 900/901; it is not implemented or qualified here. Donor interactions, protein multipoles/bonds and hydration training are not changed. CPU EXCHPOL support does not establish GPU/periodic parity or production suitability. A useful response would justify a separately declared model-development question, not establish unique missing physics or reference truth.

## Prospective numerical clarification before native results

The initial exact-invariance source/spec and preparation_v1 are preserved. Before any EXCHPOL output, root and reviewer agreed to inherited engine replay tolerances: |ΔE|≤1e−4 kcal/mol and maximum Cartesian force-component |ΔF|≤1e−3 kcal/mol/Å for both zero-q and metal-absent controls. Exact printed-equality flags remain descriptive. Identity EXCHPOL tensors can still use different floating-point solver algebra. This changes no scientific score gate, FD tolerance or induction tolerance.
