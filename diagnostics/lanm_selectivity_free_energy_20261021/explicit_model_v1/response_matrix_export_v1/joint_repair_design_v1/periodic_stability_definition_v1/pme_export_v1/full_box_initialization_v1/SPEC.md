# Full-box PME initialization: operator-only source route

Source only. No key has been derived from an actual candidate, no full box has
been initialized, and no program has been linked or executed. Final candidate,
composition, fixture selection, runtime and release bindings remain null. The
existing 49.324 Å BIO La/Dy fixtures remain unchanged.

## Minimal native path

`INIT_PME.f90` is a separate program linked to the qualified exporter library;
it does not patch the library. Its path is:

`initial → getxyz → mechanic → chkpole → rotpole('MPOLE') → localresponsesetup`

followed by the exact `dfield0c` assignments `nfft=nefft`,
`bsorder=bsporder`, `aewald=apewald`, and **one `udirect1` call**. The latter
performs native getchunk/moduli/bspline_fill/table_fill, permanent scatter,
forward FFT, qfac construction and the already qualified export, followed by
the original backward FFT/gather. It returns only a temporary reciprocal
permanent field, which is checked finite and discarded. It does not call the
real-space direct field, induce, an SCF solver, energy, gradient, force or MD.
There is no dense global response construction. One future bound fixture means
one program and one reciprocal initialization, with zero induced solves.

The normal `localresponsesetup` is still called. Both response guards remain
fatal in the unchanged library. This program rejects response records before
setup and requires lractive/lrperiodic false: the absence is explicit and is
never presented as execution of the composed response model.

## Exact derivative input and future join

The prospective composition adapter emits two translated response records,
900/901–402. Pure `DERIVE.derive` removes exactly those two full lines, retaining
their literal bytes and indices so the original key reconstructs exactly.
Everything else in the composed key, including candidate metal–O R/epsilon,
fixed POLPAIR, intrinsic alpha, masses, permanent moments and BIO terms, remains
literal. The transform appends only the declared initialization controls and
`PME-NORM-EXPORT 1`. Unknown preexisting periodic/field/response controls are
rejected rather than overridden. XYZ and the complete BIO parameter file must
remain byte-identical to the bound composition and original fixture.

Before any preparation/execution, the caller must verify the future terminal
candidate/configuration/actual-review and exact composition manifest, together
with the frozen composition source and static input inventory. The derivative
manifest must join those pins, removed records, reconstructed original key,
new key hash, unchanged XYZ/BIO hashes and separately pinned native executable.
This supplies **operator-only input**, not energy/model admission. The prepared
composition's `prepared_only_not_executable` status is not changed, and no
qualified periodic response-engine capability is fabricated.

## Frozen prospective numerical settings

Use the existing literal 49.324 Å orthogonal cube, 3,980 waters, seven Na, ten Cl
and one metal (11,958 positive-alpha sites). The source checks these counts and
the selected La/Dy type. Use the previously declared 9 Å real Ewald cutoff,
permanent/polarization alpha 0.4 Å⁻¹, 64³ FFTW grid, both orders 5, vdW cutoff
9 Å/taper 8.1 Å and no vdW correction. Conducting TINFOIL is native default;
`EWALD-BOUNDARY` is forbidden because this native parser selects VACUUM on any
such record. Retain the actual composed polar-eps 1e-8 and explicitly declare
100 iterations/MUTUAL, although **no induced solve is called**. Do not substitute
the historical unexecuted 1e-7 proposal. One thread; no coordinates are wrapped
or sampled. Runtime limits and selected fixture are unbound until root release.

## Why this yields the same reciprocal operator weights

The pinned `udirect1` expression for qfac depends on the initialized cell and
reciprocal lattice, periodic/nonprism flags, nfft, aewald, spline moduli and the
native exp cutoff. Those moduli depend on grid and spline order. It does not
read response amplitudes, beta, h, R/epsilon, POLPAIR, Thole, atom mass or induced
dipoles. Permanent multipoles affect the temporary qgrid/returned field, not
qfac. The B-spline values/indices depend on the literal positions, cell, grid and
order. Chunking is retained with the same one-thread native implementation.

The normalized scatter B adds sqrt(native polarity) and the active douind site
selection. These come from unchanged intrinsic POLARIZE records and topology;
the twelve fitted parameters do not modify intrinsic alpha. H changes the
inverse response operator separately, and must not be substituted for native
alpha in B. POLPAIR/Thole affect real-space interactions, not these PME weights.
Thus deleting only h/H response records preserves the exported reciprocal
operator for the **same future composed geometry and native alpha assignment**.
This argument does not justify using another cell, geometry, grid, alpha,
active-site subset or future parameterization that changes intrinsic alpha.
La and Dy differ at least in metal alpha; never reuse a normalized bound across
them without an explicit new exact alpha transformation or separate export.

## Numerical sign, symmetry and integration contract

Reuse admitted SNAPSHOT.py to parse all actual stored weights and replay norms.
Return qmin, qmax, negative count, zero mode, max|q(k)−q(−k)|, c1, c∞, cF,
`U_R=M*max|q|*min(c1*c∞,cF)` and
`delta_odd=M*max|q−q_conjugate|/2*min(c1*c∞,cF)`.
All are ordinary-double numerical quantities. The first finite water fixture
qualified the actual scatter/FFT/gather chain and observed small nonzero packed
leakage; it is not a measured full-box application-error bound.

For nonnegative stored q, the conceptual even multiplier
`q_even=(q(k)+q(−k))/2` defines a PSD, channel-decoupled reciprocal operator in
exact arithmetic. The original stored q is retained, with its odd part bounded
by delta_odd; no actual multiplier is clamped or symmetrized. Negative q or
nonfinite data blocks this sufficient PSD argument. Export/norm evidence alone
never declares PSD or full stability. The harness returns no dual solution.

A future sufficient stability test can charge delta_odd and a separately
qualified ordinary-double application allowance e_app against its lower margin.
For example, a positive factor of
`B_real+self+H − [1e-6*(U_B+U_R+e_app)+delta_odd+e_app] I`
is sufficient under the already declared upper/lower and arithmetic assumptions.
This preserves the 1e-6 criterion. It is a contract for the owning stability
implementation, not an executed factorization, an error estimate inferred from
PCG convergence, or a new mandatory interval-proof rule. No e_app value is
invented here. Inconclusive bounds stay inconclusive; physical accuracy,
trajectory use and full Hamiltonian forces remain unqualified.
