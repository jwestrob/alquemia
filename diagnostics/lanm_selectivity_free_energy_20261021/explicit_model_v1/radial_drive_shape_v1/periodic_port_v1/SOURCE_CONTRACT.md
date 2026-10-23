# Periodic radial port — source handoff

## Implemented and bounded

Root declaration47750f75 governs this source-only package. The copied source inventory is the immutable PME-export parent source_v3, including its prior compiler-import recovery. All615 inventoried source files are copied to the fresh workspace source_v2. Only source/polar.f and source/induce.f differ; the other613 entries retain their parent hashes. The two literal diffs and EDITS.json are retained. Reversing the four induction joins, one module addition and appended helper restores both original files byte for byte.

`lrpshape(2)` stores independent BIO e values. The dispatcher notices POLAR-RADIAL-DRIVE even without a local-response record and resets new state on every setup. The dedicated parser accepts only900/901–402, finite e∈[−1,1], unique metal records and no trailing token. A malformed or mixed new record fails. Parsing e precedes the all-zero return, so an e-only active state receives the complete existing topology, periodic-mode and stability guards. An absent/zero e performs no added field or gradient arithmetic. Supported original vacuum500/501–36 setup/geometry/gradient routines are byte-identical. This periodic library does not implement the separate vacuum POLAR-RADIAL-DRIVE500/501 keyword; that keyword is rejected here instead of being silently interpreted as BIO.

After the old periodic h loop, the new field helper adds delta h=3e Ftilde phi psi R, psi=(r²−6.25)/6.25, to each selected primary oxygen. It uses the same nearest-image geometry and compact kernel already computed by the parent. Outside21 Å it returns before forming psi. The helper uses nw-sized arguments and adds no dense tensor or other-water loop. e does not change H, beta, intrinsic alpha, damping or the original environment sums.

The existing second gradient pass passes w=u+v to the new helper. Before the parent's factor K/2, it adds

    −3e [phi psi (w·R) grad(Ftilde)
          + Ftilde {phi psi w + (psi phi′/r + 2phi/6.25)(w·R)R}].

The parent's opposite metal reaction and symmetric image virial then apply exactly once. Ftilde′ already includes the20–21 Å switch derivative; no separate switch is evaluated. At contact the phi′/r contribution is set to zero and the finite Jacobian is −3e F(0)kappa I. Exact r0 leaves delta h zero but the radial derivative nonzero. The common dfield0/ufield0 code, stationary energy hooks, PME kernels, native frame torques and both fatal periodic-guard bodies are unchanged. Total h is therefore added to both direct RHS through the original common hook; no separate scalar energy term is introduced.

## Actual controls and preserved preparation failure

The first preparer assumed every subroutine had no arguments and failed while locating periodicresponsesetup(handled), before copying source files. Its source, INPUTS.json, empty source_v1 and FAILURE.json bcc0d7a4 remain. PREPARE_SOURCE_v2.py only corrects that locator to accept routine arguments and selects the versioned inputs/output. No radial Fortran formula changed.

The successful v2 preparation took4.464096 wall/2.231557 CPU seconds on one thread. All62 focused source/manufactured controls passed in .692372 wall/.688190 CPU seconds. They establish inverse/source joins; unchanged vacuum/common fields/guards and other source hashes; explicit e-only/zero/parser guard placement; independent B/B′ versus native-factored unequal-dual contraction; Cartesian finite differences at contact, near-contact, r0, interior and switch radii; r0 slope; exact support/odd-e/zero properties; switch limiting orders; reaction and affine virial finite differences; and orthogonal-image equivalence/ties outside support.

These controls evaluate a small manufactured explicit-field contraction with fixed unequal u/v. They do not compile or execute Fortran, exercise the native parser, solve an induced system or establish full self-consistent energy/force consistency. Source parser checks are static. No physical geometry, candidate, FFT, factorization or solver was evaluated. No broad native qualification was repeated.

## Concrete future bindings; execution still disabled

Both normal periodic guards still call fatal. There is no new key/environment bypass. A later qualification source must retain the exact-state transaction and add e to the versioned state serializer, parser, raw reporting and identity comparison. Native h/H alone does not identify e at r0. The future composition adapter must bind an admitted final14 candidate, retain the original12 mapping/BIO assets and render e explicitly. Nothing here selects or composes a candidate.

Future active qualification must bind12 Å real Ewald cutoff consistently in key, gate, independent assembler and operator certificate; the old e41f9 Å state/factor is not a substitute. Other settings remain explicit. e alone changes only the RHS at fixed old coefficients, but that observation does not authorize cross-state certificate reuse. The source tree is unbuilt; a build must supply the new polar module/data object and matching induce object without replacing the parent archive. This package intentionally does not construct the69-program runtime.

Reuse the ordinary-double nominal stability contract with e_app unknown/optional. It does not certify an unsampled trajectory or full force consistency. The inactive anchored12 baseline pass remains separate; future active strain checks must use the corrected anchored derivative and prospectively bound steps. The previous compact-tail bounds omit psi and are not radial14 approximation bounds. No physical accuracy, global stability, native qualification, MD performance or production admission follows from this source handoff.
