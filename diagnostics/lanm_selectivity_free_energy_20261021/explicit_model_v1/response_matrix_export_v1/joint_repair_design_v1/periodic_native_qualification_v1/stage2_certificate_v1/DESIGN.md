# Stage2 exact-state gate and one-gradient harness

Source preparation only. The question is whether the unchanged compact periodic
response Hamiltonian has conservative native Cartesian/cell derivatives at the
old e41f La center and the 32 already prepared displacements. No build, native
initialization, assembly, factor, solve, gradient or outcome read is performed
by this package. The Stage1 actual admission and root Stage2 release are null.

## Finite inventory and source isolation

Reuse numerical contract c5b85118, 32-deck manifest dbd21ded, fixed directions
48ed837a, and root geometry review befdcb02. There are 33 active E/G programs
(one center, 32 displacements), and exactly three central bypass E/G programs:
unchanged old library/absent response, qualification clone/absent response,
qualification clone/zero response. Each calls native gradient once, returns
total energy, all 35874 positive Cartesian energy derivatives and full 3x3
native virial. No internal finite differences or optimizer is present.

The center may reuse job20053 only after exact Stage1 initialized-state admission.
Every displacement needs its own one-direct-initialization program, own PME
bound, own unshifted assembly and sufficient-factor admission before its solve.
Thus the contract retains 32 initialization + 36 E/G = 68 Stage2 native programs,
32 assemblies and 32 factors. Unadmitted slots remain untouched; no central
certificate is propagated to another state. Source files below do not launch
any of these tasks.

The generator copies the pinned parent induce.f into a separate qualification
object. It changes only the two guard call sites and adds observation/counter
calls. The two original fatal guard bodies remain bytewise unchanged. All PCG,
field, kernel, energy, torque and derivative arithmetic is unchanged. Explicit
qualification induce object + reporter/gate objects link before the immutable
archive61e2fab0; the original archive/tree is never replaced. The old-library
bypass reporter links without the qualification induce object.

## Exact certificate transaction

1. The future wrapper verifies Stage1 admission, exact row/XYZ/key/prm, candidate,
   source/build/executable, certificate producer/result and root admission pins.
   An active certificate must bind its own state file, PME snapshot/bound,
   unshifted matrix and nominal factor result, including unchanged delta_odd,
   e_app:null and the original nominal threshold. The wrapper stages an immutable
   expected_state.dat and a small input ticket in a fresh directory. The ticket
   itself is not a scientific certificate and cannot authorize a different state.
2. The early backend hook checks the explicit qualification role and supported
   fixed runtime before any induced solve. It cannot validate PME before it exists.
3. Common dfield0 saves native raw Ed/Ep immediately before the shared h addition
   and both RHS immediately afterward. Its arithmetic is unchanged. The difference
   is checked against h with 128*N*eps*max(1,raw,RHS,h), separately per channel.
4. The existing stability hook runs immediately after that same dfield0 call,
   before direct guesses, preconditioner or the first PCG mutual application.
   It writes the actual state with the **same Stage1 state serializer** and
   compares every complete line against the wrapper-pinned expected state.
   Exact mismatch, missing file or trailing data is fatal. Only this comparison
   admits the single subsequent solve; a startup ticket alone never does.
5. The initialization-only executable has no induce/gradient call. It uses the
   same checked setup, dfield0 and serializer, emits its state and terminates.
   That output becomes an expected state only after the separate external
   assembly/PME/factor/root reviews. INITIALIZE role cannot enter an induced solve.

The serialized state contains exact atom/type/order, coordinates, intrinsic
alpha, damping, thole and metal-O pair gamma, both response coefficient records,
native h/H, raw unequal Ed/Ep, PME grid/order/exponent/volume/transform, igrid,
the first two spline derivative rows and every unsymmetrized qfac value. These
are sufficient for the admitted induced operator, with topology, all other
parameters, all-unity mutual scales and mode pinned and checked separately.
The decimal tokens use the original 18-significant-digit rendering, sufficient
to round-trip binary64; comparison is stricter than a tolerance comparison.
The cell enters through the exact input pin, geometry/volume/reciprocal transform,
and explicit native mode checks. No hash implementation is added inside Fortran.

Central expected_state.dat is the future **admitted Stage1 parity_state.dat**,
not a hand-created snapshot. For displacements it is the same serializer's
initialization output. No new center initialization is counted or needed. The
saved qfac/splines and native h/H must match; equal input file hashes alone are
insufficient. Platform-dependent token differences fail and are preserved.

## Reporting and numerical accounting

The reporter runs one full gradient, then writes total E/G/virial and the native
polarization ep separately. It does not read a last finite-difference state.
For active rows, direct fields, corrected RHS, intrinsic alpha, final u/v and
h/H are retained. Central replay is -electric/(2*dielec)*sum(u*(Ep+h)); native ep
is compared under 1e-8 kcal/mol plus the inherited arithmetic allowance, not
against total energy. The three bypass E/G/virial numeric tokens must be exact.

The clone counts common direct and mutual field calls. It captures the solver's
existing pre-peek iteration/residual values, labels them as pre-peek, and records
the final post-peek u/v. There are **zero extra residual/operator readbacks**:
the number of ordinary PCG mutual applications is an actual runtime quantity,
not hidden as zero. No native fallback is enabled (ulist is already forbidden).
Native poleps=1e-8 and politer=100 remain unchanged; failures do not trigger a retry.

Future collection uses the frozen eight direction arrays and literal E/G/virial
tokens, both step sizes and all 24 original slope/step checks. Printing allowance
and 0.001 thresholds are exactly the numerical contract. It preserves all
missing/failed rows and cannot promote a partial panel. No new physical threshold,
application-error proof, changed candidate or production/trajectory claim follows.

## Remaining implementation/binding work

The supplied native sources and isolated source generator are unbuilt. The
certificate wrapper and final finite-difference collector are **not implemented
in this package**. They must verify the transaction above, including row-specific
certificate status/pins and Stage1 actual admission, before execution is possible.
The 32 assembly/PME/factor adapter bindings and resource release remain null.
The future build recipe must be source-reviewed against the one-link Stage1
compiler/module/archive recipe; no automatic compilation is supplied. Read-only
source checks here establish edits/interfaces only, not Fortran compilation or
native validity. These are concrete remaining engineering steps, not scientific
evidence of stability at changed states.

The three literal bypass decks are prepared in bypass_preparation_v1: both
absent inputs delete exactly the two response records; the zero input changes
only their eight amplitudes to zero. Geometry, topology, all other key records
and BIO parameters are identical. The unchanged old-library bypass has no new
internal counter hooks, which is reported as unavailable rather than fabricated
zero; its one gradient call is explicit and its single induction path is source
bound. The two clone bypasses retain measured internal counters and require no
active certificate. Their exact output-token equality remains untested.

Initial generated source_v1 is preserved as an unbuilt draft. Source inspection
corrected two module imports (PME variables from pme; electric from chgpot) in
source_v2. No compiler/native failure or result was involved. Both generations'
inverse source edits restore the same original induce.f exactly; this check is
not a substitute for a future compile or independent source review.

Root's scheduler design permits at most16 concurrent unchanged12-CPU factors
(192 CPUs, leaving four of the20 job slots), default48000MiB per case. This is
only a prospective resource envelope; initialization/gradient programs remain
one native thread, and no factor allocation or runtime release exists here.
