# Sparse real-operator assembly prototype

Root declares source development and tiny artificial controls only. No saved
water box, factorization backend, native program, new parameter value, or PME
calculation is executed by this task. The running component fit is unchanged.

Assemble the lower triangle of B=I-S(Treal+Tself)S+E beta E^T from the already
reviewed real-pair tensor implementation, streamed in bounded pair batches.
The prototype accepts positive-alpha sites, orthorhombic cells, a cutoff strictly
below half the shortest side, all mutual scales equal to one, per-type Thole
values and explicit symmetric pair overrides. It rejects unsupported domains;
it must not silently approximate native group exceptions or zero-alpha maps.
The fixture-specific wrapper is intentionally absent. A later actual execution
must bind the exact admitted fixture/type inventory and native-operator parity.

Use native real and self signs, coefficient normalization, and minimum images.
Store only one symmetric triangle with no duplicate pairs. Include a fixed
nonnegative diagonal shift only when explicitly supplied. Report the unshifted
site-block upper bound and the applied shift separately. A successful assembly
does not certify the complete PME operator or the original relative spectral
criterion. A standalone B factorization without the qualified reciprocal upper
bound is not that criterion either.

Before any actual large assembly, check the declared input-storage ceiling using
the exact returned pair count and preallocate numeric arrays rather than Python
objects per matrix entry. The estimate covers sparse assembly arrays and bounded
tensor batches, not CHOLMOD fill or all process memory. The backend must perform
its own symbolic/factor memory checks before numeric factorization.

Finite controls: compare a three-site assembly to a directly constructed dense
matrix using the pinned tensor, verify off-axis tensor orientation, rotation and
periodic reimage invariance, explicit pair damping, metal normalized stiffness,
diagonal shift, upper spectral bound, and rejection of unsupported/invalid inputs
and insufficient assembly storage. Ordinary double absolute/relative allowances
are fixed at 1e-13 for these tiny controls. One numerical thread; no physical
accuracy or global-stability claim is inferred from these arithmetic checks.
