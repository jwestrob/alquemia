# Independent cavity-scout source review

**Preparation, Fortran component caller and final native driver have no identified
launch blocker.**
This is static source inspection, not numerical qualification. The Python native
driver has also been inspected; its final review is recorded below. No scientific imports,
calculations, compilation, tests or submissions were performed by this reviewer.

## Reviewed snapshots

- Preparation manifest: `workspaces/native_cavity_shielding_20261014/preparation_v1/manifest.json`,
  SHA256 `e4e0a3d878a6a73ae7df7af75d79afc318c212e58ce4f9477bf38fc642c6b247`.
- `prepare/PREPARE.py`: `20ea8ee297242af4fc08f8793bb5a162fbb2f60ff6d62270de4ae3a7e58a05b0`.
- `native/CAVITY.f90`: `0a61d74b43bfb2bf9f8a3c0afb1e12416900a847ad55796ec6f9db832c5e5380`.
- Final `native/RUN.py`: `c2b8201a655dea878ebf99dbc93abbfc52a5c90e55a970a4401c0b6b605b3c81`.
- `native/RUN.sbatch`: `40a241f1a7c4d61bca24ad2e75716324d32eb2411148b870191bf46e2a21fa4f`.

The preparation matches source atom identities rather than proximity, replaces
the active source Ca exactly once, retains the spectator Ca, preserves the three
modeled physical PQQ hydrogens and distinguishes all 14 synthetic caps. The local
atoms retain their original coordinates and archived charges. Added physical
source atoms have exact zero diagnostic charges. Removing cap centers leaves all
real-atom charges unchanged. Translation, proper rotation and reversal carry the
charges, cap flags and uncertainty flags with the atoms. A preparation failure
remains unavailable; paired metal coordinates are checked before admission.

The caller agrees with the pinned native `gbsa.f90` and `kernel.f90`: `Gelec`
excludes the `2*hbw` matrix diagonal, the shape contribution contains the factor
one half, and the two physical/cap matrix orientations give one complete cross
term. The native P16 routine populates both matrix triangles. Charge-precision
bounds include only the archived population sites; appended zeros and cap sites
zeroed for the control have no artificial rounding uncertainty. The cap-removed
call uses the same real charge vector as the cap-zero union call. The wrapper
does not copy or print a second dense matrix, though native dense allocations
remain and require measured allocation-aware execution.

## Interpretation and runtime gates

This evaluates a capped-union exclusion sensitivity. It neither supplies a
physical source-only cavity nor gives the added protein atoms physical zero
charges. Quadratic cap attribution and the cap-volume control are diagnostics,
not electronic states, error bounds or score corrections. Whole-source SASA,
hydrogen bonding and reference shifts must remain separately reported.

Before the remaining six cavity updates, the sample0-Ca scout must pass the
declared local replay, repeated-charge, electrostatic invariance and quadratic
closure checks at `1e-9 Eh`, with an informative charge bound at `1e-6 Eh`.
The complete declared denominator is 12 cavity updates; a failed scout must
leave continuation cells visibly unexecuted. Numerical success cannot be inferred
from this source review.

## Native-driver review

The driver preserves the six-update scout followed by six continuation updates,
with a successful scout required before expansion. It verifies the actual prepared
coordinates, charge vectors, masks and transformations; rechecks archived native
and component-v2 acceptance; preserves the four-source denominator on failure;
and separates shape and cap-volume responses from the primary electrostatic
transfer. No classifier bands or full-source nonpolar correction enter the result.

The actual build4543 receipt reports a successful wrapper-only build. Its build
settings, native-helper pin and linked-library pins match the qualified manifest.
A copied-source path comparison and an assumption that the build manifest contained
a `native_model` field were repaired before freezing the numerical runner. These
were metadata problems, not scientific failures or reasons to rebuild the engine.

The final driver applies the existing `1e-6 Eh` ceiling to the reported summed
uncertainty of union-minus-local and cap-volume differences before continuation
and endpoint admission. Paired and source-double-difference uncertainty flags are
also explicit; a complete endpoint matrix alone does not assert that every derived
contrast is informative. Dispatch accounting distinguishes START intent records
from actual child execution receipts. Both requested changes are source-verified;
tolerances and energies are unchanged. No runtime numerical qualification is
claimed by this review.
