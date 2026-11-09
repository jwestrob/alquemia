# One bounded preparation comparison

Declared9 November2026 after audit66010, before this comparison. The first
observed damage is already in relaxation frame20: trans600 carbon7 separates
from atoms6/8/9. Frame10 was still intact. Input atom order, periodic images
and checkpoint/HDF correspondence pass. This does not by itself distinguish
an optimizer overshoot from an unfavorable potential surface.

Use exactly the original hydrated start0-A input and unchanged L potential,
head5, charge0/multiplicity6, float64/highest precision, PME and cutoffs. First
compare one ASE force evaluation at the exact final native checkpoint with
the saved native force array. JAX-MD space.transform's custom JVP returns
Cartesian forces even for fractional positions; do not multiply forces by42.
Require agreement within1e-5 eV/A, the established numerical tolerance. Save
full raw arrays. If that fails, stop this comparison as an interface mismatch.

If it passes, make one conservative optimization comparison: standard ASE
LBFGS, maximum displacement0.05 A per atom, default alpha70/memory100/damping1,
no line search, at most200 accepted steps. Stop at the original maximum-force
component criterion0.2 eV/A. All atoms remain mobile. No bond, donor, fragment
charge or association restraints are added; the potential is unchanged.
Record every accepted geometry and stop if an original ligand bond leaves
the already declared0.7–2.0 A range. Retain proton-transfer diagnostics.

This tests the numerical preparation path. It does not retune selectivity,
extend the failed FIRE trajectory, change the model, or establish chemical
accuracy. Preserve the failed FIRE/MD result. There is no optimizer setting
search or automatic extension if this single comparison fails its checks.

Allocate one H200 and12 CPUs under standard_expand, default scheduler memory
and time. Bound the static check and optimization together by an internal
600-second timeout. Count its entire allocation within the existing GPU-day.
Do not start additional MD before inspecting the actual outcome.

9 November outcome and separate diagnostic scope: job66047 stopped at the
force-bridge gate, before any optimization. Maximum force difference was
1.9643643e-5 eV/A, above the frozen1e-5 tolerance. Preserve that failure; do not
raise its tolerance or call the native interface qualified. This difference
is0.000453 kcal/mol/A and does not show a gross coordinate-unit factor.

Run the previously specified capped optimization once as an independent
ASE-only diagnostic (`--ase-only-diagnostic`), reusing the failed bridge record
without repeating its calculation. Its energy and gradient come from the
same pinned ASE potential, so the numerical preparation question can be
examined without asserting native/ASE equivalence. Keep all the declared
LBFGS settings, original input,200-step limit and chemical checks unchanged.
This is the first execution of that capped optimizer, not a setting search.
Use a fresh directory and the same600-second GPU bound. No further MD is
released by this diagnostic, whatever its outcome.
