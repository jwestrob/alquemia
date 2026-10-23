# Model-side susceptibility convention

The quantum target remains the fixed h=0.001 atomic-unit energy curvature for
each of the four states. No target, weight or physical acceptance gate changes.

At each fixed geometry the declared native family remains linear in the induced
dipoles: its geometry-dependent h and H do not depend on the external field.
With K the native electrostatic constant, B the Bohr/Angstrom conversion, and
Q the Hartree/kcal conversion, the projected model susceptibility is

    chi_zz = Q/(K B^2) z^T (A+H)^(-1) z.

This is the exact centered energy curvature at h=0.001 for the quadratic model
in exact arithmetic. The evaluator uses the qualified zero-field exported
operator, with conditioned solve checks, instead of making additional field
calls. It never substitutes a moment estimator for the quantum energy target.

The focused saved controls reconstruct all four original model curvatures from
12 pinned zero/plus/minus native energy outputs. The analytic-minus-printed
curvature differences are about -8.756e-6, -8.002e-6, -4.993e-6 and -8.433e-6 a0^3
(La dimer/trimer, Dy dimer/trimer). Each lies within the original 3.18720e-5 a0^3
energy printing bound. The raw tokens, exact propagated printing bounds and
both baseline objective values are saved in SUSCEPTIBILITY_CONVENTION.json.
Initial J changes from 2.2181741096669696 to 2.218165375513068, solely from this
model-side numerical representation apart from negligible cycle summation
roundoff. Original E/F, pair and C blocks remain the same to their arithmetic
precision. Root explicitly accepted this documented convention before fitting.

Printing bounds are not SCF/physical error bars. The prior actual 54-fixture
operator review also checks the four field aliases against saved site moments.
New native qualification must separately establish the new functional's
energy, derivative, moment and export consistency before any fit call.
