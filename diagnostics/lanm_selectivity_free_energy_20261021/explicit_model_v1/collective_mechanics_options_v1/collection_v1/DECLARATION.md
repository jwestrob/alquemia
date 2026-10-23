# Saved three-body collection

Collect exactly the eight cycles declared in ../TRIMER_DECLARATION.md. Require
all 16 new reference rows, all 16 new model rows, and all explicitly named reused
terms. This source performs no native molecular calls and no fitting.

Use electronic reference energies without D4, 627.509474 kcal/mol per hartree.
Report model total and individual component E3 separately. Retain all raw energy
and convergence metadata. Seven-term subtraction is performed with accurate
summation; no energy sign or magnitude determines admission.

For electric moments, independently read the ORCA raw output and full-precision
property file. Reconstruct the common-origin dipole from the effective nuclear
charges and full-precision Bohr coordinates, using the established nuclear
moment reconstruction of the native origin. Check it against the printed origin,
rounded-mass center of mass, and printed dipoles using the previously admitted
rounding bounds. No SCF-accuracy or physical uncertainty is inferred from these
bounds. Reused full-precision moment rows must match their admitted authority.

For the model, replay the unchanged native induced/total moments and energy
reader. Permanent dipoles cancel because complete atom identities and local
frames have zero net cycle coefficients. Check the precise induced sum against
the lower-precision total laboratory-frame moment cycle. Convert quantum moments
with the pinned native convention 0.529177210544 * 4.803206802 Debye/au.

Report every vector and model-minus-reference error, along with energy component
cycles and coefficient-propagated printing bounds. Retain near/far results even
if symmetry suppresses a dipole. Do not sum these eight sampled E3 values into B,
interpret all quantum E3 as polarization, promote a reference method, or choose a
repair coefficient from this diagnostic.
