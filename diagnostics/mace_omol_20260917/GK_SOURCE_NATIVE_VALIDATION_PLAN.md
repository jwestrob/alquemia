# Validate the consequential source-only GK component with the native solver

The analytical decomposition of the pinned GK expression attributes +18.1290
kcal/mol of the 2FW0-minus-2FVY GK contrast to the projected source interacting
with its own reaction field, offset by -5.5288 in the source/environment cross
term. This motivates checking the analytical decomposition before choosing an
improved solvation model. No predictor is changed or accepted by this test.

Use the original whole physical cavities and actual primary native inputs for
GGR1GLGextended, GGR2FW0, GGR2FVY, alpha1F6S and alpha6IP9. Five geometries,
three static evaluations each: Ca source-only, La source-only, and entirely
uncharged cavity. Exactly15native energy calls; no fields, induced solve,
quantum, charge fitting, MACE, gradients or geometry changes. Native thread
allocation64CPU, existing scheduler policy, no project time/compute cap.

Build an isolated frontend against the same pinned native library. Add only
two explicit static modes to the qualified frontend: `source_self` sets all
environment permanent moments to zero while retaining the actual source
monopoles, and `empty` zeros all permanent moments. Both disable induced
response while retaining the same native Born/cavity settings and all physical
atoms. No base parameter, cavity radius, descreening, neck or tanh change.
The existing energy routine still computes nonpolar solvation. Subtract the
empty-cavity energy from each source-only energy to remove that common term.

Verify actual source charge/zero higher moments, zero exterior moments, paired
coordinates, retained physical atom inventory, base parameters, zero induced
arrays and no extra active terms. Retain native effective radii and require
agreement with the archived full-system radii within1e-10 A. Compare each
native source-only-minus-empty GK energy and its Ca-minus-La difference to
the independently implemented analytical expression within1e-7 kcal/mol.
These tolerances address deterministic algebra/backend agreement, not affinity.

Preserve all native logs/receipts, failed attempts, measured cost and old
results. A pass validates only the component accounting. A failure must be
diagnosed without modifying old results or reporting an environmental scalar.
