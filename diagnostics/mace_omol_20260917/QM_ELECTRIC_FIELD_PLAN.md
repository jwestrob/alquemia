# Saved-density electric fields before induction scoring

Declared after Tinker framework controls and before any new field outputs.
Active-goal authorization applies. No biological classification or parameter
selection is made in this experiment.

## Question and real inputs

Does the existing native CHELPG distribution, including its fixed physical cap
projection, reproduce the electric field of the saved quantum endpoint density
at real environmental atoms? Earlier potential/coupling tests cannot establish
this. A poor source field could spoil induced response even with a fast solver.

Use all8normalized vacuum r2SCAN-3c endpoints from1200950 and the accepted
normalized_charge_report_v1. Four consumed preparations: GGR extended/connected,
alpha1F6S/6IP9; bothCa/La. No new SCF, fit, gradient, geometry, microstate, water,
assembly, MACE, force-field energy or response optimization. The QM source and
physical projection remain exactly as already recorded, with no refitting.

Probe every physical protein/water atom outside the recorded QM projection
support, in physical source order. Obtain AMOEBA2018 polarizability weights
from the real framework mapping; omit no atom based on the output field.
All sites, including hydrogens, remain recorded. Ca/La probes match exactly.
ConnectedGGR uses the same physical framework with its own recorded support.

## Eight native utility calls and fixed numerical differentiation

One ORCA6.1.1 `orca_vpot` call per saved endpoint. Pack the six Cartesian
offsets for each of two central-difference spacings,0.001 and0.0005bohr,
around every probe into that call. Use minus the spatial derivative of the
potential for electric field. This is differentiation of the saved-density
potential at observation points, **not numerical DFT nuclear gradients**.
Do not move atoms or rerun SCF. Preserve the actual native potential precision,
point ordering, units and executable/wavefunction hashes.

Compare both derivatives. Require max vector difference≤1e−6atomic field units
and change in each paired diagnostic below≤0.01kcal/mol. Use the smaller step
for the declared comparisons. Point-charge fields are evaluated analytically
at exactly the same sites, first at original QM/cap positions, then with the
existing charge/dipole-conserving physical projection. No damping or covalent
scaling is imposed: this isolates the source representation before a complete
AMOEBA boundary model exists.

## Outputs and interpretation fixed before execution

Record exact, fitted and projected vector fields; endpoint and Ca-minus-La
errors; ordinary and polarizability-weighted RMS, maximum errors, and distance
to the nearest QM/cap atom. Flag a representation as passing the field screen
when weighted RMS error≤1e−4atomic units OR relative weighted RMS≤0.10,
for every endpoint and matched pair. The10% criterion is a representation
screen, not a guarantee of sub-kcal induction accuracy or a classifier gate.

Also report the deliberately narrow diagonal response diagnostic
`U0=-0.5 sum_i(alpha_i |E_i|^2)` in consistent atomic units, converted once
to kcal/mol. It ignores mutual polarization, protein permanent fields,
solvent and native covalent/damping rules. It is **not** an environmental
correction. Flag paired projected-minus-exact U0 errors above1kcal/mol, half
the existing2kcal partition tolerance, as a reason to investigate field quality.
Do not add U0 to any score. Report all-site and ≥3Å-from-every-QM/cap-atom
strata separately; the latter is a fixed geometric diagnostic, not a rescue
subset or an alternative scientific result. Close boundary sites remain in
the all-site denominator even though eventual native scaling can alter them.

Failures distinguish native numerical resolution, charge fit and cap projection.
No post-result threshold change or charge-scheme selection per protein. A pass
cannot settle source damping, covalent exclusions, QM/FF boundary charge closure
or full-model partition/predictive accuracy. A failure motivates a new version
of the representation, retaining the inspected results as development evidence.

Use8one-thread utility workers in the established shared CPU allocation,16GB.
Eightcalls,~12times the points in the earlier8-call potential experiment whose
summed utility time was88.160987seconds. Expected seconds to minutes, no project
time/CPU limit. Actual cost and failures count. Existing pinned utilities,
execution receipts, immutable attempts, restart checks and cache hashes apply.
Preparation is read-only source extraction; every compute product goes under
workspaces/. No production baseline or reference changes.
