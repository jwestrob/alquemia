# Electronic limits of the seven completed PedE DFT endpoints

Read-only inspection of actual completed outputs,13October2026. No calculations,
new tests, acceptance changes or submissions. The unfinished sample1Ca vacuum
attempt is excluded from this inspection and contributes no energy here.

**The independent reference remains useful but is conditional.** Final orbital
tables show near-degeneracy in several water states, apparent non-Aufbau ordering
in sample3La vacuum, and a positive occupied HOMO in sample3Ca vacuum. These
features occur in completed calculations and are not merely initial guesses.
They warrant care when interpreting a solvent double difference as a test of
GFN2. They do not uniquely prove electron detachment, a preparation defect or
the physically correct electronic state.

## Exact printed final orbitals

All energies below are printed Eh, not newly computed gaps. Occupations are2.0000
and0.0000. La frontier indices are413/414; Ca indices418/419. Paths are beneath
`workspaces/pede_solvent_response_20261013/`.

| Source/metal/medium | HOMO | LUMO | Output path and exact lines |
|---|---:|---:|---|
| sample1La water | −0.164384 | −0.064384 | `reference_v1/tasks/sample1_La_water/endpoint.out`:2393–2394 |
| sample1La vacuum | −0.022786 | +0.226229 | `reference_v1/tasks/sample1_La_vacuum/endpoint.out`:2486–2487 |
| sample3La water | −0.185536 | −0.185447 | `reference_sample3_La_v1/tasks/sample3_La_water/endpoint.out`:2289–2290 |
| sample3La vacuum | −0.018918 | −0.019662 | `reference_sample3_La_v1/tasks/sample3_La_vacuum/endpoint.out`:2284–2285 |
| sample1Ca water | −0.156607 | −0.156548 | `reference_sample1_Ca_v1/tasks/sample1_Ca_water/endpoint.out`:2331–2332 |
| sample3Ca water | −0.176744 | −0.176466 | `reference_sample3_Ca_v1/tasks/sample3_Ca_water/endpoint.out`:2265–2266 |
| sample3Ca vacuum | +0.071346 | +0.171346 | `reference_sample3_Ca_v1/tasks/sample3_Ca_vacuum/endpoint.out`:2325–2326 |

The sample3La vacuum occupied level is printed above an unoccupied one. Its
SOSCF energy/gradient stopping condition therefore should not be treated as
independent evidence of a stable ground-state occupation. The positive
sample3Ca vacuum occupied level is a gas-phase electronic-reference warning;
finite-basis approximate orbital energies alone do not establish electron loss
or a unique explanation for the solvent contrast.

The large positive table separations in three TRAH cases need a separate caveat:
sample1La water has a near-final small-gap warning0.000117 at1913 and another
post-convergence warning at1929; sample1La vacuum has a near-final negative-gap
warning−0.000983 at2008 and another at2029; sample3Ca vacuum has−0.004565 at1847
and another at1863. These precede or follow the final converged macro. They are
not just startup warnings. The table/warning discrepancy is retained; no specific
level-shift or final-diagonalization explanation was established by this audit.
Do not use the printed positive table separations as proof of wide-gap stability.

## Solver residuals and state accounting

The actual final TRAH/NR residuals are below the printed1e-5 threshold:
sample1La water5.663e-7(line1922), sample1La vacuum8.692805e-7(line2022),
sample1Ca water4.453178e-7(line1854), sample3Ca vacuum8.312700e-7(line1856).
The other three end in SOSCF, where the active stopping evidence is:

| Endpoint | Native stop | Final printed orbital gradient | Lines |
|---|---|---:|---|
| sample3La vacuum | Energy check | 1.9240e-6 | 1819–1824,1862 |
| sample3La water | Energy check | 1.1701e-5 | 1817–1822,1867 |
| sample3Ca water | Gradient check | 2.3077e-6 | 1788–1793,1838 |

Sample3La water's printed gradient is slightly above its printed1e-5 tolerance,
although the native energy stopping condition passes. Several printed density
and orbital-rotation residuals also exceed their displayed tolerances. This is
not an assertion that ORCA required every criterion simultaneously; nor are
stale DIIS summaries substituted for actual TRAH residuals. The existing
diagnostic admissions are unchanged. Fine response accuracy, repeatability and
electronic-state uniqueness remain unqualified.

No total-state mismatch was found. Ca outputs print restricted singlet,charge−3,
838explicit electrons at1602–1605; La outputs print restricted singlet,charge−2,
828explicit electrons at1624–1627 and La ECP46 at1457. The differing explicit
electron counts reflect the declared DFT ECP treatment, not inconsistent paired
charge bookkeeping. Integrated alpha/beta densities match each other and stay
close to the specified counts. Printed atomic charges sum to−3/−2. The neutral
La atomic-guess doublet is not the molecular multiplicity.

No alternate spin/occupation search or explicit stability result is available.
ORCA distinguishes convergence from stability and global-solution selection;
TRAH's local-minimum property does not validate a chosen physical spin/state or
its biological transfer.
[Official SCF documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html).

## Consequence for the comparison

Use actual final endpoint totals for a conditional same-composition solvent
response comparison. Preserve these electronic warnings beside the result.
Do not equate Mulliken/Loewdin populations across DFT and xTB, infer a measured
oxidation state from their fractional values, or promote normal termination to
proof that this DFT treatment is an oracle. The present evidence justifies
finishing the declared finite comparison; it does not authorize a broad spin,
functional or basis sweep, or a retuned classifier threshold.
