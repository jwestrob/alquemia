# Electronic limits of all eight completed PedE DFT endpoints

Read-only inspection of actual completed outputs,13October2026. No calculations,
new tests, acceptance changes or submissions. Updated after primary job3372
completed sample1Ca vacuum and saved-output comparison3377 completed all eight
declared cells. The conditional seeded-restart preparation was not needed.

**The independent reference remains useful but is conditional.** Final orbital
tables show near-degeneracy in several water states, apparent non-Aufbau ordering
in sample3La vacuum, and positive occupied HOMOs in both Ca vacuum endpoints. These
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
| sample1Ca vacuum | +0.060446 | +0.310295 | `reference_sample1_Ca_v1/tasks/sample1_Ca_vacuum/endpoint.out`:3422–3423 |
| sample3Ca water | −0.176744 | −0.176466 | `reference_sample3_Ca_v1/tasks/sample3_Ca_water/endpoint.out`:2265–2266 |
| sample3Ca vacuum | +0.071346 | +0.171346 | `reference_sample3_Ca_v1/tasks/sample3_Ca_vacuum/endpoint.out`:2325–2326 |

The sample3La vacuum occupied level is printed above an unoccupied one. Its
SOSCF energy/gradient stopping condition therefore should not be treated as
independent evidence of a stable ground-state occupation. The positive
Ca vacuum occupied levels are gas-phase electronic-reference warnings;
finite-basis approximate orbital energies alone do not establish electron loss
or a unique explanation for the solvent contrast.

The large positive table separations in four TRAH cases need a separate caveat:
sample1La water has a near-final small-gap warning0.000117 at1913 and another
post-convergence warning at1929; sample1La vacuum has a near-final negative-gap
warning−0.000983 at2008 and another at2029; sample3Ca vacuum has−0.004565 at1847
and another at1863; sample1Ca vacuum has−0.000142 at2938 and a small negative
gap warning after convergence at2960. These precede or follow the final
converged macro. They are
not just startup warnings. The table/warning discrepancy is retained; no specific
level-shift or final-diagonalization explanation was established by this audit.
Do not use the printed positive table separations as proof of wide-gap stability.

## Solver residuals and state accounting

The actual final TRAH/NR residuals are below the printed1e-5 threshold:
sample1La water5.663e-7(line1922), sample1La vacuum8.692805e-7(line2022),
sample1Ca water4.453178e-7(line1854), sample3Ca vacuum8.312700e-7(line1856),
and sample1Ca vacuum6.267086e-6(line2953).
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

### What the final long attempt actually did

The sample1Ca vacuum output prints default `MaxIter125` at1656, executes50
ordinary DIIS iterations through1777, switches automatically to TRAH, reaches
macro98 at2953, and reports `SCF CONVERGED AFTER843CYCLES` at2957. The final
energy is −6718.448813177293Eh at7187; the final accepted NR residual is quoted
above. The large density and DIIS values at2994–2996 are not substituted for
that final NR residual.

The executed `endpoint.runtime.inp` requests r2SCAN-3c, NoAutostart, DefGrid3,
TightSCF,12 ranks and MaxCore3100; it contains no explicit MaxIter change.
Consequently, the earlier description of125 as a total iteration/work cap was
incorrect: native ordinary iterations, TRAH macro/micro iterations and its
reported SCF cycles are different counters. No explicit iteration, rank or
time-limit override is evidenced by these files. The exact internal limit
semantics are not inferred from the total843 count. This unchanged primary
attempt converged; it was not replaced by the prepared seeded contingency.

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
proof that this DFT treatment is an oracle.

The completed comparison reports the sample1-minus-sample3 solvent double
difference as +10.03421633832605kcal/mol for r2SCAN-3c/CPCM and
+21.256628301423902kcal/mol for corrected GFN2/ALPB. The direction survives
the independent treatment, while its magnitude is model-dependent. This is
neither zero solvent sensitivity nor proof that DFT supplies an exact target.
Both the electronic treatment and continuum model changed together.

The next useful distinction is frozen-distribution solvent response versus
additional self-consistent redistribution. Replay the existing ALPB energy
functional with each source's saved vacuum and solvent charges, first requiring
reproduction of its actual converged component totals. Preserve the geometry,
parameters and every component, including charge-dependent hydrogen bonding.
Do not replace the production correction with that frozen-charge diagnostic or
discard its intrinsic compensation. A maintained implementation/API feasibility
check is owned separately; this note does not claim that replay has run.

The historical [native-GFN2/CPCM pilot](../../solvent_cpcm_20260922/REPORT.md)
produced no complete metal pair. The [CHELPG sampling](../../mace_omol_20260917/CHELPG_SAMPLING_STABILITY_REPORT.md)
and [ddX refinement](../../mace_omol_20260917/DDX_CPCM_REFINEMENT_REPORT.md)
studies retained sampling/refinement failures. Repeating those campaigns is
not the next recommendation. A later four-cell GBSA test would be a whole-
solvent-model sensitivity test, because its fitted Born, surface and hydrogen-
bond parameters differ from ALPB; it would not by itself isolate density
response. No broad spin, functional or basis sweep, new fit, or relaxed
tolerance is proposed.
