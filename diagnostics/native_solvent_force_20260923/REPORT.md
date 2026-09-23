# Continued native solvent energies do not yet support reliable accommodation forces

**Both tested proteins fail the frozen force-consistency gates, despite passing
the energy-settling gates.** Do not use this fixed two-continuation solvent
recipe to drive composite geometry optimization. The result does not change the
released static scorer or establish a new classification failure.

All 96 new native GFN2 calls completed normally:32 initial SAD calculations and64
confirmed own-seed native restarts. Eight q0 analytic gradients were reused from
the separately owned fixed-pool continuation job. No MACE, DFT, optimization,
new structure/state or duplicate q0 call ran. All14predeclared derivative
quantities remain reported;7 pass and 7 fail. Neither complete source qualifies.
[Scope and gates](PLAN.md), [machine summary](RESULT.json), [commands](COMMANDS.md).

## Matched physical comparison

The existing source-connected maps rotate the complete terminal carboxylate of
4MAE Glu172 and Q88JH5 Glu221 about CG–CD. The4MAE mode is eligible in the full
physical map but is absent from the adaptive four selected modes; this was
explicitly resolved before execution. Both metals use identical coordinates
and the same preserved charge/protonation/context policy. The existing cap
chain rule, bonds, fixed atoms and overlap checks all pass. Maximum physical
heavy displacement is0.00110422Å. No independent cap motion is introduced.

At q=±0.001 and±0.0005rad, each cell starts fresh, then continues its own matching
GBW/xtbw twice. All native restart activations are confirmed from actual output.
The third energy is always used, even when higher. MaxIter 500, 300 K, native mixer,
Hamiltonian, parameters and solvent settings are fixed. The q0 analytical
gradient is projected through the same physical Jacobian, in kcal/mol/radian;
it is dE/dq, not force. Hartree and bohr conversions are applied once.

The gradient files' energies **exactly match** the final printed energies for
all eight q0 cells at retained precision. Native SCF/CN and appropriate ALPB
analytic drivers, atom ordering, coordinates, charge and electron counts pass.
This excludes an accidentally mismatched energy/geometry/gradient artifact.

## Numerical gates and derivative results

All 40 cells, including eight shared centers, pass the declared stage1→2 energy
change limit0.1kcal/mol. Maximum displaced-cell change is0.0643395kcal/mol.
All 10 matched Ca−La solvent contrasts pass0.2kcal/mol; maximum change0.0643872.
That scale is still far too large for a derivative across0.001rad.

Both centered-h versus centered-h/2 and analytical versus centered-h/2 must
agree within max(0.2kcal/mol/radian,5%of the finest derivative). Rounded values:

| Source | Quantity | Centered h | Centered h/2 | Analytic | Gate |
|---|---|---:|---:|---:|---|
|4MAE|Ca vacuum|0.85789|0.85789|0.85363|pass|
|4MAE|Ca ALPB|0.64575|0.64575|0.64624|pass|
|4MAE|La vacuum|9.98318|10.00507|10.00302|pass|
|4MAE|La ALPB|10.79001|−54.05485|11.24542|fail both|
|4MAE|Ca solvent transfer|−0.21214|−0.21214|−0.20739|pass|
|4MAE|La solvent transfer|0.80683|−64.05992|1.24240|fail both|
|4MAE|Ca−La solvent contrast|−1.01896|63.84778|−1.44980|fail both|
|Q88JH5|Ca vacuum|−1.25985|−1.52108|−0.80619|fail both|
|Q88JH5|Ca ALPB|0.41702|0.41703|0.41640|pass|
|Q88JH5|La vacuum|26.99745|27.16944|27.17974|pass|
|Q88JH5|La ALPB|33.30900|33.08575|33.24078|pass|
|Q88JH5|Ca solvent transfer|1.67687|1.93810|1.22258|fail both|
|Q88JH5|La solvent transfer|6.31155|5.91632|6.06104|fail refinement|
|Q88JH5|Ca−La solvent contrast|−4.63468|−3.97821|−4.83846|fail both|

All derivative values are kcal/mol/radian. The Ca−La analytical errors are
65.2976 and0.860248, against respective tolerances3.19239 and0.2. The moving
physical-heavy RMS speeds are1.088699 and1.100443Å/radian; dividing by these gives
**59.9776 and0.781729kcal/mol/Å** errors. All-heavy RMS speeds and actual atom IDs
are retained separately. Only the two terminal oxygens move in these modes.

## What the retained components show

The largest failure is 4MAE La/ALPB at −0.0005 rad. Its energies are:

| Stage | Energy (Eh) |
|---|---:|
|Fresh SAD|−352.318057069206|
|Continuation1|−352.318055466368|
|Continuation2, reported|−352.317952934843|

The final continuation raises this energy by0.0643395kcal/mol. Its printed maximum
and RMS density diagnostics remain large; maximum0.01237 versus displayed1e−5.
Nearby cells and the earlier stage give a much smoother finite difference.
They are preserved diagnostics, **not selected replacements**. Q88Ca vacuum has
smaller but still derivative-relevant sensitivity. The data show residual
solution/energy sensitivity after two continuations; they do not prove its
unique internal cause or the correct electronic root.

Fractional orbital occupations are present at the requested300K. No separate
electronic entropy/free-energy scalar is printed. Generic electronic+nuclear
entries do not sum to the native total, so their residual cannot be labeled
entropy, ALPB or another unique term. There is no justified entropy correction
here. The matching .engrad/printed energies verify artifact consistency; they
do not guarantee that an incompletely settled electronic solution is variational.
Native special-mixer convergence is not established solely by generic printed
density tolerances. All raw occupations, residuals, component energies and
input/output/gradient pins are retained in COMPONENT_AUDIT_v1.json.

## Cost, tests and disposition

Job 1210338 completed in**235s on 64 allocated CPUs = 15040 core-seconds**,128GiB
requested RAM, no GPU, node-128-512g-8gpu-1. Scheduler batch MaxRSS7908148KiB
(~7.54GiB). Stage executor wall times94.0102,61.1326,61.1823s are nested inside
that allocation and are not added again. Local preparation/testing/reporting is
unmetered. Eight shared q0 calls belong to the other experiment's recorded cost.

**Nine real-fixture tests pass**, zero skips: source/state/geometry/caps, exact
native recipes and own-seed links, real analytic-gradient parsing, signs/unit
conversion, corrupted-input rejection and explicit missing-gradient propagation.
These parser/algebra tests are distinct from the 96 executed scientific calls.
The first pytest invocation was unavailable in the pinned environment; unittest
ran the full tests without installing anything. Earlier test-harness typo and
all versioned logs remain visible.

New opt-in protocol: `native_GFN2_self_continued_solvent_force_check_v1`.
Production/default/calibration and all prior outputs remain unchanged. **Close
this continuation-as-force-remedy test with no further calls.** It establishes
that acceptable score-level settling alone does not qualify small-displacement
forces. It provides no classifier gain, geometry minimum or affinity/entropy
claim. Existing analytic-gradient capability and successful scalar-continuation
results remain useful within their separately demonstrated limits.

Full evidence: `workspaces/native_solvent_force_20260923/run_v1/COMPARISON_v3.json`,
`COMPONENT_AUDIT_v1.json`, `COSTS.json`, `RUN.json`, three stage manifests and all
receipts. The next runnable command is the read-only comparison in COMMANDS.md;
no new scientific execution is proposed by this report.
