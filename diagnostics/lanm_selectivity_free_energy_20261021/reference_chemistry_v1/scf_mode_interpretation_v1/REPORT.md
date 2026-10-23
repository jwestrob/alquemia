# What the printed SCF flags mean

**These 14 calculations converged under their declared ORCA stopping logic.**
The earlier raw residual tables remain correct, but requiring all six printed
thresholds would impose an undeclared stricter test. They do not establish a
native convergence failure or a measured numerical error.

ORCA's VeryTightSCF preset uses ConvCheckMode2: the total and one-electron energy
checks, rather than requiring every printed density/DIIS/rotation threshold.
The manual distinguishes this from mode0 and explains why density precision
can differ. SOSCF changes the iteration variable to orbital rotations and
gradients. [ORCA6.1 convergence and SOSCF documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html).

## Actual output evidence

All six canonical cells and all eight hydrated references explicitly print
`Total+1el-Energy`. All finish in SOSCF; 13 print an energy stopping signal and
the aquo8 Dy cell prints a gradient stopping signal. No TRAH recovery is used
in this set. In all 14, the reported final DIIS value matches the **first**
SOSCF gradient at the printed precision, while the final orbital-gradient
field matches the final SOSCF iteration. Thus that DIIS field is stale for
assessing the terminal solution. The density changes are current quantities,
but are not independently mandatory under the actual stopping mode.

Example: aquo8 Dy reports DIIS 1.4382e-3, matching initial SOSCF MaxGrad 1.44e-3;
the terminal MaxGrad is 2.15e-7 and triggers gradient convergence. Its energy
change being above TolE does not negate that native stopping path.

The canonical CC iterations separately meet STol1e-8: final residuals range
2e-9 to7e-9. Neither the SCF nor CC stopping tests prove accuracy of the
Hamiltonian, basis, effective-core treatment, or sampled free energy.

## Bounded precision check if a quantitative reference is needed

Do not repeatedly enumerate inactive flags. After the matched CP comparison,
test numerical sensitivity of the observable directly: a separately declared
repeat of its six CC and six PBE0 electronic CP cells with the same geometries,
bases, ECPs, grids and CC settings, but a stricter SCF orbital-gradient target
(e.g. explicit TRAH/TolG1e-8), seeded from the matching saved orbitals. Preserve
the original results and check orbital occupations/root continuity. Reuse D4.
Predeclare acceptance of <0.025kcal/mol change in each method's paired contrast
and <0.05kcal/mol in either individual interaction. These allocate a small part
of the existing0.25kcal/mol exploratory resolution, not an experimental error bar.

For the hydrated amide exchange and forces, the corresponding finite test is
four tighter energy/gradient cells: amide and aquo8, each La/Dy. Require the
paired contrast shift<0.025kcal/mol and each atom-vector force RMS change and
the paired ligand-translation projection change<0.05kcal/mol/A. The vector
RMS means sqrt(sum_i |delta F_i|^2/N). This tests numerical sensitivity only.
No such calculation is released or run by this report; root sets the priority
after the current matched comparison. Frozen-f, basis and ensemble uncertainty
cannot be resolved by tightening SCF.

Machine-readable evidence: RESULT.json (SHA256
`c2064671c80e6a6636bc728953d13076bd3e9d614cc5dac7da326491b33dc471`).
This also retains all largest-amplitude lines, including ORCA's negative-index
sentinels for singles that the earlier supplemental list omitted. Raw native
outputs and both original reviews are preserved.
