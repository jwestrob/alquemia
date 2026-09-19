# Approved PQQ contextual physical-H preparation — 2026-09-19

Root approved this exact follow-on while static context DFT remained live.
Question: can inexpensive physical-H preparation in the actual second-shell
context improve original-core PQQ discrimination, independently of adding context
energy? Cases1H4I/MxaF and4MAE/XoxF, bothCa/La, are consumeddevelopmentcontrols.

Four native unmasked OMOL-0 100M float64 constrained searches use the exact frozen
154/202atom contexts. Allheavy atoms and artificialcaps stay fixed. Physical
proteinH and canonicalPQQphysicalH may move within0.35Å of each startingcoordinate.
A single archived start, fixedSLSQP solver(max200iterations,ftol1e-9eV), and
projectedHgradient≤0.03eV/Å define the optimizer and constrainedstationarity.
Nearestheavycovalentparent/protonation mustremainunchanged. No added/removedatoms,
newwater, metalmovement, extraSCFparameter, or label-dependentselection.

The selectedproposal is the finalSLSQPiterate when coordinates/domain/chemical
identity pass. Convergence andboundaryflags are separatelyreported; an unfinished
or boundaryproposal is not claimed to be an unconstrainedminimum. The same rule
is used for everyendpoint, regardlessenergy/classification. Transfer only physical
H belongingtooriginal73/80atomcores back into their exact originalatomorder.

Four original-native-recipe r2SCAN-3c/CPCM/DefGrid3 singlepoints assess Ca−LaR and
XoxF−MxaFgap. Four additionalnativeOMOLcoreevaluations record thesameproposedcores.
ExistingoriginalcoreDFT/MACE energies are reused. No DFTOpt/gradients/Hessian.
No originalprotocolbands inherited afterHchanges; no baseline/defaultchange.

OneGPU16CPUs200000MiB, estimatedminutes; fourconcurrent16rankDFT endpoints on
64CPUs256GiB, estimatedminutes. Algorithmiterationlimit is a fixed solver setting,
not a compute/timebudget. Allattempts andfinalmeasuredcosts are retained.
