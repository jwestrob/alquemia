# Finite additive mechanics candidate: declaration before component evaluation

27September2026. Objective: make the missing protein mechanics/cross interactions
explicit for the actual1H4I source, not optimize or score biological preferences.
Use exported exact9113atom ff19SB System and PQQ3minus27atom GAFF2crossLJ table.
Physical system is protein+PQQ+exchangedmetal; electronic54atom model includes
three mapped syntheticcaps. All inputpositions and proton/redox state pinned.

Candidate energy ledger for review:
E=E_QM(physicalQM+mappedcaps;redistributedMMcharges,DoEQfalse)
 + E_FF_bonded(terms with at least one physicalMM atom)
 + E_MM-MM_LJ + E_MM-MM_Coulomb(redistributedcharges)
 + E_realQM-MM_LJ.
Exclude whollyQM classicalbonded/nonbonded terms, all syntheticcap classical
terms, QM-MM classicalCoulomb (already electronicembedding), allC4induction,
and any second solvent/dispersion term. MM-MM electrostaticexceptions must be
recomputed consistently with redistributedcharges and recordedoriginalscale;
QM-MM LJ honors actual protein1-2/1-3 exclusions and1-4exceptionparameters.
PQQ/metal are not covalentlybonded toprotein and have no inventedcovalentexceptions.
Crossbonded terms need full explicit atom supports; capdependentQM contributions
remain a boundaryapproximation, not exact physical cancellation. Review scheme
against documented additiveQM/MM; flag unsupported interactions instead of zeros.

Freeze metalcrossLJ candidate to installedpure12-6TIP3P ions234lm_126 parameterfile,
CaCM andLaIOD values documented in hybrid_feasibility/REPORT.md. This parameter
family is chosen before combinedenergy outcomes because it supplies explicit
pure12-6 coefficients, not a12-6-4 fit with induction removed. It remains an
unqualified hydration-derived approximation for embeddedQM metals; different
fitobjectives and possible effectivepolarization are limitations. Do not select
variants by desired classification or call it a validated LanM forcefield.
Use standardAmber combiningconvention fromactualforcefield, no inventedconstants.

This stage prepares ledger/sourceindices/exclusions and lowlevel componenttest
commands only. No quantumreruns, MD, optimization or newclassifier. Subsequent
finite lowlevel tests must verify actual energies/forces and sourceJmapping;
completecombinedforce qualification remains separate. NativeQM numericalgate
and exactMLcheckpointblockers persist. Production untouched.
