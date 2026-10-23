# Matched water radial force qualification

## Question and scope

Does the choice of uncorrected PBE0-D4 forces account for the apparent La-water force-field mismatch? This finite diagnostic compares methods and counterpoise at the exact existing water geometry. It does not fit parameters, establish a physical reference, inspect reserved validation outcomes or release protein sampling.

## Frozen inventory

Use the original paired central La/Dy water coordinates. Define u = (O − metal)/|O − metal|. Keep the metal fixed and translate every water atom by q u, for q = −0.010, −0.005, +0.005 and +0.010 Å. Positive q increases the metal–oxygen separation; all water internal distances remain fixed. Four displacements × two metals × three counterpoise roles × two methods gives exactly **48 new single points**. No central point is recalculated.

The three roles are real complex (+3, 18 explicit electrons), real ion with ghost water (+3, 8 electrons), and real water with ghost ion (0, 10 electrons). All have the same four centers and 90 orbital basis functions. La/Dy replace 46/55 core electrons respectively; both real ions have effective nuclear charge 11. Ghost metal retains its orbital basis and explicitly removes its ECP. Physical La singlet and Dy sextet metadata remain distinct from these effective frozen-f singlet calculations.

Canonical RHF CCSD(T), NoFrozenCore, def2-TZVP light atoms and the existing external metal ECP/basis use the unchanged proven water/donor recipe, VeryTightSCF, PModel, SCF limit 300, CC limit 100 and CC residual target 1e−8. PBE0 uses the same bases, RIJCOSX, DefGrid3 and SCF controls. PBE0 **real complexes alone include D4**; all ghost calculations omit D4. Canonical CC has no D4. The free real-water D4 energy is unchanged by translation and cancels in derivatives; ion D4 is the exact empty-pair/triple sum zero. No D4 calculation with ghosts or additional fragment calculation is needed.

## Observables and decisions

For each metal and h = 0.010 and 0.005 Å, report [E(+h) − E(−h)]/(2h) in kcal/mol/Å for:

- CC complex energy and CC counterpoise interaction E(complex) − E(ion in full basis) − E(water in full basis).
- PBE0 electronic complex and counterpoise interactions.
- The real-complex D4 term separately, then PBE0-D4 complex and counterpoise derivatives. Constant free-fragment dispersion cancels.

Compare with the saved central PBE0-D4 analytic derivative and native AMOEBA analytic derivative, both computed as **−sum(water forces)·u**. Also retain the total force sum and metal projection so translation closure is visible. Report Dy−La derivatives and model-minus-reference differences without confusing the coordinate derivative with the oppositely signed force.

The declared diagnostic derivative resolution is **0.1 kcal/mol/Å**. Two-step differences at or below this value are resolved at this diagnostic scale; larger differences remain unresolved, with both estimates retained. Compare each method/AMOEBA residual at both steps. Use the smaller-step result as the prespecified summary only when its two-step comparison passes. No extrapolation, fitted tolerance or third step is allowed. PBE0-D4 complex FD versus saved central analytic agreement is a separate sampled consistency check at the same resolution. Method disagreement does not invalidate a numerically admitted native result. No conclusion that CC is exact follows from either agreement or disagreement.

## Numerical evidence and release

Retain raw native output, all printed SCF flags, active convergence mode/stop signal, CC residual trace/T1/triples, and all energy-stage/final-integration/COSX metadata. In particular, DefGrid3 alone does not establish that the reported final energy and analytic gradient describe an identical surface. Preserve central EnGrad versus displaced SP declarations and flag changed stage signatures; the actual finite difference tests sampled consistency. Density/DIIS flags are not retroactively made mandatory when native SOSCF uses energy or orbital-gradient convergence.

Each cell must pass the existing native method/state, center/ghost/ECP, basis, electron, normal termination and CC residual checks. Prepare frozen coordinates and saved-output/failure-control tests before root submits. First scout is cc__La__plus_010__complex, the largest water electron inventory (18 electrons, 90 AOs); reuse that exact result in the 48-cell total. Root alone decides whether actual admission and memory/cost permit the remaining 47 cells. Four MPI ranks, one thread each, standard_expand, scheduler-default memory/time and MaxCore 75% of the actual memory grant. No retry, extra step, solver change or fit is implicit.
