# Saved-endpoint solvent-response comparison

Declared 2026-10-13 while the seven remaining reference endpoints run. This join uses the exact eight-cell dispatch in `REFERENCE_DISPATCH_v1.json`; no optimization, energy call, model invocation, refit or production change. Root submits all scientific execution and tests on one Slurm CPU.

## Inputs and admission

Use the accepted sample1 La/water endpoint from native3354, collected by the corrected3370 parser. Join seven new cells only from their assigned manifests and3371–3374 lanes. Each lane's final collection takes precedence; its own immutable per-task collection is usable if final collection is absent. Never borrow a same-named task from another manifest's unused copied directory. Preserve a valid first endpoint if its later lane endpoint fails.

Verify exact protocol/state/charge/multiplicity,192-atom coordinates,14caps, dry inventory, original input hashes and same source maps. Reparse raw final single-point energies and Mulliken/Loewdin populations using the frozen parser; verify receipt success, raw output hash, native version/basis/La46-electron ECP, expected molecular job and electron counts. Retain solver residuals and current diagnostic admission without claiming tighter density-response qualification. A failed required cell leaves only its dependent contrasts unavailable; preserve the complete eight-cell denominator. Future result files receive actual hashes when read, not fabricated advance pins.

## Algebra and charge interpretation

For each backend independently, U_M(s)=E_water,M(s)−E_vacuum,M(s). Report each U, U_Ca−U_La per source, and D=[U_Ca(1)−U_La(1)]−[U_Ca(3)−U_La(3)]. Also report each metal's sample1−sample3 water, vacuum and transfer response. All subtraction uses unrounded Hartree values; display conversion is applied once using the existing pinned627.509474 kcal/mol/Hartree convention. Printed CPCM dielectric energy is an auxiliary component, never U. No old aquo reference or classifier bands enter.

Compare against the actual saved standalone xTB6.7.1/GFN2 ALPB-vacuum endpoints, with independent tblite component qualification inherited from their pinned completed collection. Verify exact state/XYZ and stored endpoint identity. DFT/CPCM and GFN2/ALPB change both electronic and solvent models; agreement or disagreement cannot uniquely diagnose either component.

Partition populations by actual source-to-QM maps into target metal, complete27-atom PQQ, synthetic caps and each retained source-residue contribution. These partial-region residue groups are not full-protein formal charges. Keep Mulliken and Loewdin separate. Report within-scheme water−vacuum and source changes only when both endpoints exist. Require coordinate/element/map closure and population sum within the worst-case six-decimal printing bound192×0.5×10^-6 e +10^-9 e; this is a serialization bound, not a chemistry acceptance criterion or a cross-scheme comparison.

## Real partial-fixture test

A separate `partial-test` operation uses only the immutable3370 collection: one accepted sample1 La/water cell and seven explicitly unresolved cells. It must preserve that raw energy/population, leave every incomplete U and D unavailable, and reproduce the already-saved GFN2 double-difference algebra. No invented completed DFT outputs, dummy data or scientific executable replacement. Run the test on Slurm before the final join; AST/hash/config preparation alone is not a test pass.

The final collector runs afterany of3371–3374 so failures remain visible. It creates fresh result/receipt files, retains unrounded endpoints/group populations/source pins and actual execution receipts, and reports partial coverage honestly. No molecular retries or new work are triggered by this program.
