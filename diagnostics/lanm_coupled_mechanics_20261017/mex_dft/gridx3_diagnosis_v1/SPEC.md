# Saved Mex final-integration and field-gradient diagnosis

Declared 17 October 2026 UTC after the original scout pair failed. Root owns execution. This is a saved-data diagnosis of that failure, not a replacement qualification test or changed molecular recipe.

## Question

Does the geometry-dependent post-SCF COSX/gridx3 correction account for the embedded-PBE0 directional mismatch, and does the archived raw point-charge gradient reproduce the assembled analytic value?

Source review observes final-grid corrections of -0.000019953 Eh (plus), -0.000018673 Eh (minus), and -0.000019569 Eh (central). These are native printed observations. Their directional contributions and residuals have not been recomputed in this source-preparation phase.

## Inputs and fixed scope

Use exact collector13271 RESULT334bb6b4, its actual two endpoint results/native outputs, the admitted central origin12849/native output/pcgrad, original prepared directions and physical/field mapping. Preserve all native source, state, method, geometry and field identities. The moved A/97//O is physical index1019 and external row832; core coordinates stay identical. The original displacement is +/-0.001 Å.

One one-CPU Slurm job reads these artifacts and performs:

1. Exact replay of the original final-energy finite difference and reported analytic comparison, using native-unit subtraction before conversion.
2. Extraction of native last-converged SCF energies, final-grid correction and final energy components. Verify declared fixed-core D4 equality before treating its derivative as zero.
3. A clearly labeled decomposition into the pre-final-grid energy derivative and final-grid correction contribution. Preserve printed decimal precision and report the decomposition/rounding closure; do not pretend printed corrections are exact beyond their digits.
4. Direct replay of central raw pcgrad row832 with the frozen source direction and Hartree/bohr conversion, alongside the assembled central component gradient. Check identity/order, sign and unit application; report any mismatch rather than silently repairing it.
5. Full reporting of the residual after accounting for the final-grid contribution, whether or not it becomes smaller. Keep the original failed PBE0/total qualification and all six unrun endpoint slots explicit.

No DFT, MACE, GFN, MM, ALPB, optimization, new geometry, different step or solver call is allowed. No scientific imports or numerical tests run on login. Use the original qualified Python environment with one numerical thread and default scheduler memory/time policy. Record actual command, source/config pins, scheduler cost and completion wake.

## Interpretation

A substantial final-grid contribution would localize a numerical energy/gradient consistency problem; it would not validate an alternative energy definition, a new gradient, physical accuracy, or electronic stability. The untouched final-energy derivative gate remains failed. A remaining discrepancy must be reported, not erased by redefinition or a new tolerance. Single-step truncation uncertainty and the original central/plus/minus SCF-state/residual/runtime limitations remain.

Saved replay success permits interpretation of the existing failure only. The remaining six qualification singlepoints and eight response references remain stopped. Any subsequent molecular protocol must be separately declared with its own unchanged test before execution; this diagnostic launches none.
