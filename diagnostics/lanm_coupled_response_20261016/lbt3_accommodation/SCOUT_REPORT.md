# Intact LBT3 accommodation scout

The first 7CCO/La search lowers the complete energy by **2.880329 kcal/mol**
with unchanged covalent geometry and no newly severe contacts. It is an admitted
finite candidate, not a stationary minimum or a preference result. The unchanged
three remaining source/metal searches pass the execution/physical continuation
gate; there is no favorable-sign criterion.

Nine accepted steps used nine new MACE force evaluations and eighteen standalone
GFN2 analytic energy/gradient evaluations. The MACE contribution is −2.701812
kcal/mol; solvent adds −0.178516. The accepted geometry moves at most 0.195947 Å
from the source and the metal moves 0.114032 Å. All four final trial multipliers
exceeded the frozen 0.20 Å source displacement envelope, without new energy calls.
The final active projected maximum force remains 6.15677 kcal/mol/Å. Thus this
is an envelope-limited search, even though its final coordinates sit just inside
the numerical boundary flag. Omitted modes and full residual forces stay recorded.

Slurm12533 completed in107 s on12 CPUs and one H200:1284 allocated CPU-s and107
allocated GPU-s. Application elapsed104.760 s; model load6.465 s; nine MACE calls
3.843 s; eighteen GFN calls137.517 summed worker-seconds. These include parallel
CPU calls and are not additive latency. Scheduler TotalCPU149.492 s;
MaxRSS1507712K. Peak CUDA allocation7975476736 bytes. No speedup claim.

Exact source, state, recipe and physical guards are in the frozen
LBT3_ACCOMMODATION_PLAN.md and workspace CONFIG.json. See SCOUT_RESULT.json and
SCOUT_ACCOUNTING.tsv. No PQQ change, equilibrium inference, new experimental
label or natural-protein validation follows from this peptide scout.
