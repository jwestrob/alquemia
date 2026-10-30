# La/Dy charge-transfer diagnostic — inconclusive

**The bounded pilot did not complete. We cannot say whether suppressing
metal–shell charge transfer changes the La/Dy water-to-acetate differential.**
The remaining three systems were not run because the pilot was not admitted.

## Exact failure

The initial La-aquo attempt, Slurm **21618**, ran out of memory after 126 seconds
with 96,000 MiB allocated. One resource-only recovery, **21623**, used 24 CPUs
(12 MPI × 2 OpenMP), 256 GiB and `standard_expand`, with the same scientific input.
It reached the predeclared **7,073-second native-process limit** while still
optimizing the isolated charged-La fragment. The runner then sent SIGTERM to
its own process group. Slurm records FAILED, exit 1:0, 7,087 seconds elapsed.
This was not another OOM, a 300-iteration SCF failure, or a demonstrated
Hamiltonian incompatibility.

The first fragment had 8 electrons, 4 occupied orbitals and 47 orbital functions.
Its 43rd and last primary iteration printed residual **4e-7**, above the declared
**1e-10** threshold. No fragment converged; the shell, frozen complex, polarized
complex and full complex energies were not obtained. No R/F/P/U decomposition,
full-energy comparison or closure check can be admitted. The auxiliary atomic
initial guess and intermediate SCF energies are not completed fragment results.

The native command took 7,078.734 seconds including 5.734 seconds of termination
latency after the cap. Both molecular attempts together used 7,197.389 native
seconds; scheduler elapsed time also includes startup and cleanup. The timeout
exception occurs before the parser is called, so the missing parser `RESULT.json`
is expected. The preserved failure receipt is [RUN_RESULT.json](cp2k_results_memory_v1/aquo8_central_La/RUN_RESULT.json),
with [raw output](cp2k_results_memory_v1/aquo8_central_La/endpoint.out) and an
[independent terminal review](cp2k_memory_recovery_v1/TERMINAL_REVIEW.md).

## Scope and provenance

Exactly the four requested archived geometries were prepared. Their XYZ hashes,
unchanged coordinates and identical nonmetal coordinates within each La/Dy pair
were verified; finalization rechecked all 31 execution pins. No geometry was
optimized or regenerated. See [INPUT_AUDIT.json](INPUT_AUDIT.json),
[inputs/PROVENANCE.json](inputs/PROVENANCE.json), [FINAL_RESULT.json](FINAL_RESULT.json)
and the retained [goal wording](GOAL_SCOPE.txt).

| Geometry | Final execution status |
|---|---|
| `aquo8_central_La` | Initial OOM; unchanged-input recovery timed out; not admitted |
| `aquo8_central_Dy` | Prepared, unrun because pilot not admitted |
| `hydrated_acetate_mono_central_La` | Prepared, unrun because pilot not admitted |
| `hydrated_acetate_mono_central_Dy` | Prepared, unrun because pilot not admitted |

The archived protocol is PBE0 with La ECP46/Dy ECP55, archived orbital bases,
and effective singlets; physical Dy sextet metadata is distinct. Two fragments
are M3+ and the entire shell: 8 waters (shell charge 0), or acetate plus 7 waters
(shell charge −1). Complex charges are +3/+2 and explicit electron counts 88/110.
These states were bound in all four inputs; only the initial La fragment was
actually evaluated. Dy numerical compatibility remains untested.

The isolated, source-audited CP2K 2026.2 build uses its native ALMO implementation
and a reviewed seven-line Dy55 core-occupation reader repair. ECP coefficients,
integral routines and ALMO equations are unchanged. Startup/syntax success did
not establish scientific compatibility. [INSTALLATION.md](INSTALLATION.md),
[CONVERSION.md](cp2k_inputs_v1/CONVERSION.md), [RUN_SPEC.md](cp2k_inputs_v1/RUN_SPEC.md)
and [the resource-recovery specification](cp2k_memory_recovery_v1/RECOVERY_SPEC.md)
retain the exact implementation and settings. The archived light-atom basis is
recovered through checked named-library equivalence, not original printed
primitive coefficients; the full-energy compatibility test never finished.

## Requested results

The sign convention is
`C_X = (X_Dy − X_La)_mixed − (X_Dy − X_La)_aquo8`;
negative values favor Dy in this conditional substitution contrast.

| Quantity | kcal/mol | Status |
|---|---:|---|
| Archived full PBE0-D4 differential | +0.5750525133 | Existing reference; not a new ALMO result |
| Archived electronic differential | +0.4571786583 | Existing reference |
| Archived effective D4 differential | +0.1178738550 | Existing reference bookkeeping |
| Transferred AMOEBA differential | −6.2358150700 | Existing model |
| Model minus archived full reference | −6.8108675833 | Existing discrepancy |
| New full quantum differential | Unavailable | Pilot incomplete |
| CT-suppressed differential | Unavailable | Pilot incomplete |
| CT and polarization differentials | Unavailable | Pilot incomplete |
| Frozen-density contribution | Unavailable | Pilot incomplete |
| Native reference/preparation and reconstruction residuals | Unavailable | Required states absent |
| CT numerical uncertainty | Not established | No CT estimate |

Machine-readable values and explicit missing statuses are in
[RESULTS_TABLE.csv](RESULTS_TABLE.csv); unchanged archived values and caveats are
in [ENERGY_LEDGER.json](cp2k_inputs_v1/ENERGY_LEDGER.json). Missing values are not zero.

The planned bookkeeping remains explicit: R is the own-basis isolated metal
plus whole-shell reference; F is frozen, P polarized without metal–shell CT,
and U fully relaxed. FRZ=F−R, POL=P−F and CT=U−P. Full and suppressed totals
would add the same archived D4 scalar once to U and P. Free acetate/water terms
cancel in the metal double difference. Matched fragment references cancel
algebraically, but actual reference residuals were to be retained; none can
be measured here. No separate electrostatic/Pauli split, deformation energy,
missing absolute isolated-shell D4 reference, or numerical closure is invented.

## What this establishes

1. **Does CT materially alter or reverse the differential?** Unknown. No constrained
   or full complex state completed, so comparison with numerical uncertainty or
   the 6.81-kcal/mol model discrepancy is impossible.
2. **Which components reinforce or cancel?** Unknown. The existing full/electronic/D4
   bookkeeping is not an electronic decomposition into frozen, polarization and CT.
3. **What is established?** This exact CP2K setup did not finish its first charged
   fragment within the bounded runtime. The input/provenance work is reusable;
   scientific compatibility and the CT counterfactual remain unqualified.
4. **What remains unresolved?** Every requested new component and differential,
   full-energy agreement, decomposition closure, and numerical uncertainty.
   This neither identifies the force-field error as missing CT nor establishes
   a solution/protein binding free energy. The retained effective Dy state and
   archived convergence/integration limitations also remain.

## Actual allocation cost

| Job | Purpose / terminal state | Seconds | CPUs | Allocated CPU-seconds |
|---|---|---:|---:|---:|
| 21396 | First software build; failed | 1197 | 24 | 28728 |
| 21440 | Build recovery; completed | 700 | 24 | 16800 |
| 21618 | Initial pilot; OOM | 126 | 24 | 3024 |
| 21623 | Pilot recovery; process timeout | 7087 | 24 | 170088 |
| **Total** | | | | **218640** |

This is 60.73 allocated CPU-hours, not measured utilization. Recovery peak batch
RSS was 251,326,840 KiB, below its 262,144 MiB grant. [COSTS.json](COSTS.json) retains
parent and batch accounting separately; batch rows are not double-counted.
All owned jobs in this bounded diagnostic are terminal. No additional molecular
calls, fitting, protein work or follow-on campaign were released.

## Recommendation

**Stop this attempt as inconclusive. Do not use it to justify force-field changes
or automatically extend the CP2K campaign.** The CT hypothesis remains unresolved;
any further work requires a new, explicit decision outside this bounded diagnostic.
