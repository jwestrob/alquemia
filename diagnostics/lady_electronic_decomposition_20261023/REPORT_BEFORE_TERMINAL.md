# La/Dy electronic diagnostic — compatibility pilot running

**No charge-transfer differential is available yet.** The exact La-aquo pilot
is running as Slurm **21623**, using 24 CPUs (12 MPI ranks × 2 OpenMP threads),
256 GiB and `standard_expand`. Direct scheduler and output inspection on
23 October 2026 at approximately 21:39 UTC confirmed live charged-La fragment
iterations. Neither fragment convergence nor the full ALMO chain is established.

The initial attempt, **21618**, ended out of memory after 126 scheduler seconds
with 96,000 MiB allocated. Its raw output and failure are preserved. Recovery
21623 changes only the memory request and execution paths; the electronic deck,
binary, geometry, thresholds and admission criteria are unchanged. Its native
cap is 7,073 seconds, after charging 127 seconds against the original 7,200-second
pilot budget. No further automatic escalation or restart is authorized.

## Question and exact scope

At fixed geometry, does suppressing metal–whole-shell charge transfer materially
change the La/Dy energetic contrast for replacing coordinated water with acetate?
The two fragments are M3+ and the complete ligand shell. The aquo shell is neutral;
the acetate plus seven-water shell has charge −1.

| Archived geometry | Current result |
|---|---|
| `aquo8_central_La` | Pilot running; no admitted energy decomposition |
| `aquo8_central_Dy` | Prepared; requires pilot admission |
| `hydrated_acetate_mono_central_La` | Prepared; requires pilot admission |
| `hydrated_acetate_mono_central_Dy` | Prepared; requires pilot admission |

All four archived XYZ hashes and exact matching nonmetal coordinates were
verified in [INPUT_AUDIT.json](INPUT_AUDIT.json). Original inputs and provenance
are retained under [inputs/](inputs/). No coordinates were optimized or regenerated.
The archived PBE0 protocol uses La ECP46 and Dy ECP55 with effective singlets;
physical Dy sextet metadata is not the effective calculation state. Aquo complexes
have charge +3 and 88 explicit electrons; mixed complexes have charge +2 and 110.

## Implementation and acceptance

The isolated CP2K 2026.2 implementation uses native frozen, block-polarized and
fully relaxed electronic states. A reviewed seven-line Dy55 input-reader patch
adds the core occupation case; ECP coefficients, integrals and ALMO equations
are unchanged. Installation and syntax checks passed. These do not establish
scientific compatibility. Exact assets, source audit and invocation are recorded
in [INSTALLATION.md](INSTALLATION.md) and [cp2k_inputs_v1/](cp2k_inputs_v1/).

Every cell must pass native convergence and state checks, full electronic-energy
agreement with its archived reference within 1e-4 hartree, decomposition and
reference-sum consistency within 1e-8 hartree, and raw-output review. The three
remaining systems will run concurrently only after the pilot passes. These
thresholds and the finite runtime were declared before execution; no result has
been admitted or used to alter them.

## Energy bookkeeping and existing evidence

For each state or component X, the reported contrast will be
`C_X = (X_Dy − X_La)_mixed − (X_Dy − X_La)_aquo8`.
Negative values favor Dy in this conditional substitution comparison.

Native states are R (isolated metal plus whole-shell reference), F (frozen),
P (polarized without metal–shell CT), and U (fully relaxed). Components are
`FRZ = F − R`, `POL = P − F`, and `CT = U − P`. Full and CT-suppressed totals
will use U and P plus the same archived, geometry-dependent D4 term. Actual
reference double-difference residuals will be retained rather than forced to zero.
Free acetate and displaced-water terms cancel in the metal double difference.

| Quantity | kcal/mol | Evidence |
|---|---:|---|
| Archived full PBE0-D4 contrast | +0.5750525133 | Existing reference quartet |
| Transferred AMOEBA contrast | −6.2358150700 | Existing model comparison |
| Model minus archived reference | −6.8108675833 | Existing comparison |
| New full quantum contrast | Unavailable | Pilot incomplete |
| CT-suppressed contrast | Unavailable | Pilot incomplete |
| Polarization, CT and frozen contrasts | Unavailable | Pilot incomplete |
| New reconstruction residuals | Unavailable | Pilot incomplete |

The exact source values, D4 separation and archived convergence caveats are in
[ENERGY_LEDGER.json](cp2k_inputs_v1/ENERGY_LEDGER.json). No new numerical uncertainty
or CT magnitude can yet be inferred. Arithmetic closure alone would not establish
SCF accuracy. These are nonequilibrium cluster energies, not solution or protein
binding free energies; a future CT contribution would not by itself identify the
force-field error as missing CT.

## Execution records and next decision

The current pinned execution is
[CONFIG.json](cp2k_memory_recovery_v1/CONFIG.json); the rationale and finite limit
are in [RECOVERY_SPEC.md](cp2k_memory_recovery_v1/RECOVERY_SPEC.md). Actual live
receipts and raw output are under
[the pilot task](../../workspaces/lady_electronic_decomposition_20261023/cp2k_quartet_memory_v1/tasks/aquo8_central_La/).
[CURRENT.md](CURRENT.md) tracks the latest reviewed state. [COSTS.json](COSTS.json)
records 48,552 allocated CPU-seconds for terminal setup and the failed initial
attempt; it excludes the still-running recovery and is not measured utilization.

**Current recommendation:** finish this unchanged, bounded compatibility pilot
and use its declared acceptance checks to decide whether the three remaining
systems can run. No scientific mechanism recommendation is supported yet.
The final deliverable will replace this interim status after actual terminal
review. Earlier software-investigation and status text is preserved in
[REPORT_HISTORY_BEFORE_PILOT_RESULT.md](REPORT_HISTORY_BEFORE_PILOT_RESULT.md).
