# Existing natural-transfer inputs

19 October 2026. [INVENTORY.json](INVENTORY.json) binds existing metadata for the
three sources required by [SPEC.md](SPEC.md). No molecular evaluation, geometry
generation, scientific numerical analysis, search, submission or state selection
was performed. Metadata fidelity and readiness to score remain separate.

| Exact INPUTS source ID | Protein/water atoms | Waters | Local capped charge | Whole charge | Spectators | Deposited protein assembly |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `Hans8DQ2` | 1887 | 74 | −1 | +6 | La/La/Na | A+C or B+D dimers |
| `Hans_8FNR` | 2202 | 179 | −1 | +8 | Dy/Dy/Dy | A+B or C+D dimers |
| `Mex_8FNS` | 2045 | 167 | 0 | +2 | Nd/Nd/Nd | Monomer, chain A |

Counts and charges above are saved declarations, not recalculated molecular
properties. Each source retains its own topology, hydration, protonation and
spectator identities. Both Hans deposits are correlated sources of the same
protein; their differences cannot be described as a pure conformational effect.

## Artifact semantics that must survive reuse

- Prepared `physical_source_A.json` and the force-field export contain protein
  and water atoms only. The EF3 target and three spectators are separate INPUTS
  declarations. The bound classical-origin receipts contain the full physical
  system in their own order: protein/water, three spectators, then target.
  A core map instead places the target at QM index 0. No universal target index
  or interchangeable atom order is declared.
- INPUTS configurations A/B mean the original and recorded 2-degree carboxylate
  perturbation. They do not name biological chains. This contract selects A;
  B remains historical evidence.
- The original exported topology/charges and reviewed C-bound-H coordinates
  are separate pinned artifacts. The latter's review retains an unavailable
  optimizer-status receipt and its stated geometric-review policy. No new
  convergence or physical-state qualification follows from that admission.
- Historical INPUTS represent spectators as formal monopoles without ion LJ.
  Later classical ledgers include ion-specific parameters and a declared
  boundary policy. Both versions are identified; neither selects the energy
  model for future transfer or licenses adding classical terms to a whole
  electronic model.
- Local La singlet and Dy sextet declarations are copied exactly. Aggregate
  multiplicity, spin and physical-state qualification remain null for all
  three whole systems. The earlier Mex 10/15 and 4/9 panel is linked only as
  two scoped conditioning hypotheses, with neither selected here.

## Assembly and exchange limits

The inventory pins the original 8DQ2/8FNR/8FNS structures and copies the relevant
assembly metadata literally. All prepared INPUTS describe a conditional
extracted chain A with no dimer/folding equilibrium. Deposited partner chains
provide a possible no-new-fold route, but no partner, complex, charge/spin state
or water/ion ownership is prepared by this contract. Deposited fractional
occupancy is not a field in these INPUTS; the source structures remain pinned,
and no solution population is inferred.

The two symbolic balanced contrasts keep Hans8DQ2 and Hans_8FNR separate against
Mex_8FNS. Their fixed backgrounds are preserved within each La/Dy pair. These
are conditional EF3 exchange definitions, not affinities or a validated natural
monomer comparison. The consumed experimental-target note is bound as context,
not an EF3 label or fitting rule.

## Validation and readiness

The inventory creation verified 106 explicitly included unique file hashes and
the named INPUTS/export/core/source joins. The archived classical La/Dy records
have identical source-ID and coordinate lists within each source. This used
only standard-library metadata operations; no scientific modules or numerical
arrays were evaluated. The independent checker and review are separate files.

`scoring_ready` remains false. Mechanics accuracy/response, a common energy
recipe and explicit whole-system states remain unresolved. Contract validation
does not reopen the closed electronic-base candidate, change running reference
gates or authorize transfer scoring.
