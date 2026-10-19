# Independent contract integration review

19 October 2026. **PASS for the frozen source-only contract, checker source and actual metadata-check result.** No material semantic or field/logic error was found. This review inspected source text, saved declarations, file hashes and root's completed receipt; it did not run or rerun the checker, scientific code, molecular calls or jobs.

## Reviewed pins

- SPEC.md: `c348744d96d158b6d8f8d73274cfa4e0144767ae7b0c59550d2bb2a19962e661`.
- REVIEW_REQUIREMENTS.md: `632e06de4a19e372e7a8d146ebabed895606959ca7006f951d50380e8094aabe`.
- INVENTORY.json: `e7d7ce48a27e02e1f5fb5392682943a55821ecd8d5cfd529510bfb3afef3a72e`.
- INVENTORY_NOTES.md: `1fd86ebd1a30d924fb4cb08ab1a101ff4a6987ae93452e3def1a31add7a7a49c`.
- CONTRACT.json: `9788d2f5134b5688bdd3a0dbf94bf2f8b6c34d2e5987c363492513b86707b5a3`.
- CHECK_CONTRACT.py: `ae7c3e30d59578bfcdefb37383f57bc666ad902351fe06e33b878f72b4572eea`.

CONTRACT binds the final inventory, specification and requirements exactly. The superseded inventory `fa5e150e…` is not its authority.

## Physical scope and sign

The two separate contrasts use `Delta_s = E_s(Dy_EF3) - E_s(La_EF3)` and `Delta_Hans - Delta_Mex`. This is the reaction Hans–La + Mex–Dy → Hans–Dy + Mex–La, so positive contrast conditionally favors the starting Hans–La/Mex–Dy assignment. Both endpoints of each source retain their own background. Neither differing inventories nor the experimental target justify combining endpoints from different Hans sources, selecting a favorable source, assigning populations or treating these as aqueous binding free energies.

The inventory preserves the essential distinctions: physical_sources/export are protein and water only; classical origins add three spectators then the target; core maps put the target at QM index zero. No universal full-system target index is invented. Source target IDs retain deposited La/Dy/Nd identity while endpoint substitutions are separately La/Dy. Configuration A/B means original/perturbed geometry, not biological chains.

Local charges −1/−1/0, whole charges +6/+8/+2, source-specific water/spectator inventories and local singlet/sextet declarations remain separate. Aggregate spin, whole Hamiltonian, partner preparation and physical-state qualification remain unassigned. Historical Mex spin hypotheses are linked without selecting one. Historical formal monopoles and later classical LJ terms are labeled as different stages rather than silently combined.

Hans extracted chain A remains distinct from the deposited dimers and from an equilibrated monomer. The final Mex wording correctly preserves its deposited monomer annotation while declining a solution-population claim. Future partner cycles still require explicit complex/fragment identities, solvent treatment, water ownership, charge and electronic states. No partner coordinates or score are prepared.

## Checker coverage and limits

The checker uses only standard-library hashing, JSON/text reads and exact metadata comparisons. It verifies declared file pins without expanding unrelated ancestors, all three ordered INPUTS identities, A/B namespace, actual local endpoint declarations, whole-versus-local charge/state separation, water/spectator records, physical/export/core/boundary joins, target row zero, literal assembly excerpts and separate contrast operands. It requires null whole-state/recipe/partner fields and false scoring readiness. Errors retain a failed result; exclusive output creation prevents overwriting a prior receipt.

Symbolic sign and interpreted assembly semantics were independently checked in this review. The program checks their operand identities and literal evidence; it is not a general physical-state or geometry validator. A successful execution establishes fidelity of this reviewed metadata snapshot, not new topology, cap derivatives, convergence, force accuracy, a protein-level preference or permission to score. No tests or numerical routines are invoked by importing the checker.

## Actual metadata receipt

Root's single execution completed with exit 0 and empty stderr. METADATA_CHECK.json SHA `5946410f6653724e88c07572b7ad789a022b47d4fbc23ead7efff90c7f43bfe2` records `metadata_fidelity_passed`, all 298 checks true and 107 unique declared pins. Its inventory, contract and checker hashes exactly match this review. CHECK_EXECUTION.json SHA `90567cd8c80bd085c369528e8b435fdd25cbe93fd52eb59d39bbd6dceaa1e12c` preserves the exact command and matching stdout. The 107 pins include the inventory itself through CONTRACT, alongside the 106 inventory declarations. Scoring readiness remains false, with zero molecular calls and zero scientific allocations.

The accompanying report is consistent with this scope. The failed electronic-base candidate remains closed, and the ongoing independent reference retains its existing admission/cost gates. This receipt closes metadata preparation only.
