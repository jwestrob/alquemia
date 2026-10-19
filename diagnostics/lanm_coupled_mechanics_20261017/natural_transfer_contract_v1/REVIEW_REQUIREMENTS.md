# Natural-transfer state and assembly requirements

19 October 2026. Independent source-only conceptual review of [PLAN module 5](../PLAN.md), the [requirement matrix](../plan_evidence_20261019_v1/REQUIREMENT_MATRIX.md), [consumed experimental target](../../lanm_coupled_response_20261016/EXPERIMENTAL_TARGET.md) and [second-Mex feasibility](../../lanm_coupled_response_20261016/MEX_SECOND_SOURCE_FEASIBILITY.md). No external search, geometry generation, numerical analysis or scoring was performed. Harvey owns exact source inventory and pins.

**A faithful metadata contract is feasible now. Transfer scoring is not qualified by that contract:** mechanics/response gates remain unmet, and whole-system electronic states are not supplied by the historical local inputs.

## Within-source EF3 exchange

For a declared source state `s`, define `Delta_s = E_s(Dy at EF3) - E_s(La at EF3)`. Each pair must preserve:

| Invariant | Required binding |
|---|---|
| Physical identity | Exact deposit/model/author-chain mapping, construct, atom order/IDs, topology, terminal chemistry, protonation/hydrogens and preparation lineage. Only the declared target element/state changes. |
| Occupancy and hydration | Exact EF3 identity; every spectator's identity, species, modeled occupancy, charge and coordinates; all water identities/coordinates and their topology. Deposited fractional occupancy and the prepared integer hypothesis stay distinct. |
| Geometry and response | At a fixed point, the same physical coordinates for both metals. A later response comparison requires its own shared proposal pool and identical search/selection rules within that source; a common pool cannot mix source inventories. |
| Hamiltonian | One declared energy recipe, model/parameters, solvent convention and component ledger. Historical formal-charge spectators, later ion-specific classical terms and quantum spectators are different treatments and need explicit versions. |
| Charge and spin | Local capped charges differ from whole charges: Hans local −1 versus Mex local 0; prepared whole +6/+8/+2 respectively. Local La1/Dy6 multiplicities are not whole-system multiplicities. Aggregate spin is absent from these INPUTS and must remain null/unassigned until separately declared. Prior Mex 10/15 and 4/9 panels are scoped hypotheses, not defaults or unique coupled states. |
| Boundaries and evidence | If a regional model is used, bind exact real atoms, caps/anchors, external field, exclusions and force mapping. Preserve native/electronic caveats, unavailable endpoints and qualification status; a source match does not qualify a new Hamiltonian. |

INPUTS configuration labels `A/B` mean original versus the recorded carboxylate perturbation; they are **not** biological chain identifiers. Keep these namespaces separate in the validator.

## What can cross sources

The balanced contrast is `R(s,t) = Delta_s - Delta_t`: equivalently one EF3 Dy replaces La in `s` while La replaces Dy in `t`. All four endpoints must be available under consistent method/state conventions. Each source's background cancels only within its own pair; the two source inventories need not be identical for this conditional exchange to be stoichiometrically balanced. Different hydration, occupancy and assembly contexts nevertheless remain part of its interpretation.

Report Hans8DQ2-versus-Mex8FNS and Hans8FNR-versus-Mex8FNS separately. Do not combine La from one Hans source with Dy from the other, minimize raw energies across them, average them as independent proteins, infer populations, or select the favorable contrast. Hans8DQ2 retains La/La/Na spectators and 74 waters; Hans8FNR retains Dy/Dy/Dy and 179; Mex8FNS retains Nd/Nd/Nd and 167. Their differences cannot be labeled a pure conformational effect. Absolute single-source `Delta_s` is not an aqueous binding free energy; isolated-ion reference cancellation in the balanced cycle does not supply missing protonation, hydration, occupancy or assembly thermodynamics.

The consumed pH-5 apparent Hans-versus-Mex selectivity is a protein-level developmental target, not an EF3 label, a fitting rule or a quantitative cross-assay calibration. Two Hans deposits remain correlated states of one protein; LBT and a second Mex conformation add no independent natural protein.

## Future partner-contribution cycle

Prepared chain A is a conditional extracted subunit. Deposited Hans assemblies are 8DQ2 A+C/B+D and 8FNR A+B/C+D; Mex8FNS assembly 1 is monomeric. Neither Hans extraction demonstrates an equilibrated dissociated monomer. Assembly operators, partner identity and alternative-chain occupancy must be bound explicitly, without silently choosing a different dimer.

A later fixed-geometry cycle may define `I_s(M) = E(A_M+B) - E(A_M) - E(B)` and `I_s(Dy)-I_s(La)`, changing only the chosen A-chain EF3. It requires exact complex/fragment atom conservation, explicit water/ion ownership, unchanged protonation and partner occupancy, compatible charge/electronic states, and one solvent rule for each physical assembly. Complex solvation is not automatically the sum of fragment solvation. No water deletion, added caps or fragment charge/spin defaults may close a missing definition. Exchanging both subunits or relaxing fragments is a different cycle and must be declared separately. This conditional partner contribution is neither a dimer dissociation free energy nor a prediction of dimer populations or metal affinity.

## Validator outcome

Separate **source-contract fidelity** from **scoring readiness**. Require exact pins and semantic joins; retain unresolved aggregate spin, partner preparation, hydration ownership or protonation choices as explicit unavailable fields with reasons. A successful metadata validation must not fill them, claim physical-state qualification, authorize scoring, or reopen the failed electronic-base candidate. No independent natural monomer-only comparator is established by the present three-source contract.
