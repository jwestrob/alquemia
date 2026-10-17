# Nikasha coupled mechanics implementation plan

Prepared 17 October 2026 UTC following Jacob's request for a plan and a goal to activate. **Goal activated by Jacob on 17 October 2026 UTC; staged implementation and execution are authorized.** Root owns integration and submissions. Preserve PQQ production and its completed canonical-midpoint decision. Read CURRENT.md and actual receipts for execution status.

## Question and starting evidence

Can a coherent model of local electronic chemistry, protein mechanics and solvent produce credible La/Dy accommodation when the entire represented system can respond?

Separate two explanations: restricted motion prevents accommodation, and the energy model assigns the wrong cost to motion. Change accessible coordinates and the Hamiltonian in separate comparisons. Use existing structures and states. The rejected clamped-core/internal-strain diagnosis is not part of this plan.

The completed 48-atom LBT reference supports a large local source penalty: PBE0 gives La/Dy source strains +69.911132/+3.584798 kcal/mol, compared with MACE +77.640037/+3.950156. Flatter GFN2 energies are therefore not evidence of better physics. This result does not qualify whole-protein mechanics or affinity. Retained SCF and frozen-f limitations are in [the reference report](../lanm_coupled_response_20261016/lbt_local_reference/FINAL_REPORT.md) and [force comparison](../lanm_coupled_response_20261016/lbt_local_force_comparison/REPORT.md).

## Initial systems and denominator

| System | Physical system | Initial searches | Current energy |
| --- | --- | --- | --- |
| LBT3, 7CCO and 7CCN | 256 atoms per source, intact peptide and target | Two sources × La/Dy = four | Whole-peptide MACE vacuum + fixed-FF-charge ALPB |
| Mex-LanM, 8FNS | 2,049 atoms, including waters, EF3 target and three Nd spectators | One source × La/Dy = two | Local MACE + GFN2 embedding difference + MM ledger + fixed-FF-charge global ALPB |

7CCN is Lu-conditioned, not Dy-conditioned. LBT is an implementation/source-consistency control, not a second natural LanM. Mex's conditional EF3 exchange with Nd spectators is not a pure-La/pure-Dy protein binding experiment.

Start from the exact admitted terminal geometries of the completed restricted searches. Bind their actual manifests/results before execution. Preserve original coordinates and existing common pools; do not repeat searches or pool different source inventories. The source anchor is commit `32deb90`; pin reviewed code hashes because the shared tree contains unrelated work.

## Module 1 — Full physical energy and gradient

Expose one small contract:

```text
evaluate(physical_positions, declared_state)
    -> total_energy, physical_gradient, component_energies, provenance
```

Reuse the following implementations:

| Existing code | Reuse |
| --- | --- |
| `LBT_MATCHED_RESPONSE.py::MatchedEvaluator.evaluate` | Arbitrary 256-atom coordinates and composite gradients. |
| `MEX_RESPONSE.py::prepare_cells` | Electronic coordinates, external charges and caps from physical positions. |
| `scripts/metal_environment_force_assembly.py::scatter` | Local, cap and point-charge derivatives mapped to physical atoms. |
| `scripts/metal_environment_component_checks.py::build_system` and `MEX_CONTINUATION.py::mm_stage` | Existing MM interactions, exclusions and forces. |
| Existing L-BFGS direction helper | Optimizer kernel, without inherited rigid-bond and cumulative source-displacement restrictions. |

The diagnostic Python files are in `diagnostics/lanm_coupled_response_20261016/`. Mex retains:

```text
E = MACE_vac(Q) + GFN2_emb(Q;MM) - GFN2_vac(Q)
  + MM_ledger(R) + ALPB(fixed_FF_charges;R)
```

Every moved coordinate updates all dependent terms: caps, field, MM and solvent cavity. Spectator identities/states stay fixed; positional freedom is explicit. Moving fixed charges is not environment polarization. Cache identities include the complete physical/environment coordinates, state, model and protocol.

Keep the qualified corrected-La and Dy parameter lineage, imported files and hashes explicit in every GFN2 execution. Native executable defaults are never a fallback. Record total charge, physical multiplicity, backend spin convention and any ECP/f-in-core treatment separately.

**Gate:** exact origin replay, cache invalidation, component closure, and small directional finite differences of the complete gradient. Include newly enabled backbone, boundary-anchor, water and external-charge directions. Reuse compatible earlier checks. Freeze numerical tolerances before execution that resolve the intended convergence criterion. Numerical consistency does not establish physical accuracy.

## Module 2 — Broad response under unchanged energies

Use limited-memory optimization with analytic total gradients and adaptive per-step trust control. Permit backbone, sidechain, bond/angle, metal and represented environment motion. Remove only irrelevant global rigid motion. Keep topology, proton inventory, elements/states and cap mapping fixed.

Physical step-admission checks preserve connectivity, stereochemistry and boundary validity and reject singular contacts. Do not add springs to the source, freeze every inherited bond length, or impose the former 0.2/0.4 Å cumulative source windows. Record rejected steps; never silently clip forces.

Before execution, define source-independent criteria for target/spectator dissociation and loss of a retained coordinating water. A coordination-number change alone is not failure. Clear state-changing departures remain recorded outcomes, with their energies preserved, but cannot silently enter the same bound-state accommodation pool. Unchanged atom counts do not establish unchanged occupancy or hydration state.

Finite specification:

- One LBT scout: at most 12 accepted steps or 40 distinct total-energy/gradient evaluations.
- After its gates pass, complete all six searches under the same rule: at most 80 accepted steps or 200 distinct evaluations per search, including scout work.
- Convergence requires maximum free-atom force ≤0.10 kcal/mol/Å and atom-vector RMS ≤0.03 kcal/mol/Å over the entire movable system. The evaluator must numerically resolve these targets first.
- Budget exhaustion, persistent rejection or excess residual force means incomplete relaxation, not a minimum.
- Cross-score each new terminal candidate with both metals in its compatible source pool. At most six additional opposite-metal evaluations are needed beyond the six own-metal terminal evaluations. Preserve mathematical minima and the existing numerical selection policy.

Measure the scout before scheduling the remainder. Further continuation needs a concrete diagnosis and separately recorded finite budget. Compare additional work, component changes, full residual forces and source dependence with the restricted baseline. This isolates the effect of accessible motion, not the correctness of the energy model.

## Module 3 — Comparable scaffold mechanics

In parallel with Module 1 implementation, prepare a local-electronic/MM-scaffold model for LBT using its existing ff19SB topology and the same interaction policy as Mex. First compare this distinct Hamiltonian on the saved geometry pool; only then optimize it with the same full-coordinate procedure.

Include complete donor chemistry and relevant local hydrogen bonds in the electronic region. Account exactly once for electronic internal energy, environment internal energy, embedding, cross-region repulsion/dispersion, boundary bonded terms, exclusions/1–4 scaling, charge redistribution and cap-force projection. Never add full MM energy atop whole-peptide MACE.

The legacy LBT48 reference is not an embedded model. Cap47 depends on three physical anchors; its external-charge and cross-boundary interaction contract is incomplete. Build a valid source-derived partition and correct Jacobians. Qualify one small expanded region. Freeze the region/force-field rule before inspecting new preference results. If a valid boundary cannot be implemented in a contained effort, retain whole-peptide LBT as a control and record the blocker.

**Gate:** complete interaction ledger, state/element support, full-gradient consistency and admissible same-geometry comparisons. Distinguish coordinate-space effects from Hamiltonian effects throughout.

## Module 4 — Independent accuracy checks

### Electronic response

Reuse the working embedded PBE0 route for Mex:

```text
E_reference = PBE0_emb(Q;MM) + same_MM_ledger + same_fixed_FF_ALPB
```

PBE0 replaces both MACE vacuum and the GFN2 embedding difference. Two source-origin cells already exist; reuse compatible artifacts. Qualify its total derivative separately. The approximate-model finite differences do not qualify DFT.

Bound that qualification to two representative directions per metal: one MM-only direction and one coupled core/boundary direction. Central differences require eight displaced singlepoints, plus at most two central analytic-gradient evaluations before compatible reuse. Freeze step sizes and tolerances first. Failed derivatives stop the dependent reference comparison; further diagnostic calls require a specific recorded numerical diagnosis.

The finite response matrix is two old own-metal restricted endpoints and two new coupled endpoints, each for both metals: eight embedded cells before exact reuse. Start with one numerical/cost scout. Continue only after admission and tractability checks. Compare each metal's work, the differential and relevant force projections. Retain the existing 1 kcal/mol response resolution; freeze force-projection targets and numerical-uncertainty requirements before new reference outcomes. Preserve SCF and physical-state limitations.

### Scaffold stiffness

PBE0 plus the same MM environment cannot independently validate MM stiffness. Test that separately on one backbone/sidechain segment selected by the largest contribution to coupled motion, before new reference scoring. An expanded electronic region must include that segment and its local hydrogen-bond network. Evaluate before/after coordinates for both metals: four cells before compatible reuse. Freeze selection, boundaries and environment first.

Compare deformation work and restoring forces under matched physical states, not absolute totals across different partitions. If the important motion is too delocalized for a defensible small region, report the limit. Whole-peptide DFT is a possible separately measured alternative, not an automatic additional campaign. No full Hessian, functional sweep, fitted stiffness correction or repeated unchanged SCF failure is required. DFT may replace MACE as the research backend when accuracy and cost support it.

## Module 5 — Natural LanM transfer

After the mechanics and response gates pass, assess existing Mex8FNS and both Hans8DQ2/8FNR sources. Preserve both Hans results and their different water/spectator inventories. Use a consistent Hamiltonian and balanced exchange:

```text
R(P,Q) = [E_Dy - E_La]_P - [E_Dy - E_La]_Q
```

This is a conditional energy contrast, not a binding free energy. Experimental preferences are external developmental evidence, not a means of choosing sources, restraints or stopping rules. LBT and Mex assays differ; no quantitative cross-assay affinity fit is justified.

Hans experimental transfer requires appropriate assembly/occupancy. Extracted-subunit results remain conditional. Existing deposited dimers can support a separately declared partner-contribution comparison after monomer mechanics gates; no new folding campaign is required. Preserve reserved library outcomes. If a supported natural monomer-only comparator cannot be established, state the validation gap.

See [experimental target](../lanm_coupled_response_20261016/EXPERIMENTAL_TARGET.md), [DFT readiness](../lanm_coupled_response_20261016/DFT_BACKEND_READINESS.md), and [source feasibility](../lanm_coupled_response_20261016/MEX_SECOND_SOURCE_FEASIBILITY.md).

## Decisions and completion

| Finding | Decision |
| --- | --- |
| Broad response gives material work supported by references | Restricted motion was limiting; retain coupled optimization. |
| Same-pool model changes or expanded-region references materially alter the work | Identify and repair the specific electronic, boundary or MM interaction. |
| Response converges and is supported, but preferences remain source/state dependent | Motion alone is insufficient; identify the remaining occupancy, hydration, assembly or population question. |
| Searches or references fail their gates | Preserve unavailable results and issue a bounded diagnosis, without selecting favorable survivors. |

Deliver executable commands, complete result tables, common-pool comparisons, residual forces, component energies, reference errors, failures and actual CPU/GPU costs. Finish with an evidence-supported recommendation. A smooth trajectory, extra lowering or expected sign alone is insufficient. A supported negative result can close a branch.

## Execution and restart

After activation, delegate disjoint modules: evaluator/optimizer, interaction ledger/LBT partition, and independent review. Root integrates and submits chemistry; no duplicate calculations. All numerical science, including real-fixture tests, uses Slurm under current policies. Use allocation-derived ORCA memory and matching MPI ranks, without policy bypasses or new scheduler memory/time overrides. Orchestration stays on login.

Reuse warm workers and exact caches, checkpoint accepted steps, arm tested completion/health wakes, and avoid idle waiting collectors. On activation read this plan, the parent CURRENT and actual ownership/receipts, then bind the six starts and qualify the evaluator. Declare consequential settings before execution and retain original failures. Commit/push scoped progress under existing authorization; preserve unrelated work.

Implementation references: [SciPy L-BFGS-B](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html) supports analytic gradients and finite stopping limits. [OpenMM LocalEnergyMinimizer](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LocalEnergyMinimizer.html) minimizes the energy represented in its Context; an MM-only Context cannot substitute for the external composite evaluator. Pin installed versions rather than upgrade environments from documentation.
