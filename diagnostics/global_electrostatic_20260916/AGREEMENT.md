# Approved global electrostatic experiment

User approval, 2026-09-16: **“Approved.”** This followed the explicit
proposal for four partition-test endpoints and a conditional ten-endpoint
accuracy comparison.

Frozen proposal: `diagnostics/accuracy_strategy_20260915/GLOBAL_ELECTROSTATIC_PROPOSAL_20260916.md`
SHA256: `3a60009bfef4caa987fb27bc6096ff26eee0daf3b53839205e923c491ab8d5eb`

Approved model: native r2SCAN-3c vacuum core energy + endpoint MBIS/core-to-fixed-protein Coulomb + one whole-protein TABI-PB reaction-field energy. Fixed geometry, membership, protonation and waters. Dielectrics 1/78.54, zero salt, 298.15 K, prior Bondi/common-metal radii and 1.4-A probe. Surface/solver sequence is recorded before solver execution.

Stage 1: 1H4I qm33/qm36, La/Ca (four endpoints), ESP quality and solver refinement/rigid-transform/component/reduction checks, then partition tolerance 2 kcal/mol. Stage 2 runs only after physical feasibility: GGR1GLG v3, alpha1F6S/6IP9 v3, canonical1H4I/4MAE (ten endpoints). No new reference, threshold fitting, gradients, geometry search, default change, full-protein quantum calculation, or backend/model rescue is included. Technical recovery preserves the agreed scientific tasks and all attempts.

The user previously removed CPU-time, wall-time and spending stopping budgets; none are reinstated. The named task counts define the experiment. Cost and resource use are measured. Candidate products stay under `workspaces/global_electrostatic_20260916/`. Baseline, old experiments and unrelated dirty work remain untouched. Scoped implementation/tests, normal submissions under existing allocation rules, collection, reporting, vault updates and scoped commits are included in executing this approved plan.
