# Matched hybrid response: GGR structural transfer

Declared 2026-09-18 before new endpoint/grid outputs, under Jacob's active
MACE discriminator goal and blanket pilot authorization. Preserve production.

## Scientific question and fixed inventory

Does the actual alpha/GGR ordering improvement survive the two other consumed
GGR structures, 2FW0 and 2FVY, where direct masked MACE failed? This tests
structural robustness on the same biological pair, not a new affinity label.
Use both existing whole preparations:
`workspaces/mace_omol_20260917/ggr_source_bridge_{2fw0,2fvy}_v1/source_preparation.json`.
Source rows/evidence/exclusions are in `ggr_structure_sources_v1/`. Preserve
observed chain A, terminal chemistry, pH7 source microstates, zero waters,
explicit nonprotein exclusions, normalized H and all source heavy coordinates.
Keep the existing open/sugar-free and glucose-bound/radiation-damage caveats.
Do not infer a pure sugar or geometry effect from the difference.

For each structure generate BOTH existing source-graph policies: `alpha_caps`
(the extended 58-atom GGR core) and `connected_segment` (111 atoms). Reuse the
repaired v3 parent in each source row and the existing ggr_preparation code;
map its source hydrogens to the same normalized whole preparation. Keep caps
at their original heavy-anchor coordinates. Verify exact chemistry/counts,
formal charges (-1 Ca, 0 La), heavy coordinates, paired coordinates and source
map. Reproduce existing normalized 1GLG inputs as an implementation test.
Unsupported source chemistry remains explicit, not replaced by another core.

## Unchanged model and checks

Same energy H = DFT_vac(core) + T_mask(full) - T_mask(core); native ORCA6.1.1
r2SCAN-3c/DefGrid3/TightSCF EnGrad, exact OMOL checkpoint/mask/float64/singlet.
No solvation, geometry search, new functional, entropy or fitted coefficient.
Same physical three-coordinate metal response, Hessian(T_mask(full)), 0.20 Å
sphere, coarse/fine 0.02/0.01 Å grid, and positive curvature without clamping.
Use the exact numerical/energy/gradient/donor criteria in
MATCHED_H_RESPONSE_PLAN.md. Both response and final partition differences must
be <=2 kcal-equivalent within EACH structure. Existing 1GLG qualification
failures remain failures; this pilot cannot retrospectively qualify them.

Freeze predictions before native evaluation. Raw actual directions and qualified
scores remain separate. Reuse actual alpha1F6S/6IP9 and 1GLG from
`workspaces/mace_omol_hybrid_response_20260918/native_report_v2/result.json`.
Report all 12 alpha-minus-GGR comparisons: 2 alpha structures x 3 GGR structures
x 2 GGR representations. Raw robustness requires all >0.02; missing cases stay
in the denominator. Report eight new margins as well as the four reused ones.
No averaging away failures, best-core selection, absolute calibration or PQQ
bands. GGR remains one same-assay Ca-favoring biological group; alpha remains
one qualified La-favoring group with its assay limitations.

## Runs and cost

Initial: 8 native core DFT gradients; 8 masked core gradients, 4 whole gradients,
and 144 whole scalar grid calls (36 x 2 metals x 2 structures), total156MACE.
Compatible archived whole center energies provide independent scalar checks;
new gradients are required. No archive substitute for missing quantum gradients.
Conditional: at most8native DFT plus16matching full/core MACE gradients at the
fixed eligible points, retaining donor membership. This is two high-level
endpoints per metal state including the center; no open-ended optimization.

Existing runners: one 64CPU/4x16-rank quantum batch; GPU centers and separate
72-point grids use the existing A5000/16CPU/64474MiB policy. Previous GGR72-point
batch took1523allocatedGPU-seconds; the8native batch took614wall seconds/64CPU.
Expect development work on the scale of minutes of DFT and roughly two such
GPU grid batches; record actual time, memory and failures. No project compute
or time cap; finite inventory and scheduler policies remain. Larger benchmark
or analytic-Hessian development is outside this inventory.

## Deliverables

Source-backed preparation, frozen manifests and existing runner dispatch;
actual collection, static/response contrasts, all physical qualification fields,
measured costs, report and vault note. No default change, production rescore,
label change, push, or interference with other work.
