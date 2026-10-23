# Water interaction-strength diagnostic

**The extra parameter is useful locally; a new actual fit is required.**

Job17596 completed all four probes and32 native programs. The eight physical
midpoint checks pass, with maximum energy and force-component discrepancies about
5×10⁻⁹ in their respective kcal/mol and kcal/mol/Å units. Polarization terms
are unchanged. Only the intended metal–water-O pair well depth changed.

The existing radius/damping derivative matrix has rank4. Adding the two well-depth
columns supplies two independent directions under the same absolute singular-value
threshold0.00013718; their orthogonal singular values are0.79899 and0.33708.
Most well-depth response overlaps the old directions, but its small independent
part aligns with the remaining error.

| Local weighted residual | Loss |
|---|---:|
| Actual selected17593 model |0.419490|
| After projection onto four old parameter directions |0.419490|
| After projection including two well-depth directions |0.140344|

The66.5% reduction is an unrestricted **local linear projection**, not an achieved
fit, a bounded feasibility result or a physical accuracy guarantee. The La block
accounts for most of the projected improvement; aquo8 blocks can worsen.

Changing well depth alone by±5% worsens the actual objective in all four probes:
La−5%5.035530; La+5%5.103063; Dy−5%2.754929; Dy+5%2.763971. Each fails energy
controls. No probe is selected as a replacement model. This pattern motivates
joint radius/damping/well-depth variation, with actual evaluation of every trial.

Root declared a separate six-parameter experiment in
`diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_six_parameter_v1/ROOT_DECLARATION.md`. It retains the
same consumed data, weights and thresholds, with positive epsilon bounds0.5–1.5
original values and finite24-step/384-vector/4096-program limits. No donor transfer,
reserved validation or full-system production is admitted by this diagnostic.

## Cost and evidence

Scheduler:12s×4CPUs=48 allocated CPU-seconds, exit0:0. Adapter wall10.507s,
parentCPU4.326s, aggregate childuser8.537s/system11.063s. Exact result SHA256:
`a93f1502a82ee41f34cd9ff40999b37100de514ff3c19d07b547dfcf70396d80`.
Independent raw replay is recorded separately in execution_actual_review_v1.
